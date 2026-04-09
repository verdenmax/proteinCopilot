//! MCP Tool definitions — thin wrappers around library crate functions.
//!
//! Each tool is a `Result<Json<T>, ErrorData>` returning function that:
//! 1. Parses parameters
//! 2. Delegates to a library crate
//! 3. Returns structured JSON or a proper MCP error

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use rmcp::handler::server::tool::ToolRouter;
use rmcp::handler::server::wrapper::{Json, Parameters};
use rmcp::model::{ErrorCode, ServerInfo};
use rmcp::schemars;
use rmcp::{ErrorData, ServerHandler};
use serde::{Deserialize, Serialize};
use uuid::Uuid;

use protein_copilot_spectrum_io::reader::SpectrumReader;

use protein_copilot_core::ai_decision::AiDecision;
use protein_copilot_core::engine::{HealthStatus, SearchEngineAdapter};
use protein_copilot_core::search_params::SearchParams;
use protein_copilot_core::search_result::{SearchResult, SearchResultSummary};
use protein_copilot_core::spectrum::{AcquisitionMode, Spectrum, SpectrumSummary};
use protein_copilot_dia_extraction::{
    extract_dia_precursors as run_dia_extraction,
    extract_single_spectrum_precursors, DiaExtractionConfig, IsotopePatternExtractor,
    SingleSpectrumExtractionResult,
};
use protein_copilot_param_recommend::{ParamRecommender, SearchPreset, UserHints};
use protein_copilot_report::ReportGenerator;
use protein_copilot_search_engine::{SearchProgress, SimpleSearchEngine};

use protein_copilot_result_import::{
    custom_json::CustomJsonParser, converter::build_search_result,
    diann::DiannParser, scan_matcher::{ScanMatcherConfig, match_scans},
    unimod::UnimodDb, ImportFormat, ImportResult, ResultParser,
};

// ---------------------------------------------------------------------------
// Tool input types
// ---------------------------------------------------------------------------

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ReadSpectraInput {
    /// Path to the spectrum file (.mgf or .mzML)
    file_path: String,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct GetSpectrumInput {
    /// Path to the spectrum file (.mgf or .mzML)
    file_path: String,
    /// Scan number to retrieve (1-based)
    scan_number: u32,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct RecommendParamsInput {
    /// Spectrum summary (from read_spectra). If provided, uses this directly.
    summary: Option<SpectrumSummary>,
    /// Path to spectrum file. Used to generate summary if summary is not provided.
    file_path: Option<String>,
    /// Optional user hints (experiment_type, instrument_type, enzyme)
    #[serde(default, deserialize_with = "deserialize_hints")]
    #[schemars(with = "Option<UserHints>")]
    hints: Option<UserHints>,
    /// FASTA database path. If provided, sets database_path in the recommended params.
    database_path: Option<String>,
}

/// Deserialize hints from either a JSON object or a JSON string containing an object.
fn deserialize_hints<'de, D>(deserializer: D) -> Result<Option<UserHints>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    use serde::de::Error;

    let value: Option<serde_json::Value> = Option::deserialize(deserializer)?;
    match value {
        None => Ok(None),
        Some(serde_json::Value::Object(map)) => {
            let hints: UserHints =
                serde_json::from_value(serde_json::Value::Object(map)).map_err(D::Error::custom)?;
            Ok(Some(hints))
        }
        Some(serde_json::Value::String(s)) => {
            let hints: UserHints = serde_json::from_str(&s).map_err(D::Error::custom)?;
            Ok(Some(hints))
        }
        Some(other) => Err(D::Error::custom(format!(
            "hints must be an object or JSON string, got: {other}"
        ))),
    }
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct RunSearchInput {
    /// Search parameters (from recommend_params decision). If not provided, auto-recommends.
    #[serde(default, deserialize_with = "deserialize_params")]
    #[schemars(with = "Option<SearchParams>")]
    params: Option<SearchParams>,
    /// Paths to spectrum files
    input_files: Vec<String>,
    /// FASTA database path (used if params is not provided or params.database_path is placeholder)
    database_path: Option<String>,
    /// Optional user hints for auto-recommendation (used when params not provided)
    #[serde(default, deserialize_with = "deserialize_hints")]
    #[schemars(with = "Option<UserHints>")]
    hints: Option<UserHints>,
    /// Optional run_id from extract_dia_precursors. When provided, uses cached
    /// DIA-extracted spectra instead of reading from input_files.
    dia_run_id: Option<String>,
}

/// Deserialize params from either a JSON object or a JSON string.
fn deserialize_params<'de, D>(deserializer: D) -> Result<Option<SearchParams>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    use serde::de::Error;
    let value: Option<serde_json::Value> = Option::deserialize(deserializer)?;
    match value {
        None => Ok(None),
        Some(serde_json::Value::Object(map)) => {
            let p: SearchParams =
                serde_json::from_value(serde_json::Value::Object(map)).map_err(D::Error::custom)?;
            Ok(Some(p))
        }
        Some(serde_json::Value::String(s)) => {
            let p: SearchParams = serde_json::from_str(&s).map_err(D::Error::custom)?;
            Ok(Some(p))
        }
        Some(other) => Err(D::Error::custom(format!(
            "params must be an object or JSON string, got: {other}"
        ))),
    }
}

/// Deserialize MassTolerance from either a JSON object or a JSON string.
fn deserialize_tolerance<'de, D>(
    deserializer: D,
) -> Result<Option<protein_copilot_core::search_params::MassTolerance>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    use protein_copilot_core::search_params::MassTolerance;
    use serde::de::Error;
    let value: Option<serde_json::Value> = Option::deserialize(deserializer)?;
    match value {
        None => Ok(None),
        Some(serde_json::Value::Object(map)) => {
            let t: MassTolerance =
                serde_json::from_value(serde_json::Value::Object(map)).map_err(D::Error::custom)?;
            Ok(Some(t))
        }
        Some(serde_json::Value::String(s)) => {
            let t: MassTolerance = serde_json::from_str(&s).map_err(D::Error::custom)?;
            Ok(Some(t))
        }
        Some(other) => Err(D::Error::custom(format!(
            "fragment_tolerance must be an object or JSON string, got: {other}"
        ))),
    }
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct GenerateSummaryInput {
    /// Search result to summarize (provide either this or run_id)
    result: Option<SearchResult>,
    /// Run ID from a previous run_search call (server retrieves cached result)
    run_id: Option<String>,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ExportResultsInput {
    /// Search result to export (provide either this or run_id)
    result: Option<SearchResult>,
    /// Run ID from a previous run_search call (server retrieves cached result)
    run_id: Option<String>,
    /// Output directory path
    output_dir: String,
}

/// Engine status response
#[derive(Debug, Serialize, schemars::JsonSchema)]
struct EngineStatus {
    engine: protein_copilot_core::engine::EngineInfo,
    status: HealthStatus,
}

/// Response when search is started asynchronously
#[derive(Debug, Serialize, schemars::JsonSchema)]
struct SearchStarted {
    /// Run ID to use with get_search_status, generate_summary, export_results
    run_id: String,
    /// Current status
    status: String,
    /// Message for the user
    message: String,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct GetSearchStatusInput {
    /// Run ID from run_search
    run_id: String,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct CancelSearchInput {
    /// Run ID of the search to cancel.
    run_id: String,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ListSearchesInput {
    /// Filter by status prefix (e.g. "Completed", "Failed"). Optional.
    #[serde(default)]
    status_filter: Option<String>,
    /// Maximum results to return. Default 20.
    #[serde(default)]
    limit: Option<u32>,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct ListSearchesResponse {
    searches: Vec<crate::history::SearchHistoryEntry>,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct AnnotateSpectrumInput {
    /// Run ID — use to annotate an existing PSM from a search result.
    #[serde(default)]
    run_id: Option<String>,
    /// Spectrum file path — use for manual annotation mode.
    #[serde(default)]
    file_path: Option<String>,
    /// Scan number (1-based) to annotate.
    scan_number: u32,
    /// Peptide sequence — required for manual mode.
    #[serde(default)]
    peptide_sequence: Option<String>,
    /// Charge state — required for manual mode.
    #[serde(default)]
    charge: Option<i32>,
    /// Protein accession(s) — optional for manual mode (e.g. ["sp|P00001|TEST_HUMAN"]).
    #[serde(default)]
    protein_accessions: Option<Vec<String>>,
    /// Output HTML file path. Default: ./annotation_scan{N}.html
    #[serde(default)]
    output_path: Option<String>,
    /// Fragment mass tolerance. Default: 20 ppm.
    #[serde(default, deserialize_with = "deserialize_tolerance")]
    fragment_tolerance: Option<protein_copilot_core::search_params::MassTolerance>,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct AnnotateResult {
    output_path: String,
    scan_number: u32,
    peptide_sequence: String,
    charge: i32,
    score: f64,
    matched_ions: u32,
    total_ions: u32,
    delta_mass_ppm: f64,
    protein_accessions: Vec<String>,
    message: String,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ExtractDiaPrecursorsInput {
    /// Path to the spectrum file (.mzML)
    file_path: String,
    /// Output mode: "multi" (multiple precursors per spectrum) or "pseudo" (one precursor per spectrum). Default: "pseudo"
    #[serde(default = "default_output_mode")]
    output_mode: String,
    /// Minimum charge state to consider (default: 2)
    min_charge: Option<i32>,
    /// Maximum charge state to consider (default: 5)
    max_charge: Option<i32>,
    /// Override acquisition mode detection: "DDA" or "DIA". If not set, auto-detects.
    acquisition_mode: Option<String>,
}

fn default_output_mode() -> String {
    "pseudo".to_string()
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ExtractSpectrumPrecursorsInput {
    /// Path to the spectrum file (.mzML). The file is read to find both the
    /// target MS2 scan and nearby MS1 spectra for isotope pattern extraction.
    file_path: String,
    /// Scan number (1-based) of the MS2 spectrum to extract precursors for.
    scan_number: u32,
    /// Minimum charge state to consider (default: 2)
    min_charge: Option<i32>,
    /// Maximum charge state to consider (default: 5)
    max_charge: Option<i32>,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct DiaExtractionOutput {
    detected_mode: String,
    ms1_count: u32,
    ms2_count: u32,
    total_precursors_extracted: u32,
    avg_precursors_per_ms2: f64,
    charge_distribution: std::collections::HashMap<i32, u32>,
    output_spectra_count: u32,
    run_id: String,
    message: String,
}

/// Input for the `extract_xic` MCP tool.
#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ExtractXicInput {
    /// Spectrum file path (mzML only).
    #[schemars(description = "Path to the spectrum file (.mzML). XIC extraction requires mzML format for MS1+MS2 and isolation window data.")]
    file_path: Option<String>,
    /// Target scan number (1-based).
    #[schemars(description = "Scan number (1-based) to center the XIC around.")]
    scan_number: u32,
    /// Peptide sequence.
    #[schemars(description = "Peptide amino acid sequence (one-letter codes).")]
    peptide_sequence: Option<String>,
    /// Charge state.
    #[schemars(description = "Precursor charge state.")]
    charge: Option<i32>,
    /// Real precursor m/z (not DIA isolation window center).
    #[schemars(description = "True precursor m/z. For DIA data, use the PSM-derived value, not the isolation window center.")]
    precursor_mz: Option<f64>,
    /// Complete modifications list (fixed + applied variable).
    #[schemars(description = "Modifications applied to this peptide (fixed + variable). If omitted, uses unmodified sequence.")]
    modifications: Option<Vec<protein_copilot_core::search_params::Modification>>,
    /// Number of DIA cycles before/after target (default: 5).
    #[schemars(description = "Number of DIA cycles before and after target scan. Default: 5.")]
    n_cycles: Option<u32>,
    /// Number of top ions to display (default: 6).
    #[schemars(description = "Number of top fragment ions to display. Default: 6.")]
    top_n_ions: Option<usize>,
    /// Heavy-label type for SILAC comparison.
    #[schemars(description = "Heavy-label configuration. Use {\"Silac\": {\"heavy_k_delta\": 8.014199, \"heavy_r_delta\": 10.008269}} for standard SILAC.")]
    label_type: Option<protein_copilot_xic::LabelType>,
    /// m/z extraction tolerance (default: 20 ppm).
    #[schemars(description = "Mass tolerance for XIC peak extraction. Default: 20 ppm.")]
    extraction_tolerance: Option<protein_copilot_core::search_params::MassTolerance>,
    /// Intensity extraction rule (default: MaxInWindow).
    #[schemars(description = "How to extract intensity from peaks within tolerance. Default: MaxInWindow.")]
    intensity_rule: Option<protein_copilot_xic::IntensityRule>,
    /// Plotly loading mode (default: Cdn).
    #[schemars(description = "Plotly.js loading: 'Cdn' (default, smaller) or 'Embedded' (offline).")]
    plotly_mode: Option<protein_copilot_xic::PlotlyMode>,
    /// Output HTML file path.
    #[schemars(description = "Output HTML file path. Default: ./output/xic_scan{N}.html")]
    output_path: Option<String>,
    /// Run ID to resolve PSM context (single-file searches only).
    #[schemars(description = "Run ID from a previous search. Auto-fills peptide, charge, mods, precursor_mz. MVP: single-file searches only.")]
    run_id: Option<String>,
}

/// Result returned by `extract_xic`.
#[derive(Debug, Serialize, schemars::JsonSchema)]
struct ExtractXicResult {
    /// Path to the generated HTML file.
    output_path: String,
    /// Number of MS2 scans in the XIC window.
    ms2_scan_count: usize,
    /// Number of light fragment traces extracted.
    light_trace_count: usize,
    /// Number of heavy fragment traces extracted.
    heavy_trace_count: usize,
    /// Whether MS1 precursor XIC was found.
    has_ms1_xic: bool,
    /// Summary message.
    summary: String,
}

/// Presets list response
#[derive(Debug, Serialize, schemars::JsonSchema)]
struct PresetsResponse {
    presets: Vec<SearchPreset>,
}

/// Input for the import_search_results tool.
#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ImportSearchResultsInput {
    /// Path to external search result file (.json, .parquet, .spectra).
    result_file: String,
    /// Result file format. 'auto' detects from extension. Options: auto, custom_json, diann_parquet, pfind_spectra.
    #[serde(default = "default_import_format")]
    format: String,
    /// Directory containing mzML files. File association: raw_name + '.mzML'.
    mzml_dir: String,
    /// Path to unimod.xml. If not provided, uses builtin modification database (~22 common mods).
    #[serde(default)]
    unimod_path: Option<String>,
    /// RT tolerance in seconds for scan matching. Default: 30.
    #[serde(default = "default_rt_tolerance")]
    rt_tolerance_sec: f64,
    /// Q-value threshold for filtering (DIA-NN). Default: 0.01.
    #[serde(default = "default_filter_qvalue")]
    filter_qvalue: f64,
    /// Optional: only import PSMs from this specific run/raw_title.
    #[serde(default)]
    run_filter: Option<String>,
}

fn default_import_format() -> String {
    "auto".to_string()
}
fn default_rt_tolerance() -> f64 {
    30.0
}
fn default_filter_qvalue() -> f64 {
    0.01
}

/// Helper to create MCP error with suggestion from CoreError
fn mcp_core_err(err: protein_copilot_core::error::CoreError) -> ErrorData {
    let suggestion = err.suggestion().to_string();
    let message = format!("{err}\n\nSuggestion: {suggestion}");
    ErrorData::new(ErrorCode::INTERNAL_ERROR, message, None)
}

/// Helper to create MCP error from any Display error
fn mcp_err(code: ErrorCode, err: impl std::fmt::Display) -> ErrorData {
    ErrorData::new(code, err.to_string(), None)
}

/// State for a single search run.
struct RunState {
    progress: SearchProgress,
    result: Option<SearchResult>,
    handle: Option<tokio::task::JoinHandle<()>>,
}

/// Maximum number of cached runs before eviction.
const MAX_CACHE_SIZE: usize = 100;

/// Maximum number of cached DIA extraction runs before eviction.
/// Ordered DIA cache — insertion order tracked for FIFO eviction.
struct OrderedDiaCache {
    entries: HashMap<Uuid, Vec<Spectrum>>,
    order: Vec<Uuid>,
}

const MAX_DIA_CACHE_SIZE: usize = 10;

impl OrderedDiaCache {
    fn new() -> Self {
        Self {
            entries: HashMap::new(),
            order: Vec::new(),
        }
    }

    fn remove(&mut self, id: &Uuid) -> Option<Vec<Spectrum>> {
        if let Some(spectra) = self.entries.remove(id) {
            self.order.retain(|x| x != id);
            Some(spectra)
        } else {
            None
        }
    }

    fn insert(&mut self, id: Uuid, spectra: Vec<Spectrum>) {
        while self.order.len() >= MAX_DIA_CACHE_SIZE {
            if let Some(oldest) = self.order.first().copied() {
                self.order.remove(0);
                self.entries.remove(&oldest);
            }
        }
        self.entries.insert(id, spectra);
        self.order.push(id);
    }
}

/// Ordered run cache — insertion order tracked for FIFO eviction.
struct OrderedRunCache {
    map: HashMap<Uuid, RunState>,
    order: Vec<Uuid>,
}

impl OrderedRunCache {
    fn new() -> Self {
        Self {
            map: HashMap::new(),
            order: Vec::new(),
        }
    }

    fn insert(&mut self, id: Uuid, state: RunState) {
        if !self.map.contains_key(&id) {
            self.order.push(id);
        }
        self.map.insert(id, state);
    }

    fn get(&self, id: &Uuid) -> Option<&RunState> {
        self.map.get(id)
    }

    fn get_mut(&mut self, id: &Uuid) -> Option<&mut RunState> {
        self.map.get_mut(id)
    }

    #[allow(dead_code)]
    fn len(&self) -> usize {
        self.map.len()
    }

    /// Evict oldest non-running entries until under limit.
    fn evict_if_full(&mut self) {
        while self.map.len() >= MAX_CACHE_SIZE {
            let pos = self.order.iter().position(|id| {
                self.map
                    .get(id)
                    .is_none_or(|s| s.progress.status != "Running")
            });
            if let Some(i) = pos {
                let id = self.order.remove(i);
                self.map.remove(&id);
            } else {
                break; // all running, can't evict
            }
        }
    }
}

type RunCache = Arc<Mutex<OrderedRunCache>>;

/// Guard that sets progress to "Failed" if the task terminates abnormally.
struct PanicGuard {
    run_id: Uuid,
    cache: RunCache,
    start: Instant,
}

impl Drop for PanicGuard {
    fn drop(&mut self) {
        if let Ok(mut cache) = self.cache.lock() {
            if let Some(state) = cache.get_mut(&self.run_id) {
                // Only overwrite if still Running — don't clobber Cancelled or Failed
                if state.progress.status == "Running" {
                    state.progress.status = "Failed: task panicked".to_string();
                    state.progress.elapsed_sec = self.start.elapsed().as_secs_f64();
                    state.progress.progress_pct = None;
                    state.progress.stage = None;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Server struct
// ---------------------------------------------------------------------------

pub struct ProteinCopilotServer {
    tool_router: ToolRouter<Self>,
    registry: protein_copilot_search_engine::EngineRegistry,
    /// Unified cache for all search runs (progress + result in one lock).
    run_cache: RunCache,
    /// Cache of DIA-extracted spectra, keyed by run_id from extract_dia_precursors.
    dia_cache: Arc<Mutex<OrderedDiaCache>>,
    /// LRU cache of indexed spectrum readers for O(1) scan lookup.
    reader_cache: Arc<Mutex<lru::LruCache<PathBuf, Arc<dyn SpectrumReader>>>>,
}

impl ProteinCopilotServer {
    pub fn new() -> Self {
        let mut registry = protein_copilot_search_engine::EngineRegistry::new();
        registry.register(Box::new(SimpleSearchEngine::new()));
        Self {
            tool_router: Self::tool_router(),
            registry,
            run_cache: Arc::new(Mutex::new(OrderedRunCache::new())),
            dia_cache: Arc::new(Mutex::new(OrderedDiaCache::new())),
            reader_cache: Arc::new(Mutex::new(lru::LruCache::new(
                std::num::NonZeroUsize::new(8).unwrap(),
            ))),
        }
    }

    /// Get or create a cached indexed reader for the given file path.
    fn get_or_create_reader(
        &self,
        path: &Path,
    ) -> Result<Arc<dyn SpectrumReader>, ErrorData> {
        let canonical = std::fs::canonicalize(path).unwrap_or_else(|_| path.to_path_buf());

        // Check cache first
        {
            let mut cache = self.reader_cache.lock().unwrap_or_else(|e| e.into_inner());
            if let Some(reader) = cache.get(&canonical) {
                return Ok(Arc::clone(reader));
            }
        }

        // Create new indexed reader
        let reader: Arc<dyn SpectrumReader> = Arc::from(
            protein_copilot_spectrum_io::create_indexed_reader(path)
                .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?,
        );

        // Insert into cache
        {
            let mut cache = self.reader_cache.lock().unwrap_or_else(|e| e.into_inner());
            cache.put(canonical, Arc::clone(&reader));
        }

        Ok(reader)
    }

    /// Resolve a SearchResult from direct input or cached run_id.
    fn get_result(
        &self,
        direct: &Option<SearchResult>,
        run_id: &Option<String>,
    ) -> Result<SearchResult, ErrorData> {
        if let Some(r) = direct {
            return Ok(r.clone());
        }
        if let Some(id_str) = run_id {
            let id = Uuid::parse_str(id_str)
                .map_err(|_| mcp_err(ErrorCode::INVALID_PARAMS, "invalid run_id format"))?;
            let cache = self
                .run_cache
                .lock()
                .map_err(|_| mcp_err(ErrorCode::INTERNAL_ERROR, "cache lock failed"))?;
            let state = cache.get(&id).ok_or_else(|| {
                mcp_err(ErrorCode::INVALID_PARAMS, format!("run_id {id} not found — it may have been evicted from the cache (max 100 recent runs are kept)"))
            })?;
            return state.result.clone().ok_or_else(|| {
                mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    format!(
                        "search {id} not yet completed (status: {})",
                        state.progress.status
                    ),
                )
            });
        }
        Err(mcp_err(
            ErrorCode::INVALID_PARAMS,
            "provide either 'result' or 'run_id'",
        ))
    }
}

#[rmcp::tool_handler]
impl ServerHandler for ProteinCopilotServer {
    fn get_info(&self) -> ServerInfo {
        let mut info = ServerInfo::default();
        info.instructions = Some(
            "ProteinCopilot: AI-driven proteomics mass spectrometry search platform. \
             Use read_spectra to analyze spectrum files, recommend_params to get search \
             parameter suggestions, run_search to execute database search, and \
             generate_summary to interpret results."
                .into(),
        );
        info
    }
}

/// Validate that a file_path string is non-empty.
fn validate_file_path(path: &str) -> Result<(), ErrorData> {
    if path.trim().is_empty() {
        return Err(mcp_err(
            ErrorCode::INVALID_PARAMS,
            "file_path cannot be empty",
        ));
    }
    Ok(())
}

/// Validate that a scan_number is >= 1 (1-based indexing).
fn validate_scan_number(scan: u32) -> Result<(), ErrorData> {
    if scan == 0 {
        return Err(mcp_err(
            ErrorCode::INVALID_PARAMS,
            "scan_number must be >= 1 (1-based indexing)",
        ));
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Tool implementations
// ---------------------------------------------------------------------------

#[rmcp::tool_router]
impl ProteinCopilotServer {
    /// Read a mass spectrometry file and return a statistical summary.
    #[rmcp::tool(
        name = "read_spectra",
        description = "Read a mass spectrometry file (mgf/mzML) and return a statistical summary including spectrum count, m/z range, RT range, charge distribution, and median peaks per spectrum. Use this as the first step to understand input data."
    )]
    fn read_spectra(
        &self,
        Parameters(input): Parameters<ReadSpectraInput>,
    ) -> Result<Json<SpectrumSummary>, ErrorData> {
        validate_file_path(&input.file_path)?;
        let path = Path::new(&input.file_path);
        let info = protein_copilot_spectrum_io::detect_format(path)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        let reader = protein_copilot_spectrum_io::create_reader(&info);
        let summary = reader
            .read_summary(path)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        Ok(Json(summary))
    }

    /// Read a single spectrum by scan number.
    #[rmcp::tool(
        name = "get_spectrum",
        description = "Read a single spectrum from a file by scan number (1-based). Returns the spectrum with m/z array, intensity array, precursor info, and MS level."
    )]
    fn get_spectrum(
        &self,
        Parameters(input): Parameters<GetSpectrumInput>,
    ) -> Result<Json<Spectrum>, ErrorData> {
        validate_file_path(&input.file_path)?;
        validate_scan_number(input.scan_number)?;
        let path = Path::new(&input.file_path);
        let reader = self.get_or_create_reader(path)?;
        let spectrum = reader
            .read_spectrum(path, input.scan_number)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        Ok(Json(spectrum))
    }

    /// Recommend search parameters based on spectrum characteristics.
    #[rmcp::tool(
        name = "recommend_params",
        description = "Recommend search parameters based on spectrum file characteristics. Input: SpectrumSummary from read_spectra + optional UserHints (experiment_type, instrument_type, enzyme). Output: recommended SearchParams with confidence score and explanation. Note: set database_path in params to the FASTA file path."
    )]
    fn recommend_params(
        &self,
        Parameters(input): Parameters<RecommendParamsInput>,
    ) -> Result<Json<AiDecision<SearchParams>>, ErrorData> {
        // Get summary: use provided summary or read from file_path
        let summary = if let Some(s) = input.summary {
            s
        } else if let Some(ref fp) = input.file_path {
            validate_file_path(fp)?;
            let path = std::path::Path::new(fp);
            let info = protein_copilot_spectrum_io::detect_format(path)
                .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
            let reader = protein_copilot_spectrum_io::create_reader(&info);
            reader
                .read_summary(path)
                .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?
        } else {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                "provide either 'summary' or 'file_path'",
            ));
        };

        let recommender = ParamRecommender;
        let mut decision = recommender
            .recommend(&summary, input.hints.as_ref())
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

        // Inject database_path if provided (saves LLM from JSON manipulation)
        if let Some(ref db_path) = input.database_path {
            decision.decision.database_path = db_path.clone();
        }

        Ok(Json(decision))
    }

    /// List all available search parameter presets.
    #[rmcp::tool(
        name = "list_presets",
        description = "List all built-in search parameter presets (standard, phospho, TMT, SILAC, open search). Each preset includes name, description, parameters, and applicable scenarios."
    )]
    fn list_presets(&self) -> Json<PresetsResponse> {
        Json(PresetsResponse {
            presets: ParamRecommender::list_presets(),
        })
    }

    /// Execute a database search asynchronously.
    #[rmcp::tool(
        name = "run_search",
        description = "Run a proteomics database search. Returns immediately with a run_id. The search runs in the background. Call get_search_status(run_id) to check progress. When status is Completed, use generate_summary(run_id) and export_results(run_id)."
    )]
    async fn run_search(
        &self,
        Parameters(input): Parameters<RunSearchInput>,
    ) -> Result<Json<SearchStarted>, ErrorData> {
        // -------------------------------------------------------------------
        // DIA branch: use cached spectra from extract_dia_precursors
        // -------------------------------------------------------------------
        if let Some(ref run_id_str) = input.dia_run_id {
            let dia_uuid = Uuid::parse_str(run_id_str)
                .map_err(|_| mcp_err(ErrorCode::INVALID_PARAMS, "invalid dia_run_id format"))?;

            // Resolve and validate params BEFORE removing from cache,
            // so that a validation failure doesn't consume cached spectra.
            let mut params = if let Some(p) = input.params {
                p
            } else {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "params required when using dia_run_id \
                     (cannot auto-recommend without input files)",
                ));
            };

            if let Some(ref db_path) = input.database_path {
                params.database_path = db_path.clone();
            }

            params
                .validate()
                .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

            // Move spectra out of the DIA cache (frees the slot).
            // This happens after validation so a param error won't lose cached spectra.
            let dia_spectra = {
                let mut cache = self.dia_cache.lock().map_err(|_| {
                    mcp_err(ErrorCode::INTERNAL_ERROR, "DIA cache lock is poisoned")
                })?;
                cache.remove(&dia_uuid).ok_or_else(|| {
                    mcp_err(
                        ErrorCode::INVALID_PARAMS,
                        format!(
                            "dia_run_id '{}' not found in cache \
                             (may have been evicted or already used)",
                            run_id_str
                        ),
                    )
                })?
            };

            let run_id = Uuid::new_v4();

            // Evict + initialize in unified run cache
            {
                let mut cache = self.run_cache.lock().map_err(|_| {
                    mcp_err(ErrorCode::INTERNAL_ERROR, "run cache lock is poisoned")
                })?;
                cache.evict_if_full();
                cache.insert(
                    run_id,
                    RunState {
                        progress: SearchProgress {
                            run_id,
                            status: "Running".to_string(),
                            stage: None,
                            progress_pct: Some(0.0),
                            elapsed_sec: 0.0,
                            estimated_remaining_sec: None,
                        },
                        result: None,
                        handle: None,
                    },
                );
            }

            let run_cache_clone = Arc::clone(&self.run_cache);
            let engine = SimpleSearchEngine::new();
            let dia_source = vec![PathBuf::from(format!("dia:{}", run_id_str))];

            let progress_cache = Arc::clone(&self.run_cache);
            let progress_run_id = run_id;
            let on_progress: protein_copilot_core::progress::ProgressCallback =
                Box::new(move |p: SearchProgress| {
                    if let Ok(mut cache) = progress_cache.lock() {
                        if let Some(state) = cache.get_mut(&progress_run_id) {
                            if state.progress.status == "Running" {
                                state.progress.stage = p.stage;
                                state.progress.progress_pct = p.progress_pct;
                                state.progress.elapsed_sec = p.elapsed_sec;
                                state.progress.estimated_remaining_sec = p.estimated_remaining_sec;
                            }
                        }
                    }
                });

            let handle = tokio::spawn(async move {
                let start = Instant::now();

                // Panic guard: on abnormal exit, sets status to "Failed: task panicked"
                let _guard = PanicGuard {
                    run_id,
                    cache: Arc::clone(&run_cache_clone),
                    start,
                };

                let search_result = engine
                    .search_with_spectra(&params, dia_spectra, on_progress)
                    .await;
                let duration = start.elapsed().as_secs_f64();

                // Single lock — update progress + result atomically
                let (_updated, history_entry) = if let Ok(mut cache) = run_cache_clone.lock() {
                    if let Some(state) = cache.get_mut(&run_id) {
                        // If already cancelled, don't overwrite the status
                        if state.progress.status == "Cancelled" {
                            (true, None)
                        } else {
                            // Clear the JoinHandle — task is finishing
                            state.handle = None;
                            let entry = match search_result {
                                Ok(mut result) => {
                                    result.run_id = run_id;
                                    result.metadata.run_id = run_id;
                                    let entry = crate::history::SearchHistoryEntry {
                                        run_id,
                                        status: "Completed".to_string(),
                                        created_at: result.metadata.created_at,
                                        elapsed_sec: duration,
                                        engine_info: result.engine_info.clone(),
                                        input_files: dia_source.clone(),
                                        params_used: result.params_used.clone(),
                                        total_psms: Some(result.summary.total_psms),
                                        psms_at_1pct_fdr: Some(result.summary.psms_at_1pct_fdr),
                                        identification_rate: Some(
                                            result.summary.identification_rate,
                                        ),
                                        protein_groups: Some(
                                            result.summary.protein_groups_at_1pct_fdr,
                                        ),
                                    };
                                    state.result = Some(result);
                                    state.progress.status = "Completed".to_string();
                                    state.progress.stage = None;
                                    state.progress.progress_pct = Some(1.0);
                                    state.progress.elapsed_sec = duration;
                                    state.progress.estimated_remaining_sec = Some(0.0);
                                    Some(entry)
                                }
                                Err(e) => {
                                    let entry = crate::history::SearchHistoryEntry {
                                        run_id,
                                        status: format!("Failed: {e}"),
                                        created_at: chrono::Utc::now(),
                                        elapsed_sec: duration,
                                        engine_info: protein_copilot_core::engine::EngineInfo {
                                            name: "SimpleSearch".into(),
                                            version: "0.1.0".into(),
                                            supported_features: vec![],
                                        },
                                        input_files: dia_source.clone(),
                                        params_used: params.clone(),
                                        total_psms: None,
                                        psms_at_1pct_fdr: None,
                                        identification_rate: None,
                                        protein_groups: None,
                                    };
                                    state.progress.status = format!("Failed: {e}");
                                    state.progress.progress_pct = None;
                                    state.progress.elapsed_sec = duration;
                                    Some(entry)
                                }
                            };
                            (true, entry)
                        }
                    } else {
                        (false, None)
                    }
                } else {
                    tracing::error!("run cache lock poisoned after search {run_id}; result lost");
                    (false, None)
                };

                // Persist history to disk (outside the lock)
                if let Some(entry) = history_entry {
                    crate::history::save_entry(&entry);
                }

                // Always forget the guard — the search completed normally (success or error).
                // The guard should only trigger on unexpected task termination (panic/abort).
                // If the lock failed above, we logged the error; the guard trying to lock
                // again would also fail, leaving status as "Running" forever.
                std::mem::forget(_guard);
            });

            if let Ok(mut cache) = self.run_cache.lock() {
                if let Some(state) = cache.get_mut(&run_id) {
                    state.handle = Some(handle);
                }
            }

            return Ok(Json(SearchStarted {
                run_id: run_id.to_string(),
                status: "Running".to_string(),
                message: "DIA search started from cached spectra. \
                          Call get_search_status(run_id) to check progress."
                    .to_string(),
            }));
        }

        // -------------------------------------------------------------------
        // File-based branch: read spectra from input_files
        // -------------------------------------------------------------------
        if input.input_files.is_empty() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                "input_files is empty — provide at least one spectrum file path, \
                 or use dia_run_id",
            ));
        }

        let mut params = if let Some(p) = input.params {
            p
        } else {
            let first_file = input
                .input_files
                .first()
                .ok_or_else(|| mcp_err(ErrorCode::INVALID_PARAMS, "input_files is empty"))?;
            let path = Path::new(first_file);
            let info = protein_copilot_spectrum_io::detect_format(path)
                .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
            let reader = protein_copilot_spectrum_io::create_reader(&info);
            let summary = reader
                .read_summary(path)
                .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
            let decision = ParamRecommender
                .recommend(&summary, input.hints.as_ref())
                .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
            decision.decision
        };

        if let Some(ref db_path) = input.database_path {
            params.database_path = db_path.clone();
        }

        params
            .validate()
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

        let run_id = Uuid::new_v4();
        let files: Vec<PathBuf> = input.input_files.iter().map(PathBuf::from).collect();

        // Evict + initialize in unified cache
        {
            let mut cache = self
                .run_cache
                .lock()
                .map_err(|_| mcp_err(ErrorCode::INTERNAL_ERROR, "run cache lock is poisoned"))?;
            cache.evict_if_full();
            cache.insert(
                run_id,
                RunState {
                    progress: SearchProgress {
                        run_id,
                        status: "Running".to_string(),
                        stage: None,
                        progress_pct: Some(0.0),
                        elapsed_sec: 0.0,
                        estimated_remaining_sec: None,
                    },
                    result: None,
                    handle: None,
                },
            );
        }

        let run_cache_clone = Arc::clone(&self.run_cache);
        let engine = SimpleSearchEngine::new();

        // Construct progress callback that writes stage updates to the cache
        let progress_cache = Arc::clone(&self.run_cache);
        let progress_run_id = run_id;
        let on_progress: protein_copilot_core::progress::ProgressCallback =
            Box::new(move |p: SearchProgress| {
                if let Ok(mut cache) = progress_cache.lock() {
                    if let Some(state) = cache.get_mut(&progress_run_id) {
                        if state.progress.status == "Running" {
                            state.progress.stage = p.stage;
                            state.progress.progress_pct = p.progress_pct;
                            state.progress.elapsed_sec = p.elapsed_sec;
                            state.progress.estimated_remaining_sec = p.estimated_remaining_sec;
                        }
                    }
                }
            });

        let handle = tokio::spawn(async move {
            let start = Instant::now();

            // Panic guard: on abnormal exit, sets status to "Failed: task panicked"
            let _guard = PanicGuard {
                run_id,
                cache: Arc::clone(&run_cache_clone),
                start,
            };

            let search_result = engine.search(&params, &files, on_progress).await;
            let duration = start.elapsed().as_secs_f64();

            // Single lock — update progress + result atomically
            let (_updated, history_entry) = if let Ok(mut cache) = run_cache_clone.lock() {
                if let Some(state) = cache.get_mut(&run_id) {
                    // If already cancelled, don't overwrite the status
                    if state.progress.status == "Cancelled" {
                        (true, None)
                    } else {
                        // Clear the JoinHandle — task is finishing
                        state.handle = None;
                        let entry = match search_result {
                            Ok(mut result) => {
                                result.run_id = run_id;
                                result.metadata.run_id = run_id;
                                let entry = crate::history::SearchHistoryEntry {
                                    run_id,
                                    status: "Completed".to_string(),
                                    created_at: result.metadata.created_at,
                                    elapsed_sec: duration,
                                    engine_info: result.engine_info.clone(),
                                    input_files: result.metadata.input_files.clone(),
                                    params_used: result.params_used.clone(),
                                    total_psms: Some(result.summary.total_psms),
                                    psms_at_1pct_fdr: Some(result.summary.psms_at_1pct_fdr),
                                    identification_rate: Some(result.summary.identification_rate),
                                    protein_groups: Some(result.summary.protein_groups_at_1pct_fdr),
                                };
                                state.result = Some(result);
                                state.progress.status = "Completed".to_string();
                                state.progress.stage = None;
                                state.progress.progress_pct = Some(1.0);
                                state.progress.elapsed_sec = duration;
                                state.progress.estimated_remaining_sec = Some(0.0);
                                Some(entry)
                            }
                            Err(e) => {
                                let entry = crate::history::SearchHistoryEntry {
                                    run_id,
                                    status: format!("Failed: {e}"),
                                    created_at: chrono::Utc::now(),
                                    elapsed_sec: duration,
                                    engine_info: protein_copilot_core::engine::EngineInfo {
                                        name: "SimpleSearch".into(),
                                        version: "0.1.0".into(),
                                        supported_features: vec![],
                                    },
                                    input_files: files.clone(),
                                    params_used: params.clone(),
                                    total_psms: None,
                                    psms_at_1pct_fdr: None,
                                    identification_rate: None,
                                    protein_groups: None,
                                };
                                state.progress.status = format!("Failed: {e}");
                                state.progress.progress_pct = None;
                                state.progress.elapsed_sec = duration;
                                Some(entry)
                            }
                        };
                        (true, entry)
                    }
                } else {
                    (false, None)
                }
            } else {
                // Lock is poisoned — try to recover by replacing the entire Mutex.
                // Set a new cache with this run marked as failed.
                tracing::error!("run cache lock poisoned after search {run_id}; result lost");
                (false, None)
            };

            // Persist history to disk (outside the lock)
            if let Some(entry) = history_entry {
                crate::history::save_entry(&entry);
            }

            // Always forget the guard — the search completed normally (success or error).
            // The guard should only trigger on unexpected task termination (panic/abort).
            // If the lock failed above, we logged the error; the guard trying to lock
            // again would also fail, leaving status as "Running" forever.
            std::mem::forget(_guard);
        });

        // Store the JoinHandle so cancel_search can abort it
        if let Ok(mut cache) = self.run_cache.lock() {
            if let Some(state) = cache.get_mut(&run_id) {
                state.handle = Some(handle);
            }
        }

        Ok(Json(SearchStarted {
            run_id: run_id.to_string(),
            status: "Running".to_string(),
            message: "Search started. Call get_search_status(run_id) to check progress."
                .to_string(),
        }))
    }

    /// Check the status of a running search.
    #[rmcp::tool(
        name = "get_search_status",
        description = "Check the status of a search started by run_search. Returns progress percentage and elapsed time. When status is Completed, use generate_summary(run_id) to get results."
    )]
    fn get_search_status(
        &self,
        Parameters(input): Parameters<GetSearchStatusInput>,
    ) -> Result<Json<SearchProgress>, ErrorData> {
        let id = Uuid::parse_str(&input.run_id)
            .map_err(|_| mcp_err(ErrorCode::INVALID_PARAMS, "invalid run_id format"))?;
        let cache = self
            .run_cache
            .lock()
            .map_err(|_| mcp_err(ErrorCode::INTERNAL_ERROR, "cache lock failed"))?;
        let state = cache
            .get(&id)
            .ok_or_else(|| mcp_err(ErrorCode::INVALID_PARAMS, format!("run_id {id} not found — it may have been evicted from the cache (max 100 recent runs are kept)")))?;
        Ok(Json(state.progress.clone()))
    }

    /// Cancel a running search.
    #[rmcp::tool(
        name = "cancel_search",
        description = "Cancel a running search. The search task is immediately terminated and status is set to Cancelled."
    )]
    fn cancel_search(
        &self,
        Parameters(input): Parameters<CancelSearchInput>,
    ) -> Result<Json<SearchProgress>, ErrorData> {
        let id = Uuid::parse_str(&input.run_id)
            .map_err(|_| mcp_err(ErrorCode::INVALID_PARAMS, "invalid run_id format"))?;
        let mut cache = self
            .run_cache
            .lock()
            .map_err(|_| mcp_err(ErrorCode::INTERNAL_ERROR, "cache lock failed"))?;
        let state = cache
            .get_mut(&id)
            .ok_or_else(|| mcp_err(ErrorCode::INVALID_PARAMS, format!("run_id {id} not found — it may have been evicted from the cache (max 100 recent runs are kept)")))?;

        if state.progress.status != "Running" {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!("search is not running (status: {})", state.progress.status),
            ));
        }

        // Abort the tokio task
        if let Some(handle) = state.handle.take() {
            handle.abort();
        }

        state.progress.status = "Cancelled".to_string();
        state.progress.stage = Some("Cancelled by user".to_string());
        state.progress.progress_pct = None;

        Ok(Json(state.progress.clone()))
    }

    /// Check available search engines and their health status.
    #[rmcp::tool(
        name = "check_engine",
        description = "Check available search engines and their health status. Returns engine name, version, supported features, and availability."
    )]
    async fn check_engine(&self) -> Json<EngineStatus> {
        // Use first engine from registry
        let engines = self.registry.list_available();
        if let Some(info) = engines.first() {
            let status = if let Some(engine) = self.registry.get(&info.name) {
                engine
                    .health_check()
                    .await
                    .unwrap_or(HealthStatus::Unavailable {
                        reason: "health check failed".to_string(),
                    })
            } else {
                HealthStatus::Unavailable {
                    reason: "engine not found".to_string(),
                }
            };
            Json(EngineStatus {
                engine: info.clone(),
                status,
            })
        } else {
            Json(EngineStatus {
                engine: protein_copilot_core::engine::EngineInfo {
                    name: "none".to_string(),
                    version: "0.0.0".to_string(),
                    supported_features: vec![],
                },
                status: HealthStatus::Unavailable {
                    reason: "no engines registered".to_string(),
                },
            })
        }
    }

    /// Generate a statistical summary with FDR filtering.
    #[rmcp::tool(
        name = "generate_summary",
        description = "Generate a statistical summary from search results with 1% FDR filtering. Includes identification rate, median score, median delta ppm, modification and charge distributions. Use this after run_search to interpret results."
    )]
    fn generate_summary(
        &self,
        Parameters(input): Parameters<GenerateSummaryInput>,
    ) -> Result<Json<SearchResultSummary>, ErrorData> {
        let result = self.get_result(&input.result, &input.run_id)?;
        Ok(Json(ReportGenerator::generate_summary(&result)))
    }

    /// Export search results as TSV and JSON files.
    #[rmcp::tool(
        name = "export_results",
        description = "Export search results to files. Creates psm.tsv, peptide.tsv, protein.tsv, result.json, and run_metadata.json in the specified output directory."
    )]
    fn export_results(
        &self,
        Parameters(input): Parameters<ExportResultsInput>,
    ) -> Result<String, ErrorData> {
        let result = self.get_result(&input.result, &input.run_id)?;
        let output_dir = Path::new(&input.output_dir);

        ReportGenerator::export_tsv(&result, output_dir)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        ReportGenerator::export_json(&result, &output_dir.join("result.json"))
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        ReportGenerator::export_metadata(&result.metadata, &output_dir.join("run_metadata.json"))
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

        Ok(format!(
            "Exported to {}: psm.tsv, peptide.tsv, protein.tsv, result.json, run_metadata.json",
            output_dir.display()
        ))
    }

    /// List recent search runs with status and metrics.
    #[rmcp::tool(
        name = "list_searches",
        description = "List recent search runs with their status, duration, and key metrics. Includes both active searches and completed history."
    )]
    fn list_searches(
        &self,
        Parameters(input): Parameters<ListSearchesInput>,
    ) -> Json<ListSearchesResponse> {
        let limit = input.limit.unwrap_or(20) as usize;
        let mut entries = crate::history::load_all();

        // Merge active runs from in-memory cache
        if let Ok(cache) = self.run_cache.lock() {
            for (id, state) in &cache.map {
                if !entries.iter().any(|e| e.run_id == *id) {
                    entries.push(crate::history::SearchHistoryEntry {
                        run_id: *id,
                        status: state.progress.status.clone(),
                        created_at: chrono::Utc::now(),
                        elapsed_sec: state.progress.elapsed_sec,
                        engine_info: protein_copilot_core::engine::EngineInfo {
                            name: "SimpleSearch".into(),
                            version: "0.1.0".into(),
                            supported_features: vec![],
                        },
                        input_files: vec![],
                        params_used: protein_copilot_core::search_params::SearchParams {
                            enzyme: protein_copilot_core::search_params::Enzyme::Trypsin,
                            missed_cleavages: 0,
                            fixed_modifications: vec![],
                            variable_modifications: vec![],
                            precursor_tolerance:
                                protein_copilot_core::search_params::MassTolerance {
                                    value: 0.0,
                                    unit: protein_copilot_core::search_params::ToleranceUnit::Ppm,
                                },
                            fragment_tolerance:
                                protein_copilot_core::search_params::MassTolerance {
                                    value: 0.0,
                                    unit: protein_copilot_core::search_params::ToleranceUnit::Da,
                                },
                            database_path: String::new(),
                            decoy_strategy:
                                protein_copilot_core::search_params::DecoyStrategy::Reverse,
                            acquisition_mode: None,
            max_variable_modifications: 3,
            min_peptide_length: 7,
            max_peptide_length: 50,
                        },
                        total_psms: None,
                        psms_at_1pct_fdr: None,
                        identification_rate: None,
                        protein_groups: None,
                    });
                }
            }
        }

        if let Some(ref filter) = input.status_filter {
            entries.retain(|e| e.status.starts_with(filter.as_str()));
        }
        entries.sort_by(|a, b| b.created_at.cmp(&a.created_at));
        entries.truncate(limit);
        Json(ListSearchesResponse { searches: entries })
    }

    /// Annotate a single spectrum with peptide fragment ion matching.
    #[rmcp::tool(
        name = "annotate_spectrum",
        description = "Annotate a single spectrum with peptide fragment ion matching. Generates an interactive HTML file showing matched b/y ions. Two modes: (1) provide run_id + scan_number to annotate an existing PSM, or (2) provide file_path + scan_number + peptide_sequence + charge for manual annotation."
    )]
    fn annotate_spectrum(
        &self,
        Parameters(input): Parameters<AnnotateSpectrumInput>,
    ) -> Result<Json<AnnotateResult>, ErrorData> {
        use protein_copilot_core::search_params::{MassTolerance, Modification, ToleranceUnit};

        // Resolve mode and gather PSM info + spectrum file path
        let (spectrum_file, peptide_seq, charge, modifications, protein_accs) = if let Some(
            ref rid,
        ) = input.run_id
        {
            // Mode 1: from cached search result
            let result = self.get_result(&None, &Some(rid.clone()))?;
            let psm = result
                .psms
                .iter()
                .find(|p| p.spectrum_scan == input.scan_number)
                .ok_or_else(|| {
                    mcp_err(
                        ErrorCode::INVALID_PARAMS,
                        format!("no PSM found for scan {} in run {}", input.scan_number, rid),
                    )
                })?;
            let file = result
                .metadata
                .input_files
                .first()
                .ok_or_else(|| {
                    mcp_err(ErrorCode::INTERNAL_ERROR, "no input files in search result")
                })?
                .clone();
            (
                file,
                psm.peptide_sequence.clone(),
                psm.charge,
                psm.modifications.clone(),
                psm.protein_accessions.clone(),
            )
        } else if let (Some(ref fp), Some(ref pep), Some(ch)) =
            (&input.file_path, &input.peptide_sequence, input.charge)
        {
            // Mode 2: manual annotation
            validate_file_path(fp)?;
            validate_scan_number(input.scan_number)?;
            if pep.trim().is_empty() {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "peptide_sequence cannot be empty in manual annotation mode",
                ));
            }
            (
                PathBuf::from(fp),
                pep.clone(),
                ch,
                Vec::<Modification>::new(),
                input.protein_accessions.clone().unwrap_or_default(),
            )
        } else {
            return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "provide either 'run_id' (mode 1) or 'file_path' + 'peptide_sequence' + 'charge' (mode 2)",
                ));
        };

        // Read the spectrum
        let reader = self.get_or_create_reader(&spectrum_file)?;
        let spectrum = reader
            .read_spectrum(&spectrum_file, input.scan_number)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

        let frag_tol = input.fragment_tolerance.unwrap_or(MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        });

        // Perform annotation
        let annotation = protein_copilot_search_engine::annotate::annotate_spectrum(
            &spectrum,
            &peptide_seq,
            charge,
            &frag_tol,
            &modifications,
            protein_accs.clone(),
        )
        .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))?;

        // Render HTML
        let out_path = input.output_path.map(PathBuf::from).unwrap_or_else(|| {
            PathBuf::from(format!("output/annotation_scan{}.html", input.scan_number))
        });

        ReportGenerator::render_annotation(&annotation, &out_path)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

        Ok(Json(AnnotateResult {
            output_path: out_path.display().to_string(),
            scan_number: annotation.scan_number,
            peptide_sequence: annotation.peptide_sequence,
            charge: annotation.charge,
            score: annotation.score,
            matched_ions: annotation.matched_ions,
            total_ions: annotation.total_ions,
            delta_mass_ppm: annotation.delta_mass_ppm,
            protein_accessions: annotation.protein_accessions,
            message: format!(
                "Annotation saved to {}. Matched {}/{} ions (score: {:.3}). Open in browser to view.",
                out_path.display(),
                annotation.matched_ions,
                annotation.total_ions,
                annotation.score,
            ),
        }))
    }

    /// Extract candidate precursor ions from DIA mass spectrometry data.
    #[rmcp::tool(
        name = "extract_dia_precursors",
        description = "Extract candidate precursor ions from DIA mass spectrometry data. Reads mzML file, detects DIA mode from isolation window widths, extracts precursor candidates from MS1 isotope patterns, and caches enhanced spectra for use with run_search. Returns extraction statistics."
    )]
    fn extract_dia_precursors(
        &self,
        Parameters(input): Parameters<ExtractDiaPrecursorsInput>,
    ) -> Result<Json<DiaExtractionOutput>, ErrorData> {
        validate_file_path(&input.file_path)?;
        let path = Path::new(&input.file_path);
        let info = protein_copilot_spectrum_io::detect_format(path)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        let reader = protein_copilot_spectrum_io::create_reader(&info);
        let spectra = reader
            .read_all(path)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

        // Configure extractor
        let mut extractor = IsotopePatternExtractor::default();
        if let Some(min_c) = input.min_charge {
            if min_c < 1 {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "min_charge must be >= 1",
                ));
            }
            extractor.min_charge = min_c;
        }
        if let Some(max_c) = input.max_charge {
            if max_c < extractor.min_charge {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "max_charge must be >= min_charge",
                ));
            }
            extractor.max_charge = max_c;
        }

        // Configure extraction
        let acq_mode = input.acquisition_mode.as_deref().and_then(|m| match m {
            "DDA" | "dda" => Some(AcquisitionMode::DDA),
            "DIA" | "dia" => Some(AcquisitionMode::DIA),
            _ => None,
        });
        let config = DiaExtractionConfig {
            acquisition_mode: acq_mode,
            ..DiaExtractionConfig::default()
        };

        let result = run_dia_extraction(&spectra, &extractor, &config)
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?;

        let output_spectra = if input.output_mode == "multi" {
            result.enhanced_spectra.clone()
        } else {
            result.expand_to_pseudo_spectra()
        };

        let output_count = output_spectra.len() as u32;
        let run_id = Uuid::new_v4();

        // Cache for future use
        let mut cache = self
            .dia_cache
            .lock()
            .map_err(|_| mcp_err(ErrorCode::INTERNAL_ERROR, "DIA cache lock is poisoned"))?;
        cache.insert(run_id, output_spectra);

        let output = DiaExtractionOutput {
            detected_mode: format!("{}", result.detected_mode),
            ms1_count: result.stats.ms1_count,
            ms2_count: result.stats.ms2_count,
            total_precursors_extracted: result.stats.total_precursors_extracted,
            avg_precursors_per_ms2: result.stats.avg_precursors_per_ms2,
            charge_distribution: result.stats.charge_distribution,
            output_spectra_count: output_count,
            run_id: run_id.to_string(),
            message: format!(
                "DIA extraction complete. {} precursors extracted from {} MS2 spectra. \
                 Pass dia_run_id=\"{}\" to run_search to search these spectra.",
                result.stats.total_precursors_extracted, result.stats.ms2_count, run_id
            ),
        };

        Ok(Json(output))
    }

    /// Extract precursor candidates for a single MS2 spectrum from an mzML file.
    #[rmcp::tool(
        name = "extract_spectrum_precursors",
        description = "Extract candidate precursor ions from a single MS2 spectrum. Reads the mzML file, finds the target MS2 by scan number, correlates it to the nearest MS1, and runs isotope pattern analysis within the isolation window. Returns extracted precursor candidates with charge states and the correlation method used."
    )]
    fn extract_spectrum_precursors(
        &self,
        Parameters(input): Parameters<ExtractSpectrumPrecursorsInput>,
    ) -> Result<Json<SingleSpectrumExtractionResult>, ErrorData> {
        validate_file_path(&input.file_path)?;
        validate_scan_number(input.scan_number)?;
        let path = Path::new(&input.file_path);
        let info = protein_copilot_spectrum_io::detect_format(path)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        let reader = protein_copilot_spectrum_io::create_reader(&info);
        let spectra = reader
            .read_all(path)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

        let mut extractor = IsotopePatternExtractor::default();
        if let Some(min_c) = input.min_charge {
            if min_c < 1 {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "min_charge must be >= 1",
                ));
            }
            extractor.min_charge = min_c;
        }
        if let Some(max_c) = input.max_charge {
            if max_c < extractor.min_charge {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "max_charge must be >= min_charge",
                ));
            }
            extractor.max_charge = max_c;
        }

        let result = extract_single_spectrum_precursors(&spectra, input.scan_number, &extractor)
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))?;

        Ok(Json(result))
    }

    /// Extract XIC (Extracted Ion Chromatogram) for a peptide from an mzML file.
    #[rmcp::tool(
        name = "extract_xic",
        description = "Extract XIC (Extracted Ion Chromatogram) for a peptide from an mzML file. Generates an interactive HTML file with Plotly.js showing MS1 precursor and MS2 fragment ion chromatograms. Supports SILAC heavy-label comparison. Two modes: (1) provide run_id + scan_number to use PSM context, or (2) provide file_path + scan_number + peptide_sequence + charge + precursor_mz."
    )]
    fn extract_xic(
        &self,
        #[allow(unused_variables)] Parameters(input): Parameters<ExtractXicInput>,
    ) -> Result<Json<ExtractXicResult>, ErrorData> {
        use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};

        validate_scan_number(input.scan_number)?;

        // Resolve mode: run_id or manual
        let (file_path, peptide, charge, precursor_mz, modifications) =
            if let Some(ref rid) = input.run_id {
                let result = self.get_result(&None, &Some(rid.clone()))?;
                let psm = result
                    .psms
                    .iter()
                    .find(|p| p.spectrum_scan == input.scan_number)
                    .ok_or_else(|| {
                        mcp_err(
                            ErrorCode::INVALID_PARAMS,
                            format!("no PSM for scan {} in run {}", input.scan_number, rid),
                        )
                    })?;
                let file = result
                    .metadata
                    .input_files
                    .first()
                    .ok_or_else(|| {
                        mcp_err(ErrorCode::INTERNAL_ERROR, "no input files in search result")
                    })?
                    .clone();
                (
                    file,
                    psm.peptide_sequence.clone(),
                    psm.charge,
                    psm.precursor_mz,
                    psm.modifications.clone(),
                )
            } else if let (Some(ref fp), Some(ref pep), Some(ch), Some(mz)) = (
                &input.file_path,
                &input.peptide_sequence,
                input.charge,
                input.precursor_mz,
            ) {
                validate_file_path(fp)?;
                if pep.trim().is_empty() {
                    return Err(mcp_err(
                        ErrorCode::INVALID_PARAMS,
                        "peptide_sequence cannot be empty",
                    ));
                }
                (
                    PathBuf::from(fp),
                    pep.clone(),
                    ch,
                    mz,
                    input.modifications.clone().unwrap_or_default(),
                )
            } else {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "provide either 'run_id' or 'file_path' + 'peptide_sequence' + 'charge' + 'precursor_mz'",
                ));
            };

        let params = protein_copilot_xic::ExtractionParams {
            mz_tolerance: input.extraction_tolerance.unwrap_or(MassTolerance {
                value: 20.0,
                unit: ToleranceUnit::Ppm,
            }),
            n_cycles: input.n_cycles.unwrap_or(5),
            top_n_ions: input.top_n_ions.unwrap_or(6),
            label_type: input.label_type.clone(),
            intensity_rule: input
                .intensity_rule
                .unwrap_or(protein_copilot_xic::IntensityRule::MaxInWindow),
        };

        // Extract XIC
        let xic_data = protein_copilot_xic::extract::extract_xic(
            &file_path,
            input.scan_number,
            &peptide,
            charge,
            precursor_mz,
            &modifications,
            &params,
        )
        .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?;

        // Render HTML
        let out_path = input.output_path.map(PathBuf::from).unwrap_or_else(|| {
            PathBuf::from(format!("output/xic_scan{}.html", input.scan_number))
        });

        let plotly_mode = input
            .plotly_mode
            .unwrap_or(protein_copilot_xic::PlotlyMode::Cdn);

        protein_copilot_report::xic_visualize::render_xic_html(&xic_data, &out_path, plotly_mode)
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?;

        let ms2_count = xic_data
            .fragment_xic_traces
            .first()
            .map(|t| t.data_points.len())
            .unwrap_or(0);

        let summary = format!(
            "XIC extracted for {} ({}+) — {} light traces, {} heavy traces, {} MS2 scans",
            peptide,
            charge,
            xic_data.fragment_xic_traces.len(),
            xic_data.heavy_fragment_xic_traces.len(),
            ms2_count,
        );

        Ok(Json(ExtractXicResult {
            output_path: out_path.to_string_lossy().to_string(),
            ms2_scan_count: ms2_count,
            light_trace_count: xic_data.fragment_xic_traces.len(),
            heavy_trace_count: xic_data.heavy_fragment_xic_traces.len(),
            has_ms1_xic: xic_data.ms1_precursor_xic.is_some(),
            summary,
        }))
    }

    /// Import external search results and match to mzML scans.
    #[rmcp::tool(
        name = "import_search_results",
        description = "Import external search results (DIA-NN, custom JSON, pFind) and match to mzML scans. Returns a run_id for use with annotate_spectrum, extract_xic, and generate_summary."
    )]
    fn import_search_results(
        &self,
        Parameters(input): Parameters<ImportSearchResultsInput>,
    ) -> Result<Json<ImportResult>, ErrorData> {
        let start = Instant::now();
        let result_path = PathBuf::from(&input.result_file);
        let mzml_dir = PathBuf::from(&input.mzml_dir);

        if !result_path.exists() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!("result file not found: {}", result_path.display()),
            ));
        }
        if !mzml_dir.is_dir() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!("mzml_dir is not a directory: {}", mzml_dir.display()),
            ));
        }
        if input.rt_tolerance_sec < 0.0 || !input.rt_tolerance_sec.is_finite() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                "rt_tolerance_sec must be a non-negative finite number",
            ));
        }
        if !(0.0..=1.0).contains(&input.filter_qvalue) {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                "filter_qvalue must be between 0.0 and 1.0",
            ));
        }

        // Load Unimod database
        let unimod = if let Some(ref xml_path) = input.unimod_path {
            UnimodDb::from_xml(Path::new(xml_path)).map_err(|e| {
                mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    format!("unimod.xml error: {e}"),
                )
            })?
        } else {
            UnimodDb::builtin()
        };

        // Detect format
        let format = match input.format.as_str() {
            "auto" => protein_copilot_result_import::detect_format(&result_path)
                .map_err(|e| mcp_err(ErrorCode::INVALID_PARAMS, e.to_string()))?,
            "custom_json" => ImportFormat::CustomJson,
            "diann_parquet" => ImportFormat::DiannParquet,
            "pfind_spectra" => ImportFormat::PFindSpectra,
            other => {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    format!(
                        "unknown format: '{other}'. Supported: auto, custom_json, diann_parquet, pfind_spectra"
                    ),
                ));
            }
        };

        // Parse
        let mut psms = match format {
            ImportFormat::CustomJson => CustomJsonParser
                .parse(&result_path, &unimod)
                .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?,
            ImportFormat::DiannParquet => {
                let mut parser = DiannParser::new();
                parser.filter_qvalue = Some(input.filter_qvalue);
                parser.run_filter = input.run_filter.clone();
                parser
                    .parse(&result_path, &unimod)
                    .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?
            }
            ImportFormat::PFindSpectra => {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "pFind .spectra format import is not yet supported. \
                     Supported formats: custom_json, diann_parquet",
                ));
            }
        };

        if psms.is_empty() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                "no PSMs parsed from the result file (check format and filters)",
            ));
        }

        // Scan matching
        let config = ScanMatcherConfig {
            rt_tolerance_sec: input.rt_tolerance_sec,
            mzml_dir: mzml_dir.clone(),
        };

        let match_report = match_scans(&mut psms, &config, &|path| {
            protein_copilot_spectrum_io::create_indexed_reader(path).map_err(|e| {
                protein_copilot_result_import::ResultImportError::SpectrumIo(format!(
                    "failed to open {}: {e}",
                    path.display(),
                ))
            })
        })
        .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?;

        // Collect actual mzML paths used (for downstream annotate_spectrum/extract_xic)
        let mzml_files: Vec<PathBuf> = psms
            .iter()
            .map(|p| p.raw_name.as_str())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .map(|raw| {
                let mzml = mzml_dir.join(format!("{raw}.mzML"));
                if mzml.exists() {
                    mzml
                } else {
                    mzml_dir.join(format!("{raw}.mzml"))
                }
            })
            .collect();

        // Convert to SearchResult
        let format_name = match format {
            ImportFormat::CustomJson => "custom_json",
            ImportFormat::DiannParquet => "diann_parquet",
            ImportFormat::PFindSpectra => "pfind_spectra",
        };
        let (mut search_result, import_result) =
            build_search_result(&psms, match_report, format_name, mzml_files);

        // Fix metadata: set status to Completed and record duration
        let duration = start.elapsed().as_secs_f64();
        search_result.metadata.status =
            protein_copilot_core::run_metadata::RunStatus::Completed;
        search_result.metadata.duration_sec = Some(duration);

        // Store in run_cache
        let run_id = search_result.run_id;
        {
            let mut cache = self.run_cache.lock().map_err(|_| {
                mcp_err(ErrorCode::INTERNAL_ERROR, "run cache lock poisoned")
            })?;
            cache.evict_if_full();
            cache.insert(
                run_id,
                RunState {
                    progress: SearchProgress {
                        run_id,
                        status: "Completed".to_string(),
                        stage: Some("Imported".to_string()),
                        progress_pct: Some(1.0),
                        elapsed_sec: duration,
                        estimated_remaining_sec: None,
                    },
                    result: Some(search_result.clone()),
                    handle: None,
                },
            );
        }

        // Persist history entry
        let history_entry = crate::history::SearchHistoryEntry {
            run_id,
            status: "Completed".to_string(),
            created_at: search_result.metadata.created_at,
            elapsed_sec: duration,
            engine_info: search_result.engine_info.clone(),
            input_files: search_result.metadata.input_files.clone(),
            params_used: search_result.params_used.clone(),
            total_psms: Some(search_result.summary.total_psms),
            psms_at_1pct_fdr: Some(search_result.summary.psms_at_1pct_fdr),
            identification_rate: Some(search_result.summary.identification_rate),
            protein_groups: Some(search_result.summary.protein_groups_at_1pct_fdr),
        };
        crate::history::save_entry(&history_entry);

        Ok(Json(import_result))
    }
}
