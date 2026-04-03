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

use protein_copilot_core::ai_decision::AiDecision;
use protein_copilot_core::engine::{HealthStatus, SearchEngineAdapter};
use protein_copilot_core::search_params::SearchParams;
use protein_copilot_core::search_result::{SearchResult, SearchResultSummary};
use protein_copilot_core::spectrum::{AcquisitionMode, Spectrum, SpectrumSummary};
use protein_copilot_dia_extraction::{
    extract_dia_precursors as run_dia_extraction, DiaExtractionConfig, IsotopePatternExtractor,
};
use protein_copilot_param_recommend::{ParamRecommender, SearchPreset, UserHints};
use protein_copilot_report::ReportGenerator;
use protein_copilot_search_engine::{SearchProgress, SimpleSearchEngine};

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
    /// Fragment mass tolerance. Default: 0.02 Da.
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

/// Presets list response
#[derive(Debug, Serialize, schemars::JsonSchema)]
struct PresetsResponse {
    presets: Vec<SearchPreset>,
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
        }
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
                mcp_err(ErrorCode::INVALID_PARAMS, format!("run_id {id} not found"))
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
        let path = Path::new(&input.file_path);
        let info = protein_copilot_spectrum_io::detect_format(path)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        let reader = protein_copilot_spectrum_io::create_reader(&info);
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
            .ok_or_else(|| mcp_err(ErrorCode::INVALID_PARAMS, format!("run_id {id} not found")))?;
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
            .ok_or_else(|| mcp_err(ErrorCode::INVALID_PARAMS, format!("run_id {id} not found")))?;

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
        let info = protein_copilot_spectrum_io::detect_format(&spectrum_file)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        let reader = protein_copilot_spectrum_io::create_reader(&info);
        let spectrum = reader
            .read_spectrum(&spectrum_file, input.scan_number)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

        let frag_tol = input.fragment_tolerance.unwrap_or(MassTolerance {
            value: 0.02,
            unit: ToleranceUnit::Da,
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
            extractor.min_charge = min_c;
        }
        if let Some(max_c) = input.max_charge {
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
}
