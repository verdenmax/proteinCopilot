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
use protein_copilot_core::protein_group::InferenceResult;
use protein_copilot_core::search_params::SearchParams;
use protein_copilot_core::search_result::{SearchResult, SearchResultSummary};
#[cfg(test)]
use protein_copilot_core::spectrum::SpectrumRepresentation;
use protein_copilot_core::spectrum::{AcquisitionMode, Spectrum, SpectrumSummary};
use protein_copilot_dia_extraction::{
    extract_dia_precursors as run_dia_extraction, extract_single_spectrum_precursors,
    DiaExtractionConfig, IsotopePatternExtractor, SingleSpectrumExtractionResult,
};
use protein_copilot_param_recommend::{ParamRecommender, SearchPreset, UserHints};
use protein_copilot_report::ReportGenerator;
use protein_copilot_search_engine::{SearchProgress, SimpleSearchEngine};

use protein_copilot_result_import::{
    converter::build_search_result,
    custom_json::CustomJsonParser,
    diann::DiannParser,
    pfind_tsv::PFindTsvParser,
    scan_matcher::{match_scans, validate_raw_name, ScanMatcherConfig},
    unimod::UnimodDb,
    ImportFormat, ImportResult, ResultParser,
};

// ---------------------------------------------------------------------------
// Tool input types
// ---------------------------------------------------------------------------

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ReadSpectraInput {
    /// Path to the spectrum file (supports mzML, MGF, and PFB formats)
    file_path: String,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct GetSpectrumInput {
    /// Path to the spectrum file (supports mzML, MGF, and PFB formats)
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

/// Deserialize `Option<LabelType>` accepting both JSON object and JSON-in-string.
///
/// Some MCP clients serialize nested objects as JSON strings rather than proper
/// objects.  This handles both `{"Silac": {...}}` and `"{\"Silac\": {...}}"`.
fn deserialize_label_type<'de, D>(
    deserializer: D,
) -> Result<Option<protein_copilot_core::label::LabelType>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    use protein_copilot_core::label::LabelType;
    use serde::de::Error;
    let value: Option<serde_json::Value> = Option::deserialize(deserializer)?;
    match value {
        None => Ok(None),
        Some(serde_json::Value::Null) => Ok(None),
        Some(serde_json::Value::Object(map)) => {
            let t: LabelType =
                serde_json::from_value(serde_json::Value::Object(map)).map_err(D::Error::custom)?;
            Ok(Some(t))
        }
        Some(serde_json::Value::String(s)) => {
            let t: LabelType = serde_json::from_str(&s).map_err(D::Error::custom)?;
            Ok(Some(t))
        }
        Some(other) => Err(D::Error::custom(format!(
            "label_type must be an object or JSON string, got: {other}"
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

/// Input for the `infer_proteins` MCP tool.
#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct InferProteinsInput {
    /// Run ID from a previous search. Uses cached PSMs for inference.
    #[serde(default)]
    run_id: Option<String>,

    /// Direct SearchResult (alternative to run_id).
    #[serde(default)]
    result: Option<SearchResult>,

    /// Q-value threshold for filtering PSMs before inference (default: 0.01).
    #[serde(default = "default_qvalue_threshold")]
    q_value_threshold: Option<f64>,

    /// Path to FASTA database for sequence coverage calculation.
    /// If not provided, coverage is not calculated.
    #[serde(default)]
    fasta_path: Option<String>,
}

fn default_qvalue_threshold() -> Option<f64> {
    Some(0.01)
}

/// Engine status response
#[derive(Debug, Serialize, schemars::JsonSchema)]
struct EngineStatus {
    engine: protein_copilot_core::engine::EngineInfo,
    status: HealthStatus,
    /// All available engines (when multiple are registered)
    #[serde(skip_serializing_if = "Vec::is_empty")]
    all_engines: Vec<protein_copilot_core::engine::EngineInfo>,
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
struct DiagnoseSearchInput {
    /// The run_id to diagnose (from run_search or get_search_status)
    run_id: String,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct DiagnoseSearchOutput {
    /// The run_id
    run_id: String,
    /// Overall search status: "Completed", "Failed: ...", "Cancelled"
    overall_status: String,
    /// Error classification (only for failed searches)
    error_category: Option<protein_copilot_core::diagnostics::ErrorCategory>,
    /// Stage where failure occurred
    failure_stage: Option<String>,
    /// Error detail message
    error_detail: Option<String>,
    /// Per-stage metrics
    stages: Vec<protein_copilot_core::diagnostics::DiagnosticStage>,
    /// Detected quality anomalies
    anomalies: Vec<protein_copilot_core::diagnostics::SearchAnomaly>,
    /// Repair/optimization suggestions, sorted by priority
    suggestions: Vec<protein_copilot_core::diagnostics::DiagnosticSuggestion>,
    /// Total search duration in seconds
    total_elapsed_sec: f64,
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
    /// Spectrum file path — use for manual annotation mode (supports mzML, MGF, and PFB).
    #[serde(default)]
    file_path: Option<String>,
    /// Scan number (1-based) to annotate. Set to 0 with retention_time_min to auto-find the nearest scan by retention time.
    scan_number: u32,
    /// Peptide sequence — required for manual mode.
    #[serde(default)]
    peptide_sequence: Option<String>,
    /// Charge state — required for manual mode.
    #[serde(default)]
    charge: Option<i32>,
    /// Retention time in minutes — used when scan_number=0 to auto-find the nearest MS2 scan at this RT.
    #[serde(default)]
    retention_time_min: Option<f64>,
    /// Protein accession(s) — optional for manual mode (e.g. ["sp|P00001|TEST_HUMAN"]).
    #[serde(default)]
    protein_accessions: Option<Vec<String>>,
    /// Output HTML file path. Default: ./annotation_scan{N}.html
    #[serde(default)]
    output_path: Option<String>,
    /// Fragment mass tolerance as {"value": N, "unit": "Ppm"} or {"value": N, "unit": "Da"}. Default: {"value": 20.0, "unit": "Ppm"}.
    #[serde(default, deserialize_with = "deserialize_tolerance")]
    fragment_tolerance: Option<protein_copilot_core::search_params::MassTolerance>,
    /// Number of DIA cycles before/after target for XIC extraction (default: 5). XIC requires MS1+MS2 scans and retention time data (supported by mzML and PFB with isolation windows; MGF falls back to basic annotation).
    #[serde(default)]
    n_cycles: Option<u32>,
    /// Limit XIC traces to the top N most intense fragment ions (default: all ions; zero-intensity traces are always excluded).
    #[serde(default)]
    top_n_ions: Option<usize>,
    /// Heavy-label type for SILAC comparison. Standard SILAC: {"Silac": {"heavy_k_delta": 8.014199, "heavy_r_delta": 10.008269}}. Custom labeling: {"Custom": {"residue_deltas": [["K", 8.014199]]}}. Enables dual light/heavy annotation and XIC comparison. Requires XIC-compatible format (mzML or PFB).
    #[serde(default, deserialize_with = "deserialize_label_type")]
    label_type: Option<protein_copilot_xic::LabelType>,
    /// m/z extraction tolerance for XIC as {"value": N, "unit": "Ppm"} or {"value": N, "unit": "Da"}. Default: {"value": 20.0, "unit": "Ppm"}.
    #[serde(default, deserialize_with = "deserialize_tolerance")]
    extraction_tolerance: Option<protein_copilot_core::search_params::MassTolerance>,
    /// Plotly.js loading mode for generated HTML: "Cdn" (default, smaller file, requires internet) or "Embedded" (larger file, works offline).
    #[serde(default)]
    plotly_mode: Option<protein_copilot_xic::PlotlyMode>,
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

/// Single annotation item for batch annotation via `batch_annotate_spectra`.
#[derive(Debug, Clone, Deserialize, schemars::JsonSchema)]
struct BatchAnnotateItem {
    /// Scan number (1-based) to annotate.
    scan_number: u32,
    /// Peptide sequence (stripped, no modifications).
    peptide_sequence: String,
    /// Charge state (must be > 0).
    charge: i32,
    /// Optional modifications as (position, delta_mass) pairs.
    #[serde(default)]
    modifications: Option<Vec<(usize, f64)>>,
    /// Optional protein accession(s).
    #[serde(default)]
    protein_accessions: Option<Vec<String>>,
}

/// Input for the `batch_annotate_spectra` tool. Two input modes:
/// (1) pass `items` as a JSON array, or (2) pass `annotations_file` as a TSV path.
/// XIC and SILAC parameters are global — all items share the same settings.
#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct BatchAnnotateInput {
    /// Path to the spectrum file (supports mzML, MGF, and PFB).
    file_path: String,
    /// JSON array of annotation items.
    #[serde(default)]
    items: Option<Vec<BatchAnnotateItem>>,
    /// Path to a TSV (tab-separated) file with columns: scan_number, peptide_sequence, charge, and optional modifications.
    #[serde(default)]
    annotations_file: Option<String>,
    /// Output directory for HTML files (default: ./batch_annotations/).
    #[serde(default)]
    output_dir: Option<String>,
    /// Fragment mass tolerance. Default: {"value": 20.0, "unit": "Ppm"}.
    #[serde(default, deserialize_with = "deserialize_tolerance")]
    fragment_tolerance: Option<protein_copilot_core::search_params::MassTolerance>,
    /// Heavy-label type for SILAC comparison. Standard SILAC: {"Silac": {"heavy_k_delta": 8.014199, "heavy_r_delta": 10.008269}}. Requires XIC-compatible format (mzML or PFB).
    #[serde(default, deserialize_with = "deserialize_label_type")]
    label_type: Option<protein_copilot_xic::LabelType>,
    /// Number of DIA cycles before/after target for XIC extraction (default: 5).
    #[serde(default)]
    n_cycles: Option<u32>,
    /// Limit XIC traces to the top N most intense fragment ions (default: all).
    #[serde(default)]
    top_n_ions: Option<usize>,
    /// m/z extraction tolerance for XIC. Default: {"value": 20.0, "unit": "Ppm"}.
    #[serde(default, deserialize_with = "deserialize_tolerance")]
    extraction_tolerance: Option<protein_copilot_core::search_params::MassTolerance>,
    /// Plotly.js loading mode: "Cdn" (default, smaller, requires internet) or "Embedded" (larger, offline).
    #[serde(default)]
    plotly_mode: Option<protein_copilot_xic::PlotlyMode>,
}

/// Result for a single item within a batch annotation.
#[derive(Debug, Serialize, schemars::JsonSchema)]
struct BatchAnnotateItemResult {
    scan_number: u32,
    peptide_sequence: String,
    charge: i32,
    /// Output HTML path (None if annotation failed).
    output_path: Option<String>,
    matched_ions: Option<u32>,
    total_ions: Option<u32>,
    score: Option<f64>,
    delta_mass_ppm: Option<f64>,
    protein_accessions: Vec<String>,
    /// Error message if this item failed (None on success).
    error: Option<String>,
}

/// Output for the `batch_annotate_spectra` tool.
#[derive(Debug, Serialize, schemars::JsonSchema)]
struct BatchAnnotateOutput {
    total: usize,
    succeeded: usize,
    failed: usize,
    output_dir: String,
    results: Vec<BatchAnnotateItemResult>,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ExtractDiaPrecursorsInput {
    /// Path to the spectrum file (mzML recommended; requires MS1+MS2 and isolation window metadata)
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
    /// Path to the spectrum file (mzML recommended; requires MS1+MS2 scans for isotope pattern extraction)
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

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct GetDiaCacheStatusInput {
    /// The dia_run_id returned by extract_dia_precursors. Set to empty string to list all cached extractions.
    #[serde(default)]
    dia_run_id: String,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct DiaCacheStatusOutput {
    /// Whether the cached extraction exists
    exists: bool,
    /// Where the cache is stored: "memory", "disk", or "not_found"
    location: String,
    /// Number of spectra (only available if in memory)
    spectrum_count: Option<usize>,
    /// When the extraction was performed
    extracted_at: Option<String>,
}

/// Output for listing all cached DIA extraction runs.
#[derive(Debug, Serialize, schemars::JsonSchema)]
#[allow(dead_code)]
struct ListDiaCachesOutput {
    total: usize,
    cached: Vec<DiaCacheEntry>,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct DiaCacheEntry {
    dia_run_id: String,
    location: String,
    spectrum_count: Option<usize>,
    extracted_at: Option<String>,
}

/// View mode for the `extract_xic` tool.
#[derive(Debug, Clone, Copy, Deserialize, schemars::JsonSchema)]
enum XicView {
    /// Standard 2D XIC line chart (default).
    #[serde(rename = "standard")]
    Standard,
    /// 3D MS2 overview + per-scan b/y annotated spectra.
    #[serde(rename = "3d")]
    ThreeD,
}

/// Input for the `extract_xic` MCP tool.
#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ExtractXicInput {
    /// Spectrum file path (mzML only).
    #[schemars(
        description = "Path to the spectrum file (mzML recommended; XIC extraction requires MS1+MS2 and retention time data)"
    )]
    file_path: Option<String>,
    /// Target scan number (1-based).
    #[schemars(description = "Scan number (1-based) to center the XIC around.")]
    scan_number: u32,
    /// Retention time in minutes — alternative to scan_number.
    #[serde(default)]
    #[schemars(
        description = "Retention time in minutes. When scan_number is 0, auto-finds the closest MS2 scan matching this RT and precursor_mz."
    )]
    retention_time_min: Option<f64>,
    /// Peptide sequence.
    #[schemars(description = "Peptide amino acid sequence (one-letter codes).")]
    peptide_sequence: Option<String>,
    /// Charge state.
    #[schemars(description = "Precursor charge state.")]
    charge: Option<i32>,
    /// Real precursor m/z (not DIA isolation window center).
    #[schemars(
        description = "True precursor m/z. For DIA data, use the PSM-derived value, not the isolation window center."
    )]
    precursor_mz: Option<f64>,
    /// Complete modifications list (fixed + applied variable).
    #[schemars(
        description = "Modifications applied to this peptide (fixed + variable). If omitted, uses unmodified sequence."
    )]
    modifications: Option<Vec<protein_copilot_core::search_params::Modification>>,
    /// Number of DIA cycles before/after target (default: 5).
    #[schemars(description = "Number of DIA cycles before and after target scan. Default: 5.")]
    n_cycles: Option<u32>,
    /// Number of top ions to display (default: 6).
    #[schemars(
        description = "Number of top fragment ions to display. Default: all (zero-intensity excluded)."
    )]
    top_n_ions: Option<usize>,
    /// Heavy-label type for SILAC comparison.
    #[serde(default, deserialize_with = "deserialize_label_type")]
    #[schemars(
        description = "Heavy-label configuration. Use {\"Silac\": {\"heavy_k_delta\": 8.014199, \"heavy_r_delta\": 10.008269}} for standard SILAC."
    )]
    label_type: Option<protein_copilot_xic::LabelType>,
    /// m/z extraction tolerance (default: 20 ppm).
    #[schemars(description = "Mass tolerance for XIC peak extraction. Default: 20 ppm.")]
    extraction_tolerance: Option<protein_copilot_core::search_params::MassTolerance>,
    /// Intensity extraction rule (default: MaxInWindow).
    #[schemars(
        description = "How to extract intensity from peaks within tolerance. Default: MaxInWindow."
    )]
    intensity_rule: Option<protein_copilot_xic::IntensityRule>,
    /// Plotly loading mode (default: Cdn).
    #[schemars(
        description = "Plotly.js loading: 'Cdn' (default, smaller) or 'Embedded' (offline)."
    )]
    plotly_mode: Option<protein_copilot_xic::PlotlyMode>,
    /// Output HTML file path.
    #[schemars(description = "Output HTML file path. Default: ./output/xic_scan{N}.html")]
    output_path: Option<String>,
    /// Run ID to resolve PSM context (single-file searches only).
    #[schemars(
        description = "Run ID from a previous search. Auto-fills peptide, charge, mods, precursor_mz. MVP: single-file searches only."
    )]
    run_id: Option<String>,

    /// View mode: standard 2D XIC (default) or 3D MS2 annotation view.
    #[serde(default)]
    #[schemars(
        description = "View mode: 'standard' (default, 2D XIC line chart) or '3d' (3D MS2 overview RT x m/z x intensity + per-scan b/y annotated spectra with total peak counts)."
    )]
    view: Option<XicView>,

    /// (3D only) Max non-matched peaks per scan drawn in the 3D overview.
    #[serde(default)]
    #[schemars(
        description = "Only for view=3d: max non-matched peaks per scan drawn in the 3D overview (display declutter; matched b/y always kept). Default 200."
    )]
    max_peaks_per_scan_3d: Option<usize>,
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
    /// (3D mode only) Number of MS2 scans annotated in the window.
    #[serde(skip_serializing_if = "Option::is_none")]
    annotated_scan_count: Option<usize>,
    /// (3D mode only) Matched fragment ions in the target scan.
    #[serde(skip_serializing_if = "Option::is_none")]
    target_matched_ions: Option<u32>,
}

/// Presets list response
#[derive(Debug, Serialize, schemars::JsonSchema)]
struct PresetsResponse {
    presets: Vec<SearchPreset>,
}

// --- FASTA Database Management ---

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct ListDatabasesOutput {
    databases: Vec<protein_copilot_fasta_db::DatabaseStatus>,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct ExportResultsOutput {
    /// Directory where files were exported.
    output_dir: String,
    /// List of exported file names.
    files: Vec<String>,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ListDatabasesInput {
    /// Override cache directory. Default: .proteincopilot/databases/
    #[serde(default)]
    cache_dir: Option<String>,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct DownloadDatabaseInput {
    /// Database ID (e.g. "human_swissprot", "mouse_swissprot", "crap")
    database_id: String,
    /// Override cache directory. Default: .proteincopilot/databases/
    #[serde(default)]
    cache_dir: Option<String>,
    /// Force re-download even if already cached
    #[serde(default)]
    force: Option<bool>,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct GetDatabaseInfoInput {
    /// Database ID to get info for
    database_id: String,
    /// Override cache directory. Default: .proteincopilot/databases/
    #[serde(default)]
    cache_dir: Option<String>,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct PrepareSearchInput {
    /// Paths to spectrum files (supports mzML, MGF, and PFB formats)
    input_files: Vec<String>,
    /// Optional user hints (experiment_type, instrument_type, enzyme)
    #[serde(default, deserialize_with = "deserialize_hints")]
    #[schemars(with = "Option<UserHints>")]
    hints: Option<UserHints>,
    /// Target organism for auto database resolution (e.g. "human", "mouse", "E.coli", "小鼠").
    organism: Option<String>,
    /// Direct FASTA database path. Takes priority over organism auto-resolution.
    database_path: Option<String>,
    /// Search engine: "Sage" or "SimpleSearch". Default: "SimpleSearch".
    engine: Option<String>,
    /// Override cache directory for database downloads.
    #[serde(default)]
    cache_dir: Option<String>,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct PrepareSearchOutput {
    /// Recommended search parameters with real database_path filled in.
    params: SearchParams,
    /// Explanation of why these parameters were recommended.
    reasoning: String,
    /// Confidence score (0.0 to 1.0).
    confidence: f64,
    /// Alternative approaches the user might consider.
    alternatives: Vec<String>,
    /// Evidence supporting the recommendation.
    evidence: Vec<String>,
    /// Summary of the input spectra.
    spectra_summary: SpectrumSummary,
    /// Database info if auto-resolved.
    database_info: Option<PreparedDatabaseInfo>,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct PreparedDatabaseInfo {
    /// Database ID (e.g. "human_swissprot")
    id: String,
    /// Local file path
    path: String,
    /// Number of protein sequences
    protein_count: u64,
    /// Whether this was freshly downloaded
    freshly_downloaded: bool,
}

/// Input for the import_search_results tool.
#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ImportSearchResultsInput {
    /// Path to external search result file (.json, .parquet, .spectra, .tsv).
    result_file: String,
    /// Result file format. 'auto' detects from extension. Options: auto, custom_json, diann_parquet, pfind_tsv.
    #[serde(default = "default_import_format")]
    format: String,
    /// Directory containing mzML files. File association: raw_name + '.mzML'.
    mzml_dir: String,
    /// Path to unimod.xml. If not provided, uses builtin modification database (~22 common mods).
    #[serde(default)]
    unimod_path: Option<String>,
    /// RT tolerance in minutes for scan matching. Default: 0.5.
    #[serde(default = "default_rt_tolerance")]
    rt_tolerance_min: f64,
    /// Q-value threshold for filtering (DIA-NN). Default: 0.01.
    #[serde(default = "default_filter_qvalue")]
    filter_qvalue: f64,
    /// Optional: only import PSMs from this specific run/raw_title.
    #[serde(default)]
    run_filter: Option<String>,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ClassifyEntrapmentHitsInput {
    /// Path to search results file (.parquet for DIA-NN or .tsv)
    results_file: String,
    /// Result format override. Auto-detects from extension if omitted.
    format: Option<String>,
    /// Path to YAML config file defining target/trap rules
    config_file: String,
    /// Path to target FASTA database
    target_fasta: String,
    /// Output directory (default: ./output/entrapment/)
    output_dir: Option<String>,
    /// Directory containing mzML spectrum files for provenance tracing.
    /// When provided, runs fragment ion provenance analysis after classification.
    #[serde(default)]
    #[schemars(description = "Directory containing mzML files for provenance tracing (optional)")]
    pub mzml_dir: Option<String>,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct AnalyzeEntrapmentStatsInput {
    /// Path to classified TSV file (output from classify_entrapment_hits)
    classified_file: String,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct FindSimilarTargetsInput {
    /// Peptide sequence to look up
    peptide: String,
    /// Path to target FASTA database
    target_fasta: String,
    /// Maximum mismatches to consider (default: 2)
    max_mismatches: Option<u16>,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct AnnotateProvenanceInput {
    /// Path to the spectrum file (supports mzML, MGF, and PFB formats)
    file_path: String,
    /// Scan number (1-based)
    scan_number: u32,
    /// Trap peptide sequence (stripped, no modifications)
    trap_sequence: String,
    /// Target peptide sequence (stripped). Empty string if L4.
    #[serde(default)]
    target_sequence: String,
    /// Modifications as (position, delta_mass) pairs
    #[serde(default)]
    modifications: Vec<(usize, f64)>,
    /// Fragment mass tolerance in ppm (default: 20.0)
    #[serde(default = "default_frag_tol")]
    fragment_tolerance_ppm: f64,
    /// Maximum fragment charge state (default: 2)
    #[serde(default = "default_max_charge")]
    max_fragment_charge: i32,
    /// Chimera threshold for shared_ratio (default: 0.3)
    #[serde(default = "default_chimera_threshold")]
    chimera_threshold: f64,
    /// Output HTML file path (default: ./provenance_scan{N}.html)
    #[serde(default)]
    output_path: Option<String>,
}

fn default_frag_tol() -> f64 {
    20.0
}
fn default_max_charge() -> i32 {
    2
}
fn default_chimera_threshold() -> f64 {
    0.3
}

// --- Entrapment analysis output schemas ---
// rmcp requires outputSchema with root type "object", so we define typed
// output structs instead of using Json<serde_json::Value>.

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct ClassifyEntrapmentOutput {
    total_psms: usize,
    target_psms: usize,
    trap_psms: usize,
    ambiguous_psms: usize,
    level_counts: EntrapmentLevelCountsOutput,
    top_razor_families: Vec<EntrapmentRazorFamilyOutput>,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct EntrapmentLevelCountsOutput {
    l0: usize,
    l1: usize,
    l2: usize,
    l3: usize,
    l4: usize,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct EntrapmentRazorFamilyOutput {
    family: String,
    count: usize,
    example_peptide: String,
    example_trap_protein: String,
    example_target_protein: String,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct AnalyzeEntrapmentStatsOutput {
    total_classified: usize,
    level_distribution: std::collections::HashMap<String, usize>,
    delta_mass_stats: DeltaMassStats,
    top_protein_families: Vec<(String, usize)>,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct DeltaMassStats {
    count: usize,
    min: f64,
    max: f64,
    mean: f64,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct FindSimilarTargetsOutput {
    peptide: String,
    level: String,
    best_target_peptide: Option<String>,
    best_target_protein: Option<String>,
    mismatches: Option<u16>,
    delta_mass_da: Option<f64>,
    diff_positions: Option<String>,
    index_size: usize,
    substitution_type: Option<String>,
    edit_distance: Option<u32>,
    alignment_detail: Option<String>,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct AnnotateProvenanceOutput {
    /// Path to the generated HTML mirror plot.
    output_path: String,
    /// Trap peptide sequence.
    trap_sequence: String,
    /// Target peptide sequence (empty if L4).
    target_sequence: String,
    /// Number of peaks matching only trap ions.
    trap_matched_count: u32,
    /// Number of peaks matching only target ions.
    target_matched_count: u32,
    /// Number of peaks matching both trap and target ions.
    shared_count: u32,
    /// Number of peaks matching neither.
    unassigned_count: u32,
    /// shared / (trap_matched + target_matched + shared).
    shared_ratio: f64,
    /// Whether shared_ratio exceeds the chimera threshold.
    is_chimeric: bool,
    /// Total peaks in the spectrum.
    total_peaks: usize,
}

fn default_import_format() -> String {
    "auto".to_string()
}
fn default_rt_tolerance() -> f64 {
    0.5
}
fn default_filter_qvalue() -> f64 {
    0.01
}

/// Helper to create MCP error with suggestion from CoreError
fn mcp_core_err(err: protein_copilot_core::error::CoreError) -> ErrorData {
    let code = core_err_code(&err);
    let suggestion = err.suggestion().to_string();
    let message = format!("{err}\n\nSuggestion: {suggestion}");
    ErrorData::new(code, message, None)
}

/// Classify a [`CoreError`] into the right JSON-RPC error code: user-input /
/// data problems are `INVALID_PARAMS` (the caller can fix them), genuine
/// system/engine failures are `INTERNAL_ERROR`.
fn core_err_code(err: &protein_copilot_core::error::CoreError) -> ErrorCode {
    use protein_copilot_core::error::CoreError as E;
    match err {
        E::SpectrumParseError { .. }
        | E::InvalidSearchParams { .. }
        | E::FileNotFound { .. }
        | E::UnsupportedFormat { .. }
        | E::ResultParseError { .. }
        | E::ValidationError { .. } => ErrorCode::INVALID_PARAMS,
        E::SearchEngineError { .. } | E::SshConnectionError { .. } => ErrorCode::INTERNAL_ERROR,
    }
}

fn default_cache_dir(override_dir: &Option<String>) -> std::path::PathBuf {
    if let Some(ref dir) = override_dir {
        std::path::PathBuf::from(dir)
    } else {
        std::path::PathBuf::from(".proteincopilot/databases")
    }
}

/// Maps common organism names/keywords to database IDs.
fn organism_to_database_id(organism: &str) -> Option<&'static str> {
    let lower = organism.to_lowercase();
    // Check exact IDs first
    match lower.as_str() {
        "human_swissprot" => return Some("human_swissprot"),
        "mouse_swissprot" => return Some("mouse_swissprot"),
        "ecoli_swissprot" => return Some("ecoli_swissprot"),
        "yeast_swissprot" => return Some("yeast_swissprot"),
        "arabidopsis_swissprot" => return Some("arabidopsis_swissprot"),
        "crap" => return Some("crap"),
        _ => {}
    }
    // Fuzzy keyword matching
    if lower.contains("human")
        || lower.contains("人")
        || lower.contains("homo sapiens")
        || lower.contains("9606")
    {
        Some("human_swissprot")
    } else if lower.contains("mouse")
        || lower.contains("小鼠")
        || lower.contains("mus musculus")
        || lower.contains("10090")
    {
        Some("mouse_swissprot")
    } else if lower.contains("ecoli")
        || lower.contains("e.coli")
        || lower.contains("大肠杆菌")
        || lower.contains("escherichia")
    {
        Some("ecoli_swissprot")
    } else if lower.contains("yeast") || lower.contains("酵母") || lower.contains("saccharomyces")
    {
        Some("yeast_swissprot")
    } else if lower.contains("arabidopsis") || lower.contains("拟南芥") {
        Some("arabidopsis_swissprot")
    } else if lower.contains("contaminant") || lower.contains("污染") || lower.contains("crap") {
        Some("crap")
    } else {
        None
    }
}

/// Helper to create MCP error from any Display error
fn mcp_err(code: ErrorCode, err: impl std::fmt::Display) -> ErrorData {
    ErrorData::new(code, err.to_string(), None)
}

/// Base directory for default-generated output files. Honors the
/// `PROTEIN_OUTPUT_DIR` environment variable so clients can redirect where
/// HTML/TSV artifacts land; defaults to `output/` relative to the CWD.
fn output_base_dir() -> std::path::PathBuf {
    std::env::var_os("PROTEIN_OUTPUT_DIR")
        .filter(|v| !v.is_empty())
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::path::PathBuf::from("output"))
}

/// Default output path for a generated file, under [`output_base_dir`].
fn default_output_path(file_name: &str) -> std::path::PathBuf {
    output_base_dir().join(file_name)
}

/// Make a path absolute (join with the CWD if relative) so returned output
/// paths are unambiguous regardless of how the client launched the server.
fn absolutize(p: &std::path::Path) -> std::path::PathBuf {
    if p.is_absolute() {
        p.to_path_buf()
    } else {
        std::env::current_dir()
            .map(|d| d.join(p))
            .unwrap_or_else(|_| p.to_path_buf())
    }
}

/// Ensure the parent directory of `path` exists, creating it if needed.
fn ensure_parent_dir(path: &std::path::Path) -> Result<(), ErrorData> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() && !parent.exists() {
            std::fs::create_dir_all(parent).map_err(|e| {
                mcp_err(
                    ErrorCode::INTERNAL_ERROR,
                    format!(
                        "failed to create output directory {}: {e}",
                        parent.display()
                    ),
                )
            })?;
        }
    }
    Ok(())
}

/// Find the observed precursor m/z from the closest MS1 scan.
///
/// In DIA mode the mzML `selected ion m/z` is the isolation window center,
/// not the true observed precursor. This function searches MS1 scans near
/// the target RT for the highest-intensity peak within `tol_ppm` of the
/// theoretical precursor m/z.
fn find_precursor_in_ms1(
    ms1_scans: &[protein_copilot_xic::RawScan],
    target_rt_min: f64,
    theoretical_mz: f64,
    tol_ppm: f64,
) -> Option<f64> {
    if ms1_scans.is_empty() {
        return None;
    }
    // Find the MS1 scan closest to the target RT
    let closest = ms1_scans.iter().min_by(|a, b| {
        let da = (a.retention_time_min - target_rt_min).abs();
        let db = (b.retention_time_min - target_rt_min).abs();
        da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
    })?;

    let tol_da = theoretical_mz * tol_ppm * 1e-6;
    let mut best_mz = None;
    let mut best_int = 0.0_f64;

    for (&mz, &intensity) in closest.mz_array.iter().zip(closest.intensity_array.iter()) {
        if (mz - theoretical_mz).abs() <= tol_da && intensity > best_int {
            best_mz = Some(mz);
            best_int = intensity;
        }
    }
    best_mz
}

/// Maximum FASTA file size for in-memory loading (512 MB).
const MAX_FASTA_SIZE: u64 = 512 * 1024 * 1024;

/// Load FASTA sequences into a HashMap<accession, sequence>.
fn load_fasta_sequences(fasta_path: &str) -> Result<HashMap<String, String>, String> {
    let metadata = std::fs::metadata(fasta_path)
        .map_err(|e| format!("Cannot stat FASTA file {fasta_path}: {e}"))?;
    if metadata.len() > MAX_FASTA_SIZE {
        return Err(format!(
            "FASTA file too large ({:.0} MB, max {} MB): {fasta_path}",
            metadata.len() as f64 / 1_048_576.0,
            MAX_FASTA_SIZE / 1_048_576,
        ));
    }
    let content = std::fs::read_to_string(fasta_path)
        .map_err(|e| format!("Failed to read FASTA file {fasta_path}: {e}"))?;

    let mut sequences = HashMap::new();
    let mut current_accession = String::new();
    let mut current_sequence = String::new();

    for line in content.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if !current_accession.is_empty() {
                sequences.insert(current_accession.clone(), current_sequence.clone());
                current_sequence.clear();
            }
            let acc = header.split_whitespace().next().unwrap_or("").to_string();
            if acc.is_empty() {
                tracing::warn!("Skipping FASTA entry with empty accession");
                current_accession.clear();
                continue;
            }
            current_accession = acc;
        } else if !line.starts_with('#') && !line.starts_with(';') {
            current_sequence.push_str(line.trim());
        }
    }
    if !current_accession.is_empty() {
        sequences.insert(current_accession, current_sequence);
    }

    Ok(sequences)
}

/// State for a single search run.
struct RunState {
    progress: SearchProgress,
    result: Option<SearchResult>,
    handle: Option<tokio::task::JoinHandle<()>>,
    diagnostics: Option<protein_copilot_core::diagnostics::SearchDiagnostics>,
    params_used: Option<protein_copilot_core::search_params::SearchParams>,
}

/// Maximum number of cached runs before eviction.
const MAX_CACHE_SIZE: usize = 100;

/// Parse a TSV/CSV annotations file into a list of `BatchAnnotateItem`.
///
/// Supports tab-separated files with columns:
///   scan_number, peptide_sequence, charge, modifications (optional)
///
/// - First line is treated as a header if it starts with a non-numeric character.
/// - The modifications column uses format: `pos:delta;pos:delta` (e.g. `0:42.011;5:79.966`).
fn parse_annotations_file(path: &str) -> Result<Vec<BatchAnnotateItem>, ErrorData> {
    let content = std::fs::read_to_string(path).map_err(|e| {
        mcp_err(
            ErrorCode::INVALID_PARAMS,
            format!("cannot read annotations file '{path}': {e}"),
        )
    })?;
    let lines: Vec<&str> = content.lines().filter(|l| !l.trim().is_empty()).collect();
    if lines.is_empty() {
        return Err(mcp_err(
            ErrorCode::INVALID_PARAMS,
            "annotations file is empty",
        ));
    }

    // Auto-detect header: if first char is non-numeric, skip line 0 as header
    let start = if lines[0]
        .chars()
        .next()
        .map(|c| !c.is_ascii_digit())
        .unwrap_or(false)
    {
        1
    } else {
        0
    };
    if start >= lines.len() {
        return Err(mcp_err(
            ErrorCode::INVALID_PARAMS,
            "annotations file has header only, no data rows",
        ));
    }

    let mut items = Vec::new();
    for (i, line) in lines[start..].iter().enumerate() {
        let row_num = start + i + 1; // 1-based for error messages
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!("line {row_num}: expected at least 3 tab-separated fields (scan_number, peptide_sequence, charge), got {}", fields.len()),
            ));
        }
        let scan_number: u32 = fields[0].trim().parse().map_err(|e| {
            mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!("line {row_num}: invalid scan_number '{}': {e}", fields[0]),
            )
        })?;
        let peptide_sequence = fields[1].trim().to_string();
        if peptide_sequence.is_empty() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!("line {row_num}: peptide_sequence is empty"),
            ));
        }
        let charge: i32 = fields[2].trim().parse().map_err(|e| {
            mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!("line {row_num}: invalid charge '{}': {e}", fields[2]),
            )
        })?;
        if charge <= 0 {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!("line {row_num}: charge must be > 0, got {charge}"),
            ));
        }

        let modifications = if fields.len() >= 4 && !fields[3].trim().is_empty() {
            let mods_str = fields[3].trim();
            let mut mods = Vec::new();
            for pair in mods_str.split(';') {
                let pair = pair.trim();
                if pair.is_empty() {
                    continue;
                }
                let parts: Vec<&str> = pair.splitn(2, ':').collect();
                if parts.len() != 2 {
                    return Err(mcp_err(
                        ErrorCode::INVALID_PARAMS,
                        format!(
                            "line {row_num}: invalid modification '{}' (expected pos:delta)",
                            pair
                        ),
                    ));
                }
                let pos: usize = parts[0].trim().parse().map_err(|e| {
                    mcp_err(
                        ErrorCode::INVALID_PARAMS,
                        format!(
                            "line {row_num}: invalid modification position '{}': {e}",
                            parts[0]
                        ),
                    )
                })?;
                let delta: f64 = parts[1].trim().parse().map_err(|e| {
                    mcp_err(
                        ErrorCode::INVALID_PARAMS,
                        format!(
                            "line {row_num}: invalid modification delta '{}': {e}",
                            parts[1]
                        ),
                    )
                })?;
                mods.push((pos, delta));
            }
            Some(mods)
        } else {
            None
        };

        items.push(BatchAnnotateItem {
            scan_number,
            peptide_sequence,
            charge,
            modifications,
            protein_accessions: None,
        });
    }
    if items.is_empty() {
        return Err(mcp_err(
            ErrorCode::INVALID_PARAMS,
            "no valid annotation items found in file",
        ));
    }
    Ok(items)
}

/// Default RT tolerance (minutes) for auto-scanning MS2 lookup.
const RT_AUTO_LOOKUP_TOLERANCE_MIN: f64 = 0.5;

/// DIA isolation window detection threshold (Da).
/// Spectra with median isolation window wider than this are classified as DIA.
const DIA_ISOLATION_WINDOW_THRESHOLD_DA: f64 = 5.0;

/// Default FDR threshold (1%) for protein inference filtering.
const FDR_1PCT_THRESHOLD: f64 = 0.01;

/// Maximum number of cached DIA extraction runs before eviction.
/// Ordered DIA cache — insertion order tracked for FIFO eviction.
/// When entries exceed `MAX_DIA_CACHE_SIZE` in memory, the oldest are spilled to
/// disk under `.proteincopilot/dia_cache/` and can still be retrieved.
struct OrderedDiaCache {
    entries: HashMap<Uuid, Vec<Spectrum>>,
    order: Vec<Uuid>,
    spill_dir: PathBuf,
    extracted_at: HashMap<Uuid, chrono::DateTime<chrono::Utc>>,
}

const MAX_DIA_CACHE_SIZE: usize = 10;

impl OrderedDiaCache {
    fn new() -> Self {
        let spill_dir = PathBuf::from(".proteincopilot/dia_cache");
        Self {
            entries: HashMap::new(),
            order: Vec::new(),
            spill_dir,
            extracted_at: HashMap::new(),
        }
    }

    fn remove(&mut self, id: &Uuid) -> Option<Vec<Spectrum>> {
        if let Some(spectra) = self.entries.remove(id) {
            self.order.retain(|x| x != id);
            self.extracted_at.remove(id);
            return Some(spectra);
        }
        let path = self.spill_dir.join(format!("{}.bin", id));
        if path.exists() {
            match std::fs::read(&path) {
                Ok(data) => {
                    let _ = std::fs::remove_file(&path);
                    self.extracted_at.remove(id);
                    match bincode::deserialize(&data) {
                        Ok(spectra) => return Some(spectra),
                        Err(e) => {
                            tracing::warn!("Failed to deserialize DIA cache {}: {}", id, e);
                            return None;
                        }
                    }
                }
                Err(e) => {
                    tracing::warn!("Failed to read DIA cache file {}: {}", id, e);
                    return None;
                }
            }
        }
        None
    }

    fn insert(&mut self, id: Uuid, spectra: Vec<Spectrum>) {
        // Deduplicate: if this UUID already exists, remove old position
        if self.entries.contains_key(&id) {
            self.order.retain(|x| x != &id);
        }

        while self.order.len() >= MAX_DIA_CACHE_SIZE {
            if let Some(oldest) = self.order.first().copied() {
                self.order.remove(0);
                if let Some(old_spectra) = self.entries.remove(&oldest) {
                    if !self.spill_to_disk(oldest, &old_spectra) {
                        // Spill failed — keep entry in memory to avoid data loss
                        tracing::warn!("Keeping DIA cache {} in memory (spill failed)", oldest);
                        self.entries.insert(oldest, old_spectra);
                        self.order.insert(0, oldest);
                        break;
                    }
                }
            }
        }
        self.extracted_at.insert(id, chrono::Utc::now());
        self.entries.insert(id, spectra);
        self.order.push(id);
    }

    /// Spill spectra to disk. Returns true on success, false on failure.
    fn spill_to_disk(&self, id: Uuid, spectra: &[Spectrum]) -> bool {
        if let Err(e) = std::fs::create_dir_all(&self.spill_dir) {
            tracing::error!(
                dir = %self.spill_dir.display(),
                error = %e,
                "Cannot create DIA spill directory — all spectra kept in memory. \
                 Large DIA runs may OOM. Ensure directory is writable."
            );
            return false;
        }
        let path = self.spill_dir.join(format!("{}.bin", id));
        match bincode::serialize(spectra) {
            Ok(data) => match std::fs::write(&path, &data) {
                Ok(()) => true,
                Err(e) => {
                    tracing::warn!("Failed to write DIA cache to disk {}: {}", id, e);
                    false
                }
            },
            Err(e) => {
                tracing::warn!("Failed to serialize DIA cache {}: {}", id, e);
                false
            }
        }
    }

    fn status(&self, id: &Uuid) -> DiaCacheLocation {
        if self.entries.contains_key(id) {
            let count = self.entries[id].len();
            let ts = self.extracted_at.get(id).copied();
            DiaCacheLocation::Memory {
                spectrum_count: count,
                extracted_at: ts,
            }
        } else {
            let path = self.spill_dir.join(format!("{}.bin", id));
            if path.exists() {
                let ts = self.extracted_at.get(id).copied();
                DiaCacheLocation::Disk { extracted_at: ts }
            } else {
                DiaCacheLocation::NotFound
            }
        }
    }

    /// List all cached DIA extraction run IDs with their status.
    fn list_ids(&self) -> Vec<(Uuid, DiaCacheLocation)> {
        self.order
            .iter()
            .map(|id| {
                let location = self.status(id);
                (*id, location)
            })
            .collect()
    }
}

enum DiaCacheLocation {
    Memory {
        spectrum_count: usize,
        extracted_at: Option<chrono::DateTime<chrono::Utc>>,
    },
    Disk {
        extracted_at: Option<chrono::DateTime<chrono::Utc>>,
    },
    NotFound,
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
        registry.register(Box::new(
            protein_copilot_search_engine::adapters::sage::SageAdapter::default(),
        ));
        Self {
            tool_router: Self::tool_router(),
            registry,
            run_cache: Arc::new(Mutex::new(OrderedRunCache::new())),
            dia_cache: Arc::new(Mutex::new(OrderedDiaCache::new())),
            reader_cache: Arc::new(Mutex::new(lru::LruCache::new(
                // SAFETY: 8 is a compile-time constant, always non-zero
                std::num::NonZeroUsize::new(8).expect("reader cache size is hardcoded to 8"),
            ))),
        }
    }

    /// Returns the full tool catalog (name, description, input/output JSON
    /// Schema) without starting the stdio server. Used by the `--list-tools`
    /// CLI mode and by documentation generators.
    pub fn list_tools(&self) -> Vec<rmcp::model::Tool> {
        self.tool_router.list_all()
    }

    /// Get or create a cached indexed reader for the given file path.
    fn get_or_create_reader(&self, path: &Path) -> Result<Arc<dyn SpectrumReader>, ErrorData> {
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

    /// Validate an engine name against the registry and construct its adapter.
    ///
    /// This has no side effects, so it can be called BEFORE any cache mutation
    /// to ensure an invalid engine name fails fast without consuming cached
    /// spectra or inserting an orphaned run entry.
    fn resolve_engine(&self, engine_name: &str) -> Result<Box<dyn SearchEngineAdapter>, ErrorData> {
        if self.registry.get(engine_name).is_none()
            && !engine_name.eq_ignore_ascii_case("sage")
            && !engine_name.eq_ignore_ascii_case("simplesearch")
        {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!(
                    "Engine '{}' not registered. Available: {:?}",
                    engine_name,
                    self.registry
                        .list_available()
                        .iter()
                        .map(|e| &e.name)
                        .collect::<Vec<_>>()
                ),
            ));
        }
        Ok(if engine_name.eq_ignore_ascii_case("sage") {
            Box::new(protein_copilot_search_engine::adapters::sage::SageAdapter::default())
        } else {
            Box::new(SimpleSearchEngine::new())
        })
    }

    /// Resolve a SearchResult from direct input or cached run_id.
    fn get_result(
        &self,
        direct: &Option<SearchResult>,
        run_id: &Option<String>,
    ) -> Result<SearchResult, ErrorData> {
        // If both provided, reject ambiguity
        if direct.is_some() && run_id.is_some() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                "provide either 'result' or 'run_id', not both",
            ));
        }
        if let Some(r) = direct {
            return Ok(r.clone());
        }
        if let Some(id_str) = run_id {
            let id = Uuid::parse_str(id_str)
                .map_err(|_| mcp_err(ErrorCode::INVALID_PARAMS, "invalid run_id format"))?;
            let cache = self.run_cache.lock().unwrap_or_else(|e| e.into_inner());
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

    /// Core annotation logic shared by `annotate_spectrum` and `batch_annotate_spectra`.
    /// Takes already-resolved parameters (spectrum file, scan, peptide, etc.) and
    /// performs: read spectrum → fragment matching → SILAC heavy (optional) →
    /// XIC extraction (optional) → HTML render.
    #[allow(clippy::too_many_arguments)]
    fn annotate_one_spectrum(
        &self,
        spectrum_file: &Path,
        scan_number: u32,
        peptide_seq: &str,
        charge: i32,
        modifications: &[protein_copilot_core::search_params::Modification],
        protein_accs: Vec<String>,
        frag_tol: &protein_copilot_core::search_params::MassTolerance,
        label_type: Option<&protein_copilot_core::label::LabelType>,
        n_cycles: u32,
        top_n_ions: usize,
        extraction_tolerance: &protein_copilot_core::search_params::MassTolerance,
        plotly_mode: protein_copilot_xic::PlotlyMode,
        output_path: &Path,
    ) -> Result<AnnotateResult, ErrorData> {
        // Read the spectrum
        let reader = self.get_or_create_reader(spectrum_file)?;
        let spectrum = reader
            .read_spectrum(spectrum_file, scan_number)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

        // Perform fragment ion matching
        let annotation = protein_copilot_search_engine::annotate::annotate_spectrum(
            &spectrum,
            peptide_seq,
            charge,
            frag_tol,
            modifications,
            protein_accs.clone(),
            false,
            false,
        )
        .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))?;

        let mut annotation = annotation;
        annotation.source_file = spectrum_file
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_default();

        // ── SILAC: find and annotate heavy scan ──
        if let Some(label) = label_type {
            let heavy_delta = protein_copilot_core::label::total_heavy_delta(peptide_seq, label);

            if heavy_delta.abs() > 1e-6 {
                let is_dia = spectrum
                    .precursors
                    .first()
                    .and_then(|p| p.isolation_window.as_ref())
                    .map(|w| (w.lower_offset + w.upper_offset) > 1.0)
                    .unwrap_or(false);

                let heavy_prec_mz = protein_copilot_core::label::compute_heavy_precursor_mz(
                    annotation.theoretical_mz,
                    charge,
                    peptide_seq,
                    label,
                );

                let heavy_scan_result = if is_dia {
                    reader
                        .find_by_rt(
                            spectrum_file,
                            spectrum.retention_time_min,
                            heavy_prec_mz,
                            0.5,
                        )
                        .unwrap_or(None)
                } else {
                    let candidates = reader
                        .find_by_rt(
                            spectrum_file,
                            spectrum.retention_time_min,
                            heavy_prec_mz,
                            0.5,
                        )
                        .unwrap_or(None);
                    if let Some((scan, _delta)) = candidates {
                        match reader.read_spectrum(spectrum_file, scan) {
                            Ok(spec) => {
                                let prec_mz = spec.precursors.first().map(|p| p.mz).unwrap_or(0.0);
                                let ppm_err =
                                    ((prec_mz - heavy_prec_mz) / heavy_prec_mz * 1e6).abs();
                                if ppm_err < 20.0 {
                                    Some((scan, 0.0))
                                } else {
                                    None
                                }
                            }
                            Err(_) => None,
                        }
                    } else {
                        None
                    }
                };

                if let Some((heavy_scan_num, _)) = heavy_scan_result {
                    if let Ok(heavy_spectrum) = reader.read_spectrum(spectrum_file, heavy_scan_num)
                    {
                        match protein_copilot_search_engine::annotate::annotate_heavy_spectrum(
                            &heavy_spectrum,
                            peptide_seq,
                            charge,
                            frag_tol,
                            modifications,
                            label,
                            false,
                            false,
                        ) {
                            Ok(heavy_ann) => {
                                tracing::info!(
                                    heavy_scan = heavy_scan_num,
                                    matched = heavy_ann.matched_ions,
                                    total = heavy_ann.total_ions,
                                    is_dia = is_dia,
                                    "Heavy annotation complete"
                                );
                                annotation.heavy_annotation = Some(heavy_ann);
                            }
                            Err(e) => {
                                tracing::warn!("Heavy annotation failed: {e}");
                            }
                        }
                    }
                } else {
                    tracing::info!(
                        heavy_prec_mz = format!("{:.4}", heavy_prec_mz),
                        is_dia = is_dia,
                        "No heavy scan found"
                    );
                }
            } else {
                tracing::info!(peptide = peptide_seq, "Skipping heavy: no K/R residues");
            }
        }

        let source_file = spectrum_file
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| spectrum_file.display().to_string());

        let annotation_theo_mz = annotation.theoretical_mz;

        // ── XIC extraction (mzML / PFB) or basic render (MGF) ──
        let is_xic_format = spectrum_file
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.eq_ignore_ascii_case("mzml") || e.eq_ignore_ascii_case("pfb"))
            .unwrap_or(false);

        let mut render_mode = "annotation";

        if is_xic_format {
            let xic_params = protein_copilot_xic::ExtractionParams {
                mz_tolerance: extraction_tolerance.clone(),
                n_cycles,
                top_n_ions,
                label_type: label_type.cloned(),
                intensity_rule: protein_copilot_xic::IntensityRule::MaxInWindow,
            };

            let cached_reader = self.get_or_create_reader(spectrum_file)?;

            match protein_copilot_xic::extract::extract_xic_unified(
                cached_reader.as_ref(),
                spectrum_file,
                scan_number,
                peptide_seq,
                charge,
                annotation.theoretical_mz,
                modifications,
                &xic_params,
                20.0,
            ) {
                Ok(unified_result) => {
                    let xic_meaningful = unified_result
                        .xic_data
                        .fragment_xic_traces
                        .first()
                        .map(|t| t.data_points.len() > 1)
                        .unwrap_or(false);

                    if xic_meaningful {
                        if let Some(observed) = find_precursor_in_ms1(
                            &unified_result.raw_scans.ms1_scans,
                            annotation.retention_time_min,
                            annotation.theoretical_mz,
                            20.0,
                        ) {
                            annotation.precursor_mz = observed;
                            annotation.delta_mass_ppm = (observed - annotation.theoretical_mz)
                                / annotation.theoretical_mz
                                * 1e6;
                        }

                        if let Some(ref mut ha) = annotation.heavy_annotation {
                            let heavy_theo = ha.precursor_mz;
                            if let Some(obs_heavy) = find_precursor_in_ms1(
                                &unified_result.raw_scans.ms1_scans,
                                annotation.retention_time_min,
                                heavy_theo,
                                20.0,
                            ) {
                                ha.delta_mass_ppm =
                                    Some((obs_heavy - heavy_theo) / heavy_theo * 1e6);
                            } else {
                                ha.delta_mass_ppm = None;
                            }
                        }

                        let peptide_info = protein_copilot_report::unified_types::PeptideInfo {
                            sequence: peptide_seq.to_string(),
                            charge,
                            precursor_mz: annotation_theo_mz,
                            total_k: peptide_seq.chars().filter(|&c| c == 'K').count() as u32,
                            total_r: peptide_seq.chars().filter(|&c| c == 'R').count() as u32,
                        };
                        let unified_data = protein_copilot_report::unified_types::UnifiedViewData {
                            source_file: source_file.clone(),
                            annotation: annotation.clone(),
                            xic: Some(unified_result.xic_data),
                            raw_scans: Some(unified_result.raw_scans),
                            ion_metadata: unified_result.ion_metadata,
                            peptide_info,
                        };
                        ReportGenerator::render_unified(&unified_data, output_path, plotly_mode)
                            .map_err(|e| {
                                mcp_core_err(protein_copilot_core::error::CoreError::from(e))
                            })?;
                        render_mode = "unified+xic";
                    } else {
                        let peptide_info = protein_copilot_report::unified_types::PeptideInfo {
                            sequence: peptide_seq.to_string(),
                            charge,
                            precursor_mz: annotation_theo_mz,
                            total_k: peptide_seq.chars().filter(|&c| c == 'K').count() as u32,
                            total_r: peptide_seq.chars().filter(|&c| c == 'R').count() as u32,
                        };
                        let unified_data = protein_copilot_report::unified_types::UnifiedViewData {
                            source_file: source_file.clone(),
                            annotation: annotation.clone(),
                            xic: None,
                            raw_scans: None,
                            ion_metadata: vec![],
                            peptide_info,
                        };
                        ReportGenerator::render_unified(&unified_data, output_path, plotly_mode)
                            .map_err(|e| {
                                mcp_core_err(protein_copilot_core::error::CoreError::from(e))
                            })?;
                        render_mode = "unified";
                    }
                }
                Err(_) => {
                    let peptide_info = protein_copilot_report::unified_types::PeptideInfo {
                        sequence: peptide_seq.to_string(),
                        charge,
                        precursor_mz: annotation_theo_mz,
                        total_k: peptide_seq.chars().filter(|&c| c == 'K').count() as u32,
                        total_r: peptide_seq.chars().filter(|&c| c == 'R').count() as u32,
                    };
                    let unified_data = protein_copilot_report::unified_types::UnifiedViewData {
                        source_file: source_file.clone(),
                        annotation: annotation.clone(),
                        xic: None,
                        raw_scans: None,
                        ion_metadata: vec![],
                        peptide_info,
                    };
                    ReportGenerator::render_unified(&unified_data, output_path, plotly_mode)
                        .map_err(|e| {
                            mcp_core_err(protein_copilot_core::error::CoreError::from(e))
                        })?;
                    render_mode = "unified";
                }
            }
        } else {
            // Non-mzML/PFB file (e.g. MGF) — basic annotation
            ReportGenerator::render_annotation(&annotation, output_path)
                .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        }

        Ok(AnnotateResult {
            output_path: output_path.display().to_string(),
            scan_number: annotation.scan_number,
            peptide_sequence: annotation.peptide_sequence.clone(),
            charge: annotation.charge,
            score: annotation.score,
            matched_ions: annotation.matched_ions,
            total_ions: annotation.total_ions,
            delta_mass_ppm: annotation.delta_mass_ppm,
            protein_accessions: annotation.protein_accessions.clone(),
            message: format!(
                "Annotation saved to {}. Matched {}/{} ions (score: {:.3}). {}",
                output_path.display(),
                annotation.matched_ions,
                annotation.total_ions,
                annotation.score,
                match render_mode {
                    "unified+xic" => "Includes XIC + SILAC controls.",
                    "unified" => "Unified view (annotation only).",
                    _ => "Open in browser to view.",
                },
            ),
        })
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
    if !std::path::Path::new(path).exists() {
        return Err(mcp_err(
            ErrorCode::INVALID_PARAMS,
            format!("file does not exist: {path}"),
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
        description = "Read a mass spectrometry file (supports mzML, MGF, and PFB) and return a statistical summary including spectrum count, m/z range, RT range, charge distribution, and median peaks per spectrum. Use this as the first step to understand input data."
    )]
    fn read_spectra(
        &self,
        Parameters(input): Parameters<ReadSpectraInput>,
    ) -> Result<Json<SpectrumSummary>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "read_spectra").entered();
        tracing::info!(file = %input.file_path, "started");
        validate_file_path(&input.file_path)?;
        let path = Path::new(&input.file_path);
        let reader = self.get_or_create_reader(path)?;
        let summary = reader
            .read_summary(path)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        tracing::info!(
            ms1 = summary.ms1_count,
            ms2 = summary.ms2_count,
            total = summary.total_spectra,
            "completed"
        );
        Ok(Json(summary))
    }

    /// Read a single spectrum by scan number.
    #[rmcp::tool(
        name = "get_spectrum",
        description = "Read a single spectrum from a file by scan number (1-based). Supports mzML, MGF, and PFB formats. Returns the spectrum with m/z array, intensity array, precursor info, and MS level."
    )]
    fn get_spectrum(
        &self,
        Parameters(input): Parameters<GetSpectrumInput>,
    ) -> Result<Json<Spectrum>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "get_spectrum").entered();
        tracing::info!(file = %input.file_path, scan = input.scan_number, "started");
        validate_file_path(&input.file_path)?;
        validate_scan_number(input.scan_number)?;
        let path = Path::new(&input.file_path);
        let reader = self.get_or_create_reader(path)?;
        let spectrum = reader
            .read_spectrum(path, input.scan_number)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        tracing::info!("completed");
        Ok(Json(spectrum))
    }

    /// Recommend search parameters based on spectrum characteristics.
    #[rmcp::tool(
        name = "recommend_params",
        description = "Recommend search parameters based on spectrum file characteristics. Input: SpectrumSummary from read_spectra + optional UserHints (experiment_type, instrument_type, enzyme). Output: recommended SearchParams with confidence score and explanation. Note: set database_path in params to the FASTA file path. For a one-shot convenience tool that combines reading, recommending, and database resolution, use prepare_search."
    )]
    fn recommend_params(
        &self,
        Parameters(input): Parameters<RecommendParamsInput>,
    ) -> Result<Json<AiDecision<SearchParams>>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "recommend_params").entered();
        tracing::info!(
            file = input.file_path.as_deref().unwrap_or("none"),
            database = input.database_path.as_deref().unwrap_or("none"),
            has_summary = input.summary.is_some(),
            has_hints = input.hints.is_some(),
            "started"
        );
        // Get summary: use provided summary or read from file_path
        let summary = if let Some(s) = input.summary {
            s
        } else if let Some(ref fp) = input.file_path {
            validate_file_path(fp)?;
            let path = std::path::Path::new(fp);
            let reader = self.get_or_create_reader(path)?;
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

        tracing::info!(enzyme = ?decision.decision.enzyme, confidence = decision.confidence, "completed");
        Ok(Json(decision))
    }

    /// List all available search parameter presets.
    #[rmcp::tool(
        name = "list_presets",
        description = "List all built-in search parameter presets (standard, phospho, TMT, SILAC, open search). Each preset includes name, description, parameters, and applicable scenarios."
    )]
    fn list_presets(&self) -> Json<PresetsResponse> {
        let _span = tracing::info_span!("mcp_tool", name = "list_presets").entered();
        tracing::info!("started");
        let presets = ParamRecommender::list_presets();
        tracing::info!(count = presets.len(), "completed");
        Json(PresetsResponse { presets })
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
        let _span = tracing::info_span!("mcp_tool", name = "run_search");
        let _enter = _span.enter();
        tracing::info!(
            engine = input.params.as_ref().and_then(|p| p.engine.as_deref()).unwrap_or("auto"),
            files = input.input_files.len(),
            dia_run_id = input.dia_run_id.as_deref().unwrap_or("none"),
            database = input.params.as_ref().map(|p| p.database_path.as_str()).or(input.database_path.as_deref()).unwrap_or("auto"),
            enzyme = ?input.params.as_ref().map(|p| &p.enzyme),
            precursor_tol = ?input.params.as_ref().map(|p| &p.precursor_tolerance),
            fragment_tol = ?input.params.as_ref().map(|p| &p.fragment_tolerance),
            missed_cleavages = input.params.as_ref().map(|p| p.missed_cleavages).unwrap_or(0),
            fixed_mods = input.params.as_ref().map(|p| p.fixed_modifications.len()).unwrap_or(0),
            var_mods = input.params.as_ref().map(|p| p.variable_modifications.len()).unwrap_or(0),
            "started"
        );
        drop(_enter);
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

            let db_path = Path::new(&params.database_path);
            if !db_path.exists() {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    format!(
                        "Database file not found: {}. Use list_databases to see available databases, \
                         or download_database to fetch one.",
                        params.database_path
                    ),
                ));
            }

            // Resolve + construct the engine BEFORE any cache mutation, so an
            // invalid engine name fails fast without consuming cached spectra
            // or inserting an orphaned run entry.
            let engine_name = params.engine.as_deref().unwrap_or("SimpleSearch");
            let engine = self.resolve_engine(engine_name)?;
            let engine_info_for_error = engine.engine_info();

            // Move spectra out of the DIA cache (frees the slot).
            // This happens after engine + param validation so a failure won't
            // lose cached spectra.
            let dia_spectra = {
                let mut cache = self.dia_cache.lock().unwrap_or_else(|e| e.into_inner());
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
                let mut cache = self.run_cache.lock().unwrap_or_else(|e| e.into_inner());
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
                            error_category: None,
                            has_diagnostics: false,
                        },
                        result: None,
                        handle: None,
                        diagnostics: None,
                        params_used: None,
                    },
                );
            }

            let run_cache_clone = Arc::clone(&self.run_cache);
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
                // Keep: feeds PanicGuard + SearchProgress.elapsed_sec / MCP client
                let start = Instant::now();

                // Panic guard: on abnormal exit, sets status to "Failed: task panicked"
                let _guard = PanicGuard {
                    run_id,
                    cache: Arc::clone(&run_cache_clone),
                    start,
                };

                let mut diagnostics = protein_copilot_core::diagnostics::SearchDiagnostics::new();
                let search_result = engine
                    .search_with_spectra(&params, dia_spectra, on_progress, &mut diagnostics)
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
                                    // Populate input_files for DIA results
                                    result.metadata.input_files = dia_source.clone();
                                    // Run anomaly detection
                                    let tol_ppm = if params.precursor_tolerance.unit
                                        == protein_copilot_core::search_params::ToleranceUnit::Ppm
                                    {
                                        Some(params.precursor_tolerance.value)
                                    } else {
                                        None
                                    };
                                    let decoy_count =
                                        result.psms.iter().filter(|p| p.is_decoy).count() as u64;
                                    diagnostics.finalize(
                                        Some(result.summary.identification_rate),
                                        Some(result.summary.psms_at_1pct_fdr),
                                        Some(decoy_count),
                                        duration,
                                        tol_ppm,
                                    );
                                    state.diagnostics = Some(diagnostics);
                                    state.params_used = Some(params.clone());
                                    state.progress.has_diagnostics = true;
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
                                    state.progress.error_category =
                                        diagnostics.error_category.clone();
                                    state.diagnostics = Some(diagnostics);
                                    state.params_used = Some(params.clone());
                                    state.progress.has_diagnostics = true;
                                    let entry = crate::history::SearchHistoryEntry {
                                        run_id,
                                        status: format!("Failed: {e}"),
                                        created_at: chrono::Utc::now(),
                                        elapsed_sec: duration,
                                        engine_info: engine_info_for_error.clone(),
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

            tracing::info!(run_id = %run_id, "completed");
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

        // Validate that all input files exist before proceeding
        for file_str in &input.input_files {
            let p = Path::new(file_str);
            if !p.exists() {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    format!("input file does not exist: {file_str}"),
                ));
            }
        }

        let mut params = if let Some(p) = input.params {
            p
        } else {
            let first_file = input
                .input_files
                .first()
                .ok_or_else(|| mcp_err(ErrorCode::INVALID_PARAMS, "input_files is empty"))?;
            let path = Path::new(first_file);
            let reader = self.get_or_create_reader(path)?;
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

        let db_path = Path::new(&params.database_path);
        if !db_path.exists() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!(
                    "Database file not found: {}. Use list_databases to see available databases, \
                     or download_database to fetch one.",
                    params.database_path
                ),
            ));
        }

        // DIA guard: detect if input is DIA and reject without dia_run_id.
        // Raw DIA spectra have isolation window centers as precursor m/z, not real
        // precursors. Searching without extraction produces false positives.
        {
            let first_path = match input.input_files.first() {
                Some(p) => Path::new(p),
                None => {
                    return Err(mcp_err(
                        ErrorCode::INVALID_PARAMS,
                        "input_files is empty (internal error)",
                    ));
                }
            };
            if let Ok(reader) = self.get_or_create_reader(first_path) {
                if let Ok(summary) = reader.read_summary(first_path) {
                    if let Some(w) = summary.median_isolation_window_da {
                        if w > DIA_ISOLATION_WINDOW_THRESHOLD_DA {
                            return Err(mcp_err(
                                ErrorCode::INVALID_PARAMS,
                                format!(
                                    "Input file appears to be DIA data (median isolation window \
                                     {w:.1} Da > 5.0 Da). DIA requires precursor extraction \
                                     before searching. Please call extract_dia_precursors first, \
                                     then use the returned dia_run_id with run_search."
                                ),
                            ));
                        }
                    }
                }
            }
        }

        // Resolve + construct the engine BEFORE any cache mutation, so an
        // invalid engine name fails fast without inserting an orphaned run entry.
        let engine_name = params.engine.as_deref().unwrap_or("SimpleSearch");
        let engine = self.resolve_engine(engine_name)?;
        let engine_info_for_error = engine.engine_info();

        let run_id = Uuid::new_v4();
        let files: Vec<PathBuf> = input.input_files.iter().map(PathBuf::from).collect();

        // Evict + initialize in unified cache
        {
            let mut cache = self.run_cache.lock().unwrap_or_else(|e| e.into_inner());
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
                        error_category: None,
                        has_diagnostics: false,
                    },
                    result: None,
                    handle: None,
                    diagnostics: None,
                    params_used: None,
                },
            );
        }

        let run_cache_clone = Arc::clone(&self.run_cache);

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
            // Keep: feeds PanicGuard + SearchProgress.elapsed_sec / MCP client
            let start = Instant::now();

            // Panic guard: on abnormal exit, sets status to "Failed: task panicked"
            let _guard = PanicGuard {
                run_id,
                cache: Arc::clone(&run_cache_clone),
                start,
            };

            let mut diagnostics = protein_copilot_core::diagnostics::SearchDiagnostics::new();
            let search_result = engine
                .search(&params, &files, on_progress, &mut diagnostics)
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
                                // Run anomaly detection
                                let tol_ppm = if params.precursor_tolerance.unit
                                    == protein_copilot_core::search_params::ToleranceUnit::Ppm
                                {
                                    Some(params.precursor_tolerance.value)
                                } else {
                                    None
                                };
                                let decoy_count =
                                    result.psms.iter().filter(|p| p.is_decoy).count() as u64;
                                diagnostics.finalize(
                                    Some(result.summary.identification_rate),
                                    Some(result.summary.psms_at_1pct_fdr),
                                    Some(decoy_count),
                                    duration,
                                    tol_ppm,
                                );
                                state.diagnostics = Some(diagnostics);
                                state.params_used = Some(params.clone());
                                state.progress.has_diagnostics = true;
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
                                state.progress.error_category = diagnostics.error_category.clone();
                                state.diagnostics = Some(diagnostics);
                                state.params_used = Some(params.clone());
                                state.progress.has_diagnostics = true;
                                let entry = crate::history::SearchHistoryEntry {
                                    run_id,
                                    status: format!("Failed: {e}"),
                                    created_at: chrono::Utc::now(),
                                    elapsed_sec: duration,
                                    engine_info: engine_info_for_error.clone(),
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

        tracing::info!(run_id = %run_id, "completed");
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
        description = "Check the status of a search started by run_search. Returns progress percentage and elapsed time. When status is Completed, use generate_summary(run_id) to get results. If the search failed or has warnings, call diagnose_search(run_id) for a diagnostic report."
    )]
    fn get_search_status(
        &self,
        Parameters(input): Parameters<GetSearchStatusInput>,
    ) -> Result<Json<SearchProgress>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "get_search_status").entered();
        tracing::info!(run_id = %input.run_id, "started");
        let id = Uuid::parse_str(&input.run_id)
            .map_err(|_| mcp_err(ErrorCode::INVALID_PARAMS, "invalid run_id format"))?;
        let cache = self.run_cache.lock().unwrap_or_else(|e| e.into_inner());
        let state = cache
            .get(&id)
            .ok_or_else(|| mcp_err(ErrorCode::INVALID_PARAMS, format!("run_id {id} not found — it may have been evicted from the cache (max 100 recent runs are kept)")))?;
        tracing::info!(status = %state.progress.status, "completed");
        Ok(Json(state.progress.clone()))
    }

    /// Cancel a running search.
    #[rmcp::tool(
        name = "cancel_search",
        description = "Cancel a running search started by run_search. The search task is immediately terminated and status is set to Cancelled. Partial results may still be available via generate_summary."
    )]
    fn cancel_search(
        &self,
        Parameters(input): Parameters<CancelSearchInput>,
    ) -> Result<Json<SearchProgress>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "cancel_search").entered();
        tracing::info!(run_id = %input.run_id, "started");
        let id = Uuid::parse_str(&input.run_id)
            .map_err(|_| mcp_err(ErrorCode::INVALID_PARAMS, "invalid run_id format"))?;
        let mut cache = self.run_cache.lock().unwrap_or_else(|e| e.into_inner());
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

        tracing::info!("completed");
        Ok(Json(state.progress.clone()))
    }

    /// Check available search engines and their health status.
    #[rmcp::tool(
        name = "check_engine",
        description = "Check available search engines and their health status. Returns engine name, version, supported features, and availability."
    )]
    async fn check_engine(&self) -> Json<EngineStatus> {
        let _span = tracing::info_span!("mcp_tool", name = "check_engine");
        let _enter = _span.enter();
        tracing::info!("started");
        drop(_enter);
        let engines = self.registry.list_available();
        let all_engines = engines.clone();
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
            tracing::info!("completed");
            Json(EngineStatus {
                engine: info.clone(),
                status,
                all_engines,
            })
        } else {
            tracing::info!("completed");
            Json(EngineStatus {
                engine: protein_copilot_core::engine::EngineInfo {
                    name: "none".to_string(),
                    version: "0.0.0".to_string(),
                    supported_features: vec![],
                },
                status: HealthStatus::Unavailable {
                    reason: "no engines registered".to_string(),
                },
                all_engines: vec![],
            })
        }
    }

    /// Generate a statistical summary with FDR filtering.
    #[rmcp::tool(
        name = "generate_summary",
        description = "Generate a statistical summary from search results with 1% FDR filtering (provide either result or run_id, not both). Includes identification rate, median score, median delta ppm, modification and charge distributions. Use this after run_search to interpret results. To export raw data as TSV/JSON files, use export_results."
    )]
    fn generate_summary(
        &self,
        Parameters(input): Parameters<GenerateSummaryInput>,
    ) -> Result<Json<SearchResultSummary>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "generate_summary").entered();
        tracing::info!(
            run_id = input.run_id.as_deref().unwrap_or("direct"),
            "started"
        );
        let result = self.get_result(&input.result, &input.run_id)?;
        let summary = ReportGenerator::generate_summary(&result);
        tracing::info!(
            psms_1pct = summary.psms_at_1pct_fdr,
            id_rate = summary.identification_rate,
            "completed"
        );
        Ok(Json(summary))
    }

    /// Export search results as TSV and JSON files.
    #[rmcp::tool(
        name = "export_results",
        description = "Export search results from run_search to files (provide run_id). Creates psm.tsv, peptide.tsv, protein.tsv, result.json, and run_metadata.json in the specified output directory. For a statistical overview, call generate_summary first."
    )]
    fn export_results(
        &self,
        Parameters(input): Parameters<ExportResultsInput>,
    ) -> Result<Json<ExportResultsOutput>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "export_results").entered();
        tracing::info!(
            output_dir = %input.output_dir,
            run_id = input.run_id.as_deref().unwrap_or("direct"),
            "started"
        );
        let result = self.get_result(&input.result, &input.run_id)?;
        let output_dir = Path::new(&input.output_dir);

        ReportGenerator::export_tsv(&result, output_dir)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        ReportGenerator::export_json(&result, &output_dir.join("result.json"))
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        ReportGenerator::export_metadata(&result.metadata, &output_dir.join("run_metadata.json"))
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

        let files = vec![
            "psm.tsv".to_string(),
            "peptide.tsv".to_string(),
            "protein.tsv".to_string(),
            "result.json".to_string(),
            "run_metadata.json".to_string(),
        ];

        tracing::info!("completed");
        Ok(Json(ExportResultsOutput {
            output_dir: output_dir.display().to_string(),
            files,
        }))
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
        let _span = tracing::info_span!("mcp_tool", name = "list_searches").entered();
        tracing::info!(
            status_filter = input.status_filter.as_deref().unwrap_or("all"),
            limit = input.limit.unwrap_or(20),
            "started"
        );
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
                            engine: None,
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
        entries.sort_by_key(|e| std::cmp::Reverse(e.created_at));
        entries.truncate(limit);
        tracing::info!(count = entries.len(), "completed");
        Json(ListSearchesResponse { searches: entries })
    }

    /// Annotate a single spectrum with peptide fragment ion matching.
    #[rmcp::tool(
        name = "annotate_spectrum",
        description = "Annotate a single spectrum with peptide fragment ion matching. Generates an interactive HTML file showing matched b/y ions. Two modes (choose one): Mode 1 (from search) — provide run_id + scan_number to annotate an existing PSM. Mode 2 (manual) — provide file_path + scan_number + peptide_sequence + charge. In mode 2 only: set scan_number=0 with retention_time_min to auto-find the nearest scan. Supports mzML, MGF, and PFB formats. mzML and PFB enable XIC extraction (set n_cycles, top_n_ions, extraction_tolerance) and SILAC heavy/light comparison (set label_type). MGF provides basic fragment annotation only. For batch annotation, use batch_annotate_spectra."
    )]
    fn annotate_spectrum(
        &self,
        Parameters(input): Parameters<AnnotateSpectrumInput>,
    ) -> Result<Json<AnnotateResult>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "annotate_spectrum").entered();
        tracing::info!(
            scan = input.scan_number,
            peptide = input.peptide_sequence.as_deref().unwrap_or("from_psm"),
            charge = input.charge.unwrap_or(0),
            run_id = input.run_id.as_deref().unwrap_or("none"),
            file = input.file_path.as_deref().unwrap_or("none"),
            "started"
        );
        use protein_copilot_core::search_params::{MassTolerance, Modification, ToleranceUnit};

        // Allow scan_number=0 only when retention_time_min is provided for auto-lookup
        if input.scan_number == 0 && input.retention_time_min.is_none() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                "scan_number must be >= 1, or set scan_number=0 with retention_time_min for auto lookup",
            ));
        }
        if input.scan_number != 0 {
            validate_scan_number(input.scan_number)?;
        }
        // RT-based lookup is only for Mode 2 (manual); Mode 1 (run_id) requires valid scan_number
        if input.scan_number == 0 && input.run_id.is_some() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                "retention_time_min auto-lookup is only supported in manual mode (without run_id). With run_id, provide a valid scan_number.",
            ));
        }

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
            if pep.trim().is_empty() {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "peptide_sequence cannot be empty in manual annotation mode",
                ));
            }
            if ch <= 0 {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    format!("charge must be > 0, got {ch}"),
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

        // Defense-in-depth: validate charge > 0 before any m/z calculation.
        // Mode 2 validates at entry; this catches mode 1 PSMs with invalid charge.
        if charge <= 0 {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!("charge must be > 0 (got {charge}); PSM may have invalid charge state"),
            ));
        }

        // Resolve scan_number: if 0 and retention_time_min provided, auto-match via RT
        let resolved_scan = if input.scan_number == 0 {
            if let Some(rt) = input.retention_time_min {
                let reader = self.get_or_create_reader(&spectrum_file)?;
                use protein_copilot_search_engine::chemistry::{
                    residue_mass, PROTON_MASS, WATER_MASS,
                };
                let base_mass: f64 =
                    peptide_seq.chars().filter_map(residue_mass).sum::<f64>() + WATER_MASS;
                let mod_mass: f64 = modifications.iter().map(|m| m.mass_delta).sum();
                let precursor_mz =
                    (base_mass + mod_mass + charge as f64 * PROTON_MASS) / charge as f64;
                reader
                    .find_by_rt(
                        &spectrum_file,
                        rt,
                        precursor_mz,
                        RT_AUTO_LOOKUP_TOLERANCE_MIN,
                    )
                    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?
                    .map(|(scan, _)| scan)
                    .ok_or_else(|| {
                        mcp_err(
                            ErrorCode::INVALID_PARAMS,
                            format!("No MS2 scan near RT={rt:.2}min mz={precursor_mz:.4}"),
                        )
                    })?
            } else {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "scan_number is 0: provide retention_time_min for auto lookup",
                ));
            }
        } else {
            input.scan_number
        };

        // Read the spectrum
        let reader = self.get_or_create_reader(&spectrum_file)?;
        let spectrum = reader
            .read_spectrum(&spectrum_file, resolved_scan)
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
            false,
            false,
        )
        .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))?;

        // Set source file name on annotation for display
        let mut annotation = annotation;
        annotation.source_file = spectrum_file
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_default();

        // ── SILAC: find and annotate heavy scan (DIA or DDA) ──
        if let Some(ref label) = input.label_type {
            let core_label: &protein_copilot_core::label::LabelType = label;
            let heavy_delta =
                protein_copilot_core::label::total_heavy_delta(&peptide_seq, core_label);

            if heavy_delta.abs() < 1e-6 {
                tracing::info!(
                    peptide = peptide_seq,
                    "Skipping heavy annotation: peptide has no K/R, zero SILAC shift"
                );
            } else {
                let is_dia = spectrum
                    .precursors
                    .first()
                    .and_then(|p| p.isolation_window.as_ref())
                    .map(|w| (w.lower_offset + w.upper_offset) > 1.0)
                    .unwrap_or(false);

                let heavy_prec_mz = protein_copilot_core::label::compute_heavy_precursor_mz(
                    annotation.theoretical_mz,
                    charge,
                    &peptide_seq,
                    core_label,
                );

                // Use O(log N) binary search to find the heavy scan via RT + m/z.
                // For DIA: find_by_rt checks isolation window contains heavy_prec_mz.
                // For DDA: find_by_rt accepts scans without isolation windows (RT-only
                // fallback), then we verify the precursor m/z at 20 ppm below.
                let heavy_scan_result = if is_dia {
                    // DIA: binary search finds scan whose isolation window contains heavy m/z
                    reader
                        .find_by_rt(
                            &spectrum_file,
                            spectrum.retention_time_min,
                            heavy_prec_mz,
                            0.5,
                        )
                        .unwrap_or(None)
                } else {
                    // DDA: binary search finds RT-nearest MS2 scan, then verify precursor m/z
                    // We use a wider RT tolerance since DDA heavy scans can be several scans away
                    let candidates = reader
                        .find_by_rt(
                            &spectrum_file,
                            spectrum.retention_time_min,
                            heavy_prec_mz,
                            0.5,
                        )
                        .unwrap_or(None);
                    // Verify precursor m/z match at 20 ppm for DDA
                    if let Some((scan, delta)) = candidates {
                        match reader.read_spectrum(&spectrum_file, scan) {
                            Ok(spec) => {
                                let prec_mz = spec.precursors.first().map(|p| p.mz).unwrap_or(0.0);
                                let ppm_err =
                                    ((prec_mz - heavy_prec_mz) / heavy_prec_mz * 1e6).abs();
                                if ppm_err < 20.0 {
                                    Some((scan, delta))
                                } else {
                                    None
                                }
                            }
                            Err(_) => None,
                        }
                    } else {
                        None
                    }
                };

                let heavy_scan_num = heavy_scan_result.map(|(scan, _)| scan);

                if let Some(heavy_scan_num) = heavy_scan_num {
                    if let Ok(heavy_spectrum) = reader.read_spectrum(&spectrum_file, heavy_scan_num)
                    {
                        match protein_copilot_search_engine::annotate::annotate_heavy_spectrum(
                            &heavy_spectrum,
                            &peptide_seq,
                            charge,
                            &frag_tol,
                            &modifications,
                            core_label,
                            false,
                            false,
                        ) {
                            Ok(heavy_ann) => {
                                tracing::info!(
                                    heavy_scan = heavy_scan_num,
                                    heavy_prec_mz = format!("{:.4}", heavy_prec_mz),
                                    matched = heavy_ann.matched_ions,
                                    total = heavy_ann.total_ions,
                                    is_dia = is_dia,
                                    "Heavy annotation complete"
                                );
                                annotation.heavy_annotation = Some(heavy_ann);
                            }
                            Err(e) => {
                                tracing::warn!("Heavy annotation failed: {e}");
                            }
                        }
                    }
                } else {
                    let mode = if is_dia {
                        "DIA window"
                    } else {
                        "DDA precursor"
                    };
                    tracing::info!(
                        heavy_prec_mz = format!("{:.4}", heavy_prec_mz),
                        mode = mode,
                        "No {mode} found for heavy precursor"
                    );
                }
            }
        }

        let out_path = input.output_path.map(PathBuf::from).unwrap_or_else(|| {
            default_output_path(&format!("annotation_scan{}.html", resolved_scan))
        });
        let out_path = absolutize(&out_path);
        ensure_parent_dir(&out_path)?;

        let is_mzml = spectrum_file
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.eq_ignore_ascii_case("mzml") || e.eq_ignore_ascii_case("pfb"))
            .unwrap_or(false);

        let source_file = spectrum_file
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| spectrum_file.display().to_string());

        let make_peptide_info =
            |seq: &str, ch: i32, mz: f64| protein_copilot_report::unified_types::PeptideInfo {
                sequence: seq.to_string(),
                charge: ch,
                precursor_mz: mz,
                total_k: seq.chars().filter(|&c| c == 'K').count() as u32,
                total_r: seq.chars().filter(|&c| c == 'R').count() as u32,
            };
        let plotly_mode = input
            .plotly_mode
            .unwrap_or(protein_copilot_xic::PlotlyMode::Cdn);
        let annotation_theo_mz = annotation.theoretical_mz;
        let render_unified_without_xic = || -> Result<(), ErrorData> {
            let unified_data = protein_copilot_report::unified_types::UnifiedViewData {
                source_file: source_file.clone(),
                annotation: annotation.clone(),
                xic: None,
                raw_scans: None,
                ion_metadata: vec![],
                peptide_info: make_peptide_info(&peptide_seq, charge, annotation.theoretical_mz),
            };
            ReportGenerator::render_unified(&unified_data, &out_path, plotly_mode)
                .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
            Ok(())
        };

        let mut render_mode = "annotation";
        if is_mzml {
            let xic_params = protein_copilot_xic::ExtractionParams {
                mz_tolerance: input.extraction_tolerance.unwrap_or(MassTolerance {
                    value: 20.0,
                    unit: ToleranceUnit::Ppm,
                }),
                n_cycles: input.n_cycles.unwrap_or(5),
                top_n_ions: input.top_n_ions.unwrap_or(usize::MAX),
                label_type: input.label_type.clone(),
                intensity_rule: protein_copilot_xic::IntensityRule::MaxInWindow,
            };

            // Get cached indexed reader for O(1) scan lookups
            let cached_reader = self.get_or_create_reader(&spectrum_file)?;

            match protein_copilot_xic::extract::extract_xic_unified(
                cached_reader.as_ref(),
                &spectrum_file,
                resolved_scan,
                &peptide_seq,
                charge,
                annotation.theoretical_mz,
                &modifications,
                &xic_params,
                20.0,
            ) {
                Ok(unified_result) => {
                    let xic_meaningful = unified_result
                        .xic_data
                        .fragment_xic_traces
                        .first()
                        .map(|t| t.data_points.len() > 1)
                        .unwrap_or(false);

                    if xic_meaningful {
                        // Refine precursor_mz from MS1 scan
                        let mut annotation = annotation.clone();
                        if let Some(observed) = find_precursor_in_ms1(
                            &unified_result.raw_scans.ms1_scans,
                            annotation.retention_time_min,
                            annotation.theoretical_mz,
                            20.0,
                        ) {
                            annotation.precursor_mz = observed;
                            annotation.delta_mass_ppm = (observed - annotation.theoretical_mz)
                                / annotation.theoretical_mz
                                * 1e6;
                        }

                        // Also refine heavy precursor delta from MS1
                        if let Some(ref mut ha) = annotation.heavy_annotation {
                            let heavy_theo = ha.precursor_mz;
                            if let Some(obs_heavy) = find_precursor_in_ms1(
                                &unified_result.raw_scans.ms1_scans,
                                annotation.retention_time_min,
                                heavy_theo,
                                20.0,
                            ) {
                                ha.delta_mass_ppm =
                                    Some((obs_heavy - heavy_theo) / heavy_theo * 1e6);
                            } else {
                                ha.delta_mass_ppm = None;
                            }
                        }

                        let unified_data = protein_copilot_report::unified_types::UnifiedViewData {
                            source_file: source_file.clone(),
                            annotation,
                            xic: Some(unified_result.xic_data),
                            raw_scans: Some(unified_result.raw_scans),
                            ion_metadata: unified_result.ion_metadata,
                            peptide_info: make_peptide_info(
                                &peptide_seq,
                                charge,
                                annotation_theo_mz,
                            ),
                        };

                        ReportGenerator::render_unified(&unified_data, &out_path, plotly_mode)
                            .map_err(|e| {
                                mcp_core_err(protein_copilot_core::error::CoreError::from(e))
                            })?;
                        render_mode = "unified+xic";
                    } else {
                        render_unified_without_xic()?;
                        render_mode = "unified";
                    }
                }
                Err(_) => {
                    render_unified_without_xic()?;
                    render_mode = "unified";
                }
            }
        } else {
            // Non-mzML/PFB file (e.g. MGF) — basic annotation only, no XIC
            ReportGenerator::render_annotation(&annotation, &out_path)
                .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        }

        tracing::info!("completed");
        Ok(Json(AnnotateResult {
            output_path: out_path.display().to_string(),
            scan_number: annotation.scan_number,
            peptide_sequence: annotation.peptide_sequence.clone(),
            charge: annotation.charge,
            score: annotation.score,
            matched_ions: annotation.matched_ions,
            total_ions: annotation.total_ions,
            delta_mass_ppm: annotation.delta_mass_ppm,
            protein_accessions: annotation.protein_accessions.clone(),
            message: format!(
                "Annotation saved to {}. Matched {}/{} ions (score: {:.3}). {}",
                out_path.display(),
                annotation.matched_ions,
                annotation.total_ions,
                annotation.score,
                match render_mode {
                    "unified+xic" => "Includes XIC + SILAC controls (DIA).",
                    "unified" => "Unified view (annotation only, no XIC for DDA).",
                    _ => "Open in browser to view.",
                },
            ),
        }))
    }

    /// Batch-annotate multiple spectra with peptide fragment ion matching.
    ///
    /// Two input modes:
    /// (1) pass `items` as a JSON array of {scan_number, peptide_sequence, charge, modifications?} objects,
    /// (2) pass `annotations_file` pointing to a TSV file with the same columns.
    ///
    /// XIC and SILAC parameters are global — all items share the same settings.
    /// Generates one HTML file per spectrum in `output_dir` (default: ./batch_annotations/).
    /// Failed items are recorded in the result with error messages; the batch continues.
    #[rmcp::tool(
        name = "batch_annotate_spectra",
        description = "Annotate multiple spectra with peptide fragment ion matching — batch version of annotate_spectrum. Two input modes: (1) provide `items` as a JSON array of {scan_number, peptide_sequence, charge, modifications?}, or (2) provide `annotations_file` pointing to a TSV file with columns: scan_number, peptide_sequence, charge, modifications (optional). Supports mzML, MGF, PFB formats. Set label_type for SILAC, n_cycles for DIA XIC extraction. Output: one HTML per scan in output_dir."
    )]
    fn batch_annotate_spectra(
        &self,
        Parameters(input): Parameters<BatchAnnotateInput>,
    ) -> Result<Json<BatchAnnotateOutput>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "batch_annotate_spectra").entered();
        tracing::info!(
            file = %input.file_path,
            items = input.items.as_ref().map(|v| v.len()).unwrap_or(0),
            annotations_file = input.annotations_file.as_deref().unwrap_or("none"),
            "started"
        );

        validate_file_path(&input.file_path)?;
        let spectrum_file = PathBuf::from(&input.file_path);

        // Resolve items from JSON array or TSV file
        let items: Vec<BatchAnnotateItem> = match (&input.items, &input.annotations_file) {
            (Some(json_items), None) => {
                if json_items.is_empty() {
                    return Err(mcp_err(ErrorCode::INVALID_PARAMS, "items array is empty"));
                }
                json_items.to_vec()
            }
            (None, Some(file)) => parse_annotations_file(file)?,
            (Some(_), Some(_)) => {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "provide either 'items' or 'annotations_file', not both",
                ));
            }
            (None, None) => {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "provide either 'items' (JSON array) or 'annotations_file' (TSV path)",
                ));
            }
        };
        let total = items.len();

        // Defaults
        use protein_copilot_core::search_params::{MassTolerance, Modification, ToleranceUnit};
        let frag_tol = input.fragment_tolerance.unwrap_or(MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        });
        let extraction_tolerance = input.extraction_tolerance.unwrap_or(MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        });
        let n_cycles = input.n_cycles.unwrap_or(5);
        let top_n_ions = input.top_n_ions.unwrap_or(usize::MAX);
        let plotly_mode = input
            .plotly_mode
            .unwrap_or(protein_copilot_xic::PlotlyMode::Cdn);
        let core_label: Option<&protein_copilot_core::label::LabelType> = input
            .label_type
            .as_ref()
            .map(|l| l as &protein_copilot_core::label::LabelType);

        let output_dir = input
            .output_dir
            .map(PathBuf::from)
            .unwrap_or_else(|| PathBuf::from("./batch_annotations"));
        std::fs::create_dir_all(&output_dir).map_err(|e| {
            mcp_err(
                ErrorCode::INTERNAL_ERROR,
                format!("cannot create output_dir '{}': {e}", output_dir.display()),
            )
        })?;
        let output_dir = output_dir.canonicalize().unwrap_or(output_dir);

        // Progress tracking
        let mut results: Vec<BatchAnnotateItemResult> = Vec::with_capacity(total);
        let mut succeeded = 0usize;
        let mut failed = 0usize;

        for (_idx, item) in items.iter().enumerate() {
            // Validate per-item
            if item.peptide_sequence.trim().is_empty() {
                results.push(BatchAnnotateItemResult {
                    scan_number: item.scan_number,
                    peptide_sequence: item.peptide_sequence.clone(),
                    charge: item.charge,
                    output_path: None,
                    matched_ions: None,
                    total_ions: None,
                    score: None,
                    delta_mass_ppm: None,
                    protein_accessions: item.protein_accessions.clone().unwrap_or_default(),
                    error: Some("peptide_sequence is empty".to_string()),
                });
                failed += 1;
                continue;
            }
            if item.charge <= 0 {
                results.push(BatchAnnotateItemResult {
                    scan_number: item.scan_number,
                    peptide_sequence: item.peptide_sequence.clone(),
                    charge: item.charge,
                    output_path: None,
                    matched_ions: None,
                    total_ions: None,
                    score: None,
                    delta_mass_ppm: None,
                    protein_accessions: item.protein_accessions.clone().unwrap_or_default(),
                    error: Some(format!("charge must be > 0, got {}", item.charge)),
                });
                failed += 1;
                continue;
            }
            if item.scan_number == 0 {
                results.push(BatchAnnotateItemResult {
                    scan_number: item.scan_number,
                    peptide_sequence: item.peptide_sequence.clone(),
                    charge: item.charge,
                    output_path: None,
                    matched_ions: None,
                    total_ions: None,
                    score: None,
                    delta_mass_ppm: None,
                    protein_accessions: item.protein_accessions.clone().unwrap_or_default(),
                    error: Some("scan_number must be >= 1".to_string()),
                });
                failed += 1;
                continue;
            }

            let modifications: Vec<Modification> = item
                .modifications
                .clone()
                .unwrap_or_default()
                .into_iter()
                .map(|(pos, delta)| Modification {
                    name: format!("Mod@{}", pos),
                    mass_delta: delta,
                    residues: vec![item.peptide_sequence.chars().nth(pos).unwrap_or('*')],
                    position: protein_copilot_core::search_params::ModPosition::Anywhere,
                })
                .collect();
            let protein_accs = item.protein_accessions.clone().unwrap_or_default();

            let out_path = output_dir.join(format!("annotation_scan{}.html", item.scan_number));

            let scan = item.scan_number;
            let pep = item.peptide_sequence.clone();
            let ch = item.charge;

            match self.annotate_one_spectrum(
                &spectrum_file,
                scan,
                &pep,
                ch,
                &modifications,
                protein_accs.clone(),
                &frag_tol,
                core_label,
                n_cycles,
                top_n_ions,
                &extraction_tolerance,
                plotly_mode,
                &out_path,
            ) {
                Ok(ann) => {
                    results.push(BatchAnnotateItemResult {
                        scan_number: scan,
                        peptide_sequence: pep,
                        charge: ch,
                        output_path: Some(ann.output_path),
                        matched_ions: Some(ann.matched_ions),
                        total_ions: Some(ann.total_ions),
                        score: Some(ann.score),
                        delta_mass_ppm: Some(ann.delta_mass_ppm),
                        protein_accessions: ann.protein_accessions,
                        error: None,
                    });
                    succeeded += 1;
                }
                Err(e) => {
                    results.push(BatchAnnotateItemResult {
                        scan_number: scan,
                        peptide_sequence: pep,
                        charge: ch,
                        output_path: None,
                        matched_ions: None,
                        total_ions: None,
                        score: None,
                        delta_mass_ppm: None,
                        protein_accessions: protein_accs,
                        error: Some(format!("{e}")),
                    });
                    failed += 1;
                }
            }
        }

        tracing::info!(total, succeeded, failed, "completed");
        Ok(Json(BatchAnnotateOutput {
            total,
            succeeded,
            failed,
            output_dir: output_dir.display().to_string(),
            results,
        }))
    }

    /// Extract candidate precursor ions from DIA mass spectrometry data.
    #[rmcp::tool(
        name = "extract_dia_precursors",
        description = "Extract candidate precursor ions from DIA mass spectrometry data. Reads spectrum file (mzML recommended; requires MS1+MS2 and isolation window metadata), detects DIA mode from isolation window widths, extracts precursor candidates from MS1 isotope patterns, and caches enhanced spectra for use with run_search. Returns a dia_run_id. Check cache status with get_dia_cache_status. For single-spectrum precursor extraction, use extract_spectrum_precursors."
    )]
    fn extract_dia_precursors(
        &self,
        Parameters(input): Parameters<ExtractDiaPrecursorsInput>,
    ) -> Result<Json<DiaExtractionOutput>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "extract_dia_precursors").entered();
        tracing::info!(
            file = %input.file_path,
            output_mode = %input.output_mode,
            min_charge = input.min_charge.unwrap_or(2),
            max_charge = input.max_charge.unwrap_or(5),
            acquisition_mode = input.acquisition_mode.as_deref().unwrap_or("auto"),
            "started"
        );
        validate_file_path(&input.file_path)?;
        let path = Path::new(&input.file_path);
        let reader = self.get_or_create_reader(path)?;
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
        let acq_mode = match input.acquisition_mode.as_deref() {
            Some(m) => {
                let upper = m.to_uppercase();
                match upper.as_str() {
                    "DDA" => Some(AcquisitionMode::DDA),
                    "DIA" => Some(AcquisitionMode::DIA),
                    _ => {
                        return Err(mcp_err(
                            ErrorCode::INVALID_PARAMS,
                            format!(
                                "invalid acquisition_mode '{}': expected 'DDA' or 'DIA' (case-insensitive). \
                                 Omit the field or set to 'auto' for automatic detection.",
                                m
                            ),
                        ));
                    }
                }
            }
            None => None,
        };
        let config = DiaExtractionConfig {
            acquisition_mode: acq_mode,
            ..DiaExtractionConfig::default()
        };

        let result = run_dia_extraction(&spectra, &extractor, &config)
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?;
        let detected_mode = result.detected_mode;
        let stats = result.stats.clone();

        let output_spectra = match input.output_mode.as_str() {
            "multi" => result.into_enhanced_spectra(),
            "pseudo" => result.into_pseudo_spectra(),
            other => {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    format!(
                        "invalid output_mode '{}': expected 'multi' or 'pseudo'",
                        other
                    ),
                ));
            }
        };

        let output_count = output_spectra.len() as u32;
        let run_id = Uuid::new_v4();

        // Cache for future use
        let mut cache = self.dia_cache.lock().unwrap_or_else(|e| e.into_inner());
        cache.insert(run_id, output_spectra);

        let output = DiaExtractionOutput {
            detected_mode: format!("{}", detected_mode),
            ms1_count: stats.ms1_count,
            ms2_count: stats.ms2_count,
            total_precursors_extracted: stats.total_precursors_extracted,
            avg_precursors_per_ms2: stats.avg_precursors_per_ms2,
            charge_distribution: stats.charge_distribution,
            output_spectra_count: output_count,
            run_id: run_id.to_string(),
            message: format!(
                "DIA extraction complete. {} precursors extracted from {} MS2 spectra. \
                 Pass dia_run_id=\"{}\" to run_search to search these spectra.",
                stats.total_precursors_extracted, stats.ms2_count, run_id
            ),
        };

        tracing::info!(candidates = output.total_precursors_extracted, "completed");
        Ok(Json(output))
    }

    /// Extract precursor candidates for a single MS2 spectrum from an mzML file.
    #[rmcp::tool(
        name = "extract_spectrum_precursors",
        description = "Extract candidate precursor ions from a single MS2 spectrum — single-spectrum version of extract_dia_precursors. Reads the spectrum file (mzML recommended; requires MS1+MS2 scans), finds the target MS2 by scan number, correlates it to the nearest MS1, and runs isotope pattern analysis within the isolation window. Returns extracted precursor candidates with charge states and the correlation method used."
    )]
    fn extract_spectrum_precursors(
        &self,
        Parameters(input): Parameters<ExtractSpectrumPrecursorsInput>,
    ) -> Result<Json<SingleSpectrumExtractionResult>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "extract_spectrum_precursors").entered();
        tracing::info!(file = %input.file_path, scan = input.scan_number, "started");
        validate_file_path(&input.file_path)?;
        validate_scan_number(input.scan_number)?;
        let path = Path::new(&input.file_path);
        let reader = self.get_or_create_reader(path)?;

        // Read target MS2 via O(1) indexed seek
        let target_ms2 = reader
            .read_spectrum(path, input.scan_number)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        if target_ms2.ms_level != protein_copilot_core::spectrum::MsLevel::MS2 {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!(
                    "scan {} is not MS2 (ms_level={:?})",
                    input.scan_number, target_ms2.ms_level
                ),
            ));
        }
        let target_rt_min = target_ms2.retention_time_min;

        // Use index to find nearby MS1 scans (±1 minute RT window)
        let scan_metas = reader
            .list_scan_meta(path)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        const MS1_RT_WINDOW_MIN: f64 = 1.0;
        let nearby_ms1_scans: Vec<u32> = scan_metas
            .iter()
            .filter(|m| m.ms_level == 1 && (m.rt_min - target_rt_min).abs() <= MS1_RT_WINDOW_MIN)
            .map(|m| m.scan_number)
            .collect();

        // Read only the needed spectra
        let mut spectra = Vec::with_capacity(nearby_ms1_scans.len() + 1);
        spectra.push(target_ms2);
        for scan_no in &nearby_ms1_scans {
            if let Ok(s) = reader.read_spectrum(path, *scan_no) {
                spectra.push(s);
            }
        }

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

        tracing::info!(candidates = result.precursors.len(), "completed");
        Ok(Json(result))
    }

    /// Check if a DIA extraction result is still cached.
    #[rmcp::tool(
        name = "get_dia_cache_status",
        description = "Check if a DIA extraction result is still cached for use with run_search. Set dia_run_id='' (empty string) to list all cached DIA extractions. Call this before run_search(dia_run_id=...) to verify the extraction hasn't been evicted."
    )]
    fn get_dia_cache_status(
        &self,
        Parameters(input): Parameters<GetDiaCacheStatusInput>,
    ) -> Result<Json<DiaCacheStatusOutput>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "get_dia_cache_status").entered();
        tracing::info!(dia_run_id = %input.dia_run_id, "started");

        // If dia_run_id is empty, list all cached extractions
        if input.dia_run_id.trim().is_empty() {
            let cache = self.dia_cache.lock().unwrap_or_else(|e| e.into_inner());
            let entries: Vec<DiaCacheEntry> = cache
                .list_ids()
                .into_iter()
                .map(|(id, loc)| {
                    let (location, spectrum_count, extracted_at) = match loc {
                        DiaCacheLocation::Memory {
                            spectrum_count: sc,
                            extracted_at: ts,
                        } => ("memory".to_string(), Some(sc), ts),
                        DiaCacheLocation::Disk { extracted_at: ts } => {
                            ("disk".to_string(), None, ts)
                        }
                        DiaCacheLocation::NotFound => ("not_found".to_string(), None, None),
                    };
                    DiaCacheEntry {
                        dia_run_id: id.to_string(),
                        location,
                        spectrum_count,
                        extracted_at: extracted_at.map(|t| t.to_rfc3339()),
                    }
                })
                .collect();
            let _total = entries.len();
            return Ok(Json(DiaCacheStatusOutput {
                exists: !entries.is_empty(),
                location: if entries.is_empty() {
                    "no_entries".to_string()
                } else {
                    "list".to_string()
                },
                spectrum_count: entries.iter().filter_map(|e| e.spectrum_count).next(),
                extracted_at: None,
            }));
        }

        let dia_uuid = Uuid::parse_str(&input.dia_run_id)
            .map_err(|_| mcp_err(ErrorCode::INVALID_PARAMS, "invalid dia_run_id format"))?;

        let cache = self.dia_cache.lock().unwrap_or_else(|e| e.into_inner());

        let output = match cache.status(&dia_uuid) {
            DiaCacheLocation::Memory {
                spectrum_count,
                extracted_at,
            } => DiaCacheStatusOutput {
                exists: true,
                location: "memory".to_string(),
                spectrum_count: Some(spectrum_count),
                extracted_at: extracted_at.map(|t| t.to_rfc3339()),
            },
            DiaCacheLocation::Disk { extracted_at } => DiaCacheStatusOutput {
                exists: true,
                location: "disk".to_string(),
                spectrum_count: None,
                extracted_at: extracted_at.map(|t| t.to_rfc3339()),
            },
            DiaCacheLocation::NotFound => DiaCacheStatusOutput {
                exists: false,
                location: "not_found".to_string(),
                spectrum_count: None,
                extracted_at: None,
            },
        };

        tracing::info!(location = %output.location, "completed");
        Ok(Json(output))
    }

    /// Extract XIC (Extracted Ion Chromatogram) for a peptide from an mzML file.
    #[rmcp::tool(
        name = "extract_xic",
        description = "Extract XIC (Extracted Ion Chromatogram) for a peptide from a spectrum file (mzML recommended; requires MS1+MS2 and retention time data). Generates an interactive HTML file with Plotly.js showing MS1 precursor and MS2 fragment ion chromatograms. Supports SILAC heavy-label comparison. Two modes (choose one): Mode 1 (from search) — provide run_id + scan_number to use PSM context. Mode 2 (manual REQUIRES: file_path + scan_number + peptide_sequence + charge + precursor_mz). In mode 2 only: set scan_number=0 with retention_time_min to auto-find scan. Set view='3d' for a 3D MS2 overview (output: xic3d_scan{N}.html). For fragment ion annotation (not chromatograms), see annotate_spectrum."
    )]
    fn extract_xic(
        &self,
        #[allow(unused_variables)] Parameters(input): Parameters<ExtractXicInput>,
    ) -> Result<Json<ExtractXicResult>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "extract_xic").entered();
        tracing::info!(
            scan = input.scan_number,
            peptide = input.peptide_sequence.as_deref().unwrap_or("from_run"),
            charge = input.charge.unwrap_or(0),
            precursor_mz = input.precursor_mz.unwrap_or(0.0),
            run_id = input.run_id.as_deref().unwrap_or("none"),
            file = input.file_path.as_deref().unwrap_or("none"),
            "started"
        );
        use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};

        // Allow scan_number=0 when retention_time_min is provided
        if input.scan_number == 0 && input.retention_time_min.is_none() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                "scan_number must be >= 1, or set scan_number=0 with retention_time_min for auto lookup",
            ));
        }
        if input.scan_number != 0 {
            validate_scan_number(input.scan_number)?;
        }
        if input.scan_number == 0 && input.run_id.is_some() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                "retention_time_min auto-lookup is only supported in manual mode (without run_id). With run_id, provide a valid scan_number.",
            ));
        }

        // Resolve mode: run_id or manual
        let (file_path, peptide, charge, precursor_mz, modifications) = if let Some(ref rid) =
            input.run_id
        {
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
                    "provide either 'run_id' + 'scan_number' or 'file_path' + 'scan_number' + 'peptide_sequence' + 'charge' + 'precursor_mz'",
                ));
        };

        // Resolve scan_number: if 0 and retention_time_min provided, auto-match via RT
        let resolved_scan = if input.scan_number == 0 {
            if let Some(rt) = input.retention_time_min {
                let reader = self.get_or_create_reader(&file_path)?;
                reader
                    .find_by_rt(&file_path, rt, precursor_mz, RT_AUTO_LOOKUP_TOLERANCE_MIN)
                    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?
                    .map(|(scan, _)| scan)
                    .ok_or_else(|| {
                        mcp_err(
                            ErrorCode::INVALID_PARAMS,
                            format!("No MS2 scan near RT={rt:.2}min mz={precursor_mz:.4}"),
                        )
                    })?
            } else {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "scan_number is 0: provide retention_time_min for auto lookup",
                ));
            }
        } else {
            input.scan_number
        };

        let params = protein_copilot_xic::ExtractionParams {
            mz_tolerance: input.extraction_tolerance.unwrap_or(MassTolerance {
                value: 20.0,
                unit: ToleranceUnit::Ppm,
            }),
            n_cycles: input.n_cycles.unwrap_or(5),
            top_n_ions: input.top_n_ions.unwrap_or(usize::MAX),
            label_type: input.label_type.clone(),
            intensity_rule: input
                .intensity_rule
                .unwrap_or(protein_copilot_xic::IntensityRule::MaxInWindow),
        };

        // Extract XIC via unified (index-planned) path
        let cached_reader = self.get_or_create_reader(&file_path)?;
        let unified_result = protein_copilot_xic::extract::extract_xic_unified(
            cached_reader.as_ref(),
            &file_path,
            resolved_scan,
            &peptide,
            charge,
            precursor_mz,
            &modifications,
            &params,
            20.0,
        )
        .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?;
        // 3D view: annotate every MS2 scan in the window and render 3D HTML.
        if matches!(input.view, Some(XicView::ThreeD)) {
            let plotly_mode = input
                .plotly_mode
                .unwrap_or(protein_copilot_xic::PlotlyMode::Cdn);
            let out_path = input
                .output_path
                .clone()
                .map(PathBuf::from)
                .unwrap_or_else(|| {
                    PathBuf::from(format!("output/xic3d_scan{}.html", resolved_scan))
                });

            let data = protein_copilot_report::xic3d_build::build_xic3d_data(
                &unified_result.raw_scans.ms2_scans,
                &peptide,
                charge,
                precursor_mz,
                &modifications,
                &params.mz_tolerance,
                resolved_scan,
                &file_path.to_string_lossy(),
            )
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?;

            protein_copilot_report::xic3d_visualize::render_xic_3d(
                &data,
                &out_path,
                plotly_mode,
                input.max_peaks_per_scan_3d,
            )
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?;

            let annotated = data.scans.len();
            let target_matched = data
                .scans
                .iter()
                .find(|s| s.is_target)
                .map(|s| s.annotation.matched_ions);

            tracing::info!(scans = annotated, "completed (3d)");
            return Ok(Json(ExtractXicResult {
                output_path: out_path.to_string_lossy().to_string(),
                ms2_scan_count: annotated,
                light_trace_count: 0,
                heavy_trace_count: 0,
                has_ms1_xic: false,
                annotated_scan_count: Some(annotated),
                target_matched_ions: target_matched,
                summary: format!(
                    "3D MS2 view for {} ({}+): {} MS2 scans annotated in window",
                    peptide, charge, annotated
                ),
            }));
        }
        let xic_data = unified_result.xic_data;

        // Render HTML
        let out_path = input
            .output_path
            .map(PathBuf::from)
            .unwrap_or_else(|| default_output_path(&format!("xic_scan{}.html", resolved_scan)));
        let out_path = absolutize(&out_path);
        ensure_parent_dir(&out_path)?;

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

        tracing::info!("completed");
        Ok(Json(ExtractXicResult {
            output_path: out_path.to_string_lossy().to_string(),
            ms2_scan_count: ms2_count,
            light_trace_count: xic_data.fragment_xic_traces.len(),
            heavy_trace_count: xic_data.heavy_fragment_xic_traces.len(),
            has_ms1_xic: xic_data.ms1_precursor_xic.is_some(),
            summary,
            annotated_scan_count: None,
            target_matched_ions: None,
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
        let _span = tracing::info_span!("mcp_tool", name = "import_search_results").entered();
        tracing::info!(
            result_file = %input.result_file,
            mzml_dir = %input.mzml_dir,
            format = %input.format,
            filter_qvalue = input.filter_qvalue,
            rt_tolerance_min = input.rt_tolerance_min,
            run_filter = input.run_filter.as_deref().unwrap_or("none"),
            "started"
        );
        // Keep: feeds RunMetadata.duration_sec / MCP client
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
        if input.rt_tolerance_min < 0.0 || !input.rt_tolerance_min.is_finite() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                "rt_tolerance_min must be a non-negative finite number",
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
            UnimodDb::from_xml(Path::new(xml_path))
                .map_err(|e| mcp_err(ErrorCode::INVALID_PARAMS, format!("unimod.xml error: {e}")))?
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
            "pfind_tsv" => ImportFormat::PFindTsv,
            other => {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    format!(
                        "unknown format: '{other}'. Supported: auto, custom_json, diann_parquet, pfind_tsv"
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
                     Supported formats: custom_json, diann_parquet, pfind_tsv",
                ));
            }
            ImportFormat::PFindTsv => PFindTsvParser
                .parse(&result_path, &unimod)
                .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?,
        };

        // Apply run_filter for formats that don't handle it internally (e.g. pFind TSV)
        if let Some(ref filter) = input.run_filter {
            let before = psms.len();
            psms.retain(|p| p.raw_name == *filter);
            if before != psms.len() {
                tracing::info!(
                    run_filter = %filter,
                    before = before,
                    after = psms.len(),
                    "applied run_filter to imported PSMs"
                );
            }
        }

        if psms.is_empty() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                "no PSMs parsed from the result file (check format and filters)",
            ));
        }

        // Scan matching — skip if all PSMs already have scan numbers (e.g. pFind TSV)
        let all_scans_present = psms.iter().all(|p| p.matched_scan.is_some());

        let match_report = if all_scans_present {
            // pFind TSV already has scan numbers — build MatchReport directly
            let mut per_file = std::collections::HashMap::new();
            for psm in &psms {
                let entry = per_file.entry(psm.raw_name.clone()).or_insert(
                    protein_copilot_result_import::FileMatchStats {
                        total: 0,
                        matched: 0,
                        ms2_count: 0,
                    },
                );
                entry.total += 1;
                entry.matched += 1;
            }
            tracing::info!(
                psm_count = psms.len(),
                "all PSMs have scan numbers — skipping RT-based scan matching"
            );
            protein_copilot_result_import::MatchReport {
                total_psms: psms.len(),
                matched: psms.len(),
                unmatched: 0,
                median_rt_delta_min: 0.0,
                max_rt_delta_min: 0.0,
                per_file,
            }
        } else {
            // Normal path: RT-based scan matching
            let config = ScanMatcherConfig {
                rt_tolerance_min: input.rt_tolerance_min,
                mzml_dir: mzml_dir.clone(),
            };
            match_scans(&mut psms, &config, &|path| {
                protein_copilot_spectrum_io::create_indexed_reader(path).map_err(|e| {
                    protein_copilot_result_import::ResultImportError::SpectrumIo(format!(
                        "failed to open {}: {e}",
                        path.display(),
                    ))
                })
            })
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?
        };

        if match_report.matched == 0 {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!(
                    "parsed {} PSMs but matched 0 to mzML scans — check that mzML files \
                     correspond to the search results and RT values are correct",
                    match_report.total_psms
                ),
            ));
        }

        // Collect actual mzML paths used (for downstream annotate_spectrum/extract_xic).
        // pFind imports skip RT scan matching (which validates raw names), so
        // confine names to mzml_dir here too, rejecting path traversal. Sort the
        // unique names so mzml_files (and thus metadata.input_files) is
        // deterministic across imports rather than HashSet-ordered.
        let mut raw_names: Vec<&str> = psms
            .iter()
            .map(|p| p.raw_name.as_str())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();
        raw_names.sort_unstable();
        let mut mzml_files: Vec<PathBuf> = Vec::with_capacity(raw_names.len());
        for raw in raw_names {
            validate_raw_name(raw)
                .map_err(|e| mcp_err(ErrorCode::INVALID_PARAMS, e.to_string()))?;
            let mzml = mzml_dir.join(format!("{raw}.mzML"));
            mzml_files.push(if mzml.exists() {
                mzml
            } else {
                mzml_dir.join(format!("{raw}.mzml"))
            });
        }

        if mzml_files.len() > 1 {
            tracing::warn!(
                "imported PSMs span {} raw files — annotate_spectrum/extract_xic currently \
                 only use the first file. Use run_filter to limit to a single raw file \
                 for reliable downstream annotation.",
                mzml_files.len()
            );
        }

        // Convert to SearchResult
        let format_name = match format {
            ImportFormat::CustomJson => "custom_json",
            ImportFormat::DiannParquet => "diann_parquet",
            ImportFormat::PFindSpectra => "pfind_spectra",
            ImportFormat::PFindTsv => "pfind_tsv",
        };
        let (mut search_result, import_result) =
            build_search_result(&psms, match_report, format_name, mzml_files);

        // Fix metadata: set status to Completed and record duration
        let duration = start.elapsed().as_secs_f64();
        search_result.metadata.status = protein_copilot_core::run_metadata::RunStatus::Completed;
        search_result.metadata.duration_sec = Some(duration);

        // Store in run_cache
        let run_id = search_result.run_id;
        {
            let mut cache = self.run_cache.lock().unwrap_or_else(|e| e.into_inner());
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
                        error_category: None,
                        has_diagnostics: false,
                    },
                    result: Some(search_result.clone()),
                    handle: None,
                    diagnostics: None,
                    params_used: None,
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

        tracing::info!(psm_count = import_result.imported_psm_count, "completed");
        Ok(Json(import_result))
    }

    // -----------------------------------------------------------------------
    // Protein Inference
    // -----------------------------------------------------------------------

    /// Run protein inference on search results.
    #[rmcp::tool(
        name = "infer_proteins",
        description = "Run protein inference on search results from run_search (provide run_id). Performs parsimony analysis, razor peptide assignment, protein-level FDR, and optional sequence coverage. Input: run_id from a previous search or direct SearchResult. Returns protein groups with scores, q-values, and peptide assignments."
    )]
    fn infer_proteins(
        &self,
        Parameters(input): Parameters<InferProteinsInput>,
    ) -> Result<Json<InferenceResult>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "infer_proteins").entered();
        tracing::info!(
            run_id = input.run_id.as_deref().unwrap_or("direct"),
            q_value_threshold = input.q_value_threshold.unwrap_or(0.01),
            fasta_path = input.fasta_path.as_deref().unwrap_or("none"),
            "started"
        );
        // Validate q_value_threshold
        if let Some(q) = input.q_value_threshold {
            if q.is_nan() || q.is_infinite() || !(0.0..=1.0).contains(&q) {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    format!("q_value_threshold must be between 0.0 and 1.0, got {q}"),
                ));
            }
        }

        // Validate fasta_path upfront
        if let Some(ref fp) = input.fasta_path {
            validate_file_path(fp)?;
        }

        let result = self.get_result(&input.result, &input.run_id)?;

        // Build peptide-protein map
        let map = protein_copilot_protein_inference::mapper::build_peptide_protein_map(
            &result.psms,
            input.q_value_threshold,
        )
        .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))?;

        // Run parsimony
        let mut groups = protein_copilot_protein_inference::parsimony::run_parsimony(&map)
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))?;

        // Assign razor peptides
        let razor_map =
            protein_copilot_protein_inference::razor::assign_razor_peptides(&mut groups, &map);

        // Calculate protein-level FDR (picked-protein)
        let fdr_result = protein_copilot_fdr::protein_fdr::calculate_protein_fdr(&groups);

        let mut final_groups = match fdr_result {
            Ok(fdr) => fdr.groups,
            Err(e) => {
                tracing::warn!(
                    "Protein FDR calculation failed: {e}. Returning groups without q-values."
                );
                groups
            }
        };

        // Optional: sequence coverage
        if let Some(fasta_path) = &input.fasta_path {
            let fasta_sequences = load_fasta_sequences(fasta_path)
                .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))?;
            protein_copilot_protein_inference::coverage::calculate_coverage(
                &mut final_groups,
                &fasta_sequences,
            );
        }

        let total_target = final_groups.iter().filter(|g| !g.is_decoy).count() as u64;
        let total_decoy = final_groups.iter().filter(|g| g.is_decoy).count() as u64;
        let groups_at_1pct = final_groups
            .iter()
            .filter(|g| !g.is_decoy && g.q_value.is_some_and(|q| q <= FDR_1PCT_THRESHOLD))
            .count() as u64;
        let unique_peptides_used = map.peptide_to_proteins.len() as u64;

        tracing::info!(groups = groups_at_1pct, "completed");
        Ok(Json(InferenceResult {
            groups: final_groups,
            razor_map,
            total_target_groups: total_target,
            total_decoy_groups: total_decoy,
            groups_at_1pct_fdr: groups_at_1pct,
            unique_peptides_used,
        }))
    }

    // -----------------------------------------------------------------------
    // Composite: prepare_search
    // -----------------------------------------------------------------------

    /// One-shot preparation for a database search: read spectra, recommend
    /// parameters, and resolve a FASTA database — all in a single call.
    #[rmcp::tool(
        name = "prepare_search",
        description = "One-shot search preparation: reads spectrum files, recommends search parameters, and resolves a FASTA database. Combines read_spectra + recommend_params + download_database into a single call. Provide either 'database_path' (direct FASTA path) or 'organism' (e.g. 'human', 'mouse', 'E.coli', '小鼠') for auto-resolution. Returns ready-to-use SearchParams that can be passed directly to run_search."
    )]
    async fn prepare_search(
        &self,
        Parameters(input): Parameters<PrepareSearchInput>,
    ) -> Result<Json<PrepareSearchOutput>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "prepare_search");
        let _enter = _span.enter();
        tracing::info!(
            files = input.input_files.len(),
            organism = input.organism.as_deref().unwrap_or("none"),
            database = input.database_path.as_deref().unwrap_or("auto"),
            engine = input.engine.as_deref().unwrap_or("auto"),
            "started"
        );
        drop(_enter);
        // 1. Validate input_files
        if input.input_files.is_empty() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                "input_files is empty — provide at least one spectrum file path",
            ));
        }
        for f in &input.input_files {
            validate_file_path(f)?;
        }

        // 2. Read spectrum summary from first file
        let first_path = Path::new(&input.input_files[0]);
        let reader = self.get_or_create_reader(first_path)?;
        let summary = reader
            .read_summary(first_path)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

        // 3. Recommend parameters
        let recommender = ParamRecommender;
        let decision = recommender
            .recommend(&summary, input.hints.as_ref())
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        let mut params = decision.decision;

        // 4. Set engine if specified
        if let Some(ref engine) = input.engine {
            params.engine = Some(engine.clone());
        }

        // 5. Resolve database
        let database_info = if let Some(ref db_path) = input.database_path {
            // Direct path takes priority
            let p = Path::new(db_path);
            if !p.exists() {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    format!("database_path does not exist: {db_path}"),
                ));
            }
            params.database_path = db_path.clone();
            None
        } else if let Some(ref organism) = input.organism {
            let db_id = organism_to_database_id(organism).ok_or_else(|| {
                mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    format!(
                        "Could not resolve organism '{}' to a database. \
                         Supported: human, mouse, E.coli, yeast, arabidopsis, \
                         or use database_path for a custom FASTA.",
                        organism
                    ),
                )
            })?;

            let cache_dir = default_cache_dir(&input.cache_dir);

            // Check if already downloaded
            let databases = protein_copilot_fasta_db::list_databases(&cache_dir)
                .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))?;

            let already_downloaded = databases.iter().any(|d| {
                d.id == db_id
                    && matches!(
                        d.status,
                        protein_copilot_fasta_db::DownloadStatus::Downloaded { .. }
                    )
            });

            let dl_result = protein_copilot_fasta_db::download_database(db_id, &cache_dir, false)
                .await
                .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))?;

            params.database_path = dl_result.path.clone();

            Some(PreparedDatabaseInfo {
                id: dl_result.id,
                path: dl_result.path,
                protein_count: dl_result.protein_count,
                freshly_downloaded: !already_downloaded,
            })
        } else {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                "Provide either 'database_path' (direct FASTA path) \
                 or 'organism' (e.g. 'human', 'mouse') for database resolution.",
            ));
        };

        // 6. Validate params
        params
            .validate()
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

        // 7. Verify database file exists
        let db_path = Path::new(&params.database_path);
        if !db_path.exists() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!(
                    "Database file not found after resolution: {}. \
                     Use list_databases to see available databases, \
                     or download_database to fetch one.",
                    params.database_path
                ),
            ));
        }

        tracing::info!(engine = params.engine.as_deref().unwrap_or("auto"), database = %params.database_path, "completed");
        Ok(Json(PrepareSearchOutput {
            params,
            reasoning: decision.explanation,
            confidence: decision.confidence,
            alternatives: decision.alternatives,
            evidence: decision.evidence,
            spectra_summary: summary,
            database_info,
        }))
    }

    // -----------------------------------------------------------------------
    // FASTA Database Management
    // -----------------------------------------------------------------------

    /// List all available FASTA databases and their cache status.
    #[rmcp::tool(
        name = "list_databases",
        description = "List all built-in FASTA protein databases (Human, Mouse, E.coli, Yeast, Arabidopsis, cRAP contaminants) with download status. Shows which databases are cached locally and which are available for download."
    )]
    fn list_databases(
        &self,
        Parameters(input): Parameters<ListDatabasesInput>,
    ) -> Result<Json<ListDatabasesOutput>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "list_databases").entered();
        tracing::info!("started");
        let cache_dir = default_cache_dir(&input.cache_dir);
        protein_copilot_fasta_db::list_databases(&cache_dir)
            .map(|dbs| {
                tracing::info!(count = dbs.len(), "completed");
                Json(ListDatabasesOutput { databases: dbs })
            })
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))
    }

    /// Download a FASTA database by ID.
    #[rmcp::tool(
        name = "download_database",
        description = "Download a FASTA protein database by ID (e.g. 'human_swissprot', 'mouse_swissprot', 'ecoli_swissprot', 'yeast_swissprot', 'arabidopsis_swissprot', 'crap'). Downloads from UniProt via HTTPS and caches locally. Returns the local file path for use as database_path in search parameters. Use list_databases first to see available options."
    )]
    async fn download_database(
        &self,
        Parameters(input): Parameters<DownloadDatabaseInput>,
    ) -> Result<Json<protein_copilot_fasta_db::DownloadDatabaseResult>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "download_database");
        let _enter = _span.enter();
        tracing::info!(database_id = %input.database_id, "started");
        drop(_enter);
        let cache_dir = default_cache_dir(&input.cache_dir);
        let force = input.force.unwrap_or(false);
        let result =
            protein_copilot_fasta_db::download_database(&input.database_id, &cache_dir, force)
                .await
                .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))?;
        tracing::info!(path = %result.path, "completed");
        Ok(Json(result))
    }

    /// Get detailed info about a cached FASTA database.
    #[rmcp::tool(
        name = "get_database_info",
        description = "Get detailed information about a downloaded FASTA database: protein count, file size, SHA256 hash, download date, and first 5 protein accessions. The database must be downloaded first using download_database."
    )]
    fn get_database_info(
        &self,
        Parameters(input): Parameters<GetDatabaseInfoInput>,
    ) -> Result<Json<protein_copilot_fasta_db::DatabaseInfo>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "get_database_info").entered();
        tracing::info!(database_id = %input.database_id, "started");
        let cache_dir = default_cache_dir(&input.cache_dir);
        let info = protein_copilot_fasta_db::get_database_info(&input.database_id, &cache_dir)
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))?;
        tracing::info!("completed");
        Ok(Json(info))
    }

    /// Get diagnostic report for a search run.
    #[rmcp::tool(
        name = "diagnose_search",
        description = "Get diagnostic report for a search started by run_search. Works for both failed searches (error analysis) and completed searches (quality assessment). Returns stage metrics, detected anomalies, and repair suggestions. Call after get_search_status shows the search has finished (status is Completed, Failed, or Cancelled). Use has_diagnostics=true from get_search_status to confirm diagnostics are available."
    )]
    fn diagnose_search(
        &self,
        Parameters(input): Parameters<DiagnoseSearchInput>,
    ) -> Result<Json<DiagnoseSearchOutput>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "diagnose_search").entered();
        tracing::info!(run_id = %input.run_id, "started");
        let id = Uuid::parse_str(&input.run_id).map_err(|_| {
            mcp_err(
                ErrorCode::INVALID_PARAMS,
                "invalid run_id format — expected UUID",
            )
        })?;

        let cache = self.run_cache.lock().unwrap_or_else(|e| e.into_inner());

        let state = cache.get(&id).ok_or_else(|| {
            mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!("run_id '{}' not found in cache", input.run_id),
            )
        })?;

        if state.progress.status == "Running" {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                "Search is still running. Wait for completion before diagnosing.",
            ));
        }

        let diag = state.diagnostics.as_ref().ok_or_else(|| {
            mcp_err(
                ErrorCode::INTERNAL_ERROR,
                "No diagnostics available for this run. The search may have been started before the diagnostics feature was added.",
            )
        })?;

        tracing::info!("completed");
        Ok(Json(DiagnoseSearchOutput {
            run_id: input.run_id,
            overall_status: state.progress.status.clone(),
            error_category: diag.error_category.clone(),
            failure_stage: diag.failure_stage.clone(),
            error_detail: diag.error_detail.clone(),
            stages: diag.stages.clone(),
            anomalies: diag.anomalies.clone(),
            suggestions: diag.suggestions.clone(),
            total_elapsed_sec: diag.total_elapsed_sec,
        }))
    }

    // -----------------------------------------------------------------------
    // Entrapment analysis tools
    // -----------------------------------------------------------------------

    #[rmcp::tool(
        name = "classify_entrapment_hits",
        description = "Classify trap-database PSM hits by homology to target proteome. Reads search results, applies target/trap rules from YAML config, digests target FASTA, and classifies each trap PSM as L0-L4. Optionally traces fragment ion provenance when mzml_dir is provided. Outputs classified.tsv, razor_errors.tsv, run_metadata.json, and entrapment_report.html. Follow up with analyze_entrapment_stats for detailed statistics, find_similar_targets to investigate specific PSMs, or annotate_provenance for fragment ion mirror plots."
    )]
    fn classify_entrapment_hits(
        &self,
        Parameters(input): Parameters<ClassifyEntrapmentHitsInput>,
    ) -> Result<Json<ClassifyEntrapmentOutput>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "classify_entrapment_hits").entered();
        tracing::info!(
            results_file = %input.results_file,
            config_file = %input.config_file,
            target_fasta = %input.target_fasta,
            mzml_dir = input.mzml_dir.as_deref().unwrap_or("none"),
            output_dir = input.output_dir.as_deref().unwrap_or("default"),
            "started"
        );
        use protein_copilot_entrapment_analysis::{
            config::EntrapmentConfig,
            loader::{self, ResultFormat},
            output::{self, EntrapmentRunMetadata},
            EntrapmentAnalyzer,
        };

        let config = EntrapmentConfig::from_yaml(std::path::Path::new(&input.config_file))
            .map_err(|e| ErrorData::new(ErrorCode::INVALID_PARAMS, format!("{e}"), None))?;

        let format = match input.format.as_deref() {
            Some("diann_parquet") => ResultFormat::DiannParquet,
            Some("generic_tsv") => ResultFormat::GenericTsv,
            Some("pfind_tsv") => ResultFormat::PFindTsv,
            _ => ResultFormat::from_path(std::path::Path::new(&input.results_file))
                .map_err(|e| ErrorData::new(ErrorCode::INVALID_PARAMS, format!("{e}"), None))?,
        };

        let psms = loader::load_psms(std::path::Path::new(&input.results_file), &format, None)
            .map_err(|e| ErrorData::new(ErrorCode::INTERNAL_ERROR, format!("{e}"), None))?;

        let analyzer =
            EntrapmentAnalyzer::new(config.clone(), std::path::Path::new(&input.target_fasta))
                .map_err(|e| ErrorData::new(ErrorCode::INTERNAL_ERROR, format!("{e}"), None))?;

        let mut classified = analyzer
            .classify_all(&psms)
            .map_err(|e| ErrorData::new(ErrorCode::INTERNAL_ERROR, format!("{e}"), None))?;

        // Provenance tracing (optional)
        if let Some(mzml_dir) = &input.mzml_dir {
            use protein_copilot_entrapment_analysis::trace_provenance_batch;
            let mzml_path = std::path::Path::new(mzml_dir);
            match trace_provenance_batch(&mut classified, mzml_path, &config) {
                Ok(count) => {
                    tracing::info!("Provenance traced for {} PSMs", count);
                }
                Err(e) => {
                    tracing::warn!("Provenance tracing failed: {}", e);
                }
            }
        }

        let summary = analyzer.summary(&classified);

        let out_dir = std::path::PathBuf::from(
            input
                .output_dir
                .unwrap_or_else(|| "output/entrapment".to_string()),
        );
        std::fs::create_dir_all(&out_dir).map_err(|e| {
            ErrorData::new(
                ErrorCode::INTERNAL_ERROR,
                format!("create output dir: {e}"),
                None,
            )
        })?;

        output::write_classified_tsv(&classified, &out_dir.join("classified.tsv"))
            .map_err(|e| ErrorData::new(ErrorCode::INTERNAL_ERROR, format!("{e}"), None))?;

        output::write_razor_errors_tsv(&classified, &out_dir.join("razor_errors.tsv"))
            .map_err(|e| ErrorData::new(ErrorCode::INTERNAL_ERROR, format!("{e}"), None))?;

        let metadata = EntrapmentRunMetadata {
            tool_version: env!("CARGO_PKG_VERSION").to_string(),
            run_timestamp: chrono::Utc::now().to_rfc3339(),
            input_file: input.results_file.clone(),
            input_sha256: output::file_sha256(std::path::Path::new(&input.results_file))
                .unwrap_or_else(|_| "unknown".to_string()),
            fasta_file: input.target_fasta.clone(),
            fasta_sha256: output::file_sha256(std::path::Path::new(&input.target_fasta))
                .unwrap_or_else(|_| "unknown".to_string()),
            config_snapshot: serde_json::to_value(&config).unwrap_or_default(),
            total_psms: summary.total_psms,
            trap_psms: summary.trap_psms,
            level_counts: summary.level_counts.clone(),
        };
        output::write_run_metadata(&metadata, &out_dir.join("run_metadata.json"))
            .map_err(|e| ErrorData::new(ErrorCode::INTERNAL_ERROR, format!("{e}"), None))?;

        protein_copilot_entrapment_analysis::report::render_report(
            &summary,
            &classified,
            &out_dir.join("entrapment_report.html"),
        )
        .map_err(|e| ErrorData::new(ErrorCode::INTERNAL_ERROR, format!("{e}"), None))?;

        let output = ClassifyEntrapmentOutput {
            total_psms: summary.total_psms,
            target_psms: summary.target_psms,
            trap_psms: summary.trap_psms,
            ambiguous_psms: summary.ambiguous_psms,
            level_counts: EntrapmentLevelCountsOutput {
                l0: summary.level_counts.l0,
                l1: summary.level_counts.l1,
                l2: summary.level_counts.l2,
                l3: summary.level_counts.l3,
                l4: summary.level_counts.l4,
            },
            top_razor_families: summary
                .top_razor_families
                .iter()
                .map(|f| EntrapmentRazorFamilyOutput {
                    family: f.family.clone(),
                    count: f.count,
                    example_peptide: f.example_peptide.clone(),
                    example_trap_protein: f.example_trap_protein.clone(),
                    example_target_protein: f.example_target_protein.clone(),
                })
                .collect(),
        };

        tracing::info!("completed");
        Ok(Json(output))
    }

    #[rmcp::tool(
        name = "analyze_entrapment_stats",
        description = "Get detailed statistics from a classified entrapment TSV file. Returns level distribution, protein family clusters, and delta-mass analysis. Use after classify_entrapment_hits to interpret results."
    )]
    fn analyze_entrapment_stats(
        &self,
        Parameters(input): Parameters<AnalyzeEntrapmentStatsInput>,
    ) -> Result<Json<AnalyzeEntrapmentStatsOutput>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "analyze_entrapment_stats").entered();
        tracing::info!(classified_file = %input.classified_file, "started");
        let path = std::path::Path::new(&input.classified_file);
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(path)
            .map_err(|e| {
                ErrorData::new(
                    ErrorCode::INVALID_PARAMS,
                    format!("cannot read file: {e}"),
                    None,
                )
            })?;

        let headers = rdr
            .headers()
            .map_err(|e| {
                ErrorData::new(
                    ErrorCode::INVALID_PARAMS,
                    format!("cannot read headers: {e}"),
                    None,
                )
            })?
            .clone();

        use protein_copilot_entrapment_analysis::output::columns;

        let level_idx = headers.iter().position(|h| h == columns::LEVEL);
        let delta_idx = headers.iter().position(|h| h == columns::DELTA_MASS_DA);
        let target_protein_idx = headers
            .iter()
            .position(|h| h == columns::BEST_TARGET_PROTEIN);

        let mut level_counts = std::collections::HashMap::<String, usize>::new();
        let mut delta_masses: Vec<f64> = Vec::new();
        let mut protein_families = std::collections::HashMap::<String, usize>::new();
        let mut total = 0usize;

        for result in rdr.records() {
            let record = result.map_err(|e| {
                ErrorData::new(ErrorCode::INTERNAL_ERROR, format!("parse row: {e}"), None)
            })?;
            total += 1;

            if let Some(idx) = level_idx {
                if let Some(level) = record.get(idx) {
                    *level_counts.entry(level.to_string()).or_insert(0) += 1;
                }
            }
            if let Some(idx) = delta_idx {
                if let Some(delta_str) = record.get(idx) {
                    if let Ok(d) = delta_str.parse::<f64>() {
                        delta_masses.push(d);
                    }
                }
            }
            if let Some(idx) = target_protein_idx {
                if let Some(target_protein) = record.get(idx) {
                    if !target_protein.is_empty() {
                        let family = target_protein
                            .split('|')
                            .nth(2)
                            .and_then(|s| s.split('_').next())
                            .unwrap_or(target_protein)
                            .to_string();
                        *protein_families.entry(family).or_insert(0) += 1;
                    }
                }
            }
        }

        let mut top_families: Vec<_> = protein_families.into_iter().collect();
        top_families.sort_by_key(|b| std::cmp::Reverse(b.1));
        top_families.truncate(20);

        let stats = AnalyzeEntrapmentStatsOutput {
            total_classified: total,
            level_distribution: level_counts,
            delta_mass_stats: DeltaMassStats {
                count: delta_masses.len(),
                min: if delta_masses.is_empty() {
                    0.0
                } else {
                    delta_masses.iter().copied().fold(f64::INFINITY, f64::min)
                },
                max: if delta_masses.is_empty() {
                    0.0
                } else {
                    delta_masses
                        .iter()
                        .copied()
                        .fold(f64::NEG_INFINITY, f64::max)
                },
                mean: if delta_masses.is_empty() {
                    0.0
                } else {
                    delta_masses.iter().sum::<f64>() / delta_masses.len() as f64
                },
            },
            top_protein_families: top_families,
        };

        tracing::info!("completed");
        Ok(Json(stats))
    }

    #[rmcp::tool(
        name = "find_similar_targets",
        description = "Find similar target peptides for a given sequence — use after classify_entrapment_hits to investigate individual trap PSMs. Digests the target FASTA, compares the query peptide against target peptides using edit distance (Hamming for same-length, Levenshtein for cross-length). Returns closest matches with mass differences and substitution type annotations."
    )]
    fn find_similar_targets(
        &self,
        Parameters(input): Parameters<FindSimilarTargetsInput>,
    ) -> Result<Json<FindSimilarTargetsOutput>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "find_similar_targets").entered();
        tracing::info!(
            peptide = %input.peptide,
            target_fasta = %input.target_fasta,
            max_mismatches = input.max_mismatches.unwrap_or(2),
            "started"
        );
        use protein_copilot_entrapment_analysis::{
            config::SimilarityConfig,
            digest::TargetDigestIndex,
            similarity::classify_single,
            types::{PsmGroup, UnifiedPsm},
        };

        let max_mm = input.max_mismatches.unwrap_or(2);
        let sim_config = SimilarityConfig {
            max_mismatches: max_mm,
            ..SimilarityConfig::default()
        };

        let index = TargetDigestIndex::from_fasta(
            std::path::Path::new(&input.target_fasta),
            sim_config.max_missed_cleavages,
            sim_config.max_mismatches,
        )
        .map_err(|e| ErrorData::new(ErrorCode::INTERNAL_ERROR, format!("{e}"), None))?;

        let psm = UnifiedPsm {
            peptide: input.peptide.clone(),
            charge: None,
            precursor_mz: None,
            retention_time: None,
            rt_start: None,
            rt_stop: None,
            scan_number: None,
            spectrum_file: None,
            protein_ids: String::new(),
            q_value: None,
            modifications: Vec::new(),
        };

        let result = classify_single(&psm, PsmGroup::Trap, &index, &sim_config);

        let output = FindSimilarTargetsOutput {
            peptide: input.peptide,
            level: result.level.as_str().to_string(),
            best_target_peptide: result.best_target_peptide,
            best_target_protein: result.best_target_protein,
            mismatches: result.mismatches,
            delta_mass_da: result.delta_mass_da,
            diff_positions: result.diff_positions,
            index_size: index.len(),
            substitution_type: Some(result.substitution_type.to_string()),
            edit_distance: result.edit_distance,
            alignment_detail: result.alignment_detail,
        };

        tracing::info!(matches = output.index_size, "completed");
        Ok(Json(output))
    }

    // -----------------------------------------------------------------------
    // Provenance annotation (single PSM → mirror plot)
    // -----------------------------------------------------------------------

    #[rmcp::tool(
        name = "annotate_provenance",
        description = "Annotate a single spectrum with fragment ion provenance analysis — part of the entrapment workflow (use after classify_entrapment_hits). Generates a mirror plot HTML file showing which peaks come from the trap peptide, target peptide, both (shared), or neither (unassigned). Supports mzML, MGF, and PFB formats."
    )]
    fn annotate_provenance(
        &self,
        Parameters(input): Parameters<AnnotateProvenanceInput>,
    ) -> Result<Json<AnnotateProvenanceOutput>, ErrorData> {
        let _span = tracing::info_span!("mcp_tool", name = "annotate_provenance").entered();
        tracing::info!(
            scan = input.scan_number,
            trap_sequence = %input.trap_sequence,
            target_sequence = %input.target_sequence,
            file = %input.file_path,
            fragment_tolerance_ppm = input.fragment_tolerance_ppm,
            "started"
        );
        use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};
        use protein_copilot_entrapment_analysis::mirror_plot::render_mirror_plot;
        use protein_copilot_entrapment_analysis::provenance::trace_provenance;

        // Validate inputs
        validate_file_path(&input.file_path)?;
        validate_scan_number(input.scan_number)?;

        if input.trap_sequence.trim().is_empty() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                "trap_sequence cannot be empty",
            ));
        }
        if input.fragment_tolerance_ppm <= 0.0 || !input.fragment_tolerance_ppm.is_finite() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!(
                    "fragment_tolerance_ppm must be a positive finite number, got {}",
                    input.fragment_tolerance_ppm
                ),
            ));
        }
        if input.max_fragment_charge < 1 {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!(
                    "max_fragment_charge must be >= 1, got {}",
                    input.max_fragment_charge
                ),
            ));
        }

        // Read spectrum via cached indexed reader
        let file_path = Path::new(&input.file_path);
        let reader = self.get_or_create_reader(file_path)?;
        let spectrum = reader
            .read_spectrum(file_path, input.scan_number)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

        // Build tolerance
        let tolerance = MassTolerance {
            value: input.fragment_tolerance_ppm,
            unit: ToleranceUnit::Ppm,
        };

        // Run provenance analysis
        let provenance = trace_provenance(
            &spectrum.mz_array,
            &spectrum.intensity_array,
            &input.trap_sequence,
            &input.target_sequence,
            &input.modifications,
            &tolerance,
            input.max_fragment_charge,
        );

        // Generate output path
        let output_path = input.output_path.map(PathBuf::from).unwrap_or_else(|| {
            default_output_path(&format!("provenance_scan{}.html", input.scan_number))
        });
        let output_path = absolutize(&output_path);
        ensure_parent_dir(&output_path)?;

        // Render mirror plot
        render_mirror_plot(&provenance, &output_path).map_err(|e| {
            mcp_err(
                ErrorCode::INTERNAL_ERROR,
                format!("failed to write mirror plot: {e}"),
            )
        })?;

        tracing::info!("completed");
        Ok(Json(AnnotateProvenanceOutput {
            output_path: output_path.display().to_string(),
            trap_sequence: provenance.trap_sequence,
            target_sequence: provenance.target_sequence,
            trap_matched_count: provenance.trap_matched_count,
            target_matched_count: provenance.target_matched_count,
            shared_count: provenance.shared_count,
            unassigned_count: provenance.unassigned_count,
            shared_ratio: provenance.shared_ratio,
            is_chimeric: provenance.shared_ratio > input.chimera_threshold,
            total_peaks: provenance.annotated_peaks.len(),
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::search_params::{
        DecoyStrategy, Enzyme, MassTolerance, ToleranceUnit,
    };
    use protein_copilot_core::spectrum::MsLevel;

    /// Absolute path to a repo test fixture, independent of the working dir.
    fn fixture(name: &str) -> String {
        format!(
            "{}/../../tests/fixtures/{}",
            env!("CARGO_MANIFEST_DIR"),
            name
        )
    }

    #[test]
    fn user_input_errors_map_to_invalid_params() {
        use protein_copilot_core::error::CoreError;
        let unsupported = CoreError::UnsupportedFormat {
            format: "py".to_string(),
            supported: vec!["mzML".to_string(), "mgf".to_string()],
        };
        assert_eq!(core_err_code(&unsupported), ErrorCode::INVALID_PARAMS);
        let not_found = CoreError::FileNotFound {
            path: std::path::PathBuf::from("/no/such"),
        };
        assert_eq!(core_err_code(&not_found), ErrorCode::INVALID_PARAMS);
        let engine = CoreError::SearchEngineError {
            engine: "Sage".to_string(),
            detail: "boom".to_string(),
            suggestion: "retry".to_string(),
        };
        assert_eq!(core_err_code(&engine), ErrorCode::INTERNAL_ERROR);
    }

    #[test]
    fn output_path_respects_env() {
        std::env::remove_var("PROTEIN_OUTPUT_DIR");
        assert_eq!(
            default_output_path("a.html"),
            std::path::PathBuf::from("output").join("a.html")
        );
        std::env::set_var("PROTEIN_OUTPUT_DIR", "/tmp/pc_out");
        assert_eq!(
            default_output_path("a.html"),
            std::path::PathBuf::from("/tmp/pc_out").join("a.html")
        );
        std::env::remove_var("PROTEIN_OUTPUT_DIR");
    }

    #[test]
    fn absolutize_makes_relative_absolute() {
        let abs = absolutize(std::path::Path::new("output/x.html"));
        assert!(abs.is_absolute(), "relative path should become absolute");
        assert!(abs.ends_with("output/x.html"));
        let already = std::path::Path::new("/tmp/y.html");
        assert_eq!(absolutize(already), already.to_path_buf());
    }

    /// Minimal-but-valid `SearchParams` (passes `validate()`), with a chosen
    /// database path and engine name.
    fn valid_params(database_path: String, engine: &str) -> SearchParams {
        SearchParams {
            enzyme: Enzyme::Trypsin,
            missed_cleavages: 2,
            fixed_modifications: vec![],
            variable_modifications: vec![],
            precursor_tolerance: MassTolerance {
                value: 10.0,
                unit: ToleranceUnit::Ppm,
            },
            fragment_tolerance: MassTolerance {
                value: 0.02,
                unit: ToleranceUnit::Da,
            },
            database_path,
            decoy_strategy: DecoyStrategy::Reverse,
            acquisition_mode: None,
            max_variable_modifications: 3,
            min_peptide_length: 7,
            max_peptide_length: 50,
            engine: Some(engine.to_string()),
        }
    }

    fn dummy_spectrum() -> Spectrum {
        Spectrum {
            scan_number: 1,
            ms_level: MsLevel::MS2,
            retention_time_min: 0.0,
            precursors: vec![],
            mz_array: vec![100.0, 200.0],
            intensity_array: vec![10.0, 20.0],
            representation: SpectrumRepresentation::Centroid,
        }
    }

    /// FIX 1 (DIA path): an invalid engine name must NOT consume cached DIA
    /// spectra. The entry must survive the failed run so the user need not
    /// re-run extract_dia_precursors.
    #[tokio::test]
    async fn dia_invalid_engine_preserves_cached_spectra() {
        let server = ProteinCopilotServer::new();
        let dia_uuid = Uuid::new_v4();

        {
            let mut cache = server.dia_cache.lock().unwrap_or_else(|e| e.into_inner());
            cache.insert(dia_uuid, vec![dummy_spectrum()]);
        }

        let input = RunSearchInput {
            params: Some(valid_params(fixture("small_test.fasta"), "BogusEngine")),
            input_files: vec![],
            database_path: None,
            hints: None,
            dia_run_id: Some(dia_uuid.to_string()),
        };

        let result = server.run_search(Parameters(input)).await;
        let err = match result {
            Ok(_) => panic!("bogus engine must fail"),
            Err(e) => e,
        };
        assert!(
            err.message.contains("not registered"),
            "expected engine-validation error, got: {}",
            err.message
        );

        let cache = server.dia_cache.lock().unwrap_or_else(|e| e.into_inner());
        assert_eq!(
            cache.entries.get(&dia_uuid).map(|s| s.len()),
            Some(1),
            "DIA spectra must NOT be consumed when engine validation fails"
        );
    }

    /// FIX 1 (file path): an invalid engine name must NOT insert an orphaned
    /// "Running" entry into the run cache.
    #[tokio::test]
    async fn file_invalid_engine_leaves_no_orphan_run() {
        let server = ProteinCopilotServer::new();

        let input = RunSearchInput {
            params: Some(valid_params(fixture("small_test.fasta"), "BogusEngine")),
            input_files: vec![fixture("small_test.mgf")],
            database_path: None,
            hints: None,
            dia_run_id: None,
        };

        let result = server.run_search(Parameters(input)).await;
        let err = match result {
            Ok(_) => panic!("bogus engine must fail"),
            Err(e) => e,
        };
        assert!(
            err.message.contains("not registered"),
            "expected engine-validation error, got: {}",
            err.message
        );

        let cache = server.run_cache.lock().unwrap_or_else(|e| e.into_inner());
        assert_eq!(
            cache.len(),
            0,
            "no orphaned Running run must remain when engine validation fails"
        );
    }

    /// FIX 2: a poisoned run cache must self-heal instead of bricking the tool.
    #[test]
    fn poisoned_run_cache_self_heals_in_get_result() {
        let server = ProteinCopilotServer::new();

        let cache = Arc::clone(&server.run_cache);
        let _ = std::thread::spawn(move || {
            let _guard = cache.lock().expect("acquire lock to poison");
            panic!("intentionally poison the run cache mutex");
        })
        .join();

        let missing = Uuid::new_v4().to_string();
        let err = server
            .get_result(&None, &Some(missing))
            .expect_err("missing run_id should error");
        assert_eq!(
            err.code,
            ErrorCode::INVALID_PARAMS,
            "poisoned lock must self-heal, not return INTERNAL_ERROR"
        );
        assert!(err.message.contains("not found"));
    }

    /// Build a minimal valid pFind TSV (16 columns) from `(file_name, scan)`
    /// rows for import-path tests.
    fn write_pfind_tsv(dir: &std::path::Path, rows: &[(&str, u32)]) -> String {
        let path = dir.join("pfind.tsv");
        let mut s = String::from(
            "FileName\tPeptideSequence\tModifications\tPepMass\tPredRT\tCleavageType\tProNCTerm\tProteins\tMH+\tCharge\tScanNo\tRawScore\tDeltaMassPPM\tDeltaRT(Min)\tFinalScore\tQValue\n",
        );
        for (raw, scan) in rows {
            s.push_str(&format!(
                "{raw}\tPEPTIDEK\t\t900.4\t15.1\t3\t0\tsp|P12345|TEST/\t901.4\t2\t{scan}\t10.0\t1.0\t0.0\t0.001\t0\n"
            ));
        }
        std::fs::write(&path, s).expect("write pfind tsv");
        path.to_string_lossy().to_string()
    }

    /// ISSUE 3: pFind import bypasses RT scan matching, so a traversal raw_name
    /// must be rejected before it is joined onto mzml_dir.
    #[test]
    fn pfind_import_rejects_traversal_raw_name() {
        let server = ProteinCopilotServer::new();
        let dir = tempfile::TempDir::new().unwrap();
        let tsv = write_pfind_tsv(dir.path(), &[("../evil", 100)]);

        let input = ImportSearchResultsInput {
            result_file: tsv,
            format: "pfind_tsv".to_string(),
            mzml_dir: dir.path().to_string_lossy().to_string(),
            unimod_path: None,
            rt_tolerance_min: 0.5,
            filter_qvalue: 1.0,
            run_filter: None,
        };

        let err = match server.import_search_results(Parameters(input)) {
            Ok(_) => panic!("traversal raw_name must be rejected"),
            Err(e) => e,
        };
        assert!(
            err.message.contains("invalid raw file name") || err.message.contains("evil"),
            "expected raw-name rejection, got: {}",
            err.message
        );
    }

    /// ISSUE 1: input_files derived from pFind raw names must be deterministically
    /// sorted, not in HashSet order, so annotate_spectrum/extract_xic pick a
    /// stable first file across imports.
    #[test]
    fn pfind_import_input_files_sorted() {
        let server = ProteinCopilotServer::new();
        let dir = tempfile::TempDir::new().unwrap();
        let names: Vec<String> = (0..8).rev().map(|i| format!("run_{i:02}")).collect();
        let rows: Vec<(&str, u32)> = names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.as_str(), (i + 1) as u32))
            .collect();
        let tsv = write_pfind_tsv(dir.path(), &rows);

        let input = ImportSearchResultsInput {
            result_file: tsv,
            format: "pfind_tsv".to_string(),
            mzml_dir: dir.path().to_string_lossy().to_string(),
            unimod_path: None,
            rt_tolerance_min: 0.5,
            filter_qvalue: 1.0,
            run_filter: None,
        };

        let run_id = server
            .import_search_results(Parameters(input))
            .expect("import should succeed")
            .0
            .run_id;
        let result = server.get_result(&None, &Some(run_id)).expect("cached run");

        let mut expected = result.metadata.input_files.clone();
        expected.sort();
        assert_eq!(
            result.metadata.input_files, expected,
            "input_files must be sorted deterministically"
        );
    }
}
