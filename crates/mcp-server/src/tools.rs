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
use protein_copilot_core::spectrum::{Spectrum, SpectrumSummary};
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

// ---------------------------------------------------------------------------
// Server struct
// ---------------------------------------------------------------------------

#[allow(dead_code)]
pub struct ProteinCopilotServer {
    tool_router: ToolRouter<Self>,
    registry: protein_copilot_search_engine::EngineRegistry,
    /// Server-side cache of SearchResults keyed by run_id.
    result_cache: Arc<Mutex<HashMap<Uuid, SearchResult>>>,
    /// Progress tracking for async searches.
    progress_cache: Arc<Mutex<HashMap<Uuid, SearchProgress>>>,
}

impl ProteinCopilotServer {
    pub fn new() -> Self {
        let mut registry = protein_copilot_search_engine::EngineRegistry::new();
        registry.register(Box::new(SimpleSearchEngine::new()));
        Self {
            tool_router: Self::tool_router(),
            registry,
            result_cache: Arc::new(Mutex::new(HashMap::new())),
            progress_cache: Arc::new(Mutex::new(HashMap::new())),
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
                .result_cache
                .lock()
                .map_err(|_| mcp_err(ErrorCode::INTERNAL_ERROR, "cache lock failed"))?;
            return cache.get(&id).cloned().ok_or_else(|| {
                mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    format!("run_id {id} not found in cache"),
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

        if let Ok(mut progress) = self.progress_cache.lock() {
            progress.insert(
                run_id,
                SearchProgress {
                    run_id,
                    status: "Running".to_string(),
                    progress_pct: Some(0.0),
                    elapsed_sec: 0.0,
                    estimated_remaining_sec: None,
                },
            );
        }

        let result_cache = Arc::clone(&self.result_cache);
        let progress_cache = Arc::clone(&self.progress_cache);
        let engine = SimpleSearchEngine::new();

        tokio::spawn(async move {
            let start = Instant::now();
            match engine.search(&params, &files).await {
                Ok(mut result) => {
                    result.run_id = run_id;
                    result.metadata.run_id = run_id;
                    let duration = start.elapsed().as_secs_f64();
                    if let Ok(mut cache) = result_cache.lock() {
                        cache.insert(run_id, result);
                    }
                    if let Ok(mut progress) = progress_cache.lock() {
                        progress.insert(
                            run_id,
                            SearchProgress {
                                run_id,
                                status: "Completed".to_string(),
                                progress_pct: Some(1.0),
                                elapsed_sec: duration,
                                estimated_remaining_sec: Some(0.0),
                            },
                        );
                    }
                }
                Err(e) => {
                    let duration = start.elapsed().as_secs_f64();
                    if let Ok(mut progress) = progress_cache.lock() {
                        progress.insert(
                            run_id,
                            SearchProgress {
                                run_id,
                                status: format!("Failed: {e}"),
                                progress_pct: None,
                                elapsed_sec: duration,
                                estimated_remaining_sec: None,
                            },
                        );
                    }
                }
            }
        });

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
        let progress = self
            .progress_cache
            .lock()
            .map_err(|_| mcp_err(ErrorCode::INTERNAL_ERROR, "cache lock failed"))?;
        let p = progress
            .get(&id)
            .cloned()
            .ok_or_else(|| mcp_err(ErrorCode::INVALID_PARAMS, format!("run_id {id} not found")))?;
        Ok(Json(p))
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
}
