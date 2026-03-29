//! MCP Tool definitions — thin wrappers around library crate functions.
//!
//! Each tool is a `Result<Json<T>, ErrorData>` returning function that:
//! 1. Parses parameters
//! 2. Delegates to a library crate
//! 3. Returns structured JSON or a proper MCP error

use std::path::{Path, PathBuf};

use rmcp::handler::server::tool::ToolRouter;
use rmcp::handler::server::wrapper::{Json, Parameters};
use rmcp::model::{ErrorCode, ServerInfo};
use rmcp::schemars;
use rmcp::{ErrorData, ServerHandler};
use serde::{Deserialize, Serialize};

use protein_copilot_core::ai_decision::AiDecision;
use protein_copilot_core::engine::{HealthStatus, SearchEngineAdapter};
use protein_copilot_core::search_params::SearchParams;
use protein_copilot_core::search_result::{SearchResult, SearchResultSummary};
use protein_copilot_core::spectrum::{Spectrum, SpectrumSummary};
use protein_copilot_param_recommend::{ParamRecommender, SearchPreset, UserHints};
use protein_copilot_report::ReportGenerator;
use protein_copilot_search_engine::SimpleSearchEngine;

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
    /// Spectrum summary (from read_spectra)
    summary: SpectrumSummary,
    /// Optional user hints
    hints: Option<UserHints>,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct RunSearchInput {
    /// Search parameters
    params: SearchParams,
    /// Paths to spectrum files
    input_files: Vec<String>,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct GenerateSummaryInput {
    /// Search result to summarize
    result: SearchResult,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ExportResultsInput {
    /// Search result to export
    result: SearchResult,
    /// Output directory path
    output_dir: String,
}

/// Engine status response
#[derive(Debug, Serialize, schemars::JsonSchema)]
struct EngineStatus {
    engine: protein_copilot_core::engine::EngineInfo,
    status: HealthStatus,
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
}

impl ProteinCopilotServer {
    pub fn new() -> Self {
        Self {
            tool_router: Self::tool_router(),
        }
    }
}

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
            .map_err(|e| mcp_err(ErrorCode::INVALID_PARAMS, e))?;
        let reader = protein_copilot_spectrum_io::create_reader(&info);
        let summary = reader
            .read_summary(path)
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))?;
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
            .map_err(|e| mcp_err(ErrorCode::INVALID_PARAMS, e))?;
        let reader = protein_copilot_spectrum_io::create_reader(&info);
        let spectrum = reader
            .read_spectrum(path, input.scan_number)
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))?;
        Ok(Json(spectrum))
    }

    /// Recommend search parameters based on spectrum characteristics.
    #[rmcp::tool(
        name = "recommend_params",
        description = "Recommend search parameters based on spectrum file characteristics. Input: SpectrumSummary from read_spectra + optional UserHints (experiment_type, instrument_type, enzyme). Output: recommended SearchParams with confidence score and explanation."
    )]
    fn recommend_params(
        &self,
        Parameters(input): Parameters<RecommendParamsInput>,
    ) -> Result<Json<AiDecision<SearchParams>>, ErrorData> {
        let recommender = ParamRecommender;
        let result = recommender
            .recommend(&input.summary, input.hints.as_ref())
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))?;
        Ok(Json(result))
    }

    /// List all available search parameter presets.
    #[rmcp::tool(
        name = "list_presets",
        description = "List all built-in search parameter presets (standard, phospho, TMT, SILAC, open search). Each preset includes name, description, parameters, and applicable scenarios."
    )]
    fn list_presets(&self) -> Json<Vec<SearchPreset>> {
        Json(ParamRecommender::list_presets())
    }

    /// Execute a database search against spectrum files.
    #[rmcp::tool(
        name = "run_search",
        description = "Run a proteomics database search. Input: SearchParams (from recommend_params or manual) + spectrum file paths. Output: SearchResult with PSMs, peptides, proteins, and summary statistics. Note: set database_path in params to the FASTA file path."
    )]
    async fn run_search(
        &self,
        Parameters(input): Parameters<RunSearchInput>,
    ) -> Result<Json<SearchResult>, ErrorData> {
        let engine = SimpleSearchEngine::new();
        let files: Vec<PathBuf> = input.input_files.iter().map(PathBuf::from).collect();
        let result = engine
            .search(&input.params, &files)
            .await
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))?;
        Ok(Json(result))
    }

    /// Check available search engines and their health status.
    #[rmcp::tool(
        name = "check_engine",
        description = "Check available search engines and their health status. Returns engine name, version, supported features, and availability."
    )]
    async fn check_engine(&self) -> Json<EngineStatus> {
        let engine = SimpleSearchEngine::new();
        let info = engine.engine_info();
        let status = engine
            .health_check()
            .await
            .unwrap_or(HealthStatus::Unavailable {
                reason: "health check failed".to_string(),
            });
        Json(EngineStatus {
            engine: info,
            status,
        })
    }

    /// Generate a statistical summary with FDR filtering.
    #[rmcp::tool(
        name = "generate_summary",
        description = "Generate a statistical summary from search results with 1% FDR filtering. Includes identification rate, median score, median delta ppm, modification and charge distributions. Use this after run_search to interpret results."
    )]
    fn generate_summary(
        &self,
        Parameters(input): Parameters<GenerateSummaryInput>,
    ) -> Json<SearchResultSummary> {
        Json(ReportGenerator::generate_summary(&input.result))
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
        let output_dir = Path::new(&input.output_dir);

        ReportGenerator::export_tsv(&input.result, output_dir)
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))?;
        ReportGenerator::export_json(&input.result, &output_dir.join("result.json"))
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))?;
        ReportGenerator::export_metadata(
            &input.result.metadata,
            &output_dir.join("run_metadata.json"),
        )
        .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))?;

        Ok(format!(
            "Exported to {}: psm.tsv, peptide.tsv, protein.tsv, result.json, run_metadata.json",
            output_dir.display()
        ))
    }
}
