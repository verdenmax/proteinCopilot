//! MCP Tool definitions — thin wrappers around library crate functions.

use std::path::{Path, PathBuf};

use rmcp::handler::server::tool::ToolRouter;
use rmcp::handler::server::wrapper::{Json, Parameters};
use rmcp::model::ServerInfo;
use rmcp::schemars;
use rmcp::tool;
use rmcp::tool_router;
use rmcp::ServerHandler;
use serde::Deserialize;

use protein_copilot_core::ai_decision::AiDecision;
use protein_copilot_core::engine::{HealthStatus, SearchEngineAdapter};
use protein_copilot_core::search_params::SearchParams;
use protein_copilot_core::search_result::{SearchResult, SearchResultSummary};
use protein_copilot_core::spectrum::SpectrumSummary;
use protein_copilot_param_recommend::{ParamRecommender, UserHints};
use protein_copilot_report::ReportGenerator;
use protein_copilot_search_engine::SimpleSearchEngine;
use protein_copilot_spectrum_io::{create_reader, detect_format};

// ---------------------------------------------------------------------------
// Tool input types
// ---------------------------------------------------------------------------

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ReadSpectraInput {
    /// Path to the spectrum file (.mgf or .mzML)
    file_path: String,
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

#[tool_router]
impl ProteinCopilotServer {
    /// Read a mass spectrometry file and return a statistical summary.
    /// Use this as the first step to understand the input data characteristics.
    /// Supports .mgf and .mzML formats (DDA and DIA).
    #[tool(
        name = "read_spectra",
        description = "Read a mass spectrometry file (mgf/mzML) and return a statistical summary including spectrum count, m/z range, RT range, charge distribution, and median peaks per spectrum."
    )]
    fn read_spectra(
        &self,
        Parameters(input): Parameters<ReadSpectraInput>,
    ) -> Json<SpectrumSummary> {
        let path = Path::new(&input.file_path);
        let info = detect_format(path)
            .map_err(|e| {
                rmcp::ErrorData::new(rmcp::model::ErrorCode::INVALID_PARAMS, e.to_string(), None)
            })
            .unwrap();
        let reader = create_reader(&info);
        let summary = reader
            .read_summary(path)
            .map_err(|e| {
                rmcp::ErrorData::new(rmcp::model::ErrorCode::INTERNAL_ERROR, e.to_string(), None)
            })
            .unwrap();
        Json(summary)
    }

    /// Recommend search parameters based on spectrum characteristics.
    /// Call read_spectra first to obtain the SpectrumSummary, then pass it here.
    #[tool(
        name = "recommend_params",
        description = "Recommend search parameters based on spectrum file characteristics. Input: SpectrumSummary from read_spectra + optional UserHints. Output: recommended SearchParams with confidence score and explanation."
    )]
    fn recommend_params(
        &self,
        Parameters(input): Parameters<RecommendParamsInput>,
    ) -> Json<AiDecision<SearchParams>> {
        let recommender = ParamRecommender;
        let result = recommender
            .recommend(&input.summary, input.hints.as_ref())
            .map_err(|e| {
                rmcp::ErrorData::new(rmcp::model::ErrorCode::INTERNAL_ERROR, e.to_string(), None)
            })
            .unwrap();
        Json(result)
    }

    /// List all available search parameter presets.
    #[tool(
        name = "list_presets",
        description = "List all built-in search parameter presets (standard, phospho, TMT, SILAC, open search)."
    )]
    fn list_presets(&self) -> String {
        let presets = ParamRecommender::list_presets();
        serde_json::to_string_pretty(&presets).unwrap_or_else(|e| format!("Error: {e}"))
    }

    /// Execute a database search against spectrum files.
    /// Requires SearchParams (from recommend_params or manual) and spectrum file paths.
    #[tool(
        name = "run_search",
        description = "Run a proteomics database search. Input: SearchParams + spectrum file paths. Output: SearchResult with PSMs, peptides, proteins, and summary statistics."
    )]
    async fn run_search(
        &self,
        Parameters(input): Parameters<RunSearchInput>,
    ) -> Json<SearchResult> {
        let engine = SimpleSearchEngine::new();
        let files: Vec<PathBuf> = input.input_files.iter().map(PathBuf::from).collect();
        let result = engine
            .search(&input.params, &files)
            .await
            .map_err(|e| {
                rmcp::ErrorData::new(rmcp::model::ErrorCode::INTERNAL_ERROR, e.to_string(), None)
            })
            .unwrap();
        Json(result)
    }

    /// Check available search engines and their health status.
    #[tool(
        name = "check_engine",
        description = "Check available search engines and their health status. Returns engine name, version, and availability."
    )]
    async fn check_engine(&self) -> String {
        let engine = SimpleSearchEngine::new();
        let info = engine.engine_info();
        let status = engine
            .health_check()
            .await
            .unwrap_or(HealthStatus::Unavailable {
                reason: "health check failed".to_string(),
            });
        serde_json::to_string_pretty(&serde_json::json!({
            "engine": info,
            "status": status,
        }))
        .unwrap_or_else(|e| format!("Error: {e}"))
    }

    /// Generate a statistical summary of search results with FDR filtering.
    #[tool(
        name = "generate_summary",
        description = "Generate a statistical summary from search results with 1% FDR filtering. Includes identification rate, median score, modification/charge distributions."
    )]
    fn generate_summary(
        &self,
        Parameters(input): Parameters<GenerateSummaryInput>,
    ) -> Json<SearchResultSummary> {
        let summary = ReportGenerator::generate_summary(&input.result);
        Json(summary)
    }

    /// Export search results as TSV and JSON files.
    #[tool(
        name = "export_results",
        description = "Export search results to files. Creates psm.tsv, peptide.tsv, protein.tsv, result.json, and run_metadata.json in the specified output directory."
    )]
    fn export_results(&self, Parameters(input): Parameters<ExportResultsInput>) -> String {
        let output_dir = Path::new(&input.output_dir);
        let mut exported = Vec::new();

        if let Err(e) = ReportGenerator::export_tsv(&input.result, output_dir) {
            return format!("TSV export error: {e}");
        }
        exported.push("psm.tsv, peptide.tsv, protein.tsv");

        if let Err(e) = ReportGenerator::export_json(&input.result, &output_dir.join("result.json"))
        {
            return format!("JSON export error: {e}");
        }
        exported.push("result.json");

        if let Err(e) = ReportGenerator::export_metadata(
            &input.result.metadata,
            &output_dir.join("run_metadata.json"),
        ) {
            return format!("Metadata export error: {e}");
        }
        exported.push("run_metadata.json");

        format!(
            "Exported to {}: {}",
            output_dir.display(),
            exported.join(", ")
        )
    }
}
