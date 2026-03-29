//! # ProteinCopilot Report Generator
//!
//! Converts [`SearchResult`] into various output formats:
//! - Statistical summary ([`SearchResultSummary`]) for LLM interpretation
//! - TSV files (PSM, peptide, protein level) for researcher analysis
//! - JSON export for programmatic consumption
//! - Run metadata for auditability
//!
//! This is a pure computation crate — no MCP, network, or LLM dependencies.

pub mod error;
pub mod export;
pub mod summary;

pub use error::ReportError;

use std::path::Path;

use protein_copilot_core::run_metadata::RunMetadata;
use protein_copilot_core::search_result::{SearchResult, SearchResultSummary};

/// Report generator — stateless utility for result formatting and export.
pub struct ReportGenerator;

impl ReportGenerator {
    /// Generates a statistical summary from search results.
    ///
    /// Applies 1% FDR filtering (if q-values are available) and computes
    /// identification rate, median scores, modification/charge distributions.
    pub fn generate_summary(result: &SearchResult) -> SearchResultSummary {
        summary::generate_summary(result)
    }

    /// Exports search results as TSV files (psm.tsv, peptide.tsv, protein.tsv).
    pub fn export_tsv(result: &SearchResult, output_dir: &Path) -> Result<(), ReportError> {
        export::export_tsv(result, output_dir)
    }

    /// Exports the complete search result as JSON.
    pub fn export_json(result: &SearchResult, output_path: &Path) -> Result<(), ReportError> {
        export::export_json(result, output_path)
    }

    /// Exports run metadata as JSON.
    pub fn export_metadata(metadata: &RunMetadata, output_path: &Path) -> Result<(), ReportError> {
        export::export_metadata(metadata, output_path)
    }
}
