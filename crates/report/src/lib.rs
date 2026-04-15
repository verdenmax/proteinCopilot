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
pub mod visualize;
pub mod xic_visualize;

pub mod unified_types;
pub mod unified_visualize;

pub use error::ReportError;

/// Escapes JSON for safe embedding inside HTML `<script>` tags.
///
/// Replaces `<` with `\u003c` and `>` with `\u003e` to prevent:
/// - Premature `</script>` tag closure
/// - HTML injection via `<` and `>` in JSON string values
///
/// The escaped string remains valid JSON (parseable by `JSON.parse()`).
pub fn escape_json_for_html(json: &str) -> String {
    json.replace('<', r"\u003c").replace('>', r"\u003e")
}

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

    /// Renders a spectrum annotation as a self-contained HTML file.
    pub fn render_annotation(
        annotation: &protein_copilot_search_engine::annotate::SpectrumAnnotation,
        output_path: &Path,
    ) -> Result<(), ReportError> {
        visualize::render_annotation_html(annotation, output_path)
    }

    /// Renders XIC data as a self-contained HTML file with Plotly.js charts.
    pub fn render_xic(
        xic_data: &protein_copilot_xic::XicData,
        output_path: &Path,
        plotly_mode: protein_copilot_xic::PlotlyMode,
    ) -> Result<(), ReportError> {
        xic_visualize::render_xic_html(xic_data, output_path, plotly_mode)
    }

    /// Renders unified annotation + XIC as a self-contained HTML file.
    pub fn render_unified(
        data: &crate::unified_types::UnifiedViewData,
        output_path: &Path,
        plotly_mode: protein_copilot_xic::PlotlyMode,
    ) -> Result<(), ReportError> {
        unified_visualize::render_unified_html(data, output_path, plotly_mode)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn escape_script_close_tag() {
        let input = r#"{"val":"</script><b>"}"#;
        let escaped = escape_json_for_html(input);
        assert!(!escaped.contains("</script>"), "must not contain literal </script>");
        assert!(!escaped.contains('<'), "must not contain literal <");
        assert!(!escaped.contains('>'), "must not contain literal >");
    }

    #[test]
    fn escape_angle_brackets() {
        let input = r#"{"key":"<img onerror=alert(1)>"}"#;
        let escaped = escape_json_for_html(input);
        assert!(!escaped.contains('<'));
        assert!(!escaped.contains('>'));
        assert!(escaped.contains(r"\u003c"));
        assert!(escaped.contains(r"\u003e"));
    }

    #[test]
    fn escape_preserves_normal_json() {
        let input = r#"{"score":0.95,"peptide":"ACDK"}"#;
        let escaped = escape_json_for_html(input);
        assert_eq!(input, escaped, "no angle brackets = no changes");
    }
}
