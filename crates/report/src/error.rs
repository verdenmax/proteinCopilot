//! Error types for the report module.

use std::path::PathBuf;

use thiserror::Error;

/// Errors that can occur during report generation or export.
#[derive(Debug, Error)]
pub enum ReportError {
    /// I/O error writing output file.
    #[error("I/O error writing {path:?}: {detail}")]
    IoError {
        /// Output file path.
        path: PathBuf,
        /// What went wrong.
        detail: String,
    },

    /// Serialization error.
    #[error("serialization error: {0}")]
    SerializationError(String),

    /// The search result is empty (no PSMs).
    #[error("search result contains no PSMs")]
    EmptyResult,

    /// No MS2 scans were found in the RT window around the target scan.
    #[error("no MS2 scans found in the RT window around scan {scan}")]
    EmptyMs2Window {
        /// Target scan number.
        scan: u32,
    },

    /// Per-scan annotation failed.
    #[error("annotation failed for scan {scan}: {detail}")]
    AnnotationError {
        /// Scan number that failed.
        scan: u32,
        /// Underlying error detail.
        detail: String,
    },
}

impl From<ReportError> for protein_copilot_core::error::CoreError {
    fn from(err: ReportError) -> Self {
        protein_copilot_core::error::CoreError::ValidationError {
            context: "report".to_string(),
            detail: err.to_string(),
            suggestion: "Check search result data".to_string(),
        }
    }
}
