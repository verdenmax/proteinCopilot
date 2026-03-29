//! Error types for the search engine module.

use std::path::PathBuf;

use thiserror::Error;

/// Errors that can occur during search engine operations.
#[derive(Debug, Error)]
pub enum SearchEngineError {
    /// Search parameters failed validation.
    #[error("invalid search parameters: {detail}")]
    InvalidParams {
        /// What is wrong.
        detail: String,
    },

    /// FASTA database file could not be read.
    #[error("FASTA error ({path:?}): {detail}")]
    FastaError {
        /// Path to the FASTA file.
        path: PathBuf,
        /// What went wrong.
        detail: String,
    },

    /// I/O error during search execution.
    #[error("I/O error: {detail}")]
    IoError {
        /// What went wrong.
        detail: String,
    },

    /// The requested engine was not found in the registry.
    #[error("engine '{name}' not found in registry")]
    EngineNotFound {
        /// Engine name that was requested.
        name: String,
    },

    /// Search execution failed.
    #[error("search execution failed: {detail}")]
    ExecutionError {
        /// What went wrong.
        detail: String,
    },

    /// No spectra provided for search.
    #[error("no input spectra provided")]
    NoInputSpectra,
}

impl From<SearchEngineError> for protein_copilot_core::error::CoreError {
    fn from(err: SearchEngineError) -> Self {
        protein_copilot_core::error::CoreError::SearchEngineError {
            engine: "search-engine".to_string(),
            detail: err.to_string(),
            suggestion: match &err {
                SearchEngineError::InvalidParams { .. } => {
                    "Review and correct the search parameters".to_string()
                }
                SearchEngineError::FastaError { .. } => {
                    "Check the FASTA database file path and format".to_string()
                }
                SearchEngineError::IoError { .. } => {
                    "Check file permissions and disk availability".to_string()
                }
                SearchEngineError::EngineNotFound { .. } => {
                    "Use list_engines to see available engines".to_string()
                }
                SearchEngineError::ExecutionError { .. } => {
                    "Check engine logs for details".to_string()
                }
                SearchEngineError::NoInputSpectra => {
                    "Provide at least one spectrum file".to_string()
                }
            },
        }
    }
}
