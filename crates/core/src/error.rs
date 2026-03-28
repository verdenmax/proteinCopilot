//! Domain error types for ProteinCopilot.
//!
//! [`CoreError`] is the unified error enum used across MCP tool boundaries.
//! Each variant carries contextual information and a human-readable
//! `suggestion` to help users (and the LLM) understand how to fix the issue.
//!
//! Module-internal errors (`SpectrumError`, `SearchParamsError`, etc.) are
//! converted into `CoreError` at MCP tool boundaries via `From` impls.

use std::path::PathBuf;

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use thiserror::Error;

// ---------------------------------------------------------------------------
// CoreError
// ---------------------------------------------------------------------------

/// Unified error type for ProteinCopilot domain operations.
///
/// Every variant includes enough context for the AI layer to diagnose
/// the problem and suggest corrective actions to the user.
#[derive(Debug, Error)]
pub enum CoreError {
    /// Failed to parse a spectrum file.
    #[error("spectrum parse error ({format}): {detail}")]
    SpectrumParseError {
        /// File format that failed to parse (e.g. "mzML", "mgf").
        format: String,
        /// What went wrong.
        detail: String,
        /// Suggested corrective action.
        suggestion: String,
    },

    /// Search parameters failed validation.
    #[error("invalid search parameter '{field}': {reason}")]
    InvalidSearchParams {
        /// Which parameter field is invalid.
        field: String,
        /// Why it is invalid.
        reason: String,
        /// Suggested corrective action.
        suggestion: String,
    },

    /// Search engine execution failed.
    #[error("search engine error ({engine}): {detail}")]
    SearchEngineError {
        /// Name of the search engine.
        engine: String,
        /// What went wrong.
        detail: String,
        /// Suggested corrective action.
        suggestion: String,
    },

    /// A required file was not found.
    #[error("file not found: {path}")]
    FileNotFound {
        /// The missing file path.
        path: PathBuf,
    },

    /// The spectrum file format is not supported.
    #[error("unsupported format '{format}', supported: {}", supported.join(", "))]
    UnsupportedFormat {
        /// The unsupported format name.
        format: String,
        /// List of supported format names.
        supported: Vec<String>,
    },

    /// SSH connection to a remote search engine host failed.
    #[error("SSH connection to '{host}' failed: {detail}")]
    SshConnectionError {
        /// The remote hostname or address.
        host: String,
        /// What went wrong.
        detail: String,
    },

    /// Failed to parse search engine result output.
    #[error("result parse error ({engine}, {file:?}): {detail}")]
    ResultParseError {
        /// Name of the search engine.
        engine: String,
        /// Path to the result file.
        file: PathBuf,
        /// What went wrong.
        detail: String,
    },
}

impl CoreError {
    /// Returns a user-facing suggestion for how to fix this error.
    ///
    /// The AI layer uses this to provide actionable guidance.
    pub fn suggestion(&self) -> &str {
        match self {
            CoreError::SpectrumParseError { suggestion, .. } => suggestion,
            CoreError::InvalidSearchParams { suggestion, .. } => suggestion,
            CoreError::SearchEngineError { suggestion, .. } => suggestion,
            CoreError::FileNotFound { .. } => {
                "Check that the file path is correct and the file exists"
            }
            CoreError::UnsupportedFormat { .. } => {
                "Convert the file to a supported format (mzML or mgf)"
            }
            CoreError::SshConnectionError { .. } => {
                "Check SSH credentials, network connectivity, and that the remote host is running"
            }
            CoreError::ResultParseError { .. } => {
                "Check the search engine output for corruption or version incompatibility"
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Serializable error summary (for MCP tool responses)
// ---------------------------------------------------------------------------

/// A serializable representation of a [`CoreError`] for MCP tool responses.
///
/// MCP tools return JSON, so errors must be serializable. This struct
/// carries the error category, message, and suggestion as structured data.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct ErrorReport {
    /// Error category (e.g. "SpectrumParseError", "FileNotFound").
    pub category: String,
    /// Human-readable error message.
    pub message: String,
    /// Suggested corrective action.
    pub suggestion: String,
}

impl From<&CoreError> for ErrorReport {
    fn from(err: &CoreError) -> Self {
        let category = match err {
            CoreError::SpectrumParseError { .. } => "SpectrumParseError",
            CoreError::InvalidSearchParams { .. } => "InvalidSearchParams",
            CoreError::SearchEngineError { .. } => "SearchEngineError",
            CoreError::FileNotFound { .. } => "FileNotFound",
            CoreError::UnsupportedFormat { .. } => "UnsupportedFormat",
            CoreError::SshConnectionError { .. } => "SshConnectionError",
            CoreError::ResultParseError { .. } => "ResultParseError",
        };
        ErrorReport {
            category: category.to_string(),
            message: err.to_string(),
            suggestion: err.suggestion().to_string(),
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spectrum_parse_error_display_and_suggestion() {
        let err = CoreError::SpectrumParseError {
            format: "mzML".to_string(),
            detail: "unexpected EOF at byte 1024".to_string(),
            suggestion: "Re-download the file or check for truncation".to_string(),
        };
        assert!(err.to_string().contains("mzML"));
        assert!(err.to_string().contains("unexpected EOF"));
        assert!(err.suggestion().contains("Re-download"));
    }

    #[test]
    fn invalid_search_params_display_and_suggestion() {
        let err = CoreError::InvalidSearchParams {
            field: "precursor_tolerance".to_string(),
            reason: "value 0 is not positive".to_string(),
            suggestion: "Set precursor_tolerance to a positive value (e.g. 20 ppm)".to_string(),
        };
        assert!(err.to_string().contains("precursor_tolerance"));
        assert!(err.suggestion().contains("20 ppm"));
    }

    #[test]
    fn file_not_found_has_default_suggestion() {
        let err = CoreError::FileNotFound {
            path: PathBuf::from("/data/missing.fasta"),
        };
        assert!(err.to_string().contains("missing.fasta"));
        assert!(err.suggestion().contains("file path"));
    }

    #[test]
    fn unsupported_format_lists_supported() {
        let err = CoreError::UnsupportedFormat {
            format: "raw".to_string(),
            supported: vec!["mzML".to_string(), "mgf".to_string()],
        };
        assert!(err.to_string().contains("raw"));
        assert!(err.to_string().contains("mzML"));
        assert!(err.to_string().contains("mgf"));
    }

    #[test]
    fn ssh_connection_error_display() {
        let err = CoreError::SshConnectionError {
            host: "compute-01.lab.edu".to_string(),
            detail: "connection refused".to_string(),
        };
        assert!(err.to_string().contains("compute-01"));
        assert!(err.suggestion().contains("SSH"));
    }

    #[test]
    fn result_parse_error_display() {
        let err = CoreError::ResultParseError {
            engine: "pFind".to_string(),
            file: PathBuf::from("/results/output.spectra"),
            detail: "missing header line".to_string(),
        };
        assert!(err.to_string().contains("pFind"));
        assert!(err.to_string().contains("output.spectra"));
    }

    #[test]
    fn error_report_from_core_error() {
        let err = CoreError::SearchEngineError {
            engine: "pFind".to_string(),
            detail: "process exited with code 1".to_string(),
            suggestion: "Check pFind log files for details".to_string(),
        };
        let report = ErrorReport::from(&err);
        assert_eq!(report.category, "SearchEngineError");
        assert!(report.message.contains("pFind"));
        assert!(report.suggestion.contains("log files"));
    }

    #[test]
    fn error_report_serde_roundtrip() {
        let report = ErrorReport {
            category: "FileNotFound".to_string(),
            message: "file not found: /data/test.fasta".to_string(),
            suggestion: "Check the file path".to_string(),
        };
        let json = serde_json::to_string_pretty(&report).unwrap();
        let back: ErrorReport = serde_json::from_str(&json).unwrap();
        assert_eq!(report.category, back.category);
        assert_eq!(report.message, back.message);
        assert_eq!(report.suggestion, back.suggestion);
    }
}
