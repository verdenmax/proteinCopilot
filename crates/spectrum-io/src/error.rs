//! Error types for spectrum file I/O operations.

use std::path::PathBuf;

use thiserror::Error;

/// Errors that can occur during spectrum file parsing.
#[derive(Debug, Error)]
pub enum SpectrumIoError {
    /// File not found at the specified path.
    #[error("file not found: {path}")]
    FileNotFound {
        /// The missing file path.
        path: PathBuf,
    },

    /// File format could not be determined.
    #[error("unable to detect format for: {path}")]
    UnknownFormat {
        /// The file path.
        path: PathBuf,
    },

    /// Unsupported file format.
    #[error("unsupported format '{format}' for: {path}")]
    UnsupportedFormat {
        /// The detected format name.
        format: String,
        /// The file path.
        path: PathBuf,
    },

    /// I/O error while reading the file.
    #[error("I/O error reading {path}: {source}")]
    IoError {
        /// The file path being read.
        path: PathBuf,
        /// The underlying I/O error.
        source: std::io::Error,
    },

    /// Parse error in spectrum file content.
    #[error("parse error in {path} at line {line}: {detail}")]
    ParseError {
        /// The file path.
        path: PathBuf,
        /// Approximate line number (0 if unknown).
        line: usize,
        /// What went wrong.
        detail: String,
    },

    /// XML parsing error (mzML-specific).
    #[error("XML error in {path}: {detail}")]
    XmlError {
        /// The file path.
        path: PathBuf,
        /// What went wrong.
        detail: String,
    },

    /// Binary data decoding error (mzML-specific).
    #[error("binary decode error in {path}: {detail}")]
    BinaryDecodeError {
        /// The file path.
        path: PathBuf,
        /// What went wrong.
        detail: String,
    },

    /// Spectrum validation failed after parsing.
    #[error("spectrum validation failed (scan {scan}): {detail}")]
    ValidationError {
        /// Scan number of the invalid spectrum.
        scan: u32,
        /// What validation failed.
        detail: String,
    },

    /// Requested scan number was not found in the file.
    #[error("scan {scan} not found in {path}")]
    ScanNotFound {
        /// The file path.
        path: PathBuf,
        /// The requested scan number.
        scan: u32,
    },
}

impl From<std::io::Error> for SpectrumIoError {
    fn from(err: std::io::Error) -> Self {
        SpectrumIoError::IoError {
            path: PathBuf::from("<unknown>"),
            source: err,
        }
    }
}

impl From<SpectrumIoError> for protein_copilot_core::error::CoreError {
    fn from(err: SpectrumIoError) -> Self {
        match err {
            SpectrumIoError::FileNotFound { path } => {
                protein_copilot_core::error::CoreError::FileNotFound { path }
            }
            SpectrumIoError::UnknownFormat { path } => {
                protein_copilot_core::error::CoreError::UnsupportedFormat {
                    format: path
                        .extension()
                        .map(|e| e.to_string_lossy().to_string())
                        .unwrap_or_default(),
                    supported: vec!["mzML".to_string(), "mgf".to_string()],
                }
            }
            SpectrumIoError::UnsupportedFormat { format, .. } => {
                protein_copilot_core::error::CoreError::UnsupportedFormat {
                    format,
                    supported: vec!["mzML".to_string(), "mgf".to_string()],
                }
            }
            SpectrumIoError::IoError { path, source } => {
                protein_copilot_core::error::CoreError::SpectrumParseError {
                    format: "unknown".to_string(),
                    detail: format!("I/O error reading {}: {source}", path.display()),
                    suggestion: "Check file permissions and disk availability".to_string(),
                }
            }
            SpectrumIoError::ParseError { path, line, detail } => {
                protein_copilot_core::error::CoreError::SpectrumParseError {
                    format: path
                        .extension()
                        .map(|e| e.to_string_lossy().to_string())
                        .unwrap_or_default(),
                    detail: format!("line {line}: {detail}"),
                    suggestion: "Check spectrum file integrity".to_string(),
                }
            }
            SpectrumIoError::XmlError { path, detail } => {
                protein_copilot_core::error::CoreError::SpectrumParseError {
                    format: "mzML".to_string(),
                    detail: format!("{}: {detail}", path.display()),
                    suggestion: "Verify the mzML file is well-formed XML".to_string(),
                }
            }
            SpectrumIoError::BinaryDecodeError { path, detail } => {
                protein_copilot_core::error::CoreError::SpectrumParseError {
                    format: "mzML".to_string(),
                    detail: format!("{}: {detail}", path.display()),
                    suggestion: "Check binary data encoding in the mzML file".to_string(),
                }
            }
            SpectrumIoError::ValidationError { scan, detail } => {
                protein_copilot_core::error::CoreError::SpectrumParseError {
                    format: "spectrum".to_string(),
                    detail: format!("scan {scan}: {detail}"),
                    suggestion: "Check spectrum data integrity".to_string(),
                }
            }
            SpectrumIoError::ScanNotFound { path, scan } => {
                protein_copilot_core::error::CoreError::SpectrumParseError {
                    format: path
                        .extension()
                        .map(|e| e.to_string_lossy().to_string())
                        .unwrap_or_default(),
                    detail: format!("scan {scan} not found in {}", path.display()),
                    suggestion: "Verify the scan number exists in the file".to_string(),
                }
            }
        }
    }
}
