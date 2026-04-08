//! Error types for XIC extraction.

use std::path::PathBuf;

use thiserror::Error;

/// Errors that can occur during XIC extraction.
#[derive(Debug, Error)]
pub enum XicError {
    /// Input file format is not supported for XIC extraction.
    #[error("XIC extraction requires mzML format, got: {path}")]
    UnsupportedFormat {
        /// The file path.
        path: PathBuf,
    },

    /// Target scan not found in the spectrum file.
    #[error("target scan {scan} not found in {path}")]
    ScanNotFound {
        /// The file path.
        path: PathBuf,
        /// The requested scan number.
        scan: u32,
    },

    /// Target scan has no isolation window (required for DIA XIC).
    #[error("scan {scan} has no isolation window — XIC requires MS2 with isolation window")]
    NoIsolationWindow {
        /// The scan number.
        scan: u32,
    },

    /// No matching MS2 scans found in the same isolation window.
    #[error("no MS2 scans found in same isolation window as scan {scan}")]
    NoMatchingCycles {
        /// The target scan number.
        scan: u32,
    },

    /// Peptide sequence is empty or contains invalid residues.
    #[error("invalid peptide sequence: {detail}")]
    InvalidPeptide {
        /// What went wrong.
        detail: String,
    },

    /// Spectrum I/O error.
    #[error("spectrum I/O error: {0}")]
    SpectrumIo(#[from] protein_copilot_spectrum_io::SpectrumIoError),
}
