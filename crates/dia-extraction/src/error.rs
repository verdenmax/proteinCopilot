//! Error types for DIA precursor extraction.

use thiserror::Error;

/// Errors that can occur during DIA precursor extraction.
#[derive(Debug, Error)]
pub enum DiaExtractionError {
    /// DIA data detected but no MS1 spectra found for precursor extraction.
    #[error("DIA data detected but no MS1 spectra found for precursor extraction")]
    NoMs1Spectra,

    /// No MS2 spectra found in input data.
    #[error("No MS2 spectra found in input data")]
    NoMs2Spectra,

    /// Invalid isolation window configuration.
    #[error("Invalid isolation window: {detail}")]
    InvalidIsolationWindow {
        /// Description of the isolation window problem.
        detail: String,
    },

    /// General extraction failure.
    #[error("Precursor extraction failed: {detail}")]
    ExtractionFailed {
        /// Description of the extraction failure.
        detail: String,
    },
}
