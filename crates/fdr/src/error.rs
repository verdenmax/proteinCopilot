//! Error types for FDR calculation.

use thiserror::Error;

/// Errors from FDR operations.
#[derive(Debug, Error)]
pub enum FdrError {
    /// No PSMs provided for FDR calculation.
    #[error("no PSMs provided for FDR calculation")]
    NoPsms,

    /// No decoy PSMs found — FDR cannot be estimated.
    #[error("no decoy PSMs found; target-decoy FDR requires decoy hits")]
    NoDecoyHits,

    /// Invalid score values.
    #[error("non-finite score value encountered")]
    InvalidScore,
}

impl From<FdrError> for protein_copilot_core::error::CoreError {
    fn from(err: FdrError) -> Self {
        protein_copilot_core::error::CoreError::ValidationError {
            context: "FDR calculation".to_string(),
            detail: err.to_string(),
            suggestion: match &err {
                FdrError::NoPsms => "Ensure search produced PSM results".to_string(),
                FdrError::NoDecoyHits => {
                    "Check decoy strategy; try Reverse if using Shuffle".to_string()
                }
                FdrError::InvalidScore => "Check search engine scoring".to_string(),
            },
        }
    }
}
