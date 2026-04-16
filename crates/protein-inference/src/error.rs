//! Error types for protein inference.

use thiserror::Error;

/// Errors that can occur during protein inference.
#[derive(Debug, Error)]
pub enum InferenceError {
    /// No PSMs provided for inference.
    #[error("no PSMs provided for protein inference")]
    NoPsms,

    /// No target proteins found in PSMs.
    #[error("no target proteins found in PSMs")]
    NoTargetProteins,

    /// FDR calculation failed.
    #[error("FDR calculation failed: {0}")]
    FdrFailed(String),

    /// FASTA file required but not provided.
    #[error("FASTA file required for coverage calculation but not provided")]
    NoFastaProvided,

    /// Invalid peptide sequence encountered.
    #[error("invalid peptide sequence: {0}")]
    InvalidPeptide(String),
}
