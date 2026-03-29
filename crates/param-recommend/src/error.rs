//! Error types for the parameter recommendation engine.

use thiserror::Error;

/// Errors that can occur during parameter recommendation.
#[derive(Debug, Error)]
pub enum ParamRecommendError {
    /// The input spectrum summary is empty (no spectra).
    #[error("spectrum summary is empty — cannot infer parameters from an empty file")]
    EmptySummary,

    /// A required field in the summary is invalid.
    #[error("invalid summary field '{field}': {detail}")]
    InvalidSummary {
        /// Which field is problematic.
        field: &'static str,
        /// What is wrong with it.
        detail: String,
    },
}
