//! # ProteinCopilot Parameter Recommendation Engine
//!
//! Deterministic rule-based engine that recommends search parameters
//! based on spectrum file characteristics. This is a pure library crate
//! with no MCP, network, or LLM dependencies.
//!
//! # Architecture Role
//!
//! ```text
//! SpectrumSummary ──┐
//!                   ├──▶ ParamRecommender::recommend() ──▶ AiDecision<SearchParams>
//! UserHints ────────┘
//! ```
//!
//! The LLM layer calls `recommend_params` MCP tool, which delegates to
//! this engine. The output `AiDecision<SearchParams>` contains the
//! recommended parameters along with confidence, explanation, and
//! alternatives — ready for the LLM to present to the user.

pub mod error;
pub mod hints;
pub mod preset;
pub mod rules;

pub use error::ParamRecommendError;
pub use hints::UserHints;
pub use preset::SearchPreset;

use protein_copilot_core::ai_decision::AiDecision;
use protein_copilot_core::search_params::SearchParams;
use protein_copilot_core::spectrum::SpectrumSummary;

/// Deterministic parameter recommendation engine.
///
/// All recommendation logic is pure functions: same input → same output.
/// No LLM calls, no network I/O, no randomness.
pub struct ParamRecommender;

impl ParamRecommender {
    /// Recommends search parameters based on spectrum characteristics.
    ///
    /// # Arguments
    /// - `summary` — Statistical summary of the spectrum file
    /// - `hints` — Optional user hints (experiment type, instrument, notes)
    ///
    /// # Returns
    /// An `AiDecision<SearchParams>` containing the recommendation with
    /// confidence, explanation, evidence, and alternatives.
    pub fn recommend(
        &self,
        summary: &SpectrumSummary,
        hints: Option<&UserHints>,
    ) -> Result<AiDecision<SearchParams>, ParamRecommendError> {
        rules::recommend(summary, hints)
    }

    /// Returns all built-in search parameter presets.
    pub fn list_presets() -> Vec<SearchPreset> {
        preset::all_presets()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn list_presets_returns_4() {
        let presets = ParamRecommender::list_presets();
        assert_eq!(presets.len(), 4);
    }
}
