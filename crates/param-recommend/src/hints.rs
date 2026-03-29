//! User hints — structured input from the LLM layer.
//!
//! The LLM translates user natural language intent (e.g., "search this
//! phospho dataset") into a [`UserHints`] struct, which the rule engine
//! uses to adjust its recommendations.

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

/// Structured hints from the user (via LLM translation).
///
/// All fields are optional — the rule engine has sensible defaults.
/// When provided, hints override or adjust the automatic recommendation.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct UserHints {
    /// Experiment type (e.g., "phosphorylation", "TMT", "SILAC", "standard").
    pub experiment_type: Option<String>,
    /// Instrument type (e.g., "Orbitrap", "TOF", "QExactive").
    pub instrument_type: Option<String>,
    /// Free-form notes from the user (e.g., "use 5ppm tolerance").
    pub custom_notes: Option<String>,
}
