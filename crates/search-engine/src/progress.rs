//! Search progress tracking.

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use uuid::Uuid;

/// Progress information for a running search.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SearchProgress {
    /// Unique identifier for the search run.
    pub run_id: Uuid,
    /// Current status description.
    pub status: String,
    /// Progress percentage (0.0 to 1.0), `None` if indeterminate.
    pub progress_pct: Option<f64>,
    /// Elapsed time in seconds.
    pub elapsed_sec: f64,
    /// Estimated remaining time in seconds, `None` if unknown.
    pub estimated_remaining_sec: Option<f64>,
}
