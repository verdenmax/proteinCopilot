//! Search progress tracking (shared data structure).
//!
//! Defines [`SearchProgress`] for reporting search status and stage,
//! plus [`ProgressCallback`] for engines to report progress.

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use uuid::Uuid;

/// Progress information for a running search.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SearchProgress {
    /// Unique identifier for the search run.
    pub run_id: Uuid,
    /// Current status: `"Running"`, `"Completed"`, `"Failed: ..."`, or `"Cancelled"`.
    pub status: String,
    /// Current search stage, e.g. `"Matching spectra (300/1000)"`.
    pub stage: Option<String>,
    /// Progress percentage (0.0 to 1.0), `None` if indeterminate.
    pub progress_pct: Option<f64>,
    /// Elapsed time in seconds.
    pub elapsed_sec: f64,
    /// Estimated remaining time in seconds, `None` if unknown.
    pub estimated_remaining_sec: Option<f64>,
}

/// Callback type for progress reporting from search engines.
///
/// Search engine adapters call this at each stage to report progress.
/// The MCP server layer captures these updates and writes them into
/// the run cache for `get_search_status` queries.
pub type ProgressCallback = Box<dyn Fn(SearchProgress) + Send + Sync>;

/// A no-op progress callback for cases where progress reporting is not needed.
pub fn noop_progress() -> ProgressCallback {
    Box::new(|_| {})
}
