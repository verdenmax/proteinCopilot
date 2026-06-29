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

    /// Error classification (set only when status starts with "Failed").
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub error_category: Option<crate::diagnostics::ErrorCategory>,

    /// Whether detailed diagnostics are available via diagnose_search.
    #[serde(default)]
    pub has_diagnostics: bool,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn search_progress_roundtrip_without_error_category() {
        let progress = SearchProgress {
            run_id: Uuid::new_v4(),
            status: "Running".to_string(),
            stage: Some("Matching spectra".to_string()),
            progress_pct: Some(0.5),
            elapsed_sec: 12.3,
            estimated_remaining_sec: None,
            error_category: None,
            has_diagnostics: false,
        };

        let json = serde_json::to_string(&progress).unwrap();
        assert!(
            !json.contains("error_category"),
            "None field should be skipped"
        );

        let restored: SearchProgress = serde_json::from_str(&json).unwrap();
        assert_eq!(restored.status, "Running");
        assert!(restored.error_category.is_none());
    }
}
