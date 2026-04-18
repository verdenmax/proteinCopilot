use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::time::Instant;

/// Error classification for failed searches.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub enum ErrorCategory {
    /// Input data problem: file corrupt, wrong format, no MS2 spectra
    InputData,
    /// Search parameter mismatch: tolerance wrong, enzyme mismatch
    Parameters,
    /// Database problem: species mismatch, FASTA format error, too small
    Database,
    /// Search engine internal error or resource exhaustion (OOM)
    Engine,
}

/// Quality anomaly classification for completed searches.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub enum AnomalyCategory {
    /// PSM identification rate < 10%
    LowIdentificationRate,
    /// FDR distribution abnormal (too many decoys)
    HighFdr,
    /// No decoy hits (FDR cannot be calculated)
    NoDecoyHits,
    /// Precursor tolerance too narrow, few candidates
    NarrowTolerance,
    /// Precursor tolerance too wide, FDR unreliable
    WideTolerance,
    /// Median fragment ion count too low
    LowSpectraQuality,
    /// Species likely mismatched with database
    DatabaseMismatch,
    /// Matching phase consumed > 90% of total time
    SlowSearch,
}

/// Metrics captured for one search phase.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct DiagnosticStage {
    /// Stage name: "file_reading", "fasta_parsing", "digestion", "matching", "fdr_calculation"
    pub name: String,
    /// "completed", "failed", or "skipped"
    pub status: String,
    /// Wall-clock seconds for this stage
    pub elapsed_sec: f64,
    /// Items processed (spectra, proteins, peptides, PSMs — context-dependent)
    pub items_processed: Option<u64>,
    /// Total items expected (if known)
    pub items_total: Option<u64>,
    /// Extra detail (error message on failure, summary on completion)
    pub detail: Option<String>,
}

/// A detected quality anomaly.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SearchAnomaly {
    /// "warning" or "error"
    pub severity: String,
    /// Anomaly classification
    pub category: AnomalyCategory,
    /// Human-readable description
    pub message: String,
    /// Metric name (e.g. "identification_rate")
    pub metric_name: Option<String>,
    /// Actual metric value
    pub metric_value: Option<f64>,
    /// Expected range (e.g. "10-40%")
    pub expected_range: Option<String>,
}

/// A repair/optimization suggestion.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct DiagnosticSuggestion {
    /// Priority 1 (highest) to 5 (lowest)
    pub priority: u8,
    /// What the user should do
    pub action: String,
    /// Why this is suggested
    pub reason: String,
    /// Concrete parameter changes (field name → new value)
    pub param_changes: Option<HashMap<String, serde_json::Value>>,
}

/// Complete diagnostics report for a search run.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SearchDiagnostics {
    /// Error classification (only set when search failed)
    pub error_category: Option<ErrorCategory>,
    /// Stage where failure occurred
    pub failure_stage: Option<String>,
    /// Error detail string
    pub error_detail: Option<String>,
    /// Per-stage metrics (always collected)
    pub stages: Vec<DiagnosticStage>,
    /// Detected anomalies (populated by finalize)
    pub anomalies: Vec<SearchAnomaly>,
    /// Repair suggestions (populated by finalize)
    pub suggestions: Vec<DiagnosticSuggestion>,
    /// Total search duration in seconds
    pub total_elapsed_sec: f64,

    /// Internal: timestamp when current stage began (not serialized)
    #[serde(skip)]
    #[schemars(skip)]
    current_stage_start: Option<Instant>,
    /// Internal: name of current stage
    #[serde(skip)]
    #[schemars(skip)]
    current_stage_name: Option<String>,
}

impl SearchDiagnostics {
    /// Create a new empty diagnostics collector.
    pub fn new() -> Self {
        Self {
            error_category: None,
            failure_stage: None,
            error_detail: None,
            stages: Vec::new(),
            anomalies: Vec::new(),
            suggestions: Vec::new(),
            total_elapsed_sec: 0.0,
            current_stage_start: None,
            current_stage_name: None,
        }
    }

    /// Mark the beginning of a named stage.
    pub fn begin_stage(&mut self, name: &str) {
        self.current_stage_start = Some(Instant::now());
        self.current_stage_name = Some(name.to_string());
    }

    /// Mark the current stage as completed.
    pub fn end_stage(&mut self, items_processed: Option<u64>) {
        let elapsed = self
            .current_stage_start
            .map(|s| s.elapsed().as_secs_f64())
            .unwrap_or(0.0);
        let name = self
            .current_stage_name
            .take()
            .unwrap_or_else(|| "unknown".to_string());
        self.current_stage_start = None;
        self.stages.push(DiagnosticStage {
            name,
            status: "completed".to_string(),
            elapsed_sec: elapsed,
            items_processed,
            items_total: None,
            detail: None,
        });
    }

    /// Mark the current stage as failed.
    pub fn fail_stage(&mut self, detail: &str) {
        let elapsed = self
            .current_stage_start
            .map(|s| s.elapsed().as_secs_f64())
            .unwrap_or(0.0);
        let name = self
            .current_stage_name
            .take()
            .unwrap_or_else(|| "unknown".to_string());
        self.current_stage_start = None;
        self.failure_stage = Some(name.clone());
        self.stages.push(DiagnosticStage {
            name,
            status: "failed".to_string(),
            elapsed_sec: elapsed,
            items_processed: None,
            items_total: None,
            detail: Some(detail.to_string()),
        });
    }

    /// Set the top-level error classification.
    pub fn set_error(&mut self, category: ErrorCategory, detail: &str) {
        self.error_category = Some(category);
        self.error_detail = Some(detail.to_string());
    }
}

impl Default for SearchDiagnostics {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_diagnostics_is_empty() {
        let d = SearchDiagnostics::new();
        assert!(d.error_category.is_none());
        assert!(d.failure_stage.is_none());
        assert!(d.stages.is_empty());
        assert!(d.anomalies.is_empty());
        assert!(d.suggestions.is_empty());
    }

    #[test]
    fn test_stage_tracking() {
        let mut d = SearchDiagnostics::new();
        d.begin_stage("file_reading");
        std::thread::sleep(std::time::Duration::from_millis(10));
        d.end_stage(Some(1000));

        assert_eq!(d.stages.len(), 1);
        assert_eq!(d.stages[0].name, "file_reading");
        assert_eq!(d.stages[0].status, "completed");
        assert_eq!(d.stages[0].items_processed, Some(1000));
        assert!(d.stages[0].elapsed_sec >= 0.01);
    }

    #[test]
    fn test_fail_stage() {
        let mut d = SearchDiagnostics::new();
        d.begin_stage("fasta_parsing");
        d.fail_stage("Invalid header at line 42");

        assert_eq!(d.stages.len(), 1);
        assert_eq!(d.stages[0].status, "failed");
        assert_eq!(
            d.stages[0].detail.as_deref(),
            Some("Invalid header at line 42")
        );
    }

    #[test]
    fn test_set_error() {
        let mut d = SearchDiagnostics::new();
        d.set_error(ErrorCategory::Database, "FASTA parse failure");

        assert_eq!(d.error_category, Some(ErrorCategory::Database));
        assert_eq!(d.error_detail.as_deref(), Some("FASTA parse failure"));
    }

    #[test]
    fn test_multiple_stages_tracked_in_order() {
        let mut d = SearchDiagnostics::new();
        d.begin_stage("file_reading");
        d.end_stage(Some(500));
        d.begin_stage("fasta_parsing");
        d.end_stage(Some(20000));
        d.begin_stage("matching");
        d.end_stage(Some(500));

        assert_eq!(d.stages.len(), 3);
        assert_eq!(d.stages[0].name, "file_reading");
        assert_eq!(d.stages[1].name, "fasta_parsing");
        assert_eq!(d.stages[2].name, "matching");
    }

    #[test]
    fn test_serialization_roundtrip() {
        let mut d = SearchDiagnostics::new();
        d.begin_stage("matching");
        d.end_stage(Some(100));
        d.set_error(ErrorCategory::InputData, "No MS2 spectra");

        let json = serde_json::to_string(&d).unwrap();
        let d2: SearchDiagnostics = serde_json::from_str(&json).unwrap();
        assert_eq!(d2.error_category, Some(ErrorCategory::InputData));
        assert_eq!(d2.stages.len(), 1);
    }
}