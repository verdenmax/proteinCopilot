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
    ///
    /// This is a terminal operation — no `end_stage()` needed after this.
    /// Idempotent: if no stage is active, this is a no-op.
    pub fn fail_stage(&mut self, detail: &str) {
        if self.current_stage_name.is_none() {
            return;
        }
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

    /// Run anomaly detection rules after search completes.
    ///
    /// Arguments come from `SearchResultSummary` and `SearchParams`:
    /// - `identification_rate`: PSMs at 1% FDR / total spectra (0.0 to 1.0)
    /// - `psms_at_1pct_fdr`: absolute count of PSMs passing 1% FDR
    /// - `decoy_count`: number of decoy PSM hits
    /// - `total_elapsed_sec`: total wall-clock seconds
    /// - `precursor_tolerance_ppm`: precursor mass tolerance in ppm (None if Da units)
    pub fn finalize(
        &mut self,
        identification_rate: Option<f64>,
        psms_at_1pct_fdr: Option<u64>,
        decoy_count: Option<u64>,
        total_elapsed_sec: f64,
        precursor_tolerance_ppm: Option<f64>,
    ) {
        // Guard against double-finalize: anomalies would accumulate
        if !self.anomalies.is_empty() || !self.suggestions.is_empty() {
            return;
        }

        self.total_elapsed_sec = total_elapsed_sec;

        // Rule 1: Low identification rate
        if let Some(rate) = identification_rate {
            if rate < 0.10 {
                self.anomalies.push(SearchAnomaly {
                    severity: "warning".to_string(),
                    category: AnomalyCategory::LowIdentificationRate,
                    message: format!(
                        "PSM identification rate is {:.1}%, well below typical 10-40%",
                        rate * 100.0
                    ),
                    metric_name: Some("identification_rate".to_string()),
                    metric_value: Some(rate),
                    expected_range: Some("10-40%".to_string()),
                });
                self.suggestions.push(DiagnosticSuggestion {
                    priority: 1,
                    action: "Confirm sample species matches database organism".to_string(),
                    reason: "Species mismatch is the most common cause of low identification rate"
                        .to_string(),
                    param_changes: None,
                });
                self.suggestions.push(DiagnosticSuggestion {
                    priority: 2,
                    action: "Check enzyme setting matches sample preparation".to_string(),
                    reason: "Wrong enzyme causes systematic peptide mismatch".to_string(),
                    param_changes: None,
                });
            }
        }

        // Rule 2: No decoy hits
        if let Some(0) = decoy_count {
            self.anomalies.push(SearchAnomaly {
                severity: "error".to_string(),
                category: AnomalyCategory::NoDecoyHits,
                message: "No decoy PSM hits detected — FDR calculation is unreliable".to_string(),
                metric_name: Some("decoy_count".to_string()),
                metric_value: Some(0.0),
                expected_range: Some("> 0".to_string()),
            });
            self.suggestions.push(DiagnosticSuggestion {
                priority: 1,
                action: "Verify FASTA database contains decoy sequences".to_string(),
                reason: "Target-decoy FDR requires decoy sequences in the database".to_string(),
                param_changes: None,
            });
        }

        // Rule 3: Very few PSMs at 1% FDR
        if let Some(count) = psms_at_1pct_fdr {
            if count == 0 {
                self.anomalies.push(SearchAnomaly {
                    severity: "error".to_string(),
                    category: AnomalyCategory::HighFdr,
                    message: "No PSMs pass 1% FDR — search produced no reliable identifications"
                        .to_string(),
                    metric_name: Some("psms_at_1pct_fdr".to_string()),
                    metric_value: Some(0.0),
                    expected_range: Some("> 100".to_string()),
                });
                self.suggestions.push(DiagnosticSuggestion {
                    priority: 1,
                    action: "Check database species, enzyme, and modification settings".to_string(),
                    reason: "Zero identifications usually indicates fundamental parameter mismatch"
                        .to_string(),
                    param_changes: None,
                });
            } else if count < 50 {
                self.anomalies.push(SearchAnomaly {
                    severity: "warning".to_string(),
                    category: AnomalyCategory::HighFdr,
                    message: format!(
                        "Only {} PSMs pass 1% FDR — results may not be statistically reliable",
                        count
                    ),
                    metric_name: Some("psms_at_1pct_fdr".to_string()),
                    metric_value: Some(count as f64),
                    expected_range: Some("> 100".to_string()),
                });
            }
        }

        // Rule 4: Narrow precursor tolerance + low identification
        if let (Some(tol), Some(rate)) = (precursor_tolerance_ppm, identification_rate) {
            if tol < 5.0 && rate < 0.10 {
                self.anomalies.push(SearchAnomaly {
                    severity: "warning".to_string(),
                    category: AnomalyCategory::NarrowTolerance,
                    message: format!(
                        "Precursor tolerance {:.1} ppm is narrow and identification rate is low ({:.1}%)",
                        tol, rate * 100.0
                    ),
                    metric_name: Some("precursor_tolerance_ppm".to_string()),
                    metric_value: Some(tol),
                    expected_range: Some("5-20 ppm".to_string()),
                });
                self.suggestions.push(DiagnosticSuggestion {
                    priority: 3,
                    action: "Widen precursor mass tolerance to 10-20 ppm".to_string(),
                    reason: "Narrow tolerance may exclude correct matches on uncalibrated instruments".to_string(),
                    param_changes: Some(HashMap::from([(
                        "precursor_tolerance".to_string(),
                        serde_json::json!({"value": 20.0, "unit": "ppm"}),
                    )])),
                });
            }
        }

        // Rule 5: Wide precursor tolerance
        if let Some(tol) = precursor_tolerance_ppm {
            if tol > 50.0 {
                self.anomalies.push(SearchAnomaly {
                    severity: "warning".to_string(),
                    category: AnomalyCategory::WideTolerance,
                    message: format!(
                        "Precursor tolerance {:.1} ppm is very wide — FDR may be unreliable",
                        tol
                    ),
                    metric_name: Some("precursor_tolerance_ppm".to_string()),
                    metric_value: Some(tol),
                    expected_range: Some("5-20 ppm".to_string()),
                });
            }
        }

        // Rule 6: Database mismatch heuristic
        if let Some(rate) = identification_rate {
            if rate < 0.05 && self.error_category.is_none() {
                let already_has = self
                    .anomalies
                    .iter()
                    .any(|a| a.category == AnomalyCategory::DatabaseMismatch);
                if !already_has {
                    self.anomalies.push(SearchAnomaly {
                        severity: "warning".to_string(),
                        category: AnomalyCategory::DatabaseMismatch,
                        message: "Identification rate < 5% — database species may not match sample"
                            .to_string(),
                        metric_name: Some("identification_rate".to_string()),
                        metric_value: Some(rate),
                        expected_range: Some("> 10%".to_string()),
                    });
                }
            }
        }

        // Rule 7: Slow search (matching > 90% of total time)
        if total_elapsed_sec > 10.0 {
            if let Some(matching) = self.stages.iter().find(|s| s.name == "matching") {
                if matching.elapsed_sec / total_elapsed_sec > 0.90 {
                    self.anomalies.push(SearchAnomaly {
                        severity: "warning".to_string(),
                        category: AnomalyCategory::SlowSearch,
                        message: format!(
                            "Matching phase took {:.0}s ({:.0}% of total) — consider reducing search space",
                            matching.elapsed_sec,
                            matching.elapsed_sec / total_elapsed_sec * 100.0
                        ),
                        metric_name: Some("matching_time_pct".to_string()),
                        metric_value: Some(matching.elapsed_sec / total_elapsed_sec),
                        expected_range: Some("< 90%".to_string()),
                    });
                    self.suggestions.push(DiagnosticSuggestion {
                        priority: 4,
                        action: "Use Sage engine for parallel matching or reduce variable modifications".to_string(),
                        reason: "Matching is the bottleneck; parallelism or smaller search space helps".to_string(),
                        param_changes: Some(HashMap::from([(
                            "engine".to_string(),
                            serde_json::json!("Sage"),
                        )])),
                    });
                }
            }
        }

        // Sort suggestions by priority
        self.suggestions.sort_by_key(|s| s.priority);
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

    #[test]
    fn test_finalize_detects_low_identification_rate() {
        let mut d = SearchDiagnostics::new();
        d.begin_stage("matching");
        d.end_stage(Some(10000));
        // 5% identification rate
        d.finalize(Some(0.05), Some(500), Some(50), 10.0, None);

        assert!(d.anomalies.iter().any(|a| a.category == AnomalyCategory::LowIdentificationRate));
        assert!(!d.suggestions.is_empty());
    }

    #[test]
    fn test_finalize_detects_no_decoy_hits() {
        let mut d = SearchDiagnostics::new();
        d.begin_stage("fdr_calculation");
        d.end_stage(None);
        // 0 decoy hits
        d.finalize(Some(0.30), Some(3000), Some(0), 5.0, None);

        assert!(d.anomalies.iter().any(|a| a.category == AnomalyCategory::NoDecoyHits));
    }

    #[test]
    fn test_finalize_detects_narrow_tolerance() {
        let mut d = SearchDiagnostics::new();
        d.begin_stage("matching");
        d.end_stage(Some(5000));
        // 3 ppm + low rate
        d.finalize(Some(0.02), Some(100), Some(10), 8.0, Some(3.0));

        assert!(d.anomalies.iter().any(|a| a.category == AnomalyCategory::NarrowTolerance));
    }

    #[test]
    fn test_finalize_no_anomalies_for_good_search() {
        let mut d = SearchDiagnostics::new();
        d.begin_stage("matching");
        d.end_stage(Some(10000));
        // 25% rate, plenty of decoys, 10 ppm
        d.finalize(Some(0.25), Some(2500), Some(200), 15.0, Some(10.0));

        assert!(d.anomalies.is_empty());
        assert!(d.suggestions.is_empty());
    }

    #[test]
    fn test_finalize_detects_slow_search() {
        let mut d = SearchDiagnostics::new();
        // Manually push stages (not using begin/end to avoid timing issues)
        d.stages.push(DiagnosticStage {
            name: "file_reading".to_string(),
            status: "completed".to_string(),
            elapsed_sec: 1.0,
            items_processed: Some(1000),
            items_total: None,
            detail: None,
        });
        d.stages.push(DiagnosticStage {
            name: "matching".to_string(),
            status: "completed".to_string(),
            elapsed_sec: 95.0,
            items_processed: Some(1000),
            items_total: None,
            detail: None,
        });
        d.total_elapsed_sec = 100.0;
        d.finalize(Some(0.25), Some(250), Some(30), 100.0, Some(10.0));

        assert!(d.anomalies.iter().any(|a| a.category == AnomalyCategory::SlowSearch));
    }
}