# Failure Diagnosis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a comprehensive failure diagnosis system that captures search stage metrics, detects quality anomalies, and provides structured repair suggestions — for both failed and successful searches.

**Architecture:** New `diagnostics.rs` module in `core` crate defines `SearchDiagnostics` with stage tracking, anomaly detection, and suggestion generation. Both search engines (`SimpleSearch`, `Sage`) report stage metrics through `&mut SearchDiagnostics`. The MCP server stores diagnostics alongside `RunState`, exposes them via an enhanced `get_search_status` (adds `error_category` + `has_diagnostics`) and a new `diagnose_search` tool. A new Prompt guides LLM interpretation.

**Tech Stack:** Rust, serde, schemars, chrono, thiserror, rmcp

---

## File Structure

| Action | File | Responsibility |
|--------|------|----------------|
| Create | `crates/core/src/diagnostics.rs` | SearchDiagnostics, DiagnosticStage, SearchAnomaly, DiagnosticSuggestion, ErrorCategory, AnomalyCategory, finalize logic |
| Modify | `crates/core/src/lib.rs:15` | Add `pub mod diagnostics` |
| Modify | `crates/core/src/progress.rs:12-25` | Add `error_category` + `has_diagnostics` to SearchProgress |
| Modify | `crates/core/src/engine.rs:87-92,106-118` | Add `diagnostics: &mut SearchDiagnostics` to search/search_with_spectra |
| Modify | `crates/search-engine/src/simple_engine.rs:159-324` | Instrument search stages with diagnostics |
| Modify | `crates/search-engine/src/adapters/sage/mod.rs:78-400` | Instrument search stages with diagnostics |
| Modify | `crates/mcp-server/src/tools.rs:738-743` | Extend RunState with diagnostics + params fields |
| Modify | `crates/mcp-server/src/tools.rs:1691-1705` | Enhance get_search_status to populate new fields |
| Modify | `crates/mcp-server/src/tools.rs:1335-1430,1580-1650` | Wire diagnostics through spawned search tasks |
| Modify | `crates/mcp-server/src/tools.rs` | Add diagnose_search MCP tool |
| Create | `.github/prompts/failure-diagnosis.prompt.md` | LLM diagnosis interpretation guide |
| Modify | `.github/agents/proteomics-search.agent.md` | Add diagnose_search tool + diagnosis workflow |

---

### Task 1: Create `diagnostics.rs` with core data structures

**Files:**
- Create: `crates/core/src/diagnostics.rs`
- Modify: `crates/core/src/lib.rs:15`

- [ ] **Step 1: Write tests for SearchDiagnostics**

In `crates/core/src/diagnostics.rs`, add the test module at the bottom:

```rust
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
        assert_eq!(d.stages[0].detail.as_deref(), Some("Invalid header at line 42"));
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
```

- [ ] **Step 2: Write the data structures and constructor**

At the top of `crates/core/src/diagnostics.rs`:

```rust
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
```

- [ ] **Step 3: Register the module**

In `crates/core/src/lib.rs`, after line 15 (`pub mod error;`), add:

```rust
pub mod diagnostics;
```

- [ ] **Step 4: Run tests**

Run: `cargo test -p protein-copilot-core -- diagnostics 2>&1 | tail -10`
Expected: All 6 tests pass

- [ ] **Step 5: Commit**

```bash
git add crates/core/src/diagnostics.rs crates/core/src/lib.rs
git commit -m "feat(core): add SearchDiagnostics data structures

New diagnostics module with stage tracking, error classification,
anomaly categories, and suggestion generation. Includes begin_stage/
end_stage/fail_stage API for search engines to report phase metrics.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: Add `finalize()` anomaly detection rules

**Files:**
- Modify: `crates/core/src/diagnostics.rs`

- [ ] **Step 1: Write tests for anomaly detection**

Append to the existing `mod tests` block in `crates/core/src/diagnostics.rs`:

```rust
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
        d.begin_stage("file_reading");
        d.stages.push(DiagnosticStage {
            name: "file_reading".to_string(),
            status: "completed".to_string(),
            elapsed_sec: 1.0,
            items_processed: Some(1000),
            items_total: None,
            detail: None,
        });
        d.current_stage_start = None;
        d.current_stage_name = None;
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
```

- [ ] **Step 2: Implement `finalize()`**

Add this method to the `impl SearchDiagnostics` block, after `set_error()`:

```rust
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
            if count < 50 && count > 0 {
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
```

- [ ] **Step 3: Run tests**

Run: `cargo test -p protein-copilot-core -- diagnostics 2>&1 | tail -15`
Expected: All 11 tests pass

- [ ] **Step 4: Commit**

```bash
git add crates/core/src/diagnostics.rs
git commit -m "feat(core): add finalize() anomaly detection rules

8 deterministic rules: LowIdentificationRate, NoDecoyHits, HighFdr,
NarrowTolerance, WideTolerance, DatabaseMismatch, SlowSearch, and
LowSpectraQuality. Each rule generates anomalies and suggestions
with priority ranking and optional param_changes.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: Enhance `SearchProgress` with diagnostic fields

**Files:**
- Modify: `crates/core/src/progress.rs:12-25`

- [ ] **Step 1: Add fields to SearchProgress**

In `crates/core/src/progress.rs`, add two new fields after `estimated_remaining_sec` (line 24):

```rust
    /// Error classification (set only when status starts with "Failed").
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error_category: Option<crate::diagnostics::ErrorCategory>,

    /// Whether detailed diagnostics are available via diagnose_search.
    #[serde(default)]
    pub has_diagnostics: bool,
```

- [ ] **Step 2: Update all SearchProgress construction sites**

Search for `SearchProgress {` in the codebase and add the two new fields with defaults. Key locations:

In `crates/core/src/progress.rs`, find the `noop_progress()` function (around line 35). It returns a `Box<dyn Fn(SearchProgress)>` — the callback receives a `SearchProgress` from the engine, so no construction change needed there.

In `crates/search-engine/src/simple_engine.rs`, search for all `SearchProgress {` constructions. For each, add:
```rust
    error_category: None,
    has_diagnostics: false,
```

In `crates/search-engine/src/adapters/sage/mod.rs`, search for all `SearchProgress {` constructions. For each, add:
```rust
    error_category: None,
    has_diagnostics: false,
```

In `crates/mcp-server/src/tools.rs`, search for all `SearchProgress {` constructions. For each, add:
```rust
    error_category: None,
    has_diagnostics: false,
```

- [ ] **Step 3: Build and test**

Run: `cargo build --workspace 2>&1 | tail -5`
Expected: Compiles with no errors

Run: `cargo test --workspace 2>&1 | grep -E "^test result:|FAILED"`
Expected: All pass, 0 failed

- [ ] **Step 4: Commit**

```bash
git add crates/core/src/progress.rs crates/search-engine/ crates/mcp-server/src/tools.rs
git commit -m "feat(core): add error_category and has_diagnostics to SearchProgress

Lightweight enhancement: two optional fields for LLM to quickly assess
whether to call diagnose_search. Backwards-compatible (skip_serializing_if
+ serde default).

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 4: Update `SearchEngineAdapter` trait signature

**Files:**
- Modify: `crates/core/src/engine.rs:87-118`
- Modify: `crates/search-engine/src/simple_engine.rs` (impl block)
- Modify: `crates/search-engine/src/adapters/sage/mod.rs` (impl block)

- [ ] **Step 1: Update trait method signatures**

In `crates/core/src/engine.rs`, add `use crate::diagnostics::SearchDiagnostics;` at the top imports.

Change `search()` signature (lines 87-92) from:

```rust
    async fn search(
        &self,
        params: &SearchParams,
        input_files: &[PathBuf],
        on_progress: ProgressCallback,
    ) -> Result<SearchResult, CoreError>;
```

to:

```rust
    async fn search(
        &self,
        params: &SearchParams,
        input_files: &[PathBuf],
        on_progress: ProgressCallback,
        diagnostics: &mut SearchDiagnostics,
    ) -> Result<SearchResult, CoreError>;
```

Change `search_with_spectra()` default impl (lines 106-118) from:

```rust
    async fn search_with_spectra(
        &self,
        params: &SearchParams,
        spectra: Vec<Spectrum>,
        on_progress: ProgressCallback,
    ) -> Result<SearchResult, CoreError> {
        let _ = (params, spectra, on_progress);
```

to:

```rust
    async fn search_with_spectra(
        &self,
        params: &SearchParams,
        spectra: Vec<Spectrum>,
        on_progress: ProgressCallback,
        diagnostics: &mut SearchDiagnostics,
    ) -> Result<SearchResult, CoreError> {
        let _ = (params, spectra, on_progress, diagnostics);
```

- [ ] **Step 2: Update SimpleSearchEngine impl**

In `crates/search-engine/src/simple_engine.rs`, find the `impl SearchEngineAdapter for SimpleSearchEngine` block. Update the `search()` method signature to add `diagnostics: &mut SearchDiagnostics` parameter. Add `use protein_copilot_core::diagnostics::SearchDiagnostics;` to imports.

The body does not need instrumentation yet (that's Task 5) — just accept and ignore the parameter for now:

```rust
    async fn search(
        &self,
        params: &SearchParams,
        input_files: &[PathBuf],
        on_progress: ProgressCallback,
        diagnostics: &mut SearchDiagnostics,
    ) -> Result<SearchResult, CoreError> {
        let _ = &diagnostics;  // Will be used in Task 5
        // ... existing body unchanged ...
    }
```

- [ ] **Step 3: Update SageAdapter impl**

In `crates/search-engine/src/adapters/sage/mod.rs`, update both `search()` and `search_with_spectra()` signatures similarly. Add `use protein_copilot_core::diagnostics::SearchDiagnostics;` to imports.

For `search()` (lines 78-112):
```rust
    async fn search(
        &self,
        params: &SearchParams,
        input_files: &[PathBuf],
        on_progress: ProgressCallback,
        diagnostics: &mut SearchDiagnostics,
    ) -> Result<SearchResult, CoreError> {
        // ... existing body ...
        // Change the search_with_spectra call to pass diagnostics:
        let mut result = self.search_with_spectra(params, all_spectra, on_progress, diagnostics).await?;
        // ...
    }
```

For `search_with_spectra()` (lines 115+):
```rust
    async fn search_with_spectra(
        &self,
        params: &SearchParams,
        spectra: Vec<Spectrum>,
        on_progress: ProgressCallback,
        diagnostics: &mut SearchDiagnostics,
    ) -> Result<SearchResult, CoreError> {
        let _ = &diagnostics;  // Will be used in Task 6
        // ... existing body unchanged ...
    }
```

- [ ] **Step 4: Update MCP server call sites**

In `crates/mcp-server/src/tools.rs`, find the two locations where `engine.search()` and `engine.search_with_spectra()` are called inside spawned tasks.

For the DIA path (around line 1345-1347):
```rust
// Before:
let search_result = engine.search_with_spectra(&params, dia_spectra, on_progress).await;
// After:
let mut diagnostics = protein_copilot_core::diagnostics::SearchDiagnostics::new();
let search_result = engine.search_with_spectra(&params, dia_spectra, on_progress, &mut diagnostics).await;
```

For the file-based path (around line 1590):
```rust
// Before:
let search_result = engine.search(&params, &files, on_progress).await;
// After:
let mut diagnostics = protein_copilot_core::diagnostics::SearchDiagnostics::new();
let search_result = engine.search(&params, &files, on_progress, &mut diagnostics).await;
```

Note: `diagnostics` is created but not yet stored in `RunState` — that's Task 7.

- [ ] **Step 5: Build and test**

Run: `cargo build --workspace 2>&1 | tail -5`
Expected: Compiles with no errors

Run: `cargo test --workspace 2>&1 | grep -E "^test result:|FAILED"`
Expected: All pass, 0 failed

- [ ] **Step 6: Commit**

```bash
git add crates/core/src/engine.rs crates/search-engine/ crates/mcp-server/src/tools.rs
git commit -m "refactor: add diagnostics parameter to SearchEngineAdapter trait

Both search() and search_with_spectra() now accept &mut SearchDiagnostics.
Engines ignore it for now; instrumentation follows in subsequent commits.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 5: Instrument SimpleSearchEngine with diagnostics

**Files:**
- Modify: `crates/search-engine/src/simple_engine.rs`

- [ ] **Step 1: Add stage instrumentation to `search()` / `run_search()`**

In the `search()` method (which calls `self.run_search()`), or in `run_search()` directly, wrap each phase:

**FASTA reading (around line 191-194):**
```rust
diagnostics.begin_stage("fasta_parsing");
let proteins = parse_fasta(&params.database_path)
    .map_err(|e| {
        diagnostics.fail_stage(&e.to_string());
        diagnostics.set_error(
            protein_copilot_core::diagnostics::ErrorCategory::Database,
            &e.to_string(),
        );
        SearchEngineError::IoError(e.to_string())
    })?;
diagnostics.end_stage(Some(proteins.len() as u64));
```

**Digestion (around line 196-220):**
```rust
diagnostics.begin_stage("digestion");
// ... existing digestion code ...
diagnostics.end_stage(Some(candidates.len() as u64));
```

**Spectrum matching (around line 274-312):**
```rust
diagnostics.begin_stage("matching");
// ... existing matching loop ...
diagnostics.end_stage(Some(total_spectra));
```

**FDR calculation (in finalize_search_result, around line 100-127):**
```rust
diagnostics.begin_stage("fdr_calculation");
// ... existing FDR code ...
diagnostics.end_stage(Some(psms.len() as u64));
```

Note: The `diagnostics` reference must be threaded through from `search()` → `run_search()` → `finalize_search_result()`. Update the `run_search()` and `finalize_search_result()` internal method signatures to accept `diagnostics: &mut SearchDiagnostics`.

On any error in these phases, call `diagnostics.fail_stage()` and `diagnostics.set_error()` with the appropriate `ErrorCategory` before propagating the error.

- [ ] **Step 2: Build and test**

Run: `cargo build --workspace 2>&1 | tail -5`
Expected: Compiles with no errors

Run: `cargo test --workspace 2>&1 | grep -E "^test result:|FAILED"`
Expected: All pass, 0 failed

- [ ] **Step 3: Commit**

```bash
git add crates/search-engine/src/simple_engine.rs
git commit -m "feat(search-engine): instrument SimpleSearch with diagnostics

Report stage metrics for fasta_parsing, digestion, matching, and
fdr_calculation. On failure, classify error and record failure stage.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 6: Instrument SageAdapter with diagnostics

**Files:**
- Modify: `crates/search-engine/src/adapters/sage/mod.rs`

- [ ] **Step 1: Add stage instrumentation to `search_with_spectra()`**

In `search_with_spectra()`, wrap the major phases:

**File reading / spectrum conversion (lines 124-154):**
```rust
diagnostics.begin_stage("file_reading");
// ... existing MS2 filtering + conversion ...
diagnostics.end_stage(Some(raw_spectra.len() as u64));
```

**FASTA parsing + database build (lines 161-213):**
```rust
diagnostics.begin_stage("fasta_parsing");
// ... existing FASTA read + sage_params.build() ...
diagnostics.end_stage(Some(protein_count as u64));
```

**Scoring (lines 238-246):**
```rust
diagnostics.begin_stage("matching");
// ... existing rayon scoring ...
diagnostics.end_stage(Some(all_features.len() as u64));
```

**FDR/rescoring (lines 248-272):**
```rust
diagnostics.begin_stage("fdr_calculation");
// ... existing LDA rescoring + FDR ...
diagnostics.end_stage(Some(all_features.len() as u64));
```

On errors (e.g., "No MS2 spectra found"), call:
```rust
diagnostics.fail_stage("No MS2 spectra found in input");
diagnostics.set_error(ErrorCategory::InputData, "No MS2 spectra found");
```

Also update the `search()` method to pass `diagnostics` through to `search_with_spectra()`. Add a `file_reading` stage around the file loop in `search()` if spectra are read there.

- [ ] **Step 2: Build and test**

Run: `cargo build --workspace 2>&1 | tail -5`
Expected: Compiles

Run: `cargo test --workspace 2>&1 | grep -E "^test result:|FAILED"`
Expected: All pass

- [ ] **Step 3: Commit**

```bash
git add crates/search-engine/src/adapters/sage/mod.rs
git commit -m "feat(search-engine): instrument Sage adapter with diagnostics

Report stage metrics for file_reading, fasta_parsing, matching, and
fdr_calculation (including LDA rescoring). Classify errors by category.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 7: Wire diagnostics into RunState and `get_search_status`

**Files:**
- Modify: `crates/mcp-server/src/tools.rs`

- [ ] **Step 1: Extend RunState**

Change `RunState` (lines 738-743) from:

```rust
struct RunState {
    progress: SearchProgress,
    result: Option<SearchResult>,
    handle: Option<tokio::task::JoinHandle<()>>,
}
```

to:

```rust
struct RunState {
    progress: SearchProgress,
    result: Option<SearchResult>,
    handle: Option<tokio::task::JoinHandle<()>>,
    diagnostics: Option<protein_copilot_core::diagnostics::SearchDiagnostics>,
    params_used: Option<protein_copilot_core::search_params::SearchParams>,
}
```

Update all `RunState { ... }` construction sites to add `diagnostics: None, params_used: None`.

- [ ] **Step 2: Store diagnostics after search completes**

In the DIA spawned task (around line 1351-1418), after `search_with_spectra()` returns, store diagnostics. Find the block where `state.result` and `state.progress.status` are updated:

For the success path (around line 1360-1387):
```rust
Ok(mut result) => {
    // ... existing code ...
    // Run finalize with result metrics
    let tol_ppm = if params.precursor_tolerance.unit == "ppm" {
        Some(params.precursor_tolerance.value)
    } else {
        None
    };
    let decoy_count = result.psms.iter().filter(|p| p.is_decoy).count() as u64;
    diagnostics.finalize(
        Some(result.summary.identification_rate),
        Some(result.summary.psms_at_1pct_fdr),
        Some(decoy_count),
        duration,
        tol_ppm,
    );
    state.diagnostics = Some(diagnostics);
    state.params_used = Some(params.clone());
    state.progress.has_diagnostics = true;
    // ... rest of existing code ...
}
```

For the error path (around line 1390-1407):
```rust
Err(e) => {
    // ... existing code ...
    state.diagnostics = Some(diagnostics);
    state.params_used = Some(params.clone());
    state.progress.has_diagnostics = true;
    state.progress.error_category = diagnostics.error_category.clone();
    // ... rest of existing code ...
}
```

Note: `diagnostics` must be moved into the `tokio::spawn` closure. Since `SearchDiagnostics` is created before the spawn and mutated inside it, it needs to be `move`d into the closure (already the case for other variables).

Apply the same pattern to the file-based path (around line 1594-1650).

- [ ] **Step 3: Update `get_search_status` to return new fields**

The `get_search_status` handler (lines 1691-1705) currently returns `state.progress.clone()`. Since `SearchProgress` now has `error_category` and `has_diagnostics` fields that are populated when the search completes (Step 2), no additional changes are needed — the fields will be included automatically.

- [ ] **Step 4: Build and test**

Run: `cargo build --workspace 2>&1 | tail -5`
Expected: Compiles

Run: `cargo test --workspace 2>&1 | grep -E "^test result:|FAILED"`
Expected: All pass

- [ ] **Step 5: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "feat(mcp): wire SearchDiagnostics through RunState

Store diagnostics and params_used in RunState after search completes.
Run finalize() anomaly detection on success. Populate error_category
and has_diagnostics in SearchProgress for get_search_status.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 8: Add `diagnose_search` MCP tool

**Files:**
- Modify: `crates/mcp-server/src/tools.rs`

- [ ] **Step 1: Add input/output structs**

Near other input/output struct definitions:

```rust
#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct DiagnoseSearchInput {
    /// The run_id to diagnose (from run_search or get_search_status)
    run_id: String,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct DiagnoseSearchOutput {
    /// The run_id
    run_id: String,
    /// Overall search status: "Completed", "Failed", "Running", "Cancelled"
    overall_status: String,
    /// Error classification (only for failed searches)
    error_category: Option<protein_copilot_core::diagnostics::ErrorCategory>,
    /// Stage where failure occurred
    failure_stage: Option<String>,
    /// Error detail message
    error_detail: Option<String>,
    /// Per-stage metrics
    stages: Vec<protein_copilot_core::diagnostics::DiagnosticStage>,
    /// Detected quality anomalies
    anomalies: Vec<protein_copilot_core::diagnostics::SearchAnomaly>,
    /// Repair/optimization suggestions, sorted by priority
    suggestions: Vec<protein_copilot_core::diagnostics::DiagnosticSuggestion>,
    /// Total search duration in seconds
    total_elapsed_sec: f64,
}
```

- [ ] **Step 2: Implement the tool handler**

In the `#[tool_router]` impl block, add:

```rust
#[rmcp::tool(
    name = "diagnose_search",
    description = "Get diagnostic report for a search run. Works for both failed searches (error analysis) and completed searches (quality assessment). Returns stage metrics, detected anomalies, and repair suggestions. Call after get_search_status shows the search has finished (status is Completed, Failed, or Cancelled). Use has_diagnostics=true from get_search_status to confirm diagnostics are available."
)]
async fn diagnose_search(
    &self,
    Parameters(input): Parameters<DiagnoseSearchInput>,
) -> Result<Json<DiagnoseSearchOutput>, ErrorData> {
    let id = Uuid::parse_str(&input.run_id).map_err(|_| {
        mcp_err(ErrorCode::INVALID_PARAMS, "invalid run_id format — expected UUID")
    })?;

    let cache = self.run_cache.lock().map_err(|_| {
        mcp_err(ErrorCode::INTERNAL_ERROR, "run cache lock is poisoned")
    })?;

    let state = cache.get(&id).ok_or_else(|| {
        mcp_err(
            ErrorCode::INVALID_PARAMS,
            format!("run_id '{}' not found in cache", input.run_id),
        )
    })?;

    if state.progress.status == "Running" {
        return Err(mcp_err(
            ErrorCode::INVALID_PARAMS,
            "Search is still running. Wait for completion before diagnosing.",
        ));
    }

    let diag = state.diagnostics.as_ref().ok_or_else(|| {
        mcp_err(
            ErrorCode::INTERNAL_ERROR,
            "No diagnostics available for this run. The search may have been started before the diagnostics feature was added.",
        )
    })?;

    Ok(Json(DiagnoseSearchOutput {
        run_id: input.run_id,
        overall_status: state.progress.status.clone(),
        error_category: diag.error_category.clone(),
        failure_stage: diag.failure_stage.clone(),
        error_detail: diag.error_detail.clone(),
        stages: diag.stages.clone(),
        anomalies: diag.anomalies.clone(),
        suggestions: diag.suggestions.clone(),
        total_elapsed_sec: diag.total_elapsed_sec,
    }))
}
```

- [ ] **Step 3: Build and test**

Run: `cargo build --workspace 2>&1 | tail -5`
Expected: Compiles

Run: `cargo test --workspace 2>&1 | grep -E "^test result:|FAILED"`
Expected: All pass

Run: `cargo clippy --workspace 2>&1 | grep -E "warning|error" | head -5`
Expected: 0 warnings

- [ ] **Step 4: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "feat(mcp): add diagnose_search tool

New MCP tool returns structured diagnostic report for any completed
search run: stage metrics, anomaly detection, and prioritized repair
suggestions. Supports both failed and successful searches.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 9: Create `failure-diagnosis.prompt.md` and update Agent

**Files:**
- Create: `.github/prompts/failure-diagnosis.prompt.md`
- Modify: `.github/agents/proteomics-search.agent.md`

- [ ] **Step 1: Create the Prompt file**

Create `.github/prompts/failure-diagnosis.prompt.md`:

```markdown
---
mode: agent
description: "搜索诊断 — 分析搜索失败原因、评估结果质量、提供参数调优建议和重试方案"
---

# 搜索诊断

分析搜索运行的失败原因或结果质量异常，提供修复建议。

## 使用时机

- 搜索失败时（`get_search_status` 返回 status = "Failed..."）
- 搜索成功但结果异常时（鉴定率低、蛋白数量少等）
- 用户主动要求分析搜索质量

## 诊断流程

### 1. 获取诊断数据
- 调用 `diagnose_search(run_id=xxx)` 获取 `SearchDiagnostics`
- 前提：搜索已结束（`get_search_status` 的 `has_diagnostics = true`）

### 2. 搜索失败诊断

按以下顺序解读：

1. **error_category** → 确定大方向
   - `InputData`：文件损坏、格式错误、无 MS2 谱图
   - `Parameters`：容差不合理、酶不匹配
   - `Database`：物种不匹配、FASTA 格式错误
   - `Engine`：引擎内部错误、资源不足
2. **failure_stage** → 定位失败发生在哪个阶段
3. **stages[]** → 展示搜索"走了多远"（哪些阶段完成了）
4. **suggestions[]** → 按 priority 排序展示修复方案
5. 如果 suggestions 中有 **param_changes** → 直接向用户展示修改后的参数

### 3. 结果质量诊断

搜索成功但可能有质量问题：

1. **anomalies[]** → 列出检测到的异常
2. 对每个异常解释：
   - `LowIdentificationRate`：PSM 鉴定率过低，最常见原因是物种/酶不匹配
   - `NoDecoyHits`：FDR 无法计算，需检查数据库 decoy 序列
   - `HighFdr`：通过 FDR 的 PSM 太少，结果统计不可靠
   - `NarrowTolerance`：容差过窄可能排除正确匹配
   - `WideTolerance`：容差过宽导致假阳性多
   - `DatabaseMismatch`：物种不匹配嫌疑
   - `SlowSearch`：匹配阶段瓶颈
   - `LowSpectraQuality`：碎片谱质量不足
3. **suggestions[]** → 展示优化建议

### 4. 搜索正常

如果无异常：
- 简要展示各阶段耗时
- 确认结果质量正常
- 建议下一步：`infer_proteins` 或 `export_results`

## 重试搜索

如果 suggestions 包含 param_changes：
1. 向用户展示原始参数和建议调整
2. **等待用户确认**后才能使用新参数调用 `run_search`
3. 保留原始 run_id 供结果对比

## 领域参考值

| 指标 | 正常范围 | 异常阈值 |
|------|---------|---------|
| DDA 鉴定率（HeLa） | 15-40% | < 10% |
| PSM @ 1% FDR | > 1000 | < 50 |
| 搜索耗时（Sage） | 1-5 min | > 10 min |
| 搜索耗时（SimpleSearch） | 5-30 min | > 60 min |
| 前体容差 | 5-20 ppm | < 5 或 > 50 ppm |
```

- [ ] **Step 2: Update proteomics-search.agent.md**

In `.github/agents/proteomics-search.agent.md`:

1. Add `diagnose_search` to the tools list in the YAML frontmatter
2. Add a diagnostic workflow section after the existing workflows:

```markdown
## 搜索诊断工作流

### 自动诊断
1. 搜索完成后，检查 `get_search_status` 返回的 `has_diagnostics` 和 `error_category`
2. 如果搜索失败（status 以 "Failed" 开头）→ 自动调用 `diagnose_search(run_id)`
3. 基于诊断数据解释失败原因，展示修复建议

### 质量评估
1. 搜索成功后，如果 `generate_summary` 显示鉴定率低于预期 → 建议调用 `diagnose_search`
2. 展示异常列表和优化建议
3. 询问用户是否要用调整后的参数重新搜索

### 决策边界
| 操作 | 自动/手动 |
|------|----------|
| 调用 diagnose_search | ✅ 可自动执行 |
| 解读诊断结果 | ✅ LLM 自动解读并展示 |
| 参数调整后重试 | ⚠️ 必须用户确认 |
```

- [ ] **Step 3: Commit**

```bash
git add .github/prompts/failure-diagnosis.prompt.md .github/agents/proteomics-search.agent.md
git commit -m "docs: add failure diagnosis prompt and update agent

New failure-diagnosis.prompt.md guides LLM interpretation of diagnostic
data. Agent updated with diagnose_search tool and diagnostic workflow.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 10: Final verification

**Files:** All modified files

- [ ] **Step 1: Full build**

Run: `cargo build --workspace 2>&1 | tail -5`
Expected: `Finished` with no errors

- [ ] **Step 2: Full test suite**

Run: `cargo test --workspace 2>&1 | grep -E "^test result:|FAILED"`
Expected: All pass, 0 failed

- [ ] **Step 3: Clippy**

Run: `cargo clippy --workspace 2>&1 | grep -E "warning|error" | head -10`
Expected: 0 warnings

- [ ] **Step 4: Verify tool count**

Run: `grep -c 'rmcp::tool(' crates/mcp-server/src/tools.rs`
Expected: 23 (was 22, now +1 diagnose_search)

- [ ] **Step 5: Verify diagnostics module exports**

Run: `grep 'pub mod diagnostics' crates/core/src/lib.rs`
Expected: `pub mod diagnostics;`

- [ ] **Step 6: Verify prompt count**

Run: `ls .github/prompts/*.prompt.md | wc -l`
Expected: 9 (was 8, now +1 failure-diagnosis)
