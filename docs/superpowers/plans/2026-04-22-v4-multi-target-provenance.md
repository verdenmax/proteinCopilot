# v4 Multi-Target Fragment Ion Provenance — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** For each L2/L3 trap PSM, find all co-eluting target peptides (light + heavy SILAC), match every observed fragment ion to specific target sources, and generate per-PSM HTML reports with mirror spectra and attribution tables.

**Architecture:** Extends entrapment-analysis crate with 3 new modules: `coelution.rs` (co-elution index build + query), `multi_provenance.rs` (multi-target fragment matching), `multi_report.rs` (per-PSM HTML rendering). Reuses v3's `generate_theoretical_ions()`, RT-based scan lookup, and IndexedMzMLReader. Config extends `ProvenanceConfig` with `SilacConfig`.

**Tech Stack:** Rust, serde, Plotly.js (embedded HTML), protein_copilot_spectrum_io (mzML reading), protein_copilot_core (MassTolerance)

**Spec:** `docs/superpowers/specs/2026-04-22-v4-multi-target-provenance-design.md`

---

### Task 1: Config — SilacConfig + ProvenanceConfig Extension

**Files:**
- Modify: `crates/entrapment-analysis/src/config.rs`
- Test: `crates/entrapment-analysis/src/config.rs` (inline `#[cfg(test)]` module)

- [ ] **Step 1: Write failing test for SilacConfig deserialization**

Add to the existing `#[cfg(test)]` module at the bottom of `config.rs`:

```rust
#[test]
fn test_silac_config_default() {
    let silac = SilacConfig::default();
    assert!((silac.heavy_k_delta - 8.014199).abs() < 1e-6);
    assert!((silac.heavy_r_delta - 10.008269).abs() < 1e-6);
    assert!(silac.enable_heavy_search);
}

#[test]
fn test_provenance_config_with_silac() {
    let yaml = r#"
fragment_tolerance_ppm: 20.0
max_fragment_charge: 2
silac:
  heavy_k_delta: 8.014199
  heavy_r_delta: 10.008269
  enable_heavy_search: true
generate_per_psm_reports: true
max_co_eluting_candidates: 15
"#;
    let config: ProvenanceConfig = serde_yaml::from_str(yaml).unwrap();
    assert!(config.silac.is_some());
    let silac = config.silac.unwrap();
    assert!((silac.heavy_k_delta - 8.014199).abs() < 1e-6);
    assert!(config.generate_per_psm_reports);
    assert_eq!(config.max_co_eluting_candidates, 15);
}

#[test]
fn test_provenance_config_backward_compat() {
    // v3 config without silac fields should still parse
    let yaml = r#"
fragment_tolerance_ppm: 20.0
max_fragment_charge: 2
chimera_threshold: 0.3
"#;
    let config: ProvenanceConfig = serde_yaml::from_str(yaml).unwrap();
    assert!(config.silac.is_none());
    assert!(config.generate_per_psm_reports);
    assert_eq!(config.max_co_eluting_candidates, 20);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p protein-copilot-entrapment-analysis test_silac_config -- --nocapture`
Expected: FAIL — `SilacConfig` not defined, `silac`/`generate_per_psm_reports`/`max_co_eluting_candidates` fields don't exist.

- [ ] **Step 3: Implement SilacConfig and extend ProvenanceConfig**

Add after the existing `ProvenanceConfig` struct (before `impl Default for ProvenanceConfig`):

```rust
/// SILAC heavy-label configuration for co-eluting target search.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SilacConfig {
    /// Delta mass for heavy Lysine (¹³C₆¹⁵N₂-Lys).
    #[serde(default = "default_heavy_k_delta")]
    pub heavy_k_delta: f64,

    /// Delta mass for heavy Arginine (¹³C₆¹⁵N₄-Arg).
    #[serde(default = "default_heavy_r_delta")]
    pub heavy_r_delta: f64,

    /// Whether to search for heavy-labeled co-eluting targets.
    #[serde(default = "default_true")]
    pub enable_heavy_search: bool,
}

impl Default for SilacConfig {
    fn default() -> Self {
        Self {
            heavy_k_delta: default_heavy_k_delta(),
            heavy_r_delta: default_heavy_r_delta(),
            enable_heavy_search: true,
        }
    }
}

fn default_heavy_k_delta() -> f64 {
    8.014199 // ¹³C₆¹⁵N₂-Lys
}
fn default_heavy_r_delta() -> f64 {
    10.008269 // ¹³C₆¹⁵N₄-Arg
}
```

Add 3 new fields to `ProvenanceConfig`:

```rust
    /// SILAC heavy-label configuration. If present, enables heavy co-eluting target search.
    #[serde(default)]
    pub silac: Option<SilacConfig>,

    /// Generate per-PSM HTML provenance reports for L2/L3 traps.
    #[serde(default = "default_true")]
    pub generate_per_psm_reports: bool,

    /// Maximum number of co-eluting candidates per trap PSM (prevents explosion in dense regions).
    #[serde(default = "default_max_co_eluting_candidates")]
    pub max_co_eluting_candidates: usize,
```

Add default function:

```rust
fn default_max_co_eluting_candidates() -> usize {
    20
}
```

Update `impl Default for ProvenanceConfig` to include the 3 new fields:

```rust
impl Default for ProvenanceConfig {
    fn default() -> Self {
        Self {
            fragment_tolerance_ppm: default_fragment_tolerance_ppm(),
            max_fragment_charge: default_max_fragment_charge(),
            chimera_threshold: default_chimera_threshold(),
            min_peaks_for_analysis: default_min_peaks_for_analysis(),
            levels_to_trace: default_levels_to_trace(),
            rt_tolerance_min: default_rt_tolerance_min(),
            silac: None,
            generate_per_psm_reports: true,
            max_co_eluting_candidates: default_max_co_eluting_candidates(),
        }
    }
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test -p protein-copilot-entrapment-analysis test_silac_config test_provenance_config -- --nocapture`
Expected: 3 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add crates/entrapment-analysis/src/config.rs
git commit -m "feat(entrapment): add SilacConfig and extend ProvenanceConfig for v4"
```

---

### Task 2: Core Types — CoElutingCandidate, LabelForm, MultiAnnotatedPeak, MultiTargetProvenance

**Files:**
- Modify: `crates/entrapment-analysis/src/types.rs`
- Test: `crates/entrapment-analysis/src/types.rs` (inline tests)

- [ ] **Step 1: Write failing test for new types**

Add to the `#[cfg(test)]` module in `types.rs`:

```rust
#[test]
fn test_label_form_light() {
    let form = LabelForm::Light;
    assert!(matches!(form, LabelForm::Light));
}

#[test]
fn test_label_form_heavy() {
    let form = LabelForm::Heavy {
        precursor_mz_heavy: 556.13,
        residue_deltas: vec![(9, 8.014199)], // K at position 9
    };
    if let LabelForm::Heavy { precursor_mz_heavy, residue_deltas } = form {
        assert!((precursor_mz_heavy - 556.13).abs() < 0.01);
        assert_eq!(residue_deltas.len(), 1);
    } else {
        panic!("expected Heavy");
    }
}

#[test]
fn test_co_eluting_candidate() {
    let candidate = CoElutingCandidate {
        peptide: "STTSGHLVYK".to_string(),
        protein_ids: vec!["sp|P12345|EF1A_HUMAN".to_string()],
        precursor_mz: 548.12,
        charge: 2,
        rt_start: 34.5,
        rt_stop: 35.8,
        label_form: LabelForm::Light,
        modifications: vec![],
    };
    assert_eq!(candidate.peptide, "STTSGHLVYK");
    assert_eq!(candidate.charge, 2);
}

#[test]
fn test_multi_annotated_peak() {
    let peak = MultiAnnotatedPeak {
        mz_observed: 285.155,
        intensity: 45230.0,
        trap_ion: Some("b3+1".to_string()),
        target_matches: vec![
            TargetIonMatch {
                candidate_index: 0,
                ion_label: "b3+1".to_string(),
                delta_ppm: -2.1,
            },
        ],
    };
    assert!(peak.trap_ion.is_some());
    assert_eq!(peak.target_matches.len(), 1);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p protein-copilot-entrapment-analysis test_label_form test_co_eluting test_multi_annotated -- --nocapture`
Expected: FAIL — types not defined.

- [ ] **Step 3: Implement the v4 types**

Add at the end of `types.rs` (before the `#[cfg(test)]` module):

```rust
// ---------------------------------------------------------------------------
// v4 Multi-Target Provenance Types
// ---------------------------------------------------------------------------

/// Whether a co-eluting candidate is a light or heavy (SILAC) form.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum LabelForm {
    /// Light (unlabeled) form.
    Light,
    /// Heavy (SILAC-labeled) form with shifted precursor and residue deltas.
    Heavy {
        /// Heavy precursor m/z.
        precursor_mz_heavy: f64,
        /// (0-based position, delta_Da) for each labeled residue (K or R).
        residue_deltas: Vec<(usize, f64)>,
    },
}

/// A co-eluting target peptide candidate found within the same DIA window and RT range.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoElutingCandidate {
    /// Target peptide sequence.
    pub peptide: String,
    /// Protein accession(s).
    pub protein_ids: Vec<String>,
    /// Precursor m/z (light form; for Heavy, use `label_form` field).
    pub precursor_mz: f64,
    /// Charge state.
    pub charge: i32,
    /// Elution window start (minutes).
    pub rt_start: f64,
    /// Elution window stop (minutes).
    pub rt_stop: f64,
    /// Light or Heavy form.
    pub label_form: LabelForm,
    /// Modifications: (0-based position, delta_mass_Da).
    pub modifications: Vec<(usize, f64)>,
}

/// A single target ion match for a multi-target annotated peak.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TargetIonMatch {
    /// Index into `MultiTargetProvenance::candidates`.
    pub candidate_index: usize,
    /// Ion label, e.g. "b3+1", "y5+2".
    pub ion_label: String,
    /// Matching error in ppm.
    pub delta_ppm: f64,
}

/// An observed peak annotated with multi-target provenance information.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiAnnotatedPeak {
    /// Observed m/z.
    pub mz_observed: f64,
    /// Observed intensity.
    pub intensity: f64,
    /// Ion label from the trap peptide (if matched).
    pub trap_ion: Option<String>,
    /// All target ion matches (may be 0, 1, or many).
    pub target_matches: Vec<TargetIonMatch>,
}

/// Complete multi-target provenance result for one trap PSM.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiTargetProvenance {
    /// The trap PSM being analyzed.
    pub trap_peptide: String,
    /// Scan number of the analyzed spectrum.
    pub scan_number: u32,
    /// All co-eluting target candidates (light + heavy).
    pub candidates: Vec<CoElutingCandidate>,
    /// Per-peak multi-target annotation.
    pub annotated_peaks: Vec<MultiAnnotatedPeak>,
    /// Summary: count of peaks matching only trap ions.
    pub trap_only_count: u32,
    /// Summary: count of peaks matching at least one target (not trap).
    pub target_only_count: u32,
    /// Summary: count of peaks matching both trap and at least one target.
    pub shared_count: u32,
    /// Summary: count of peaks matching nothing.
    pub unassigned_count: u32,
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test -p protein-copilot-entrapment-analysis test_label_form test_co_eluting test_multi_annotated -- --nocapture`
Expected: 4 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add crates/entrapment-analysis/src/types.rs
git commit -m "feat(entrapment): add v4 multi-target provenance types"
```

---

### Task 3: CoElutionIndex — Build and Query

**Files:**
- Create: `crates/entrapment-analysis/src/coelution.rs`
- Modify: `crates/entrapment-analysis/src/lib.rs` (add `pub mod coelution;`)

- [ ] **Step 1: Write failing test for CoElutionIndex**

Create `crates/entrapment-analysis/src/coelution.rs` with test module:

```rust
//! Co-elution index for finding target peptides that share an MS2 scan with a trap PSM.
//!
//! Builds a per-run, RT-sorted index of target PSMs from DIA-NN results.
//! Queries find all targets whose elution window overlaps the trap's RT range
//! and whose precursor m/z falls in the same DIA isolation window.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::config::SilacConfig;
use crate::types::{CoElutingCandidate, LabelForm, PsmGroup, UnifiedPsm};

/// An entry in the co-elution index (one per target PSM per run).
#[derive(Debug, Clone)]
struct TargetEntry {
    peptide: String,
    protein_ids: Vec<String>,
    precursor_mz: f64,
    charge: i32,
    rt_start: f64,
    rt_stop: f64,
    modifications: Vec<(usize, f64)>,
}

/// DIA isolation window definition.
#[derive(Debug, Clone, Copy)]
pub struct DiaWindow {
    pub center: f64,
    pub low: f64,
    pub high: f64,
}

/// Index of co-eluting target PSMs, organized by run name.
pub struct CoElutionIndex {
    by_run: HashMap<String, Vec<TargetEntry>>,
    dia_windows: Vec<DiaWindow>,
    silac: Option<SilacConfig>,
    max_candidates: usize,
}

// Implementation will go here...

#[cfg(test)]
mod tests {
    use super::*;

    fn make_target_psm(peptide: &str, mz: f64, rt_start: f64, rt_stop: f64, run: &str) -> (UnifiedPsm, PsmGroup) {
        let psm = UnifiedPsm {
            peptide: peptide.to_string(),
            charge: Some(2),
            precursor_mz: Some(mz),
            retention_time: Some((rt_start + rt_stop) / 2.0),
            rt_start: Some(rt_start),
            rt_stop: Some(rt_stop),
            scan_number: None,
            spectrum_file: Some(run.to_string()),
            protein_ids: "sp|P12345|TEST_HUMAN".to_string(),
            q_value: Some(0.001),
            modifications: vec![],
        };
        (psm, PsmGroup::Target)
    }

    fn make_trap_psm(peptide: &str, mz: f64, rt_start: f64, rt_stop: f64, run: &str) -> UnifiedPsm {
        UnifiedPsm {
            peptide: peptide.to_string(),
            charge: Some(2),
            precursor_mz: Some(mz),
            retention_time: Some((rt_start + rt_stop) / 2.0),
            rt_start: Some(rt_start),
            rt_stop: Some(rt_stop),
            scan_number: None,
            spectrum_file: Some(run.to_string()),
            protein_ids: "sp|P99999|TRAP_YEAST".to_string(),
            q_value: Some(0.005),
            modifications: vec![],
        }
    }

    #[test]
    fn test_build_index_and_query_basic() {
        // Two targets in same run, same DIA window, overlapping RT
        let targets = vec![
            make_target_psm("PEPTIDEA", 548.1, 34.5, 35.8, "Rep1"),
            make_target_psm("PEPTIDEB", 547.9, 35.0, 35.4, "Rep1"),
            make_target_psm("FARAWAY", 700.0, 50.0, 51.0, "Rep1"), // different RT
        ];
        let psms: Vec<UnifiedPsm> = targets.iter().map(|(p, _)| p.clone()).collect();
        let groups: Vec<PsmGroup> = targets.iter().map(|(_, g)| *g).collect();

        let windows = vec![DiaWindow { center: 548.0, low: 546.0, high: 550.0 }];
        let index = CoElutionIndex::build(&psms, &groups, &windows, None, 20);

        let trap = make_trap_psm("TRAPPEP", 548.3, 34.8, 35.6, "Rep1");
        let results = index.find_co_eluting(&trap, "Rep1");

        // Should find PEPTIDEA and PEPTIDEB (RT overlap + same window), not FARAWAY
        assert_eq!(results.len(), 2);
        let peptides: Vec<&str> = results.iter().map(|c| c.peptide.as_str()).collect();
        assert!(peptides.contains(&"PEPTIDEA"));
        assert!(peptides.contains(&"PEPTIDEB"));
    }

    #[test]
    fn test_no_rt_overlap() {
        let targets = vec![
            make_target_psm("PEPTIDEA", 548.1, 30.0, 31.0, "Rep1"), // RT too early
        ];
        let psms: Vec<UnifiedPsm> = targets.iter().map(|(p, _)| p.clone()).collect();
        let groups: Vec<PsmGroup> = targets.iter().map(|(_, g)| *g).collect();

        let windows = vec![DiaWindow { center: 548.0, low: 546.0, high: 550.0 }];
        let index = CoElutionIndex::build(&psms, &groups, &windows, None, 20);

        let trap = make_trap_psm("TRAPPEP", 548.3, 34.8, 35.6, "Rep1");
        let results = index.find_co_eluting(&trap, "Rep1");
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_different_dia_window() {
        let targets = vec![
            make_target_psm("PEPTIDEA", 600.0, 34.5, 35.8, "Rep1"), // different DIA window
        ];
        let psms: Vec<UnifiedPsm> = targets.iter().map(|(p, _)| p.clone()).collect();
        let groups: Vec<PsmGroup> = targets.iter().map(|(_, g)| *g).collect();

        let windows = vec![
            DiaWindow { center: 548.0, low: 546.0, high: 550.0 },
            DiaWindow { center: 600.0, low: 598.0, high: 602.0 },
        ];
        let index = CoElutionIndex::build(&psms, &groups, &windows, None, 20);

        let trap = make_trap_psm("TRAPPEP", 548.3, 34.8, 35.6, "Rep1");
        let results = index.find_co_eluting(&trap, "Rep1");
        assert_eq!(results.len(), 0); // different window
    }

    #[test]
    fn test_silac_heavy_pairing() {
        // Target with K at end → heavy should be shifted by 8.014/2 = 4.007 Da
        let targets = vec![
            make_target_psm("PEPTIDEK", 548.1, 34.5, 35.8, "Rep1"),
        ];
        let psms: Vec<UnifiedPsm> = targets.iter().map(|(p, _)| p.clone()).collect();
        let groups: Vec<PsmGroup> = targets.iter().map(|(_, g)| *g).collect();

        let silac = SilacConfig {
            heavy_k_delta: 8.014199,
            heavy_r_delta: 10.008269,
            enable_heavy_search: true,
        };

        // Window wide enough to include both light (548.1) and heavy (548.1 + 4.007 = 552.1)
        let windows = vec![
            DiaWindow { center: 548.0, low: 546.0, high: 550.0 },
            DiaWindow { center: 552.0, low: 550.0, high: 554.0 },
        ];
        let index = CoElutionIndex::build(&psms, &groups, &windows, Some(&silac), 20);

        let trap = make_trap_psm("TRAPPEP", 548.3, 34.8, 35.6, "Rep1");
        let results = index.find_co_eluting(&trap, "Rep1");

        // Should find light (same window) and heavy (different window)
        assert_eq!(results.len(), 2);
        let light = results.iter().find(|c| matches!(c.label_form, LabelForm::Light)).unwrap();
        assert_eq!(light.peptide, "PEPTIDEK");
        let heavy = results.iter().find(|c| matches!(c.label_form, LabelForm::Heavy { .. })).unwrap();
        assert_eq!(heavy.peptide, "PEPTIDEK");
    }

    #[test]
    fn test_max_candidates_cap() {
        // Create 25 overlapping targets, cap at 20
        let targets: Vec<(UnifiedPsm, PsmGroup)> = (0..25)
            .map(|i| make_target_psm(&format!("PEP{:02}", i), 548.0 + i as f64 * 0.01, 34.5, 35.8, "Rep1"))
            .collect();
        let psms: Vec<UnifiedPsm> = targets.iter().map(|(p, _)| p.clone()).collect();
        let groups: Vec<PsmGroup> = targets.iter().map(|(_, g)| *g).collect();

        let windows = vec![DiaWindow { center: 548.0, low: 546.0, high: 550.0 }];
        let index = CoElutionIndex::build(&psms, &groups, &windows, None, 20);

        let trap = make_trap_psm("TRAPPEP", 548.3, 34.8, 35.6, "Rep1");
        let results = index.find_co_eluting(&trap, "Rep1");
        assert!(results.len() <= 20);
    }
}
```

- [ ] **Step 2: Register module and run tests to verify they fail**

Add to `lib.rs` after the existing `pub mod provenance;` line:

```rust
pub mod coelution;
```

Run: `cargo test -p protein-copilot-entrapment-analysis coelution -- --nocapture`
Expected: FAIL — `CoElutionIndex::build` and `find_co_eluting` not implemented.

- [ ] **Step 3: Implement CoElutionIndex**

Add the implementation to `coelution.rs` between the struct definitions and the test module:

```rust
impl CoElutionIndex {
    /// Build from all PSMs + group assignments.
    ///
    /// Filters to Target PSMs with valid RT window + precursor_mz + spectrum_file.
    /// Sorts each run's entries by rt_start for efficient binary search.
    pub fn build(
        psms: &[UnifiedPsm],
        groups: &[PsmGroup],
        dia_windows: &[DiaWindow],
        silac: Option<&SilacConfig>,
        max_candidates: usize,
    ) -> Self {
        let mut by_run: HashMap<String, Vec<TargetEntry>> = HashMap::new();

        for (psm, group) in psms.iter().zip(groups.iter()) {
            if *group != PsmGroup::Target {
                continue;
            }
            let (Some(rt_start), Some(rt_stop)) = (psm.rt_start, psm.rt_stop) else {
                continue;
            };
            let Some(precursor_mz) = psm.precursor_mz else { continue };
            let Some(ref spectrum_file) = psm.spectrum_file else { continue };

            let entry = TargetEntry {
                peptide: psm.peptide.clone(),
                protein_ids: psm.protein_ids.split(';').map(|s| s.trim().to_string()).collect(),
                precursor_mz,
                charge: psm.charge.unwrap_or(2),
                rt_start,
                rt_stop,
                modifications: psm.modifications.clone(),
            };
            by_run.entry(spectrum_file.clone()).or_default().push(entry);
        }

        // Sort each run by rt_start for binary search
        for entries in by_run.values_mut() {
            entries.sort_by(|a, b| a.rt_start.partial_cmp(&b.rt_start).unwrap_or(std::cmp::Ordering::Equal));
        }

        Self {
            by_run,
            dia_windows: dia_windows.to_vec(),
            silac: silac.cloned(),
            max_candidates,
        }
    }

    /// Find all co-eluting targets for a trap PSM in the given run.
    ///
    /// Criteria:
    /// 1. RT windows overlap: target.[rt_start, rt_stop] ∩ trap.[rt_start, rt_stop] ≠ ∅
    /// 2. Same DIA isolation window: both precursor_mz values fall in the same window
    /// 3. For each light target, optionally generate a Heavy candidate if SILAC is configured
    pub fn find_co_eluting(&self, trap: &UnifiedPsm, run: &str) -> Vec<CoElutingCandidate> {
        let entries = match self.by_run.get(run) {
            Some(e) => e,
            None => return vec![],
        };

        let trap_rt_start = match trap.rt_start {
            Some(v) => v,
            None => return vec![],
        };
        let trap_rt_stop = match trap.rt_stop {
            Some(v) => v,
            None => return vec![],
        };
        let trap_mz = match trap.precursor_mz {
            Some(v) => v,
            None => return vec![],
        };

        let trap_window = self.find_dia_window(trap_mz);

        let mut candidates = Vec::new();

        // Binary search: find first entry where rt_start could overlap
        // An entry overlaps if entry.rt_stop >= trap_rt_start AND entry.rt_start <= trap_rt_stop
        // Since entries are sorted by rt_start, find starting position
        let start = entries.partition_point(|e| e.rt_stop < trap_rt_start);

        for entry in entries[start..].iter() {
            if entry.rt_start > trap_rt_stop {
                break; // past the overlap window
            }

            // Confirm RT overlap: entry.rt_stop >= trap_rt_start (guaranteed by start index)
            // AND entry.rt_start <= trap_rt_stop (guaranteed by break above)

            // Check same DIA window
            let entry_window = self.find_dia_window(entry.precursor_mz);
            let same_window = match (trap_window, entry_window) {
                (Some(tw), Some(ew)) => (tw.center - ew.center).abs() < 0.01,
                _ => false,
            };
            if !same_window {
                continue;
            }

            // Skip if it's the same peptide as the trap
            if entry.peptide == trap.peptide {
                continue;
            }

            // Add light candidate
            candidates.push(CoElutingCandidate {
                peptide: entry.peptide.clone(),
                protein_ids: entry.protein_ids.clone(),
                precursor_mz: entry.precursor_mz,
                charge: entry.charge,
                rt_start: entry.rt_start,
                rt_stop: entry.rt_stop,
                label_form: LabelForm::Light,
                modifications: entry.modifications.clone(),
            });

            // Add heavy candidate if SILAC configured
            if let Some(ref silac) = self.silac {
                if silac.enable_heavy_search {
                    let heavy_delta = compute_heavy_delta(&entry.peptide, silac);
                    if heavy_delta > 0.0 {
                        let heavy_mz = entry.precursor_mz + heavy_delta / entry.charge as f64;
                        let residue_deltas = compute_residue_deltas(&entry.peptide, silac);

                        candidates.push(CoElutingCandidate {
                            peptide: entry.peptide.clone(),
                            protein_ids: entry.protein_ids.clone(),
                            precursor_mz: entry.precursor_mz,
                            charge: entry.charge,
                            rt_start: entry.rt_start,
                            rt_stop: entry.rt_stop,
                            label_form: LabelForm::Heavy {
                                precursor_mz_heavy: heavy_mz,
                                residue_deltas,
                            },
                            modifications: entry.modifications.clone(),
                        });
                    }
                }
            }

            if candidates.len() >= self.max_candidates {
                break;
            }
        }

        candidates
    }

    /// Find which DIA window a precursor m/z falls into.
    fn find_dia_window(&self, mz: f64) -> Option<&DiaWindow> {
        self.dia_windows.iter().find(|w| mz >= w.low && mz <= w.high)
    }
}

/// Compute total heavy delta mass for a peptide sequence.
fn compute_heavy_delta(sequence: &str, silac: &SilacConfig) -> f64 {
    sequence.chars().fold(0.0, |acc, c| {
        acc + match c {
            'K' => silac.heavy_k_delta,
            'R' => silac.heavy_r_delta,
            _ => 0.0,
        }
    })
}

/// Compute per-residue heavy delta positions.
fn compute_residue_deltas(sequence: &str, silac: &SilacConfig) -> Vec<(usize, f64)> {
    sequence
        .chars()
        .enumerate()
        .filter_map(|(i, c)| match c {
            'K' => Some((i, silac.heavy_k_delta)),
            'R' => Some((i, silac.heavy_r_delta)),
            _ => None,
        })
        .collect()
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test -p protein-copilot-entrapment-analysis coelution -- --nocapture`
Expected: 5 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add crates/entrapment-analysis/src/coelution.rs crates/entrapment-analysis/src/lib.rs
git commit -m "feat(entrapment): implement CoElutionIndex for multi-target co-elution search"
```

---

### Task 4: Multi-Target Fragment Matching

**Files:**
- Create: `crates/entrapment-analysis/src/multi_provenance.rs`
- Modify: `crates/entrapment-analysis/src/lib.rs` (add `pub mod multi_provenance;`)
- Modify: `crates/entrapment-analysis/src/provenance.rs` (make `generate_theoretical_ions` and `amino_acid_mass` pub(crate))

- [ ] **Step 1: Make v3 helpers accessible**

In `provenance.rs`, change `generate_theoretical_ions` and `amino_acid_mass` from `fn` to `pub(crate) fn`, and `TheoreticalIon` from private to `pub(crate)`:

```rust
// Change:  struct TheoreticalIon {
// To:
pub(crate) struct TheoreticalIon {
```

```rust
// Change:  fn generate_theoretical_ions(
// To:
pub(crate) fn generate_theoretical_ions(
```

```rust
// Change:  fn amino_acid_mass(aa: char) -> f64 {
// To:
pub(crate) fn amino_acid_mass(aa: char) -> f64 {
```

- [ ] **Step 2: Write failing test for multi-target matching**

Create `crates/entrapment-analysis/src/multi_provenance.rs`:

```rust
//! Multi-target fragment ion matching.
//!
//! Extends v3's single-target provenance to match each observed MS2 peak
//! against theoretical ions from multiple co-eluting target peptides.

use protein_copilot_core::search_params::MassTolerance;

use crate::config::SilacConfig;
use crate::provenance::{generate_theoretical_ions, TheoreticalIon};
use crate::types::{
    CoElutingCandidate, LabelForm, MultiAnnotatedPeak, MultiTargetProvenance, TargetIonMatch,
};

// Implementation will go here...

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::search_params::ToleranceUnit;

    fn tolerance_20ppm() -> MassTolerance {
        MassTolerance { value: 20.0, unit: ToleranceUnit::Ppm }
    }

    #[test]
    fn test_match_single_target_light() {
        // Trap: ABCDE, Target_A: ABCDF (differ at position 4)
        // They share b1-b4 ions, differ in y-ions from y1
        let trap_seq = "STTTG";
        let target_seq = "STTSG"; // differ at pos 3: T→S

        let candidates = vec![CoElutingCandidate {
            peptide: target_seq.to_string(),
            protein_ids: vec!["P12345".to_string()],
            precursor_mz: 450.0,
            charge: 2,
            rt_start: 34.0,
            rt_stop: 36.0,
            label_form: LabelForm::Light,
            modifications: vec![],
        }];

        // Generate trap theoretical ions to use as "observed" peaks
        let trap_ions = generate_theoretical_ions(trap_seq, &[], 1);
        let observed_mz: Vec<f64> = trap_ions.iter().map(|i| i.mz).collect();
        let observed_int: Vec<f64> = vec![1000.0; observed_mz.len()];

        let result = trace_multi_target(
            &observed_mz,
            &observed_int,
            trap_seq,
            &[],
            &candidates,
            &tolerance_20ppm(),
            1,
        );

        assert_eq!(result.trap_peptide, "STTTG");
        assert_eq!(result.candidates.len(), 1);
        // All peaks should match trap; some should also match target (shared)
        assert!(result.shared_count > 0 || result.trap_only_count > 0);
        assert_eq!(result.annotated_peaks.len(), observed_mz.len());
    }

    #[test]
    fn test_match_multiple_targets() {
        let trap_seq = "STTTGHLIYK";
        let candidates = vec![
            CoElutingCandidate {
                peptide: "STTSGHLVYK".to_string(),
                protein_ids: vec!["P12345".to_string()],
                precursor_mz: 548.1,
                charge: 2,
                rt_start: 34.5,
                rt_stop: 35.8,
                label_form: LabelForm::Light,
                modifications: vec![],
            },
            CoElutingCandidate {
                peptide: "AETFGHLK".to_string(),
                protein_ids: vec!["P67890".to_string()],
                precursor_mz: 547.9,
                charge: 2,
                rt_start: 35.0,
                rt_stop: 35.4,
                label_form: LabelForm::Light,
                modifications: vec![],
            },
        ];

        let trap_ions = generate_theoretical_ions(trap_seq, &[], 1);
        let observed_mz: Vec<f64> = trap_ions.iter().map(|i| i.mz).collect();
        let observed_int: Vec<f64> = vec![1000.0; observed_mz.len()];

        let result = trace_multi_target(
            &observed_mz,
            &observed_int,
            trap_seq,
            &[],
            &candidates,
            &tolerance_20ppm(),
            1,
        );

        assert_eq!(result.candidates.len(), 2);
        assert_eq!(result.annotated_peaks.len(), observed_mz.len());
        // Every peak should have trap_ion set (since observed = trap theoretical)
        for peak in &result.annotated_peaks {
            assert!(peak.trap_ion.is_some());
        }
    }

    #[test]
    fn test_heavy_label_ion_generation() {
        let candidates = vec![CoElutingCandidate {
            peptide: "PEPTIDEK".to_string(),
            protein_ids: vec!["P11111".to_string()],
            precursor_mz: 450.0,
            charge: 2,
            rt_start: 34.0,
            rt_stop: 36.0,
            label_form: LabelForm::Heavy {
                precursor_mz_heavy: 454.007,
                residue_deltas: vec![(7, 8.014199)], // K at position 7
            },
            modifications: vec![],
        }];

        // Use a known m/z that should match the heavy y1 ion (K + heavy delta)
        // Light y1 for K: (128.09496 + 18.010565 + 1.007276) / 1 = 147.113
        // Heavy y1 for K: (128.09496 + 8.014199 + 18.010565 + 1.007276) / 1 = 155.127
        let observed_mz = vec![155.127];
        let observed_int = vec![1000.0];

        let result = trace_multi_target(
            &observed_mz,
            &observed_int,
            "ABCDEFGH", // trap (won't match 155.127)
            &[],
            &candidates,
            &MassTolerance { value: 50.0, unit: ToleranceUnit::Ppm }, // wider tolerance for test
            1,
        );

        // The heavy y1 ion should match
        assert_eq!(result.target_only_count, 1);
        assert_eq!(result.annotated_peaks[0].target_matches.len(), 1);
    }
}
```

- [ ] **Step 3: Register module and run tests to verify they fail**

Add to `lib.rs`:

```rust
pub mod multi_provenance;
```

Run: `cargo test -p protein-copilot-entrapment-analysis multi_provenance -- --nocapture`
Expected: FAIL — `trace_multi_target` not defined.

- [ ] **Step 4: Implement trace_multi_target**

Add implementation to `multi_provenance.rs` (before the test module):

```rust
/// Trace fragment ion provenance against multiple co-eluting targets.
///
/// For each observed peak, matches against:
/// 1. Trap peptide theoretical b/y ions
/// 2. Each candidate's theoretical b/y ions (with heavy-label shifts if applicable)
///
/// Returns per-peak attribution with all matching targets identified.
pub fn trace_multi_target(
    observed_mz: &[f64],
    observed_intensity: &[f64],
    trap_sequence: &str,
    trap_modifications: &[(usize, f64)],
    candidates: &[CoElutingCandidate],
    fragment_tolerance: &MassTolerance,
    max_fragment_charge: i32,
) -> MultiTargetProvenance {
    use protein_copilot_search_engine::matching::within_tolerance;

    // 1. Generate trap theoretical ions
    let trap_ions = generate_theoretical_ions(trap_sequence, trap_modifications, max_fragment_charge);

    // 2. Generate theoretical ions for each candidate
    let candidate_ions: Vec<Vec<TheoreticalIon>> = candidates
        .iter()
        .map(|c| generate_candidate_ions(c, max_fragment_charge))
        .collect();

    // 3. Match each observed peak
    let mut annotated_peaks = Vec::with_capacity(observed_mz.len());
    let mut trap_only = 0u32;
    let mut target_only = 0u32;
    let mut shared = 0u32;
    let mut unassigned = 0u32;

    for (i, &mz) in observed_mz.iter().enumerate() {
        let intensity = observed_intensity.get(i).copied().unwrap_or(0.0);

        // Match against trap ions
        let trap_ion = find_best_match(mz, &trap_ions, fragment_tolerance);

        // Match against each candidate's ions
        let mut target_matches = Vec::new();
        for (ci, ions) in candidate_ions.iter().enumerate() {
            if let Some((label, ppm)) = find_best_match_with_ppm(mz, ions, fragment_tolerance) {
                target_matches.push(TargetIonMatch {
                    candidate_index: ci,
                    ion_label: label,
                    delta_ppm: ppm,
                });
            }
        }

        // Classify
        match (trap_ion.is_some(), !target_matches.is_empty()) {
            (true, true) => shared += 1,
            (true, false) => trap_only += 1,
            (false, true) => target_only += 1,
            (false, false) => unassigned += 1,
        }

        annotated_peaks.push(MultiAnnotatedPeak {
            mz_observed: mz,
            intensity,
            trap_ion: trap_ion.map(|(l, _)| l),
            target_matches,
        });
    }

    MultiTargetProvenance {
        trap_peptide: trap_sequence.to_string(),
        scan_number: 0, // set by caller
        candidates: candidates.to_vec(),
        annotated_peaks,
        trap_only_count: trap_only,
        target_only_count: target_only,
        shared_count: shared,
        unassigned_count: unassigned,
    }
}

/// Generate theoretical ions for a candidate, applying heavy-label shifts if applicable.
fn generate_candidate_ions(
    candidate: &CoElutingCandidate,
    max_charge: i32,
) -> Vec<TheoreticalIon> {
    let light_ions = generate_theoretical_ions(
        &candidate.peptide,
        &candidate.modifications,
        max_charge,
    );

    match &candidate.label_form {
        LabelForm::Light => light_ions,
        LabelForm::Heavy { residue_deltas, .. } => {
            // Shift each ion's m/z based on heavy residue deltas
            shift_ions_heavy(&candidate.peptide, &light_ions, residue_deltas, max_charge)
        }
    }
}

/// Shift theoretical ions for heavy-labeled peptide.
///
/// For b_n ions: add cumulative heavy delta of residues [0..n]
/// For y_n ions: add cumulative heavy delta of residues [len-n..len]
fn shift_ions_heavy(
    sequence: &str,
    light_ions: &[TheoreticalIon],
    residue_deltas: &[(usize, f64)],
    _max_charge: i32,
) -> Vec<TheoreticalIon> {
    let n = sequence.len();
    // Build a map: position → delta
    let delta_map: std::collections::HashMap<usize, f64> = residue_deltas.iter().cloned().collect();

    light_ions
        .iter()
        .map(|ion| {
            // Parse ion label to determine type and number
            let (ion_type, ion_num, charge) = parse_ion_label(&ion.label);
            let cumulative_delta: f64 = match ion_type {
                'b' => {
                    // b_n covers residues 0..n
                    (0..ion_num).filter_map(|pos| delta_map.get(&pos)).sum()
                }
                'y' => {
                    // y_n covers residues (n - ion_num)..n
                    ((n - ion_num)..n).filter_map(|pos| delta_map.get(&pos)).sum()
                }
                _ => 0.0,
            };
            TheoreticalIon {
                mz: ion.mz + cumulative_delta / charge as f64,
                label: format!("{}(H)", ion.label), // mark as heavy
            }
        })
        .collect()
}

/// Parse an ion label like "b3+1" into (type='b', number=3, charge=1).
fn parse_ion_label(label: &str) -> (char, usize, i32) {
    let ion_type = label.chars().next().unwrap_or('?');
    let rest = &label[1..];
    let parts: Vec<&str> = rest.split('+').collect();
    let ion_num: usize = parts[0].parse().unwrap_or(0);
    let charge: i32 = parts.get(1).and_then(|s| s.parse().ok()).unwrap_or(1);
    (ion_type, ion_num, charge)
}

/// Find best matching ion, returning (label, delta_ppm).
fn find_best_match_with_ppm(
    observed_mz: f64,
    ions: &[TheoreticalIon],
    tolerance: &MassTolerance,
) -> Option<(String, f64)> {
    use protein_copilot_search_engine::matching::within_tolerance;

    let mut best: Option<(f64, String, f64)> = None;
    for ion in ions {
        if within_tolerance(observed_mz, ion.mz, tolerance) {
            let error = (observed_mz - ion.mz).abs();
            let ppm = (observed_mz - ion.mz) / ion.mz * 1e6;
            match best {
                Some((best_err, _, _)) if error < best_err => {
                    best = Some((error, ion.label.clone(), ppm));
                }
                None => {
                    best = Some((error, ion.label.clone(), ppm));
                }
                _ => {}
            }
        }
    }
    best.map(|(_, label, ppm)| (label, ppm))
}

/// Find best matching ion, returning (label, delta_ppm) tuple.
fn find_best_match(
    observed_mz: f64,
    ions: &[TheoreticalIon],
    tolerance: &MassTolerance,
) -> Option<(String, f64)> {
    find_best_match_with_ppm(observed_mz, ions, tolerance)
}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `cargo test -p protein-copilot-entrapment-analysis multi_provenance -- --nocapture`
Expected: 3 tests PASS.

- [ ] **Step 6: Commit**

```bash
git add crates/entrapment-analysis/src/multi_provenance.rs crates/entrapment-analysis/src/provenance.rs crates/entrapment-analysis/src/lib.rs
git commit -m "feat(entrapment): implement multi-target fragment ion matching"
```

---

### Task 5: Per-PSM HTML Report Renderer

**Files:**
- Create: `crates/entrapment-analysis/src/multi_report.rs`
- Modify: `crates/entrapment-analysis/src/lib.rs` (add `pub mod multi_report;`)

- [ ] **Step 1: Write failing test**

Create `crates/entrapment-analysis/src/multi_report.rs` with test:

```rust
//! Per-PSM HTML report renderer for multi-target provenance.
//!
//! Generates a self-contained HTML file with:
//! 1. Co-eluting target candidate table
//! 2. Mirror spectrum (Plotly.js) — trap peaks up, target peaks down
//! 3. Fragment ion attribution table

use crate::types::{
    CoElutingCandidate, LabelForm, MultiAnnotatedPeak, MultiTargetProvenance, TargetIonMatch,
};

// Implementation will go here...

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_provenance() -> MultiTargetProvenance {
        MultiTargetProvenance {
            trap_peptide: "STTTGHLIYK".to_string(),
            scan_number: 12345,
            candidates: vec![
                CoElutingCandidate {
                    peptide: "STTSGHLVYK".to_string(),
                    protein_ids: vec!["sp|P12345|EF1A_HUMAN".to_string()],
                    precursor_mz: 548.12,
                    charge: 2,
                    rt_start: 34.5,
                    rt_stop: 35.8,
                    label_form: LabelForm::Light,
                    modifications: vec![],
                },
            ],
            annotated_peaks: vec![
                MultiAnnotatedPeak {
                    mz_observed: 285.155,
                    intensity: 45230.0,
                    trap_ion: Some("b3+1".to_string()),
                    target_matches: vec![TargetIonMatch {
                        candidate_index: 0,
                        ion_label: "b3+1".to_string(),
                        delta_ppm: -2.1,
                    }],
                },
                MultiAnnotatedPeak {
                    mz_observed: 386.203,
                    intensity: 72100.0,
                    trap_ion: Some("b4+1".to_string()),
                    target_matches: vec![],
                },
                MultiAnnotatedPeak {
                    mz_observed: 512.334,
                    intensity: 8200.0,
                    trap_ion: None,
                    target_matches: vec![],
                },
            ],
            trap_only_count: 1,
            target_only_count: 0,
            shared_count: 1,
            unassigned_count: 1,
        }
    }

    #[test]
    fn test_generate_html_contains_sections() {
        let prov = make_test_provenance();
        let html = generate_multi_provenance_html(&prov);

        assert!(html.contains("Fragment Ion Provenance Report"));
        assert!(html.contains("STTTGHLIYK"));
        assert!(html.contains("STTSGHLVYK"));
        assert!(html.contains("285.155")); // observed m/z in table
        assert!(html.contains("plotly")); // Plotly.js reference
        assert!(html.contains("b3+1")); // ion label
        assert!(html.contains("TrapOnly")); // origin classification
        assert!(html.contains("Shared")); // origin classification
    }

    #[test]
    fn test_write_html_file() {
        let prov = make_test_provenance();
        let dir = std::env::temp_dir().join("test_multi_report");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test_report.html");

        render_multi_provenance_report(&prov, &path).unwrap();

        assert!(path.exists());
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("<!DOCTYPE html>"));
        assert!(content.contains("STTTGHLIYK"));

        // Cleanup
        std::fs::remove_dir_all(&dir).ok();
    }
}
```

- [ ] **Step 2: Register module and run tests to verify they fail**

Add to `lib.rs`:

```rust
pub mod multi_report;
```

Run: `cargo test -p protein-copilot-entrapment-analysis multi_report -- --nocapture`
Expected: FAIL — functions not defined.

- [ ] **Step 3: Implement HTML renderer**

Add implementation to `multi_report.rs` (before test module):

```rust
use std::fmt::Write;
use std::path::Path;

/// Color palette for candidates (up to 10 distinct colors).
const CANDIDATE_COLORS: &[&str] = &[
    "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b",
    "#e377c2", "#bcbd22", "#17becf", "#1f77b4", "#aec7e8",
];

const TRAP_COLOR: &str = "#1f77b4";
const SHARED_COLOR: &str = "#9467bd";
const UNASSIGNED_COLOR: &str = "#7f7f7f";

/// Generate complete HTML string for a multi-target provenance report.
pub fn generate_multi_provenance_html(prov: &MultiTargetProvenance) -> String {
    let mut html = String::with_capacity(16_000);

    // Header
    write_header(&mut html, prov);

    // Section 1: Candidate table
    write_candidate_table(&mut html, prov);

    // Section 2: Mirror spectrum (Plotly.js)
    write_mirror_spectrum(&mut html, prov);

    // Section 3: Attribution table
    write_attribution_table(&mut html, prov);

    // Footer summary
    write_footer(&mut html, prov);

    // Close HTML
    html.push_str("</body></html>");
    html
}

/// Write the report to an HTML file.
pub fn render_multi_provenance_report(
    prov: &MultiTargetProvenance,
    output_path: &Path,
) -> Result<(), std::io::Error> {
    let html = generate_multi_provenance_html(prov);
    std::fs::write(output_path, html)
}

fn write_header(html: &mut String, prov: &MultiTargetProvenance) {
    let total_peaks = prov.annotated_peaks.len();
    write!(html, r#"<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Provenance: {trap}</title>
<script src="https://cdn.plot.ly/plotly-2.35.0.min.js"></script>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 0; color: #333; }}
.header {{ padding: 15px 20px; background: #2c3e50; color: white; }}
.header h1 {{ margin: 0; font-size: 18px; }}
.header .meta {{ font-size: 13px; margin-top: 5px; opacity: 0.85; }}
.section {{ padding: 15px 20px; border-bottom: 1px solid #ddd; }}
.section h2 {{ font-size: 15px; color: #2c3e50; margin: 0 0 10px 0; }}
table {{ width: 100%; border-collapse: collapse; font-size: 12px; }}
th {{ padding: 6px 8px; text-align: left; background: #ecf0f1; }}
td {{ padding: 5px 8px; }}
tr:nth-child(even) {{ background: #f9f9f9; }}
.footer {{ padding: 12px 20px; background: #ecf0f1; font-size: 12px; }}
.color-dot {{ display: inline-block; width: 10px; height: 10px; border-radius: 2px; margin-right: 4px; vertical-align: middle; }}
.origin-trap {{ color: {trap_color}; font-weight: bold; }}
.origin-shared {{ color: {shared_color}; font-weight: bold; }}
.origin-target {{ color: #d62728; font-weight: bold; }}
.origin-unassigned {{ color: {unassigned_color}; }}
</style>
</head><body>
<div class="header">
<h1>🔬 Fragment Ion Provenance Report</h1>
<div class="meta">Trap: <b>{trap}</b> | Scan: {scan} | Candidates: {ncand} | Peaks: {npeaks}</div>
</div>
"#,
        trap = prov.trap_peptide,
        scan = prov.scan_number,
        ncand = prov.candidates.len(),
        npeaks = total_peaks,
        trap_color = TRAP_COLOR,
        shared_color = SHARED_COLOR,
        unassigned_color = UNASSIGNED_COLOR,
    ).ok();
}

fn write_candidate_table(html: &mut String, prov: &MultiTargetProvenance) {
    html.push_str(r#"<div class="section"><h2>📋 Co-eluting Target Candidates</h2><table>
<tr><th>#</th><th>Peptide</th><th>Protein</th><th>Label</th><th>m/z</th><th>z</th><th>RT range</th><th>Matched ions</th></tr>"#);

    for (i, c) in prov.candidates.iter().enumerate() {
        let color = CANDIDATE_COLORS.get(i).unwrap_or(&"#333");
        let label_str = match &c.label_form {
            LabelForm::Light => "🔵 Light".to_string(),
            LabelForm::Heavy { precursor_mz_heavy, .. } => format!("🔴 Heavy ({:.1})", precursor_mz_heavy),
        };
        let matched = prov.annotated_peaks.iter()
            .filter(|p| p.target_matches.iter().any(|m| m.candidate_index == i))
            .count();
        write!(html, r#"<tr><td><span class="color-dot" style="background:{color}"></span></td>
<td style="font-family:monospace"><b>{pep}</b></td>
<td style="font-size:11px">{prot}</td>
<td>{label}</td>
<td>{mz:.2}</td><td>{z}</td>
<td>{rt_s:.1} – {rt_e:.1}</td>
<td style="font-weight:bold;color:{color}">{matched}</td></tr>"#,
            color = color,
            pep = c.peptide,
            prot = c.protein_ids.first().map(|s| s.as_str()).unwrap_or(""),
            label = label_str,
            mz = c.precursor_mz,
            z = c.charge,
            rt_s = c.rt_start,
            rt_e = c.rt_stop,
            matched = matched,
        ).ok();
    }
    html.push_str("</table></div>\n");
}

fn write_mirror_spectrum(html: &mut String, prov: &MultiTargetProvenance) {
    html.push_str(r#"<div class="section"><h2>📊 Mirror Spectrum (Trap ↑ vs Targets ↓)</h2>
<div id="mirror-plot" style="width:100%;height:400px;"></div>
<script>
"#);

    // Build trap traces (positive y)
    let mut trap_mz = Vec::new();
    let mut trap_int = Vec::new();
    let mut trap_text = Vec::new();
    let mut trap_colors = Vec::new();

    // Build target traces (negative y) grouped by candidate
    let mut target_mz = Vec::new();
    let mut target_int = Vec::new();
    let mut target_text = Vec::new();
    let mut target_colors = Vec::new();

    for peak in &prov.annotated_peaks {
        let has_trap = peak.trap_ion.is_some();
        let has_target = !peak.target_matches.is_empty();

        if has_trap {
            trap_mz.push(peak.mz_observed);
            trap_int.push(peak.intensity);
            let label = peak.trap_ion.as_deref().unwrap_or("");
            trap_text.push(label.to_string());
            if has_target {
                trap_colors.push(SHARED_COLOR.to_string());
            } else {
                trap_colors.push(TRAP_COLOR.to_string());
            }
        }

        for tm in &peak.target_matches {
            target_mz.push(peak.mz_observed);
            target_int.push(-peak.intensity); // negative for mirror
            target_text.push(format!("{} ({})", tm.ion_label, prov.candidates.get(tm.candidate_index).map(|c| c.peptide.as_str()).unwrap_or("?")));
            let color = CANDIDATE_COLORS.get(tm.candidate_index).unwrap_or(&"#333");
            target_colors.push(color.to_string());
        }
    }

    // Write Plotly traces
    write!(html, "var trapTrace = {{x:{trap_mz:?},y:{trap_int:?},text:{trap_text:?},type:'bar',name:'Trap',marker:{{color:{trap_colors:?}}},hovertemplate:'%{{text}}<br>m/z: %{{x:.3f}}<br>Int: %{{y:.0f}}'}};\n").ok();
    write!(html, "var targetTrace = {{x:{target_mz:?},y:{target_int:?},text:{target_text:?},type:'bar',name:'Targets',marker:{{color:{target_colors:?}}},hovertemplate:'%{{text}}<br>m/z: %{{x:.3f}}<br>Int: %{{y:.0f}}'}};\n").ok();

    write!(html, r#"var layout = {{
title: 'Trap: {trap}',
xaxis: {{title: 'm/z'}},
yaxis: {{title: 'Intensity', zeroline: true}},
barmode: 'overlay',
bargap: 0,
hovermode: 'closest'
}};
Plotly.newPlot('mirror-plot', [trapTrace, targetTrace], layout);
</script></div>
"#, trap = prov.trap_peptide).ok();
}

fn write_attribution_table(html: &mut String, prov: &MultiTargetProvenance) {
    html.push_str(r#"<div class="section"><h2>📝 Fragment Ion Attribution Table</h2><table>
<tr><th>Obs. m/z</th><th>Intensity</th><th>Trap Ion</th><th>Origin</th><th>Target Matches</th><th>Δppm</th></tr>"#);

    for peak in &prov.annotated_peaks {
        let has_trap = peak.trap_ion.is_some();
        let has_target = !peak.target_matches.is_empty();

        let origin = match (has_trap, has_target) {
            (true, true) => format!(r#"<span class="origin-shared">Shared×{}</span>"#, peak.target_matches.len()),
            (true, false) => r#"<span class="origin-trap">TrapOnly</span>"#.to_string(),
            (false, true) => r#"<span class="origin-target">TargetOnly</span>"#.to_string(),
            (false, false) => r#"<span class="origin-unassigned">Unassigned</span>"#.to_string(),
        };

        let trap_label = peak.trap_ion.as_deref().unwrap_or("—");

        let mut targets_html = String::new();
        let mut ppm_html = String::new();
        if peak.target_matches.is_empty() {
            targets_html.push_str("—");
            ppm_html.push_str("—");
        } else {
            for (j, tm) in peak.target_matches.iter().enumerate() {
                if j > 0 { targets_html.push_str("<br>"); ppm_html.push_str("<br>"); }
                let color = CANDIDATE_COLORS.get(tm.candidate_index).unwrap_or(&"#333");
                let pep = prov.candidates.get(tm.candidate_index).map(|c| c.peptide.as_str()).unwrap_or("?");
                write!(targets_html, r#"<span class="color-dot" style="background:{}"></span>{} {}"#, color, pep, tm.ion_label).ok();
                write!(ppm_html, "{:+.1}", tm.delta_ppm).ok();
            }
        }

        write!(html, "<tr><td>{:.3}</td><td>{:.0}</td><td style=\"font-family:monospace\">{}</td><td>{}</td><td>{}</td><td>{}</td></tr>\n",
            peak.mz_observed, peak.intensity, trap_label, origin, targets_html, ppm_html
        ).ok();
    }
    html.push_str("</table></div>\n");
}

fn write_footer(html: &mut String, prov: &MultiTargetProvenance) {
    let total = prov.annotated_peaks.len();
    write!(html, r#"<div class="footer"><b>Summary:</b> {} peaks | {} TrapOnly | {} Shared | {} TargetOnly | {} Unassigned | {} candidates</div>"#,
        total, prov.trap_only_count, prov.shared_count, prov.target_only_count, prov.unassigned_count, prov.candidates.len()
    ).ok();
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test -p protein-copilot-entrapment-analysis multi_report -- --nocapture`
Expected: 2 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add crates/entrapment-analysis/src/multi_report.rs crates/entrapment-analysis/src/lib.rs
git commit -m "feat(entrapment): implement per-PSM HTML provenance report renderer"
```

---

### Task 6: Batch Pipeline Integration

**Files:**
- Modify: `crates/entrapment-analysis/src/lib.rs`

- [ ] **Step 1: Write failing test**

Add to the integration tests file `crates/entrapment-analysis/tests/v3_e2e_provenance.rs` (or create a new `v4_multi_target.rs`):

Create `crates/entrapment-analysis/tests/v4_multi_target.rs`:

```rust
//! v4 multi-target provenance integration tests.

use protein_copilot_entrapment_analysis::coelution::{CoElutionIndex, DiaWindow};
use protein_copilot_entrapment_analysis::multi_provenance::trace_multi_target;
use protein_copilot_entrapment_analysis::types::*;

use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};

#[test]
fn test_full_pipeline_mock() {
    // 1. Create mock target PSMs with RT windows
    let targets: Vec<(UnifiedPsm, PsmGroup)> = vec![
        (UnifiedPsm {
            peptide: "STTSGHLVYK".to_string(),
            charge: Some(2),
            precursor_mz: Some(548.12),
            retention_time: Some(35.15),
            rt_start: Some(34.5),
            rt_stop: Some(35.8),
            scan_number: None,
            spectrum_file: Some("Rep1".to_string()),
            protein_ids: "sp|P12345|EF1A_HUMAN".to_string(),
            q_value: Some(0.001),
            modifications: vec![],
        }, PsmGroup::Target),
    ];

    let psms: Vec<UnifiedPsm> = targets.iter().map(|(p, _)| p.clone()).collect();
    let groups: Vec<PsmGroup> = targets.iter().map(|(_, g)| *g).collect();

    // 2. Build index
    let windows = vec![DiaWindow { center: 548.0, low: 546.0, high: 550.0 }];
    let index = CoElutionIndex::build(&psms, &groups, &windows, None, 20);

    // 3. Query for a trap PSM
    let trap = UnifiedPsm {
        peptide: "STTTGHLIYK".to_string(),
        charge: Some(2),
        precursor_mz: Some(548.30),
        retention_time: Some(35.2),
        rt_start: Some(34.8),
        rt_stop: Some(35.6),
        scan_number: Some(12345),
        spectrum_file: Some("Rep1".to_string()),
        protein_ids: "sp|P99999|TRAP_YEAST".to_string(),
        q_value: Some(0.005),
        modifications: vec![],
    };

    let candidates = index.find_co_eluting(&trap, "Rep1");
    assert!(!candidates.is_empty(), "should find co-eluting targets");
    assert_eq!(candidates[0].peptide, "STTSGHLVYK");

    // 4. Run multi-target matching with synthetic peaks
    let tolerance = MassTolerance { value: 20.0, unit: ToleranceUnit::Ppm };
    // Use a few known m/z values that should match
    let observed_mz = vec![200.0, 300.0, 400.0, 500.0, 600.0, 700.0];
    let observed_int = vec![1000.0; 6];

    let result = trace_multi_target(
        &observed_mz,
        &observed_int,
        &trap.peptide,
        &trap.modifications,
        &candidates,
        &tolerance,
        2,
    );

    assert_eq!(result.trap_peptide, "STTTGHLIYK");
    assert_eq!(result.annotated_peaks.len(), 6);
}
```

- [ ] **Step 2: Run test to verify it passes (basic integration)**

Run: `cargo test -p protein-copilot-entrapment-analysis v4_multi_target -- --nocapture`
Expected: PASS — this uses already-implemented components.

- [ ] **Step 3: Implement `trace_multi_target_provenance` batch function in lib.rs**

Add after the existing `trace_provenance_batch` function in `lib.rs`:

```rust
// ---------------------------------------------------------------------------
// v4 Multi-target provenance batch tracing
// ---------------------------------------------------------------------------

/// Extract DIA isolation windows from an mzML file via its spectrum index.
pub fn extract_dia_windows(
    reader: &dyn protein_copilot_spectrum_io::reader::SpectrumReader,
    mzml_path: &std::path::Path,
) -> Vec<coelution::DiaWindow> {
    let mut seen = std::collections::HashSet::new();
    let mut windows = Vec::new();

    if let Ok(index) = reader.get_index(mzml_path) {
        for meta in index.entries() {
            if meta.ms_level == 2 {
                if let Some((center, lower, upper)) = meta.isolation_window {
                    let key = ((center * 100.0) as i64, (lower * 100.0) as i64, (upper * 100.0) as i64);
                    if seen.insert(key) {
                        windows.push(coelution::DiaWindow {
                            center,
                            low: center - lower,
                            high: center + upper,
                        });
                    }
                }
            }
        }
    }

    windows
}

/// Run multi-target provenance tracing on classified PSMs.
///
/// For each L2/L3 trap PSM:
/// 1. Find co-eluting targets from the full PSM list
/// 2. Read the MS2 spectrum from mzML
/// 3. Match observed peaks against all candidate theoretical ions
/// 4. Generate per-PSM HTML report
///
/// Returns the number of PSMs successfully traced and the provenance results.
pub fn trace_multi_target_provenance(
    classified: &[ClassifiedPsm],
    all_psms: &[UnifiedPsm],
    all_groups: &[PsmGroup],
    mzml_dir: &Path,
    config: &EntrapmentConfig,
    output_dir: &Path,
) -> Result<(u32, Vec<types::MultiTargetProvenance>), EntrapmentError> {
    use std::collections::HashSet;
    use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};

    let tolerance = MassTolerance {
        value: config.provenance.fragment_tolerance_ppm,
        unit: ToleranceUnit::Ppm,
    };

    let levels_to_trace: HashSet<&str> = config
        .provenance
        .levels_to_trace
        .iter()
        .map(|s| s.as_str())
        .collect();

    // Collect eligible trap PSMs
    let mut eligible: Vec<(usize, String)> = Vec::new(); // (index, run_name)
    for (idx, cpsm) in classified.iter().enumerate() {
        if cpsm.group != PsmGroup::Trap { continue; }
        if !levels_to_trace.contains(cpsm.level.as_str()) { continue; }

        let has_scan = cpsm.psm.scan_number.map_or(false, |s| s > 0);
        let has_rt_mz = cpsm.psm.retention_time.is_some() && cpsm.psm.precursor_mz.is_some();
        if !has_scan && !has_rt_mz { continue; }

        if let Some(ref run) = cpsm.psm.spectrum_file {
            eligible.push((idx, run.clone()));
        }
    }

    if eligible.is_empty() {
        return Ok((0, vec![]));
    }

    // Get unique runs
    let runs: HashSet<String> = eligible.iter().map(|(_, r)| r.clone()).collect();

    // For each run: build reader, extract DIA windows, build index, trace
    let mut results = Vec::new();
    let mut traced_count = 0u32;

    // Create provenance output directory
    let prov_dir = output_dir.join("provenance");
    if config.provenance.generate_per_psm_reports {
        std::fs::create_dir_all(&prov_dir).map_err(|e| {
            EntrapmentError::IoError(format!("failed to create provenance dir: {}", e))
        })?;
    }

    for run in &runs {
        let mzml_path = match find_mzml_file(mzml_dir, run) {
            Ok(p) => p,
            Err(_) => {
                tracing::warn!(run = %run, "mzML not found, skipping multi-target provenance");
                continue;
            }
        };

        let reader = match protein_copilot_spectrum_io::create_indexed_reader(&mzml_path) {
            Ok(r) => r,
            Err(e) => {
                tracing::warn!(file = %mzml_path.display(), error = %e, "failed to create reader");
                continue;
            }
        };

        // Extract DIA windows from this mzML
        let dia_windows = extract_dia_windows(reader.as_ref(), &mzml_path);

        // Build co-elution index for this run
        let index = coelution::CoElutionIndex::build(
            all_psms,
            all_groups,
            &dia_windows,
            config.provenance.silac.as_ref(),
            config.provenance.max_co_eluting_candidates,
        );

        // Process eligible PSMs for this run
        for &(idx, ref psm_run) in &eligible {
            if psm_run != run { continue; }

            let cpsm = &classified[idx];

            // Resolve scan number
            let scan_number = if let Some(s) = cpsm.psm.scan_number.filter(|&s| s > 0) {
                s
            } else if let (Some(rt), Some(mz)) = (cpsm.psm.retention_time, cpsm.psm.precursor_mz) {
                let tol = match (cpsm.psm.rt_start, cpsm.psm.rt_stop) {
                    (Some(start), Some(stop)) if stop > start => (stop - start) / 2.0,
                    _ => config.provenance.rt_tolerance_min,
                };
                match reader.find_by_rt(&mzml_path, rt, mz, tol) {
                    Ok(Some((scan, _))) => scan,
                    _ => continue,
                }
            } else {
                continue;
            };

            // Read spectrum
            let spectrum = match reader.read_spectrum(&mzml_path, scan_number) {
                Ok(s) => s,
                Err(_) => continue,
            };

            if (spectrum.mz_array.len() as u32) < config.provenance.min_peaks_for_analysis {
                continue;
            }

            // Find co-eluting targets
            let candidates = index.find_co_eluting(&cpsm.psm, run);
            if candidates.is_empty() {
                continue;
            }

            // Multi-target matching
            let mut prov = multi_provenance::trace_multi_target(
                &spectrum.mz_array,
                &spectrum.intensity_array,
                &cpsm.psm.peptide,
                &cpsm.psm.modifications,
                &candidates,
                &tolerance,
                config.provenance.max_fragment_charge,
            );
            prov.scan_number = scan_number;

            // Generate per-PSM HTML report
            if config.provenance.generate_per_psm_reports {
                let filename = format!(
                    "{}_z{}_{}.html",
                    cpsm.psm.peptide,
                    cpsm.psm.charge.unwrap_or(0),
                    run
                );
                let report_path = prov_dir.join(&filename);
                if let Err(e) = multi_report::render_multi_provenance_report(&prov, &report_path) {
                    tracing::warn!(file = %report_path.display(), error = %e, "failed to write per-PSM report");
                }
            }

            results.push(prov);
            traced_count += 1;
        }
    }

    Ok((traced_count, results))
}
```

- [ ] **Step 4: Build and verify compilation**

Run: `cargo build -p protein-copilot-entrapment-analysis`
Expected: BUILD SUCCESS.

- [ ] **Step 5: Commit**

```bash
git add crates/entrapment-analysis/src/lib.rs crates/entrapment-analysis/tests/v4_multi_target.rs
git commit -m "feat(entrapment): implement multi-target provenance batch pipeline"
```

---

### Task 7: CLI Integration

**Files:**
- Modify: `crates/entrapment-cli/src/main.rs`

- [ ] **Step 1: Add multi-target provenance call in run_analyze**

In `crates/entrapment-cli/src/main.rs`, after the existing v3 provenance tracing block (around line 196), add:

```rust
    // 5c. Multi-target provenance (v4) — if mzml_dir + silac config present
    if let Some(mzml_dir) = mzml_dir {
        use protein_copilot_entrapment_analysis::trace_multi_target_provenance;
        use protein_copilot_entrapment_analysis::tagger::Tagger;

        // Need group assignments for all PSMs to build co-elution index
        let tagger = Tagger::new(&config)?;
        let groups: Vec<protein_copilot_entrapment_analysis::PsmGroup> = psms
            .iter()
            .map(|psm| tagger.tag(&psm.protein_ids).unwrap_or(protein_copilot_entrapment_analysis::PsmGroup::Target))
            .collect();

        println!("Running multi-target provenance tracing...");
        match trace_multi_target_provenance(
            &classified,
            &psms,
            &groups,
            mzml_dir,
            &config,
            &out_dir,
        ) {
            Ok((count, _results)) => {
                println!("Multi-target provenance traced for {} PSMs", count);
                println!("Per-PSM reports written to: {}/provenance/", out_dir.display());
            }
            Err(e) => eprintln!("Warning: multi-target provenance failed: {}", e),
        }
    }
```

- [ ] **Step 2: Build and verify**

Run: `cargo build -p entrapment-cli`
Expected: BUILD SUCCESS.

- [ ] **Step 3: Commit**

```bash
git add crates/entrapment-cli/src/main.rs
git commit -m "feat(entrapment-cli): integrate multi-target provenance pipeline"
```

---

### Task 8: Summary Report Generator

**Files:**
- Modify: `crates/entrapment-analysis/src/multi_report.rs`

- [ ] **Step 1: Write failing test for summary report**

Add to the test module in `multi_report.rs`:

```rust
    #[test]
    fn test_summary_report() {
        let provs = vec![make_test_provenance()];
        let html = generate_provenance_summary_html(&provs);
        assert!(html.contains("Provenance Summary"));
        assert!(html.contains("STTTGHLIYK"));
        assert!(html.contains("<!DOCTYPE html>"));
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p protein-copilot-entrapment-analysis test_summary_report -- --nocapture`
Expected: FAIL — `generate_provenance_summary_html` not defined.

- [ ] **Step 3: Implement summary report generator**

Add to `multi_report.rs`:

```rust
/// Generate an HTML summary report for all multi-target provenance results.
pub fn generate_provenance_summary_html(results: &[MultiTargetProvenance]) -> String {
    let mut html = String::with_capacity(8_000);

    write!(html, r#"<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Provenance Summary</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 20px; color: #333; }}
h1 {{ color: #2c3e50; }}
table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
th {{ padding: 8px; text-align: left; background: #2c3e50; color: white; }}
td {{ padding: 6px 8px; }}
tr:nth-child(even) {{ background: #f5f5f5; }}
tr:hover {{ background: #ebf5fb; }}
a {{ color: #2980b9; text-decoration: none; }}
a:hover {{ text-decoration: underline; }}
.chimeric {{ color: #e74c3c; font-weight: bold; }}
</style>
</head><body>
<h1>🔬 Provenance Summary — {} PSMs Traced</h1>
<table>
<tr><th>Peptide</th><th>Scan</th><th>#Candidates</th><th>TrapOnly</th><th>Shared</th><th>TargetOnly</th><th>Unassigned</th><th>Report</th></tr>
"#, results.len()).ok();

    for prov in results {
        let filename = format!("{}_z0_{}.html", prov.trap_peptide, prov.scan_number);
        let total = prov.annotated_peaks.len();
        let shared_pct = if total > 0 { prov.shared_count as f64 / total as f64 * 100.0 } else { 0.0 };

        write!(html, r#"<tr>
<td style="font-family:monospace"><b>{pep}</b></td>
<td>{scan}</td>
<td>{ncand}</td>
<td>{trap}</td>
<td>{shared}{chimeric}</td>
<td>{target}</td>
<td>{unassigned}</td>
<td><a href="provenance/{file}">View →</a></td>
</tr>"#,
            pep = prov.trap_peptide,
            scan = prov.scan_number,
            ncand = prov.candidates.len(),
            trap = prov.trap_only_count,
            shared = prov.shared_count,
            chimeric = if shared_pct > 30.0 { format!(r#" <span class="chimeric">({:.0}%)</span>"#, shared_pct) } else { String::new() },
            target = prov.target_only_count,
            unassigned = prov.unassigned_count,
            file = filename,
        ).ok();
    }

    html.push_str("</table></body></html>");
    html
}

/// Write summary report to file.
pub fn render_provenance_summary(
    results: &[MultiTargetProvenance],
    output_path: &Path,
) -> Result<(), std::io::Error> {
    let html = generate_provenance_summary_html(results);
    std::fs::write(output_path, html)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p protein-copilot-entrapment-analysis test_summary_report -- --nocapture`
Expected: PASS.

- [ ] **Step 5: Wire summary report into batch pipeline**

In `lib.rs`, add at the end of `trace_multi_target_provenance`, before `Ok((traced_count, results))`:

```rust
    // Generate summary report
    if !results.is_empty() && config.provenance.generate_per_psm_reports {
        let summary_path = output_dir.join("provenance_summary.html");
        if let Err(e) = multi_report::render_provenance_summary(&results, &summary_path) {
            tracing::warn!(error = %e, "failed to write provenance summary report");
        } else {
            tracing::info!(
                path = %summary_path.display(),
                count = results.len(),
                "wrote provenance summary report"
            );
        }
    }
```

- [ ] **Step 6: Commit**

```bash
git add crates/entrapment-analysis/src/multi_report.rs crates/entrapment-analysis/src/lib.rs
git commit -m "feat(entrapment): add provenance summary report with per-PSM links"
```

---

### Task 9: Update Example Config + Documentation

**Files:**
- Modify: `examples/hela-mix-2da-entrapment.yaml`
- Modify: `docs/entrapment-analysis.md`

- [ ] **Step 1: Add silac config to example YAML**

In `examples/hela-mix-2da-entrapment.yaml`, add under the `provenance:` section:

```yaml
  silac:
    heavy_k_delta: 8.014199
    heavy_r_delta: 10.008269
    enable_heavy_search: true
  generate_per_psm_reports: true
  max_co_eluting_candidates: 20
```

- [ ] **Step 2: Update user documentation**

In `docs/entrapment-analysis.md`, add a new section after the "版本历程" section:

```markdown
### v4 使用指南

#### 多目标碎片溯源

v4 对每个 L2/L3 trap PSM 自动查找所有共洗脱的 target 肽段（轻标 + 重标 SILAC），将每个观测碎片离子归属到具体的 target 来源。

**前置条件：**
- `--mzml-dir` 参数指向 mzML 文件目录
- config 中配置 `provenance.silac` 块（可选，启用重标搜索）

**输出：**
- `provenance_summary.html` — 所有溯源 PSMs 的汇总表
- `provenance/` 目录 — 每个 PSM 一份 HTML 报告

**配置参数：**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `silac.heavy_k_delta` | 8.014199 | 重标 K delta mass |
| `silac.heavy_r_delta` | 10.008269 | 重标 R delta mass |
| `silac.enable_heavy_search` | true | 是否搜索重标候选 |
| `generate_per_psm_reports` | true | 是否生成 per-PSM HTML |
| `max_co_eluting_candidates` | 20 | 每个 trap 的最大候选数 |
```

- [ ] **Step 3: Commit**

```bash
git add examples/hela-mix-2da-entrapment.yaml docs/entrapment-analysis.md
git commit -m "docs: add v4 multi-target provenance config and user guide"
```

---

### Task 10: Full Build + Test Verification

- [ ] **Step 1: Run full workspace build**

Run: `cargo build --workspace`
Expected: BUILD SUCCESS with no errors.

- [ ] **Step 2: Run full test suite**

Run: `cargo test --workspace`
Expected: All tests PASS.

- [ ] **Step 3: Run clippy**

Run: `cargo clippy --workspace -- -D warnings`
Expected: No warnings.

- [ ] **Step 4: Fix any issues found**

Address any compiler warnings, test failures, or clippy lints.

- [ ] **Step 5: Final commit if fixes needed**

```bash
git add -A
git commit -m "fix: address clippy warnings and test issues"
```

---

### Task 11: Real Data End-to-End Test

- [ ] **Step 1: Run CLI with real hela-mix-2da data**

```bash
cd /home/verden/pfind/2026-spring/code/proteinCopilot
cargo run -p entrapment-cli -- analyze \
  --results output/hela-mix-2da_report.parquet \
  --config examples/hela-mix-2da-entrapment.yaml \
  --target-fasta .proteincopilot/databases/human_swissprot.fasta \
  --mzml-dir /home/verden/pfind/2025-fall/code/2da/ \
  --out output/entrapment-v4
```

- [ ] **Step 2: Verify output files**

Check that:
- `output/entrapment-v4/provenance_summary.html` exists and opens in browser
- `output/entrapment-v4/provenance/` contains per-PSM HTML files
- Per-PSM reports contain candidate tables, mirror spectra, and attribution tables
- Light + heavy candidates are shown where applicable

- [ ] **Step 3: Spot-check a specific report**

Open one per-PSM HTML and verify:
- Candidates have correct DIA window info
- Mirror spectrum is interactive (Plotly hover works)
- Attribution table shows correct origin classifications

- [ ] **Step 4: Commit if all good**

```bash
git add -A
git commit -m "test: verify v4 multi-target provenance with real hela-mix-2da data"
```
