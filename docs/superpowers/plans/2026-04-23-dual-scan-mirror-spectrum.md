# Dual-Scan Mirror Spectrum Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Light and Heavy mirror spectra read from separate DIA MS2 scans (different isolation windows) instead of sharing one scan, with the trap's heavy ions generated for the heavy mirror.

**Architecture:** In DIA, light and heavy SILAC precursors have different m/z and fall into different isolation windows → different MS2 scans at the same RT. We split `MultiTargetProvenance` into `light: MirrorData` + `heavy: Option<MirrorData>`, each holding its own scan number and annotated peaks. The pipeline reads two scans via `find_by_rt` (one for light precursor m/z, one for heavy), runs matching separately, and the report renders each mirror from its own observed spectrum.

**Tech Stack:** Rust, protein-copilot-entrapment-analysis crate, protein-copilot-spectrum-io (find_by_rt)

**Key insight:** `find_by_rt(rt, precursor_mz, tolerance)` already filters by isolation window. Calling it with heavy precursor m/z naturally finds the scan whose DIA window covers the heavy m/z.

**Trap heavy ions:** For the heavy mirror, we also generate the trap peptide's theoretical heavy ions (shifted by SILAC config), so we can see if any observed peak in the heavy scan matches what the trap *would* look like in heavy form.

---

### Task 1: Add MirrorData struct and restructure MultiTargetProvenance

**Files:**
- Modify: `crates/entrapment-analysis/src/types.rs`

This task replaces the flat `annotated_peaks` + `scan_number` + counters on `MultiTargetProvenance` with a `MirrorData` sub-struct, duplicated for `light` and `heavy`.

- [ ] **Step 1: Add MirrorData struct to types.rs**

Add after `MultiAnnotatedPeak` (around line 309):

```rust
/// Data for one mirror spectrum (light or heavy scan).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MirrorData {
    /// Scan number of the observed MS2 spectrum.
    pub scan_number: u32,
    /// Per-peak multi-target annotation.
    pub annotated_peaks: Vec<MultiAnnotatedPeak>,
    /// Count of peaks matching only trap ions.
    pub trap_only_count: u32,
    /// Count of peaks matching at least one target (not trap).
    pub target_only_count: u32,
    /// Count of peaks matching both trap and at least one target.
    pub shared_count: u32,
    /// Count of peaks matching nothing.
    pub unassigned_count: u32,
}
```

- [ ] **Step 2: Restructure MultiTargetProvenance**

Replace the current `MultiTargetProvenance` struct with:

```rust
/// Complete multi-target provenance result for one trap PSM.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiTargetProvenance {
    /// The trap PSM being analyzed.
    pub trap_peptide: String,
    /// Precursor m/z (light form) of the trap PSM.
    pub trap_precursor_mz: f64,
    /// Precursor m/z (heavy SILAC form) of the trap PSM, if applicable.
    pub trap_precursor_mz_heavy: Option<f64>,
    /// Charge state of the trap PSM.
    pub trap_charge: i32,
    /// Spectrum/raw file name.
    pub spectrum_file: String,
    /// All co-eluting target candidates (light + heavy).
    /// Candidate indices in TargetIonMatch reference this vec.
    pub candidates: Vec<CoElutingCandidate>,
    /// Light mirror data (observed from the light-precursor DIA scan).
    pub light: MirrorData,
    /// Heavy mirror data (observed from the heavy-precursor DIA scan).
    /// None if no heavy precursor or no matching scan found.
    pub heavy: Option<MirrorData>,
}
```

Remove the old flat fields: `scan_number`, `annotated_peaks`, `trap_only_count`, `target_only_count`, `shared_count`, `unassigned_count`.

- [ ] **Step 3: Verify types.rs compiles in isolation**

Run: `cargo check -p protein-copilot-entrapment-analysis 2>&1 | head -5`
Expected: compilation errors in dependent files (multi_provenance.rs, multi_report.rs, lib.rs) — that's OK, we fix them in subsequent tasks.

- [ ] **Step 4: Commit**

```bash
git add crates/entrapment-analysis/src/types.rs
git commit -m "refactor(types): add MirrorData struct, restructure MultiTargetProvenance for dual-scan"
```

---

### Task 2: Update trace_multi_target to return MirrorData

**Files:**
- Modify: `crates/entrapment-analysis/src/multi_provenance.rs`

The core matching function now returns `MirrorData` instead of `MultiTargetProvenance`. The caller assembles the full provenance. Also make `shift_ions_heavy` and `generate_candidate_ions` pub(crate) so the caller can generate trap heavy ions.

- [ ] **Step 1: Change trace_multi_target return type to MirrorData**

```rust
use crate::types::{
    CoElutingCandidate, LabelForm, MirrorData, MultiAnnotatedPeak, TargetIonMatch,
};

/// Trace fragment ion provenance for one mirror (one observed spectrum).
///
/// Matches observed peaks against trap theoretical ions and candidate
/// theoretical ions. Returns `MirrorData` with per-peak annotations.
///
/// `scan_number` is set by the caller after this function returns.
pub fn trace_multi_target(
    observed_mz: &[f64],
    observed_intensity: &[f64],
    trap_sequence: &str,
    trap_modifications: &[(usize, f64)],
    candidates: &[CoElutingCandidate],
    fragment_tolerance: &MassTolerance,
    max_fragment_charge: i32,
) -> MirrorData {
    // 1. Generate theoretical ions for the trap peptide.
    let trap_ions =
        generate_theoretical_ions(trap_sequence, trap_modifications, max_fragment_charge);

    // 2. Generate theoretical ions for each candidate.
    let candidate_ions: Vec<Vec<TheoreticalIon>> = candidates
        .iter()
        .map(|c| generate_candidate_ions(c, max_fragment_charge))
        .collect();

    // 3. Classify each observed peak.
    let mut annotated_peaks = Vec::with_capacity(observed_mz.len());
    let mut trap_only_count = 0u32;
    let mut target_only_count = 0u32;
    let mut shared_count = 0u32;
    let mut unassigned_count = 0u32;

    for (i, &mz) in observed_mz.iter().enumerate() {
        let intensity = observed_intensity.get(i).copied().unwrap_or(0.0);

        let trap_match = find_best_match(mz, &trap_ions, fragment_tolerance);
        let trap_ion = trap_match.map(|(label, _ppm)| label);

        let mut target_matches = Vec::new();
        for (ci, ions) in candidate_ions.iter().enumerate() {
            if let Some((label, delta_ppm)) =
                find_best_match_with_ppm(mz, ions, fragment_tolerance)
            {
                target_matches.push(TargetIonMatch {
                    candidate_index: ci,
                    ion_label: label,
                    delta_ppm,
                });
            }
        }

        match (trap_ion.is_some(), !target_matches.is_empty()) {
            (true, true) => shared_count += 1,
            (true, false) => trap_only_count += 1,
            (false, true) => target_only_count += 1,
            (false, false) => unassigned_count += 1,
        }

        annotated_peaks.push(MultiAnnotatedPeak {
            mz_observed: mz,
            intensity,
            trap_ion,
            target_matches,
        });
    }

    MirrorData {
        scan_number: 0, // caller sets
        annotated_peaks,
        trap_only_count,
        target_only_count,
        shared_count,
        unassigned_count,
    }
}
```

- [ ] **Step 2: Add trace_mirror_with_trap_ions for heavy mirror**

Add a new public function that takes pre-generated trap ions instead of generating them internally. This allows the caller to pass heavy-shifted trap ions for the heavy mirror.

```rust
/// Trace fragment ion provenance using pre-generated trap ions.
///
/// Used for the heavy mirror: the caller generates the trap's heavy
/// theoretical ions (shifted by SILAC deltas) and passes them directly.
pub fn trace_mirror_with_trap_ions(
    observed_mz: &[f64],
    observed_intensity: &[f64],
    trap_ions: &[TheoreticalIon],
    candidates: &[CoElutingCandidate],
    fragment_tolerance: &MassTolerance,
    max_fragment_charge: i32,
) -> MirrorData {
    let candidate_ions: Vec<Vec<TheoreticalIon>> = candidates
        .iter()
        .map(|c| generate_candidate_ions(c, max_fragment_charge))
        .collect();

    let mut annotated_peaks = Vec::with_capacity(observed_mz.len());
    let mut trap_only_count = 0u32;
    let mut target_only_count = 0u32;
    let mut shared_count = 0u32;
    let mut unassigned_count = 0u32;

    for (i, &mz) in observed_mz.iter().enumerate() {
        let intensity = observed_intensity.get(i).copied().unwrap_or(0.0);

        let trap_match = find_best_match(mz, trap_ions, fragment_tolerance);
        let trap_ion = trap_match.map(|(label, _ppm)| label);

        let mut target_matches = Vec::new();
        for (ci, ions) in candidate_ions.iter().enumerate() {
            if let Some((label, delta_ppm)) =
                find_best_match_with_ppm(mz, ions, fragment_tolerance)
            {
                target_matches.push(TargetIonMatch {
                    candidate_index: ci,
                    ion_label: label,
                    delta_ppm,
                });
            }
        }

        match (trap_ion.is_some(), !target_matches.is_empty()) {
            (true, true) => shared_count += 1,
            (true, false) => trap_only_count += 1,
            (false, true) => target_only_count += 1,
            (false, false) => unassigned_count += 1,
        }

        annotated_peaks.push(MultiAnnotatedPeak {
            mz_observed: mz,
            intensity,
            trap_ion,
            target_matches,
        });
    }

    MirrorData {
        scan_number: 0,
        annotated_peaks,
        trap_only_count,
        target_only_count,
        shared_count,
        unassigned_count,
    }
}
```

- [ ] **Step 3: Make shift_ions_heavy and generate_theoretical_ions accessible**

`shift_ions_heavy` is currently private. The pipeline in `lib.rs` needs it to generate trap heavy ions. Make it `pub(crate)`:

```rust
pub(crate) fn shift_ions_heavy(...)  // was: fn shift_ions_heavy(...)
```

Also import and re-export `TheoreticalIon` from `provenance.rs` if not already pub(crate). Check that `generate_theoretical_ions` in `provenance.rs` is `pub(crate)`.

- [ ] **Step 4: Update existing tests**

The tests in `multi_provenance.rs` reference `MultiTargetProvenance` — update them to use `MirrorData` instead. For example:

```rust
#[test]
fn test_match_single_target_light() {
    // ... (same setup) ...
    let result = trace_multi_target(
        &observed_mz, &observed_int, trap_seq, &[],
        &candidates, &tolerance_20ppm(), 1,
    );

    // result is now MirrorData, not MultiTargetProvenance
    assert!(result.shared_count > 0 || result.trap_only_count > 0);
    assert_eq!(result.annotated_peaks.len(), observed_mz.len());
}

#[test]
fn test_scan_number_is_zero() {
    let result = trace_multi_target(
        &[100.0], &[1000.0], "AG", &[], &[],
        &tolerance_20ppm(), 1,
    );
    assert_eq!(result.scan_number, 0);
}
```

Update all 11 tests to expect `MirrorData` return type. Remove references to `result.trap_peptide`, `result.candidates`.

- [ ] **Step 5: Verify tests compile and pass**

Run: `cargo test -p protein-copilot-entrapment-analysis -- multi_provenance`
Expected: all tests pass.

- [ ] **Step 6: Commit**

```bash
git add crates/entrapment-analysis/src/multi_provenance.rs
git commit -m "refactor(multi_provenance): return MirrorData, add trace_mirror_with_trap_ions"
```

---

### Task 3: Update pipeline to read dual scans

**Files:**
- Modify: `crates/entrapment-analysis/src/lib.rs`

This is the core pipeline change. For each eligible trap PSM:
1. Read the light scan (existing logic, using trap's light precursor m/z)
2. Compute trap heavy precursor m/z from SILAC config
3. Read the heavy scan via `find_by_rt(rt, heavy_mz, tol)`
4. Run `trace_multi_target` for light candidates against light spectrum
5. Generate trap heavy ions, run `trace_mirror_with_trap_ions` for heavy candidates against heavy spectrum
6. Assemble `MultiTargetProvenance { light, heavy, ... }`

- [ ] **Step 1: Update the per-PSM processing loop**

In `trace_multi_target_provenance()`, replace the section from "Trace multi-target provenance" to "Generate per-PSM HTML report" (approximately lines 603–631) with:

```rust
            // Partition candidates into light and heavy groups.
            // Candidate indices remain relative to the full `candidates` vec.
            let light_indices: Vec<usize> = candidates
                .iter()
                .enumerate()
                .filter(|(_, c)| matches!(c.label_form, LabelForm::Light))
                .map(|(i, _)| i)
                .collect();
            let heavy_indices: Vec<usize> = candidates
                .iter()
                .enumerate()
                .filter(|(_, c)| !matches!(c.label_form, LabelForm::Light))
                .map(|(i, _)| i)
                .collect();

            // --- Light mirror ---
            let mut light_mirror = multi_provenance::trace_multi_target(
                &spectrum.mz_array,
                &spectrum.intensity_array,
                &cpsm.psm.peptide,
                &cpsm.psm.modifications,
                &candidates,
                &tolerance,
                config.provenance.max_fragment_charge,
            );
            light_mirror.scan_number = scan_number;

            // --- Heavy mirror ---
            let heavy_mirror = if let Some(silac) = &config.provenance.silac {
                // Compute trap heavy precursor m/z.
                let seq_chars: Vec<char> = cpsm.psm.peptide.chars().collect();
                let total_delta: f64 = seq_chars
                    .iter()
                    .map(|&c| match c {
                        'K' => silac.heavy_k_delta,
                        'R' => silac.heavy_r_delta,
                        _ => 0.0,
                    })
                    .sum();

                if total_delta > 0.0 && !heavy_indices.is_empty() {
                    let charge = cpsm.psm.charge.unwrap_or(2) as f64;
                    let heavy_mz = cpsm.psm.precursor_mz.unwrap_or(0.0) + total_delta / charge;

                    // Find the heavy scan at the same RT.
                    let heavy_tol = match (cpsm.psm.rt_start, cpsm.psm.rt_stop) {
                        (Some(start), Some(stop)) if stop > start => (stop - start) / 2.0,
                        _ => config.provenance.rt_tolerance_min,
                    };
                    let heavy_rt = cpsm.psm.retention_time.unwrap_or(0.0);

                    match reader.find_by_rt(&mzml_path, heavy_rt, heavy_mz, heavy_tol) {
                        Ok(Some((heavy_scan, _delta))) => {
                            // Read the heavy spectrum.
                            match reader.read_spectrum(&mzml_path, heavy_scan) {
                                Ok(heavy_spec) => {
                                    // Generate trap heavy ions.
                                    let residue_deltas: Vec<(usize, f64)> = seq_chars
                                        .iter()
                                        .enumerate()
                                        .filter_map(|(i, &c)| match c {
                                            'K' => Some((i, silac.heavy_k_delta)),
                                            'R' => Some((i, silac.heavy_r_delta)),
                                            _ => None,
                                        })
                                        .collect();

                                    let trap_light_ions =
                                        crate::provenance::generate_theoretical_ions(
                                            &cpsm.psm.peptide,
                                            &cpsm.psm.modifications,
                                            config.provenance.max_fragment_charge,
                                        );
                                    let trap_heavy_ions =
                                        multi_provenance::shift_ions_heavy(
                                            &cpsm.psm.peptide,
                                            &trap_light_ions,
                                            &residue_deltas,
                                            config.provenance.max_fragment_charge,
                                        );

                                    let mut heavy_data =
                                        multi_provenance::trace_mirror_with_trap_ions(
                                            &heavy_spec.mz_array,
                                            &heavy_spec.intensity_array,
                                            &trap_heavy_ions,
                                            &candidates,
                                            &tolerance,
                                            config.provenance.max_fragment_charge,
                                        );
                                    heavy_data.scan_number = heavy_scan;

                                    Some(heavy_data)
                                }
                                Err(err) => {
                                    tracing::debug!(
                                        scan = heavy_scan,
                                        error = %err,
                                        "could not read heavy scan, skipping heavy mirror"
                                    );
                                    None
                                }
                            }
                        }
                        Ok(None) => {
                            tracing::debug!(
                                peptide = %cpsm.psm.peptide,
                                heavy_mz = heavy_mz,
                                "no heavy scan found, skipping heavy mirror"
                            );
                            None
                        }
                        Err(err) => {
                            tracing::debug!(error = %err, "heavy scan lookup failed");
                            None
                        }
                    }
                } else {
                    None
                }
            } else {
                None
            };

            // Compute trap heavy precursor m/z for metadata.
            let trap_precursor_mz_heavy = config.provenance.silac.as_ref().and_then(|silac| {
                let seq_chars: Vec<char> = cpsm.psm.peptide.chars().collect();
                let total_delta: f64 = seq_chars
                    .iter()
                    .map(|&c| match c {
                        'K' => silac.heavy_k_delta,
                        'R' => silac.heavy_r_delta,
                        _ => 0.0,
                    })
                    .sum();
                if total_delta > 0.0 {
                    let charge = cpsm.psm.charge.unwrap_or(2) as f64;
                    Some(cpsm.psm.precursor_mz.unwrap_or(0.0) + total_delta / charge)
                } else {
                    None
                }
            });

            let prov = MultiTargetProvenance {
                trap_peptide: cpsm.psm.peptide.clone(),
                trap_precursor_mz: cpsm.psm.precursor_mz.unwrap_or(0.0),
                trap_precursor_mz_heavy,
                trap_charge: cpsm.psm.charge.unwrap_or(0),
                spectrum_file: run_name.clone(),
                candidates,
                light: light_mirror,
                heavy: heavy_mirror,
            };
```

- [ ] **Step 2: Remove the old flat-field assignments**

Delete the old lines that set `prov.scan_number`, `prov.trap_precursor_mz`, etc. since they're now set in the struct constructor above. Also remove the old `trace_multi_target` call and the SILAC heavy precursor computation block that was added in the previous commit (they're now integrated into the new code).

- [ ] **Step 3: Verify it compiles**

Run: `cargo check -p protein-copilot-entrapment-analysis 2>&1 | head -20`
Expected: may have errors in multi_report.rs (fixed in Task 4).

- [ ] **Step 4: Commit**

```bash
git add crates/entrapment-analysis/src/lib.rs
git commit -m "feat(pipeline): read dual DIA scans for light/heavy mirrors"
```

---

### Task 4: Update report renderer to use MirrorData

**Files:**
- Modify: `crates/entrapment-analysis/src/multi_report.rs`

The report renderer now reads from `prov.light` and `prov.heavy` instead of `prov.annotated_peaks`.

- [ ] **Step 1: Update generate_multi_provenance_html**

Replace references to `prov.annotated_peaks` with `prov.light.annotated_peaks`. Update scan number references:

In `write_header`: display both scan numbers:
```rust
// Scan numbers
let light_scan = prov.light.scan_number;
let heavy_scan_str = match &prov.heavy {
    Some(h) => format!("{}", h.scan_number),
    None => "N/A".to_string(),
};
// In the info-grid:
// <div><span class="label">Scan (Light):</span> <span class="value">{light_scan}</span></div>
// <div><span class="label">Scan (Heavy):</span> <span class="value">{heavy_scan_str}</span></div>
```

For mirror rendering:
```rust
// Light mirror — always rendered
write_mirror_spectrum_from_data(html, prov, &prov.light, MirrorKind::Light);

// Heavy mirror — only if heavy data exists
if let Some(heavy) = &prov.heavy {
    write_mirror_spectrum_from_data(html, prov, heavy, MirrorKind::Heavy);
}
```

- [ ] **Step 2: Refactor write_mirror_spectrum to take MirrorData**

Rename to `write_mirror_spectrum_from_data` with signature:

```rust
fn write_mirror_spectrum_from_data(
    html: &mut String,
    prov: &MultiTargetProvenance,  // for candidates list
    mirror: &MirrorData,           // annotated peaks + scan number
    kind: MirrorKind,
)
```

Change all references from `prov.annotated_peaks` to `mirror.annotated_peaks` and `prov.scan_number` to `mirror.scan_number` inside this function.

The candidate filtering by Light/Heavy still uses `prov.candidates` and `kind`.

- [ ] **Step 3: Update attribution table to use light mirror**

Change `write_attribution_table` to reference `prov.light.annotated_peaks` (the primary mirror for trap ion attribution):

```rust
fn write_attribution_table(html: &mut String, prov: &MultiTargetProvenance) {
    // ...
    let max_intensity = prov.light.annotated_peaks
        .iter()
        .map(|p| p.intensity)
        .fold(0.0_f64, f64::max)
        .max(1.0);

    for peak in &prov.light.annotated_peaks {
        // ... same logic, reading from light mirror
    }
}
```

- [ ] **Step 4: Update footer counts**

Footer should show combined or per-mirror counts:

```rust
let light = &prov.light;
// Footer shows light mirror stats (primary)
write!(html, "Light: TrapOnly={} Shared={} TargetOnly={} Unassigned={}", 
    light.trap_only_count, light.shared_count, 
    light.target_only_count, light.unassigned_count);
if let Some(heavy) = &prov.heavy {
    write!(html, " | Heavy: TrapOnly={} Shared={} TargetOnly={} Unassigned={}",
        heavy.trap_only_count, heavy.shared_count,
        heavy.target_only_count, heavy.unassigned_count);
}
```

- [ ] **Step 5: Update summary report**

In `generate_provenance_summary_html`, change `prov.scan_number` to `prov.light.scan_number` and counts to `prov.light.*`:

```rust
let scan = prov.light.scan_number;
// Use light mirror counts for summary
let total = prov.light.trap_only_count + prov.light.shared_count 
    + prov.light.target_only_count + prov.light.unassigned_count;
let shared_frac = if total > 0 {
    f64::from(prov.light.shared_count) / f64::from(total)
} else { 0.0 };
```

- [ ] **Step 6: Lower bold line width**

Change line widths from `2.5` to `1.5` in both `write_normalized_trace` and `write_target_trace_with_bold`:

```rust
let line_width = if bold { 1.5 } else { 0.0 };
// ...
let line_widths: Vec<String> = bold_flags.iter()
    .map(|&b| if b { "1.5".to_string() } else { "0".to_string() })
    .collect();
```

- [ ] **Step 7: Commit**

```bash
git add crates/entrapment-analysis/src/multi_report.rs
git commit -m "feat(report): render light/heavy mirrors from separate scans, lower bold width"
```

---

### Task 5: Update tests and integration test

**Files:**
- Modify: `crates/entrapment-analysis/src/multi_report.rs` (test section)
- Modify: `crates/entrapment-analysis/tests/v4_multi_target.rs`

- [ ] **Step 1: Update make_test_provenance in multi_report.rs tests**

```rust
fn make_test_provenance() -> MultiTargetProvenance {
    MultiTargetProvenance {
        trap_peptide: "STTTGHLIYK".to_string(),
        trap_precursor_mz: 547.789,
        trap_precursor_mz_heavy: Some(556.803),
        trap_charge: 2,
        spectrum_file: "550_600_2Da_Rep1".to_string(),
        candidates: vec![CoElutingCandidate {
            peptide: "STTSGHLVYK".to_string(),
            protein_ids: vec!["sp|P12345|EF1A_HUMAN".to_string()],
            precursor_mz: 548.12,
            charge: 2,
            rt_start: 34.5,
            rt_stop: 35.8,
            label_form: LabelForm::Light,
            modifications: vec![],
        }],
        light: MirrorData {
            scan_number: 12345,
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
        },
        heavy: None,
    }
}
```

- [ ] **Step 2: Update test assertions**

Tests that referenced `prov.scan_number` → `prov.light.scan_number`.
Tests that checked `prov.trap_only_count` → `prov.light.trap_only_count`.

For the chimeric test:
```rust
let prov = MultiTargetProvenance {
    trap_peptide: "CHIMERIC".to_string(),
    trap_precursor_mz: 500.0,
    trap_precursor_mz_heavy: None,
    trap_charge: 2,
    spectrum_file: "test_run".to_string(),
    candidates: vec![],
    light: MirrorData {
        scan_number: 999,
        annotated_peaks: vec![],
        trap_only_count: 1,
        target_only_count: 0,
        shared_count: 5,
        unassigned_count: 4,
    },
    heavy: None,
};
```

Similarly update the empty candidates test.

- [ ] **Step 3: Add test for dual-scan provenance**

Add a test that constructs `MultiTargetProvenance` with both `light` and `heavy` `MirrorData`, and verifies the HTML contains both mirror sections:

```rust
#[test]
fn test_dual_scan_both_mirrors() {
    let prov = MultiTargetProvenance {
        trap_peptide: "PEPTIDEK".to_string(),
        trap_precursor_mz: 450.0,
        trap_precursor_mz_heavy: Some(454.007),
        trap_charge: 2,
        spectrum_file: "test_run".to_string(),
        candidates: vec![
            CoElutingCandidate {
                peptide: "PEPTIDER".to_string(),
                protein_ids: vec!["P1".to_string()],
                precursor_mz: 451.0,
                charge: 2,
                rt_start: 30.0,
                rt_stop: 32.0,
                label_form: LabelForm::Light,
                modifications: vec![],
            },
            CoElutingCandidate {
                peptide: "PEPTIDER".to_string(),
                protein_ids: vec!["P1".to_string()],
                precursor_mz: 451.0,
                charge: 2,
                rt_start: 30.0,
                rt_stop: 32.0,
                label_form: LabelForm::Heavy {
                    precursor_mz_heavy: 456.0,
                    residue_deltas: vec![(7, 10.008269)],
                },
                modifications: vec![],
            },
        ],
        light: MirrorData {
            scan_number: 100,
            annotated_peaks: vec![MultiAnnotatedPeak {
                mz_observed: 300.0,
                intensity: 1000.0,
                trap_ion: Some("b3+1".to_string()),
                target_matches: vec![TargetIonMatch {
                    candidate_index: 0,
                    ion_label: "b3+1".to_string(),
                    delta_ppm: 1.0,
                }],
            }],
            trap_only_count: 0,
            target_only_count: 0,
            shared_count: 1,
            unassigned_count: 0,
        },
        heavy: Some(MirrorData {
            scan_number: 105,
            annotated_peaks: vec![MultiAnnotatedPeak {
                mz_observed: 305.0,
                intensity: 800.0,
                trap_ion: Some("b3+1(H)".to_string()),
                target_matches: vec![TargetIonMatch {
                    candidate_index: 1,
                    ion_label: "b3+1(H)".to_string(),
                    delta_ppm: 2.0,
                }],
            }],
            trap_only_count: 0,
            target_only_count: 0,
            shared_count: 1,
            unassigned_count: 0,
        }),
    };

    let html = generate_multi_provenance_html(&prov);
    // Both mirrors present
    assert!(html.contains("Light Targets (Scan 100)"));
    assert!(html.contains("Heavy Targets (Scan 105)"));
    // Both scan numbers in header
    assert!(html.contains("100"));
    assert!(html.contains("105"));
}
```

- [ ] **Step 4: Update v4_multi_target integration test**

View `crates/entrapment-analysis/tests/v4_multi_target.rs` and update any references to the old struct layout.

- [ ] **Step 5: Run all tests**

Run: `cargo test -p protein-copilot-entrapment-analysis 2>&1 | tail -20`
Expected: all tests pass.

- [ ] **Step 6: Commit**

```bash
git add crates/entrapment-analysis/
git commit -m "test: update all tests for dual-scan MirrorData structure"
```

---

### Task 6: Build, clippy, and full test verification

**Files:** None (verification only)

- [ ] **Step 1: Full workspace build**

Run: `cargo build --workspace 2>&1 | tail -5`
Expected: builds successfully.

- [ ] **Step 2: Clippy**

Run: `cargo clippy --workspace 2>&1 | grep "^warning:" | grep -v "generated"`
Expected: only pre-existing `too_many_arguments` warnings.

- [ ] **Step 3: Full test suite**

Run: `cargo test --workspace 2>&1 | tail -20`
Expected: all tests pass.

- [ ] **Step 4: Commit any remaining fixes**

If clippy or tests revealed issues, fix and commit.

---

### Task 7: Update visual preview and reduce bold

**Files:**
- Create: `/tmp/v4-report-preview-v2.html` (temporary, not committed)

- [ ] **Step 1: Update visual companion**

Regenerate the preview HTML with:
- Two mirrors showing different scan numbers (Light: Scan 48721, Heavy: Scan 48725)
- Lower bold width (1.5px instead of 2.5px)
- Header showing both scan numbers
- Heavy mirror showing fewer shared peaks (demonstrating the diagnostic value)

- [ ] **Step 2: Verify in browser**

Serve and verify at `http://localhost:8765/v4-report-preview-v2.html`.
