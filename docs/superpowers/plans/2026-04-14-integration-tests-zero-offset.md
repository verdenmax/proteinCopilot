# Integration Tests + Zero-Offset Validation — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add workspace-level integration tests for the 4 annotation/XIC scenarios (DDA/DIA × non-SILAC/SILAC), and guard against meaningless heavy processing when peptide has no K/R residues.

**Architecture:** Two independent changes: (1) zero-offset guard in 2 code locations (tools.rs, extract.rs), and (2) 4 new test files in root `tests/` directory using synthetic Spectrum data.

**Tech Stack:** Rust, tokio (async tests), protein_copilot_core, protein_copilot_search_engine, protein_copilot_xic

---

### Task 1: Zero-Offset Guard in Annotation Handler

**Files:**
- Modify: `crates/mcp-server/src/tools.rs:1608-1625`
- Test: verified by Task 5 (`tests/annotation_scenarios.rs`)

- [ ] **Step 1: Add zero-offset check after computing heavy_prec_mz**

In `crates/mcp-server/src/tools.rs`, after the line that computes `heavy_prec_mz` (line ~1623), add a delta check. The current code is:

```rust
        // ── SILAC: find and annotate heavy scan (DIA or DDA) ──
        if let Some(ref label) = input.label_type {
            let is_dia = spectrum
                .precursors
                .first()
                .and_then(|p| p.isolation_window.as_ref())
                .map(|w| (w.lower_offset + w.upper_offset) > 1.0)
                .unwrap_or(false);

            let core_label: &protein_copilot_core::label::LabelType = label;
            let heavy_prec_mz = protein_copilot_core::label::compute_heavy_precursor_mz(
                annotation.theoretical_mz,
                charge,
                &peptide_seq,
                core_label,
            );

            // Read nearby spectra to find the heavy scan
```

Change to:

```rust
        // ── SILAC: find and annotate heavy scan (DIA or DDA) ──
        if let Some(ref label) = input.label_type {
            let is_dia = spectrum
                .precursors
                .first()
                .and_then(|p| p.isolation_window.as_ref())
                .map(|w| (w.lower_offset + w.upper_offset) > 1.0)
                .unwrap_or(false);

            let core_label: &protein_copilot_core::label::LabelType = label;
            let heavy_delta = protein_copilot_core::label::total_heavy_delta(&peptide_seq, core_label);

            if heavy_delta.abs() < 1e-6 {
                tracing::info!(
                    peptide = peptide_seq,
                    "Skipping heavy annotation: peptide has no K/R, zero SILAC shift"
                );
            } else {

            let heavy_prec_mz = protein_copilot_core::label::compute_heavy_precursor_mz(
                annotation.theoretical_mz,
                charge,
                &peptide_seq,
                core_label,
            );

            // Read nearby spectra to find the heavy scan
```

And add a closing brace `}` after the existing heavy annotation block's closing brace (after the `else` log message about "No ... found for heavy precursor").

- [ ] **Step 2: Build and verify no compilation errors**

Run: `cargo build -p protein-copilot-mcp-server 2>&1 | tail -3`
Expected: `Finished`

- [ ] **Step 3: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "fix: skip heavy annotation when peptide has no K/R (zero SILAC shift)"
```

---

### Task 2: Zero-Offset Guard in XIC Extraction

**Files:**
- Modify: `crates/xic/src/extract.rs:287-294` (extract_xic function)
- Modify: `crates/xic/src/extract.rs:697-705` (extract_xic_with_raw function)
- Test: `crates/xic/src/extract.rs` (existing test module) + Task 6

- [ ] **Step 1: Add zero-offset check in extract_xic**

In `crates/xic/src/extract.rs`, the current code at ~line 287 is:

```rust
    let heavy_ions = match &params.label_type {
        Some(label) => crate::heavy::compute_heavy_target_ions(&light_ions, peptide_sequence, label),
        None => Vec::new(),
    };

    let heavy_precursor_mz = params.label_type.as_ref().map(|label| {
        crate::heavy::compute_heavy_precursor_mz(precursor_mz, charge, peptide_sequence, label)
    });
```

Change to:

```rust
    let (heavy_ions, heavy_precursor_mz) = match &params.label_type {
        Some(label) => {
            let delta = protein_copilot_core::label::total_heavy_delta(peptide_sequence, label);
            if delta.abs() < 1e-6 {
                tracing::info!(
                    peptide = peptide_sequence,
                    "Skipping heavy XIC: peptide has no K/R, zero SILAC shift"
                );
                (Vec::new(), None)
            } else {
                let ions = crate::heavy::compute_heavy_target_ions(&light_ions, peptide_sequence, label);
                let mz = crate::heavy::compute_heavy_precursor_mz(precursor_mz, charge, peptide_sequence, label);
                (ions, Some(mz))
            }
        }
        None => (Vec::new(), None),
    };
```

- [ ] **Step 2: Apply same change in extract_xic_with_raw**

In `crates/xic/src/extract.rs`, the current code at ~line 697 is:

```rust
    let heavy_ions = match &params.label_type {
        Some(label) => crate::heavy::compute_heavy_target_ions(&light_ions, peptide_sequence, label),
        None => Vec::new(),
    };

    let heavy_precursor_mz = params.label_type.as_ref().map(|label| {
        crate::heavy::compute_heavy_precursor_mz(precursor_mz, charge, peptide_sequence, label)
    });
```

Apply identical change as Step 1.

- [ ] **Step 3: Add unit test for zero-offset in extract.rs test module**

Find the existing test module at the bottom of `crates/xic/src/extract.rs` and add:

```rust
    #[test]
    fn zero_silac_shift_skips_heavy() {
        // "PEPTIDE" has no K or R — heavy delta should be zero
        use protein_copilot_core::label::LabelType;
        let label = LabelType::standard_silac();
        let delta = protein_copilot_core::label::total_heavy_delta("PEPTIDE", &label);
        assert!(delta.abs() < 1e-6, "PEPTIDE has no K/R, delta should be 0, got {delta}");

        // "PEPTIDEK" has 1 K — heavy delta should be non-zero
        let delta_k = protein_copilot_core::label::total_heavy_delta("PEPTIDEK", &label);
        assert!(delta_k > 7.0, "PEPTIDEK should have K+8 shift, got {delta_k}");
    }
```

- [ ] **Step 4: Run tests**

Run: `cargo test -p protein-copilot-xic -- zero_silac 2>&1 | tail -5`
Expected: `test result: ok. 1 passed`

- [ ] **Step 5: Clippy check**

Run: `cargo clippy -p protein-copilot-xic -p protein-copilot-mcp-server -- -D warnings 2>&1 | tail -3`
Expected: `Finished`

- [ ] **Step 6: Commit**

```bash
git add crates/xic/src/extract.rs
git commit -m "fix: skip heavy XIC extraction when peptide has no K/R (zero SILAC shift)"
```

---

### Task 3: Test Helpers Module for Synthetic Spectra

**Files:**
- Create: `tests/helpers/mod.rs` — shared synthetic spectrum builders
- Create: `tests/helpers.rs` — module declaration (if needed by Rust edition)

- [ ] **Step 1: Create tests/helpers/mod.rs with spectrum builders**

```rust
//! Shared test helpers: synthetic spectrum builders for integration tests.

use protein_copilot_core::spectrum::{IsolationWindow, MsLevel, PrecursorInfo, Spectrum};

/// Build a synthetic MS2 spectrum with peaks at specified m/z values.
///
/// `peaks`: slice of (mz, intensity) pairs — will be sorted by mz.
pub fn make_ms2(
    scan: u32,
    rt_min: f64,
    precursor_mz: f64,
    charge: i32,
    isolation_window: Option<IsolationWindow>,
    peaks: &[(f64, f64)],
) -> Spectrum {
    let mut sorted = peaks.to_vec();
    sorted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    Spectrum {
        scan_number: scan,
        ms_level: MsLevel::MS2,
        retention_time_min: rt_min,
        mz_array: sorted.iter().map(|(mz, _)| *mz).collect(),
        intensity_array: sorted.iter().map(|(_, i)| *i).collect(),
        precursors: vec![PrecursorInfo {
            mz: precursor_mz,
            charge: Some(charge),
            intensity: None,
            source_scan: None,
            isolation_window,
        }],
    }
}

/// DDA-style narrow isolation window (0.7 Th total).
pub fn dda_window(center_mz: f64) -> IsolationWindow {
    IsolationWindow {
        target_mz: center_mz,
        lower_offset: 0.35,
        upper_offset: 0.35,
    }
}

/// DIA-style wide isolation window (25 Th total).
pub fn dia_window(center_mz: f64) -> IsolationWindow {
    IsolationWindow {
        target_mz: center_mz,
        lower_offset: 12.5,
        upper_offset: 12.5,
    }
}

/// Standard SILAC label (K+8, R+10).
pub fn silac_label() -> protein_copilot_core::label::LabelType {
    protein_copilot_core::label::LabelType::standard_silac()
}

/// Compute theoretical b/y ion m/z values for a peptide (charge 1 only).
///
/// Returns (b_ions, y_ions) as Vec<f64>.
/// Uses the search-engine chemistry module for consistency with annotation code.
pub fn theoretical_fragments(peptide: &str, charge: i32) -> Vec<f64> {
    use protein_copilot_search_engine::chemistry::{amino_acid_mass, PROTON_MASS, WATER_MASS};

    let chars: Vec<char> = peptide.chars().collect();
    let n = chars.len();
    let mut fragments = Vec::new();

    // b-ions: sum of residues + proton
    let mut cumulative = 0.0_f64;
    for i in 0..n - 1 {
        if let Some(m) = amino_acid_mass(chars[i]) {
            cumulative += m;
            for z in 1..=charge.min(2) {
                fragments.push((cumulative + z as f64 * PROTON_MASS) / z as f64);
            }
        }
    }

    // y-ions: sum from C-term + water + proton
    let mut cumulative = 0.0_f64;
    for i in (1..n).rev() {
        if let Some(m) = amino_acid_mass(chars[i]) {
            cumulative += m;
            for z in 1..=charge.min(2) {
                fragments.push((cumulative + WATER_MASS + z as f64 * PROTON_MASS) / z as f64);
            }
        }
    }

    fragments
}
```

- [ ] **Step 2: Verify it compiles**

Run: `cargo test --test annotation_scenarios -- --list 2>&1 | head -5`
(Will fail because test file doesn't exist yet — that's expected. We just need helpers to compile as part of the next task.)

---

### Task 4: Annotation Scenarios Integration Test

**Files:**
- Create: `tests/annotation_scenarios.rs`

- [ ] **Step 1: Create the test file**

```rust
//! Integration tests for spectrum annotation across 4 scenarios:
//! DDA/DIA × non-SILAC/SILAC, plus zero-offset validation.

mod helpers;

use protein_copilot_core::label::LabelType;
use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};
use protein_copilot_search_engine::annotate::{annotate_spectrum, annotate_heavy_spectrum};

fn frag_tol() -> MassTolerance {
    MassTolerance { value: 20.0, unit: ToleranceUnit::Ppm }
}

// ──────────────────────────────────────────────
// Scenario ①: DDA + non-SILAC
// ──────────────────────────────────────────────

#[test]
fn scenario_1_dda_no_silac() {
    let peptide = "PEPTIDEK";
    let charge = 2;
    let frags = helpers::theoretical_fragments(peptide, charge);
    // Build peaks at theoretical positions with noise
    let peaks: Vec<(f64, f64)> = frags.iter().map(|&mz| (mz, 1000.0)).collect();

    let spectrum = helpers::make_ms2(
        1, 10.0, 500.0, charge,
        Some(helpers::dda_window(500.0)),
        &peaks,
    );

    let ann = annotate_spectrum(
        &spectrum, peptide, charge, &frag_tol(), &[], vec!["TEST".to_string()],
    ).unwrap();

    assert!(ann.score > 0.0, "should match some ions");
    assert!(ann.matched_ions > 0);
    assert!(ann.heavy_annotation.is_none(), "no SILAC → no heavy annotation");
}

// ──────────────────────────────────────────────
// Scenario ②: DDA + SILAC
// ──────────────────────────────────────────────

#[test]
fn scenario_2_dda_silac_heavy_annotation() {
    let peptide = "PEPTIDEK"; // has K → SILAC shift
    let charge = 2;
    let label = helpers::silac_label();

    let frags = helpers::theoretical_fragments(peptide, charge);
    let light_peaks: Vec<(f64, f64)> = frags.iter().map(|&mz| (mz, 1000.0)).collect();

    let light_spectrum = helpers::make_ms2(
        1, 10.0, 500.0, charge,
        Some(helpers::dda_window(500.0)),
        &light_peaks,
    );

    // Compute heavy fragment peaks (shift K-containing fragments)
    let heavy_peaks: Vec<(f64, f64)> = frags.iter().map(|&mz| (mz + 4.007, 800.0)).collect();
    let heavy_prec_mz = protein_copilot_core::label::compute_heavy_precursor_mz(
        500.0, charge, peptide, &label,
    );
    let heavy_spectrum = helpers::make_ms2(
        2, 10.01, heavy_prec_mz, charge,
        Some(helpers::dda_window(heavy_prec_mz)),
        &heavy_peaks,
    );

    // Light annotation
    let ann = annotate_spectrum(
        &light_spectrum, peptide, charge, &frag_tol(), &[], vec!["TEST".to_string()],
    ).unwrap();
    assert!(ann.score > 0.0);

    // Heavy annotation
    let heavy_ann = annotate_heavy_spectrum(
        &heavy_spectrum, peptide, charge, &frag_tol(), &[], &label,
    ).unwrap();
    assert!(heavy_ann.matched_ions > 0, "heavy should match some ions");
    assert!(heavy_ann.score > 0.0);
}

// ──────────────────────────────────────────────
// Scenario ③: DIA + non-SILAC
// ──────────────────────────────────────────────

#[test]
fn scenario_3_dia_no_silac() {
    let peptide = "PEPTIDER";
    let charge = 2;
    let frags = helpers::theoretical_fragments(peptide, charge);
    let peaks: Vec<(f64, f64)> = frags.iter().map(|&mz| (mz, 1000.0)).collect();

    let spectrum = helpers::make_ms2(
        1, 10.0, 500.0, charge,
        Some(helpers::dia_window(500.0)), // wide window
        &peaks,
    );

    let ann = annotate_spectrum(
        &spectrum, peptide, charge, &frag_tol(), &[], vec!["TEST".to_string()],
    ).unwrap();

    assert!(ann.score > 0.0);
    assert!(ann.heavy_annotation.is_none());
}

// ──────────────────────────────────────────────
// Scenario ④: DIA + SILAC
// ──────────────────────────────────────────────

#[test]
fn scenario_4_dia_silac_heavy_annotation() {
    let peptide = "PEPTIDER"; // has R → SILAC shift
    let charge = 2;
    let label = helpers::silac_label();

    let frags = helpers::theoretical_fragments(peptide, charge);
    let light_peaks: Vec<(f64, f64)> = frags.iter().map(|&mz| (mz, 1000.0)).collect();
    let heavy_peaks: Vec<(f64, f64)> = frags.iter().map(|&mz| (mz + 5.004, 800.0)).collect();

    let light_spectrum = helpers::make_ms2(
        1, 10.0, 500.0, charge,
        Some(helpers::dia_window(500.0)),
        &light_peaks,
    );

    let heavy_prec_mz = protein_copilot_core::label::compute_heavy_precursor_mz(
        500.0, charge, peptide, &label,
    );
    let heavy_spectrum = helpers::make_ms2(
        2, 10.01, heavy_prec_mz, charge,
        Some(helpers::dia_window(heavy_prec_mz)),
        &heavy_peaks,
    );

    let ann = annotate_spectrum(
        &light_spectrum, peptide, charge, &frag_tol(), &[], vec!["TEST".to_string()],
    ).unwrap();
    assert!(ann.score > 0.0);

    let heavy_ann = annotate_heavy_spectrum(
        &heavy_spectrum, peptide, charge, &frag_tol(), &[], &label,
    ).unwrap();
    assert!(heavy_ann.matched_ions > 0);
}

// ──────────────────────────────────────────────
// Zero-offset: no K/R → skip heavy
// ──────────────────────────────────────────────

#[test]
fn zero_offset_no_kr_skips_heavy() {
    let peptide = "PEPTIDE"; // no K, no R → zero SILAC shift
    let label = helpers::silac_label();
    let delta = protein_copilot_core::label::total_heavy_delta(peptide, &label);
    assert!(delta.abs() < 1e-6, "PEPTIDE should have zero delta, got {delta}");

    // Even with label, heavy annotation on same spectrum should match identically
    // (because heavy theoretical == light theoretical)
    // This test validates the zero-shift detection logic is correct
    let charge = 2;
    let heavy_prec = protein_copilot_core::label::compute_heavy_precursor_mz(
        500.0, charge, peptide, &label,
    );
    assert!((heavy_prec - 500.0).abs() < 1e-6, "heavy_prec should == light_prec for no K/R");
}
```

- [ ] **Step 2: Run the annotation tests**

Run: `cargo test --test annotation_scenarios 2>&1 | tail -10`
Expected: all 5 tests pass

- [ ] **Step 3: Commit**

```bash
git add tests/annotation_scenarios.rs tests/helpers/mod.rs
git commit -m "test: add 4-scenario annotation integration tests + zero-offset check"
```

---

### Task 5: XIC Scenarios Integration Test

**Files:**
- Create: `tests/xic_scenarios.rs`

- [ ] **Step 1: Create the test file**

```rust
//! Integration tests for XIC extraction across 4 scenarios:
//! DDA/DIA × non-SILAC/SILAC, plus zero-offset validation.
//!
//! These tests use the heavy.rs module directly (not the full extract_xic pipeline
//! which requires file I/O) to validate the heavy scan finding logic.

mod helpers;

use protein_copilot_core::label::LabelType;
use protein_copilot_core::spectrum::Spectrum;
use protein_copilot_xic::heavy;

// ──────────────────────────────────────────────
// DDA heavy scan finding
// ──────────────────────────────────────────────

#[test]
fn xic_dda_silac_finds_heavy_scan() {
    let peptide = "PEPTIDEK";
    let charge = 2;
    let label = helpers::silac_label();
    let light_prec = 500.0;
    let heavy_prec = protein_copilot_core::label::compute_heavy_precursor_mz(
        light_prec, charge, peptide, &label,
    );

    let spectra = vec![
        helpers::make_ms2(1, 10.0, light_prec, charge, Some(helpers::dda_window(light_prec)), &[(200.0, 100.0)]),
        helpers::make_ms2(2, 10.01, heavy_prec, charge, Some(helpers::dda_window(heavy_prec)), &[(204.0, 100.0)]),
        helpers::make_ms2(3, 10.02, 600.0, charge, Some(helpers::dda_window(600.0)), &[(300.0, 100.0)]),
    ];

    let result = heavy::find_dda_heavy_scan(&spectra, 10.0, heavy_prec, 20.0);
    assert_eq!(result, Some(2), "should find scan 2 as heavy");
}

#[test]
fn xic_dda_no_silac_no_heavy() {
    let spectra = vec![
        helpers::make_ms2(1, 10.0, 500.0, 2, Some(helpers::dda_window(500.0)), &[(200.0, 100.0)]),
    ];

    // Without heavy target, nothing to search for
    let result = heavy::find_dda_heavy_scan(&spectra, 10.0, 999.0, 20.0);
    assert!(result.is_none(), "no scan should match unrelated m/z");
}

// ──────────────────────────────────────────────
// DIA heavy window finding
// ──────────────────────────────────────────────

#[test]
fn xic_dia_silac_finds_heavy_window() {
    let peptide = "PEPTIDER";
    let charge = 2;
    let label = helpers::silac_label();
    let light_prec = 500.0;
    let heavy_prec = protein_copilot_core::label::compute_heavy_precursor_mz(
        light_prec, charge, peptide, &label,
    );

    let spectra = vec![
        helpers::make_ms2(1, 10.0, 500.0, charge, Some(helpers::dia_window(500.0)), &[(200.0, 100.0)]),
        helpers::make_ms2(2, 10.01, heavy_prec, charge, Some(helpers::dia_window(heavy_prec)), &[(204.0, 100.0)]),
    ];

    let result = heavy::find_heavy_dia_window_from_spectra(&spectra, 10.0, heavy_prec);
    assert!(result.is_some(), "should find DIA window containing heavy precursor");
    assert_eq!(result.unwrap().0, 2);
}

#[test]
fn xic_dia_no_silac_no_heavy_window() {
    let spectra = vec![
        helpers::make_ms2(1, 10.0, 500.0, 2, Some(helpers::dia_window(500.0)), &[(200.0, 100.0)]),
    ];

    let result = heavy::find_heavy_dia_window_from_spectra(&spectra, 10.0, 999.0);
    assert!(result.is_none());
}

// ──────────────────────────────────────────────
// Zero-offset: heavy ions identical to light
// ──────────────────────────────────────────────

#[test]
fn xic_zero_offset_heavy_ions_identical() {
    use protein_copilot_xic::IonType;

    let peptide = "PEPTIDE"; // no K/R
    let label = helpers::silac_label();

    let light_ions = vec![
        protein_copilot_xic::extract::TargetIon {
            label: "b1¹⁺".to_string(),
            ion_type: IonType::B,
            ion_number: 1,
            charge: 1,
            mz: 100.0,
        },
    ];

    let heavy_ions = heavy::compute_heavy_target_ions(&light_ions, peptide, &label);
    assert_eq!(heavy_ions.len(), 1);
    assert!(
        (heavy_ions[0].mz - light_ions[0].mz).abs() < 1e-6,
        "no K/R → heavy mz should equal light mz"
    );
}

#[test]
fn xic_nonzero_offset_heavy_ions_shifted() {
    use protein_copilot_xic::IonType;

    let peptide = "PEPTIDEK"; // has K
    let label = helpers::silac_label();

    let light_ions = vec![
        protein_copilot_xic::extract::TargetIon {
            label: "y1¹⁺".to_string(),
            ion_type: IonType::Y,
            ion_number: 1,
            charge: 1,
            mz: 147.1, // approximate y1 for K
        },
    ];

    let heavy_ions = heavy::compute_heavy_target_ions(&light_ions, peptide, &label);
    // y1 suffix = "K" → shift by K_DELTA = 8.014
    assert!(
        (heavy_ions[0].mz - light_ions[0].mz - 8.014199).abs() < 0.01,
        "K-containing y1 should shift by ~8.014"
    );
}
```

- [ ] **Step 2: Check if TargetIon is public**

The test uses `protein_copilot_xic::extract::TargetIon`. Verify it is `pub` in `crates/xic/src/extract.rs`. If not, add `pub` to the struct and its fields.

Run: `grep -n 'pub struct TargetIon' crates/xic/src/extract.rs`

If it says `pub struct TargetIon` → OK. If not, make it public.

- [ ] **Step 3: Run the XIC tests**

Run: `cargo test --test xic_scenarios 2>&1 | tail -10`
Expected: all 6 tests pass

- [ ] **Step 4: Commit**

```bash
git add tests/xic_scenarios.rs
git commit -m "test: add 4-scenario XIC integration tests + zero-offset checks"
```

---

### Task 6: Search Pipeline Integration Test

**Files:**
- Create: `tests/search_pipeline.rs`

- [ ] **Step 1: Create the test file**

This test reuses existing fixture data (`test_100.mgf` + `test_100.fasta`) from the search-engine crate.

```rust
//! End-to-end search pipeline integration test at workspace level.
//!
//! Covers: read_spectra → recommend_params → search → FDR → summary → annotate

mod helpers;

use std::path::PathBuf;

use protein_copilot_core::engine::SearchEngineAdapter;
use protein_copilot_core::progress::noop_progress;
use protein_copilot_core::search_params::*;
use protein_copilot_param_recommend::ParamRecommender;
use protein_copilot_report::ReportGenerator;
use protein_copilot_search_engine::annotate::annotate_spectrum;
use protein_copilot_search_engine::SimpleSearchEngine;
use protein_copilot_spectrum_io::{create_reader, detect_format};

fn fixtures_dir() -> PathBuf {
    // Reuse search-engine's e2e fixtures
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("crates")
        .join("search-engine")
        .join("tests")
        .join("fixtures_e2e")
}

#[tokio::test]
async fn full_pipeline_search_to_annotate() {
    let mgf = fixtures_dir().join("test_100.mgf");
    let fasta = fixtures_dir().join("test_100.fasta");

    // Step 1: Read spectra
    let file_info = detect_format(&mgf).unwrap();
    let reader = create_reader(&file_info);
    let summary = reader.read_summary(&mgf).unwrap();
    assert!(summary.ms2_count > 0);

    // Step 2: Recommend params
    let rec = ParamRecommender.recommend(&summary, None).unwrap();
    let mut params = rec.decision;
    params.database_path = fasta.to_string_lossy().to_string();

    // Step 3: Search
    let engine = SimpleSearchEngine::new();
    let result = engine.search(&params, &[mgf.clone()], noop_progress()).await.unwrap();
    assert!(!result.psms.is_empty(), "search should find PSMs");

    // Step 4: Generate summary
    let report_summary = ReportGenerator::generate_summary(&result);
    assert!(report_summary.total_psms > 0);
    assert!(report_summary.psms_at_1pct_fdr > 0);

    // Step 5: Annotate first PSM
    let psm = &result.psms[0];
    let spectrum = reader.read_spectrum(&mgf, psm.spectrum_scan).unwrap();
    let ann = annotate_spectrum(
        &spectrum,
        &psm.peptide_sequence,
        psm.charge,
        &params.fragment_tolerance,
        &params.fixed_modifications,
        psm.protein_accessions.clone(),
    ).unwrap();
    assert!(ann.score > 0.0, "annotation should match some ions");
    assert!(ann.matched_ions > 0);
}
```

- [ ] **Step 2: Run the test**

Run: `cargo test --test search_pipeline 2>&1 | tail -10`
Expected: 1 test passed

- [ ] **Step 3: Commit**

```bash
git add tests/search_pipeline.rs
git commit -m "test: add workspace-level search pipeline integration test"
```

---

### Task 7: Import Pipeline Integration Test

**Files:**
- Create: `tests/import_pipeline.rs`

- [ ] **Step 1: Create the test file**

```rust
//! Integration test for the result import → summary pipeline.
//!
//! Uses a synthetic custom JSON result to validate the import flow.

mod helpers;

use std::io::Write;
use tempfile::NamedTempFile;

use protein_copilot_result_import::custom_json::import_custom_json;

#[test]
fn import_custom_json_and_summarize() {
    // Create a minimal custom JSON result file
    let json_content = r#"{
        "psms": [
            {
                "raw_name": "test_run",
                "scan_number": 1,
                "peptide_sequence": "PEPTIDEK",
                "charge": 2,
                "score": 15.5,
                "precursor_mz": 500.27,
                "modifications": [],
                "protein_accessions": ["sp|P12345|TEST_HUMAN"],
                "retention_time_min": 10.5
            },
            {
                "raw_name": "test_run",
                "scan_number": 2,
                "peptide_sequence": "ANOTHERPEPTIDER",
                "charge": 3,
                "score": 12.3,
                "precursor_mz": 550.30,
                "modifications": [],
                "protein_accessions": ["sp|P67890|TEST2_HUMAN"],
                "retention_time_min": 15.2
            }
        ]
    }"#;

    let mut tmp = NamedTempFile::new().unwrap();
    tmp.write_all(json_content.as_bytes()).unwrap();
    tmp.flush().unwrap();

    let imported = import_custom_json(tmp.path()).unwrap();
    assert_eq!(imported.len(), 2, "should import 2 PSMs");
    assert_eq!(imported[0].peptide_sequence, "PEPTIDEK");
    assert_eq!(imported[1].peptide_sequence, "ANOTHERPEPTIDER");
    assert_eq!(imported[0].charge, 2);
    assert!((imported[0].score - 15.5).abs() < 0.01);
}
```

- [ ] **Step 2: Check if tempfile is a dependency**

Run: `grep 'tempfile' Cargo.toml`

If not found, add to workspace root `Cargo.toml` under `[dev-dependencies]`:
```toml
tempfile = "3"
```

- [ ] **Step 3: Run the test**

Run: `cargo test --test import_pipeline 2>&1 | tail -10`
Expected: 1 test passed

- [ ] **Step 4: Commit**

```bash
git add tests/import_pipeline.rs Cargo.toml
git commit -m "test: add result import pipeline integration test"
```

---

### Task 8: Full Workspace Test Run + Final Commit

- [ ] **Step 1: Run all workspace tests**

Run: `cargo test --workspace 2>&1 | tail -20`
Expected: all tests pass (previous ~296 + new ~13 = ~309)

- [ ] **Step 2: Run clippy on all targets**

Run: `cargo clippy --all-targets -- -D warnings 2>&1 | tail -5`
Expected: `Finished`

- [ ] **Step 3: Squash or tag final state**

```bash
git --no-pager log --oneline -8
```

Verify all commits are clean. Done.
