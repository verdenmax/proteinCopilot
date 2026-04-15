//! Integration tests for spectrum annotation across 4 scenarios:
//! DDA/DIA × non-SILAC/SILAC + zero-offset edge case.

use protein_copilot_core::search_params::Modification;
use protein_copilot_core::spectrum::Spectrum;
use protein_copilot_search_engine::annotate::{annotate_heavy_spectrum, annotate_spectrum};
use test_helpers::*;

/// Helper: annotate a spectrum and verify basic invariants.
fn assert_annotation_ok(spectrum: &Spectrum, peptide: &str, charge: i32) {
    let tol = default_frag_tolerance();
    let mods: Vec<Modification> = vec![];
    let result = annotate_spectrum(spectrum, peptide, charge, &tol, &mods, vec![], false, false);
    assert!(
        result.is_ok(),
        "annotation should succeed: {:?}",
        result.err()
    );
    let ann = result.unwrap();
    assert_eq!(ann.peptide_sequence, peptide);
    assert_eq!(ann.charge, charge);
    assert!(ann.total_ions > 0, "should have theoretical ions");
}

// ─── Scenario ①: DDA + non-SILAC ──────────────────────────────────────────

#[test]
fn scenario_1_dda_no_silac_annotates_light_only() {
    let peptide = "PEPTIDEK";
    let precursor_mz = 458.24; // approximate [M+2H]²⁺
    let peaks = synthetic_peaks_for_peptide(peptide, 1000.0);
    let spectrum = make_ms2(
        1,
        30.0,
        precursor_mz,
        2,
        Some(dda_window(precursor_mz)),
        peaks,
    );

    let tol = default_frag_tolerance();
    let mods: Vec<Modification> = vec![];
    let ann = annotate_spectrum(&spectrum, peptide, 2, &tol, &mods, vec![], false, false).unwrap();

    assert_eq!(ann.peptide_sequence, peptide);
    assert!(ann.matched_ions > 0, "should match some fragment ions");
    // No heavy annotation expected
    assert!(
        ann.heavy_annotation.is_none(),
        "DDA+non-SILAC should have no heavy annotation"
    );
}

// ─── Scenario ②: DDA + SILAC ─────────────────────────────────────────────

#[test]
fn scenario_2_dda_silac_heavy_annotation_succeeds() {
    let peptide = "PEPTIDEK"; // has K → non-zero SILAC shift
    let precursor_mz = 458.24;
    let label = silac_label();

    // Create light spectrum
    let light_peaks = synthetic_peaks_for_peptide(peptide, 1000.0);
    let light_spectrum = make_ms2(
        10,
        30.0,
        precursor_mz,
        2,
        Some(dda_window(precursor_mz)),
        light_peaks.clone(),
    );

    // Light annotation should work
    assert_annotation_ok(&light_spectrum, peptide, 2);

    // Heavy annotation should work on a spectrum with shifted peaks
    let heavy_peaks = heavy_shifted_peaks(peptide, &light_peaks, &label);
    let heavy_prec_mz =
        protein_copilot_core::label::compute_heavy_precursor_mz(precursor_mz, 2, peptide, &label);
    let heavy_spectrum = make_ms2(
        11,
        30.1,
        heavy_prec_mz,
        2,
        Some(dda_window(heavy_prec_mz)),
        heavy_peaks,
    );

    let tol = default_frag_tolerance();
    let mods: Vec<Modification> = vec![];
    let heavy_ann = annotate_heavy_spectrum(
        &heavy_spectrum,
        peptide,
        2,
        &tol,
        &mods,
        &label,
        false,
        false,
    )
    .unwrap();

    assert!(
        heavy_ann.total_ions > 0,
        "should have heavy theoretical ions"
    );
    // heavy_ann.precursor_mz is computed from peptide mass, not from input spectrum
    assert!(
        heavy_ann.precursor_mz > 0.0,
        "should have valid heavy precursor m/z"
    );
    // The delta should be positive (K adds mass)
    assert!(
        heavy_ann.precursor_mz > precursor_mz,
        "heavy m/z ({}) should exceed light m/z ({})",
        heavy_ann.precursor_mz,
        precursor_mz
    );
}

// ─── Scenario ③: DIA + non-SILAC ──────────────────────────────────────────

#[test]
fn scenario_3_dia_no_silac_annotates_light_only() {
    let peptide = "DGFLLDGFPR"; // has R
    let precursor_mz = 547.28;
    let peaks = synthetic_peaks_for_peptide(peptide, 1000.0);
    let spectrum = make_ms2(
        1,
        45.0,
        precursor_mz,
        2,
        Some(dia_window(precursor_mz)),
        peaks,
    );

    let tol = default_frag_tolerance();
    let mods: Vec<Modification> = vec![];
    let ann = annotate_spectrum(&spectrum, peptide, 2, &tol, &mods, vec![], false, false).unwrap();

    assert_eq!(ann.peptide_sequence, peptide);
    assert!(ann.matched_ions > 0, "should match some fragment ions");
    assert!(
        ann.heavy_annotation.is_none(),
        "DIA+non-SILAC should have no heavy annotation"
    );
}

// ─── Scenario ④: DIA + SILAC ──────────────────────────────────────────────

#[test]
fn scenario_4_dia_silac_heavy_annotation_succeeds() {
    let peptide = "DGFLLDGFPR"; // has R → non-zero SILAC shift
    let precursor_mz = 547.28;
    let label = silac_label();

    let light_peaks = synthetic_peaks_for_peptide(peptide, 1000.0);
    let light_spectrum = make_ms2(
        20,
        45.0,
        precursor_mz,
        2,
        Some(dia_window(precursor_mz)),
        light_peaks.clone(),
    );

    // Light annotation
    assert_annotation_ok(&light_spectrum, peptide, 2);

    // Heavy annotation
    let heavy_peaks = heavy_shifted_peaks(peptide, &light_peaks, &label);
    let heavy_prec_mz =
        protein_copilot_core::label::compute_heavy_precursor_mz(precursor_mz, 2, peptide, &label);
    let heavy_spectrum = make_ms2(
        21,
        45.1,
        heavy_prec_mz,
        2,
        Some(dia_window(heavy_prec_mz)),
        heavy_peaks,
    );

    let tol = default_frag_tolerance();
    let mods: Vec<Modification> = vec![];
    let heavy_ann = annotate_heavy_spectrum(
        &heavy_spectrum,
        peptide,
        2,
        &tol,
        &mods,
        &label,
        false,
        false,
    )
    .unwrap();

    assert!(heavy_ann.total_ions > 0);
}

// ─── Zero-offset edge case ────────────────────────────────────────────────

#[test]
fn zero_offset_peptide_no_kr_skips_heavy() {
    let peptide = "PEPTIDE"; // no K or R → zero SILAC delta
    let label = silac_label();

    let delta = protein_copilot_core::label::total_heavy_delta(peptide, &label);
    assert!(
        delta.abs() < 1e-6,
        "PEPTIDE without K/R should have zero delta, got {delta}"
    );

    // Heavy precursor m/z should equal light
    let light_mz = 400.19;
    let heavy_mz =
        protein_copilot_core::label::compute_heavy_precursor_mz(light_mz, 2, peptide, &label);
    assert!(
        (heavy_mz - light_mz).abs() < 0.001,
        "zero-offset: heavy m/z ({heavy_mz}) should equal light ({light_mz})"
    );

    // annotate_heavy_spectrum still works on a zero-offset peptide
    // (it's the caller's job to skip, but the function itself should not crash)
    let peaks = synthetic_peaks_for_peptide(peptide, 1000.0);
    let spectrum = make_ms2(1, 30.0, light_mz, 2, Some(dda_window(light_mz)), peaks);
    let tol = default_frag_tolerance();
    let mods: Vec<Modification> = vec![];
    let result = annotate_heavy_spectrum(&spectrum, peptide, 2, &tol, &mods, &label, false, false);
    assert!(
        result.is_ok(),
        "annotate_heavy_spectrum should not crash on zero-offset peptide"
    );
}
