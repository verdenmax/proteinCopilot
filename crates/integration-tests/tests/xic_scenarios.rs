//! Integration tests for XIC-related functions across DDA/DIA × non-SILAC/SILAC.
//!
//! Tests heavy scan finding logic and zero-offset behavior.

use protein_copilot_core::label::{compute_heavy_precursor_mz, total_heavy_delta};
use protein_copilot_xic::heavy::{
    compute_heavy_target_ions, find_dda_heavy_scan, find_heavy_dia_window_from_spectra,
    window_contains_mz,
};
use protein_copilot_xic::extract::build_target_ions;
use test_helpers::*;

// ─── DDA heavy scan finding ───────────────────────────────────────────────

#[test]
fn dda_find_heavy_scan_by_precursor_mz() {
    let peptide = "PEPTIDEK";
    let label = silac_label();
    let light_mz = 458.24;
    let heavy_mz = compute_heavy_precursor_mz(light_mz, 2, peptide, &label);

    // Create a set of nearby DDA MS2 spectra
    let spectra = vec![
        make_ms2(10, 30.0, light_mz, 2, Some(dda_window(light_mz)), vec![(200.0, 100.0)]),
        make_ms2(11, 30.05, heavy_mz, 2, Some(dda_window(heavy_mz)), vec![(200.0, 100.0)]),
        make_ms2(12, 30.1, 600.0, 2, Some(dda_window(600.0)), vec![(200.0, 100.0)]),
    ];

    let result = find_dda_heavy_scan(&spectra, 30.0, heavy_mz, 20.0);
    assert_eq!(result, Some(11), "should find scan 11 matching heavy m/z");
}

#[test]
fn dda_no_match_returns_none() {
    let spectra = vec![
        make_ms2(10, 30.0, 500.0, 2, Some(dda_window(500.0)), vec![(200.0, 100.0)]),
        make_ms2(11, 30.1, 501.0, 2, Some(dda_window(501.0)), vec![(200.0, 100.0)]),
    ];

    // Look for m/z far from any spectrum
    let result = find_dda_heavy_scan(&spectra, 30.0, 700.0, 20.0);
    assert_eq!(result, None, "no spectrum near 700 m/z");
}

// ─── DIA heavy scan finding ──────────────────────────────────────────────

#[test]
fn dia_find_heavy_scan_by_window_containment() {
    let peptide = "DGFLLDGFPR";
    let label = silac_label();
    let light_mz = 547.28;
    let heavy_mz = compute_heavy_precursor_mz(light_mz, 2, peptide, &label);

    // Create DIA windows — use different center m/z so only one contains heavy m/z
    // Heavy m/z is ~552.3 (light 547.28 + ~5.00 for R/2)
    // Window at 400 does NOT contain heavy (range 387.5-412.5)
    // Window at 552 DOES contain heavy (range 539.5-564.5)
    let spectra = vec![
        make_ms2(20, 45.0, 400.0, 2, Some(dia_window(400.0)), vec![(200.0, 100.0)]),
        make_ms2(21, 45.0, heavy_mz, 2, Some(dia_window(heavy_mz)), vec![(200.0, 100.0)]),
        make_ms2(22, 45.0, 900.0, 2, Some(dia_window(900.0)), vec![(200.0, 100.0)]),
    ];

    let result = find_heavy_dia_window_from_spectra(&spectra, 45.0, heavy_mz);
    assert!(result.is_some(), "should find DIA window containing heavy m/z");
    let (scan, window) = result.unwrap();
    assert_eq!(scan, 21);
    assert!(window_contains_mz(&window, heavy_mz));
}

#[test]
fn dia_no_window_contains_heavy_mz() {
    // Create windows that don't cover the target
    let spectra = vec![
        make_ms2(20, 45.0, 400.0, 2, Some(dia_window(400.0)), vec![(200.0, 100.0)]),
        make_ms2(21, 45.0, 500.0, 2, Some(dia_window(500.0)), vec![(200.0, 100.0)]),
    ];

    // Target heavy m/z is far outside both windows
    let result = find_heavy_dia_window_from_spectra(&spectra, 45.0, 800.0);
    assert!(result.is_none(), "no window should contain 800 m/z");
}

// ─── Zero-offset: heavy target ions ──────────────────────────────────────

#[test]
fn zero_offset_no_heavy_target_ions_shift() {
    let peptide = "PEPTIDE"; // no K or R
    let label = silac_label();

    let delta = total_heavy_delta(peptide, &label);
    assert!(delta.abs() < 1e-6, "PEPTIDE has no K/R, zero delta");

    let light_ions = build_target_ions(peptide, &[], 2);
    let heavy_ions = compute_heavy_target_ions(&light_ions, peptide, &label);

    // Heavy ions should exist but have identical m/z to light ions
    assert_eq!(light_ions.len(), heavy_ions.len());
    for (light, heavy) in light_ions.iter().zip(heavy_ions.iter()) {
        assert!(
            (light.mz - heavy.mz).abs() < 1e-6,
            "zero-offset: light ({}) and heavy ({}) m/z should be identical for {}",
            light.mz,
            heavy.mz,
            light.label
        );
    }
}

#[test]
fn nonzero_offset_heavy_ions_shifted() {
    let peptide = "PEPTIDEK"; // has K
    let label = silac_label();

    let delta = total_heavy_delta(peptide, &label);
    assert!(delta > 1.0, "PEPTIDEK should have ~8 Da delta");

    let light_ions = build_target_ions(peptide, &[], 2);
    let heavy_ions = compute_heavy_target_ions(&light_ions, peptide, &label);

    assert_eq!(light_ions.len(), heavy_ions.len());

    // At least some ions should be shifted (those covering the K residue)
    let shifted_count = light_ions
        .iter()
        .zip(heavy_ions.iter())
        .filter(|(l, h)| (l.mz - h.mz).abs() > 0.01)
        .count();
    assert!(
        shifted_count > 0,
        "at least some ions should be shifted for PEPTIDEK"
    );
}
