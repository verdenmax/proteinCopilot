//! SILAC heavy-label m/z calculation for XIC traces.
//!
//! Computes heavy-label m/z shifts for fragment ions and precursors
//! based on the number of K (Lysine) and R (Arginine) residues in
//! each fragment.

use crate::extract::TargetIon;
use crate::IonType;
use protein_copilot_core::label::{residue_heavy_delta, total_heavy_delta, LabelType};
use protein_copilot_core::spectrum::{IsolationWindow, Spectrum, SpectrumRepresentation};

/// Compute heavy-label m/z for a precursor ion.
///
/// Delegates to [`protein_copilot_core::label::compute_heavy_precursor_mz`].
pub fn compute_heavy_precursor_mz(
    light_mz: f64,
    charge: i32,
    peptide_sequence: &str,
    label: &LabelType,
) -> f64 {
    protein_copilot_core::label::compute_heavy_precursor_mz(
        light_mz,
        charge,
        peptide_sequence,
        label,
    )
}

/// Compute heavy-label target ions from light target ions.
///
/// For each light fragment ion, counts the K/R residues in the
/// corresponding fragment (prefix for b-ions, suffix for y-ions)
/// and applies the mass shift divided by the fragment charge.
pub fn compute_heavy_target_ions(
    light_ions: &[TargetIon],
    peptide_sequence: &str,
    label: &LabelType,
) -> Vec<TargetIon> {
    let chars: Vec<char> = peptide_sequence.chars().collect();
    let n = chars.len();

    light_ions
        .iter()
        .map(|ion| {
            let fragment_delta = match ion.ion_type {
                IonType::B => {
                    // b-ion covers prefix[0..ion_number]
                    let prefix = &chars[..((ion.ion_number as usize).min(n))];
                    residue_heavy_delta(prefix, label)
                }
                IonType::Y => {
                    // y-ion covers suffix[n - ion_number..]
                    let start = n.saturating_sub(ion.ion_number as usize);
                    let suffix = &chars[start..];
                    residue_heavy_delta(suffix, label)
                }
                IonType::Precursor => total_heavy_delta(peptide_sequence, label),
            };

            let heavy_mz = ion.mz + fragment_delta / ion.charge.max(1) as f64;

            TargetIon {
                label: ion.label.clone(),
                ion_type: ion.ion_type,
                ion_number: ion.ion_number,
                charge: ion.charge,
                mz: heavy_mz,
            }
        })
        .collect()
}

/// Check if an isolation window contains a given m/z value.
///
/// Returns `true` when `mz` falls within `[target_mz - lower_offset, target_mz + upper_offset]`.
pub fn window_contains_mz(window: &IsolationWindow, mz: f64) -> bool {
    let lo = window.target_mz - window.lower_offset;
    let hi = window.target_mz + window.upper_offset;
    mz >= lo && mz <= hi
}

/// Search DDA MS2 spectra for a scan whose precursor m/z matches the heavy target.
///
/// In DDA, each MS2 scan selects a specific precursor. The heavy peptide has a
/// different precursor m/z than the light one, so it will be in a different scan.
/// Finds the MS2 scan closest to `reference_rt_min` whose precursor m/z is within
/// `tolerance_ppm` of `target_mz`. Returns `scan_number` if found.
pub fn find_dda_heavy_scan(
    spectra: &[Spectrum],
    reference_rt_min: f64,
    target_mz: f64,
    tolerance_ppm: f64,
) -> Option<u32> {
    if spectra.is_empty() || tolerance_ppm <= 0.0 || !target_mz.is_finite() || target_mz <= 0.0 {
        return None;
    }
    spectra
        .iter()
        .filter_map(|spec| {
            let prec_mz = spec.precursors.first()?.mz;
            let ppm_err = ((prec_mz - target_mz) / target_mz * 1e6).abs();
            if ppm_err < tolerance_ppm {
                let rt_delta = (spec.retention_time_min - reference_rt_min).abs();
                Some((spec.scan_number, rt_delta))
            } else {
                None
            }
        })
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(scan, _)| scan)
}

/// Search a slice of MS2 spectra for the isolation window covering `target_mz`.
///
/// Finds the MS2 scan closest to `reference_rt_min` whose isolation window
/// contains `target_mz`. Returns `(scan_number, isolation_window)` if found.
pub fn find_heavy_dia_window_from_spectra(
    spectra: &[Spectrum],
    reference_rt_min: f64,
    target_mz: f64,
) -> Option<(u32, IsolationWindow)> {
    if spectra.is_empty() || !target_mz.is_finite() || target_mz <= 0.0 {
        return None;
    }
    spectra
        .iter()
        .filter_map(|spec| {
            let win = spec
                .precursors
                .first()
                .and_then(|p| p.isolation_window.as_ref())?;
            if window_contains_mz(win, target_mz) {
                let rt_delta = (spec.retention_time_min - reference_rt_min).abs();
                Some((spec.scan_number, win.clone(), rt_delta))
            } else {
                None
            }
        })
        .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(scan, win, _)| (scan, win))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn silac() -> LabelType {
        LabelType::standard_silac()
    }

    #[test]
    fn heavy_precursor_single_k() {
        // "PEPTIDEK" has 1 K, 0 R
        let heavy = compute_heavy_precursor_mz(500.0, 2, "PEPTIDEK", &silac());
        let expected = 500.0 + 8.014199 / 2.0;
        assert!(
            (heavy - expected).abs() < 1e-4,
            "got {heavy}, expected {expected}"
        );
    }

    #[test]
    fn heavy_precursor_k_and_r() {
        // "PEPTIDER" has 0 K, 1 R
        let heavy = compute_heavy_precursor_mz(500.0, 2, "PEPTIDER", &silac());
        let expected = 500.0 + 10.008269 / 2.0;
        assert!((heavy - expected).abs() < 1e-4);
    }

    #[test]
    fn heavy_precursor_no_label_residues() {
        // "PEPTIDE" has 0 K, 0 R
        let heavy = compute_heavy_precursor_mz(500.0, 2, "PEPTIDE", &silac());
        assert!((heavy - 500.0).abs() < 1e-6, "no K/R means no shift");
    }

    #[test]
    fn heavy_fragment_b_ion_prefix() {
        // "KPEPTIDE" — b1 covers "K" (1 K), b2 covers "KP" (1 K)
        let light = vec![
            TargetIon {
                label: "b1¹⁺".to_string(),
                ion_type: IonType::B,
                ion_number: 1,
                charge: 1,
                mz: 100.0,
            },
            TargetIon {
                label: "b2¹⁺".to_string(),
                ion_type: IonType::B,
                ion_number: 2,
                charge: 1,
                mz: 200.0,
            },
        ];
        let heavy = compute_heavy_target_ions(&light, "KPEPTIDE", &silac());
        assert_eq!(heavy.len(), 2);
        // b1 includes K → +8.014
        assert!((heavy[0].mz - (100.0 + 8.014199)).abs() < 1e-4);
        // b2 includes K → +8.014
        assert!((heavy[1].mz - (200.0 + 8.014199)).abs() < 1e-4);
    }

    #[test]
    fn heavy_fragment_y_ion_suffix() {
        // "PEPTIDEK" — y1 covers "K" (1 K), y2 covers "EK" (1 K)
        let light = vec![
            TargetIon {
                label: "y1¹⁺".to_string(),
                ion_type: IonType::Y,
                ion_number: 1,
                charge: 1,
                mz: 100.0,
            },
            TargetIon {
                label: "y2¹⁺".to_string(),
                ion_type: IonType::Y,
                ion_number: 2,
                charge: 1,
                mz: 200.0,
            },
        ];
        let heavy = compute_heavy_target_ions(&light, "PEPTIDEK", &silac());
        // y1 suffix = "K" → +8.014
        assert!((heavy[0].mz - (100.0 + 8.014199)).abs() < 1e-4);
        // y2 suffix = "EK" → +8.014
        assert!((heavy[1].mz - (200.0 + 8.014199)).abs() < 1e-4);
    }

    #[test]
    fn heavy_fragment_doubly_charged() {
        // b1 at charge 2: delta should be halved
        let light = vec![TargetIon {
            label: "b1²⁺".to_string(),
            ion_type: IonType::B,
            ion_number: 1,
            charge: 2,
            mz: 100.0,
        }];
        let heavy = compute_heavy_target_ions(&light, "KPEPTIDE", &silac());
        let expected = 100.0 + 8.014199 / 2.0;
        assert!((heavy[0].mz - expected).abs() < 1e-4);
    }

    #[test]
    fn window_contains_mz_inside() {
        let w = IsolationWindow {
            target_mz: 500.0,
            lower_offset: 12.5,
            upper_offset: 12.5,
        };
        assert!(window_contains_mz(&w, 500.0));
        assert!(window_contains_mz(&w, 487.5)); // exactly at lower edge
        assert!(window_contains_mz(&w, 512.5)); // exactly at upper edge
        assert!(window_contains_mz(&w, 505.0)); // inside
    }

    #[test]
    fn window_contains_mz_outside() {
        let w = IsolationWindow {
            target_mz: 500.0,
            lower_offset: 12.5,
            upper_offset: 12.5,
        };
        assert!(!window_contains_mz(&w, 487.4)); // just below
        assert!(!window_contains_mz(&w, 512.6)); // just above
        assert!(!window_contains_mz(&w, 600.0)); // far away
    }

    #[test]
    fn find_heavy_dia_window_basic() {
        let spectra = vec![
            make_ms2_spec(100, 10.0, 400.0, 12.5),
            make_ms2_spec(101, 10.01, 425.0, 12.5),
            make_ms2_spec(102, 10.02, 450.0, 12.5),
        ];

        // heavy m/z = 440.0 should match window #3 (450 ±12.5 = [437.5, 462.5])
        let result = find_heavy_dia_window_from_spectra(&spectra, 10.0, 440.0);
        assert!(result.is_some(), "should find window containing 440.0");
        let (scan, win) = result.unwrap();
        assert_eq!(scan, 102);
        assert!((win.target_mz - 450.0).abs() < 0.01);
    }

    #[test]
    fn find_heavy_dia_window_not_found() {
        let spectra = vec![
            make_ms2_spec(100, 10.0, 400.0, 12.5),
            make_ms2_spec(101, 10.01, 425.0, 12.5),
        ];
        let result = find_heavy_dia_window_from_spectra(&spectra, 10.0, 600.0);
        assert!(result.is_none());
    }

    #[test]
    fn find_dda_heavy_scan_basic() {
        // DDA scans with narrow windows, different precursor m/z values
        let spectra = vec![
            make_dda_spec(100, 10.0, 500.0),   // light
            make_dda_spec(101, 10.01, 504.01), // heavy (K+8 at charge 2 ≈ +4.007)
            make_dda_spec(102, 10.02, 600.0),  // unrelated
        ];
        // heavy target m/z = 504.007 (light 500 + 8.014/2)
        let result = find_dda_heavy_scan(&spectra, 10.0, 504.007, 20.0);
        assert!(result.is_some());
        assert_eq!(result.unwrap(), 101);
    }

    #[test]
    fn find_dda_heavy_scan_picks_closest_rt() {
        let spectra = vec![
            make_dda_spec(100, 10.0, 504.01), // matches but far in RT
            make_dda_spec(200, 15.0, 504.01), // matches and closer to ref RT=14.0
            make_dda_spec(300, 20.0, 504.01), // matches but farther
        ];
        let result = find_dda_heavy_scan(&spectra, 14.0, 504.007, 20.0);
        assert_eq!(result.unwrap(), 200);
    }

    #[test]
    fn find_dda_heavy_scan_not_found() {
        let spectra = vec![
            make_dda_spec(100, 10.0, 500.0),
            make_dda_spec(101, 10.01, 510.0), // too far from target
        ];
        let result = find_dda_heavy_scan(&spectra, 10.0, 504.007, 20.0);
        assert!(result.is_none());
    }

    fn make_dda_spec(scan: u32, rt_min: f64, precursor_mz: f64) -> Spectrum {
        use protein_copilot_core::spectrum::{MsLevel, PrecursorInfo};
        Spectrum {
            scan_number: scan,
            ms_level: MsLevel::MS2,
            retention_time_min: rt_min,
            mz_array: vec![],
            intensity_array: vec![],
            precursors: vec![PrecursorInfo {
                mz: precursor_mz,
                charge: Some(2),
                intensity: None,
                source_scan: None,
                isolation_window: Some(IsolationWindow {
                    target_mz: precursor_mz,
                    lower_offset: 0.35,
                    upper_offset: 0.35,
                }),
            }],
            representation: SpectrumRepresentation::Centroid,
        }
    }

    fn make_ms2_spec(scan: u32, rt_min: f64, center_mz: f64, half_width: f64) -> Spectrum {
        use protein_copilot_core::spectrum::{MsLevel, PrecursorInfo};
        Spectrum {
            scan_number: scan,
            ms_level: MsLevel::MS2,
            retention_time_min: rt_min,
            mz_array: vec![],
            intensity_array: vec![],
            precursors: vec![PrecursorInfo {
                mz: center_mz,
                charge: Some(2),
                intensity: None,
                source_scan: None,
                isolation_window: Some(IsolationWindow {
                    target_mz: center_mz,
                    lower_offset: half_width,
                    upper_offset: half_width,
                }),
            }],
            representation: SpectrumRepresentation::Centroid,
        }
    }
}
