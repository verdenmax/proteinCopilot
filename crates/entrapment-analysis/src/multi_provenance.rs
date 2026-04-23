//! Multi-Target Fragment Ion Provenance Engine
//!
//! Extends the single-target provenance engine (v3) to match each observed
//! MS2 peak against theoretical b/y ions from *N* co-eluting target peptides.
//! Supports SILAC heavy-label ion shifting for Heavy candidates.

use protein_copilot_core::search_params::MassTolerance;
use protein_copilot_search_engine::matching::within_tolerance;

use crate::provenance::{generate_theoretical_ions, TheoreticalIon};
use crate::types::{
    CoElutingCandidate, LabelForm, MirrorData, MultiAnnotatedPeak, TargetIonMatch,
};

// ---------------------------------------------------------------------------
// Core public function
// ---------------------------------------------------------------------------

/// Trace fragment ion provenance against multiple co-eluting targets.
///
/// For each observed peak, matches against theoretical ions from the trap
/// peptide and each co-eluting candidate (light or heavy), then classifies
/// the peak as trap-only, target-only, shared, or unassigned.
///
/// # Arguments
///
/// * `observed_mz` — m/z values of observed MS2 peaks.
/// * `observed_intensity` — intensities corresponding to each observed peak.
/// * `trap_sequence` — amino-acid sequence of the trap peptide.
/// * `trap_modifications` — position-delta pairs applied to the trap sequence.
/// * `candidates` — co-eluting target peptide candidates.
/// * `fragment_tolerance` — mass tolerance for matching.
/// * `max_fragment_charge` — maximum charge state for theoretical ions.
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

        // Match against trap ions.
        let trap_match = find_best_match(mz, &trap_ions, fragment_tolerance);
        let trap_ion = trap_match.map(|(label, _ppm)| label);

        // Match against each candidate's ions.
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

        // Classify.
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
        scan_number: 0, // caller sets the actual scan number
        annotated_peaks,
        trap_only_count,
        target_only_count,
        shared_count,
        unassigned_count,
    }
}

// ---------------------------------------------------------------------------
// Mirror with pre-generated trap ions
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Generate theoretical ions for a single co-eluting candidate.
///
/// For `Light` candidates, delegates directly to `generate_theoretical_ions`.
/// For `Heavy` candidates, generates light ions first and then shifts them
/// by the cumulative heavy-residue delta masses.
fn generate_candidate_ions(
    candidate: &CoElutingCandidate,
    max_charge: i32,
) -> Vec<TheoreticalIon> {
    match &candidate.label_form {
        LabelForm::Light => {
            generate_theoretical_ions(&candidate.peptide, &candidate.modifications, max_charge)
        }
        LabelForm::Heavy {
            residue_deltas, ..
        } => {
            let light_ions = generate_theoretical_ions(
                &candidate.peptide,
                &candidate.modifications,
                max_charge,
            );
            shift_ions_heavy(&candidate.peptide, &light_ions, residue_deltas, max_charge)
        }
    }
}

/// Shift light theoretical ions to heavy by adding cumulative residue deltas.
///
/// For b_n ions the delta comes from residues `[0..n]`;
/// for y_n ions the delta comes from residues `[len-n..len]`.
/// Appends `"(H)"` to each ion label to mark it as heavy.
pub(crate) fn shift_ions_heavy(
    sequence: &str,
    light_ions: &[TheoreticalIon],
    residue_deltas: &[(usize, f64)],
    _max_charge: i32,
) -> Vec<TheoreticalIon> {
    let seq_len = sequence.chars().count();

    light_ions
        .iter()
        .map(|ion| {
            let (ion_type, ion_number, charge) = parse_ion_label(&ion.label);

            // Compute cumulative delta for the residues covered by this ion.
            let cumulative_delta: f64 = residue_deltas
                .iter()
                .filter(|&&(pos, _)| match ion_type {
                    'b' => pos < ion_number,
                    'y' => pos >= seq_len.saturating_sub(ion_number),
                    _ => false,
                })
                .map(|&(_, delta)| delta)
                .sum();

            // Shift the m/z: delta is in Da, divide by charge.
            let shifted_mz = ion.mz + cumulative_delta / charge as f64;

            TheoreticalIon {
                mz: shifted_mz,
                label: format!("{}{}+{}(H)", ion_type, ion_number, charge),
            }
        })
        .collect()
}

/// Parse an ion label like `"b3+1"` or `"y5+2"` into `(type, number, charge)`.
///
/// Returns `('b', 3, 1)` for `"b3+1"` and `('y', 5, 2)` for `"y5+2"`.
fn parse_ion_label(label: &str) -> (char, usize, i32) {
    let chars: Vec<char> = label.chars().collect();
    let ion_type = chars[0]; // 'b' or 'y'

    // Find the '+' separator.
    let plus_pos = label.find('+').unwrap_or(label.len());
    let number: usize = label[1..plus_pos].parse().unwrap_or(0);
    let charge: i32 = if plus_pos < label.len() {
        // Take everything after '+', stripping any suffix like "(H)"
        let charge_str: String = label[plus_pos + 1..]
            .chars()
            .take_while(|c| c.is_ascii_digit())
            .collect();
        charge_str.parse().unwrap_or(1)
    } else {
        1
    };

    (ion_type, number, charge)
}

/// Find the best matching theoretical ion within tolerance and return its
/// label and the mass error in ppm.
fn find_best_match_with_ppm(
    observed_mz: f64,
    ions: &[TheoreticalIon],
    tolerance: &MassTolerance,
) -> Option<(String, f64)> {
    let mut best: Option<(f64, &str, f64)> = None; // (abs_error, label, ppm)

    for ion in ions {
        if within_tolerance(observed_mz, ion.mz, tolerance) {
            let error = (observed_mz - ion.mz).abs();
            let ppm = (observed_mz - ion.mz) / ion.mz * 1e6;
            match best {
                Some((best_err, _, _)) if error < best_err => {
                    best = Some((error, &ion.label, ppm));
                }
                None => {
                    best = Some((error, &ion.label, ppm));
                }
                _ => {}
            }
        }
    }

    best.map(|(_, label, ppm)| (label.to_string(), ppm))
}

/// Find the best matching theoretical ion within tolerance.
///
/// Convenience wrapper around [`find_best_match_with_ppm`].
fn find_best_match(
    observed_mz: f64,
    ions: &[TheoreticalIon],
    tolerance: &MassTolerance,
) -> Option<(String, f64)> {
    find_best_match_with_ppm(observed_mz, ions, tolerance)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::search_params::ToleranceUnit;

    fn tolerance_20ppm() -> MassTolerance {
        MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        }
    }

    #[test]
    fn test_parse_ion_label_b() {
        let (ion_type, number, charge) = parse_ion_label("b3+1");
        assert_eq!(ion_type, 'b');
        assert_eq!(number, 3);
        assert_eq!(charge, 1);
    }

    #[test]
    fn test_parse_ion_label_y() {
        let (ion_type, number, charge) = parse_ion_label("y5+2");
        assert_eq!(ion_type, 'y');
        assert_eq!(number, 5);
        assert_eq!(charge, 2);
    }

    #[test]
    fn test_parse_ion_label_heavy() {
        let (ion_type, number, charge) = parse_ion_label("b3+1(H)");
        assert_eq!(ion_type, 'b');
        assert_eq!(number, 3);
        assert_eq!(charge, 1);
    }

    #[test]
    fn test_match_single_target_light() {
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

        // Use trap theoretical ions as "observed" peaks.
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

        assert_eq!(result.annotated_peaks.len(), observed_mz.len());
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

        // Heavy y1 for K: (128.09496 + 8.014199 + 18.010565 + 1.007276) / 1 = 155.127
        let observed_mz = vec![155.127];
        let observed_int = vec![1000.0];

        let result = trace_multi_target(
            &observed_mz,
            &observed_int,
            "ABCDEFGH",
            &[],
            &candidates,
            &MassTolerance {
                value: 50.0,
                unit: ToleranceUnit::Ppm,
            },
            1,
        );

        assert_eq!(result.target_only_count, 1);
        assert_eq!(result.annotated_peaks[0].target_matches.len(), 1);
    }

    #[test]
    fn test_empty_observed_peaks() {
        let candidates = vec![CoElutingCandidate {
            peptide: "PEPTIDE".to_string(),
            protein_ids: vec!["P1".to_string()],
            precursor_mz: 400.0,
            charge: 2,
            rt_start: 30.0,
            rt_stop: 32.0,
            label_form: LabelForm::Light,
            modifications: vec![],
        }];

        let result = trace_multi_target(
            &[],
            &[],
            "STTTG",
            &[],
            &candidates,
            &tolerance_20ppm(),
            1,
        );

        assert_eq!(result.annotated_peaks.len(), 0);
        assert_eq!(result.trap_only_count, 0);
        assert_eq!(result.target_only_count, 0);
        assert_eq!(result.shared_count, 0);
        assert_eq!(result.unassigned_count, 0);
    }

    #[test]
    fn test_no_candidates() {
        let trap_seq = "STTTG";
        let trap_ions = generate_theoretical_ions(trap_seq, &[], 1);
        let observed_mz: Vec<f64> = trap_ions.iter().map(|i| i.mz).collect();
        let observed_int: Vec<f64> = vec![1000.0; observed_mz.len()];

        let result = trace_multi_target(
            &observed_mz,
            &observed_int,
            trap_seq,
            &[],
            &[],
            &tolerance_20ppm(),
            1,
        );

        // All peaks should be trap-only since there are no candidates.
        assert!(result.trap_only_count > 0);
        assert_eq!(result.target_only_count, 0);
        assert_eq!(result.shared_count, 0);
    }

    #[test]
    fn test_scan_number_is_zero() {
        let result = trace_multi_target(
            &[100.0],
            &[1000.0],
            "AG",
            &[],
            &[],
            &tolerance_20ppm(),
            1,
        );
        assert_eq!(result.scan_number, 0);
    }

    #[test]
    fn test_shift_ions_heavy_b_ion() {
        // Sequence "AK" with heavy K at position 1: delta = 8.014199
        // Light b1 for 'A' = 71.03711 + 1.007276 = 72.044386
        // b1 covers residues [0..1] = [A], K is at pos 1 → not covered → no shift
        let light_ions = generate_theoretical_ions("AK", &[], 1);
        let heavy_ions = shift_ions_heavy("AK", &light_ions, &[(1, 8.014199)], 1);

        // b1 should not be shifted (only covers pos 0)
        let b1_heavy = heavy_ions
            .iter()
            .find(|i| i.label.starts_with("b1+1"))
            .expect("b1+1(H)");
        let b1_light = light_ions
            .iter()
            .find(|i| i.label == "b1+1")
            .expect("b1+1");
        assert!(
            (b1_heavy.mz - b1_light.mz).abs() < 0.001,
            "b1 should not be shifted; heavy={}, light={}",
            b1_heavy.mz,
            b1_light.mz
        );

        // y1 covers residue [1] = K → should be shifted by 8.014199
        let y1_heavy = heavy_ions
            .iter()
            .find(|i| i.label.starts_with("y1+1"))
            .expect("y1+1(H)");
        let y1_light = light_ions
            .iter()
            .find(|i| i.label == "y1+1")
            .expect("y1+1");
        let expected_shift = 8.014199;
        assert!(
            (y1_heavy.mz - y1_light.mz - expected_shift).abs() < 0.001,
            "y1 should be shifted by ~8.014; heavy={}, light={}, diff={}",
            y1_heavy.mz,
            y1_light.mz,
            y1_heavy.mz - y1_light.mz
        );
    }

    #[test]
    fn test_heavy_label_marks_h_suffix() {
        let light_ions = generate_theoretical_ions("AK", &[], 1);
        let heavy_ions = shift_ions_heavy("AK", &light_ions, &[(1, 8.014199)], 1);

        for ion in &heavy_ions {
            assert!(
                ion.label.contains("(H)"),
                "Heavy ion label '{}' should contain '(H)'",
                ion.label
            );
        }
    }
}
