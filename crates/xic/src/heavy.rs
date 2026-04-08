//! SILAC heavy-label m/z calculation for XIC traces.
//!
//! Computes heavy-label m/z shifts for fragment ions and precursors
//! based on the number of K (Lysine) and R (Arginine) residues in
//! each fragment.

use crate::extract::TargetIon;
use crate::{IonType, LabelType};

/// Compute heavy-label m/z for a precursor ion.
///
/// Adds `(count_K × heavy_k_delta + count_R × heavy_r_delta) / charge`
/// to the light precursor m/z.
pub fn compute_heavy_precursor_mz(
    light_mz: f64,
    charge: i32,
    peptide_sequence: &str,
    label: &LabelType,
) -> f64 {
    let total_delta = total_heavy_delta(peptide_sequence, label);
    light_mz + total_delta / charge.abs().max(1) as f64
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

/// Total heavy mass delta for an entire peptide.
fn total_heavy_delta(peptide_sequence: &str, label: &LabelType) -> f64 {
    let chars: Vec<char> = peptide_sequence.chars().collect();
    residue_heavy_delta(&chars, label)
}

/// Heavy mass delta for a set of residues.
fn residue_heavy_delta(residues: &[char], label: &LabelType) -> f64 {
    match label {
        LabelType::Silac {
            heavy_k_delta,
            heavy_r_delta,
        } => {
            let count_k = residues.iter().filter(|&&c| c == 'K').count() as f64;
            let count_r = residues.iter().filter(|&&c| c == 'R').count() as f64;
            count_k * heavy_k_delta + count_r * heavy_r_delta
        }
        LabelType::Custom { residue_deltas } => {
            let mut delta = 0.0;
            for &(res, d) in residue_deltas {
                let count = residues.iter().filter(|&&c| c == res).count() as f64;
                delta += count * d;
            }
            delta
        }
    }
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
        assert!((heavy - expected).abs() < 1e-4, "got {heavy}, expected {expected}");
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
}
