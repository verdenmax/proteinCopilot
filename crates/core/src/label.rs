//! SILAC heavy-label types and mass delta computation.
//!
//! Shared between `xic` (XIC extraction) and `search-engine` (annotation)
//! crates to avoid circular dependencies.

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

/// Heavy-label type for SILAC or custom isotope labeling.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub enum LabelType {
    /// SILAC heavy amino acids.
    Silac {
        /// Mass shift for heavy Lysine (default: 8.014199 Da for ¹³C₆¹⁵N₂-Lys).
        heavy_k_delta: f64,
        /// Mass shift for heavy Arginine (default: 10.008269 Da for ¹³C₆¹⁵N₄-Arg).
        heavy_r_delta: f64,
    },
    /// Custom residue mass shifts.
    Custom {
        /// (residue, mass_delta) pairs.
        residue_deltas: Vec<(char, f64)>,
    },
}

impl LabelType {
    /// Standard SILAC: K+8, R+10.
    pub fn standard_silac() -> Self {
        LabelType::Silac {
            heavy_k_delta: 8.014199,
            heavy_r_delta: 10.008269,
        }
    }
}

/// Heavy mass delta for a slice of residues given a label type.
pub fn residue_heavy_delta(residues: &[char], label: &LabelType) -> f64 {
    match label {
        LabelType::Silac {
            heavy_k_delta,
            heavy_r_delta,
        } => {
            let count_k = residues.iter().filter(|&&c| c == 'K' || c == 'k').count() as f64;
            let count_r = residues.iter().filter(|&&c| c == 'R' || c == 'r').count() as f64;
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

/// Total heavy mass delta for an entire peptide sequence.
pub fn total_heavy_delta(peptide_sequence: &str, label: &LabelType) -> f64 {
    let chars: Vec<char> = peptide_sequence.chars().collect();
    residue_heavy_delta(&chars, label)
}

/// Compute heavy precursor m/z from light precursor m/z.
pub fn compute_heavy_precursor_mz(
    light_mz: f64,
    charge: i32,
    peptide_sequence: &str,
    label: &LabelType,
) -> f64 {
    let delta = total_heavy_delta(peptide_sequence, label);
    light_mz + delta / charge.abs().max(1) as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_silac() {
        let label = LabelType::standard_silac();
        match &label {
            LabelType::Silac {
                heavy_k_delta,
                heavy_r_delta,
            } => {
                assert!((heavy_k_delta - 8.014199).abs() < 1e-6);
                assert!((heavy_r_delta - 10.008269).abs() < 1e-6);
            }
            _ => panic!("expected Silac variant"),
        }
    }

    #[test]
    fn test_residue_heavy_delta_silac() {
        let label = LabelType::standard_silac();
        let chars: Vec<char> = "KR".chars().collect();
        let delta = residue_heavy_delta(&chars, &label);
        assert!((delta - 18.022468).abs() < 1e-4);
    }

    #[test]
    fn test_no_kr_no_delta() {
        let label = LabelType::standard_silac();
        let delta = total_heavy_delta("AGLDEF", &label);
        assert!((delta - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_heavy_precursor_mz() {
        let label = LabelType::standard_silac();
        // DGFLLDGFPR: 1 R → delta = 10.008269, charge 2 → +5.004
        let heavy = compute_heavy_precursor_mz(569.3, 2, "DGFLLDGFPR", &label);
        assert!((heavy - 569.3 - 10.008269 / 2.0).abs() < 1e-4);
    }
}
