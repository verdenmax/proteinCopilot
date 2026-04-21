//! Fragment Ion Provenance Engine
//!
//! Classifies each observed MS2 peak as originating from the trap peptide,
//! the target peptide, both (shared), or neither (unassigned) by matching
//! against theoretical b/y ions.

use protein_copilot_core::search_params::MassTolerance;
use protein_copilot_search_engine::matching::within_tolerance;
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

const PROTON_MASS: f64 = 1.007276;
const WATER_MASS: f64 = 18.010565;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Classification of a single observed MS2 peak's origin.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum IonOrigin {
    /// Peak matches only theoretical ions from the trap peptide.
    TrapOnly,
    /// Peak matches only theoretical ions from the target peptide.
    TargetOnly,
    /// Peak matches theoretical ions from both trap and target peptides.
    Shared,
    /// Peak does not match any theoretical ion.
    Unassigned,
}

/// A single annotated peak with its provenance classification.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnnotatedPeak {
    /// Observed m/z value.
    pub mz_observed: f64,
    /// Observed intensity.
    pub intensity: f64,
    /// Provenance classification.
    pub origin: IonOrigin,
    /// Ion label from the trap peptide, e.g. "b3+1", "y5+2".
    pub trap_ion_label: Option<String>,
    /// Ion label from the target peptide, e.g. "b4+1", "y7+1".
    pub target_ion_label: Option<String>,
}

/// Complete provenance analysis result for one PSM.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FragmentProvenance {
    /// Trap peptide sequence.
    pub trap_sequence: String,
    /// Target peptide sequence (empty if L4 — no target match).
    pub target_sequence: String,
    /// Per-peak annotation results.
    pub annotated_peaks: Vec<AnnotatedPeak>,
    /// Number of peaks matching only trap ions.
    pub trap_matched_count: u32,
    /// Number of peaks matching only target ions.
    pub target_matched_count: u32,
    /// Number of peaks matching both trap and target ions.
    pub shared_count: u32,
    /// Number of peaks matching neither.
    pub unassigned_count: u32,
    /// shared / (trap_matched + target_matched + shared), 0.0 if denominator is 0.
    pub shared_ratio: f64,
    /// Whether `shared_ratio` exceeds the chimera threshold (set by caller).
    pub is_chimeric: bool,
}

// ---------------------------------------------------------------------------
// Internal types
// ---------------------------------------------------------------------------

/// A theoretical ion with its m/z and human-readable label.
struct TheoreticalIon {
    mz: f64,
    label: String,
}

// ---------------------------------------------------------------------------
// Core public function
// ---------------------------------------------------------------------------

/// Trace fragment ion provenance for a single PSM.
///
/// Given observed MS2 peaks, generates theoretical b/y ions for both trap and
/// target peptides, then classifies each observed peak as `TrapOnly`,
/// `TargetOnly`, `Shared`, or `Unassigned`.
///
/// # Arguments
///
/// * `observed_mz` — m/z values of observed MS2 peaks.
/// * `observed_intensity` — intensities corresponding to each observed peak.
/// * `trap_sequence` — amino-acid sequence of the trap peptide.
/// * `target_sequence` — amino-acid sequence of the target peptide (empty for L4).
/// * `trap_modifications` — position-delta pairs applied to the trap sequence.
/// * `fragment_tolerance` — mass tolerance for matching.
/// * `max_fragment_charge` — maximum charge state for theoretical ions.
pub fn trace_provenance(
    observed_mz: &[f64],
    observed_intensity: &[f64],
    trap_sequence: &str,
    target_sequence: &str,
    trap_modifications: &[(usize, f64)],
    fragment_tolerance: &MassTolerance,
    max_fragment_charge: i32,
) -> FragmentProvenance {
    // 1. Generate theoretical ions for the trap peptide.
    let trap_ions =
        generate_theoretical_ions(trap_sequence, trap_modifications, max_fragment_charge);

    // 2. Generate theoretical ions for the target peptide (empty if L4).
    let target_ions = if target_sequence.is_empty() {
        Vec::new()
    } else {
        generate_theoretical_ions(target_sequence, &[], max_fragment_charge)
    };

    // 3. Classify each observed peak.
    let mut annotated_peaks = Vec::with_capacity(observed_mz.len());
    let mut trap_matched = 0u32;
    let mut target_matched = 0u32;
    let mut shared = 0u32;
    let mut unassigned = 0u32;

    for (i, &mz) in observed_mz.iter().enumerate() {
        let intensity = observed_intensity.get(i).copied().unwrap_or(0.0);

        let trap_match = find_matching_ion(mz, &trap_ions, fragment_tolerance);
        let target_match = find_matching_ion(mz, &target_ions, fragment_tolerance);

        let (origin, trap_label, target_label) = match (trap_match, target_match) {
            (Some(tl), Some(gl)) => {
                shared += 1;
                (IonOrigin::Shared, Some(tl), Some(gl))
            }
            (Some(tl), None) => {
                trap_matched += 1;
                (IonOrigin::TrapOnly, Some(tl), None)
            }
            (None, Some(gl)) => {
                target_matched += 1;
                (IonOrigin::TargetOnly, None, Some(gl))
            }
            (None, None) => {
                unassigned += 1;
                (IonOrigin::Unassigned, None, None)
            }
        };

        annotated_peaks.push(AnnotatedPeak {
            mz_observed: mz,
            intensity,
            origin,
            trap_ion_label: trap_label,
            target_ion_label: target_label,
        });
    }

    let total_matched = trap_matched + target_matched + shared;
    let shared_ratio = if total_matched == 0 {
        0.0
    } else {
        shared as f64 / total_matched as f64
    };

    FragmentProvenance {
        trap_sequence: trap_sequence.to_string(),
        target_sequence: target_sequence.to_string(),
        annotated_peaks,
        trap_matched_count: trap_matched,
        target_matched_count: target_matched,
        shared_count: shared,
        unassigned_count: unassigned,
        shared_ratio,
        is_chimeric: false, // caller sets based on config.chimera_threshold
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Generate all theoretical b/y ions for a peptide with optional modifications.
///
/// Modifications are position-delta pairs where position is a 0-based residue
/// index and delta is the mass shift in Daltons.
fn generate_theoretical_ions(
    sequence: &str,
    modifications: &[(usize, f64)],
    max_charge: i32,
) -> Vec<TheoreticalIon> {
    let residues: Vec<char> = sequence.chars().collect();
    let n = residues.len();

    if n == 0 {
        return Vec::new();
    }

    // Build residue masses with modifications applied.
    let mut residue_masses: Vec<f64> = residues.iter().map(|&c| amino_acid_mass(c)).collect();
    for &(pos, delta) in modifications {
        if pos < residue_masses.len() {
            residue_masses[pos] += delta;
        }
    }

    let mut ions = Vec::new();

    // b-ions: cumulative sum from N-terminus (b1 .. b_{n-1}).
    for charge in 1..=max_charge {
        let mut cumulative = 0.0;
        for (i, &mass) in residue_masses.iter().enumerate().take(n - 1) {
            cumulative += mass;
            let mz = (cumulative + PROTON_MASS * charge as f64) / charge as f64;
            ions.push(TheoreticalIon {
                mz,
                label: format!("b{}+{}", i + 1, charge),
            });
        }
    }

    // y-ions: cumulative sum from C-terminus (y1 .. y_{n-1}).
    for charge in 1..=max_charge {
        let mut cumulative = WATER_MASS;
        for i in (1..n).rev() {
            cumulative += residue_masses[i];
            let mz = (cumulative + PROTON_MASS * charge as f64) / charge as f64;
            ions.push(TheoreticalIon {
                mz,
                label: format!("y{}+{}", n - i, charge),
            });
        }
    }

    ions
}

/// Standard monoisotopic amino acid residue masses (Da).
fn amino_acid_mass(aa: char) -> f64 {
    match aa {
        'G' => 57.02146,
        'A' => 71.03711,
        'V' => 99.06841,
        'L' => 113.08406,
        'I' => 113.08406,
        'P' => 97.05276,
        'F' => 147.06841,
        'W' => 186.07931,
        'M' => 131.04049,
        'S' => 87.03203,
        'T' => 101.04768,
        'C' => 103.00919,
        'Y' => 163.06333,
        'H' => 137.05891,
        'D' => 115.02694,
        'E' => 129.04259,
        'N' => 114.04293,
        'Q' => 128.05858,
        'K' => 128.09496,
        'R' => 156.10111,
        _ => 0.0, // unknown amino acid
    }
}

/// Find the best matching theoretical ion for an observed m/z value.
///
/// Returns the label of the closest matching ion within tolerance, or `None`.
fn find_matching_ion(
    observed_mz: f64,
    theoretical_ions: &[TheoreticalIon],
    tolerance: &MassTolerance,
) -> Option<String> {
    let mut best_match: Option<(f64, &str)> = None;

    for ion in theoretical_ions {
        if within_tolerance(observed_mz, ion.mz, tolerance) {
            let error = (observed_mz - ion.mz).abs();
            match best_match {
                Some((best_err, _)) if error < best_err => {
                    best_match = Some((error, &ion.label));
                }
                None => {
                    best_match = Some((error, &ion.label));
                }
                _ => {}
            }
        }
    }

    best_match.map(|(_, label)| label.to_string())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};

    fn ppm_tolerance() -> MassTolerance {
        MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        }
    }

    #[test]
    fn test_amino_acid_mass_glycine() {
        assert!((amino_acid_mass('G') - 57.02146).abs() < 1e-4);
    }

    #[test]
    fn test_amino_acid_mass_unknown() {
        assert!((amino_acid_mass('X') - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_generate_theoretical_ions_simple() {
        let ions = generate_theoretical_ions("AG", &[], 1);
        // b1 = A mass + proton = 71.03711 + 1.007276 = 72.044386
        assert!(!ions.is_empty());
        let b1 = ions.iter().find(|i| i.label == "b1+1").expect("b1+1 ion");
        assert!((b1.mz - 72.044386).abs() < 0.01);
    }

    #[test]
    fn test_generate_theoretical_ions_with_mod() {
        // C with carbamidomethyl: 103.00919 + 57.021464 = 160.030654
        let mods = vec![(0, 57.021464)];
        let ions = generate_theoretical_ions("CA", &mods, 1);
        let b1 = ions.iter().find(|i| i.label == "b1+1").expect("b1+1 ion");
        // b1 = modified C mass + proton = 160.030654 + 1.007276 = 161.03793
        assert!((b1.mz - 161.03793).abs() < 0.01);
    }

    #[test]
    fn test_generate_theoretical_ions_empty_sequence() {
        let ions = generate_theoretical_ions("", &[], 1);
        assert!(ions.is_empty());
    }

    #[test]
    fn test_trace_provenance_all_trap_only() {
        // Trap=AGK, Target=FWR — completely different theoretical ions.
        // Observed peaks match only trap ions.
        let trap_ions = generate_theoretical_ions("AGK", &[], 1);
        let observed: Vec<f64> = trap_ions.iter().map(|i| i.mz).collect();
        let intensities: Vec<f64> = vec![100.0; observed.len()];

        let result = trace_provenance(
            &observed,
            &intensities,
            "AGK",
            "FWR",
            &[],
            &ppm_tolerance(),
            1,
        );

        assert!(result.trap_matched_count > 0);
        assert_eq!(result.target_matched_count, 0);
        assert_eq!(result.shared_count, 0);
        assert!((result.shared_ratio - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_trace_provenance_l4_no_target() {
        // L4: target_sequence is empty.
        let trap_ions = generate_theoretical_ions("AGK", &[], 1);
        let observed: Vec<f64> = trap_ions.iter().map(|i| i.mz).collect();
        let intensities: Vec<f64> = vec![100.0; observed.len()];

        let result = trace_provenance(
            &observed,
            &intensities,
            "AGK",
            "",
            &[],
            &ppm_tolerance(),
            1,
        );

        assert!(result.trap_matched_count > 0);
        assert_eq!(result.target_matched_count, 0);
        assert_eq!(result.shared_count, 0);
        assert!(result.target_sequence.is_empty());
    }

    #[test]
    fn test_trace_provenance_shared_ions() {
        // Same sequence for trap and target → all peaks should be Shared.
        let ions = generate_theoretical_ions("AGK", &[], 1);
        let observed: Vec<f64> = ions.iter().map(|i| i.mz).collect();
        let intensities: Vec<f64> = vec![100.0; observed.len()];

        let result = trace_provenance(
            &observed,
            &intensities,
            "AGK",
            "AGK",
            &[],
            &ppm_tolerance(),
            1,
        );

        assert_eq!(result.trap_matched_count, 0);
        assert_eq!(result.target_matched_count, 0);
        assert!(result.shared_count > 0);
        assert!((result.shared_ratio - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_trace_provenance_unassigned_peaks() {
        // Observed peaks that don't match any theoretical ion.
        let observed = vec![999.999, 1234.567];
        let intensities = vec![100.0, 200.0];

        let result = trace_provenance(
            &observed,
            &intensities,
            "AGK",
            "FWR",
            &[],
            &ppm_tolerance(),
            1,
        );

        assert_eq!(result.unassigned_count, 2);
        assert_eq!(result.trap_matched_count, 0);
        assert_eq!(result.target_matched_count, 0);
    }

    #[test]
    fn test_trace_provenance_empty_spectra() {
        let result = trace_provenance(&[], &[], "AGK", "FWR", &[], &ppm_tolerance(), 1);

        assert_eq!(result.trap_matched_count, 0);
        assert_eq!(result.target_matched_count, 0);
        assert_eq!(result.shared_count, 0);
        assert_eq!(result.unassigned_count, 0);
        assert!((result.shared_ratio - 0.0).abs() < 1e-10);
    }
}
