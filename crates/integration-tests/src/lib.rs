//! Shared test helpers for ProteinCopilot integration tests.
//!
//! Provides synthetic `Spectrum` builders and common constants for testing
//! annotation, XIC extraction, and search pipeline scenarios.

use protein_copilot_core::label::LabelType;
use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};
use protein_copilot_core::spectrum::{IsolationWindow, MsLevel, PrecursorInfo, Spectrum};

/// Standard SILAC label: ¹³C₆¹⁵N₂-Lys (+8.014 Da), ¹³C₆¹⁵N₄-Arg (+10.008 Da)
pub fn silac_label() -> LabelType {
    LabelType::Silac {
        heavy_k_delta: 8.014199,
        heavy_r_delta: 10.008269,
    }
}

/// Default fragment tolerance: 20 ppm
pub fn default_frag_tolerance() -> MassTolerance {
    MassTolerance {
        value: 20.0,
        unit: ToleranceUnit::Ppm,
    }
}

/// Build a DDA-style isolation window (narrow, 0.35 Th each side = 0.7 Th total)
pub fn dda_window(center_mz: f64) -> IsolationWindow {
    IsolationWindow {
        target_mz: center_mz,
        lower_offset: 0.35,
        upper_offset: 0.35,
    }
}

/// Build a DIA-style isolation window (wide, default 12.5 Th each side = 25 Th total)
pub fn dia_window(center_mz: f64) -> IsolationWindow {
    IsolationWindow {
        target_mz: center_mz,
        lower_offset: 12.5,
        upper_offset: 12.5,
    }
}

/// Build a synthetic MS2 spectrum with given peaks and precursor.
///
/// `scan_number`: 1-based scan number
/// `rt_min`: retention time in minutes
/// `precursor_mz`: selected precursor m/z
/// `charge`: precursor charge state
/// `isolation`: optional isolation window
/// `peaks`: list of (mz, intensity) pairs, will be sorted by m/z
pub fn make_ms2(
    scan_number: u32,
    rt_min: f64,
    precursor_mz: f64,
    charge: i32,
    isolation: Option<IsolationWindow>,
    peaks: Vec<(f64, f64)>,
) -> Spectrum {
    let mut sorted_peaks = peaks;
    sorted_peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let mz_array: Vec<f64> = sorted_peaks.iter().map(|p| p.0).collect();
    let intensity_array: Vec<f64> = sorted_peaks.iter().map(|p| p.1).collect();

    Spectrum {
        scan_number,
        ms_level: MsLevel::MS2,
        retention_time_min: rt_min,
        precursors: vec![PrecursorInfo {
            mz: precursor_mz,
            charge: Some(charge),
            intensity: Some(1e6),
            isolation_window: isolation,
            source_scan: None,
        }],
        mz_array,
        intensity_array,
    }
}

/// Build a synthetic MS1 spectrum with given peaks.
pub fn make_ms1(scan_number: u32, rt_min: f64, peaks: Vec<(f64, f64)>) -> Spectrum {
    let mut sorted_peaks = peaks;
    sorted_peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let mz_array: Vec<f64> = sorted_peaks.iter().map(|p| p.0).collect();
    let intensity_array: Vec<f64> = sorted_peaks.iter().map(|p| p.1).collect();

    Spectrum {
        scan_number,
        ms_level: MsLevel::MS1,
        retention_time_min: rt_min,
        precursors: vec![],
        mz_array,
        intensity_array,
    }
}

/// Ion entry: (ion_number, mz)
pub type IonEntry = (u32, f64);

/// Compute theoretical b/y ion m/z values for a peptide at charge 1+.
///
/// Returns `(b_ions, y_ions)` where each is a Vec of (ion_number, mz).
/// Uses monoisotopic residue masses. Does not include modifications.
pub fn theoretical_fragments(peptide: &str) -> (Vec<IonEntry>, Vec<IonEntry>) {
    use protein_copilot_search_engine::chemistry::{residue_mass, PROTON_MASS, WATER_MASS};

    let residues: Vec<char> = peptide.chars().collect();
    let n = residues.len();
    let mut b_ions = Vec::new();
    let mut y_ions = Vec::new();

    // b ions: cumulative mass from N-terminus + proton
    let mut cumulative = 0.0;
    for (i, &aa) in residues[..n - 1].iter().enumerate() {
        if let Some(mass) = residue_mass(aa) {
            cumulative += mass;
            let b_mz = cumulative + PROTON_MASS;
            b_ions.push(((i + 1) as u32, b_mz));
        }
    }

    // y ions: cumulative mass from C-terminus + water + proton
    let mut cumulative = 0.0;
    for i in (1..n).rev() {
        if let Some(mass) = residue_mass(residues[i]) {
            cumulative += mass;
            let y_mz = cumulative + WATER_MASS + PROTON_MASS;
            let y_num = (n - i) as u32;
            y_ions.push((y_num, y_mz));
        }
    }

    (b_ions, y_ions)
}

/// Generate a synthetic MS2 peak list containing some b/y ions for a peptide.
///
/// Returns peaks that include ~50% of theoretical fragments with realistic intensities,
/// plus some noise peaks. Useful for annotation tests.
pub fn synthetic_peaks_for_peptide(peptide: &str, base_intensity: f64) -> Vec<(f64, f64)> {
    let (b_ions, y_ions) = theoretical_fragments(peptide);
    let mut peaks = Vec::new();

    // Add ~half of b ions
    for (i, (_, mz)) in b_ions.iter().enumerate() {
        if i % 2 == 0 {
            peaks.push((
                *mz,
                base_intensity * (0.3 + 0.7 * (i as f64 / b_ions.len() as f64)),
            ));
        }
    }

    // Add ~half of y ions (y ions are typically stronger)
    for (i, (_, mz)) in y_ions.iter().enumerate() {
        if i % 2 == 0 {
            peaks.push((
                *mz,
                base_intensity * (0.5 + 0.5 * (i as f64 / y_ions.len() as f64)),
            ));
        }
    }

    // Add some noise peaks
    peaks.push((150.5, base_intensity * 0.05));
    peaks.push((250.3, base_intensity * 0.08));
    peaks.push((350.1, base_intensity * 0.03));

    peaks
}

/// Generate heavy-shifted peaks from light peaks using SILAC label.
///
/// Shifts b/y ion m/z values by the appropriate heavy delta for each ion position.
/// For simplicity, applies a uniform shift based on the peptide's total heavy delta
/// distributed across fragment positions.
pub fn heavy_shifted_peaks(
    peptide: &str,
    light_peaks: &[(f64, f64)],
    label: &LabelType,
) -> Vec<(f64, f64)> {
    let total_delta = protein_copilot_core::label::total_heavy_delta(peptide, label);
    let n_residues = peptide.len() as f64;

    // Simple approach: shift each peak proportionally
    light_peaks
        .iter()
        .map(|(mz, int)| {
            // Approximate: shift by fraction of total delta
            let shift = total_delta * (mz / (mz + 500.0)); // rough proxy
            if shift.abs() < 1e-6 {
                (*mz, *int * 0.8) // no shift for noise peaks
            } else {
                (mz + shift / n_residues.max(1.0), *int * 0.8)
            }
        })
        .collect()
}

/// Path to test fixtures in the search-engine crate.
pub fn search_engine_fixtures() -> std::path::PathBuf {
    let manifest = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    manifest
        .parent()
        .unwrap()
        .join("search-engine")
        .join("tests")
        .join("fixtures")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn theoretical_fragments_peptider() {
        let (b, y) = theoretical_fragments("PEPTIDE");
        assert_eq!(b.len(), 6, "PEPTIDE has 7 residues → 6 b ions");
        assert_eq!(y.len(), 6, "PEPTIDE has 7 residues → 6 y ions");
        // b1 = P mass + proton ≈ 97.053 + 1.007 ≈ 98.06
        assert!(
            (b[0].1 - 98.06).abs() < 0.1,
            "b1 should be ~98.06, got {}",
            b[0].1
        );
    }

    #[test]
    fn synthetic_peaks_non_empty() {
        let peaks = synthetic_peaks_for_peptide("PEPTIDEK", 1000.0);
        assert!(peaks.len() >= 5, "should have at least 5 peaks");
        assert!(
            peaks.iter().all(|(_, int)| *int > 0.0),
            "all intensities positive"
        );
    }

    #[test]
    fn silac_label_values() {
        let label = silac_label();
        match label {
            LabelType::Silac {
                heavy_k_delta,
                heavy_r_delta,
            } => {
                assert!((heavy_k_delta - 8.014199).abs() < 1e-4);
                assert!((heavy_r_delta - 10.008269).abs() < 1e-4);
            }
            _ => panic!("expected Silac variant"),
        }
    }

    #[test]
    fn dda_window_is_narrow() {
        let w = dda_window(500.0);
        assert!(
            (w.lower_offset + w.upper_offset) <= 1.0,
            "DDA window should be ≤ 1 Th"
        );
    }

    #[test]
    fn dia_window_is_wide() {
        let w = dia_window(500.0);
        assert!(
            (w.lower_offset + w.upper_offset) > 1.0,
            "DIA window should be > 1 Th"
        );
    }
}
