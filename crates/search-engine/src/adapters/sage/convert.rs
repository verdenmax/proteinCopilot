//! Type conversion between ProteinCopilot and sage-core types.

use protein_copilot_core::search_params::{MassTolerance, ModPosition, Modification, ToleranceUnit};
use protein_copilot_core::spectrum::{MsLevel, Spectrum};
use sage_core::mass::Tolerance as SageTolerance;
use sage_core::modification::ModificationSpecificity;
use sage_core::spectrum::{Precursor as SagePrecursor, RawSpectrum, Representation};

/// Convert our Spectrum to sage RawSpectrum.
///
/// `file_id` is an index identifying which input file this spectrum came from.
pub fn spectrum_to_raw(spec: &Spectrum, file_id: usize) -> RawSpectrum {
    let ms_level = match spec.ms_level {
        MsLevel::MS1 => 1u8,
        MsLevel::MS2 => 2u8,
        MsLevel::Other(n) => n,
    };

    let precursors = spec
        .precursors
        .iter()
        .map(|p| {
            let charge = p.charge.map(|c| c.unsigned_abs() as u8);
            let isolation_window = p.isolation_window.as_ref().map(|iw| {
                let half_width = ((iw.lower_offset + iw.upper_offset) / 2.0) as f32;
                SageTolerance::Da(-half_width, half_width)
            });
            SagePrecursor {
                mz: p.mz as f32,
                intensity: p.intensity.map(|i| i as f32),
                charge,
                spectrum_ref: None,
                isolation_window,
                inverse_ion_mobility: None,
            }
        })
        .collect();

    RawSpectrum {
        file_id,
        ms_level,
        id: spec.scan_number.to_string(),
        precursors,
        representation: Representation::Centroid,
        scan_start_time: (spec.retention_time_min / 60.0) as f32,
        ion_injection_time: 0.0,
        total_ion_current: spec.intensity_array.iter().sum::<f64>() as f32,
        mz: spec.mz_array.iter().map(|&v| v as f32).collect(),
        intensity: spec.intensity_array.iter().map(|&v| v as f32).collect(),
        mobility: None,
    }
}

/// Convert our MassTolerance to sage Tolerance.
pub fn mass_tolerance_to_sage(tol: &MassTolerance) -> SageTolerance {
    let v = tol.value as f32;
    match tol.unit {
        ToleranceUnit::Ppm => SageTolerance::Ppm(-v, v),
        ToleranceUnit::Da => SageTolerance::Da(-v, v),
    }
}

/// Convert a single fixed modification to sage static_mods entries.
/// Returns one entry per target residue.
pub fn fixed_mod_to_sage(m: &Modification) -> Vec<(ModificationSpecificity, f32)> {
    let mass = m.mass_delta as f32;
    specificity_targets(m)
        .into_iter()
        .map(|spec| (spec, mass))
        .collect()
}

/// Convert a single variable modification to sage variable_mods entries.
/// Returns one entry per target residue.
pub fn variable_mod_to_sage(m: &Modification) -> Vec<(ModificationSpecificity, Vec<f32>)> {
    let mass = m.mass_delta as f32;
    specificity_targets(m)
        .into_iter()
        .map(|spec| (spec, vec![mass]))
        .collect()
}

/// Map our ModPosition + residues to one or more sage ModificationSpecificity values.
fn specificity_targets(m: &Modification) -> Vec<ModificationSpecificity> {
    match m.position {
        ModPosition::Anywhere => {
            m.residues
                .iter()
                .map(|&r| ModificationSpecificity::Residue(r as u8))
                .collect()
        }
        ModPosition::AnyNTerm => {
            if m.residues.is_empty() {
                vec![ModificationSpecificity::PeptideN(None)]
            } else {
                m.residues
                    .iter()
                    .map(|&r| ModificationSpecificity::PeptideN(Some(r as u8)))
                    .collect()
            }
        }
        ModPosition::AnyCTerm => {
            if m.residues.is_empty() {
                vec![ModificationSpecificity::PeptideC(None)]
            } else {
                m.residues
                    .iter()
                    .map(|&r| ModificationSpecificity::PeptideC(Some(r as u8)))
                    .collect()
            }
        }
        ModPosition::ProteinNTerm => {
            if m.residues.is_empty() {
                vec![ModificationSpecificity::ProteinN(None)]
            } else {
                m.residues
                    .iter()
                    .map(|&r| ModificationSpecificity::ProteinN(Some(r as u8)))
                    .collect()
            }
        }
        ModPosition::ProteinCTerm => {
            if m.residues.is_empty() {
                vec![ModificationSpecificity::ProteinC(None)]
            } else {
                m.residues
                    .iter()
                    .map(|&r| ModificationSpecificity::ProteinC(Some(r as u8)))
                    .collect()
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::spectrum::{IsolationWindow, PrecursorInfo};

    fn sample_ms2_spectrum() -> Spectrum {
        Spectrum {
            scan_number: 42,
            ms_level: MsLevel::MS2,
            retention_time_min: 120.5, // seconds (confusing name but actually stores seconds)
            precursors: vec![PrecursorInfo {
                mz: 500.25,
                charge: Some(2),
                intensity: Some(1e6),
                isolation_window: Some(IsolationWindow {
                    target_mz: 500.25,
                    lower_offset: 1.0,
                    upper_offset: 1.0,
                }),
                source_scan: Some(40),
            }],
            mz_array: vec![100.0, 200.0, 300.5, 400.75],
            intensity_array: vec![1000.0, 2000.0, 500.0, 3000.0],
        }
    }

    #[test]
    fn spectrum_to_raw_basic_fields() {
        let spec = sample_ms2_spectrum();
        let raw = spectrum_to_raw(&spec, 3);

        assert_eq!(raw.file_id, 3);
        assert_eq!(raw.ms_level, 2);
        assert_eq!(raw.id, "42");
        // retention_time_min is in seconds, sage wants minutes
        let expected_rt = (120.5 / 60.0) as f32;
        assert!((raw.scan_start_time - expected_rt).abs() < 1e-5);
        assert_eq!(raw.representation, Representation::Centroid);
    }

    #[test]
    fn spectrum_to_raw_mz_and_intensity() {
        let spec = sample_ms2_spectrum();
        let raw = spectrum_to_raw(&spec, 0);

        assert_eq!(raw.mz.len(), 4);
        assert_eq!(raw.intensity.len(), 4);
        assert!((raw.mz[0] - 100.0f32).abs() < 1e-4);
        assert!((raw.mz[2] - 300.5f32).abs() < 1e-4);
        assert!((raw.intensity[1] - 2000.0f32).abs() < 1e-2);
    }

    #[test]
    fn spectrum_to_raw_precursor() {
        let spec = sample_ms2_spectrum();
        let raw = spectrum_to_raw(&spec, 0);

        assert_eq!(raw.precursors.len(), 1);
        let p = &raw.precursors[0];
        assert!((p.mz - 500.25f32).abs() < 1e-3);
        assert_eq!(p.charge, Some(2));
        assert!(p.intensity.is_some());
        assert!(p.isolation_window.is_some());
    }

    #[test]
    fn spectrum_to_raw_negative_charge_clamped() {
        let mut spec = sample_ms2_spectrum();
        spec.precursors[0].charge = Some(-2);
        let raw = spectrum_to_raw(&spec, 0);
        assert_eq!(raw.precursors[0].charge, Some(2));
    }

    #[test]
    fn spectrum_to_raw_no_precursor() {
        let spec = Spectrum {
            scan_number: 1,
            ms_level: MsLevel::MS1,
            retention_time_min: 60.0,
            precursors: vec![],
            mz_array: vec![100.0],
            intensity_array: vec![1000.0],
        };
        let raw = spectrum_to_raw(&spec, 0);
        assert_eq!(raw.ms_level, 1);
        assert!(raw.precursors.is_empty());
    }

    #[test]
    fn mass_tolerance_ppm() {
        let tol = MassTolerance { value: 20.0, unit: ToleranceUnit::Ppm };
        let sage_tol = mass_tolerance_to_sage(&tol);
        assert_eq!(sage_tol, SageTolerance::Ppm(-20.0, 20.0));
    }

    #[test]
    fn mass_tolerance_da() {
        let tol = MassTolerance { value: 0.5, unit: ToleranceUnit::Da };
        let sage_tol = mass_tolerance_to_sage(&tol);
        assert_eq!(sage_tol, SageTolerance::Da(-0.5, 0.5));
    }

    #[test]
    fn fixed_mod_single_residue() {
        let m = Modification {
            name: "Carbamidomethyl".into(),
            mass_delta: 57.021464,
            residues: vec!['C'],
            position: ModPosition::Anywhere,
        };
        let result = fixed_mod_to_sage(&m);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0, ModificationSpecificity::Residue(b'C'));
        assert!((result[0].1 - 57.021464f32).abs() < 1e-4);
    }

    #[test]
    fn fixed_mod_multiple_residues() {
        let m = Modification {
            name: "Phospho".into(),
            mass_delta: 79.966331,
            residues: vec!['S', 'T', 'Y'],
            position: ModPosition::Anywhere,
        };
        let result = fixed_mod_to_sage(&m);
        assert_eq!(result.len(), 3);
        assert_eq!(result[0].0, ModificationSpecificity::Residue(b'S'));
        assert_eq!(result[1].0, ModificationSpecificity::Residue(b'T'));
        assert_eq!(result[2].0, ModificationSpecificity::Residue(b'Y'));
    }

    #[test]
    fn fixed_mod_nterm() {
        let m = Modification {
            name: "Acetyl".into(),
            mass_delta: 42.010565,
            residues: vec![],
            position: ModPosition::AnyNTerm,
        };
        let result = fixed_mod_to_sage(&m);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0, ModificationSpecificity::PeptideN(None));
    }

    #[test]
    fn fixed_mod_protein_nterm_with_residue() {
        let m = Modification {
            name: "Acetyl".into(),
            mass_delta: 42.010565,
            residues: vec!['M'],
            position: ModPosition::ProteinNTerm,
        };
        let result = fixed_mod_to_sage(&m);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0, ModificationSpecificity::ProteinN(Some(b'M')));
    }

    #[test]
    fn variable_mod_single() {
        let m = Modification {
            name: "Oxidation".into(),
            mass_delta: 15.994915,
            residues: vec!['M'],
            position: ModPosition::Anywhere,
        };
        let result = variable_mod_to_sage(&m);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0, ModificationSpecificity::Residue(b'M'));
        assert_eq!(result[0].1.len(), 1);
        assert!((result[0].1[0] - 15.994915f32).abs() < 1e-4);
    }
}
