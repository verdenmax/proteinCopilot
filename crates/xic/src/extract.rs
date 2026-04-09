//! XIC extraction core logic.
//!
//! Implements the 1.5-pass extraction algorithm:
//! - Pass 0: read target scan to get RT and isolation window
//! - Pass 1: stream all spectra, extracting intensities for target ions

use std::path::Path;

use protein_copilot_core::search_params::{MassTolerance, Modification, ToleranceUnit};
use protein_copilot_core::spectrum::{IsolationWindow, MsLevel};
use protein_copilot_search_engine::matching::{
    generate_b_ions_with_charge, generate_y_ions_with_charge, within_tolerance,
};

use crate::{ExtractionParams, IntensityRule, XicData, XicDataPoint, XicError, XicTrace};
use crate::IonType;

/// A target ion for XIC extraction.
#[derive(Debug, Clone)]
pub struct TargetIon {
    pub label: String,
    pub ion_type: IonType,
    pub ion_number: u32,
    pub charge: u32,
    pub mz: f64,
}

/// Extract intensity for a target m/z from a spectrum's peak list.
///
/// Uses binary search for efficiency (mzML peaks are sorted by m/z).
/// Returns 0.0 if no peak is found within tolerance.
pub fn extract_intensity(
    target_mz: f64,
    exp_mz: &[f64],
    exp_int: &[f64],
    tolerance: &MassTolerance,
    rule: IntensityRule,
) -> f64 {
    if exp_mz.is_empty() || exp_mz.len() != exp_int.len() {
        return 0.0;
    }

    let pos = exp_mz.partition_point(|&m| m < target_mz);

    let max_da = match tolerance.unit {
        ToleranceUnit::Ppm => target_mz * tolerance.value * 1e-6,
        ToleranceUnit::Da => tolerance.value,
    };
    let scan_start = match exp_mz[..pos].iter().rposition(|&m| target_mz - m > max_da * 1.5) {
        Some(i) => i + 1,
        None => 0,
    };
    let scan_end = match exp_mz[pos..].iter().position(|&m| m - target_mz > max_da * 1.5) {
        Some(i) => pos + i,
        None => exp_mz.len(),
    };

    let mut best_intensity = 0.0;
    let mut sum_intensity = 0.0;
    let mut nearest_dist = f64::MAX;
    let mut nearest_intensity = 0.0;
    let mut found = false;

    for i in scan_start..scan_end {
        if within_tolerance(exp_mz[i], target_mz, tolerance) {
            found = true;
            let intensity = exp_int[i];
            if intensity > best_intensity {
                best_intensity = intensity;
            }
            sum_intensity += intensity;
            let dist = (exp_mz[i] - target_mz).abs();
            if dist < nearest_dist {
                nearest_dist = dist;
                nearest_intensity = intensity;
            }
        }
    }

    if !found {
        return 0.0;
    }

    match rule {
        IntensityRule::MaxInWindow => best_intensity,
        IntensityRule::SumInWindow => sum_intensity,
        IntensityRule::NearestPeak => nearest_intensity,
    }
}

/// Check if two isolation windows cover the same DIA region.
///
/// Compares full window bounds (not just center m/z) to avoid
/// conflating adjacent or overlapping windows.
pub fn same_isolation_window(a: &IsolationWindow, b: &IsolationWindow) -> bool {
    let a_lo = a.target_mz - a.lower_offset;
    let a_hi = a.target_mz + a.upper_offset;
    let b_lo = b.target_mz - b.lower_offset;
    let b_hi = b.target_mz + b.upper_offset;

    let center_a = (a_lo + a_hi) / 2.0;
    let center_b = (b_lo + b_hi) / 2.0;
    let center_close = (center_a - center_b).abs() < 1.0;

    let width_a = a_hi - a_lo;
    let width_b = b_hi - b_lo;
    let width_close = width_a > 0.0 && ((width_a - width_b).abs() / width_a) < 0.2;

    center_close && width_close
}

/// Build the list of target fragment ions from a peptide sequence.
///
/// Reuses `matching.rs` fragment ion generators. Does NOT duplicate ion
/// calculation logic — only wraps the results with metadata.
pub fn build_target_ions(
    sequence: &str,
    modifications: &[Modification],
    precursor_charge: i32,
) -> Vec<TargetIon> {
    let max_frag_charge = if precursor_charge >= 3 { 2 } else { 1 };
    let _n = sequence.chars().count().saturating_sub(1);

    let b_mz = generate_b_ions_with_charge(sequence, modifications, max_frag_charge);
    let y_mz = generate_y_ions_with_charge(sequence, modifications, max_frag_charge);

    let max_z = max_frag_charge.max(1) as usize;
    let mut ions = Vec::with_capacity(b_mz.len() + y_mz.len());

    for (idx, &mz) in b_mz.iter().enumerate() {
        let frag_idx = idx / max_z;
        let z = (idx % max_z) + 1;
        let ion_number = (frag_idx + 1) as u32;
        let superscript = if z == 1 { "¹⁺" } else { "²⁺" };
        ions.push(TargetIon {
            label: format!("b{ion_number}{superscript}"),
            ion_type: IonType::B,
            ion_number,
            charge: z as u32,
            mz,
        });
    }

    for (idx, &mz) in y_mz.iter().enumerate() {
        let frag_idx = idx / max_z;
        let z = (idx % max_z) + 1;
        let ion_number = (frag_idx + 1) as u32;
        let superscript = if z == 1 { "¹⁺" } else { "²⁺" };
        ions.push(TargetIon {
            label: format!("y{ion_number}{superscript}"),
            ion_type: IonType::Y,
            ion_number,
            charge: z as u32,
            mz,
        });
    }

    ions
}

/// Extract XIC data for a peptide from an mzML file.
///
/// Uses the 1.5-pass strategy:
/// - Pass 0: `read_spectrum(target_scan)` → get RT and isolation window
/// - Pass 1: `for_each_spectrum()` → stream all spectra, extract intensities
pub fn extract_xic(
    file_path: &Path,
    target_scan: u32,
    peptide_sequence: &str,
    charge: i32,
    precursor_mz: f64,
    modifications: &[Modification],
    params: &ExtractionParams,
) -> Result<XicData, XicError> {
    if peptide_sequence.is_empty() {
        return Err(XicError::InvalidPeptide {
            detail: "peptide sequence is empty".to_string(),
        });
    }

    let info = protein_copilot_spectrum_io::detect_format(file_path)?;
    if info.format != protein_copilot_core::spectrum::SpectrumFormat::MzML {
        return Err(XicError::UnsupportedFormat {
            path: file_path.to_path_buf(),
        });
    }

    let reader = protein_copilot_spectrum_io::create_reader(&info);

    // --- Pass 0: Get target scan info ---
    let target_spectrum = reader.read_spectrum(file_path, target_scan)?;
    let target_rt = target_spectrum.retention_time_sec;
    let target_window = target_spectrum
        .precursors
        .first()
        .and_then(|p| p.isolation_window.as_ref())
        .cloned();

    // --- Build target ion list ---
    let light_ions = build_target_ions(peptide_sequence, modifications, charge);
    let heavy_ions = match &params.label_type {
        Some(label) => crate::heavy::compute_heavy_target_ions(&light_ions, peptide_sequence, label),
        None => Vec::new(),
    };

    let heavy_precursor_mz = params.label_type.as_ref().map(|label| {
        crate::heavy::compute_heavy_precursor_mz(precursor_mz, charge, peptide_sequence, label)
    });

    // --- Pass 1: Stream spectra and extract intensities ---
    let mut ms2_points: Vec<(u32, f64, Vec<f64>, Vec<f64>)> = Vec::new();
    let mut ms1_light_points: Vec<XicDataPoint> = Vec::new();
    let mut ms1_heavy_points: Vec<XicDataPoint> = Vec::new();

    reader.for_each_spectrum(file_path, &mut |spec| {
        let rt = spec.retention_time_sec;

        match spec.ms_level {
            MsLevel::MS1 => {
                let light_int = extract_intensity(
                    precursor_mz,
                    &spec.mz_array,
                    &spec.intensity_array,
                    &params.mz_tolerance,
                    params.intensity_rule,
                );
                ms1_light_points.push(XicDataPoint {
                    retention_time_sec: rt,
                    scan_number: spec.scan_number,
                    intensity: light_int,
                });

                if let Some(heavy_mz) = heavy_precursor_mz {
                    let heavy_int = extract_intensity(
                        heavy_mz,
                        &spec.mz_array,
                        &spec.intensity_array,
                        &params.mz_tolerance,
                        params.intensity_rule,
                    );
                    ms1_heavy_points.push(XicDataPoint {
                        retention_time_sec: rt,
                        scan_number: spec.scan_number,
                        intensity: heavy_int,
                    });
                }
            }
            MsLevel::MS2 => {
                let matches_window = match (&target_window, spec.precursors.first()) {
                    (Some(tw), Some(prec)) => prec
                        .isolation_window
                        .as_ref()
                        .is_some_and(|w| same_isolation_window(tw, w)),
                    (None, _) => true, // DDA: no window filtering
                    _ => false,
                };

                if matches_window {
                    let light_intensities: Vec<f64> = light_ions
                        .iter()
                        .map(|ion| {
                            extract_intensity(
                                ion.mz,
                                &spec.mz_array,
                                &spec.intensity_array,
                                &params.mz_tolerance,
                                params.intensity_rule,
                            )
                        })
                        .collect();

                    let heavy_intensities: Vec<f64> = heavy_ions
                        .iter()
                        .map(|ion| {
                            extract_intensity(
                                ion.mz,
                                &spec.mz_array,
                                &spec.intensity_array,
                                &params.mz_tolerance,
                                params.intensity_rule,
                            )
                        })
                        .collect();

                    ms2_points.push((spec.scan_number, rt, light_intensities, heavy_intensities));
                }
            }
            _ => {}
        }
        Ok(true)
    })?;

    // --- Post-processing ---
    ms2_points.sort_by_key(|(scan, _, _, _)| *scan);

    let target_pos = ms2_points
        .iter()
        .position(|(scan, _, _, _)| *scan == target_scan);
    let (start, end) = match target_pos {
        Some(pos) => {
            let n = params.n_cycles as usize;
            let start = pos.saturating_sub(n);
            let end = (pos + n + 1).min(ms2_points.len());
            (start, end)
        }
        None => (0, ms2_points.len()),
    };
    let windowed = &ms2_points[start..end];

    // Build fragment XIC traces
    let mut fragment_traces: Vec<XicTrace> = light_ions
        .iter()
        .enumerate()
        .map(|(i, ion)| XicTrace {
            ion_label: ion.label.clone(),
            ion_type: ion.ion_type,
            ion_number: ion.ion_number,
            charge: ion.charge,
            theoretical_mz: ion.mz,
            data_points: windowed
                .iter()
                .map(|(scan, rt, ints, _)| XicDataPoint {
                    retention_time_sec: *rt,
                    scan_number: *scan,
                    intensity: ints.get(i).copied().unwrap_or(0.0),
                })
                .collect(),
            is_heavy: false,
        })
        .collect();

    // Top-N selection by total intensity
    fragment_traces.sort_by(|a, b| {
        let a_total: f64 = a.data_points.iter().map(|p| p.intensity).sum();
        let b_total: f64 = b.data_points.iter().map(|p| p.intensity).sum();
        b_total
            .partial_cmp(&a_total)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let top_n = params.top_n_ions.min(fragment_traces.len());
    fragment_traces.truncate(top_n);

    // Build heavy traces (matching top-N selection)
    let heavy_traces: Vec<XicTrace> = if heavy_ions.is_empty() {
        Vec::new()
    } else {
        let top_labels: Vec<String> = fragment_traces.iter().map(|t| t.ion_label.clone()).collect();
        heavy_ions
            .iter()
            .enumerate()
            .filter(|(_, ion)| top_labels.contains(&ion.label))
            .map(|(i, ion)| XicTrace {
                ion_label: ion.label.clone(),
                ion_type: ion.ion_type,
                ion_number: ion.ion_number,
                charge: ion.charge,
                theoretical_mz: ion.mz,
                data_points: windowed
                    .iter()
                    .map(|(scan, rt, _, heavy_ints)| XicDataPoint {
                        retention_time_sec: *rt,
                        scan_number: *scan,
                        intensity: heavy_ints.get(i).copied().unwrap_or(0.0),
                    })
                    .collect(),
                is_heavy: true,
            })
            .collect()
    };

    // Trim MS1 to same RT range as MS2
    let rt_range = if let (Some(first), Some(last)) = (windowed.first(), windowed.last()) {
        Some((first.1, last.1))
    } else {
        None
    };

    let ms1_precursor_xic = if ms1_light_points.is_empty() {
        None
    } else {
        let filtered: Vec<XicDataPoint> = match rt_range {
            Some((lo, hi)) => ms1_light_points
                .into_iter()
                .filter(|p| p.retention_time_sec >= lo && p.retention_time_sec <= hi)
                .collect(),
            None => ms1_light_points,
        };
        if filtered.is_empty() {
            None
        } else {
            Some(XicTrace {
                ion_label: "precursor".to_string(),
                ion_type: IonType::Precursor,
                ion_number: 0,
                charge: charge as u32,
                theoretical_mz: precursor_mz,
                data_points: filtered,
                is_heavy: false,
            })
        }
    };

    let ms1_heavy_precursor_xic = if ms1_heavy_points.is_empty() {
        None
    } else {
        let filtered: Vec<XicDataPoint> = match rt_range {
            Some((lo, hi)) => ms1_heavy_points
                .into_iter()
                .filter(|p| p.retention_time_sec >= lo && p.retention_time_sec <= hi)
                .collect(),
            None => ms1_heavy_points,
        };
        if filtered.is_empty() {
            None
        } else {
            Some(XicTrace {
                ion_label: "precursor (heavy)".to_string(),
                ion_type: IonType::Precursor,
                ion_number: 0,
                charge: charge as u32,
                theoretical_mz: heavy_precursor_mz.unwrap_or(precursor_mz),
                data_points: filtered,
                is_heavy: true,
            })
        }
    };

    Ok(XicData {
        peptide_sequence: peptide_sequence.to_string(),
        target_rt_sec: target_rt,
        target_scan,
        charge,
        precursor_mz,
        ms1_precursor_xic,
        ms1_heavy_precursor_xic,
        fragment_xic_traces: fragment_traces,
        heavy_fragment_xic_traces: heavy_traces,
        extraction_params: params.clone(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extract_intensity_max_in_window() {
        let mz = vec![100.0, 200.0, 200.01, 200.02, 300.0];
        let int = vec![50.0, 100.0, 200.0, 150.0, 75.0];
        let tol = MassTolerance {
            value: 0.05,
            unit: ToleranceUnit::Da,
        };
        let result = extract_intensity(200.01, &mz, &int, &tol, IntensityRule::MaxInWindow);
        assert!((result - 200.0).abs() < 0.01, "expected 200.0, got {result}");
    }

    #[test]
    fn extract_intensity_sum_in_window() {
        let mz = vec![100.0, 200.0, 200.01, 200.02, 300.0];
        let int = vec![50.0, 100.0, 200.0, 150.0, 75.0];
        let tol = MassTolerance {
            value: 0.05,
            unit: ToleranceUnit::Da,
        };
        let result = extract_intensity(200.01, &mz, &int, &tol, IntensityRule::SumInWindow);
        assert!((result - 450.0).abs() < 0.01, "expected 450.0, got {result}");
    }

    #[test]
    fn extract_intensity_nearest_peak() {
        let mz = vec![100.0, 200.0, 200.01, 200.02, 300.0];
        let int = vec![50.0, 100.0, 200.0, 150.0, 75.0];
        let tol = MassTolerance {
            value: 0.05,
            unit: ToleranceUnit::Da,
        };
        let result = extract_intensity(200.01, &mz, &int, &tol, IntensityRule::NearestPeak);
        assert!((result - 200.0).abs() < 0.01, "expected 200.0, got {result}");
    }

    #[test]
    fn extract_intensity_no_match_returns_zero() {
        let mz = vec![100.0, 200.0, 300.0];
        let int = vec![50.0, 100.0, 75.0];
        let tol = MassTolerance {
            value: 0.01,
            unit: ToleranceUnit::Da,
        };
        let result = extract_intensity(250.0, &mz, &int, &tol, IntensityRule::MaxInWindow);
        assert!((result - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn extract_intensity_empty_spectrum() {
        let tol = MassTolerance {
            value: 0.05,
            unit: ToleranceUnit::Da,
        };
        assert_eq!(
            extract_intensity(200.0, &[], &[], &tol, IntensityRule::MaxInWindow),
            0.0
        );
    }

    #[test]
    fn extract_intensity_ppm_tolerance() {
        let mz = vec![500.0, 500.005, 500.01];
        let int = vec![100.0, 200.0, 300.0];
        let tol = MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        };
        let result = extract_intensity(500.005, &mz, &int, &tol, IntensityRule::MaxInWindow);
        assert!(result > 0.0, "should find a peak within 20 ppm");
    }

    #[test]
    fn same_isolation_window_identical() {
        let w = IsolationWindow {
            target_mz: 500.0,
            lower_offset: 12.5,
            upper_offset: 12.5,
        };
        assert!(same_isolation_window(&w, &w));
    }

    #[test]
    fn same_isolation_window_different_center() {
        let a = IsolationWindow {
            target_mz: 500.0,
            lower_offset: 12.5,
            upper_offset: 12.5,
        };
        let b = IsolationWindow {
            target_mz: 525.0,
            lower_offset: 12.5,
            upper_offset: 12.5,
        };
        assert!(!same_isolation_window(&a, &b));
    }

    #[test]
    fn same_isolation_window_different_width() {
        let a = IsolationWindow {
            target_mz: 500.0,
            lower_offset: 12.5,
            upper_offset: 12.5,
        };
        let b = IsolationWindow {
            target_mz: 500.0,
            lower_offset: 5.0,
            upper_offset: 5.0,
        };
        assert!(!same_isolation_window(&a, &b));
    }

    #[test]
    fn build_target_ions_simple_peptide() {
        let ions = build_target_ions("PEPTIDE", &[], 2);
        assert_eq!(ions.len(), 12); // 6 b + 6 y at charge 1
        assert!(ions[0].label.starts_with('b'));
        assert!(ions.last().unwrap().label.starts_with('y'));
    }

    #[test]
    fn build_target_ions_high_charge_gets_doubly_charged() {
        let ions = build_target_ions("PEPTIDE", &[], 3);
        assert_eq!(ions.len(), 24); // 6 b × 2 charges + 6 y × 2 charges
        let has_double = ions.iter().any(|i| i.charge == 2);
        assert!(has_double, "charge 3 precursor should produce doubly-charged fragments");
    }

    #[test]
    fn build_target_ions_empty_sequence() {
        let ions = build_target_ions("", &[], 2);
        assert!(ions.is_empty());
    }
}
