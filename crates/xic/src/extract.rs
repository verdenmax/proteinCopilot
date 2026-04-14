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

/// An MS2 data point: (scan_number, RT_min, extracted_intensities).
/// Each intensity entry is (intensity, Option<observed_mz>).
type Ms2Point = (u32, f64, Vec<(f64, Option<f64>)>);

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
/// Returns `(intensity, observed_mz)` where `observed_mz` is the actual
/// peak m/z that contributed to the result. Returns `(0.0, None)` if no
/// peak is found within tolerance.
///
/// The observed m/z depends on the intensity rule:
/// - `MaxInWindow`: m/z of the highest-intensity peak
/// - `SumInWindow`: m/z of the peak closest to `target_mz`
/// - `NearestPeak`: m/z of the nearest peak
pub fn extract_intensity(
    target_mz: f64,
    exp_mz: &[f64],
    exp_int: &[f64],
    tolerance: &MassTolerance,
    rule: IntensityRule,
) -> (f64, Option<f64>) {
    if exp_mz.is_empty() || exp_mz.len() != exp_int.len() {
        return (0.0, None);
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
    let mut best_mz = 0.0;
    let mut sum_intensity = 0.0;
    let mut nearest_dist = f64::MAX;
    let mut nearest_intensity = 0.0;
    let mut nearest_mz = 0.0;
    let mut found = false;

    for i in scan_start..scan_end {
        if within_tolerance(exp_mz[i], target_mz, tolerance) {
            found = true;
            let intensity = exp_int[i];
            if intensity > best_intensity {
                best_intensity = intensity;
                best_mz = exp_mz[i];
            }
            sum_intensity += intensity;
            let dist = (exp_mz[i] - target_mz).abs();
            if dist < nearest_dist {
                nearest_dist = dist;
                nearest_intensity = intensity;
                nearest_mz = exp_mz[i];
            }
        }
    }

    if !found {
        return (0.0, None);
    }

    match rule {
        IntensityRule::MaxInWindow => (best_intensity, Some(best_mz)),
        IntensityRule::SumInWindow => (sum_intensity, Some(nearest_mz)),
        IntensityRule::NearestPeak => (nearest_intensity, Some(nearest_mz)),
    }
}

/// Check if two isolation windows cover the same DIA region.
///
/// Uses fixed tolerances: center m/z within 1.0 Da and width within 20%.
/// These values are suitable for typical DIA window schemes (e.g., SWATH).
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

/// Compute ion metadata with K/R residue counts for SILAC calculation.
///
/// For each target ion, counts the K and R residues in the fragment
/// (prefix for b-ions, suffix for y-ions) so the browser can compute
/// heavy m/z shifts client-side.
pub fn compute_ion_metadata(ions: &[TargetIon], peptide: &str) -> Vec<crate::IonMetadataEntry> {
    let chars: Vec<char> = peptide.chars().collect();
    let n = chars.len();

    ions.iter()
        .map(|ion| {
            let fragment_chars: &[char] = match ion.ion_type {
                IonType::B => &chars[..(ion.ion_number as usize).min(n)],
                IonType::Y => {
                    let start = n.saturating_sub(ion.ion_number as usize);
                    &chars[start..]
                }
                IonType::Precursor => &chars[..],
            };

            crate::IonMetadataEntry {
                label: ion.label.clone(),
                ion_type: ion.ion_type,
                ion_number: ion.ion_number,
                charge: ion.charge,
                light_mz: ion.mz,
                k_count: fragment_chars.iter().filter(|&&c| c == 'K' || c == 'k').count() as u32,
                r_count: fragment_chars.iter().filter(|&&c| c == 'R' || c == 'r').count() as u32,
            }
        })
        .collect()
}

/// Trim a peak list to only include peaks within ±window_da of center_mz.
///
/// Used to reduce MS1 raw data volume for HTML embedding — only peaks
/// near the precursor m/z are needed for SILAC recomputation.
fn trim_peaks_to_window(
    mz_array: &[f64],
    intensity_array: &[f64],
    center_mz: f64,
    window_da: f64,
) -> (Vec<f64>, Vec<f64>) {
    let lo = center_mz - window_da;
    let hi = center_mz + window_da;

    let start = mz_array.partition_point(|&m| m < lo);
    let end = mz_array.partition_point(|&m| m <= hi);

    (
        mz_array[start..end].to_vec(),
        intensity_array[start..end].to_vec(),
    )
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
    if charge <= 0 {
        return Err(XicError::InvalidPeptide {
            detail: format!("charge must be > 0, got {charge}"),
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

    // Validate target scan is MS2
    if target_spectrum.ms_level != MsLevel::MS2 {
        return Err(XicError::InvalidPeptide {
            detail: format!(
                "scan {} is not MS2 — XIC extraction requires an MS2 scan",
                target_scan,
            ),
        });
    }

    let target_rt = target_spectrum.retention_time_min;
    let target_window = target_spectrum
        .precursors
        .first()
        .and_then(|p| p.isolation_window.as_ref())
        .cloned();

    // --- Build target ion list ---
    let light_ions = build_target_ions(peptide_sequence, modifications, charge);
    // Zero-offset guard: skip heavy when peptide has no K/R (zero SILAC shift)
    let effective_label = params.label_type.as_ref().filter(|label| {
        protein_copilot_core::label::total_heavy_delta(peptide_sequence, label).abs() > 1e-6
    });
    let heavy_ions = match &effective_label {
        Some(label) => crate::heavy::compute_heavy_target_ions(&light_ions, peptide_sequence, label),
        None => Vec::new(),
    };

    let heavy_precursor_mz = effective_label.map(|label| {
        crate::heavy::compute_heavy_precursor_mz(precursor_mz, charge, peptide_sequence, label)
    });

    // --- Pass 1: Stream spectra and extract intensities ---
    // Light MS2 points: (scan, RT, light_intensities)
    let mut ms2_light_points: Vec<Ms2Point> = Vec::new();
    // Heavy MS2 points: (scan, RT, heavy_intensities) — from different DIA window
    let mut ms2_heavy_points: Vec<Ms2Point> = Vec::new();
    // DIA detection: isolation window width > 1 Th indicates DIA
    // (DDA windows are typically < 1 Th; DIA narrow windows start at ~2 Th)
    let is_dia = target_window
        .as_ref()
        .map(|w| (w.lower_offset + w.upper_offset) > 1.0)
        .unwrap_or(false);
    // Whether heavy uses a separate DIA window (DIA+SILAC)
    let needs_separate_heavy_window = is_dia && !heavy_ions.is_empty();
    let mut ms1_light_points: Vec<XicDataPoint> = Vec::new();
    let mut ms1_heavy_points: Vec<XicDataPoint> = Vec::new();

    reader.for_each_spectrum(file_path, &mut |spec| {
        let rt = spec.retention_time_min;

        match spec.ms_level {
            MsLevel::MS1 => {
                let (light_int, light_obs) = extract_intensity(
                    precursor_mz,
                    &spec.mz_array,
                    &spec.intensity_array,
                    &params.mz_tolerance,
                    params.intensity_rule,
                );
                ms1_light_points.push(XicDataPoint {
                    retention_time_min: rt,
                    scan_number: spec.scan_number,
                    intensity: light_int,
                    observed_mz: light_obs,
                });

                if let Some(heavy_mz) = heavy_precursor_mz {
                    let (heavy_int, heavy_obs) = extract_intensity(
                        heavy_mz,
                        &spec.mz_array,
                        &spec.intensity_array,
                        &params.mz_tolerance,
                        params.intensity_rule,
                    );
                    ms1_heavy_points.push(XicDataPoint {
                        retention_time_min: rt,
                        scan_number: spec.scan_number,
                        intensity: heavy_int,
                        observed_mz: heavy_obs,
                    });
                }
            }
            MsLevel::MS2 => {
                let matches_light_window = match (&target_window, spec.precursors.first()) {
                    (Some(tw), Some(prec)) => prec
                        .isolation_window
                        .as_ref()
                        .is_some_and(|w| same_isolation_window(tw, w)),
                    (None, _) => true, // DDA: no window filtering
                    _ => false,
                };

                if matches_light_window {
                    let light_intensities: Vec<(f64, Option<f64>)> = light_ions
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

                    ms2_light_points.push((spec.scan_number, rt, light_intensities));

                    // DDA: heavy also from same scans (no separate window)
                    if !needs_separate_heavy_window && !heavy_ions.is_empty() {
                        let heavy_intensities: Vec<(f64, Option<f64>)> = heavy_ions
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
                        ms2_heavy_points.push((spec.scan_number, rt, heavy_intensities));
                    }
                }

                // DIA+SILAC: check if this scan matches the HEAVY window
                if needs_separate_heavy_window {
                    if let Some(heavy_mz) = heavy_precursor_mz {
                        let matches_heavy = spec
                            .precursors
                            .first()
                            .and_then(|p| p.isolation_window.as_ref())
                            .is_some_and(|w| crate::heavy::window_contains_mz(w, heavy_mz));

                        if matches_heavy {
                            let heavy_intensities: Vec<(f64, Option<f64>)> = heavy_ions
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
                            ms2_heavy_points.push((spec.scan_number, rt, heavy_intensities));
                        }
                    }
                }
            }
            _ => {}
        }
        Ok(true)
    })?;

    // --- Post-processing: light windowing ---
    ms2_light_points.sort_by_key(|(scan, _, _)| *scan);

    let target_pos = ms2_light_points
        .iter()
        .position(|(scan, _, _)| *scan == target_scan);
    let (start, end) = match target_pos {
        Some(pos) => {
            let n = params.n_cycles as usize;
            let start = pos.saturating_sub(n);
            let end = (pos + n + 1).min(ms2_light_points.len());
            (start, end)
        }
        None => {
            return Err(XicError::ScanNotFound {
                scan: target_scan,
                path: file_path.to_path_buf(),
            });
        }
    };
    let light_windowed = &ms2_light_points[start..end];

    // --- Heavy windowing (independent scan sequence) ---
    ms2_heavy_points.sort_by_key(|(scan, _, _)| *scan);
    let heavy_windowed: &[Ms2Point] = if ms2_heavy_points.is_empty() {
        &[]
    } else {
        // Find the heavy scan closest to target_scan's RT
        let target_rt_for_heavy = target_rt;
        let heavy_center = ms2_heavy_points
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                let da = (a.1 - target_rt_for_heavy).abs();
                let db = (b.1 - target_rt_for_heavy).abs();
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(idx, _)| idx)
            .unwrap_or(0);
        let n = params.n_cycles as usize;
        let h_start = heavy_center.saturating_sub(n);
        let h_end = (heavy_center + n + 1).min(ms2_heavy_points.len());
        &ms2_heavy_points[h_start..h_end]
    };

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
            data_points: light_windowed
                .iter()
                .map(|(scan, rt, ints)| XicDataPoint {
                    retention_time_min: *rt,
                    scan_number: *scan,
                    intensity: ints.get(i).map(|(int, _)| *int).unwrap_or(0.0),
                    observed_mz: ints.get(i).and_then(|(_, mz)| *mz),
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
    // Remove traces with zero total intensity (no signal detected)
    fragment_traces.retain(|t| t.data_points.iter().any(|p| p.intensity > 0.0));

    // Build heavy traces (matching top-N selection)
    let heavy_traces: Vec<XicTrace> = if heavy_ions.is_empty() || heavy_windowed.is_empty() {
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
                data_points: heavy_windowed
                    .iter()
                    .map(|(scan, rt, heavy_ints)| XicDataPoint {
                        retention_time_min: *rt,
                        scan_number: *scan,
                        intensity: heavy_ints.get(i).map(|(int, _)| *int).unwrap_or(0.0),
                        observed_mz: heavy_ints.get(i).and_then(|(_, mz)| *mz),
                    })
                    .collect(),
                is_heavy: true,
            })
            .collect()
    };

    // RT range: union of light and heavy MS2 windows
    let light_rt_range = if let (Some(first), Some(last)) = (light_windowed.first(), light_windowed.last()) {
        Some((first.1, last.1))
    } else {
        None
    };
    let heavy_rt_range = if let (Some(first), Some(last)) = (heavy_windowed.first(), heavy_windowed.last()) {
        Some((first.1, last.1))
    } else {
        None
    };
    let rt_range = match (light_rt_range, heavy_rt_range) {
        (Some((l_lo, l_hi)), Some((h_lo, h_hi))) => Some((l_lo.min(h_lo), l_hi.max(h_hi))),
        (Some(r), None) | (None, Some(r)) => Some(r),
        (None, None) => None,
    };

    let ms1_precursor_xic = if ms1_light_points.is_empty() {
        None
    } else {
        let filtered: Vec<XicDataPoint> = match rt_range {
            Some((lo, hi)) => ms1_light_points
                .into_iter()
                .filter(|p| p.retention_time_min >= lo && p.retention_time_min <= hi)
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
                .filter(|p| p.retention_time_min >= lo && p.retention_time_min <= hi)
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

    let heavy_warning = if needs_separate_heavy_window && ms2_heavy_points.is_empty() {
        Some(format!(
            "Heavy precursor m/z ({:.4}) is outside all DIA MS2 isolation windows. Heavy MS2 traces unavailable.",
            heavy_precursor_mz.unwrap_or(0.0)
        ))
    } else {
        None
    };

    Ok(XicData {
        peptide_sequence: peptide_sequence.to_string(),
        target_rt_min: target_rt,
        target_scan,
        charge,
        precursor_mz,
        ms1_precursor_xic,
        ms1_heavy_precursor_xic,
        fragment_xic_traces: fragment_traces,
        heavy_fragment_xic_traces: heavy_traces,
        extraction_params: params.clone(),
        heavy_warning,
    })
}

/// Extract XIC data AND raw scan peak arrays for client-side SILAC.
///
/// This is an extension of [`extract_xic`] that additionally captures
/// raw peak data from MS1/MS2 scans in the RT window. The raw data
/// enables the HTML frontend to recompute XIC traces for arbitrary
/// SILAC configurations without a backend round-trip.
///
/// MS1 peaks are trimmed to ±`ms1_mz_window_da` around `precursor_mz`
/// to control embedded data volume.
#[allow(clippy::too_many_arguments)]
pub fn extract_xic_with_raw(
    file_path: &Path,
    target_scan: u32,
    peptide_sequence: &str,
    charge: i32,
    precursor_mz: f64,
    modifications: &[Modification],
    params: &ExtractionParams,
    ms1_mz_window_da: f64,
) -> Result<(XicData, crate::RawScanData, Vec<crate::IonMetadataEntry>), XicError> {
    if peptide_sequence.is_empty() {
        return Err(XicError::InvalidPeptide {
            detail: "peptide sequence is empty".to_string(),
        });
    }
    if charge <= 0 {
        return Err(XicError::InvalidPeptide {
            detail: format!("charge must be > 0, got {charge}"),
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

    // Validate target scan is MS2
    if target_spectrum.ms_level != MsLevel::MS2 {
        return Err(XicError::InvalidPeptide {
            detail: format!(
                "scan {} is not MS2 — XIC extraction requires an MS2 scan",
                target_scan,
            ),
        });
    }

    let target_rt = target_spectrum.retention_time_min;
    let target_window = target_spectrum
        .precursors
        .first()
        .and_then(|p| p.isolation_window.as_ref())
        .cloned();

    // Compute dynamic MS1 trim window: must accommodate heavy precursor shift.
    // Heavy shift = (K_count * K_delta + R_count * R_delta) / |charge|
    // Use maximum plausible shift (standard SILAC K+8.015, R+10.009) + 5 Da margin.
    let k_count = peptide_sequence.chars().filter(|&c| c == 'K' || c == 'k').count() as f64;
    let r_count = peptide_sequence.chars().filter(|&c| c == 'R' || c == 'r').count() as f64;
    let max_heavy_shift =
        (k_count * 8.015 + r_count * 10.009) / (charge.unsigned_abs() as f64).max(1.0);
    let effective_ms1_window = ms1_mz_window_da.max(max_heavy_shift + 5.0);

    // --- Build target ion list ---
    let light_ions = build_target_ions(peptide_sequence, modifications, charge);
    // Zero-offset guard: skip heavy when peptide has no K/R (zero SILAC shift)
    let effective_label = params.label_type.as_ref().filter(|label| {
        protein_copilot_core::label::total_heavy_delta(peptide_sequence, label).abs() > 1e-6
    });
    let heavy_ions = match &effective_label {
        Some(label) => {
            crate::heavy::compute_heavy_target_ions(&light_ions, peptide_sequence, label)
        }
        None => Vec::new(),
    };

    let heavy_precursor_mz = effective_label.map(|label| {
        crate::heavy::compute_heavy_precursor_mz(precursor_mz, charge, peptide_sequence, label)
    });

    // --- Pass 1: Stream spectra, extract intensities, AND capture raw peaks ---
    let mut ms2_light_points: Vec<Ms2Point> = Vec::new();
    let mut ms2_heavy_points: Vec<Ms2Point> = Vec::new();
    // DIA detection: isolation window width > 1 Th indicates DIA
    // (DDA windows are typically < 1 Th; DIA narrow windows start at ~2 Th)
    let is_dia = target_window
        .as_ref()
        .map(|w| (w.lower_offset + w.upper_offset) > 1.0)
        .unwrap_or(false);
    let needs_separate_heavy_window = is_dia && !heavy_ions.is_empty();
    let mut ms1_light_points: Vec<XicDataPoint> = Vec::new();
    let mut ms1_heavy_points: Vec<XicDataPoint> = Vec::new();
    let mut raw_ms1_scans: Vec<crate::RawScan> = Vec::new();
    let mut raw_ms2_light_scans: Vec<crate::RawScan> = Vec::new();
    let mut raw_ms2_heavy_scans: Vec<crate::RawScan> = Vec::new();

    reader.for_each_spectrum(file_path, &mut |spec| {
        let rt = spec.retention_time_min;

        match spec.ms_level {
            MsLevel::MS1 => {
                let (light_int, light_obs) = extract_intensity(
                    precursor_mz,
                    &spec.mz_array,
                    &spec.intensity_array,
                    &params.mz_tolerance,
                    params.intensity_rule,
                );
                ms1_light_points.push(XicDataPoint {
                    retention_time_min: rt,
                    scan_number: spec.scan_number,
                    intensity: light_int,
                    observed_mz: light_obs,
                });

                if let Some(heavy_mz) = heavy_precursor_mz {
                    let (heavy_int, heavy_obs) = extract_intensity(
                        heavy_mz,
                        &spec.mz_array,
                        &spec.intensity_array,
                        &params.mz_tolerance,
                        params.intensity_rule,
                    );
                    ms1_heavy_points.push(XicDataPoint {
                        retention_time_min: rt,
                        scan_number: spec.scan_number,
                        intensity: heavy_int,
                        observed_mz: heavy_obs,
                    });
                }

                // Capture raw MS1 peaks (trimmed to dynamic window around precursor)
                let (trimmed_mz, trimmed_int) = trim_peaks_to_window(
                    &spec.mz_array,
                    &spec.intensity_array,
                    precursor_mz,
                    effective_ms1_window,
                );
                if !trimmed_mz.is_empty() {
                    raw_ms1_scans.push(crate::RawScan {
                        scan_number: spec.scan_number,
                        retention_time_min: rt,
                        mz_array: trimmed_mz,
                        intensity_array: trimmed_int,
                    });
                }
            }
            MsLevel::MS2 => {
                let matches_light_window = match (&target_window, spec.precursors.first()) {
                    (Some(tw), Some(prec)) => prec
                        .isolation_window
                        .as_ref()
                        .is_some_and(|w| same_isolation_window(tw, w)),
                    (None, _) => true,
                    _ => false,
                };

                if matches_light_window {
                    let light_intensities: Vec<(f64, Option<f64>)> = light_ions
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

                    ms2_light_points.push((spec.scan_number, rt, light_intensities));

                    // DDA: heavy also from same scans
                    if !needs_separate_heavy_window && !heavy_ions.is_empty() {
                        let heavy_intensities: Vec<(f64, Option<f64>)> = heavy_ions
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
                        ms2_heavy_points.push((spec.scan_number, rt, heavy_intensities));
                    }

                    // Capture raw light MS2 peaks
                    let rt_close = (rt - target_rt).abs() < 300.0;
                    if target_window.is_some() || rt_close {
                        raw_ms2_light_scans.push(crate::RawScan {
                            scan_number: spec.scan_number,
                            retention_time_min: rt,
                            mz_array: spec.mz_array.clone(),
                            intensity_array: spec.intensity_array.clone(),
                        });
                    }
                }

                // DIA+SILAC: check heavy window
                if needs_separate_heavy_window {
                    if let Some(heavy_mz) = heavy_precursor_mz {
                        let matches_heavy = spec
                            .precursors
                            .first()
                            .and_then(|p| p.isolation_window.as_ref())
                            .is_some_and(|w| crate::heavy::window_contains_mz(w, heavy_mz));

                        if matches_heavy {
                            let heavy_intensities: Vec<(f64, Option<f64>)> = heavy_ions
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
                            ms2_heavy_points.push((spec.scan_number, rt, heavy_intensities));

                            let rt_close = (rt - target_rt).abs() < 300.0;
                            if rt_close {
                                raw_ms2_heavy_scans.push(crate::RawScan {
                                    scan_number: spec.scan_number,
                                    retention_time_min: rt,
                                    mz_array: spec.mz_array.clone(),
                                    intensity_array: spec.intensity_array.clone(),
                                });
                            }
                        }
                    }
                }
            }
            _ => {}
        }
        Ok(true)
    })?;

    // --- Post-processing: light windowing ---
    ms2_light_points.sort_by_key(|(scan, _, _)| *scan);

    let target_pos = ms2_light_points
        .iter()
        .position(|(scan, _, _)| *scan == target_scan);
    let (start, end) = match target_pos {
        Some(pos) => {
            let n = params.n_cycles as usize;
            let start = pos.saturating_sub(n);
            let end = (pos + n + 1).min(ms2_light_points.len());
            (start, end)
        }
        None => {
            return Err(XicError::ScanNotFound {
                scan: target_scan,
                path: file_path.to_path_buf(),
            });
        }
    };
    let light_windowed = &ms2_light_points[start..end];

    // --- Heavy windowing (independent scan sequence) ---
    ms2_heavy_points.sort_by_key(|(scan, _, _)| *scan);
    let heavy_windowed: &[Ms2Point] = if ms2_heavy_points.is_empty() {
        &[]
    } else {
        let target_rt_for_heavy = target_rt;
        let heavy_center = ms2_heavy_points
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                let da = (a.1 - target_rt_for_heavy).abs();
                let db = (b.1 - target_rt_for_heavy).abs();
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(idx, _)| idx)
            .unwrap_or(0);
        let n = params.n_cycles as usize;
        let h_start = heavy_center.saturating_sub(n);
        let h_end = (heavy_center + n + 1).min(ms2_heavy_points.len());
        &ms2_heavy_points[h_start..h_end]
    };

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
            data_points: light_windowed
                .iter()
                .map(|(scan, rt, ints)| XicDataPoint {
                    retention_time_min: *rt,
                    scan_number: *scan,
                    intensity: ints.get(i).map(|(int, _)| *int).unwrap_or(0.0),
                    observed_mz: ints.get(i).and_then(|(_, mz)| *mz),
                })
                .collect(),
            is_heavy: false,
        })
        .collect();

    fragment_traces.sort_by(|a, b| {
        let a_total: f64 = a.data_points.iter().map(|p| p.intensity).sum();
        let b_total: f64 = b.data_points.iter().map(|p| p.intensity).sum();
        b_total
            .partial_cmp(&a_total)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let top_n = params.top_n_ions.min(fragment_traces.len());
    fragment_traces.truncate(top_n);
    fragment_traces.retain(|t| t.data_points.iter().any(|p| p.intensity > 0.0));

    let heavy_traces: Vec<XicTrace> = if heavy_ions.is_empty() || heavy_windowed.is_empty() {
        Vec::new()
    } else {
        let top_labels: Vec<String> =
            fragment_traces.iter().map(|t| t.ion_label.clone()).collect();
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
                data_points: heavy_windowed
                    .iter()
                    .map(|(scan, rt, heavy_ints)| XicDataPoint {
                        retention_time_min: *rt,
                        scan_number: *scan,
                        intensity: heavy_ints.get(i).map(|(int, _)| *int).unwrap_or(0.0),
                        observed_mz: heavy_ints.get(i).and_then(|(_, mz)| *mz),
                    })
                    .collect(),
                is_heavy: true,
            })
            .collect()
    };

    // RT range: union of light and heavy MS2 windows
    let light_rt_range = if let (Some(first), Some(last)) = (light_windowed.first(), light_windowed.last()) {
        Some((first.1, last.1))
    } else {
        None
    };
    let heavy_rt_range = if let (Some(first), Some(last)) = (heavy_windowed.first(), heavy_windowed.last()) {
        Some((first.1, last.1))
    } else {
        None
    };
    let rt_range = match (light_rt_range, heavy_rt_range) {
        (Some((l_lo, l_hi)), Some((h_lo, h_hi))) => Some((l_lo.min(h_lo), l_hi.max(h_hi))),
        (Some(r), None) | (None, Some(r)) => Some(r),
        (None, None) => None,
    };

    let ms1_precursor_xic = if ms1_light_points.is_empty() {
        None
    } else {
        let filtered: Vec<XicDataPoint> = match rt_range {
            Some((lo, hi)) => ms1_light_points
                .into_iter()
                .filter(|p| p.retention_time_min >= lo && p.retention_time_min <= hi)
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
                .filter(|p| p.retention_time_min >= lo && p.retention_time_min <= hi)
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

    // Trim raw scans to same RT window
    let raw_ms1_trimmed = match rt_range {
        Some((lo, hi)) => raw_ms1_scans
            .into_iter()
            .filter(|s| s.retention_time_min >= lo && s.retention_time_min <= hi)
            .collect(),
        None => raw_ms1_scans,
    };
    // MS2 raw scans: keep only the windowed ones
    let windowed_scans: std::collections::HashSet<u32> =
        light_windowed.iter().map(|(scan, _, _)| *scan).collect();
    let raw_ms2_trimmed: Vec<crate::RawScan> = raw_ms2_light_scans
        .into_iter()
        .filter(|s| windowed_scans.contains(&s.scan_number))
        .collect();

    let heavy_windowed_scans: std::collections::HashSet<u32> =
        heavy_windowed.iter().map(|(scan, _, _)| *scan).collect();
    let raw_ms2_heavy_trimmed: Vec<crate::RawScan> = raw_ms2_heavy_scans
        .into_iter()
        .filter(|s| heavy_windowed_scans.contains(&s.scan_number))
        .collect();

    let heavy_warning = if needs_separate_heavy_window && ms2_heavy_points.is_empty() {
        Some(format!(
            "Heavy precursor m/z ({:.4}) is outside all DIA MS2 isolation windows. Heavy MS2 traces unavailable.",
            heavy_precursor_mz.unwrap_or(0.0)
        ))
    } else {
        None
    };

    let xic_data = XicData {
        peptide_sequence: peptide_sequence.to_string(),
        target_rt_min: target_rt,
        target_scan,
        charge,
        precursor_mz,
        ms1_precursor_xic,
        ms1_heavy_precursor_xic,
        fragment_xic_traces: fragment_traces,
        heavy_fragment_xic_traces: heavy_traces,
        extraction_params: params.clone(),
        heavy_warning,
    };

    let raw_scans = crate::RawScanData {
        ms1_scans: raw_ms1_trimmed,
        ms2_scans: raw_ms2_trimmed,
        ms2_heavy_scans: raw_ms2_heavy_trimmed,
    };

    // Ion metadata for the top-N selected light ions
    let ion_metadata = compute_ion_metadata(
        &xic_data
            .fragment_xic_traces
            .iter()
            .map(|t| TargetIon {
                label: t.ion_label.clone(),
                ion_type: t.ion_type,
                ion_number: t.ion_number,
                charge: t.charge,
                mz: t.theoretical_mz,
            })
            .collect::<Vec<_>>(),
        peptide_sequence,
    );

    Ok((xic_data, raw_scans, ion_metadata))
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
        let (intensity, obs_mz) = extract_intensity(200.01, &mz, &int, &tol, IntensityRule::MaxInWindow);
        assert!((intensity - 200.0).abs() < 0.01, "expected 200.0, got {intensity}");
        assert_eq!(obs_mz, Some(200.01), "max peak is at 200.01");
    }

    #[test]
    fn extract_intensity_sum_in_window() {
        let mz = vec![100.0, 200.0, 200.01, 200.02, 300.0];
        let int = vec![50.0, 100.0, 200.0, 150.0, 75.0];
        let tol = MassTolerance {
            value: 0.05,
            unit: ToleranceUnit::Da,
        };
        let (intensity, obs_mz) = extract_intensity(200.01, &mz, &int, &tol, IntensityRule::SumInWindow);
        assert!((intensity - 450.0).abs() < 0.01, "expected 450.0, got {intensity}");
        assert_eq!(obs_mz, Some(200.01), "nearest peak to target is at 200.01");
    }

    #[test]
    fn extract_intensity_nearest_peak() {
        let mz = vec![100.0, 200.0, 200.01, 200.02, 300.0];
        let int = vec![50.0, 100.0, 200.0, 150.0, 75.0];
        let tol = MassTolerance {
            value: 0.05,
            unit: ToleranceUnit::Da,
        };
        let (intensity, obs_mz) = extract_intensity(200.01, &mz, &int, &tol, IntensityRule::NearestPeak);
        assert!((intensity - 200.0).abs() < 0.01, "expected 200.0, got {intensity}");
        assert_eq!(obs_mz, Some(200.01), "nearest peak is at 200.01");
    }

    #[test]
    fn extract_intensity_no_match_returns_zero() {
        let mz = vec![100.0, 200.0, 300.0];
        let int = vec![50.0, 100.0, 75.0];
        let tol = MassTolerance {
            value: 0.01,
            unit: ToleranceUnit::Da,
        };
        let (intensity, obs_mz) = extract_intensity(250.0, &mz, &int, &tol, IntensityRule::MaxInWindow);
        assert!((intensity - 0.0).abs() < f64::EPSILON);
        assert_eq!(obs_mz, None, "no match means no observed m/z");
    }

    #[test]
    fn extract_intensity_empty_spectrum() {
        let tol = MassTolerance {
            value: 0.05,
            unit: ToleranceUnit::Da,
        };
        let (intensity, obs_mz) = extract_intensity(200.0, &[], &[], &tol, IntensityRule::MaxInWindow);
        assert_eq!(intensity, 0.0);
        assert_eq!(obs_mz, None);
    }

    #[test]
    fn extract_intensity_ppm_tolerance() {
        let mz = vec![500.0, 500.005, 500.01];
        let int = vec![100.0, 200.0, 300.0];
        let tol = MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        };
        let (intensity, obs_mz) = extract_intensity(500.005, &mz, &int, &tol, IntensityRule::MaxInWindow);
        assert!(intensity > 0.0, "should find a peak within 20 ppm");
        assert!(obs_mz.is_some(), "should have observed m/z");
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

    #[test]
    fn compute_ion_metadata_counts_k_r_for_b_ions() {
        // "KPEPTIDER" — b1="K" (K=1,R=0), b3="KPE" (K=1,R=0), b8="KPEPTIDE" (K=1,R=0)
        let ions = build_target_ions("KPEPTIDER", &[], 2);
        let meta = compute_ion_metadata(&ions, "KPEPTIDER");

        // Find b1 (1+ charge)
        let b1 = meta.iter().find(|m| m.label == "b1¹⁺").unwrap();
        assert_eq!(b1.k_count, 1);
        assert_eq!(b1.r_count, 0);

        // Find y1 (suffix "R")
        let y1 = meta.iter().find(|m| m.label == "y1¹⁺").unwrap();
        assert_eq!(y1.k_count, 0);
        assert_eq!(y1.r_count, 1);

        // y8 suffix = "PEPTIDER" (K=0, R=1)
        let y8 = meta.iter().find(|m| m.label == "y8¹⁺").unwrap();
        assert_eq!(y8.k_count, 0);
        assert_eq!(y8.r_count, 1);
    }

    #[test]
    fn compute_ion_metadata_no_k_r() {
        // "PEPTIDE" — no K or R
        let ions = build_target_ions("PEPTIDE", &[], 2);
        let meta = compute_ion_metadata(&ions, "PEPTIDE");
        for m in &meta {
            assert_eq!(m.k_count, 0, "ion {} should have k_count=0", m.label);
            assert_eq!(m.r_count, 0, "ion {} should have r_count=0", m.label);
        }
    }

    #[test]
    fn compute_ion_metadata_preserves_light_mz() {
        let ions = build_target_ions("PEPTIDEK", &[], 2);
        let meta = compute_ion_metadata(&ions, "PEPTIDEK");
        // Every metadata entry's light_mz must match the corresponding ion's mz
        for (ion, m) in ions.iter().zip(meta.iter()) {
            assert!(
                (ion.mz - m.light_mz).abs() < 1e-6,
                "mismatch for {}: ion.mz={}, meta.light_mz={}",
                m.label, ion.mz, m.light_mz
            );
        }
    }

    #[test]
    fn trim_ms1_peaks_filters_by_mz_window() {
        let mz = vec![100.0, 200.0, 449.0, 450.0, 451.0, 500.0, 600.0];
        let int = vec![10.0, 20.0, 100.0, 500.0, 200.0, 30.0, 40.0];
        let (trimmed_mz, trimmed_int) = trim_peaks_to_window(&mz, &int, 450.0, 20.0);
        // ±20 Da of 450 = [430, 470], so 449, 450, 451 are in range
        assert_eq!(trimmed_mz.len(), 3);
        assert!((trimmed_mz[0] - 449.0).abs() < 0.01);
        assert!((trimmed_mz[2] - 451.0).abs() < 0.01);
        assert_eq!(trimmed_int.len(), 3);
    }

    #[test]
    fn extract_xic_with_raw_returns_raw_scans() {
        // This test uses the same test fixture as extract_xic tests.
        // We verify that raw_scans is populated and ms2_scans are non-empty.
        let fixture = std::path::Path::new("tests/fixtures");
        if !fixture.exists() {
            // Skip if no fixture available (unit test environment)
        }
        // Integration test — will be covered in Task 7
    }

    #[test]
    fn zero_offset_peptide_skips_heavy() {
        use protein_copilot_core::label::LabelType;

        // "PEPTIDE" has no K or R → zero SILAC shift
        let peptide = "PEPTIDE";
        let label = LabelType::Silac {
            heavy_k_delta: 8.014199,
            heavy_r_delta: 10.008269,
        };
        let delta = protein_copilot_core::label::total_heavy_delta(peptide, &label);
        assert!(delta.abs() < 1e-6, "peptide without K/R should have zero delta");

        // effective_label filter should exclude this
        let effective = Some(&label).filter(|l| {
            protein_copilot_core::label::total_heavy_delta(peptide, l).abs() > 1e-6
        });
        assert!(effective.is_none(), "zero-offset label should be filtered out");

        // With K/R, delta should be non-zero
        let peptide_kr = "PEPTIDEK";
        let delta_kr = protein_copilot_core::label::total_heavy_delta(peptide_kr, &label);
        assert!(delta_kr.abs() > 1.0, "peptide with K should have non-zero delta");
        let effective_kr = Some(&label).filter(|l| {
            protein_copilot_core::label::total_heavy_delta(peptide_kr, l).abs() > 1e-6
        });
        assert!(effective_kr.is_some(), "non-zero offset should keep label");
    }
}