//! Single spectrum annotation for visualization and quality inspection.
//!
//! Given a spectrum and a peptide identification, this module annotates each
//! experimental peak with its matching theoretical b/y ion (if any) and
//! provides a complete [`SpectrumAnnotation`] suitable for downstream
//! visualization or report generation.

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use protein_copilot_core::search_params::{MassTolerance, Modification};
use protein_copilot_core::spectrum::Spectrum;

use crate::chemistry::{peptide_mass, peptide_mz, PROTON_MASS, WATER_MASS};
use crate::error::SearchEngineError;
use crate::matching::{mod_delta_fragment, within_tolerance};

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

/// Type of fragment ion.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub enum IonType {
    /// b-ion (N-terminal fragment).
    B,
    /// y-ion (C-terminal fragment).
    Y,
}

impl std::fmt::Display for IonType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IonType::B => write!(f, "b"),
            IonType::Y => write!(f, "y"),
        }
    }
}

/// Annotation linking an experimental peak to a theoretical fragment ion.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct IonAnnotation {
    /// Type of the matched ion.
    pub ion_type: IonType,
    /// Ion number (1-based; e.g. b3 → 3).
    pub ion_number: u32,
    /// Charge state of the fragment ion (1 for b3¹⁺, 2 for b3²⁺, etc.).
    pub charge: u32,
    /// Theoretical m/z of the ion.
    pub theoretical_mz: f64,
    /// Absolute mass deviation (Da) between experimental and theoretical m/z.
    pub delta_mz: f64,
    /// Relative mass deviation in ppm.
    pub delta_ppm: f64,
}

/// An experimental peak with an optional ion annotation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct AnnotatedPeak {
    /// Observed m/z.
    pub mz: f64,
    /// Observed intensity.
    pub intensity: f64,
    /// Annotation if this peak matched a theoretical ion.
    pub annotation: Option<IonAnnotation>,
}

/// A theoretical ion with match status against the experimental spectrum.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct TheoreticalIon {
    /// Type of the ion (b or y).
    pub ion_type: IonType,
    /// Ion number (1-based).
    pub number: u32,
    /// Charge state (1 for singly-charged, 2 for doubly-charged, etc.).
    pub charge: u32,
    /// Theoretical m/z.
    pub theoretical_mz: f64,
    /// Whether this ion was matched in the experimental spectrum.
    pub matched: bool,
    /// Observed m/z of the matched peak (if matched).
    pub matched_mz: Option<f64>,
    /// Mass deviation in ppm (if matched).
    pub delta_ppm: Option<f64>,
}

/// Complete annotation of a single spectrum against a peptide identification.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct SpectrumAnnotation {
    /// Scan number (1-based).
    pub scan_number: u32,
    /// Retention time in minutes.
    pub retention_time_min: f64,
    /// Source spectrum file name (e.g. "sample.mzML").
    #[serde(default)]
    pub source_file: String,
    /// Identified peptide sequence.
    pub peptide_sequence: String,
    /// Charge state.
    pub charge: i32,
    /// Observed precursor m/z.
    pub precursor_mz: f64,
    /// Theoretical precursor m/z.
    pub theoretical_mz: f64,
    /// Precursor mass deviation in ppm.
    pub delta_mass_ppm: f64,
    /// Match score (matched_ions / total_ions).
    pub score: f64,
    /// Number of matched fragment ions.
    pub matched_ions: u32,
    /// Total number of theoretical fragment ions.
    pub total_ions: u32,
    /// Protein accessions associated with this peptide.
    pub protein_accessions: Vec<String>,
    /// All experimental peaks with optional annotations.
    pub peaks: Vec<AnnotatedPeak>,
    /// Theoretical b-ions with match status.
    pub b_ions: Vec<TheoreticalIon>,
    /// Theoretical y-ions with match status.
    pub y_ions: Vec<TheoreticalIon>,
    /// Fixed modifications applied.
    pub modifications: Vec<Modification>,
    /// Heavy-label annotation from a separate DIA scan (DIA+SILAC only).
    /// `None` for DDA or when no heavy label is configured.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub heavy_annotation: Option<HeavyAnnotation>,
}

/// Annotation data from a heavy-label DIA scan.
///
/// When SILAC + DIA is active, the heavy precursor falls in a different
/// isolation window, so its fragments appear in a different MS2 scan.
/// This struct holds the annotation of that separate scan.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct HeavyAnnotation {
    /// Scan number of the heavy MS2 scan.
    pub scan_number: u32,
    /// Retention time of the heavy scan in minutes.
    pub retention_time_min: f64,
    /// Theoretical heavy precursor m/z (computed from theoretical light + SILAC delta).
    pub precursor_mz: f64,
    /// Precursor mass deviation in ppm (heavy observed vs heavy theoretical).
    /// `None` when observed heavy precursor m/z is unavailable (DIA wide-window).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub delta_mass_ppm: Option<f64>,
    /// All experimental peaks in the heavy scan with optional annotations.
    pub peaks: Vec<AnnotatedPeak>,
    /// Theoretical heavy b-ions with match status.
    pub b_ions: Vec<TheoreticalIon>,
    /// Theoretical heavy y-ions with match status.
    pub y_ions: Vec<TheoreticalIon>,
    /// Match score (matched_ions / total_ions).
    pub score: f64,
    /// Number of matched heavy fragment ions.
    pub matched_ions: u32,
    /// Total number of theoretical heavy fragment ions.
    pub total_ions: u32,
}

// ---------------------------------------------------------------------------
// Helper: binary-search best match
// ---------------------------------------------------------------------------

/// Finds the best matching experimental peak for a theoretical m/z.
///
/// Uses `partition_point` (binary search) on the sorted experimental m/z
/// array and checks neighbors. Returns `(matched, Option<best_index>)`.
fn find_best_match(
    theoretical_mz: f64,
    exp_mz_slice: &[f64],
    tolerance: &MassTolerance,
) -> (bool, Option<usize>) {
    if exp_mz_slice.is_empty() {
        return (false, None);
    }

    // partition_point returns the index of the first element >= theoretical_mz
    let pos = exp_mz_slice.partition_point(|&x| x < theoretical_mz);

    let mut best_idx: Option<usize> = None;
    let mut best_diff = f64::MAX;

    // Check the element at pos and the one before it
    for &candidate in &[pos.wrapping_sub(1), pos] {
        if candidate < exp_mz_slice.len() {
            let diff = (exp_mz_slice[candidate] - theoretical_mz).abs();
            if diff < best_diff
                && within_tolerance(exp_mz_slice[candidate], theoretical_mz, tolerance)
            {
                best_diff = diff;
                best_idx = Some(candidate);
            }
        }
    }

    (best_idx.is_some(), best_idx)
}

// ---------------------------------------------------------------------------
// Fixed-modification mass adjustment
// ---------------------------------------------------------------------------

/// Calculates the total mass shift from fixed modifications for a peptide.
/// Applies fixed modifications to compute total mass delta.
///
/// When `is_protein_nterm` / `is_protein_cterm` are `true`, the corresponding
/// `ProteinNTerm` / `ProteinCTerm` mods are applied; otherwise they are skipped.
fn apply_fixed_mod_mass(
    sequence: &str,
    fixed_mods: &[Modification],
    is_protein_nterm: bool,
    is_protein_cterm: bool,
) -> f64 {
    use protein_copilot_core::search_params::ModPosition;
    let mut delta = 0.0;
    for m in fixed_mods {
        if m.residues.is_empty() {
            match m.position {
                ModPosition::AnyNTerm | ModPosition::AnyCTerm | ModPosition::Anywhere => {
                    delta += m.mass_delta;
                }
                ModPosition::ProteinNTerm => {
                    if is_protein_nterm {
                        delta += m.mass_delta;
                    }
                }
                ModPosition::ProteinCTerm => {
                    if is_protein_cterm {
                        delta += m.mass_delta;
                    }
                }
            }
        } else {
            for ch in sequence.chars() {
                if m.residues.contains(&ch) {
                    delta += m.mass_delta;
                }
            }
        }
    }
    delta
}

// ---------------------------------------------------------------------------
// Fragment ion generation with modification support (multi-charge aware)
// ---------------------------------------------------------------------------

/// A theoretical fragment ion entry before matching.
struct FragmentEntry {
    ion_type: IonType,
    ion_number: u32,
    charge: u32,
    mz: f64,
}

/// Generates theoretical b-ion entries at multiple charge states.
fn generate_b_entries(
    sequence: &str,
    fixed_mods: &[Modification],
    max_charge: u32,
) -> Option<Vec<FragmentEntry>> {
    let chars: Vec<char> = sequence.chars().collect();
    let n = chars.len();
    if n < 2 {
        return Some(Vec::new());
    }
    let max_z = max_charge.max(1) as usize;
    let mut entries = Vec::with_capacity((n - 1) * max_z);
    let mut cumulative = 0.0;

    for (frag_idx, &aa) in chars[..n - 1].iter().enumerate() {
        let mass = crate::chemistry::residue_mass(aa)?;
        cumulative += mass;
        let prefix_len = frag_idx + 1;
        let mod_delta = mod_delta_fragment(&chars[..prefix_len], fixed_mods, true);
        let neutral = cumulative + mod_delta;

        for z in 1..=max_z {
            entries.push(FragmentEntry {
                ion_type: IonType::B,
                ion_number: (frag_idx + 1) as u32,
                charge: z as u32,
                mz: (neutral + z as f64 * PROTON_MASS) / z as f64,
            });
        }
    }
    Some(entries)
}

/// Generates theoretical y-ion entries at multiple charge states.
fn generate_y_entries(
    sequence: &str,
    fixed_mods: &[Modification],
    max_charge: u32,
) -> Option<Vec<FragmentEntry>> {
    let chars: Vec<char> = sequence.chars().collect();
    let n = chars.len();
    if n < 2 {
        return Some(Vec::new());
    }
    let max_z = max_charge.max(1) as usize;
    let mut entries = Vec::with_capacity((n - 1) * max_z);
    let mut cumulative = WATER_MASS;

    for (i, &aa) in chars.iter().rev().enumerate() {
        if i >= n - 1 {
            break;
        }
        let mass = crate::chemistry::residue_mass(aa)?;
        cumulative += mass;
        let suffix_start = n - 1 - i;
        let mod_delta = mod_delta_fragment(&chars[suffix_start..], fixed_mods, false);
        let neutral = cumulative + mod_delta;

        for z in 1..=max_z {
            entries.push(FragmentEntry {
                ion_type: IonType::Y,
                ion_number: (i + 1) as u32,
                charge: z as u32,
                mz: (neutral + z as f64 * PROTON_MASS) / z as f64,
            });
        }
    }
    Some(entries)
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Annotates a spectrum with b/y ion matches against a peptide.
///
/// Validates that the spectrum has peaks and precursor info, generates
/// theoretical fragment ions, matches them against experimental peaks,
/// and returns a complete [`SpectrumAnnotation`].
#[allow(clippy::too_many_arguments)]
pub fn annotate_spectrum(
    spectrum: &Spectrum,
    peptide_sequence: &str,
    charge: i32,
    fragment_tolerance: &MassTolerance,
    fixed_modifications: &[Modification],
    protein_accessions: Vec<String>,
    is_protein_nterm: bool,
    is_protein_cterm: bool,
) -> Result<SpectrumAnnotation, SearchEngineError> {
    // --- Validation ---
    if charge <= 0 {
        return Err(SearchEngineError::ExecutionError {
            detail: format!("charge must be >= 1, got {charge}"),
        });
    }
    if spectrum.mz_array.is_empty() {
        return Err(SearchEngineError::ExecutionError {
            detail: "spectrum has no peaks".to_string(),
        });
    }
    let precursor = spectrum
        .precursors
        .iter()
        .find(|p| p.charge == Some(charge))
        .or(spectrum.precursors.first())
        .ok_or_else(|| SearchEngineError::ExecutionError {
            detail: "spectrum has no precursor information".to_string(),
        })?;

    // --- Theoretical precursor m/z ---
    let neutral_mass =
        peptide_mass(peptide_sequence).ok_or_else(|| SearchEngineError::ExecutionError {
            detail: format!(
                "cannot compute mass for sequence '{}': contains non-standard residue",
                peptide_sequence
            ),
        })?;
    let mod_delta = apply_fixed_mod_mass(
        peptide_sequence,
        fixed_modifications,
        is_protein_nterm,
        is_protein_cterm,
    );
    let modified_mass = neutral_mass + mod_delta;
    let theoretical_precursor_mz = peptide_mz(modified_mass, charge);
    let observed_precursor_mz = precursor.mz;
    let precursor_delta_ppm =
        (observed_precursor_mz - theoretical_precursor_mz) / theoretical_precursor_mz * 1e6;

    // --- Generate theoretical fragment ions (with modifications, multi-charge) ---
    let max_frag_charge: u32 = if charge >= 3 { 2 } else { 1 };

    let b_entries =
        generate_b_entries(peptide_sequence, fixed_modifications, max_frag_charge).ok_or_else(
            || SearchEngineError::ExecutionError {
                detail: format!(
                    "cannot generate b-ions for '{}': non-standard residue",
                    peptide_sequence
                ),
            },
        )?;
    let y_entries =
        generate_y_entries(peptide_sequence, fixed_modifications, max_frag_charge).ok_or_else(
            || SearchEngineError::ExecutionError {
                detail: format!(
                    "cannot generate y-ions for '{}': non-standard residue",
                    peptide_sequence
                ),
            },
        )?;

    let exp_mz = &spectrum.mz_array;
    let exp_int = &spectrum.intensity_array;

    // --- Match theoretical ions against experimental peaks ---
    let mut peak_annotations: Vec<Option<IonAnnotation>> = vec![None; exp_mz.len()];

    // Helper closure: match entries and build TheoreticalIon list
    let mut match_entries =
        |entries: &[FragmentEntry]| -> Vec<TheoreticalIon> {
            let mut ions_out = Vec::with_capacity(entries.len());
            for entry in entries {
                let (_, best_idx) = find_best_match(entry.mz, exp_mz, fragment_tolerance);

                if let Some(idx) = best_idx {
                    let obs_mz = exp_mz[idx];
                    let dppm = (obs_mz - entry.mz) / entry.mz * 1e6;

                    let new_annotation = IonAnnotation {
                        ion_type: entry.ion_type,
                        ion_number: entry.ion_number,
                        charge: entry.charge,
                        theoretical_mz: entry.mz,
                        delta_mz: obs_mz - entry.mz,
                        delta_ppm: dppm,
                    };
                    match &peak_annotations[idx] {
                        Some(existing)
                            if existing.delta_mz.abs() <= new_annotation.delta_mz.abs() => {}
                        _ => {
                            peak_annotations[idx] = Some(new_annotation);
                        }
                    }

                    ions_out.push(TheoreticalIon {
                        ion_type: entry.ion_type,
                        number: entry.ion_number,
                        charge: entry.charge,
                        theoretical_mz: entry.mz,
                        matched: true,
                        matched_mz: Some(obs_mz),
                        delta_ppm: Some(dppm),
                    });
                } else {
                    ions_out.push(TheoreticalIon {
                        ion_type: entry.ion_type,
                        number: entry.ion_number,
                        charge: entry.charge,
                        theoretical_mz: entry.mz,
                        matched: false,
                        matched_mz: None,
                        delta_ppm: None,
                    });
                }
            }
            ions_out
        };

    let b_ions_out = match_entries(&b_entries);
    let y_ions_out = match_entries(&y_entries);

    // --- Build annotated peaks ---
    let peaks: Vec<AnnotatedPeak> = exp_mz
        .iter()
        .zip(exp_int.iter())
        .zip(peak_annotations)
        .map(|((&mz, &intensity), annotation)| AnnotatedPeak {
            mz,
            intensity,
            annotation,
        })
        .collect();

    // Count matched ions from the TheoreticalIon lists (avoids double-counting
    // when both a b-ion and y-ion match the same experimental peak).
    let matched_b = b_ions_out.iter().filter(|i| i.matched).count() as u32;
    let matched_y = y_ions_out.iter().filter(|i| i.matched).count() as u32;
    let matched_count = matched_b + matched_y;
    let total_ions = (b_entries.len() + y_entries.len()) as u32;
    let score = if total_ions > 0 {
        matched_count as f64 / total_ions as f64
    } else {
        0.0
    };

    Ok(SpectrumAnnotation {
        scan_number: spectrum.scan_number,
        retention_time_min: spectrum.retention_time_min,
        source_file: String::new(),
        peptide_sequence: peptide_sequence.to_string(),
        charge,
        precursor_mz: observed_precursor_mz,
        theoretical_mz: theoretical_precursor_mz,
        delta_mass_ppm: precursor_delta_ppm,
        score,
        matched_ions: matched_count,
        total_ions,
        protein_accessions,
        peaks,
        b_ions: b_ions_out,
        y_ions: y_ions_out,
        modifications: fixed_modifications.to_vec(),
        heavy_annotation: None,
    })
}

/// Annotate a heavy-label spectrum for mirror plot display.
///
/// Given a heavy MS2 spectrum and the peptide sequence + label type,
/// computes heavy theoretical fragment m/z values and matches them
/// against the heavy spectrum's peaks.
///
/// This is called separately from `annotate_spectrum` because in DIA mode,
/// the heavy fragments appear in a different MS2 scan.
#[allow(clippy::too_many_arguments)]
pub fn annotate_heavy_spectrum(
    heavy_spectrum: &Spectrum,
    peptide_sequence: &str,
    charge: i32,
    fragment_tolerance: &MassTolerance,
    fixed_modifications: &[Modification],
    label: &protein_copilot_core::label::LabelType,
    is_protein_nterm: bool,
    is_protein_cterm: bool,
) -> Result<HeavyAnnotation, SearchEngineError> {
    use protein_copilot_core::label::{compute_heavy_precursor_mz, residue_heavy_delta};

    if charge <= 0 {
        return Err(SearchEngineError::ExecutionError {
            detail: format!("charge must be >= 1, got {charge}"),
        });
    }
    if heavy_spectrum.mz_array.is_empty() {
        return Err(SearchEngineError::ExecutionError {
            detail: "heavy spectrum has no peaks".to_string(),
        });
    }

    // Compute heavy precursor m/z
    let light_precursor_mz = {
        let neutral = peptide_mass(peptide_sequence).ok_or_else(|| {
            SearchEngineError::ExecutionError {
                detail: format!("cannot compute mass for '{}'", peptide_sequence),
            }
        })?;
        let mod_delta = apply_fixed_mod_mass(
            peptide_sequence,
            fixed_modifications,
            is_protein_nterm,
            is_protein_cterm,
        );
        peptide_mz(neutral + mod_delta, charge)
    };
    let heavy_precursor_mz =
        compute_heavy_precursor_mz(light_precursor_mz, charge, peptide_sequence, label);

    // Compute delta_mass_ppm if the heavy spectrum has a precursor m/z
    let delta_mass_ppm = heavy_spectrum
        .precursors
        .first()
        .map(|p| (p.mz - heavy_precursor_mz) / heavy_precursor_mz * 1e6);

    // Generate heavy theoretical fragments
    let max_frag_charge: u32 = if charge >= 3 { 2 } else { 1 };
    let chars: Vec<char> = peptide_sequence.chars().collect();
    let n = chars.len();

    // Heavy b-ions
    let light_b = generate_b_entries(peptide_sequence, fixed_modifications, max_frag_charge)
        .ok_or_else(|| SearchEngineError::ExecutionError {
            detail: format!("cannot generate b-ions for '{}'", peptide_sequence),
        })?;
    let heavy_b_entries: Vec<FragmentEntry> = light_b
        .iter()
        .map(|e| {
            let prefix = &chars[..(e.ion_number as usize).min(n)];
            let delta = residue_heavy_delta(prefix, label);
            FragmentEntry {
                ion_type: e.ion_type,
                ion_number: e.ion_number,
                charge: e.charge,
                mz: e.mz + delta / e.charge as f64,
            }
        })
        .collect();

    // Heavy y-ions
    let light_y = generate_y_entries(peptide_sequence, fixed_modifications, max_frag_charge)
        .ok_or_else(|| SearchEngineError::ExecutionError {
            detail: format!("cannot generate y-ions for '{}'", peptide_sequence),
        })?;
    let heavy_y_entries: Vec<FragmentEntry> = light_y
        .iter()
        .map(|e| {
            let start = n.saturating_sub(e.ion_number as usize);
            let suffix = &chars[start..];
            let delta = residue_heavy_delta(suffix, label);
            FragmentEntry {
                ion_type: e.ion_type,
                ion_number: e.ion_number,
                charge: e.charge,
                mz: e.mz + delta / e.charge as f64,
            }
        })
        .collect();

    // Match against heavy spectrum peaks
    let exp_mz = &heavy_spectrum.mz_array;
    let exp_int = &heavy_spectrum.intensity_array;
    let mut peak_annotations: Vec<Option<IonAnnotation>> = vec![None; exp_mz.len()];
    let total_count = (heavy_b_entries.len() + heavy_y_entries.len()) as u32;

    // Reuse the same matching approach as annotate_spectrum
    let mut match_heavy_entries = |entries: &[FragmentEntry]| -> Vec<TheoreticalIon> {
        let mut ions_out = Vec::with_capacity(entries.len());
        for entry in entries {
            let (_, best_idx) = find_best_match(entry.mz, exp_mz, fragment_tolerance);
            if let Some(idx) = best_idx {
                let obs_mz = exp_mz[idx];
                let dppm = (obs_mz - entry.mz) / entry.mz * 1e6;
                let new_ann = IonAnnotation {
                    ion_type: entry.ion_type,
                    ion_number: entry.ion_number,
                    charge: entry.charge,
                    theoretical_mz: entry.mz,
                    delta_mz: obs_mz - entry.mz,
                    delta_ppm: dppm,
                };
                match &peak_annotations[idx] {
                    Some(existing) if existing.delta_mz.abs() <= new_ann.delta_mz.abs() => {}
                    _ => {
                        peak_annotations[idx] = Some(new_ann);
                    }
                }
                ions_out.push(TheoreticalIon {
                    ion_type: entry.ion_type,
                    number: entry.ion_number,
                    charge: entry.charge,
                    theoretical_mz: entry.mz,
                    matched: true,
                    matched_mz: Some(obs_mz),
                    delta_ppm: Some(dppm),
                });
            } else {
                ions_out.push(TheoreticalIon {
                    ion_type: entry.ion_type,
                    number: entry.ion_number,
                    charge: entry.charge,
                    theoretical_mz: entry.mz,
                    matched: false,
                    matched_mz: None,
                    delta_ppm: None,
                });
            }
        }
        ions_out
    };

    let heavy_b_ions = match_heavy_entries(&heavy_b_entries);
    let heavy_y_ions = match_heavy_entries(&heavy_y_entries);

    // Build annotated peaks
    let peaks: Vec<AnnotatedPeak> = exp_mz
        .iter()
        .zip(exp_int.iter())
        .zip(peak_annotations)
        .map(|((&mz, &intensity), annotation)| AnnotatedPeak {
            mz,
            intensity,
            annotation,
        })
        .collect();

    let matched_b = heavy_b_ions.iter().filter(|i| i.matched).count() as u32;
    let matched_y = heavy_y_ions.iter().filter(|i| i.matched).count() as u32;
    let matched_count = matched_b + matched_y;
    let score = if total_count > 0 {
        matched_count as f64 / total_count as f64
    } else {
        0.0
    };

    Ok(HeavyAnnotation {
        scan_number: heavy_spectrum.scan_number,
        retention_time_min: heavy_spectrum.retention_time_min,
        precursor_mz: heavy_precursor_mz,
        delta_mass_ppm,
        peaks,
        b_ions: heavy_b_ions,
        y_ions: heavy_y_ions,
        score,
        matched_ions: matched_count,
        total_ions: total_count,
    })
}
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matching::{generate_b_ions, generate_y_ions};
    use protein_copilot_core::search_params::ToleranceUnit;
    use protein_copilot_core::spectrum::{MsLevel, PrecursorInfo, Spectrum};

    fn fragment_tolerance_da() -> MassTolerance {
        MassTolerance {
            value: 0.02,
            unit: ToleranceUnit::Da,
        }
    }

    fn make_spectrum(precursor_mz: f64, charge: Option<i32>, peaks_mz: Vec<f64>) -> Spectrum {
        let intensities = vec![1000.0; peaks_mz.len()];
        Spectrum::new(
            1,
            MsLevel::MS2,
            120.5,
            vec![PrecursorInfo {
                mz: precursor_mz,
                charge,
                intensity: None,
                isolation_window: None,
                source_scan: None,
            }],
            peaks_mz,
            intensities,
        )
        .expect("test spectrum should be valid")
    }

    #[test]
    fn annotate_matches_all_ions() {
        let seq = "PEPTIDER";
        let charge = 2;
        let mass = peptide_mass(seq).unwrap();
        let prec_mz = peptide_mz(mass, charge);

        // Build experimental peaks from exact b/y ion positions
        let b = generate_b_ions(seq, &[]);
        let y = generate_y_ions(seq, &[]);
        let mut peaks: Vec<f64> = b.iter().chain(y.iter()).copied().collect();
        peaks.sort_by(|a, b| a.total_cmp(b));
        peaks.dedup();

        let spectrum = make_spectrum(prec_mz, Some(charge), peaks);
        let ann = annotate_spectrum(
            &spectrum,
            seq,
            charge,
            &fragment_tolerance_da(),
            &[],
            vec!["P001".to_string()],
            false,
            false,
        )
        .unwrap();

        assert_eq!(ann.peptide_sequence, seq);
        assert_eq!(ann.charge, charge);
        assert!(
            ann.score > 0.9,
            "expected high score with perfect peaks, got {}",
            ann.score
        );
        assert_eq!(ann.matched_ions, ann.total_ions);
        assert!(ann.delta_mass_ppm.abs() < 1.0);

        // All b-ions matched
        for bi in &ann.b_ions {
            assert!(bi.matched, "b{} should be matched", bi.number);
        }
        // All y-ions matched
        for yi in &ann.y_ions {
            assert!(yi.matched, "y{} should be matched", yi.number);
        }
    }

    #[test]
    fn annotate_no_match() {
        let wrong_seq = "AAAAAAAAK";
        let charge = 2;

        // Annotate against a completely different peptide with random peaks
        let wrong_mass = peptide_mass(wrong_seq).unwrap();
        let wrong_prec_mz = peptide_mz(wrong_mass, charge);
        let spectrum_wrong = make_spectrum(wrong_prec_mz, Some(charge), vec![100.0, 200.0, 300.0]);

        let ann = annotate_spectrum(
            &spectrum_wrong,
            wrong_seq,
            charge,
            &fragment_tolerance_da(),
            &[],
            vec![],
            false,
            false,
        )
        .unwrap();

        // Very few or no matches with random peaks
        assert!(
            ann.score < 0.3,
            "expected low score for mismatched peaks, got {}",
            ann.score
        );
    }

    #[test]
    fn annotate_nonstandard_residue_errors() {
        let seq = "PEPT*DE";
        let charge = 2;
        // Use a dummy precursor mz
        let spectrum = make_spectrum(500.0, Some(charge), vec![100.0, 200.0, 300.0]);

        let result = annotate_spectrum(
            &spectrum,
            seq,
            charge,
            &fragment_tolerance_da(),
            &[],
            vec![],
            false,
            false,
        );
        let err_msg = result.unwrap_err().to_string();
        assert!(
            err_msg.contains("non-standard"),
            "error should mention non-standard residue: {}",
            err_msg
        );
    }

    #[test]
    fn annotated_peaks_include_noise() {
        let seq = "PEPTIDER";
        let charge = 2;
        let mass = peptide_mass(seq).unwrap();
        let prec_mz = peptide_mz(mass, charge);

        let b = generate_b_ions(seq, &[]);
        // Add noise peaks that don't match any theoretical ion
        let noise_peaks = [50.0, 150.5, 850.0, 1200.0];
        let mut peaks: Vec<f64> = b
            .iter()
            .copied()
            .chain(noise_peaks.iter().copied())
            .collect();
        peaks.sort_by(|a, b| a.total_cmp(b));

        let spectrum = make_spectrum(prec_mz, Some(charge), peaks.clone());
        let ann = annotate_spectrum(
            &spectrum,
            seq,
            charge,
            &fragment_tolerance_da(),
            &[],
            vec![],
            false,
            false,
        )
        .unwrap();

        // Total peaks should include noise
        assert_eq!(ann.peaks.len(), peaks.len());

        // Count unannotated peaks — should be at least the noise peaks
        let unannotated = ann.peaks.iter().filter(|p| p.annotation.is_none()).count();
        assert!(
            unannotated >= noise_peaks.len(),
            "expected at least {} unannotated noise peaks, got {}",
            noise_peaks.len(),
            unannotated
        );
    }

    #[test]
    fn annotation_serde_roundtrip() {
        let seq = "PEPTIDER";
        let charge = 2;
        let mass = peptide_mass(seq).unwrap();
        let prec_mz = peptide_mz(mass, charge);

        let b = generate_b_ions(seq, &[]);
        let y = generate_y_ions(seq, &[]);
        let mut peaks: Vec<f64> = b.iter().chain(y.iter()).copied().collect();
        peaks.sort_by(|a, b| a.total_cmp(b));
        peaks.dedup();

        let spectrum = make_spectrum(prec_mz, Some(charge), peaks);
        let ann = annotate_spectrum(
            &spectrum,
            seq,
            charge,
            &fragment_tolerance_da(),
            &[],
            vec!["PROT1".to_string()],
            false,
            false,
        )
        .unwrap();

        let json = serde_json::to_string_pretty(&ann).unwrap();
        let deserialized: SpectrumAnnotation = serde_json::from_str(&json).unwrap();

        assert_eq!(ann, deserialized);
    }

    #[test]
    fn annotate_selects_precursor_by_charge() {
        // Two precursors: charge 2 at 500.0, charge 3 at 400.0
        let spectrum = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            vec![
                PrecursorInfo {
                    mz: 500.0,
                    charge: Some(2),
                    intensity: None,
                    isolation_window: None,
                    source_scan: None,
                },
                PrecursorInfo {
                    mz: 400.0,
                    charge: Some(3),
                    intensity: None,
                    isolation_window: None,
                    source_scan: None,
                },
            ],
            vec![100.0, 200.0, 300.0],
            vec![1000.0, 2000.0, 500.0],
        )
        .unwrap();

        let tol = MassTolerance {
            value: 0.5,
            unit: ToleranceUnit::Da,
        };
        // Annotate with charge 3 — should pick precursor at 400.0, not 500.0
        let result = annotate_spectrum(&spectrum, "GK", 3, &tol, &[], vec![], false, false);
        assert!(result.is_ok());
        let ann = result.unwrap();
        assert!(
            (ann.precursor_mz - 400.0).abs() < 0.01,
            "should select precursor with matching charge 3, got {}",
            ann.precursor_mz
        );
    }

    #[test]
    fn ion_type_display() {
        assert_eq!(IonType::B.to_string(), "b");
        assert_eq!(IonType::Y.to_string(), "y");
    }

    #[test]
    fn annotate_empty_spectrum_errors() {
        let spectrum = Spectrum::new(
            1,
            MsLevel::MS2,
            100.0,
            vec![PrecursorInfo {
                mz: 500.0,
                charge: Some(2),
                intensity: None,
                isolation_window: None,
                source_scan: None,
            }],
            vec![],
            vec![],
        );
        // Spectrum::new may reject empty arrays; if so, build manually
        // or test with a spectrum that has no precursor.
        if let Ok(spec) = spectrum {
            let result =
                annotate_spectrum(&spec, "PEPTIDER", 2, &fragment_tolerance_da(), &[], vec![], false, false);
            assert!(result.is_err());
        }
        // Also test missing precursor
    }

    #[test]
    fn nterm_modification_affects_b_ions_not_y_ions() {
        use protein_copilot_core::search_params::ModPosition;

        let seq = "PEPTIDE";
        // TMT6plex-like N-term modification (229.163 Da)
        let tmt_mod = Modification {
            name: "TMT6plex".to_string(),
            mass_delta: 229.162932,
            residues: vec![],
            position: ModPosition::AnyNTerm,
        };

        // Generate b-ions WITH and WITHOUT the N-term mod
        let b_no_mod: Vec<f64> = generate_b_entries(seq, &[], 1)
            .unwrap()
            .iter()
            .map(|e| e.mz)
            .collect();
        let b_with_mod: Vec<f64> =
            generate_b_entries(seq, std::slice::from_ref(&tmt_mod), 1)
                .unwrap()
                .iter()
                .map(|e| e.mz)
                .collect();

        // All b-ions should be shifted by the TMT mass
        assert_eq!(b_no_mod.len(), b_with_mod.len());
        for (no_mod, with_mod) in b_no_mod.iter().zip(b_with_mod.iter()) {
            let diff = with_mod - no_mod;
            assert!(
                (diff - 229.162932).abs() < 0.001,
                "b-ion should be shifted by TMT mass, got diff {diff:.4}"
            );
        }

        // y-ions should NOT be affected by N-term mod
        let y_no_mod: Vec<f64> = generate_y_entries(seq, &[], 1)
            .unwrap()
            .iter()
            .map(|e| e.mz)
            .collect();
        let y_with_mod: Vec<f64> = generate_y_entries(seq, &[tmt_mod], 1)
            .unwrap()
            .iter()
            .map(|e| e.mz)
            .collect();

        assert_eq!(y_no_mod.len(), y_with_mod.len());
        for (no_mod, with_mod) in y_no_mod.iter().zip(y_with_mod.iter()) {
            let diff = (with_mod - no_mod).abs();
            assert!(
                diff < 0.001,
                "y-ion should NOT be shifted by N-term mod, got diff {diff:.4}"
            );
        }
    }

    #[test]
    fn test_annotate_heavy_spectrum_basic() {
        use protein_copilot_core::label::LabelType;
        use protein_copilot_core::search_params::ToleranceUnit;

        // Heavy spectrum with a peak near heavy b2+ m/z for "AK" prefix
        // AK: A=71.03711 + K=128.09496 = 199.13207 neutral b2
        // Heavy K delta = 8.014199 → heavy b2 neutral = 207.146269
        // b2+ m/z = (207.146269 + 1.007276) / 1 = 208.153545
        let heavy_b2_mz = 208.1535;

        let spectrum = Spectrum::new(
            100,
            MsLevel::MS2,
            25.0,
            vec![PrecursorInfo {
                mz: 300.0,
                charge: Some(2),
                intensity: Some(5000.0),
                isolation_window: None,
                source_scan: None,
            }],
            vec![heavy_b2_mz, 300.0, 500.0],
            vec![1000.0, 200.0, 150.0],
        )
        .expect("test spectrum");

        let label = LabelType::standard_silac();
        let tolerance = MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        };

        let result =
            annotate_heavy_spectrum(&spectrum, "AKDEF", 2, &tolerance, &[], &label, false, false).unwrap();

        assert_eq!(result.scan_number, 100);
        assert!(result.matched_ions >= 1, "should match at least heavy b2+");
        assert!(result.total_ions > 0);
        assert!(result.score > 0.0);
        let matched_b = result.b_ions.iter().filter(|i| i.matched).count();
        assert!(
            matched_b >= 1,
            "should have at least one matched heavy b-ion"
        );
    }

    #[test]
    fn test_terminal_mod_applied_when_context_true() {
        use protein_copilot_core::search_params::{ModPosition, Modification};
        let nterm_mod = Modification {
            name: "Acetyl".to_string(),
            mass_delta: 42.010565,
            residues: vec![],
            position: ModPosition::ProteinNTerm,
        };
        // With is_protein_nterm=true, mod should be applied
        let delta_applied = apply_fixed_mod_mass("PEPTIDE", &[nterm_mod.clone()], true, false);
        assert!((delta_applied - 42.010565).abs() < 1e-6);

        // With is_protein_nterm=false, mod should be skipped
        let delta_skipped = apply_fixed_mod_mass("PEPTIDE", &[nterm_mod], false, false);
        assert!(delta_skipped.abs() < 1e-6);
    }
}
