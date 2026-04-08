//! Spectrum-to-peptide matching and scoring.
//!
//! Implements a simplified search algorithm:
//! 1. **Precursor matching**: compare observed precursor m/z with theoretical
//!    peptide m/z within a given tolerance (ppm or Da)
//! 2. **Fragment matching**: generate theoretical b/y ions and count matches
//!    against the experimental peak list
//! 3. **Scoring**: matched fragments / total theoretical fragments
//!
//! This is a simplified scoring model for MVP validation. Production engines
//! (pFind, MSFragger) use statistical scoring models (e.g., hyperscore,
//! binomial probability).

use protein_copilot_core::search_params::{MassTolerance, Modification, ToleranceUnit};
use protein_copilot_core::spectrum::Spectrum;

use crate::chemistry::{residue_mass, PROTON_MASS, WATER_MASS};
use crate::digest::{peptide_mz, DigestedPeptide};

/// A candidate match between a spectrum and a peptide.
#[derive(Debug, Clone)]
pub struct PeptideMatch {
    /// The matched peptide.
    pub peptide: DigestedPeptide,
    /// Charge state used for matching.
    pub charge: i32,
    /// Observed precursor m/z from the spectrum.
    pub observed_mz: f64,
    /// Theoretical m/z for this peptide at this charge.
    pub theoretical_mz: f64,
    /// Mass deviation in ppm.
    pub delta_mass_ppm: f64,
    /// Match score (0.0 to 1.0): matched fragments / total fragments.
    pub score: f64,
    /// Number of matched b/y ions.
    pub matched_ions: u32,
    /// Total theoretical b/y ions.
    pub total_ions: u32,
    /// Variable modifications applied in this match.
    /// Each entry is (Modification, 0-based residue position or usize::MAX for terminal).
    pub applied_variable_mods: Vec<(Modification, usize)>,
}

// ---------------------------------------------------------------------------
// Precursor matching
// ---------------------------------------------------------------------------

/// Checks if observed and theoretical m/z are within tolerance.
pub fn within_tolerance(observed: f64, theoretical: f64, tolerance: &MassTolerance) -> bool {
    match tolerance.unit {
        ToleranceUnit::Ppm => {
            let ppm_diff = ((observed - theoretical) / theoretical).abs() * 1e6;
            ppm_diff <= tolerance.value
        }
        ToleranceUnit::Da => (observed - theoretical).abs() <= tolerance.value,
    }
}

/// Calculates mass deviation in ppm.
fn calc_delta_ppm(observed: f64, theoretical: f64) -> f64 {
    (observed - theoretical) / theoretical * 1e6
}

// ---------------------------------------------------------------------------
// Fragment ion generation (b/y ions)
// ---------------------------------------------------------------------------

/// Cumulative fixed-modification mass delta for a fragment prefix or suffix.
///
/// For residue-specific mods: adds `mass_delta` for each matching residue.
/// For N-term mods: adds `mass_delta` to b-ions (which contain the N-terminus).
/// For C-term mods: adds `mass_delta` to y-ions (which contain the C-terminus).
pub(crate) fn mod_delta_fragment(
    residues: &[char],
    fixed_mods: &[Modification],
    is_b_ion: bool,
) -> f64 {
    use protein_copilot_core::search_params::ModPosition;
    let mut delta = 0.0;
    for m in fixed_mods {
        if m.residues.is_empty() {
            match m.position {
                ModPosition::AnyNTerm | ModPosition::ProteinNTerm => {
                    if is_b_ion {
                        delta += m.mass_delta;
                    }
                }
                ModPosition::AnyCTerm | ModPosition::ProteinCTerm => {
                    if !is_b_ion {
                        delta += m.mass_delta;
                    }
                }
                ModPosition::Anywhere => {
                    delta += m.mass_delta;
                }
            }
        } else {
            for &ch in residues {
                if m.residues.contains(&ch) {
                    delta += m.mass_delta;
                }
            }
        }
    }
    delta
}

/// Generates theoretical singly-charged b-ion m/z values for a peptide,
/// applying fixed modifications to each fragment.
pub fn generate_b_ions(sequence: &str, fixed_mods: &[Modification]) -> Vec<f64> {
    let chars: Vec<char> = sequence.chars().collect();
    let mut ions = Vec::with_capacity(chars.len().saturating_sub(1));
    let mut cumulative = 0.0;

    for &aa in &chars[..chars.len().saturating_sub(1)] {
        // Skip entire generation if non-standard residue encountered
        let mass = match residue_mass(aa) {
            Some(m) => m,
            None => return Vec::new(),
        };
        cumulative += mass;
        let prefix_len = ions.len() + 1;
        let mod_delta = mod_delta_fragment(&chars[..prefix_len], fixed_mods, true);
        // b-ion m/z (singly charged) = cumulative residues + mod delta + proton
        ions.push(cumulative + mod_delta + PROTON_MASS);
    }

    ions
}

/// Generates theoretical singly-charged y-ion m/z values for a peptide,
/// applying fixed modifications to each fragment.
pub fn generate_y_ions(sequence: &str, fixed_mods: &[Modification]) -> Vec<f64> {
    let chars: Vec<char> = sequence.chars().collect();
    let n = chars.len();
    let mut ions = Vec::with_capacity(n.saturating_sub(1));
    let mut cumulative = WATER_MASS;

    for (i, &aa) in chars.iter().rev().enumerate() {
        if i >= n.saturating_sub(1) {
            break;
        }
        let mass = match residue_mass(aa) {
            Some(m) => m,
            None => return Vec::new(),
        };
        cumulative += mass;
        let suffix_start = n - 1 - i;
        let mod_delta = mod_delta_fragment(&chars[suffix_start..], fixed_mods, false);
        // y-ion m/z (singly charged) = cumulative residues + water + mod delta + proton
        ions.push(cumulative + mod_delta + PROTON_MASS);
    }

    ions
}

// ---------------------------------------------------------------------------
// Fragment matching
// ---------------------------------------------------------------------------

/// Counts how many theoretical ions match peaks in the experimental spectrum.
fn count_matched_ions(
    theoretical_ions: &[f64],
    experimental_mz: &[f64],
    tolerance: &MassTolerance,
) -> u32 {
    let mut matched = 0u32;

    for &theo in theoretical_ions {
        // Binary search for closest peak in sorted experimental array.
        // total_cmp handles NaN deterministically (NaN sorts after all values).
        let idx = experimental_mz.binary_search_by(|probe| probe.total_cmp(&theo));

        let candidates = match idx {
            Ok(i) => vec![i],
            Err(i) => {
                let mut c = Vec::new();
                if i > 0 {
                    c.push(i - 1);
                }
                if i < experimental_mz.len() {
                    c.push(i);
                }
                c
            }
        };

        for &ci in &candidates {
            if within_tolerance(experimental_mz[ci], theo, tolerance) {
                matched += 1;
                break;
            }
        }
    }

    matched
}

// ---------------------------------------------------------------------------
// Public API: match spectrum against candidate peptides
// ---------------------------------------------------------------------------

/// Applies fixed modifications to the theoretical peptide mass.
///
/// Terminal mods (empty residues) are applied based on their `ModPosition`:
/// - `AnyNTerm` / `ProteinNTerm`: apply once (N-terminal)
/// - `AnyCTerm` / `ProteinCTerm`: apply once (C-terminal)
/// - `Anywhere` with empty residues: treated as a global mod, applied once
///
/// Residue-specific mods: apply once per matching residue in the sequence.
fn apply_fixed_mods(sequence: &str, mods: &[Modification]) -> f64 {
    use protein_copilot_core::search_params::ModPosition;
    let mut delta = 0.0;
    for m in mods {
        if m.residues.is_empty() {
            match m.position {
                ModPosition::AnyNTerm | ModPosition::ProteinNTerm |
                ModPosition::AnyCTerm | ModPosition::ProteinCTerm |
                ModPosition::Anywhere => {
                    // Terminal or global mod with no specific residue: apply once
                    delta += m.mass_delta;
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

/// Builds a combined modification list merging fixed mods with applied variable sites.
///
/// Variable mods at residue positions are converted to `Anywhere` single-residue
/// mods for fragment ion mass calculation via `mod_delta_fragment`.
fn build_combined_mods(
    fixed_mods: &[Modification],
    variable_mods: &[Modification],
    applied_sites: &[crate::varmod::ModSite],
    sequence: &str,
) -> Vec<Modification> {
    use protein_copilot_core::search_params::ModPosition;
    let chars: Vec<char> = sequence.chars().collect();
    let mut combined: Vec<Modification> = fixed_mods.to_vec();

    for site in applied_sites {
        let base_mod = &variable_mods[site.mod_index];
        if site.residue_pos == usize::MAX {
            combined.push(base_mod.clone());
        } else {
            let ch = chars[site.residue_pos];
            combined.push(Modification {
                name: base_mod.name.clone(),
                mass_delta: base_mod.mass_delta,
                residues: vec![ch],
                position: ModPosition::Anywhere,
            });
        }
    }

    combined
}

/// Matches a single spectrum against all candidate peptides.
///
/// Returns the best match (highest score) if any peptide matches
/// within the precursor tolerance.
pub fn match_spectrum(
    spectrum: &Spectrum,
    candidates: &[DigestedPeptide],
    precursor_tolerance: &MassTolerance,
    fragment_tolerance: &MassTolerance,
    fixed_mods: &[Modification],
    variable_mods: &[Modification],
    max_variable_mods: u32,
) -> Option<PeptideMatch> {
    // Need at least one precursor to match against
    let precursor = spectrum.precursors.first()?;
    let observed_mz = precursor.mz;

    // Try charge states 1-4 (or use observed charge if available)
    let charge_states: Vec<i32> = if let Some(c) = precursor.charge {
        vec![c]
    } else {
        vec![2, 3, 1, 4] // common charge states, ordered by frequency
    };

    let mut best_match: Option<PeptideMatch> = None;

    for peptide in candidates {
        let fixed_delta = apply_fixed_mods(&peptide.sequence, fixed_mods);

        // Discover applicable variable mod sites for this peptide
        let sites = crate::varmod::find_applicable_sites(
            &peptide.sequence,
            variable_mods,
            peptide.is_protein_nterm,
            peptide.is_protein_cterm,
        );

        let combinations =
            crate::varmod::enumerate_combinations(variable_mods, &sites, max_variable_mods);

        for combo in &combinations {
            let modified_mass = peptide.neutral_mass + fixed_delta + combo.mass_delta;

            for &charge in &charge_states {
                if charge == 0 {
                    continue;
                }
                let theoretical_mz = peptide_mz(modified_mass, charge);

                if within_tolerance(observed_mz, theoretical_mz, precursor_tolerance) {
                    let combined_mods = build_combined_mods(
                        fixed_mods,
                        variable_mods,
                        &combo.sites,
                        &peptide.sequence,
                    );

                    let b_ions = generate_b_ions(&peptide.sequence, &combined_mods);
                    let y_ions = generate_y_ions(&peptide.sequence, &combined_mods);

                    let total_theoretical = (b_ions.len() + y_ions.len()) as u32;
                    if total_theoretical == 0 {
                        continue;
                    }

                    let all_ions: Vec<f64> = b_ions.into_iter().chain(y_ions).collect();
                    let matched =
                        count_matched_ions(&all_ions, &spectrum.mz_array, fragment_tolerance);

                    let score = matched as f64 / total_theoretical as f64;
                    let delta_ppm = calc_delta_ppm(observed_mz, theoretical_mz);

                    if !score.is_finite() || !delta_ppm.is_finite() {
                        continue;
                    }

                    let is_better = match &best_match {
                        None => true,
                        Some(prev) => score > prev.score,
                    };

                    if is_better {
                        let applied = combo
                            .sites
                            .iter()
                            .map(|s| (variable_mods[s.mod_index].clone(), s.residue_pos))
                            .collect();

                        best_match = Some(PeptideMatch {
                            peptide: peptide.clone(),
                            charge,
                            observed_mz,
                            theoretical_mz,
                            delta_mass_ppm: delta_ppm,
                            score,
                            matched_ions: matched,
                            total_ions: total_theoretical,
                            applied_variable_mods: applied,
                        });
                    }
                }
            }
        }
    }

    best_match
}

/// Matches a spectrum against candidate peptides using ALL precursors.
///
/// Unlike [`match_spectrum`] which only uses the first precursor,
/// this function iterates over every precursor in `spectrum.precursors`
/// and returns the best match for each precursor that produces a hit.
///
/// Designed for DIA data where a spectrum can have multiple candidate
/// precursors after extraction.
pub fn match_spectrum_all(
    spectrum: &Spectrum,
    candidates: &[DigestedPeptide],
    precursor_tolerance: &MassTolerance,
    fragment_tolerance: &MassTolerance,
    fixed_mods: &[Modification],
    variable_mods: &[Modification],
    max_variable_mods: u32,
) -> Vec<PeptideMatch> {
    let mut results = Vec::new();

    for precursor in &spectrum.precursors {
        let observed_mz = precursor.mz;

        let charge_states: Vec<i32> = if let Some(c) = precursor.charge {
            vec![c]
        } else {
            vec![2, 3, 1, 4]
        };

        let mut best_match: Option<PeptideMatch> = None;

        for peptide in candidates {
            let fixed_delta = apply_fixed_mods(&peptide.sequence, fixed_mods);

            let sites = crate::varmod::find_applicable_sites(
                &peptide.sequence,
                variable_mods,
                peptide.is_protein_nterm,
                peptide.is_protein_cterm,
            );

            let combinations =
                crate::varmod::enumerate_combinations(variable_mods, &sites, max_variable_mods);

            for combo in &combinations {
                let modified_mass = peptide.neutral_mass + fixed_delta + combo.mass_delta;

                for &charge in &charge_states {
                    if charge == 0 {
                        continue;
                    }
                    let theoretical_mz = peptide_mz(modified_mass, charge);

                    if within_tolerance(observed_mz, theoretical_mz, precursor_tolerance) {
                        let combined_mods = build_combined_mods(
                            fixed_mods,
                            variable_mods,
                            &combo.sites,
                            &peptide.sequence,
                        );

                        let b_ions = generate_b_ions(&peptide.sequence, &combined_mods);
                        let y_ions = generate_y_ions(&peptide.sequence, &combined_mods);

                        let total_theoretical = (b_ions.len() + y_ions.len()) as u32;
                        if total_theoretical == 0 {
                            continue;
                        }

                        let all_ions: Vec<f64> = b_ions.into_iter().chain(y_ions).collect();
                        let matched =
                            count_matched_ions(&all_ions, &spectrum.mz_array, fragment_tolerance);

                        let score = matched as f64 / total_theoretical as f64;
                        let delta_ppm = calc_delta_ppm(observed_mz, theoretical_mz);

                        if !score.is_finite() || !delta_ppm.is_finite() {
                            continue;
                        }

                        let is_better = match &best_match {
                            None => true,
                            Some(prev) => score > prev.score,
                        };

                        if is_better {
                            let applied = combo
                                .sites
                                .iter()
                                .map(|s| (variable_mods[s.mod_index].clone(), s.residue_pos))
                                .collect();

                            best_match = Some(PeptideMatch {
                                peptide: peptide.clone(),
                                charge,
                                observed_mz,
                                theoretical_mz,
                                delta_mass_ppm: delta_ppm,
                                score,
                                matched_ions: matched,
                                total_ions: total_theoretical,
                                applied_variable_mods: applied,
                            });
                        }
                    }
                }
            }
        }

        if let Some(m) = best_match {
            results.push(m);
        }
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::digest::peptide_mass;
    use protein_copilot_core::spectrum::{MsLevel, PrecursorInfo, Spectrum};

    fn make_spectrum(precursor_mz: f64, charge: Option<i32>, peaks_mz: Vec<f64>) -> Spectrum {
        let intensities = vec![1000.0; peaks_mz.len()];
        Spectrum::new(
            1,
            MsLevel::MS2,
            100.0,
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

    fn make_peptide(sequence: &str, accession: &str) -> DigestedPeptide {
        DigestedPeptide {
            sequence: sequence.to_string(),
            protein_accession: accession.to_string(),
            neutral_mass: peptide_mass(sequence).expect("test uses standard residues"),
            is_protein_nterm: false,
            is_protein_cterm: false,
        }
    }

    fn default_tolerance() -> MassTolerance {
        MassTolerance {
            value: 10.0,
            unit: ToleranceUnit::Ppm,
        }
    }

    fn fragment_tolerance() -> MassTolerance {
        MassTolerance {
            value: 0.02,
            unit: ToleranceUnit::Da,
        }
    }

    #[test]
    fn within_tolerance_ppm() {
        let tol = MassTolerance {
            value: 10.0,
            unit: ToleranceUnit::Ppm,
        };
        assert!(within_tolerance(500.001, 500.0, &tol)); // 2 ppm
        assert!(!within_tolerance(500.01, 500.0, &tol)); // 20 ppm
    }

    #[test]
    fn within_tolerance_da() {
        let tol = MassTolerance {
            value: 0.5,
            unit: ToleranceUnit::Da,
        };
        assert!(within_tolerance(500.3, 500.0, &tol));
        assert!(!within_tolerance(500.6, 500.0, &tol));
    }

    #[test]
    fn b_ion_generation() {
        let ions = generate_b_ions("GG", &[]);
        // b1 for "GG" = G mass + proton = 57.021464 + 1.007276
        assert_eq!(ions.len(), 1);
        assert!((ions[0] - 58.02874).abs() < 0.001);
    }

    #[test]
    fn y_ion_generation() {
        let ions = generate_y_ions("GG", &[]);
        // y1 for "GG" = G mass + water + proton
        assert_eq!(ions.len(), 1);
        assert!((ions[0] - (57.021464 + WATER_MASS + PROTON_MASS)).abs() < 0.001);
    }

    #[test]
    fn match_correct_peptide() {
        let peptide = make_peptide("PEPTIDER", "P001");
        let mz_z2 = peptide_mz(peptide.neutral_mass, 2);

        // Generate b/y ions to use as "experimental" peaks
        let b = generate_b_ions("PEPTIDER", &[]);
        let y = generate_y_ions("PEPTIDER", &[]);
        let mut peaks: Vec<f64> = b.iter().chain(y.iter()).copied().collect();
        peaks.sort_by(|a, b| a.total_cmp(b));

        let spectrum = make_spectrum(mz_z2, Some(2), peaks);
        let result = match_spectrum(
            &spectrum,
            &[peptide],
            &default_tolerance(),
            &fragment_tolerance(),
            &[],
            &[],
            0,
        );

        assert!(result.is_some());
        let m = result.unwrap();
        assert_eq!(m.peptide.sequence, "PEPTIDER");
        assert_eq!(m.charge, 2);
        assert!(m.score > 0.5, "score should be high: {}", m.score);
        assert!(m.delta_mass_ppm.abs() < 1.0);
    }

    #[test]
    fn no_match_outside_tolerance() {
        let peptide = make_peptide("PEPTIDER", "P001");
        // Use a precursor m/z far from any theoretical value
        let spectrum = make_spectrum(999.999, Some(2), vec![100.0, 200.0, 300.0]);
        let result = match_spectrum(
            &spectrum,
            &[peptide],
            &default_tolerance(),
            &fragment_tolerance(),
            &[],
            &[],
            0,
        );

        assert!(result.is_none());
    }

    #[test]
    fn fixed_mods_adjust_mass() {
        let peptide = make_peptide("PEPTIDCK", "P001");
        let cam = Modification {
            name: "Carbamidomethyl".to_string(),
            mass_delta: 57.021464,
            residues: vec!['C'],
            position: protein_copilot_core::search_params::ModPosition::Anywhere,
        };
        let modified_mass =
            peptide.neutral_mass + apply_fixed_mods("PEPTIDCK", std::slice::from_ref(&cam));
        assert!((modified_mass - peptide.neutral_mass - 57.021464).abs() < 0.001);
    }

    #[test]
    fn unknown_charge_tries_multiple() {
        let peptide = make_peptide("PEPTIDER", "P001");
        let mz_z3 = peptide_mz(peptide.neutral_mass, 3);

        let b = generate_b_ions("PEPTIDER", &[]);
        let mut peaks: Vec<f64> = b;
        peaks.sort_by(|a, b| a.total_cmp(b));

        // No charge specified — should try 2, 3, 1, 4
        let spectrum = make_spectrum(mz_z3, None, peaks);
        let result = match_spectrum(
            &spectrum,
            &[peptide],
            &default_tolerance(),
            &fragment_tolerance(),
            &[],
            &[],
            0,
        );

        assert!(result.is_some());
        assert_eq!(result.unwrap().charge, 3);
    }

    #[test]
    fn best_score_wins() {
        let pep1 = make_peptide("PEPTIDER", "P001");
        let pep2 = make_peptide("AVCDEFGK", "P002");
        let mz_z2 = peptide_mz(pep1.neutral_mass, 2);

        // Peaks matching PEPTIDER's b/y ions
        let b = generate_b_ions("PEPTIDER", &[]);
        let y = generate_y_ions("PEPTIDER", &[]);
        let mut peaks: Vec<f64> = b.iter().chain(y.iter()).copied().collect();
        peaks.sort_by(|a, b| a.total_cmp(b));

        let spectrum = make_spectrum(mz_z2, Some(2), peaks);
        let result = match_spectrum(
            &spectrum,
            &[pep1.clone(), pep2],
            &MassTolerance {
                value: 500.0,
                unit: ToleranceUnit::Da,
            }, // wide tolerance to let both match precursor
            &fragment_tolerance(),
            &[],
            &[],
            0,
        );

        assert!(result.is_some());
        assert_eq!(result.unwrap().peptide.sequence, "PEPTIDER");
    }

    #[test]
    fn test_match_spectrum_multiple_precursors() {
        let spectrum = Spectrum::new(
            1,
            MsLevel::MS2,
            100.0,
            vec![
                PrecursorInfo {
                    mz: 500.0,
                    charge: Some(2),
                    intensity: None,
                    isolation_window: None,
                    source_scan: None,
                },
                PrecursorInfo {
                    mz: 600.0,
                    charge: Some(2),
                    intensity: None,
                    isolation_window: None,
                    source_scan: None,
                },
            ],
            vec![175.119, 262.151, 276.134, 363.166],
            vec![100.0, 80.0, 90.0, 70.0],
        )
        .unwrap();

        let peptides = vec![DigestedPeptide {
            sequence: "PEPTIDE".to_string(),
            protein_accession: "test".to_string(),
            neutral_mass: 1197.985448, // (600.0 - 1.007276) * 2
            is_protein_nterm: false,
            is_protein_cterm: false,
        }];

        let tol = MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        };
        let frag_tol = MassTolerance {
            value: 0.02,
            unit: ToleranceUnit::Da,
        };

        let results = match_spectrum_all(&spectrum, &peptides, &tol, &frag_tol, &[], &[], 0);
        assert!(!results.is_empty());
    }

    #[test]
    fn test_match_spectrum_all_empty_precursors() {
        let spectrum = Spectrum::new(
            1,
            MsLevel::MS2,
            100.0,
            vec![],
            vec![100.0, 200.0],
            vec![50.0, 50.0],
        )
        .unwrap();

        let results = match_spectrum_all(
            &spectrum,
            &[],
            &MassTolerance {
                value: 10.0,
                unit: ToleranceUnit::Ppm,
            },
            &MassTolerance {
                value: 0.02,
                unit: ToleranceUnit::Da,
            },
            &[],
            &[],
            0,
        );
        assert!(results.is_empty());
    }

    #[test]
    fn b_ions_include_fixed_modification() {
        use protein_copilot_core::search_params::ModPosition;
        // Peptide "ACK": b2 covers A+C.
        // Without CAM: b2 = M(A) + M(C) + proton = 71.03711 + 103.00919 + 1.00728 = 175.05358
        // With CAM(C): b2 = 71.03711 + 103.00919 + 57.02146 + 1.00728 = 232.07504
        let cam = Modification {
            name: "Carbamidomethyl".to_string(),
            mass_delta: 57.021464,
            residues: vec!['C'],
            position: ModPosition::Anywhere,
        };

        let b_no_mod = generate_b_ions("ACK", &[]);
        let b_with_mod = generate_b_ions("ACK", &[cam]);

        assert_eq!(b_no_mod.len(), 2); // b1, b2
        assert_eq!(b_with_mod.len(), 2);

        // b1 (A only) should be identical — no C residue
        assert!(
            (b_no_mod[0] - b_with_mod[0]).abs() < 0.001,
            "b1 should be unchanged: no_mod={:.4} vs with_mod={:.4}",
            b_no_mod[0],
            b_with_mod[0]
        );

        // b2 (A+C) should differ by CAM mass
        let diff = b_with_mod[1] - b_no_mod[1];
        assert!(
            (diff - 57.021464).abs() < 0.001,
            "b2 should be shifted by 57.02 Da: diff={diff:.4}"
        );
    }

    #[test]
    fn y_ions_include_fixed_modification() {
        use protein_copilot_core::search_params::ModPosition;
        // Peptide "ACK": y1 covers K, y2 covers C+K
        // y1 has no C → unchanged
        // y2 has C → shifted by +57.02
        let cam = Modification {
            name: "Carbamidomethyl".to_string(),
            mass_delta: 57.021464,
            residues: vec!['C'],
            position: ModPosition::Anywhere,
        };

        let y_no_mod = generate_y_ions("ACK", &[]);
        let y_with_mod = generate_y_ions("ACK", &[cam]);

        assert_eq!(y_no_mod.len(), 2);
        assert_eq!(y_with_mod.len(), 2);

        // y1 (K only) — unchanged
        assert!(
            (y_no_mod[0] - y_with_mod[0]).abs() < 0.001,
            "y1 should be unchanged: no_mod={:.4} vs with_mod={:.4}",
            y_no_mod[0],
            y_with_mod[0]
        );

        // y2 (C+K) — shifted by CAM
        let diff = y_with_mod[1] - y_no_mod[1];
        assert!(
            (diff - 57.021464).abs() < 0.001,
            "y2 should be shifted by 57.02 Da: diff={diff:.4}"
        );
    }

    #[test]
    fn modified_fragment_scoring_works() {
        use protein_copilot_core::search_params::ModPosition;
        // Create a peptide with C and generate peaks using modified fragment ions
        let peptide = make_peptide("PEPTIDCK", "P001");
        let cam = Modification {
            name: "Carbamidomethyl".to_string(),
            mass_delta: 57.021464,
            residues: vec!['C'],
            position: ModPosition::Anywhere,
        };
        let modified_mass = peptide.neutral_mass + 57.021464; // one C
        let mz_z2 = peptide_mz(modified_mass, 2);

        // Use modified fragment ions as experimental peaks
        let b = generate_b_ions("PEPTIDCK", std::slice::from_ref(&cam));
        let y = generate_y_ions("PEPTIDCK", std::slice::from_ref(&cam));
        let mut peaks: Vec<f64> = b.iter().chain(y.iter()).copied().collect();
        peaks.sort_by(|a, b| a.total_cmp(b));

        let spectrum = make_spectrum(mz_z2, Some(2), peaks);
        let result = match_spectrum(
            &spectrum,
            &[peptide],
            &default_tolerance(),
            &fragment_tolerance(),
            &[cam],
            &[],
            0,
        );

        assert!(result.is_some());
        let m = result.unwrap();
        // Score should be perfect (1.0) since peaks match exactly
        assert!(
            m.score > 0.99,
            "score should be ~1.0 with correct mod handling: {}",
            m.score
        );
    }
}
