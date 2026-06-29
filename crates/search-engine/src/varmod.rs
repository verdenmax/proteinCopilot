//! Variable modification enumeration.
//!
//! Discovers applicable modification sites on a peptide and generates
//! all valid combinations respecting `max_variable_modifications`.

use protein_copilot_core::search_params::{ModPosition, Modification};

/// Sentinel residue position for a peptide N-terminal (or global) modification.
pub const NTERM_POS: usize = usize::MAX;
/// Sentinel residue position for a peptide C-terminal modification.
///
/// Distinct from [`NTERM_POS`] so an N-terminal and a C-terminal mod are treated
/// as separate chemical sites and can co-occur on the same peptide.
pub const CTERM_POS: usize = usize::MAX - 1;

/// Returns true if `pos` is a terminal sentinel rather than a 0-based residue index.
pub fn is_terminal_pos(pos: usize) -> bool {
    pos == NTERM_POS || pos == CTERM_POS
}

/// A single modification applied at a specific site.
#[derive(Debug, Clone, PartialEq)]
pub struct ModSite {
    /// Index of the modification in the input `variable_mods` slice.
    pub mod_index: usize,
    /// 0-based residue position in the peptide, or [`NTERM_POS`]/[`CTERM_POS`] for terminal.
    pub residue_pos: usize,
}

/// A complete set of variable modifications to apply to a peptide.
#[derive(Debug, Clone)]
pub struct ModCombination {
    /// Which mods are applied and where.
    pub sites: Vec<ModSite>,
    /// Total mass delta from all applied variable mods.
    pub mass_delta: f64,
}

/// Finds all positions in a peptide where a variable modification can apply.
///
/// Respects `ModPosition` filtering:
/// - `Anywhere`: every residue matching `mod.residues`
/// - `AnyNTerm`: position 0 only (if residue matches or residues is empty)
/// - `AnyCTerm`: last position only
/// - `ProteinNTerm`: position 0 only if `is_protein_nterm`
/// - `ProteinCTerm`: last position only if `is_protein_cterm`
///
/// Returns: `Vec<(mod_index, Vec<residue_positions>)>` where [`NTERM_POS`]
/// / [`CTERM_POS`] represent the corresponding peptide terminus.
pub fn find_applicable_sites(
    sequence: &str,
    variable_mods: &[Modification],
    is_protein_nterm: bool,
    is_protein_cterm: bool,
) -> Vec<(usize, Vec<usize>)> {
    let chars: Vec<char> = sequence.chars().collect();
    let last_pos = chars.len().saturating_sub(1);
    let mut result = Vec::new();

    for (mod_idx, m) in variable_mods.iter().enumerate() {
        let mut positions = Vec::new();

        match m.position {
            ModPosition::Anywhere => {
                if m.residues.is_empty() {
                    // Global variable mod — applies once as terminal
                    positions.push(NTERM_POS);
                } else {
                    for (i, &ch) in chars.iter().enumerate() {
                        if m.residues.contains(&ch) {
                            positions.push(i);
                        }
                    }
                }
            }
            ModPosition::AnyNTerm => {
                if m.residues.is_empty() || (!chars.is_empty() && m.residues.contains(&chars[0])) {
                    positions.push(NTERM_POS);
                }
            }
            ModPosition::AnyCTerm => {
                if m.residues.is_empty()
                    || (!chars.is_empty() && m.residues.contains(&chars[last_pos]))
                {
                    positions.push(CTERM_POS);
                }
            }
            ModPosition::ProteinNTerm => {
                if is_protein_nterm
                    && (m.residues.is_empty()
                        || (!chars.is_empty() && m.residues.contains(&chars[0])))
                {
                    positions.push(NTERM_POS);
                }
            }
            ModPosition::ProteinCTerm => {
                if is_protein_cterm
                    && (m.residues.is_empty()
                        || (!chars.is_empty() && m.residues.contains(&chars[last_pos])))
                {
                    positions.push(CTERM_POS);
                }
            }
        }

        if !positions.is_empty() {
            result.push((mod_idx, positions));
        }
    }

    result
}

/// Generates all valid variable modification combinations for a peptide.
///
/// Uses iterative backtracking to enumerate combinations respecting:
/// - `max_mods`: maximum total variable modifications per peptide
/// - Each residue position can only be modified once
/// - Terminal sites ([`NTERM_POS`] / [`CTERM_POS`]) are each a single chemical
///   site, so two mods on the same terminus are mutually exclusive while an
///   N-term and a C-term mod can co-occur
///
/// Always includes the empty combination (no variable mods, mass_delta=0).
pub fn enumerate_combinations(
    variable_mods: &[Modification],
    sites: &[(usize, Vec<usize>)],
    max_mods: u32,
) -> Vec<ModCombination> {
    let mut results = vec![ModCombination {
        sites: Vec::new(),
        mass_delta: 0.0,
    }];

    // Flatten all (mod_index, position, mass_delta) as choosable items
    let mut items: Vec<(usize, usize, f64)> = Vec::new();
    for &(mod_idx, ref positions) in sites {
        let mass = variable_mods[mod_idx].mass_delta;
        for &pos in positions {
            items.push((mod_idx, pos, mass));
        }
    }

    if items.is_empty() || max_mods == 0 {
        return results;
    }

    let max_k = (max_mods as usize).min(items.len());
    let mut stack: Vec<(usize, Vec<ModSite>, f64)> = Vec::new();

    // Seed: start from each item
    for (i, &(mod_idx, pos, mass)) in items.iter().enumerate() {
        stack.push((
            i,
            vec![ModSite {
                mod_index: mod_idx,
                residue_pos: pos,
            }],
            mass,
        ));
    }

    while let Some((last_idx, current_sites, current_mass)) = stack.pop() {
        results.push(ModCombination {
            sites: current_sites.clone(),
            mass_delta: current_mass,
        });

        if current_sites.len() < max_k {
            for (j, &(mod_idx, pos, mass)) in items.iter().enumerate().skip(last_idx + 1) {
                // Skip if this site is already occupied. Residue positions conflict
                // by index; the N- and C-termini (NTERM_POS / CTERM_POS) are distinct
                // single chemical sites, so each conflicts only with its own kind —
                // letting one N-term and one C-term mod co-occur.
                let pos_conflict = current_sites.iter().any(|s| s.residue_pos == pos);
                if pos_conflict {
                    continue;
                }

                let mut new_sites = current_sites.clone();
                new_sites.push(ModSite {
                    mod_index: mod_idx,
                    residue_pos: pos,
                });
                stack.push((j, new_sites, current_mass + mass));
            }
        }
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;

    fn oxidation_m() -> Modification {
        Modification {
            name: "Oxidation".to_string(),
            mass_delta: 15.994915,
            residues: vec!['M'],
            position: ModPosition::Anywhere,
        }
    }

    fn phospho_sty() -> Modification {
        Modification {
            name: "Phospho".to_string(),
            mass_delta: 79.966331,
            residues: vec!['S', 'T', 'Y'],
            position: ModPosition::Anywhere,
        }
    }

    fn acetyl_protein_nterm() -> Modification {
        Modification {
            name: "Acetyl".to_string(),
            mass_delta: 42.010565,
            residues: vec![],
            position: ModPosition::ProteinNTerm,
        }
    }

    // ---- find_applicable_sites tests ----

    #[test]
    fn find_sites_oxidation_on_peptide_with_two_m() {
        let sites = find_applicable_sites("PEPTMDEMK", &[oxidation_m()], false, false);
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].0, 0);
        assert_eq!(sites[0].1, vec![4, 7]);
    }

    #[test]
    fn find_sites_no_matching_residues() {
        let sites = find_applicable_sites("PEPTIDEK", &[oxidation_m()], false, false);
        assert_eq!(sites.len(), 0);
    }

    #[test]
    fn find_sites_phospho_on_sty() {
        // PEPTSDYEK: T@3, S@4, Y@6
        let sites = find_applicable_sites("PEPTSDYEK", &[phospho_sty()], false, false);
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].0, 0);
        assert_eq!(sites[0].1, vec![3, 4, 6]);
    }

    #[test]
    fn find_sites_protein_nterm_on_nterm_peptide() {
        let sites = find_applicable_sites("PEPTIDEK", &[acetyl_protein_nterm()], true, false);
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].1, vec![usize::MAX]);
    }

    #[test]
    fn find_sites_protein_nterm_on_internal_peptide() {
        let sites = find_applicable_sites("PEPTIDEK", &[acetyl_protein_nterm()], false, false);
        assert_eq!(sites.len(), 0);
    }

    #[test]
    fn find_sites_multiple_mods() {
        // PEPTSMEYK: M@5 for Oxidation; T@3, S@4, Y@7 for Phospho
        let mods = vec![oxidation_m(), phospho_sty()];
        let sites = find_applicable_sites("PEPTSMEYK", &mods, false, false);
        assert_eq!(sites.len(), 2);
        assert_eq!(sites[0].1, vec![5]);
        assert_eq!(sites[1].1, vec![3, 4, 7]);
    }

    // ---- enumerate_combinations tests ----

    #[test]
    fn enumerate_empty_when_no_sites() {
        let combos = enumerate_combinations(&[], &[], 3);
        assert_eq!(combos.len(), 1);
        assert_eq!(combos[0].sites.len(), 0);
        assert!(combos[0].mass_delta.abs() < 1e-9);
    }

    #[test]
    fn enumerate_single_mod_single_site() {
        let mods = vec![oxidation_m()];
        let sites = vec![(0_usize, vec![5_usize])];
        let combos = enumerate_combinations(&mods, &sites, 3);
        assert_eq!(combos.len(), 2);
    }

    #[test]
    fn enumerate_single_mod_two_sites() {
        let mods = vec![oxidation_m()];
        let sites = vec![(0, vec![4, 7])];
        let combos = enumerate_combinations(&mods, &sites, 3);
        assert_eq!(combos.len(), 4);
    }

    #[test]
    fn enumerate_respects_max_mods() {
        let mods = vec![oxidation_m()];
        let sites = vec![(0, vec![1, 2, 3])];
        let combos = enumerate_combinations(&mods, &sites, 1);
        assert_eq!(combos.len(), 4);
    }

    #[test]
    fn enumerate_two_mod_types() {
        let mods = vec![oxidation_m(), phospho_sty()];
        let sites = vec![(0, vec![5]), (1, vec![4])];
        let combos = enumerate_combinations(&mods, &sites, 3);
        assert_eq!(combos.len(), 4);
    }

    #[test]
    fn enumerate_mass_delta_correct() {
        let mods = vec![oxidation_m()];
        let sites = vec![(0, vec![4, 7])];
        let combos = enumerate_combinations(&mods, &sites, 3);
        let double = combos.iter().find(|c| c.sites.len() == 2).unwrap();
        let expected = 15.994915 * 2.0;
        assert!((double.mass_delta - expected).abs() < 1e-6);
    }

    #[test]
    fn enumerate_terminal_mod() {
        let mods = vec![acetyl_protein_nterm()];
        let sites = vec![(0, vec![usize::MAX])];
        let combos = enumerate_combinations(&mods, &sites, 3);
        assert_eq!(combos.len(), 2);
        let acetyl = combos.iter().find(|c| c.sites.len() == 1).unwrap();
        assert!((acetyl.mass_delta - 42.010565).abs() < 1e-6);
    }

    #[test]
    fn enumerate_two_terminal_mods_mutually_exclusive() {
        // Two different N-term mods should NOT both apply — a terminus is one site
        let formyl = Modification {
            name: "Formyl".to_string(),
            mass_delta: 27.994915,
            residues: vec![],
            position: ModPosition::AnyNTerm,
        };
        let acetyl_nterm = Modification {
            name: "Acetyl".to_string(),
            mass_delta: 42.010565,
            residues: vec![],
            position: ModPosition::AnyNTerm,
        };
        let mods = vec![formyl, acetyl_nterm];
        let sites = find_applicable_sites("PEPTIDEK", &mods, false, false);
        // Both mods apply to AnyNTerm → both have usize::MAX positions
        assert_eq!(sites.len(), 2);

        let combos = enumerate_combinations(&mods, &sites, 3);
        // Expected: empty + Formyl alone + Acetyl alone = 3 (NOT 4 with both)
        assert_eq!(combos.len(), 3);
        // No combination should have 2 sites (both on same terminus)
        assert!(
            combos.iter().all(|c| c.sites.len() <= 1),
            "two terminal mods must not co-occur on the same terminus"
        );
    }

    #[test]
    fn enumerate_nterm_and_cterm_mods_coexist() {
        // Distinct termini are distinct chemical sites: an N-term and a C-term
        // variable mod must be enumerable together (doubly-modified form).
        let acetyl_nterm = Modification {
            name: "Acetyl".to_string(),
            mass_delta: 42.010565,
            residues: vec![],
            position: ModPosition::AnyNTerm,
        };
        let amide_cterm = Modification {
            name: "Amidation".to_string(),
            mass_delta: -0.984016,
            residues: vec![],
            position: ModPosition::AnyCTerm,
        };
        let mods = vec![acetyl_nterm, amide_cterm];
        let sites = find_applicable_sites("PEPTIDEK", &mods, false, false);
        assert_eq!(sites.len(), 2);

        let combos = enumerate_combinations(&mods, &sites, 3);
        // empty + N alone + C alone + both = 4
        assert_eq!(combos.len(), 4, "N-term and C-term mods must be combinable");
        let both = combos.iter().find(|c| c.sites.len() == 2);
        assert!(
            both.is_some(),
            "doubly-modified (N-term + C-term) form must exist"
        );
        let both = both.unwrap();
        assert!((both.mass_delta - (42.010565 - 0.984016)).abs() < 1e-6);
    }
}
