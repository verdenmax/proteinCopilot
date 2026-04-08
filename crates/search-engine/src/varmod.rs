//! Variable modification enumeration.
//!
//! Discovers applicable modification sites on a peptide and generates
//! all valid combinations respecting `max_variable_modifications`.

use protein_copilot_core::search_params::{ModPosition, Modification};

/// A single modification applied at a specific site.
#[derive(Debug, Clone, PartialEq)]
pub struct ModSite {
    /// Index of the modification in the input `variable_mods` slice.
    pub mod_index: usize,
    /// 0-based residue position in the peptide, or `usize::MAX` for terminal.
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
/// Returns: `Vec<(mod_index, Vec<residue_positions>)>` where `usize::MAX`
/// represents a terminal position.
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
                    positions.push(usize::MAX);
                } else {
                    for (i, &ch) in chars.iter().enumerate() {
                        if m.residues.contains(&ch) {
                            positions.push(i);
                        }
                    }
                }
            }
            ModPosition::AnyNTerm => {
                if m.residues.is_empty()
                    || (!chars.is_empty() && m.residues.contains(&chars[0]))
                {
                    positions.push(usize::MAX);
                }
            }
            ModPosition::AnyCTerm => {
                if m.residues.is_empty()
                    || (!chars.is_empty() && m.residues.contains(&chars[last_pos]))
                {
                    positions.push(usize::MAX);
                }
            }
            ModPosition::ProteinNTerm => {
                if is_protein_nterm
                    && (m.residues.is_empty()
                        || (!chars.is_empty() && m.residues.contains(&chars[0])))
                {
                    positions.push(usize::MAX);
                }
            }
            ModPosition::ProteinCTerm => {
                if is_protein_cterm
                    && (m.residues.is_empty()
                        || (!chars.is_empty() && m.residues.contains(&chars[last_pos])))
                {
                    positions.push(usize::MAX);
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
/// - Terminal sites (`usize::MAX`) are exclusive per modification index
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
            for j in (last_idx + 1)..items.len() {
                let (mod_idx, pos, mass) = items[j];
                // Skip if same residue position already modified
                let pos_conflict = if pos == usize::MAX {
                    current_sites
                        .iter()
                        .any(|s| s.residue_pos == usize::MAX && s.mod_index == mod_idx)
                } else {
                    current_sites.iter().any(|s| s.residue_pos == pos)
                };
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
}
