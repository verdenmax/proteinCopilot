//! Parsimony (greedy set cover) algorithm for protein inference.
//!
//! Given a peptide-to-protein mapping, finds the minimal set of proteins
//! that explains all observed peptides. Groups indistinguishable proteins
//! (those with identical peptide sets) into [`ProteinGroup`]s.

use std::collections::{BTreeSet, HashMap, HashSet};

use protein_copilot_core::protein_group::ProteinGroup;
use protein_copilot_core::util::is_decoy_accession;
use tracing::{debug, info};

use crate::error::InferenceError;
use crate::mapper::PeptideProteinMap;

/// An intermediate group of indistinguishable proteins (identical peptide sets).
#[derive(Debug)]
struct IndistinguishableGroup {
    /// Sorted member accessions (alphabetical for determinism).
    members: Vec<String>,
    /// The peptide set shared by all members.
    peptides: HashSet<String>,
    /// Aggregated score: best peptide score among all peptides.
    score: f64,
}

/// Run the parsimony (greedy set cover) algorithm for protein inference.
///
/// # Steps
/// 1. Group indistinguishable proteins (identical peptide sets).
/// 2. Remove subset proteins (strict subsets are subsumed into the superset group).
/// 3. Greedy set cover: iteratively pick the group covering the most unexplained peptides.
/// 4. Classify peptides as unique (only in one selected group) vs shared.
///
/// # Errors
/// Returns [`InferenceError::NoPsms`] if the input map has no proteins.
pub fn run_parsimony(map: &PeptideProteinMap) -> Result<Vec<ProteinGroup>, InferenceError> {
    if map.protein_to_peptides.is_empty() {
        return Err(InferenceError::NoPsms);
    }

    // Step 1: Group indistinguishable proteins (identical peptide sets).
    let mut indistinguishable = group_indistinguishable(map);
    info!(
        groups = indistinguishable.len(),
        "grouped indistinguishable proteins"
    );

    // Step 2: Remove subset proteins — subsume strict subsets into supersets.
    remove_subsets(&mut indistinguishable);
    info!(groups = indistinguishable.len(), "after subset removal");

    // Step 3: Greedy set cover.
    let selected = greedy_set_cover(&indistinguishable);
    info!(
        selected = selected.len(),
        "greedy set cover selected groups"
    );

    // Step 4: Build ProteinGroup output with unique/shared classification.
    let mut result = build_protein_groups(&selected, map);

    // Sort by score descending, then by leader accession for determinism.
    result.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.leader_accession.cmp(&b.leader_accession))
    });

    debug!(groups = result.len(), "parsimony complete");
    Ok(result)
}

/// Group proteins with identical peptide sets into [`IndistinguishableGroup`]s.
///
/// Uses a sorted `BTreeSet` of peptides as the hash key for grouping.
fn group_indistinguishable(map: &PeptideProteinMap) -> Vec<IndistinguishableGroup> {
    let mut key_to_proteins: HashMap<BTreeSet<String>, Vec<String>> = HashMap::new();

    for (protein, peptides) in &map.protein_to_peptides {
        let key: BTreeSet<String> = peptides.iter().cloned().collect();
        key_to_proteins
            .entry(key)
            .or_default()
            .push(protein.clone());
    }

    key_to_proteins
        .into_iter()
        .map(|(pep_key, mut members)| {
            // Sort members alphabetically for deterministic leader selection.
            members.sort();
            let peptides: HashSet<String> = pep_key.into_iter().collect();
            let score = best_peptide_score(&peptides, &map.peptide_best_score);
            IndistinguishableGroup {
                members,
                peptides,
                score,
            }
        })
        .collect()
}

/// Remove subset proteins: if group A's peptides ⊂ group B's peptides,
/// merge A's members into B and drop A.
fn remove_subsets(groups: &mut Vec<IndistinguishableGroup>) {
    // Sort by peptide count descending so larger sets come first.
    groups.sort_by(|a, b| b.peptides.len().cmp(&a.peptides.len()));

    let n = groups.len();
    let mut subsumed = vec![false; n];
    // For each subsumed group, track which superset absorbs it.
    let mut merge_into: Vec<Option<usize>> = vec![None; n];

    for i in 0..n {
        if subsumed[i] {
            continue;
        }
        for j in (i + 1)..n {
            if subsumed[j] {
                continue;
            }
            // Check if j is a strict subset of i.
            if groups[j].peptides.len() < groups[i].peptides.len()
                && groups[j].peptides.is_subset(&groups[i].peptides)
            {
                subsumed[j] = true;
                merge_into[j] = Some(i);
            }
        }
    }

    // Merge subsumed members into their supersets.
    // Collect merge operations first to avoid borrow issues.
    let merges: Vec<(usize, Vec<String>)> = merge_into
        .iter()
        .enumerate()
        .filter_map(|(j, target)| target.map(|i| (i, groups[j].members.clone())))
        .collect();

    for (target_idx, members) in merges {
        groups[target_idx].members.extend(members);
    }

    // Re-sort merged members for determinism.
    for g in groups.iter_mut() {
        g.members.sort();
    }

    // Remove subsumed groups (iterate in reverse to preserve indices).
    let mut idx = n;
    while idx > 0 {
        idx -= 1;
        if subsumed[idx] {
            groups.remove(idx);
        }
    }
}

/// Greedy set cover: iteratively pick the group that explains the most
/// uncovered peptides until all peptides are covered.
fn greedy_set_cover(groups: &[IndistinguishableGroup]) -> Vec<&IndistinguishableGroup> {
    let mut uncovered: HashSet<&str> = groups
        .iter()
        .flat_map(|g| g.peptides.iter().map(String::as_str))
        .collect();

    let mut selected: Vec<&IndistinguishableGroup> = Vec::new();
    let mut used = vec![false; groups.len()];

    while !uncovered.is_empty() {
        // Find the group covering the most uncovered peptides.
        // Ties broken by: (1) highest score, (2) first member alphabetically.
        let best_idx = groups
            .iter()
            .enumerate()
            .filter(|(i, _)| !used[*i])
            .max_by(|(_, a), (_, b)| {
                let a_cover = a
                    .peptides
                    .iter()
                    .filter(|p| uncovered.contains(p.as_str()))
                    .count();
                let b_cover = b
                    .peptides
                    .iter()
                    .filter(|p| uncovered.contains(p.as_str()))
                    .count();
                a_cover
                    .cmp(&b_cover)
                    .then_with(|| {
                        a.score
                            .partial_cmp(&b.score)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .then_with(|| b.members[0].cmp(&a.members[0])) // alphabetically first = "less" = better
            })
            .map(|(i, _)| i);

        match best_idx {
            Some(idx) => {
                used[idx] = true;
                for pep in &groups[idx].peptides {
                    uncovered.remove(pep.as_str());
                }
                selected.push(&groups[idx]);
            }
            None => break, // No more groups available (shouldn't happen with valid input).
        }
    }

    selected
}

/// Build the final [`ProteinGroup`] output from selected groups.
///
/// Classifies each peptide as unique (appears in exactly one selected group)
/// or shared (appears in multiple selected groups).
fn build_protein_groups(
    selected: &[&IndistinguishableGroup],
    _map: &PeptideProteinMap,
) -> Vec<ProteinGroup> {
    // Count how many selected groups each peptide appears in.
    let mut peptide_group_count: HashMap<&str, usize> = HashMap::new();
    for group in selected {
        for pep in &group.peptides {
            *peptide_group_count.entry(pep.as_str()).or_insert(0) += 1;
        }
    }

    selected
        .iter()
        .map(|group| {
            let leader_accession = group.members[0].clone(); // alphabetically first

            let mut peptides: Vec<String> = group.peptides.iter().cloned().collect();
            peptides.sort();

            let mut unique_peptides: Vec<String> = group
                .peptides
                .iter()
                .filter(|p| peptide_group_count.get(p.as_str()).copied().unwrap_or(0) == 1)
                .cloned()
                .collect();
            unique_peptides.sort();

            let is_decoy = group
                .members
                .iter()
                .all(|acc| is_decoy_accession(acc));

            ProteinGroup {
                leader_accession,
                leader_description: String::new(),
                member_accessions: group.members.clone(),
                peptides,
                unique_peptides,
                razor_peptides: Vec::new(),
                score: group.score,
                q_value: None,
                coverage: None,
                is_decoy,
            }
        })
        .collect()
}

/// Get the best peptide score from a set of peptides.
fn best_peptide_score(peptides: &HashSet<String>, scores: &HashMap<String, f64>) -> f64 {
    peptides
        .iter()
        .filter_map(|p| scores.get(p))
        .copied()
        .fold(f64::NEG_INFINITY, f64::max)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: build a PeptideProteinMap from a list of (peptide, proteins, score) tuples.
    fn build_map(entries: &[(&str, &[&str], f64)]) -> PeptideProteinMap {
        let mut peptide_to_proteins: HashMap<String, HashSet<String>> = HashMap::new();
        let mut protein_to_peptides: HashMap<String, HashSet<String>> = HashMap::new();
        let mut peptide_best_score: HashMap<String, f64> = HashMap::new();

        for &(pep, prots, score) in entries {
            let pep = pep.to_string();
            for &prot in prots {
                peptide_to_proteins
                    .entry(pep.clone())
                    .or_default()
                    .insert(prot.to_string());
                protein_to_peptides
                    .entry(prot.to_string())
                    .or_default()
                    .insert(pep.clone());
            }
            peptide_best_score
                .entry(pep.clone())
                .and_modify(|s| {
                    if score > *s {
                        *s = score;
                    }
                })
                .or_insert(score);
        }

        let peptide_is_decoy: HashMap<String, bool> = peptide_to_proteins
            .iter()
            .map(|(pep, prots)| {
                let all_decoy = prots.iter().all(|acc| is_decoy_accession(acc));
                (pep.clone(), all_decoy)
            })
            .collect();

        PeptideProteinMap {
            peptide_to_proteins,
            protein_to_peptides,
            peptide_best_score,
            peptide_is_decoy,
        }
    }

    #[test]
    fn test_single_protein_unique_peptides() {
        let map = build_map(&[
            ("PEP1", &["PROT_A"], 10.0),
            ("PEP2", &["PROT_A"], 8.0),
            ("PEP3", &["PROT_A"], 6.0),
        ]);
        let groups = run_parsimony(&map).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].leader_accession, "PROT_A");
        assert_eq!(groups[0].member_accessions, vec!["PROT_A"]);
        assert_eq!(groups[0].peptides.len(), 3);
        assert_eq!(
            groups[0].unique_peptides.len(),
            3,
            "all peptides should be unique"
        );
        assert_eq!(groups[0].score, 10.0);
        assert!(!groups[0].is_decoy);
    }

    #[test]
    fn test_two_indistinguishable_proteins() {
        // PROT_A and PROT_B share identical peptide sets → merged into one group.
        let map = build_map(&[
            ("PEP1", &["PROT_A", "PROT_B"], 10.0),
            ("PEP2", &["PROT_A", "PROT_B"], 8.0),
        ]);
        let groups = run_parsimony(&map).unwrap();
        assert_eq!(groups.len(), 1);
        // Leader is alphabetically first.
        assert_eq!(groups[0].leader_accession, "PROT_A");
        assert_eq!(groups[0].member_accessions, vec!["PROT_A", "PROT_B"]);
        assert_eq!(groups[0].peptides.len(), 2);
        assert!(!groups[0].is_decoy);
    }

    #[test]
    fn test_subset_protein() {
        // PROT_A has {PEP1, PEP2, PEP3}, PROT_B has {PEP1, PEP2} ⊂ PROT_A.
        // PROT_B should be subsumed into PROT_A's group.
        let map = build_map(&[
            ("PEP1", &["PROT_A", "PROT_B"], 10.0),
            ("PEP2", &["PROT_A", "PROT_B"], 8.0),
            ("PEP3", &["PROT_A"], 6.0),
        ]);
        let groups = run_parsimony(&map).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].leader_accession, "PROT_A");
        assert!(
            groups[0].member_accessions.contains(&"PROT_B".to_string()),
            "subset protein should be a member"
        );
        assert_eq!(groups[0].peptides.len(), 3);
    }

    #[test]
    fn test_greedy_cover_overlapping() {
        // PROT_A: {PEP1, PEP2, PEP3}  (covers 3)
        // PROT_B: {PEP3, PEP4}         (overlaps on PEP3)
        // PROT_C: {PEP5}               (unique peptide)
        // Greedy should pick PROT_A first (3 peptides), then PROT_B (covers PEP4),
        // then PROT_C (covers PEP5).
        let map = build_map(&[
            ("PEP1", &["PROT_A"], 10.0),
            ("PEP2", &["PROT_A"], 9.0),
            ("PEP3", &["PROT_A", "PROT_B"], 8.0),
            ("PEP4", &["PROT_B"], 7.0),
            ("PEP5", &["PROT_C"], 6.0),
        ]);
        let groups = run_parsimony(&map).unwrap();
        assert_eq!(groups.len(), 3);

        // Find the groups by leader.
        let group_a = groups
            .iter()
            .find(|g| g.leader_accession == "PROT_A")
            .unwrap();
        let group_b = groups
            .iter()
            .find(|g| g.leader_accession == "PROT_B")
            .unwrap();
        let group_c = groups
            .iter()
            .find(|g| g.leader_accession == "PROT_C")
            .unwrap();

        // PROT_A: PEP1, PEP2 unique; PEP3 shared with PROT_B.
        assert!(group_a.unique_peptides.contains(&"PEP1".to_string()));
        assert!(group_a.unique_peptides.contains(&"PEP2".to_string()));
        assert!(
            !group_a.unique_peptides.contains(&"PEP3".to_string()),
            "PEP3 is shared, not unique to PROT_A"
        );
        assert!(group_a.peptides.contains(&"PEP3".to_string()));

        // PROT_B: PEP4 unique; PEP3 shared.
        assert!(group_b.unique_peptides.contains(&"PEP4".to_string()));
        assert!(!group_b.unique_peptides.contains(&"PEP3".to_string()));

        // PROT_C: PEP5 unique.
        assert_eq!(group_c.unique_peptides, vec!["PEP5"]);
    }

    #[test]
    fn test_decoy_group() {
        let map = build_map(&[
            ("PEP1", &["REV_PROT_A"], 5.0),
            ("PEP2", &["REV_PROT_A"], 3.0),
            // Need at least one target protein for the map to be valid in production,
            // but our test helper doesn't enforce that, so add one.
            ("PEP3", &["PROT_B"], 10.0),
        ]);
        let groups = run_parsimony(&map).unwrap();
        let decoy_group = groups
            .iter()
            .find(|g| g.leader_accession == "REV_PROT_A")
            .unwrap();
        assert!(
            decoy_group.is_decoy,
            "group with all REV_ members should be decoy"
        );

        let target_group = groups
            .iter()
            .find(|g| g.leader_accession == "PROT_B")
            .unwrap();
        assert!(!target_group.is_decoy);
    }

    #[test]
    fn test_mixed_target_decoy_group() {
        // Target and decoy protein share peptides → they become indistinguishable →
        // group is NOT decoy (has target member).
        let map = build_map(&[
            ("PEP1", &["PROT_A", "REV_PROT_B"], 10.0),
            ("PEP2", &["PROT_A", "REV_PROT_B"], 8.0),
        ]);
        let groups = run_parsimony(&map).unwrap();
        assert_eq!(groups.len(), 1);
        assert!(
            !groups[0].is_decoy,
            "group with target+decoy members should NOT be decoy"
        );
    }

    #[test]
    fn test_empty_map_returns_error() {
        let map = PeptideProteinMap {
            peptide_to_proteins: HashMap::new(),
            protein_to_peptides: HashMap::new(),
            peptide_best_score: HashMap::new(),
            peptide_is_decoy: HashMap::new(),
        };
        let result = run_parsimony(&map);
        assert!(
            matches!(result, Err(InferenceError::NoPsms)),
            "empty map should return NoPsms error"
        );
    }

    #[test]
    fn test_score_selection_and_leader() {
        // Two proteins with different peptides but we verify scoring and leader selection.
        // PROT_B has a higher-scoring peptide, but leader should be alphabetically first within group.
        let map = build_map(&[("PEP1", &["PROT_A"], 5.0), ("PEP2", &["PROT_B"], 15.0)]);
        let groups = run_parsimony(&map).unwrap();
        assert_eq!(groups.len(), 2);

        // Sorted by score descending: PROT_B (15) first, PROT_A (5) second.
        assert_eq!(groups[0].leader_accession, "PROT_B");
        assert_eq!(groups[0].score, 15.0);
        assert_eq!(groups[1].leader_accession, "PROT_A");
        assert_eq!(groups[1].score, 5.0);
    }

    #[test]
    fn test_indistinguishable_leader_is_alphabetically_first() {
        // Three indistinguishable proteins: leader should be "ALPHA" (first alphabetically).
        let map = build_map(&[
            ("PEP1", &["GAMMA", "ALPHA", "BETA"], 10.0),
            ("PEP2", &["GAMMA", "ALPHA", "BETA"], 8.0),
        ]);
        let groups = run_parsimony(&map).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].leader_accession, "ALPHA");
        assert_eq!(groups[0].member_accessions, vec!["ALPHA", "BETA", "GAMMA"]);
    }

    #[test]
    fn test_subset_chain() {
        // PROT_A: {PEP1, PEP2, PEP3}
        // PROT_B: {PEP1, PEP2}  — subset of A
        // PROT_C: {PEP1}        — subset of B and A
        // All should merge into one group led by PROT_A.
        let map = build_map(&[
            ("PEP1", &["PROT_A", "PROT_B", "PROT_C"], 10.0),
            ("PEP2", &["PROT_A", "PROT_B"], 8.0),
            ("PEP3", &["PROT_A"], 6.0),
        ]);
        let groups = run_parsimony(&map).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].leader_accession, "PROT_A");
        assert_eq!(groups[0].member_accessions.len(), 3);
        assert!(groups[0].member_accessions.contains(&"PROT_B".to_string()));
        assert!(groups[0].member_accessions.contains(&"PROT_C".to_string()));
    }

    #[test]
    fn test_all_peptides_covered() {
        // Verify that every peptide in the input appears in at least one output group.
        let map = build_map(&[
            ("PEP1", &["PROT_A"], 10.0),
            ("PEP2", &["PROT_A", "PROT_B"], 8.0),
            ("PEP3", &["PROT_B"], 6.0),
            ("PEP4", &["PROT_C"], 4.0),
        ]);
        let groups = run_parsimony(&map).unwrap();

        let covered: HashSet<&str> = groups
            .iter()
            .flat_map(|g| g.peptides.iter().map(String::as_str))
            .collect();

        for pep in map.peptide_to_proteins.keys() {
            assert!(covered.contains(pep.as_str()), "peptide {pep} not covered");
        }
    }
}
