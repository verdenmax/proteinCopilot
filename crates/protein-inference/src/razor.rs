//! Razor peptide assignment.
//!
//! Shared peptides (those mapping to multiple protein groups) are assigned
//! to the group with the most unique peptide evidence ("razor" logic).
//! This maximizes quantitative accuracy by avoiding double-counting.

use std::collections::HashMap;

use protein_copilot_core::protein_group::ProteinGroup;
use tracing::debug;

use crate::mapper::PeptideProteinMap;

/// Assign shared peptides to protein groups using razor logic.
///
/// For each shared peptide (appearing in multiple groups), assigns it to the group
/// with the most unique peptides. Ties broken by score, then leader accession.
///
/// Returns a map: normalized_peptide_sequence → leader_accession of assigned group.
/// Also mutates groups in-place to populate their `razor_peptides` field.
pub fn assign_razor_peptides(
    groups: &mut [ProteinGroup],
    _map: &PeptideProteinMap,
) -> HashMap<String, String> {
    let mut razor_map: HashMap<String, String> = HashMap::new();

    if groups.len() <= 1 {
        debug!(groups = groups.len(), "skipping razor assignment (≤1 group)");
        return razor_map;
    }

    // Step 1: Build a lookup from peptide → list of group indices that contain it.
    let mut peptide_to_group_indices: HashMap<&str, Vec<usize>> = HashMap::new();
    for (idx, group) in groups.iter().enumerate() {
        for pep in &group.peptides {
            peptide_to_group_indices
                .entry(pep.as_str())
                .or_default()
                .push(idx);
        }
    }

    // Step 2: Identify shared peptides (appear in >1 group) that are not unique to any group.
    let shared_peptides: Vec<String> = peptide_to_group_indices
        .iter()
        .filter(|(_, indices)| indices.len() > 1)
        .map(|(&pep, _)| pep.to_string())
        .collect();

    // Build a snapshot of each group's unique peptide count, score, and leader
    // so we can compare without borrowing `groups` mutably.
    let group_keys: Vec<(usize, usize, f64, String)> = groups
        .iter()
        .enumerate()
        .map(|(idx, g)| {
            (
                idx,
                g.unique_peptides.len(),
                g.score,
                g.leader_accession.clone(),
            )
        })
        .collect();

    // Step 3: For each shared peptide, assign it to the best group.
    for pep in &shared_peptides {
        let indices = &peptide_to_group_indices[pep.as_str()];

        // Skip peptides that are already unique to one group.
        if indices.iter().any(|&i| groups[i].unique_peptides.contains(pep)) {
            continue;
        }

        // Pick the best group among those containing this peptide.
        let best_idx = indices
            .iter()
            .copied()
            .max_by(|&a, &b| {
                let (_, a_uniq, a_score, ref a_leader) = group_keys[a];
                let (_, b_uniq, b_score, ref b_leader) = group_keys[b];
                a_uniq
                    .cmp(&b_uniq)
                    .then_with(|| {
                        a_score
                            .partial_cmp(&b_score)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .then_with(|| b_leader.cmp(a_leader)) // alphabetically first wins
            })
            .expect("indices is non-empty for shared peptides");

        razor_map.insert(pep.clone(), groups[best_idx].leader_accession.clone());
    }

    // Step 4: Populate razor_peptides on each group from the razor_map.
    for group in groups.iter_mut() {
        group.razor_peptides.clear();
    }
    for (pep, leader) in &razor_map {
        if let Some(group) = groups.iter_mut().find(|g| g.leader_accession == *leader) {
            group.razor_peptides.push(pep.clone());
        }
    }
    // Sort for determinism.
    for group in groups.iter_mut() {
        group.razor_peptides.sort();
    }

    debug!(
        razor_count = razor_map.len(),
        "razor peptide assignment complete"
    );
    razor_map
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: create a ProteinGroup with given fields.
    fn make_group(
        leader: &str,
        peptides: &[&str],
        unique_peptides: &[&str],
        score: f64,
    ) -> ProteinGroup {
        ProteinGroup {
            leader_accession: leader.to_string(),
            leader_description: String::new(),
            member_accessions: vec![leader.to_string()],
            peptides: peptides.iter().map(|s| s.to_string()).collect(),
            unique_peptides: unique_peptides.iter().map(|s| s.to_string()).collect(),
            razor_peptides: Vec::new(),
            score,
            q_value: None,
            coverage: None,
            is_decoy: false,
        }
    }

    /// Helper: build a minimal PeptideProteinMap (not used by razor logic
    /// but required by the function signature).
    fn empty_map() -> PeptideProteinMap {
        PeptideProteinMap {
            peptide_to_proteins: HashMap::new(),
            protein_to_peptides: HashMap::new(),
            peptide_best_score: HashMap::new(),
            peptide_is_decoy: HashMap::new(),
        }
    }

    #[test]
    fn test_no_shared_peptides() {
        // All peptides unique to their groups → empty razor map.
        let mut groups = vec![
            make_group("PROT_A", &["PEP1", "PEP2"], &["PEP1", "PEP2"], 10.0),
            make_group("PROT_B", &["PEP3", "PEP4"], &["PEP3", "PEP4"], 8.0),
        ];
        let map = empty_map();
        let razor = assign_razor_peptides(&mut groups, &map);

        assert!(razor.is_empty(), "no shared peptides → empty razor map");
        assert!(groups[0].razor_peptides.is_empty());
        assert!(groups[1].razor_peptides.is_empty());
    }

    #[test]
    fn test_single_shared_peptide() {
        // PEP_S is shared; PROT_A has 2 unique, PROT_B has 1 unique → assign to PROT_A.
        let mut groups = vec![
            make_group(
                "PROT_A",
                &["PEP1", "PEP2", "PEP_S"],
                &["PEP1", "PEP2"],
                10.0,
            ),
            make_group("PROT_B", &["PEP3", "PEP_S"], &["PEP3"], 8.0),
        ];
        let map = empty_map();
        let razor = assign_razor_peptides(&mut groups, &map);

        assert_eq!(razor.len(), 1);
        assert_eq!(razor["PEP_S"], "PROT_A");
        assert_eq!(groups[0].razor_peptides, vec!["PEP_S"]);
        assert!(groups[1].razor_peptides.is_empty());
    }

    #[test]
    fn test_tiebreak_by_score() {
        // Same unique peptide count (1 each), different scores → higher score wins.
        let mut groups = vec![
            make_group("PROT_A", &["PEP1", "PEP_S"], &["PEP1"], 5.0),
            make_group("PROT_B", &["PEP2", "PEP_S"], &["PEP2"], 15.0),
        ];
        let map = empty_map();
        let razor = assign_razor_peptides(&mut groups, &map);

        assert_eq!(razor.len(), 1);
        assert_eq!(razor["PEP_S"], "PROT_B", "higher score should win tie-break");
        assert!(groups[0].razor_peptides.is_empty());
        assert_eq!(groups[1].razor_peptides, vec!["PEP_S"]);
    }

    #[test]
    fn test_tiebreak_by_accession() {
        // Same unique count, same score → alphabetically first leader wins.
        let mut groups = vec![
            make_group("PROT_B", &["PEP1", "PEP_S"], &["PEP1"], 10.0),
            make_group("PROT_A", &["PEP2", "PEP_S"], &["PEP2"], 10.0),
        ];
        let map = empty_map();
        let razor = assign_razor_peptides(&mut groups, &map);

        assert_eq!(razor.len(), 1);
        assert_eq!(
            razor["PEP_S"], "PROT_A",
            "alphabetically first leader should win"
        );
        assert!(groups[0].razor_peptides.is_empty());
        assert_eq!(groups[1].razor_peptides, vec!["PEP_S"]);
    }

    #[test]
    fn test_multiple_shared_peptides() {
        // PROT_A: 3 unique, PROT_B: 1 unique; two shared peptides → both go to PROT_A.
        let mut groups = vec![
            make_group(
                "PROT_A",
                &["PEP1", "PEP2", "PEP3", "PEP_S1", "PEP_S2"],
                &["PEP1", "PEP2", "PEP3"],
                10.0,
            ),
            make_group(
                "PROT_B",
                &["PEP4", "PEP_S1", "PEP_S2"],
                &["PEP4"],
                8.0,
            ),
        ];
        let map = empty_map();
        let razor = assign_razor_peptides(&mut groups, &map);

        assert_eq!(razor.len(), 2);
        assert_eq!(razor["PEP_S1"], "PROT_A");
        assert_eq!(razor["PEP_S2"], "PROT_A");
        assert_eq!(groups[0].razor_peptides, vec!["PEP_S1", "PEP_S2"]);
        assert!(groups[1].razor_peptides.is_empty());
    }

    #[test]
    fn test_single_group() {
        // Only one group → no sharing possible → empty razor map.
        let mut groups = vec![make_group(
            "PROT_A",
            &["PEP1", "PEP2"],
            &["PEP1", "PEP2"],
            10.0,
        )];
        let map = empty_map();
        let razor = assign_razor_peptides(&mut groups, &map);

        assert!(razor.is_empty());
        assert!(groups[0].razor_peptides.is_empty());
    }

    #[test]
    fn test_all_shared_no_unique() {
        // Both groups have 0 unique peptides → tie-break by score, then alphabetical.
        let mut groups = vec![
            make_group("PROT_B", &["PEP_S1", "PEP_S2"], &[], 10.0),
            make_group("PROT_A", &["PEP_S1", "PEP_S2"], &[], 10.0),
        ];
        let map = empty_map();
        let razor = assign_razor_peptides(&mut groups, &map);

        assert_eq!(razor.len(), 2);
        // Same unique count (0), same score (10.0) → alphabetically first leader "PROT_A" wins.
        assert_eq!(razor["PEP_S1"], "PROT_A");
        assert_eq!(razor["PEP_S2"], "PROT_A");
        assert!(groups[0].razor_peptides.is_empty());
        assert_eq!(groups[1].razor_peptides, vec!["PEP_S1", "PEP_S2"]);
    }

    #[test]
    fn test_shared_peptide_assigned_to_different_groups() {
        // Three groups with different shared peptides going to different winners.
        // PROT_A: 2 unique; PROT_B: 3 unique; PROT_C: 1 unique
        // PEP_AB shared between A and B → B wins (3 > 2)
        // PEP_AC shared between A and C → A wins (2 > 1)
        let mut groups = vec![
            make_group(
                "PROT_A",
                &["PEP1", "PEP2", "PEP_AB", "PEP_AC"],
                &["PEP1", "PEP2"],
                10.0,
            ),
            make_group(
                "PROT_B",
                &["PEP3", "PEP4", "PEP5", "PEP_AB"],
                &["PEP3", "PEP4", "PEP5"],
                8.0,
            ),
            make_group("PROT_C", &["PEP6", "PEP_AC"], &["PEP6"], 6.0),
        ];
        let map = empty_map();
        let razor = assign_razor_peptides(&mut groups, &map);

        assert_eq!(razor.len(), 2);
        assert_eq!(razor["PEP_AB"], "PROT_B");
        assert_eq!(razor["PEP_AC"], "PROT_A");
        assert_eq!(groups[0].razor_peptides, vec!["PEP_AC"]);
        assert_eq!(groups[1].razor_peptides, vec!["PEP_AB"]);
        assert!(groups[2].razor_peptides.is_empty());
    }

    #[test]
    fn test_unique_peptides_not_in_razor_map() {
        // Unique peptides must never appear in the razor map.
        let mut groups = vec![
            make_group(
                "PROT_A",
                &["PEP1", "PEP2", "PEP_S"],
                &["PEP1", "PEP2"],
                10.0,
            ),
            make_group("PROT_B", &["PEP3", "PEP_S"], &["PEP3"], 8.0),
        ];
        let map = empty_map();
        let razor = assign_razor_peptides(&mut groups, &map);

        assert!(!razor.contains_key("PEP1"));
        assert!(!razor.contains_key("PEP2"));
        assert!(!razor.contains_key("PEP3"));
        assert!(razor.contains_key("PEP_S"));
    }

    #[test]
    fn test_empty_groups_slice() {
        let mut groups: Vec<ProteinGroup> = Vec::new();
        let map = empty_map();
        let razor = assign_razor_peptides(&mut groups, &map);
        assert!(razor.is_empty());
    }
}
