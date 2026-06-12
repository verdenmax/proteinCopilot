//! Protein-level FDR using the picked-protein (picked target-decoy) approach.
//!
//! Algorithm:
//! 1. Pair each target protein group with its corresponding decoy (REV_ prefix)
//! 2. Compete within pairs: the higher-scoring group wins (ties → target wins)
//! 3. Unpaired targets win automatically; unpaired decoys enter as decoy winners
//! 4. Apply standard target-decoy FDR on the winner list
//! 5. Assign q-values to winning target groups

use std::collections::{HashMap, HashSet};

use crate::calculation::{calculate_fdr, ScoredPsm};
use crate::error::FdrError;
use protein_copilot_core::protein_group::ProteinGroup;

/// Decoy prefix used to pair target ↔ decoy protein groups.
const DECOY_PREFIX: &str = "REV_";

/// Result of protein-level FDR calculation.
#[derive(Debug, Clone)]
pub struct ProteinFdrResult {
    /// Winning target groups with q-values assigned.
    pub groups: Vec<ProteinGroup>,
    /// Number of target groups at 1% FDR.
    pub target_groups_at_1pct: u64,
    /// Total target groups (before FDR filtering).
    pub total_target_groups: u64,
    /// Total decoy groups.
    pub total_decoy_groups: u64,
}

/// Calculate protein-level FDR using picked-protein approach.
///
/// Pairs target/decoy protein groups, competes within pairs, then calculates
/// FDR on winners. Returns updated groups with q-values assigned.
pub fn calculate_protein_fdr(groups: &[ProteinGroup]) -> Result<ProteinFdrResult, FdrError> {
    if groups.is_empty() {
        return Err(FdrError::NoPsms);
    }

    let total_target_groups = groups.iter().filter(|g| !g.is_decoy).count() as u64;
    let total_decoy_groups = groups.iter().filter(|g| g.is_decoy).count() as u64;

    // Index target and decoy groups by their "base" accession for pairing.
    // Target key = leader_accession, Decoy key = leader_accession without REV_ prefix.
    let mut target_map: HashMap<&str, usize> = HashMap::new();
    let mut decoy_map: HashMap<String, usize> = HashMap::new();

    for (i, g) in groups.iter().enumerate() {
        if g.is_decoy {
            if let Some(base) = g.leader_accession.strip_prefix(DECOY_PREFIX) {
                decoy_map.insert(base.to_string(), i);
            } else {
                // Decoy without REV_ prefix — still treat as decoy, will be unpaired
                decoy_map.insert(g.leader_accession.clone(), i);
            }
        } else {
            target_map.insert(&g.leader_accession, i);
        }
    }

    // Compete within pairs and collect winners.
    let mut winners: Vec<(usize, bool)> = Vec::new(); // (index into groups, is_decoy)
    let mut paired_decoy_keys: HashSet<String> = HashSet::new();

    for (&target_acc, &target_idx) in &target_map {
        if let Some(&decoy_idx) = decoy_map.get(target_acc) {
            // Paired: compete. Ties → target wins.
            paired_decoy_keys.insert(target_acc.to_string());
            if groups[target_idx].score >= groups[decoy_idx].score {
                winners.push((target_idx, false));
            } else {
                winners.push((decoy_idx, true));
            }
        } else {
            // Unpaired target: wins automatically
            winners.push((target_idx, false));
        }
    }

    // Unpaired decoys: decoys whose base accession has no matching target
    for (base_acc, &decoy_idx) in &decoy_map {
        if !paired_decoy_keys.contains(base_acc) {
            winners.push((decoy_idx, true));
        }
    }

    let has_decoy_winners = winners.iter().any(|(_, is_decoy)| *is_decoy);
    if !has_decoy_winners {
        // No decoy winners — cannot estimate FDR via TDC.
        // All targets win with q_value = 0.0.
        let mut result_groups: Vec<ProteinGroup> = winners
            .iter()
            .filter(|(_, is_decoy)| !is_decoy)
            .map(|(idx, _)| {
                let mut g = groups[*idx].clone();
                g.q_value = Some(0.0);
                g
            })
            .collect();
        result_groups.sort_by(|a, b| {
            b.score
                .total_cmp(&a.score)
                .then_with(|| a.leader_accession.cmp(&b.leader_accession))
        });

        let target_groups_at_1pct = result_groups.len() as u64;
        return Ok(ProteinFdrResult {
            groups: result_groups,
            target_groups_at_1pct,
            total_target_groups,
            total_decoy_groups,
        });
    }

    // Build ScoredPsm entries for the FDR engine (reuses PSM-level TDC).
    let scored: Vec<ScoredPsm> = winners
        .iter()
        .enumerate()
        .map(|(i, (idx, is_decoy))| ScoredPsm {
            index: i,
            score: groups[*idx].score,
            is_decoy: *is_decoy,
        })
        .collect();

    let fdr_results = calculate_fdr(&scored)?;

    // Map q-values back to winning target groups.
    let mut q_map: HashMap<usize, f64> = HashMap::new();
    for (winner_i, q) in &fdr_results {
        q_map.insert(*winner_i, *q);
    }

    let mut result_groups: Vec<ProteinGroup> = Vec::new();
    for (winner_i, (group_idx, is_decoy)) in winners.iter().enumerate() {
        if *is_decoy {
            continue;
        }
        let mut g = groups[*group_idx].clone();
        g.q_value = q_map.get(&winner_i).copied();
        result_groups.push(g);
    }

    // Sort by score descending, then leader accession ascending for a fully
    // deterministic output order even when scores (and q-values) are tied.
    result_groups.sort_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then_with(|| a.leader_accession.cmp(&b.leader_accession))
    });

    let target_groups_at_1pct = result_groups
        .iter()
        .filter(|g| g.q_value.is_some_and(|q| q <= 0.01))
        .count() as u64;

    Ok(ProteinFdrResult {
        groups: result_groups,
        target_groups_at_1pct,
        total_target_groups,
        total_decoy_groups,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_group(accession: &str, score: f64, is_decoy: bool) -> ProteinGroup {
        let leader = accession.to_string();
        ProteinGroup {
            leader_accession: leader.clone(),
            leader_description: format!("{accession} protein"),
            member_accessions: vec![leader],
            peptides: vec!["PEPTIDE".to_string()],
            unique_peptides: vec!["PEPTIDE".to_string()],
            razor_peptides: vec![],
            score,
            q_value: None,
            coverage: None,
            is_decoy,
        }
    }

    #[test]
    fn perfect_competition_targets_win() {
        let groups = vec![
            make_group("P001", 20.0, false),
            make_group("P002", 15.0, false),
            make_group("REV_P001", 5.0, true),
            make_group("REV_P002", 3.0, true),
        ];
        let result = calculate_protein_fdr(&groups).unwrap();

        assert_eq!(result.total_target_groups, 2);
        assert_eq!(result.total_decoy_groups, 2);
        // Both targets win — with no decoy winners, all get q=0
        assert_eq!(result.groups.len(), 2);
        for g in &result.groups {
            assert!(!g.is_decoy);
            assert_eq!(g.q_value, Some(0.0));
        }
        assert_eq!(result.target_groups_at_1pct, 2);
    }

    #[test]
    fn decoy_wins_pair() {
        // P001 target=10, REV_P001 decoy=20 → decoy wins
        // P002 target=15, REV_P002 decoy=3 → target wins
        // Winners: REV_P001 (decoy), P002 (target)
        // FDR at P002: 1 decoy / 1 target = 1.0
        let groups = vec![
            make_group("P001", 10.0, false),
            make_group("P002", 15.0, false),
            make_group("REV_P001", 20.0, true),
            make_group("REV_P002", 3.0, true),
        ];
        let result = calculate_protein_fdr(&groups).unwrap();

        // Only P002 should be in the result (winning target)
        assert_eq!(result.groups.len(), 1);
        assert_eq!(result.groups[0].leader_accession, "P002");
        assert!(result.groups[0].q_value.is_some());
    }

    #[test]
    fn unpaired_targets_win_automatically() {
        // P001 paired with REV_P001 → target wins
        // P003 has no matching decoy → wins automatically
        let groups = vec![
            make_group("P001", 20.0, false),
            make_group("P003", 12.0, false),
            make_group("REV_P001", 5.0, true),
        ];
        let result = calculate_protein_fdr(&groups).unwrap();

        assert_eq!(result.total_target_groups, 2);
        assert_eq!(result.total_decoy_groups, 1);
        // Both targets win, no decoy winners → all q=0
        assert_eq!(result.groups.len(), 2);
        let accessions: Vec<&str> = result
            .groups
            .iter()
            .map(|g| g.leader_accession.as_str())
            .collect();
        assert!(accessions.contains(&"P001"));
        assert!(accessions.contains(&"P003"));
        for g in &result.groups {
            assert_eq!(g.q_value, Some(0.0));
        }
    }

    #[test]
    fn unpaired_decoys_included_in_fdr() {
        // P001 paired with REV_P001 → target wins (score 20 > 5)
        // REV_P099 has no matching target → unpaired decoy winner
        // Winners: P001 (target, score 20), REV_P099 (decoy, score 18)
        // Sorted: P001(20,T), REV_P099(18,D)
        // FDR walk: pos1: 0/1=0, pos2: 1/1=1.0
        // Monotonicity: q[0]=min(0, 1.0)=0, q[1]=1.0
        // P001 gets q=0
        let groups = vec![
            make_group("P001", 20.0, false),
            make_group("REV_P001", 5.0, true),
            make_group("REV_P099", 18.0, true),
        ];
        let result = calculate_protein_fdr(&groups).unwrap();

        assert_eq!(result.groups.len(), 1);
        assert_eq!(result.groups[0].leader_accession, "P001");
        // The unpaired decoy influences FDR but target still has good q-value
        assert!(result.groups[0].q_value.unwrap() <= 0.01);
    }

    #[test]
    fn all_targets_no_decoys_returns_all_q_zero() {
        let groups = vec![
            make_group("P001", 20.0, false),
            make_group("P002", 15.0, false),
        ];
        // No decoys at all → no decoy winners → all targets get q=0
        let result = calculate_protein_fdr(&groups).unwrap();

        assert_eq!(result.groups.len(), 2);
        assert_eq!(result.total_decoy_groups, 0);
        for g in &result.groups {
            assert_eq!(g.q_value, Some(0.0));
        }
    }

    #[test]
    fn single_pair_target_wins() {
        let groups = vec![
            make_group("P001", 15.0, false),
            make_group("REV_P001", 8.0, true),
        ];
        let result = calculate_protein_fdr(&groups).unwrap();

        assert_eq!(result.groups.len(), 1);
        assert_eq!(result.groups[0].leader_accession, "P001");
        // Single target winner with no decoy winners → q=0
        assert_eq!(result.groups[0].q_value, Some(0.0));
    }

    #[test]
    fn single_pair_decoy_wins() {
        let groups = vec![
            make_group("P001", 5.0, false),
            make_group("REV_P001", 15.0, true),
        ];
        let result = calculate_protein_fdr(&groups).unwrap();

        // Decoy won, no target winners → empty result
        assert_eq!(result.groups.len(), 0);
        assert_eq!(result.target_groups_at_1pct, 0);
    }

    #[test]
    fn score_tie_target_wins() {
        let groups = vec![
            make_group("P001", 10.0, false),
            make_group("REV_P001", 10.0, true),
        ];
        let result = calculate_protein_fdr(&groups).unwrap();

        // Tie → target wins
        assert_eq!(result.groups.len(), 1);
        assert_eq!(result.groups[0].leader_accession, "P001");
    }

    #[test]
    fn q_value_monotonicity() {
        // Create a realistic scenario with mixed winners.
        // Multiple targets and decoys so that FDR changes across the list.
        let groups = vec![
            make_group("P001", 50.0, false),
            make_group("P002", 40.0, false),
            make_group("P003", 30.0, false),
            make_group("P004", 20.0, false),
            make_group("P005", 10.0, false),
            make_group("REV_P001", 1.0, true),  // target wins
            make_group("REV_P002", 1.0, true),  // target wins
            make_group("REV_P003", 35.0, true), // decoy wins (35 > 30)
            make_group("REV_P004", 25.0, true), // decoy wins (25 > 20)
            make_group("REV_P005", 1.0, true),  // target wins
        ];
        let result = calculate_protein_fdr(&groups).unwrap();

        // Result groups should be sorted by score descending.
        // Verify q-values are monotonically non-decreasing as score decreases.
        for window in result.groups.windows(2) {
            let q_a = window[0].q_value.unwrap();
            let q_b = window[1].q_value.unwrap();
            assert!(
                q_a <= q_b + 1e-12,
                "q-values must be monotonically non-decreasing: {} (score {}) > {} (score {})",
                q_a,
                window[0].score,
                q_b,
                window[1].score,
            );
        }

        // Winners: P001(50,T), P002(40,T), REV_P003(35,D), REV_P004(25,D), P005(10,T)
        // Only target winners in output: P001, P002, P005
        assert_eq!(result.groups.len(), 3);
    }

    #[test]
    fn empty_input_returns_error() {
        let result = calculate_protein_fdr(&[]);
        assert!(matches!(result, Err(FdrError::NoPsms)));
    }

    #[test]
    fn tied_target_winners_deterministic_order() {
        // Five target winners with identical scores, plus an unpaired decoy winner
        // (so the FDR/TDC path runs). Tied groups must be emitted in a deterministic
        // order — sorted by score descending, then leader_accession ascending.
        let build = || {
            vec![
                make_group("P001", 10.0, false),
                make_group("REV_P001", 1.0, true),
                make_group("P002", 10.0, false),
                make_group("REV_P002", 1.0, true),
                make_group("P003", 10.0, false),
                make_group("REV_P003", 1.0, true),
                make_group("P004", 10.0, false),
                make_group("REV_P004", 1.0, true),
                make_group("P005", 10.0, false),
                make_group("REV_P005", 1.0, true),
                make_group("REV_P999", 0.5, true), // unpaired decoy → decoy winner
            ]
        };

        let expected_order = vec!["P001", "P002", "P003", "P004", "P005"];

        // Run many times: each call builds fresh HashMaps (different seeds), so a
        // non-deterministic ordering would surface quickly.
        for _ in 0..16 {
            let result = calculate_protein_fdr(&build()).unwrap();
            let order: Vec<&str> = result
                .groups
                .iter()
                .map(|g| g.leader_accession.as_str())
                .collect();
            assert_eq!(
                order, expected_order,
                "tied target winners must be emitted in deterministic (score desc, leader asc) order"
            );
            // Tied scores must share one q-value (see FIX 1).
            let qs: Vec<f64> = result.groups.iter().map(|g| g.q_value.unwrap()).collect();
            for q in &qs {
                assert!(
                    (q - qs[0]).abs() < f64::EPSILON,
                    "tied winners must share one q-value: {qs:?}"
                );
            }
        }
    }

    #[test]
    fn many_unpaired_decoys_no_panic_and_correct() {
        // Exercises the unpaired-decoy membership check (Vec→HashSet change):
        // mix of paired (target wins) and several unpaired decoys.
        let groups = vec![
            make_group("P001", 30.0, false),
            make_group("REV_P001", 2.0, true), // paired, target wins
            make_group("P002", 28.0, false),
            make_group("REV_P002", 2.0, true), // paired, target wins
            make_group("REV_X1", 1.0, true),   // unpaired decoy
            make_group("REV_X2", 1.0, true),   // unpaired decoy
            make_group("REV_X3", 1.0, true),   // unpaired decoy
        ];
        let result = calculate_protein_fdr(&groups).unwrap();
        let accessions: Vec<&str> = result
            .groups
            .iter()
            .map(|g| g.leader_accession.as_str())
            .collect();
        assert!(accessions.contains(&"P001"));
        assert!(accessions.contains(&"P002"));
        // Only the two winning targets are returned.
        assert_eq!(result.groups.len(), 2);
    }

    #[test]
    fn output_sorted_by_score_descending() {
        let groups = vec![
            make_group("P002", 10.0, false),
            make_group("P001", 30.0, false),
            make_group("P003", 20.0, false),
            make_group("REV_P001", 1.0, true),
            make_group("REV_P002", 1.0, true),
            make_group("REV_P003", 1.0, true),
        ];
        let result = calculate_protein_fdr(&groups).unwrap();

        for window in result.groups.windows(2) {
            assert!(
                window[0].score >= window[1].score,
                "groups should be sorted by score descending: {} < {}",
                window[0].score,
                window[1].score,
            );
        }
    }
}
