//! L0–L4 similarity classification for trap-database PSM hits.
//!
//! Compares each trap PSM against the target digest index to determine
//! how similar the trap peptide is to known target peptides, producing
//! a discriminability level from L0 (exact match / razor error) to
//! L4 (no close match / true trap hit).

use protein_copilot_search_engine::digest::residue_mass;

use crate::config::SimilarityConfig;
use crate::digest::TargetDigestIndex;
use crate::types::{ClassifiedPsm, DiscriminabilityLevel, PsmGroup, SubstitutionType, UnifiedPsm};

/// Compute Hamming-style character differences between two equal-length sequences.
///
/// Returns `None` if `a` and `b` have different lengths.
/// Otherwise returns `Some((mismatches, delta_mass, diff_positions_string))`:
///
/// - `mismatches`: number of positions where the characters differ.
/// - `delta_mass`: sum of `residue_mass(b_char) - residue_mass(a_char)` for each
///   differing position.  If either character has no known mass the contribution
///   for that position is 0.0.
/// - `diff_positions_string`: formatted as `"[pos:X->Y,pos2:A->B]"` (0-indexed).
///   Empty string `""` when there are no mismatches.
pub fn hamming_diff(a: &str, b: &str) -> Option<(u16, f64, String)> {
    if a.len() != b.len() {
        return None;
    }

    let mut mismatches: u16 = 0;
    let mut delta_mass: f64 = 0.0;
    let mut diffs: Vec<String> = Vec::new();

    for (i, (ca, cb)) in a.chars().zip(b.chars()).enumerate() {
        if ca != cb {
            mismatches = mismatches.saturating_add(1);

            let dm = match (residue_mass(ca), residue_mass(cb)) {
                (Some(ma), Some(mb)) => mb - ma,
                _ => 0.0,
            };
            delta_mass += dm;

            diffs.push(format!("{i}:{ca}->{cb}"));
        }
    }

    let diff_str = if diffs.is_empty() {
        String::new()
    } else {
        format!("[{}]", diffs.join(","))
    };

    Some((mismatches, delta_mass, diff_str))
}

/// Returns `true` when every mismatch position in `a` vs `b` is a Leu / Ile swap.
///
/// Assumes `a` and `b` have the same length.
fn is_only_li_substitution(a: &str, b: &str) -> bool {
    a.chars().zip(b.chars()).all(|(ca, cb)| {
        ca == cb || {
            let pair = (ca, cb);
            pair == ('L', 'I') || pair == ('I', 'L')
        }
    })
}

/// Classify a single PSM against the target digest index.
///
/// Non-`Trap` PSMs are returned unchanged with level L4 and no match info.
/// Trap PSMs are classified through the L0 → L1 → L2/L3/L4 priority chain.
pub fn classify_single(
    psm: &UnifiedPsm,
    group: PsmGroup,
    index: &TargetDigestIndex,
    config: &SimilarityConfig,
) -> ClassifiedPsm {
    // Non-trap PSMs always get L4 with no match info.
    if group != PsmGroup::Trap {
        return ClassifiedPsm {
            psm: psm.clone(),
            group,
            level: DiscriminabilityLevel::L4,
            best_target_peptide: None,
            best_target_protein: None,
            mismatches: None,
            delta_mass_da: None,
            diff_positions: None,
            substitution_type: SubstitutionType::None,
            edit_distance: None,
            alignment_detail: None,
        };
    }

    // --- L0: exact match -------------------------------------------------
    if index.has_exact(&psm.peptide) {
        let protein = index.exact_protein(&psm.peptide).map(|s| s.to_owned());
        return ClassifiedPsm {
            psm: psm.clone(),
            group,
            level: DiscriminabilityLevel::L0,
            best_target_peptide: Some(psm.peptide.clone()),
            best_target_protein: protein,
            mismatches: Some(0),
            delta_mass_da: Some(0.0),
            diff_positions: Some(String::new()),
            substitution_type: SubstitutionType::None,
            edit_distance: None,
            alignment_detail: None,
        };
    }

    // --- L1: L/I-normalized match (but NOT exact) -------------------------
    if index.has_normalized(&psm.peptide) {
        if let Some((orig, prot)) = index.normalized_match(&psm.peptide) {
            // Defensive: if the original target == trap peptide, this is actually L0
            // (should have been caught above, but guard against index inconsistencies)
            if orig == psm.peptide {
                return ClassifiedPsm {
                    psm: psm.clone(),
                    group,
                    level: DiscriminabilityLevel::L0,
                    best_target_peptide: Some(psm.peptide.clone()),
                    best_target_protein: Some(prot.to_owned()),
                    mismatches: Some(0),
                    delta_mass_da: Some(0.0),
                    diff_positions: Some(String::new()),
                    substitution_type: SubstitutionType::None,
                    edit_distance: None,
                    alignment_detail: None,
                };
            }
            let (mm, dm, dp) = hamming_diff(&psm.peptide, orig).unwrap_or((0, 0.0, String::new()));
            return ClassifiedPsm {
                psm: psm.clone(),
                group,
                level: DiscriminabilityLevel::L1,
                best_target_peptide: Some(orig.to_owned()),
                best_target_protein: Some(prot.to_owned()),
                mismatches: Some(mm),
                delta_mass_da: Some(dm),
                diff_positions: Some(dp),
                substitution_type: SubstitutionType::LIIsomer,
                edit_distance: None,
                alignment_detail: None,
            };
        }
    }

    // --- L2/L3/L4: brute-force hamming scan -------------------------------
    let candidates = index.peptides_of_length(psm.peptide.len());

    let mut best_mm: u16 = u16::MAX;
    let mut best_dm: f64 = f64::MAX;
    let mut best_dp = String::new();
    let mut best_seq: Option<&str> = None;
    let mut best_prot: Option<&str> = None;

    for target in candidates {
        let (mm, dm, dp) = match hamming_diff(&psm.peptide, &target.sequence) {
            Some(v) => v,
            None => continue, // shouldn't happen – same length
        };

        if mm == 0 {
            // exact match already handled above
            continue;
        }
        if mm > config.max_mismatches {
            continue;
        }
        // Skip pure L/I substitutions – those are L1 territory.
        if is_only_li_substitution(&psm.peptide, &target.sequence) {
            continue;
        }

        let abs_dm = dm.abs();
        if mm < best_mm || (mm == best_mm && abs_dm < best_dm) {
            best_mm = mm;
            best_dm = abs_dm;
            best_dp = dp;
            best_seq = Some(&target.sequence);
            best_prot = Some(&target.protein_accession);
        }
    }

    // Decide level from best match
    if best_mm == u16::MAX {
        // No close match found → L4
        return ClassifiedPsm {
            psm: psm.clone(),
            group,
            level: DiscriminabilityLevel::L4,
            best_target_peptide: None,
            best_target_protein: None,
            mismatches: None,
            delta_mass_da: None,
            diff_positions: None,
            substitution_type: SubstitutionType::None,
            edit_distance: None,
            alignment_detail: None,
        };
    }

    let level = if best_mm == 1 && best_dm < config.delta_mass_threshold_da {
        DiscriminabilityLevel::L2
    } else {
        DiscriminabilityLevel::L3
    };

    ClassifiedPsm {
        psm: psm.clone(),
        group,
        level,
        best_target_peptide: best_seq.map(|s| s.to_owned()),
        best_target_protein: best_prot.map(|s| s.to_owned()),
        mismatches: Some(best_mm),
        delta_mass_da: Some(best_dm),
        diff_positions: Some(best_dp),
        substitution_type: SubstitutionType::None,
        edit_distance: None,
        alignment_detail: None,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::digest::TargetPeptide;

    /// Helper to create a minimal [`UnifiedPsm`] for testing.
    fn make_psm(peptide: &str) -> UnifiedPsm {
        UnifiedPsm {
            peptide: peptide.to_owned(),
            charge: None,
            precursor_mz: None,
            retention_time: None,
            scan_number: None,
            spectrum_file: None,
            protein_ids: String::new(),
            q_value: None,
        }
    }

    // --- hamming_diff tests -----------------------------------------------

    #[test]
    fn test_hamming_identical() {
        let (mm, dm, dp) = hamming_diff("AGCDEK", "AGCDEK").expect("same length");
        assert_eq!(mm, 0);
        assert!((dm - 0.0).abs() < f64::EPSILON);
        assert_eq!(dp, "");
    }

    #[test]
    fn test_hamming_one_mismatch() {
        // G -> N substitution at position 1
        let (mm, _dm, dp) = hamming_diff("AGCDEK", "ANCDEK").expect("same length");
        assert_eq!(mm, 1);
        assert_eq!(dp, "[1:G->N]");
    }

    #[test]
    fn test_hamming_two_mismatches() {
        // G -> N at pos 1, D -> E at pos 3
        let (mm, _dm, dp) = hamming_diff("AGCDEK", "ANCEER").expect("same length");
        // positions 1 (G->N), 3 (D->E), 5 (K->R)
        assert_eq!(mm, 3);
        assert!(dp.starts_with('['));
        assert!(dp.contains("1:G->N"));
        assert!(dp.contains("3:D->E"));
        assert!(dp.contains("5:K->R"));

        // Two-mismatch variant
        let (mm, _dm, dp) = hamming_diff("AGCDEK", "ANCDER").expect("same length");
        assert_eq!(mm, 2);
        assert!(dp.contains("1:G->N"));
        assert!(dp.contains("5:K->R"));
    }

    #[test]
    fn test_hamming_different_length() {
        assert!(hamming_diff("ABC", "ABCD").is_none());
        assert!(hamming_diff("ABCD", "AB").is_none());
    }

    #[test]
    fn test_hamming_mass_difference() {
        // D mass = 115.026943 Da, N mass = 114.042927 Da
        // delta = 114.042927 - 115.026943 = −0.984016
        let (mm, dm, dp) = hamming_diff("DGFLLDGFPR", "NGFLLDGFPR").expect("same length");
        assert_eq!(mm, 1);
        assert!((dm.abs() - 0.984016).abs() < 0.001, "delta = {dm}");
        assert!(dm < 0.0, "N is lighter than D, so delta should be negative");
        assert_eq!(dp, "[0:D->N]");
    }

    // --- classify_single tests --------------------------------------------

    #[test]
    fn test_classify_target_psm_gets_l4() {
        let psm = make_psm("PEPTIDEK");
        let index = TargetDigestIndex::empty_for_test();
        let config = SimilarityConfig::default();

        let result = classify_single(&psm, PsmGroup::Target, &index, &config);
        assert_eq!(result.level, DiscriminabilityLevel::L4);
        assert_eq!(result.group, PsmGroup::Target);
        assert!(result.best_target_peptide.is_none());
    }

    #[test]
    fn test_classify_trap_l0_exact_match() {
        let psm = make_psm("PEPTIDEK");
        let mut index = TargetDigestIndex::empty_for_test();
        index.exact_set.insert("PEPTIDEK".to_owned());
        index
            .exact_to_protein
            .insert("PEPTIDEK".to_owned(), "P00001".to_owned());
        let config = SimilarityConfig::default();

        let result = classify_single(&psm, PsmGroup::Trap, &index, &config);
        assert_eq!(result.level, DiscriminabilityLevel::L0);
        assert_eq!(result.best_target_peptide.as_deref(), Some("PEPTIDEK"));
        assert_eq!(result.best_target_protein.as_deref(), Some("P00001"));
        assert_eq!(result.mismatches, Some(0));
    }

    #[test]
    fn test_classify_trap_l1_li_isomer() {
        // Trap peptide has I where target has L
        let psm = make_psm("PEPTIDEK");
        let mut index = TargetDigestIndex::empty_for_test();
        // Not an exact match – "PEPTLDEK" is in the target DB, not "PEPTIDEK"
        let norm = "PEPTLDEK"; // normalize_li("PEPTIDEK") = "PEPTLDEK"
        index.normalized_set.insert(norm.to_owned());
        index.normalized_to_original.insert(
            norm.to_owned(),
            ("PEPTLDEK".to_owned(), "P00002".to_owned()),
        );
        let config = SimilarityConfig::default();

        let result = classify_single(&psm, PsmGroup::Trap, &index, &config);
        assert_eq!(result.level, DiscriminabilityLevel::L1);
        assert_eq!(result.best_target_peptide.as_deref(), Some("PEPTLDEK"));
        assert_eq!(result.best_target_protein.as_deref(), Some("P00002"));
    }

    #[test]
    fn test_classify_trap_l2_near_isobaric() {
        // One mismatch with small delta mass → L2
        // D→N has delta ≈ 0.984 Da which is < 1.0 threshold
        let psm = make_psm("DGFLLDGFPR");
        let mut index = TargetDigestIndex::empty_for_test();
        index.by_length.insert(
            10,
            vec![TargetPeptide {
                sequence: "NGFLLDGFPR".to_owned(),
                protein_accession: "P00003".to_owned(),
                neutral_mass: 0.0,
            }],
        );
        let config = SimilarityConfig::default();

        let result = classify_single(&psm, PsmGroup::Trap, &index, &config);
        assert_eq!(result.level, DiscriminabilityLevel::L2);
        assert_eq!(result.mismatches, Some(1));
    }

    #[test]
    fn test_classify_trap_l3_distinguishable() {
        // One mismatch with large delta mass → L3
        // G→W: 186.079313 − 57.021464 = 129.057849 Da >> 1.0
        let psm = make_psm("AGFLLDGFPR");
        let mut index = TargetDigestIndex::empty_for_test();
        index.by_length.insert(
            10,
            vec![TargetPeptide {
                sequence: "WGFLLDGFPR".to_owned(),
                protein_accession: "P00004".to_owned(),
                neutral_mass: 0.0,
            }],
        );
        let config = SimilarityConfig::default();

        let result = classify_single(&psm, PsmGroup::Trap, &index, &config);
        assert_eq!(result.level, DiscriminabilityLevel::L3);
        assert_eq!(result.mismatches, Some(1));
    }

    #[test]
    fn test_classify_trap_l4_no_match() {
        // Empty index → L4
        let psm = make_psm("DGFLLDGFPR");
        let index = TargetDigestIndex::empty_for_test();
        let config = SimilarityConfig::default();

        let result = classify_single(&psm, PsmGroup::Trap, &index, &config);
        assert_eq!(result.level, DiscriminabilityLevel::L4);
        assert!(result.best_target_peptide.is_none());
    }

    #[test]
    fn test_classify_skips_pure_li_in_hamming_scan() {
        // If the only candidates in by_length differ only by L/I, they should
        // be skipped in the hamming scan, resulting in L4 (since there's no
        // L1 match set up via normalized_set either).
        let psm = make_psm("PEPTIDEK");
        let mut index = TargetDigestIndex::empty_for_test();
        index.by_length.insert(
            8,
            vec![TargetPeptide {
                sequence: "PEPTLDEK".to_owned(),
                protein_accession: "P00005".to_owned(),
                neutral_mass: 0.0,
            }],
        );
        let config = SimilarityConfig::default();

        let result = classify_single(&psm, PsmGroup::Trap, &index, &config);
        // Pure L/I substitution is skipped in hamming scan, no L1 set up → L4
        assert_eq!(result.level, DiscriminabilityLevel::L4);
    }
}
