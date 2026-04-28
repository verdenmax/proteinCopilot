//! Sequence coverage calculation.
//!
//! Computes the fraction of a protein's amino acid sequence that is
//! covered by identified peptides. Requires the FASTA database for
//! full protein sequences.

use std::collections::{HashMap, HashSet};

use protein_copilot_core::protein_group::ProteinGroup;

/// Normalize I/L equivalence: replace all 'I' with 'L'.
fn normalize_il(s: &str) -> String {
    s.replace('I', "L")
}

/// Calculate sequence coverage for protein groups.
///
/// For each group, finds where its peptides map in the leader protein's FASTA sequence,
/// then computes coverage = (covered residues) / (total residues).
///
/// `fasta_sequences`: protein_accession → amino acid sequence (from FASTA file)
///
/// Mutates groups in-place to set the `coverage` field. Proteins not found
/// in the FASTA map get coverage = None.
pub fn calculate_coverage(groups: &mut [ProteinGroup], fasta_sequences: &HashMap<String, String>) {
    let _span = tracing::info_span!("calculate_coverage", group_count = groups.len()).entered();

    let total = groups.len();
    let progress_interval: usize = 500;
    let loop_start = std::time::Instant::now();

    for (i, group) in groups.iter_mut().enumerate() {
        let Some(fasta_seq) = fasta_sequences.get(&group.leader_accession) else {
            tracing::warn!(
                accession = %group.leader_accession,
                "Protein not found in FASTA database; coverage set to None"
            );
            group.coverage = None;
            continue;
        };

        let seq_len = fasta_seq.len();
        if seq_len == 0 {
            group.coverage = None;
            continue;
        }

        let normalized_fasta = normalize_il(fasta_seq);
        let mut covered = vec![false; seq_len];

        // Combine peptides and razor_peptides, deduplicated
        let all_peptides: HashSet<&str> = group
            .peptides
            .iter()
            .chain(group.razor_peptides.iter())
            .map(String::as_str)
            .collect();

        for peptide in &all_peptides {
            let normalized_pep = normalize_il(peptide);
            // Find all occurrences of this peptide in the normalized FASTA sequence
            let mut start = 0;
            while let Some(pos) = normalized_fasta[start..].find(&normalized_pep) {
                let abs_pos = start + pos;
                for flag in covered.iter_mut().skip(abs_pos).take(normalized_pep.len()) {
                    *flag = true;
                }
                start = abs_pos + 1;
            }
        }

        let covered_count = covered.iter().filter(|&&b| b).count();
        group.coverage = Some(covered_count as f64 / seq_len as f64);

        if (i + 1) % progress_interval == 0 || i + 1 == total {
            let elapsed = loop_start.elapsed().as_secs_f64();
            let rate = if elapsed > 0.0 { (i + 1) as f64 / elapsed } else { 0.0 };
            let eta = if rate > 0.0 { (total - i - 1) as f64 / rate } else { 0.0 };
            tracing::info!(
                progress = i + 1,
                total = total,
                rate = format!("{:.0}/s", rate),
                eta_sec = format!("{:.1}", eta),
                "calculating coverage"
            );
        }
    }

    let with_coverage = groups.iter().filter(|g| g.coverage.is_some()).count();
    tracing::info!(groups_with_coverage = with_coverage, "coverage calculation complete");
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_group(accession: &str, peptides: Vec<&str>, razor_peptides: Vec<&str>) -> ProteinGroup {
        ProteinGroup {
            leader_accession: accession.to_string(),
            leader_description: String::new(),
            member_accessions: vec![accession.to_string()],
            peptides: peptides.into_iter().map(String::from).collect(),
            unique_peptides: Vec::new(),
            razor_peptides: razor_peptides.into_iter().map(String::from).collect(),
            score: 0.0,
            q_value: None,
            coverage: None,
            is_decoy: false,
        }
    }

    #[test]
    fn full_coverage() {
        let mut groups = vec![make_group("P1", vec!["ACDEFG"], vec![])];
        let fasta: HashMap<String, String> = [("P1".into(), "ACDEFG".into())].into_iter().collect();
        calculate_coverage(&mut groups, &fasta);
        assert_eq!(groups[0].coverage, Some(1.0));
    }

    #[test]
    fn partial_coverage() {
        let mut groups = vec![make_group("P1", vec!["ACDE"], vec![])];
        let fasta: HashMap<String, String> =
            [("P1".into(), "ACDEFGHI".into())].into_iter().collect();
        calculate_coverage(&mut groups, &fasta);
        // 4 out of 8 = 0.5
        assert_eq!(groups[0].coverage, Some(0.5));
    }

    #[test]
    fn overlapping_peptides() {
        // Protein: ABCDEFGH (len 8)
        // Peptide1: ABCDE (pos 0..5)
        // Peptide2: CDEFGH (pos 2..8)
        // Union covers all 8 positions → 1.0
        let mut groups = vec![make_group("P1", vec!["ABCDE", "CDEFGH"], vec![])];
        let fasta: HashMap<String, String> =
            [("P1".into(), "ABCDEFGH".into())].into_iter().collect();
        calculate_coverage(&mut groups, &fasta);
        assert_eq!(groups[0].coverage, Some(1.0));
    }

    #[test]
    fn non_overlapping_peptides() {
        // Protein: ABCDEFGHIJ (len 10)
        // Peptide1: ABC (pos 0..3)
        // Peptide2: HIJ (pos 7..10)
        // Covered: 6 / 10 = 0.6
        let mut groups = vec![make_group("P1", vec!["ABC", "HIJ"], vec![])];
        let fasta: HashMap<String, String> =
            [("P1".into(), "ABCDEFGHIJ".into())].into_iter().collect();
        calculate_coverage(&mut groups, &fasta);
        assert_eq!(groups[0].coverage, Some(0.6));
    }

    #[test]
    fn peptide_not_found_in_sequence() {
        // Peptide "XYZ" doesn't exist in protein "ABCDEF"
        let mut groups = vec![make_group("P1", vec!["ABC", "XYZ"], vec![])];
        let fasta: HashMap<String, String> = [("P1".into(), "ABCDEF".into())].into_iter().collect();
        calculate_coverage(&mut groups, &fasta);
        // Only ABC matches → 3/6 = 0.5
        assert_eq!(groups[0].coverage, Some(0.5));
    }

    #[test]
    fn protein_not_in_fasta() {
        let mut groups = vec![make_group("MISSING", vec!["ABC"], vec![])];
        let fasta: HashMap<String, String> = HashMap::new();
        calculate_coverage(&mut groups, &fasta);
        assert_eq!(groups[0].coverage, None);
    }

    #[test]
    fn il_equivalence() {
        // Peptide has L, FASTA has I → should still match via normalization
        let mut groups = vec![make_group("P1", vec!["ALCD"], vec![])];
        let fasta: HashMap<String, String> = [("P1".into(), "AICD".into())].into_iter().collect();
        calculate_coverage(&mut groups, &fasta);
        assert_eq!(groups[0].coverage, Some(1.0));
    }

    #[test]
    fn empty_groups() {
        let mut groups: Vec<ProteinGroup> = Vec::new();
        let fasta: HashMap<String, String> = HashMap::new();
        calculate_coverage(&mut groups, &fasta);
        assert!(groups.is_empty());
    }

    #[test]
    fn multiple_occurrences() {
        // Protein: ABCABC (len 6), peptide ABC occurs at pos 0 and 3
        let mut groups = vec![make_group("P1", vec!["ABC"], vec![])];
        let fasta: HashMap<String, String> = [("P1".into(), "ABCABC".into())].into_iter().collect();
        calculate_coverage(&mut groups, &fasta);
        // Both occurrences covered → 6/6 = 1.0
        assert_eq!(groups[0].coverage, Some(1.0));
    }

    #[test]
    fn razor_peptides_counted() {
        // Protein: ABCDEFGHIJ (len 10)
        // peptides: ABC (0..3)
        // razor_peptides: HIJ (7..10)
        let mut groups = vec![make_group("P1", vec!["ABC"], vec!["HIJ"])];
        let fasta: HashMap<String, String> =
            [("P1".into(), "ABCDEFGHIJ".into())].into_iter().collect();
        calculate_coverage(&mut groups, &fasta);
        // 6/10 = 0.6
        assert_eq!(groups[0].coverage, Some(0.6));
    }

    #[test]
    fn empty_protein_sequence() {
        let mut groups = vec![make_group("P1", vec!["ABC"], vec![])];
        let fasta: HashMap<String, String> = [("P1".into(), String::new())].into_iter().collect();
        calculate_coverage(&mut groups, &fasta);
        assert_eq!(groups[0].coverage, None);
    }

    #[test]
    fn duplicate_peptide_in_both_vecs() {
        // Same peptide in peptides and razor_peptides should not cause double-counting
        let mut groups = vec![make_group("P1", vec!["ABCD"], vec!["ABCD"])];
        let fasta: HashMap<String, String> =
            [("P1".into(), "ABCDEFGH".into())].into_iter().collect();
        calculate_coverage(&mut groups, &fasta);
        // 4/8 = 0.5
        assert_eq!(groups[0].coverage, Some(0.5));
    }
}
