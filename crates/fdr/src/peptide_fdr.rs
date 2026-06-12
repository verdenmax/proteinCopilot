//! Peptide-level FDR calculation using target-decoy competition.
//!
//! Groups PSMs by peptide sequence, keeps the best score per sequence,
//! then applies the same TDC algorithm used for PSM-level FDR.

use std::collections::HashMap;

use crate::calculation::{calculate_fdr, ScoredPsm};
use crate::error::FdrError;

/// Result of peptide-level FDR calculation.
#[derive(Debug, Clone)]
pub struct PeptideFdrResult {
    /// peptide_sequence → q_value
    pub peptide_q_values: HashMap<String, f64>,
    /// Number of target peptides at 1% FDR
    pub target_peptides_at_1pct: u64,
    /// Number of decoy peptides
    pub total_decoy_peptides: u64,
    /// Total unique peptide sequences
    pub total_peptides: u64,
}

/// Input: a peptide with its best PSM score and decoy status.
#[derive(Debug, Clone)]
pub struct PeptideScore {
    /// Peptide amino acid sequence (should be I/L normalized if desired).
    pub sequence: String,
    /// Best PSM score for this peptide (highest score among all PSMs).
    pub best_score: f64,
    /// Whether this peptide is from the decoy database.
    pub is_decoy: bool,
}

/// Calculate peptide-level FDR using target-decoy competition.
///
/// Algorithm:
/// 1. Each unique peptide is represented by its best-scoring PSM
/// 2. Apply TDC (target-decoy competition) at peptide level
/// 3. Compute q-values with monotonicity enforcement
///
/// This reuses the existing [`calculate_fdr`] function which handles
/// the core TDC algorithm.
pub fn calculate_peptide_fdr(peptides: &[PeptideScore]) -> Result<PeptideFdrResult, FdrError> {
    let scored: Vec<ScoredPsm> = peptides
        .iter()
        .enumerate()
        .map(|(i, p)| ScoredPsm {
            index: i,
            score: p.best_score,
            is_decoy: p.is_decoy,
        })
        .collect();

    let fdr_results = calculate_fdr(&scored)?;

    let mut peptide_q_values = HashMap::with_capacity(fdr_results.len());
    for (idx, q) in &fdr_results {
        peptide_q_values.insert(peptides[*idx].sequence.clone(), *q);
    }

    let total_decoy_peptides = peptides.iter().filter(|p| p.is_decoy).count() as u64;
    let total_peptides = peptides.len() as u64;

    let target_peptides_at_1pct = fdr_results
        .iter()
        .filter(|(idx, q)| !peptides[*idx].is_decoy && *q <= 0.01)
        .count() as u64;

    Ok(PeptideFdrResult {
        peptide_q_values,
        target_peptides_at_1pct,
        total_decoy_peptides,
        total_peptides,
    })
}

/// Extract unique peptides from PSMs, keeping best score per sequence.
///
/// For each unique `peptide_sequence`:
/// - Keep the highest score
/// - Use the `is_decoy` flag from the best-scoring PSM
///
/// Note: This function does NOT normalize I/L — the caller should
/// normalize first if desired.
pub fn extract_unique_peptides(
    psms: &[protein_copilot_core::search_result::Psm],
) -> Vec<PeptideScore> {
    let mut best_per_seq: HashMap<&str, (f64, bool)> = HashMap::new();

    for psm in psms {
        let entry = best_per_seq
            .entry(psm.peptide_sequence.as_str())
            .or_insert((psm.score, psm.is_decoy));
        if psm.score > entry.0 {
            *entry = (psm.score, psm.is_decoy);
        }
    }

    let mut result: Vec<PeptideScore> = best_per_seq
        .into_iter()
        .map(|(seq, (score, is_decoy))| PeptideScore {
            sequence: seq.to_string(),
            best_score: score,
            is_decoy,
        })
        .collect();
    // Deterministic ordering so downstream tie-breaking is reproducible.
    result.sort_by(|a, b| a.sequence.cmp(&b.sequence));
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_peptide(seq: &str, score: f64, is_decoy: bool) -> PeptideScore {
        PeptideScore {
            sequence: seq.to_string(),
            best_score: score,
            is_decoy,
        }
    }

    fn make_psm(seq: &str, score: f64, is_decoy: bool) -> protein_copilot_core::search_result::Psm {
        protein_copilot_core::search_result::Psm {
            spectrum_scan: 1,
            peptide_sequence: seq.to_string(),
            modifications: vec![],
            charge: 2,
            precursor_mz: 500.0,
            calculated_mz: 500.0,
            delta_mass_ppm: 0.0,
            score,
            q_value: None,
            protein_accessions: vec!["P12345".to_string()],
            is_decoy,
            extra: None,
        }
    }

    #[test]
    fn test_empty_peptides() {
        let result = calculate_peptide_fdr(&[]);
        assert!(matches!(result, Err(FdrError::NoPsms)));
    }

    #[test]
    fn test_single_target_peptide() {
        let peptides = vec![make_peptide("PEPTIDER", 10.0, false)];
        let result = calculate_peptide_fdr(&peptides);
        assert!(matches!(result, Err(FdrError::NoDecoyHits)));
    }

    #[test]
    fn test_basic_peptide_fdr() {
        let peptides = vec![
            make_peptide("PEPTIDER", 10.0, false),
            make_peptide("ANOTHERSEQ", 8.0, false),
            make_peptide("THIRDPEP", 6.0, false),
            make_peptide("REV_DECOY", 4.0, true),
        ];
        let result = calculate_peptide_fdr(&peptides).unwrap();

        assert_eq!(result.total_peptides, 4);
        assert_eq!(result.total_decoy_peptides, 1);
        assert_eq!(result.peptide_q_values.len(), 4);

        // All 4 peptides should have a q-value assigned
        for p in &peptides {
            assert!(
                result.peptide_q_values.contains_key(&p.sequence),
                "missing q-value for {}",
                p.sequence
            );
        }
    }

    #[test]
    fn test_best_score_selection() {
        let psms = vec![
            make_psm("PEPTIDER", 5.0, false),
            make_psm("PEPTIDER", 10.0, false),
            make_psm("PEPTIDER", 7.0, false),
            make_psm("OTHER", 3.0, true),
        ];
        let unique = extract_unique_peptides(&psms);

        let peptider = unique.iter().find(|p| p.sequence == "PEPTIDER").unwrap();
        assert!(
            (peptider.best_score - 10.0).abs() < f64::EPSILON,
            "should pick max score 10.0, got {}",
            peptider.best_score
        );
        assert_eq!(unique.len(), 2);
    }

    #[test]
    fn test_decoy_peptide_handling() {
        let peptides = vec![
            make_peptide("TARGET1", 10.0, false),
            make_peptide("DECOY1", 9.0, true),
            make_peptide("TARGET2", 5.0, false),
            make_peptide("DECOY2", 3.0, true),
        ];
        let result = calculate_peptide_fdr(&peptides).unwrap();

        assert_eq!(result.total_decoy_peptides, 2);
        // Decoys should have q-values assigned
        assert!(result.peptide_q_values.contains_key("DECOY1"));
        assert!(result.peptide_q_values.contains_key("DECOY2"));

        // All q-values should be in [0, 1]
        for q in result.peptide_q_values.values() {
            assert!(*q >= 0.0 && *q <= 1.0, "q-value out of range: {q}");
        }
    }

    #[test]
    fn test_target_count_at_1pct() {
        // 3 high-scoring targets, then 1 decoy at lowest score
        // With TDC: first 3 positions are targets → FDR=0/1, 0/2, 0/3 → q=0
        // Position 4 is decoy → FDR=1/3 ≈ 0.33
        // After monotonicity: targets all have q ≤ 0.01
        let peptides = vec![
            make_peptide("TARGET1", 10.0, false),
            make_peptide("TARGET2", 9.0, false),
            make_peptide("TARGET3", 8.0, false),
            make_peptide("DECOY1", 2.0, true),
        ];
        let result = calculate_peptide_fdr(&peptides).unwrap();
        assert_eq!(
            result.target_peptides_at_1pct, 3,
            "expected 3 targets at 1% FDR, got {}",
            result.target_peptides_at_1pct
        );
    }

    #[test]
    fn test_monotonic_q_values() {
        let peptides = vec![
            make_peptide("A", 10.0, false),
            make_peptide("B", 9.0, true),
            make_peptide("C", 8.0, false),
            make_peptide("D", 7.0, false),
            make_peptide("E", 6.0, true),
            make_peptide("F", 5.0, false),
        ];
        let result = calculate_peptide_fdr(&peptides).unwrap();

        // Collect q-values sorted by score descending
        let mut scored_q: Vec<(f64, f64)> = peptides
            .iter()
            .map(|p| (p.best_score, result.peptide_q_values[&p.sequence]))
            .collect();
        scored_q.sort_by(|a, b| b.0.total_cmp(&a.0));

        for window in scored_q.windows(2) {
            assert!(
                window[0].1 <= window[1].1 + 1e-12,
                "q-values must be monotonically non-decreasing with decreasing score: {:?}",
                scored_q
            );
        }
    }

    #[test]
    fn test_extract_unique_peptides_sorted_and_stable() {
        // Sequences inserted in non-alphabetical order, with a duplicate ("AAAA").
        // extract_unique_peptides must return them sorted by sequence (deterministic),
        // while preserving best-score-per-sequence semantics.
        let psms = vec![
            make_psm("DDDD", 5.0, false),
            make_psm("AAAA", 3.0, false),
            make_psm("CCCC", 7.0, false),
            make_psm("AAAA", 9.0, false), // duplicate, higher score wins
            make_psm("BBBB", 4.0, true),
            make_psm("EEEE", 2.0, false),
        ];

        let unique = extract_unique_peptides(&psms);
        let seqs: Vec<String> = unique.iter().map(|p| p.sequence.clone()).collect();
        assert_eq!(
            seqs,
            vec!["AAAA", "BBBB", "CCCC", "DDDD", "EEEE"],
            "output must be deterministically sorted by sequence"
        );

        // Stable across repeated calls.
        let again: Vec<String> = extract_unique_peptides(&psms)
            .iter()
            .map(|p| p.sequence.clone())
            .collect();
        assert_eq!(seqs, again, "ordering must be stable across runs");

        // Best-score-per-sequence semantics intact.
        let aaaa = unique.iter().find(|p| p.sequence == "AAAA").unwrap();
        assert!((aaaa.best_score - 9.0).abs() < f64::EPSILON);
        let bbbb = unique.iter().find(|p| p.sequence == "BBBB").unwrap();
        assert!(bbbb.is_decoy);
        assert_eq!(unique.len(), 5);
    }

    #[test]
    fn test_extract_from_psms() {
        let psms = vec![
            make_psm("AAAA", 10.0, false),
            make_psm("AAAA", 8.0, false),
            make_psm("BBBB", 6.0, true),
            make_psm("CCCC", 4.0, false),
            make_psm("CCCC", 9.0, false),
        ];
        let unique = extract_unique_peptides(&psms);
        assert_eq!(unique.len(), 3);

        let aaaa = unique.iter().find(|p| p.sequence == "AAAA").unwrap();
        assert!((aaaa.best_score - 10.0).abs() < f64::EPSILON);
        assert!(!aaaa.is_decoy);

        let bbbb = unique.iter().find(|p| p.sequence == "BBBB").unwrap();
        assert!((bbbb.best_score - 6.0).abs() < f64::EPSILON);
        assert!(bbbb.is_decoy);

        let cccc = unique.iter().find(|p| p.sequence == "CCCC").unwrap();
        assert!((cccc.best_score - 9.0).abs() < f64::EPSILON);
        assert!(!cccc.is_decoy);
    }
}
