//! Peptide-to-protein mapping.
//!
//! Builds a bidirectional peptide↔protein bipartite graph from PSMs.
//! Handles I/L equivalence, best-score tracking, decoy classification,
//! and optional q-value filtering.

use std::collections::{HashMap, HashSet};

use protein_copilot_core::search_result::Psm;

use crate::error::InferenceError;

/// Decoy protein accession prefix.
use protein_copilot_core::util::is_decoy_accession;

/// Result of peptide-to-protein mapping.
#[derive(Debug, Clone)]
pub struct PeptideProteinMap {
    /// peptide_sequence → set of protein accessions
    pub peptide_to_proteins: HashMap<String, HashSet<String>>,
    /// protein_accession → set of peptide sequences
    pub protein_to_peptides: HashMap<String, HashSet<String>>,
    /// Best PSM score for each peptide sequence (for scoring protein groups later)
    pub peptide_best_score: HashMap<String, f64>,
    /// Whether each peptide is from decoy (true only if ALL its proteins are decoy)
    pub peptide_is_decoy: HashMap<String, bool>,
}

/// Normalize peptide sequence for I/L equivalence.
/// Replaces all 'I' with 'L' for comparison purposes.
pub fn normalize_il(sequence: &str) -> String {
    sequence.replace('I', "L")
}

/// Build a bidirectional peptide-protein map from PSMs.
///
/// - Filters PSMs by `q_value_threshold` when q-values are available
/// - Handles I/L equivalence via [`normalize_il`]
/// - Tracks best score per peptide for later protein scoring
/// - Classifies each peptide as decoy only when ALL mapped proteins are decoy
pub fn build_peptide_protein_map(
    psms: &[Psm],
    q_value_threshold: Option<f64>,
) -> Result<PeptideProteinMap, InferenceError> {
    let _span = tracing::info_span!("build_peptide_protein_map", psm_count = psms.len()).entered();

    if psms.is_empty() {
        return Err(InferenceError::NoPsms);
    }

    let mut peptide_to_proteins: HashMap<String, HashSet<String>> = HashMap::new();
    let mut protein_to_peptides: HashMap<String, HashSet<String>> = HashMap::new();
    let mut peptide_best_score: HashMap<String, f64> = HashMap::new();

    for psm in psms {
        // Filter by q-value threshold when provided.
        // PSMs without a q-value are kept (conservative: don't discard unscored PSMs).
        if let Some(threshold) = q_value_threshold {
            if let Some(q) = psm.q_value {
                if q > threshold {
                    continue;
                }
            }
        }

        let norm_seq = normalize_il(&psm.peptide_sequence);

        // Update peptide → proteins
        let proteins = peptide_to_proteins.entry(norm_seq.clone()).or_default();
        for acc in &psm.protein_accessions {
            proteins.insert(acc.clone());
        }

        // Update protein → peptides
        for acc in &psm.protein_accessions {
            protein_to_peptides
                .entry(acc.clone())
                .or_default()
                .insert(norm_seq.clone());
        }

        // Track best score per peptide (highest wins)
        peptide_best_score
            .entry(norm_seq)
            .and_modify(|best| {
                if psm.score > *best {
                    *best = psm.score;
                }
            })
            .or_insert(psm.score);
    }

    // Classify decoy status: a peptide is decoy iff ALL its proteins are decoy
    let peptide_is_decoy: HashMap<String, bool> = peptide_to_proteins
        .iter()
        .map(|(pep, prots)| {
            let all_decoy = prots.iter().all(|acc| is_decoy_accession(acc));
            (pep.clone(), all_decoy)
        })
        .collect();

    // Validate: at least one target protein must exist
    let has_target = protein_to_peptides
        .keys()
        .any(|acc| !is_decoy_accession(acc));
    if !has_target {
        return Err(InferenceError::NoTargetProteins);
    }

    tracing::info!(
        peptides = peptide_to_proteins.len(),
        proteins = protein_to_peptides.len(),
        "peptide-protein map built"
    );

    Ok(PeptideProteinMap {
        peptide_to_proteins,
        protein_to_peptides,
        peptide_best_score,
        peptide_is_decoy,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::search_result::Psm;

    /// Helper to create a minimal PSM for testing.
    fn make_psm(
        sequence: &str,
        proteins: &[&str],
        score: f64,
        q_value: Option<f64>,
        is_decoy: bool,
    ) -> Psm {
        Psm {
            spectrum_scan: 1,
            peptide_sequence: sequence.to_string(),
            modifications: vec![],
            charge: 2,
            precursor_mz: 500.0,
            calculated_mz: 500.0,
            delta_mass_ppm: 0.0,
            score,
            q_value,
            protein_accessions: proteins.iter().map(|s| s.to_string()).collect(),
            is_decoy,
            extra: None,
        }
    }

    #[test]
    fn test_empty_psms() {
        let result = build_peptide_protein_map(&[], None);
        assert!(
            matches!(result, Err(InferenceError::NoPsms)),
            "expected NoPsms error for empty input"
        );
    }

    #[test]
    fn test_single_psm() {
        let psms = vec![make_psm("PEPTIDE", &["P001"], 10.0, None, false)];
        let map = build_peptide_protein_map(&psms, None).unwrap();

        let norm = normalize_il("PEPTIDE");
        assert!(map.peptide_to_proteins.contains_key(&norm));
        assert!(map.peptide_to_proteins[&norm].contains("P001"));
        assert!(map.protein_to_peptides["P001"].contains(&norm));
        assert_eq!(map.peptide_best_score[&norm], 10.0);
        assert!(!map.peptide_is_decoy[&norm]);
    }

    #[test]
    fn test_il_equivalence() {
        // "PEPTIDE" has I at position 5; "PEPTLDE" has L at position 5.
        // After I→L normalization, both become "PEPTLDE".
        let psms = vec![
            make_psm("PEPTIDE", &["P001"], 10.0, None, false),
            make_psm("PEPTLDE", &["P002"], 8.0, None, false),
        ];

        let map = build_peptide_protein_map(&psms, None).unwrap();

        let norm = normalize_il("PEPTIDE");
        assert_eq!(norm, "PEPTLDE");
        assert_eq!(norm, normalize_il("PEPTLDE"));
        assert_eq!(
            map.peptide_to_proteins[&norm].len(),
            2,
            "I/L-equivalent peptides should map to the same entry"
        );
        assert!(map.peptide_to_proteins[&norm].contains("P001"));
        assert!(map.peptide_to_proteins[&norm].contains("P002"));
    }

    #[test]
    fn test_shared_peptide() {
        let psms = vec![make_psm("SHAREDPEP", &["P001", "P002"], 12.0, None, false)];
        let map = build_peptide_protein_map(&psms, None).unwrap();

        let norm = normalize_il("SHAREDPEP");
        assert_eq!(map.peptide_to_proteins[&norm].len(), 2);
        assert!(map.protein_to_peptides["P001"].contains(&norm));
        assert!(map.protein_to_peptides["P002"].contains(&norm));
    }

    #[test]
    fn test_decoy_classification() {
        let psms = vec![
            make_psm("TARGETPEP", &["P001"], 10.0, None, false),
            make_psm("DECOYPEP", &["REV_P002"], 8.0, None, true),
        ];
        let map = build_peptide_protein_map(&psms, None).unwrap();

        assert!(!map.peptide_is_decoy[&normalize_il("TARGETPEP")]);
        assert!(map.peptide_is_decoy[&normalize_il("DECOYPEP")]);
    }

    #[test]
    fn test_mixed_target_decoy_peptide() {
        // A peptide mapping to both a target and a decoy protein is NOT decoy
        let psms = vec![make_psm(
            "MIXEDPEP",
            &["P001", "REV_P002"],
            10.0,
            None,
            false,
        )];
        let map = build_peptide_protein_map(&psms, None).unwrap();

        assert!(
            !map.peptide_is_decoy[&normalize_il("MIXEDPEP")],
            "peptide mapping to both target and decoy should NOT be classified as decoy"
        );
    }

    #[test]
    fn test_best_score_tracking() {
        let psms = vec![
            make_psm("SCOREPEP", &["P001"], 5.0, None, false),
            make_psm("SCOREPEP", &["P001"], 15.0, None, false),
            make_psm("SCOREPEP", &["P001"], 10.0, None, false),
        ];
        let map = build_peptide_protein_map(&psms, None).unwrap();

        assert_eq!(
            map.peptide_best_score[&normalize_il("SCOREPEP")],
            15.0,
            "should track the highest score"
        );
    }

    #[test]
    fn test_q_value_filtering() {
        let psms = vec![
            make_psm("GOODPEP", &["P001"], 15.0, Some(0.005), false),
            make_psm("BADPEP", &["P002"], 3.0, Some(0.05), false),
            make_psm("NONEQPEP", &["P003"], 8.0, None, false),
        ];
        let map = build_peptide_protein_map(&psms, Some(0.01)).unwrap();

        assert!(
            map.peptide_to_proteins
                .contains_key(&normalize_il("GOODPEP")),
            "PSM with q=0.005 should pass threshold 0.01"
        );
        assert!(
            !map.peptide_to_proteins
                .contains_key(&normalize_il("BADPEP")),
            "PSM with q=0.05 should be filtered at threshold 0.01"
        );
        assert!(
            map.peptide_to_proteins
                .contains_key(&normalize_il("NONEQPEP")),
            "PSM with no q-value should be kept"
        );
    }

    #[test]
    fn test_no_target_proteins() {
        let psms = vec![
            make_psm("DECOYONE", &["REV_P001"], 10.0, None, true),
            make_psm("DECOYTWO", &["REV_P002"], 8.0, None, true),
        ];
        let result = build_peptide_protein_map(&psms, None);
        assert!(
            matches!(result, Err(InferenceError::NoTargetProteins)),
            "expected NoTargetProteins when all proteins are decoy"
        );
    }
}
