//! Target FASTA digest index for entrapment classification.
//!
//! Wraps the search-engine crate's FASTA parsing and in-silico digestion to build
//! efficient lookup structures used by the L0/L1/L2–L4 classification pipeline.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use tracing::info;

use protein_copilot_core::search_params::Enzyme;
use protein_copilot_search_engine::digest::{digest_with_length, is_standard_sequence};
use protein_copilot_search_engine::fasta::parse_fasta;

use crate::error::EntrapmentError;

/// A single digested target peptide with its source protein.
#[derive(Debug, Clone)]
pub struct TargetPeptide {
    /// Amino acid sequence.
    pub sequence: String,
    /// Protein accession this peptide came from.
    pub protein_accession: String,
    /// Monoisotopic neutral mass (Da).
    pub neutral_mass: f64,
}

/// Replace all 'I' (isoleucine) with 'L' (leucine) for L/I-normalized comparison.
///
/// Both residues have the same monoisotopic mass (113.084064 Da) and are
/// indistinguishable by mass spectrometry.
pub fn normalize_li(seq: &str) -> String {
    seq.replace('I', "L")
}

/// Lookup index built from in-silico tryptic digestion of a target FASTA database.
///
/// Supports exact (L0), L/I-normalized (L1), and length-based (L2–L4) lookups.
#[derive(Debug)]
pub struct TargetDigestIndex {
    /// Peptides grouped by sequence length for hamming-distance scanning.
    pub by_length: HashMap<usize, Vec<TargetPeptide>>,
    /// All unique target peptide sequences for exact (L0) lookup.
    pub exact_set: HashSet<String>,
    /// L/I-normalized sequences for L1 lookup.
    pub normalized_set: HashSet<String>,
    /// Sequence → protein accession for exact (L0) matches.
    pub exact_to_protein: HashMap<String, String>,
    /// Normalized sequence → (original sequence, protein accession) for L1 matches.
    pub normalized_to_original: HashMap<String, (String, String)>,
}

impl TargetDigestIndex {
    /// Build the digest index from a FASTA file.
    ///
    /// Parses the FASTA, performs tryptic digestion with the given number of
    /// missed cleavages (length range 6–50), filters to standard-residue
    /// peptides, and constructs all lookup structures.
    pub fn from_fasta(path: &Path, max_missed_cleavages: u32) -> Result<Self, EntrapmentError> {
        let entries = parse_fasta(path).map_err(|e| EntrapmentError::FastaError {
            path: path.to_path_buf(),
            detail: e.to_string(),
        })?;

        let mut exact_set = HashSet::new();
        let mut normalized_set = HashSet::new();
        let mut exact_to_protein = HashMap::new();
        let mut normalized_to_original: HashMap<String, (String, String)> = HashMap::new();
        let mut by_length: HashMap<usize, Vec<TargetPeptide>> = HashMap::new();

        let enzyme = Enzyme::Trypsin;

        for entry in &entries {
            let peptides = digest_with_length(
                &entry.sequence,
                &entry.accession,
                &enzyme,
                max_missed_cleavages,
                6,
                50,
            );

            for dp in peptides {
                if !is_standard_sequence(&dp.sequence) {
                    continue;
                }

                let seq = &dp.sequence;

                // Build exact structures (first protein wins for duplicates)
                if exact_set.insert(seq.clone()) {
                    exact_to_protein.insert(seq.clone(), dp.protein_accession.clone());
                }

                // Build normalized structures
                let norm = normalize_li(seq);
                if !normalized_set.contains(&norm) {
                    normalized_set.insert(norm.clone());
                    normalized_to_original
                        .insert(norm, (seq.clone(), dp.protein_accession.clone()));
                }

                // Build by-length index
                let target_peptide = TargetPeptide {
                    sequence: seq.clone(),
                    protein_accession: dp.protein_accession.clone(),
                    neutral_mass: dp.neutral_mass,
                };
                by_length
                    .entry(seq.len())
                    .or_default()
                    .push(target_peptide);
            }
        }

        info!(
            path = %path.display(),
            proteins = entries.len(),
            unique_peptides = exact_set.len(),
            "built target digest index"
        );

        Ok(Self {
            by_length,
            exact_set,
            normalized_set,
            exact_to_protein,
            normalized_to_original,
        })
    }

    /// Check whether the sequence is an exact match to a target peptide (L0).
    pub fn has_exact(&self, seq: &str) -> bool {
        self.exact_set.contains(seq)
    }

    /// Return the protein accession for an exact-match sequence, if present.
    pub fn exact_protein(&self, seq: &str) -> Option<&str> {
        self.exact_to_protein.get(seq).map(|s| s.as_str())
    }

    /// Check whether the L/I-normalized sequence matches a target peptide (L1).
    pub fn has_normalized(&self, seq: &str) -> bool {
        let norm = normalize_li(seq);
        self.normalized_set.contains(&norm)
    }

    /// Return the original target sequence and protein accession for an L1 match.
    pub fn normalized_match(&self, seq: &str) -> Option<(&str, &str)> {
        let norm = normalize_li(seq);
        self.normalized_to_original
            .get(&norm)
            .map(|(orig, prot)| (orig.as_str(), prot.as_str()))
    }

    /// Return all target peptides of a given sequence length (for hamming scanning).
    ///
    /// Returns an empty slice if no peptides of this length exist.
    pub fn peptides_of_length(&self, len: usize) -> &[TargetPeptide] {
        self.by_length
            .get(&len)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Total number of unique target peptides in the index.
    pub fn len(&self) -> usize {
        self.exact_set.len()
    }

    /// Whether the index contains no peptides.
    pub fn is_empty(&self) -> bool {
        self.exact_set.is_empty()
    }

    /// Construct an empty index for unit tests.
    #[cfg(test)]
    pub fn empty_for_test() -> Self {
        Self {
            by_length: HashMap::new(),
            exact_set: HashSet::new(),
            normalized_set: HashSet::new(),
            exact_to_protein: HashMap::new(),
            normalized_to_original: HashMap::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalize_li() {
        // No I in the sequence → unchanged
        assert_eq!(normalize_li("ELTALAPSTMK"), "ELTALAPSTMK");
        // I replaced with L
        assert_eq!(normalize_li("HPFPGPGIAIR"), "HPFPGPGLALR");
    }

    #[test]
    fn test_normalize_li_empty() {
        assert_eq!(normalize_li(""), "");
    }

    #[test]
    fn test_empty_index() {
        let idx = TargetDigestIndex::empty_for_test();
        assert!(idx.is_empty());
        assert_eq!(idx.len(), 0);
        assert!(!idx.has_exact("PEPTIDEK"));
        assert!(!idx.has_normalized("PEPTIDEK"));
        assert!(idx.exact_protein("PEPTIDEK").is_none());
        assert!(idx.normalized_match("PEPTIDEK").is_none());
        assert!(idx.peptides_of_length(8).is_empty());
    }

    #[test]
    fn test_from_fasta_roundtrip() {
        use std::io::Write;

        let mut f = tempfile::NamedTempFile::new().expect("create temp file");
        write!(
            f,
            ">sp|P001|TEST_HUMAN Test protein\n\
             PEPTIDEKANSTHERPEPTIDERLASTPART\n"
        )
        .expect("write fasta");

        let idx = TargetDigestIndex::from_fasta(f.path(), 0).expect("build index");

        // Should contain the tryptic peptides >= 6 aa
        assert!(idx.has_exact("PEPTIDEK"));
        assert_eq!(idx.exact_protein("PEPTIDEK"), Some("sp|P001|TEST_HUMAN"));

        // L/I-normalised lookup: "PEPTIDEK" normalises to "PEPTLDEK"
        // Same sequence should also match via normalised path
        assert!(idx.has_normalized("PEPTIDEK"));

        // length-based lookup
        assert!(!idx.peptides_of_length(8).is_empty()); // "PEPTIDEK" is 8 chars

        assert!(!idx.is_empty());
        assert!(idx.len() > 0);
    }

    #[test]
    fn test_normalized_li_lookup() {
        use std::io::Write;

        // Create a FASTA with a peptide containing 'I'
        let mut f = tempfile::NamedTempFile::new().expect("create temp file");
        // "ELVISLSK" is 8 chars, trypsin cuts after K at end
        write!(
            f,
            ">sp|P002|LI_TEST LI test protein\n\
             ELVISLSKAGATHER\n"
        )
        .expect("write fasta");

        let idx = TargetDigestIndex::from_fasta(f.path(), 0).expect("build index");

        // "ELVISLSK" should be in the index (8 chars, ends with K)
        assert!(idx.has_exact("ELVISLSK"));

        // Normalised: "ELVISLSK" → "ELVLSLSK"
        // A query "ELVLSLSK" (with L instead of I) should match via L1
        assert!(idx.has_normalized("ELVLSLSK"));
        let (orig, _prot) = idx.normalized_match("ELVLSLSK").expect("should find L1 match");
        assert_eq!(orig, "ELVISLSK");
    }

    #[test]
    fn test_from_fasta_error_on_missing_file() {
        let result = TargetDigestIndex::from_fasta(Path::new("/nonexistent/db.fasta"), 0);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(
            err.to_string().contains("FASTA error"),
            "unexpected error: {err}"
        );
    }
}
