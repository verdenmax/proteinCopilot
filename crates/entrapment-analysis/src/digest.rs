//! Target FASTA digest index for entrapment classification.
//!
//! Wraps the search-engine crate's FASTA parsing and in-silico digestion to build
//! efficient lookup structures used by the L0/L1/L2–L4 classification pipeline.

use std::collections::hash_map::DefaultHasher;
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::path::Path;

use tracing::info;

use protein_copilot_core::search_params::Enzyme;
use protein_copilot_search_engine::digest::{digest_with_length, is_standard_sequence};
use protein_copilot_search_engine::fasta::parse_fasta;

use crate::config::SimilarityConfig;
use crate::error::EntrapmentError;
use crate::levenshtein;
use crate::types::SubstitutionType;

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

/// A match found by the similarity search.
#[derive(Debug, Clone)]
pub struct SimilarityMatch {
    /// The matching target peptide sequence.
    pub target_peptide: String,
    /// Protein accession the target peptide belongs to.
    pub target_protein: String,
    /// Edit distance between query and target peptide.
    pub edit_distance: u32,
    /// Mass difference in Da (target − query).
    pub delta_mass_da: f64,
    /// Human-readable alignment detail, e.g. "D0→N" or "ins:G@5".
    pub alignment_detail: String,
    /// Substitution type annotation (categorized later by classify_single).
    pub substitution_type: SubstitutionType,
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

    // v2: k-mer inverted index for cross-length similarity search
    /// Inverted index: k-mer hash → list of peptide IDs in `all_peptides`.
    pub kmer_index: HashMap<u64, Vec<u32>>,
    /// Flat array of all target peptides (referenced by k-mer posting lists).
    pub all_peptides: Vec<TargetPeptide>,
    /// The k-mer length used for the index (pigeonhole guarantee).
    pub kmer_k: usize,
}

/// Hash a k-mer byte slice to u64.
fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut hasher = DefaultHasher::new();
    kmer.hash(&mut hasher);
    hasher.finish()
}

/// Extract all k-mers from a sequence and return their hashes.
fn extract_kmers(seq: &[u8], k: usize) -> Vec<u64> {
    if seq.len() < k {
        return Vec::new();
    }
    (0..=(seq.len() - k))
        .map(|i| hash_kmer(&seq[i..i + k]))
        .collect()
}

impl TargetDigestIndex {
    /// Build the digest index from a FASTA file.
    ///
    /// Parses the FASTA, performs tryptic digestion with the given number of
    /// missed cleavages (length range 6–50), filters to standard-residue
    /// peptides, and constructs all lookup structures.
    pub fn from_fasta(
        path: &Path,
        max_missed_cleavages: u32,
        max_edit_distance: u16,
    ) -> Result<Self, EntrapmentError> {
        let entries = parse_fasta(path).map_err(|e| EntrapmentError::FastaError {
            path: path.to_path_buf(),
            detail: e.to_string(),
        })?;

        if entries.is_empty() {
            info!(path = %path.display(), "FASTA file contains no protein sequences");
        }

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
                by_length.entry(seq.len()).or_default().push(target_peptide);
            }
        }

        info!(
            path = %path.display(),
            proteins = entries.len(),
            unique_peptides = exact_set.len(),
            "built target digest index"
        );

        // Determine k for pigeonhole guarantee: k = min_len / (max_edit + 1)
        // Must use the actual configured max_edit_distance to ensure no false negatives.
        let min_peptide_len = 6usize; // our min digest length
        let kmer_k = (min_peptide_len / (max_edit_distance as usize + 1)).max(1);

        // Build flat peptide array and k-mer inverted index
        let mut all_peptides: Vec<TargetPeptide> = Vec::new();
        let mut kmer_index: HashMap<u64, Vec<u32>> = HashMap::new();

        for peptides in by_length.values() {
            for tp in peptides {
                let pid = all_peptides.len() as u32;
                let kmers = extract_kmers(tp.sequence.as_bytes(), kmer_k);
                for kh in kmers {
                    kmer_index.entry(kh).or_default().push(pid);
                }
                all_peptides.push(tp.clone());
            }
        }

        // Deduplicate each posting list
        for list in kmer_index.values_mut() {
            list.sort_unstable();
            list.dedup();
        }

        Ok(Self {
            by_length,
            exact_set,
            normalized_set,
            exact_to_protein,
            normalized_to_original,
            kmer_index,
            all_peptides,
            kmer_k,
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

    /// Find target peptides within `max_edit_dist` edit distance of `query`.
    ///
    /// Uses k-mer pre-filtering (pigeonhole principle) to reduce candidates,
    /// then computes full Levenshtein distance + alignment for survivors.
    pub fn find_similar(
        &self,
        query: &str,
        max_edit_dist: u16,
        len_tolerance: usize,
        _config: &SimilarityConfig,
    ) -> Vec<SimilarityMatch> {
        let query_len = query.len();
        let min_len = query_len.saturating_sub(len_tolerance);
        let max_len = query_len + len_tolerance;

        // Extract k-mers from query
        let query_kmers = extract_kmers(query.as_bytes(), self.kmer_k);
        if query_kmers.is_empty() {
            return Vec::new();
        }

        // Collect candidate peptide IDs from k-mer hits
        let mut candidate_ids: HashSet<u32> = HashSet::new();
        for kh in &query_kmers {
            if let Some(ids) = self.kmer_index.get(kh) {
                for &id in ids {
                    candidate_ids.insert(id);
                }
            }
        }

        // Filter and compute edit distance
        let mut results = Vec::new();
        for &pid in &candidate_ids {
            debug_assert!(
                (pid as usize) < self.all_peptides.len(),
                "kmer_index contains out-of-bounds peptide ID {pid}"
            );
            let tp = &self.all_peptides[pid as usize];

            // Length filter
            if tp.sequence.len() < min_len || tp.sequence.len() > max_len {
                continue;
            }

            // Skip exact matches (handled by L0 path)
            if tp.sequence == query {
                continue;
            }

            // Quick edit distance check
            let dist = levenshtein::edit_distance(query, &tp.sequence);
            if dist > max_edit_dist as u32 {
                continue;
            }

            // Full alignment for survivors
            let alignment = levenshtein::align(query, &tp.sequence);

            results.push(SimilarityMatch {
                target_peptide: tp.sequence.clone(),
                target_protein: tp.protein_accession.clone(),
                edit_distance: alignment.edit_distance,
                delta_mass_da: alignment.delta_mass_da,
                alignment_detail: alignment.alignment_detail,
                substitution_type: SubstitutionType::None, // categorized later by classify_single
            });
        }

        // Sort by edit distance, then by |delta_mass|
        results.sort_by(|a, b| {
            a.edit_distance.cmp(&b.edit_distance).then_with(|| {
                a.delta_mass_da
                    .abs()
                    .partial_cmp(&b.delta_mass_da.abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
        });

        results
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
            kmer_index: HashMap::new(),
            all_peptides: Vec::new(),
            kmer_k: 2,
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

        let idx = TargetDigestIndex::from_fasta(f.path(), 0, 2).expect("build index");

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

        let idx = TargetDigestIndex::from_fasta(f.path(), 0, 2).expect("build index");

        // "ELVISLSK" should be in the index (8 chars, ends with K)
        assert!(idx.has_exact("ELVISLSK"));

        // Normalised: "ELVISLSK" → "ELVLSLSK"
        // A query "ELVLSLSK" (with L instead of I) should match via L1
        assert!(idx.has_normalized("ELVLSLSK"));
        let (orig, _prot) = idx
            .normalized_match("ELVLSLSK")
            .expect("should find L1 match");
        assert_eq!(orig, "ELVISLSK");
    }

    #[test]
    fn test_from_fasta_error_on_missing_file() {
        let result = TargetDigestIndex::from_fasta(Path::new("/nonexistent/db.fasta"), 0, 2);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(
            err.to_string().contains("FASTA error"),
            "unexpected error: {err}"
        );
    }

    // --- k-mer index and find_similar() tests (Task 4) ---

    use crate::config::SimilarityConfig;

    #[test]
    fn test_kmer_index_built_on_from_fasta() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().expect("create temp file");
        write!(
            f,
            ">sp|P001|TEST_HUMAN Test\nPEPTIDEKANSTHERPEPTIDERLASTPART\n"
        )
        .unwrap();
        let idx = TargetDigestIndex::from_fasta(f.path(), 2, 2).expect("build index");
        // k-mer index should be populated
        assert!(!idx.kmer_index.is_empty());
        assert!(!idx.all_peptides.is_empty());
        assert!(idx.kmer_k >= 2);
    }

    #[test]
    fn test_find_similar_exact_match_excluded() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().expect("create temp file");
        write!(f, ">sp|P001|TEST_HUMAN Test\nPEPTIDEKANSTHERR\n").unwrap();
        let config = SimilarityConfig::default();
        let idx = TargetDigestIndex::from_fasta(
            f.path(),
            config.max_missed_cleavages,
            config.max_mismatches,
        )
        .unwrap();
        // "PEPTIDEK" is in the index; find_similar should not return exact matches
        let matches = idx.find_similar("PEPTIDEK", 2, 2, &config);
        for m in &matches {
            assert_ne!(m.target_peptide, "PEPTIDEK");
        }
    }

    #[test]
    fn test_find_similar_one_substitution() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().expect("create temp file");
        // Target has "NGFLLDGFPR", query is "DGFLLDGFPR" (D→N, edit=1)
        write!(f, ">sp|P001|TEST_HUMAN Test\nNGFLLDGFPR\n").unwrap();
        let config = SimilarityConfig::default();
        let idx = TargetDigestIndex::from_fasta(
            f.path(),
            config.max_missed_cleavages,
            config.max_mismatches,
        )
        .unwrap();
        let matches = idx.find_similar("DGFLLDGFPR", 2, 2, &config);
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].target_peptide, "NGFLLDGFPR");
        assert_eq!(matches[0].edit_distance, 1);
        assert!(
            (matches[0].delta_mass_da.abs() - 0.984016).abs() < 0.001,
            "unexpected delta_mass_da: {}",
            matches[0].delta_mass_da
        );
    }

    #[test]
    fn test_find_similar_different_length_indel() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().expect("create temp file");
        // Target: "PEPGGDEK" (8 chars), Query: "PEPNDEK" (7 chars)
        write!(f, ">sp|P001|TEST_HUMAN Test\nPEPGGDEKKAAAAAR\n").unwrap();
        let config = SimilarityConfig::default();
        let idx = TargetDigestIndex::from_fasta(
            f.path(),
            config.max_missed_cleavages,
            config.max_mismatches,
        )
        .unwrap();
        let matches = idx.find_similar("PEPNDEK", 2, 2, &config);
        // Should find "PEPGGDEK" as a candidate (len 8 is within len_tolerance=2 of len 7)
        let found = matches.iter().any(|m| m.target_peptide == "PEPGGDEK");
        assert!(
            found,
            "should find PEPGGDEK as similar to PEPNDEK; found: {:?}",
            matches
        );
    }

    #[test]
    fn test_find_similar_no_match_beyond_tolerance() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().expect("create temp file");
        write!(f, ">sp|P001|TEST_HUMAN Test\nWWWWWWWWR\n").unwrap();
        let config = SimilarityConfig::default();
        let idx = TargetDigestIndex::from_fasta(
            f.path(),
            config.max_missed_cleavages,
            config.max_mismatches,
        )
        .unwrap();
        // Completely different peptide
        let matches = idx.find_similar("PEPTIDEK", 2, 2, &config);
        assert!(matches.is_empty());
    }
}
