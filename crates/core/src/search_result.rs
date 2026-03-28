//! Standardized search result structures across all search engines.
//!
//! This module defines the types for representing proteomics search results
//! at three levels of granularity:
//! - [`Psm`] — Peptide-Spectrum Match (individual spectrum-to-peptide assignment)
//! - [`PeptideResult`] — Peptide-level aggregation
//! - [`ProteinResult`] — Protein-level aggregation
//! - [`SearchResultSummary`] — Statistical overview for LLM consumption
//! - [`SearchResult`] — Complete search output combining all levels
//!
//! All search engines (pFind, MSFragger, Comet, etc.) normalize their
//! output into these shared types via the adapter layer.

use std::collections::HashMap;

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use uuid::Uuid;

use crate::search_params::Modification;

// ---------------------------------------------------------------------------
// PSM (Peptide-Spectrum Match)
// ---------------------------------------------------------------------------

/// A single Peptide-Spectrum Match — the fundamental unit of a database search.
///
/// Each PSM represents the assignment of a peptide sequence to a specific
/// MS2 spectrum, along with scoring and quality metrics.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct Psm {
    /// Scan number of the matched spectrum (1-based).
    pub spectrum_scan: u32,
    /// Identified peptide sequence (one-letter amino acid codes).
    pub peptide_sequence: String,
    /// Modifications identified on this peptide.
    pub modifications: Vec<Modification>,
    /// Precursor charge state.
    pub charge: i32,
    /// Observed precursor m/z (Da).
    pub precursor_mz: f64,
    /// Theoretical precursor m/z calculated from the peptide (Da).
    pub calculated_mz: f64,
    /// Mass deviation between observed and calculated (ppm).
    pub delta_mass_ppm: f64,
    /// Search engine score (higher = better match, engine-dependent).
    pub score: f64,
    /// False discovery rate q-value (`None` if FDR not yet calculated).
    pub q_value: Option<f64>,
    /// Protein accessions this peptide maps to.
    pub protein_accessions: Vec<String>,
    /// Whether this PSM is from the decoy database.
    pub is_decoy: bool,
}

// ---------------------------------------------------------------------------
// PeptideResult
// ---------------------------------------------------------------------------

/// Peptide-level search result, aggregated from PSMs.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct PeptideResult {
    /// Peptide amino acid sequence.
    pub sequence: String,
    /// Protein accessions containing this peptide.
    pub protein_accessions: Vec<String>,
    /// Best score among all PSMs for this peptide.
    pub best_score: f64,
    /// q-value at peptide level (`None` if not calculated).
    pub q_value: Option<f64>,
    /// Number of PSMs supporting this peptide.
    pub psm_count: u64,
}

// ---------------------------------------------------------------------------
// ProteinResult
// ---------------------------------------------------------------------------

/// Protein-level search result, aggregated from peptides.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct ProteinResult {
    /// Protein accession (e.g. UniProt ID like "P12345").
    pub accession: String,
    /// Protein description / name.
    pub description: String,
    /// Sequence coverage (0.0–1.0).
    pub coverage: f64,
    /// Total number of peptides mapped to this protein.
    pub peptide_count: u64,
    /// Number of unique (non-shared) peptides.
    pub unique_peptide_count: u64,
}

// ---------------------------------------------------------------------------
// SearchResultSummary
// ---------------------------------------------------------------------------

/// Statistical summary of search results for LLM-driven interpretation.
///
/// This is the primary input for the AI layer to understand and explain
/// search quality. All fields are deterministically computed by Rust.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SearchResultSummary {
    /// Total number of spectra submitted to the search.
    pub total_spectra_searched: u64,
    /// Total number of PSMs returned by the engine (before FDR filtering).
    pub total_psms: u64,
    /// PSMs passing 1% FDR threshold.
    pub psms_at_1pct_fdr: u64,
    /// Unique peptide sequences at 1% FDR.
    pub unique_peptides_at_1pct_fdr: u64,
    /// Protein groups at 1% FDR.
    pub protein_groups_at_1pct_fdr: u64,
    /// Median search engine score across all PSMs.
    pub median_score: f64,
    /// Median mass deviation (ppm) across all PSMs.
    pub median_delta_mass_ppm: f64,
    /// Identification rate: `psms_at_1pct_fdr / total_spectra_searched`.
    pub identification_rate: f64,
    /// Modification frequency distribution (modification name → count).
    pub modification_distribution: HashMap<String, u64>,
    /// Charge state distribution (charge → count).
    pub charge_distribution: HashMap<i32, u64>,
    /// Total search duration in seconds.
    pub search_duration_sec: f64,
}

// ---------------------------------------------------------------------------
// SearchResult
// ---------------------------------------------------------------------------

/// Complete search result from one engine run.
///
/// This is the top-level output of a search engine invocation, combining
/// all PSMs, peptide/protein rollups, statistical summary, and metadata.
/// The `run_id` links this result to the originating run context.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SearchResult {
    /// Unique identifier for this analysis run.
    pub run_id: Uuid,
    /// Name and version of the search engine that produced this result.
    pub engine_name: String,
    /// Version of the search engine.
    pub engine_version: String,
    /// All PSMs returned by the search.
    pub psms: Vec<Psm>,
    /// Peptide-level aggregation.
    pub peptides: Vec<PeptideResult>,
    /// Protein-level aggregation.
    pub proteins: Vec<ProteinResult>,
    /// Statistical summary.
    pub summary: SearchResultSummary,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::search_params::{ModPosition, Modification};

    fn sample_psm() -> Psm {
        Psm {
            spectrum_scan: 42,
            peptide_sequence: "PEPTIDER".to_string(),
            modifications: vec![],
            charge: 2,
            precursor_mz: 471.2561,
            calculated_mz: 471.2555,
            delta_mass_ppm: 1.27,
            score: 45.6,
            q_value: Some(0.001),
            protein_accessions: vec!["P12345".to_string()],
            is_decoy: false,
        }
    }

    fn sample_psm_with_mod() -> Psm {
        Psm {
            spectrum_scan: 100,
            peptide_sequence: "MYPEPTIDE".to_string(),
            modifications: vec![Modification {
                name: "Oxidation".to_string(),
                mass_delta: 15.994915,
                residues: vec!['M'],
                position: ModPosition::Anywhere,
            }],
            charge: 3,
            precursor_mz: 362.1789,
            calculated_mz: 362.1785,
            delta_mass_ppm: 1.1,
            score: 38.2,
            q_value: Some(0.005),
            protein_accessions: vec!["P12345".to_string(), "Q67890".to_string()],
            is_decoy: false,
        }
    }

    fn sample_peptide_result() -> PeptideResult {
        PeptideResult {
            sequence: "PEPTIDER".to_string(),
            protein_accessions: vec!["P12345".to_string()],
            best_score: 45.6,
            q_value: Some(0.001),
            psm_count: 3,
        }
    }

    fn sample_protein_result() -> ProteinResult {
        ProteinResult {
            accession: "P12345".to_string(),
            description: "Serum albumin OS=Homo sapiens".to_string(),
            coverage: 0.35,
            peptide_count: 12,
            unique_peptide_count: 10,
        }
    }

    fn sample_summary() -> SearchResultSummary {
        let mut mod_dist = HashMap::new();
        mod_dist.insert("Oxidation".to_string(), 1200);
        mod_dist.insert("Carbamidomethyl".to_string(), 8500);

        let mut charge_dist = HashMap::new();
        charge_dist.insert(2, 15000);
        charge_dist.insert(3, 10000);

        SearchResultSummary {
            total_spectra_searched: 83000,
            total_psms: 35000,
            psms_at_1pct_fdr: 28000,
            unique_peptides_at_1pct_fdr: 12000,
            protein_groups_at_1pct_fdr: 3500,
            median_score: 32.5,
            median_delta_mass_ppm: 0.8,
            identification_rate: 0.337,
            modification_distribution: mod_dist,
            charge_distribution: charge_dist,
            search_duration_sec: 245.0,
        }
    }

    fn sample_search_result() -> SearchResult {
        SearchResult {
            run_id: Uuid::nil(),
            engine_name: "pFind".to_string(),
            engine_version: "3.1.0".to_string(),
            psms: vec![sample_psm(), sample_psm_with_mod()],
            peptides: vec![sample_peptide_result()],
            proteins: vec![sample_protein_result()],
            summary: sample_summary(),
        }
    }

    // -- PSM ------------------------------------------------------------

    #[test]
    fn psm_serde_roundtrip() {
        let psm = sample_psm();
        let json = serde_json::to_string_pretty(&psm).unwrap();
        let back: Psm = serde_json::from_str(&json).unwrap();
        assert_eq!(psm, back);
    }

    #[test]
    fn psm_with_modifications_roundtrip() {
        let psm = sample_psm_with_mod();
        let json = serde_json::to_string_pretty(&psm).unwrap();
        let back: Psm = serde_json::from_str(&json).unwrap();
        assert_eq!(psm, back);
    }

    #[test]
    fn psm_decoy_flag() {
        let mut psm = sample_psm();
        psm.is_decoy = true;
        psm.protein_accessions = vec!["REV_P12345".to_string()];
        let json = serde_json::to_string(&psm).unwrap();
        let back: Psm = serde_json::from_str(&json).unwrap();
        assert!(back.is_decoy);
    }

    #[test]
    fn psm_no_q_value() {
        let mut psm = sample_psm();
        psm.q_value = None;
        let json = serde_json::to_string(&psm).unwrap();
        let back: Psm = serde_json::from_str(&json).unwrap();
        assert!(back.q_value.is_none());
    }

    // -- PeptideResult --------------------------------------------------

    #[test]
    fn peptide_result_serde_roundtrip() {
        let pr = sample_peptide_result();
        let json = serde_json::to_string_pretty(&pr).unwrap();
        let back: PeptideResult = serde_json::from_str(&json).unwrap();
        assert_eq!(pr, back);
    }

    // -- ProteinResult --------------------------------------------------

    #[test]
    fn protein_result_serde_roundtrip() {
        let pr = sample_protein_result();
        let json = serde_json::to_string_pretty(&pr).unwrap();
        let back: ProteinResult = serde_json::from_str(&json).unwrap();
        assert_eq!(pr, back);
    }

    // -- SearchResultSummary --------------------------------------------

    #[test]
    fn summary_serde_roundtrip() {
        let s = sample_summary();
        let json = serde_json::to_string_pretty(&s).unwrap();
        let back: SearchResultSummary = serde_json::from_str(&json).unwrap();
        assert_eq!(s.total_spectra_searched, back.total_spectra_searched);
        assert_eq!(s.psms_at_1pct_fdr, back.psms_at_1pct_fdr);
        assert_eq!(
            s.modification_distribution.len(),
            back.modification_distribution.len()
        );
        assert_eq!(s.charge_distribution.len(), back.charge_distribution.len());
    }

    // -- SearchResult ---------------------------------------------------

    #[test]
    fn search_result_full_roundtrip() {
        let result = sample_search_result();
        let json = serde_json::to_string_pretty(&result).unwrap();
        let back: SearchResult = serde_json::from_str(&json).unwrap();
        assert_eq!(result.run_id, back.run_id);
        assert_eq!(result.engine_name, back.engine_name);
        assert_eq!(result.psms.len(), back.psms.len());
        assert_eq!(result.peptides.len(), back.peptides.len());
        assert_eq!(result.proteins.len(), back.proteins.len());
        assert_eq!(
            result.summary.total_spectra_searched,
            back.summary.total_spectra_searched
        );
    }

    #[test]
    fn search_result_empty_collections() {
        let result = SearchResult {
            run_id: Uuid::nil(),
            engine_name: "test".to_string(),
            engine_version: "1.0".to_string(),
            psms: vec![],
            peptides: vec![],
            proteins: vec![],
            summary: SearchResultSummary {
                total_spectra_searched: 0,
                total_psms: 0,
                psms_at_1pct_fdr: 0,
                unique_peptides_at_1pct_fdr: 0,
                protein_groups_at_1pct_fdr: 0,
                median_score: 0.0,
                median_delta_mass_ppm: 0.0,
                identification_rate: 0.0,
                modification_distribution: HashMap::new(),
                charge_distribution: HashMap::new(),
                search_duration_sec: 0.0,
            },
        };
        let json = serde_json::to_string(&result).unwrap();
        let back: SearchResult = serde_json::from_str(&json).unwrap();
        assert!(back.psms.is_empty());
        assert!(back.peptides.is_empty());
        assert!(back.proteins.is_empty());
    }
}
