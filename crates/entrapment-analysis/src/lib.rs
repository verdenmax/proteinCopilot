//! Entrapment analysis - classify trap-database PSM hits by homology to target proteome.
//!
//! Provides L0-L4 discriminability levels for each trap PSM, identifying
//! razor attribution errors, L/I isomers, near-identical homologs, and true trap hits.

pub mod config;
pub mod digest;
pub mod error;
pub mod levenshtein;
pub mod loader;
pub mod mirror_plot;
pub mod output;
pub mod mod_parser;
pub mod provenance;
pub mod report;
pub mod similarity;
pub mod tagger;
pub mod types;

pub use error::EntrapmentError;
pub use types::{
    ClassifiedPsm, DiscriminabilityLevel, EntrapmentSummary, LevelCounts, PsmGroup,
    SubstitutionType, UnifiedPsm,
};

use std::path::Path;

use config::EntrapmentConfig;
use digest::TargetDigestIndex;
use similarity::classify_single;
use tagger::Tagger;

// ---------------------------------------------------------------------------
// EntrapmentAnalyzer – public API
// ---------------------------------------------------------------------------

/// High-level API that ties together configuration, tagging, digest index,
/// and similarity classification into a single entry point.
pub struct EntrapmentAnalyzer {
    config: EntrapmentConfig,
    tagger: Tagger,
    index: TargetDigestIndex,
}

impl EntrapmentAnalyzer {
    /// Build a new analyser from a config and a target FASTA database.
    ///
    /// Compiles all tagging rules and builds the in-silico tryptic digest
    /// index from the FASTA file.
    pub fn new(config: EntrapmentConfig, fasta_path: &Path) -> Result<Self, EntrapmentError> {
        let tagger = Tagger::new(&config)?;
        let index = TargetDigestIndex::from_fasta(
            fasta_path,
            config.similarity.max_missed_cleavages,
            config.similarity.max_mismatches,
        )?;
        Ok(Self {
            config,
            tagger,
            index,
        })
    }

    /// Classify a single PSM: tag its group and determine similarity level.
    pub fn classify(&self, psm: &UnifiedPsm) -> Result<ClassifiedPsm, EntrapmentError> {
        let group = self.tagger.tag(&psm.protein_ids)?;
        Ok(classify_single(
            psm,
            group,
            &self.index,
            &self.config.similarity,
        ))
    }

    /// Classify a batch of PSMs.
    pub fn classify_all(&self, psms: &[UnifiedPsm]) -> Result<Vec<ClassifiedPsm>, EntrapmentError> {
        psms.iter().map(|psm| self.classify(psm)).collect()
    }

    /// Compute summary statistics from a set of classified PSMs.
    pub fn summary(&self, classified: &[ClassifiedPsm]) -> EntrapmentSummary {
        let mut level_counts = LevelCounts::default();
        let mut target_psms = 0usize;
        let mut trap_psms = 0usize;
        let mut ambiguous_psms = 0usize;

        // Track L0 razor families
        let mut razor_families: std::collections::HashMap<String, (usize, String, String, String)> =
            std::collections::HashMap::new();

        for cp in classified {
            match cp.group {
                PsmGroup::Target => target_psms += 1,
                PsmGroup::Trap => {
                    trap_psms += 1;
                    level_counts.increment(cp.level);

                    if cp.level == DiscriminabilityLevel::L0 {
                        let family = extract_family_name(&cp.psm.protein_ids);
                        let entry = razor_families.entry(family.clone()).or_insert_with(|| {
                            (
                                0,
                                cp.psm.peptide.clone(),
                                cp.psm.protein_ids.clone(),
                                cp.best_target_protein.clone().unwrap_or_default(),
                            )
                        });
                        entry.0 += 1;
                    }
                }
                PsmGroup::Ambiguous => ambiguous_psms += 1,
            }
        }

        let mut top_razor_families: Vec<types::RazorFamily> = razor_families
            .into_iter()
            .map(|(family, (count, pep, trap, target))| types::RazorFamily {
                family,
                count,
                example_peptide: pep,
                example_trap_protein: trap,
                example_target_protein: target,
            })
            .collect();
        top_razor_families.sort_by(|a, b| b.count.cmp(&a.count));
        top_razor_families.truncate(10);

        EntrapmentSummary {
            total_psms: classified.len(),
            target_psms,
            trap_psms,
            ambiguous_psms,
            level_counts,
            top_razor_families,
        }
    }
}

/// Extract the protein-family name from a UniProt-format accession.
///
/// E.g. `"sp|P12345|EF1A1_HUMAN"` → `"EF1A1"`.
/// Falls back to returning the first semicolon-delimited accession when
/// the format is not recognised.
fn extract_family_name(accession: &str) -> String {
    let first = accession.split(';').next().unwrap_or(accession);
    if let Some(bar_part) = first.split('|').nth(2) {
        if let Some(name) = bar_part.split('_').next() {
            return name.to_string();
        }
    }
    first.to_string()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_family_name_uniprot() {
        assert_eq!(extract_family_name("sp|P12345|EF1A1_HUMAN"), "EF1A1");
    }

    #[test]
    fn test_extract_family_name_multi_accession() {
        assert_eq!(
            extract_family_name("sp|P12345|EF1A1_HUMAN;sp|Q99999|EF1A2_HUMAN"),
            "EF1A1"
        );
    }

    #[test]
    fn test_extract_family_name_non_uniprot() {
        assert_eq!(extract_family_name("SOME_RANDOM_ID"), "SOME_RANDOM_ID");
    }
}
