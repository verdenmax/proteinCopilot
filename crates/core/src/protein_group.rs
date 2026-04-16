//! Protein group and inference result data structures.
//!
//! Used by the protein-inference crate and exported via MCP tools.

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// A group of indistinguishable or related proteins identified by shared peptide evidence.
///
/// In proteomics, proteins sharing the same set of identified peptides cannot be
/// distinguished and are grouped together. The "leader" is the representative protein
/// (typically the one with the most unique peptides or the canonical accession).
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct ProteinGroup {
    /// Representative protein accession (the "leader" of this group).
    pub leader_accession: String,

    /// Description of the leader protein.
    pub leader_description: String,

    /// All protein accessions in this group (including the leader).
    /// Proteins with identical peptide sets are grouped together.
    pub member_accessions: Vec<String>,

    /// All peptide sequences mapped to this group (unique + shared).
    pub peptides: Vec<String>,

    /// Peptide sequences unique to this group (not shared with other groups).
    pub unique_peptides: Vec<String>,

    /// Razor peptides assigned to this group (shared peptides allocated by Razor logic).
    pub razor_peptides: Vec<String>,

    /// Group-level score (e.g., best peptide score or sum of top-3).
    pub score: f64,

    /// Protein-level q-value (`None` if not yet calculated).
    pub q_value: Option<f64>,

    /// Sequence coverage (0.0–1.0). `None` if not calculated.
    pub coverage: Option<f64>,

    /// Whether this group is from the decoy database.
    pub is_decoy: bool,
}

/// Complete result of protein inference.
///
/// Produced by running the parsimony algorithm, razor peptide assignment,
/// and protein-level FDR on a set of PSMs.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct InferenceResult {
    /// Protein groups after parsimony and FDR filtering.
    pub groups: Vec<ProteinGroup>,

    /// Razor peptide mapping: peptide_sequence → assigned leader_accession.
    /// Shared peptides are assigned to the protein group with the most unique peptides.
    pub razor_map: HashMap<String, String>,

    /// Total number of target protein groups (before FDR filtering).
    pub total_target_groups: u64,

    /// Total number of decoy protein groups.
    pub total_decoy_groups: u64,

    /// Number of protein groups passing the FDR threshold (typically 1%).
    pub groups_at_1pct_fdr: u64,

    /// Number of unique peptides used for inference.
    pub unique_peptides_used: u64,
}

impl ProteinGroup {
    /// Validates the protein group for correctness.
    pub fn validate(&self) -> Result<(), ProteinGroupError> {
        if self.leader_accession.trim().is_empty() {
            return Err(ProteinGroupError::EmptyLeader);
        }
        if self.member_accessions.is_empty() {
            return Err(ProteinGroupError::NoMembers);
        }
        if !self.member_accessions.contains(&self.leader_accession) {
            return Err(ProteinGroupError::LeaderNotInMembers);
        }
        if self.peptides.is_empty() {
            return Err(ProteinGroupError::NoPeptides);
        }
        if let Some(cov) = self.coverage {
            if !(0.0..=1.0).contains(&cov) {
                return Err(ProteinGroupError::InvalidCoverage(cov));
            }
        }
        if let Some(q) = self.q_value {
            if !(0.0..=1.0).contains(&q) {
                return Err(ProteinGroupError::InvalidQValue(q));
            }
        }
        Ok(())
    }
}

/// Errors in protein group validation.
#[derive(Debug, thiserror::Error)]
pub enum ProteinGroupError {
    #[error("leader_accession must not be empty")]
    EmptyLeader,
    #[error("member_accessions must not be empty")]
    NoMembers,
    #[error("leader_accession must be in member_accessions")]
    LeaderNotInMembers,
    #[error("peptides must not be empty")]
    NoPeptides,
    #[error("coverage must be in [0.0, 1.0], got {0}")]
    InvalidCoverage(f64),
    #[error("q_value must be in [0.0, 1.0], got {0}")]
    InvalidQValue(f64),
}
