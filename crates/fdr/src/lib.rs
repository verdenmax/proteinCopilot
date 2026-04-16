//! # ProteinCopilot FDR
//!
//! Native target-decoy FDR calculation for proteomics search results.
//!
//! Provides:
//! - Decoy database generation (reverse/shuffle protein sequences)
//! - Target-decoy approach FDR calculation
//! - q-value assignment with monotonicity enforcement

pub mod calculation;
pub mod decoy;
pub mod error;
pub mod peptide_fdr;
pub mod protein_fdr;

pub use calculation::calculate_fdr;
pub use decoy::{generate_decoys, DecoyProtein};
pub use error::FdrError;
pub use peptide_fdr::{calculate_peptide_fdr, extract_unique_peptides, PeptideFdrResult, PeptideScore};
pub use protein_fdr::{calculate_protein_fdr, ProteinFdrResult};
