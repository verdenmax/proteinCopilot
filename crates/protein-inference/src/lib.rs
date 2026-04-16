//! # ProteinCopilot Protein Inference
//!
//! Implements protein inference from peptide-spectrum matches:
//! - Peptide-to-protein mapping
//! - Parsimony (greedy set cover) algorithm
//! - Razor peptide assignment
//! - Sequence coverage calculation
//! - Protein-level FDR (picked-protein approach)

pub mod coverage;
pub mod error;
pub mod mapper;
pub mod parsimony;
pub mod razor;
