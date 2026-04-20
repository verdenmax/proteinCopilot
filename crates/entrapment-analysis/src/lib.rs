//! Entrapment analysis - classify trap-database PSM hits by homology to target proteome.
//!
//! Provides L0-L4 discriminability levels for each trap PSM, identifying
//! razor attribution errors, L/I isomers, near-identical homologs, and true trap hits.

pub mod config;
pub mod error;

pub use error::EntrapmentError;
