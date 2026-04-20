//! Entrapment analysis - classify trap-database PSM hits by homology to target proteome.
//!
//! Provides L0-L4 discriminability levels for each trap PSM, identifying
//! razor attribution errors, L/I isomers, near-identical homologs, and true trap hits.

pub mod config;
pub mod digest;
pub mod error;
pub mod loader;
pub mod similarity;
pub mod tagger;
pub mod types;

pub use error::EntrapmentError;
pub use types::{
    ClassifiedPsm, DiscriminabilityLevel, EntrapmentSummary, LevelCounts, PsmGroup, UnifiedPsm,
};
