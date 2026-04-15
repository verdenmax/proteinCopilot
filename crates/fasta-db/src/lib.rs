//! FASTA database registry, download, and cache management.
//!
//! Provides a built-in registry of common proteomics databases (UniProt Swiss-Prot
//! for major model organisms + cRAP contaminants), HTTPS download with streaming,
//! and local file caching with metadata tracking.

pub mod cache;
pub mod download;
pub mod error;
pub mod registry;

pub use error::FastaDbError;
