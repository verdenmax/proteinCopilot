//! # ProteinCopilot Search Engine
//!
//! Library crate for search engine orchestration and execution.
//! Contains a simplified built-in search engine for MVP validation,
//! with architecture ready for pFind/MSFragger/Comet adapters.
//!
//! # Architecture
//!
//! ```text
//! SearchParams + input_files
//!        │
//!        ▼
//!  EngineRegistry.get("SimpleSearch")
//!        │
//!        ▼
//!  SearchEngineAdapter::search()
//!        │
//!        ├── Read FASTA → digest proteins → theoretical peptides
//!        ├── Match spectra precursors → candidate peptides
//!        ├── Score matches (b/y ion counting)
//!        └── Build SearchResult
//! ```

pub mod digest;
pub mod error;
pub mod fasta;
pub mod progress;
pub mod registry;

pub use error::SearchEngineError;
pub use progress::SearchProgress;
pub use registry::EngineRegistry;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crate_compiles() {
        let registry = EngineRegistry::new();
        assert!(registry.is_empty());
    }
}
