//! # ProteinCopilot Search Engine
//!
//! Library crate for search engine orchestration and execution.
//! Contains a simplified built-in search engine for MVP validation,
//! a Sage adapter for production-grade proteomics search via sage-core,
//! and architecture ready for pFind/MSFragger/Comet adapters.
//!
//! # Architecture
//!
//! ```text
//! SearchParams + input_files
//!        │
//!        ▼
//!  EngineRegistry.get("SimpleSearch" | "Sage")
//!        │
//!        ▼
//!  SearchEngineAdapter::search()
//!        │
//!        ├── SimpleSearch: Read FASTA → digest → b/y ion scoring
//!        └── Sage: sage-core library → rayon parallel scoring → LDA rescoring
//! ```

pub mod adapters;
pub mod annotate;
pub mod chemistry;
pub mod digest;
pub mod error;
pub mod fasta;
pub mod matching;
pub mod progress;
pub mod registry;
pub mod simple_engine;
pub mod varmod;

pub use error::SearchEngineError;
pub use progress::SearchProgress;
pub use registry::EngineRegistry;
pub use simple_engine::SimpleSearchEngine;
pub use adapters::sage::SageAdapter;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crate_compiles() {
        let registry = EngineRegistry::new();
        assert!(registry.is_empty());
    }
}
