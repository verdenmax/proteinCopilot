//! Search engine adapter trait and related types.
//!
//! This module defines the [`SearchEngineAdapter`] trait that all search
//! engine integrations (pFind, MSFragger, Comet, etc.) must implement.
//!
//! The trait enforces a uniform interface for:
//! - Running a search against input spectrum files
//! - Querying engine availability and version
//! - Health-checking the engine (local binary or remote SSH host)
//!
//! Concrete adapter implementations live in the `mcp-search-engine` crate
//! under `src/adapters/`.

use std::path::PathBuf;

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use uuid::Uuid;

use crate::error::CoreError;
use crate::progress::ProgressCallback;
use crate::search_params::SearchParams;
use crate::search_result::SearchResult;

// ---------------------------------------------------------------------------
// EngineInfo
// ---------------------------------------------------------------------------

/// Static metadata about a search engine.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct EngineInfo {
    /// Engine name (e.g. "pFind", "MSFragger", "Comet").
    pub name: String,
    /// Engine version string (e.g. "3.1.0").
    pub version: String,
    /// Features supported by this engine (e.g. "open_search", "glyco").
    pub supported_features: Vec<String>,
}

// ---------------------------------------------------------------------------
// HealthStatus
// ---------------------------------------------------------------------------

/// Health status of a search engine installation.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub enum HealthStatus {
    /// Engine is fully operational.
    Healthy,
    /// Engine is running but with issues (e.g. low disk space).
    Degraded {
        /// Description of the degradation.
        reason: String,
    },
    /// Engine is not available.
    Unavailable {
        /// Description of why the engine is unavailable.
        reason: String,
    },
}

// ---------------------------------------------------------------------------
// SearchEngineAdapter trait
// ---------------------------------------------------------------------------

/// Trait that all search engine integrations must implement.
///
/// This is the core abstraction that allows ProteinCopilot to support
/// multiple search engines through a uniform interface. Each adapter
/// handles engine-specific details (binary invocation, parameter file
/// generation, result parsing) internally.
///
/// # Implementors
///
/// - `PFindAdapter` — pFind 3.x via SSH (MVP)
/// - `MSFraggerAdapter` — planned
/// - `CometAdapter` — planned
#[async_trait::async_trait]
pub trait SearchEngineAdapter: Send + Sync {
    /// Execute a search with the given parameters against input spectrum files.
    ///
    /// The `on_progress` callback is invoked at each search stage to report
    /// progress. Callers that don't need progress can pass [`crate::progress::noop_progress()`].
    ///
    /// Returns a standardized [`SearchResult`] regardless of the underlying
    /// engine's native output format.
    async fn search(
        &self,
        params: &SearchParams,
        input_files: &[PathBuf],
        on_progress: ProgressCallback,
    ) -> Result<SearchResult, CoreError>;

    /// Returns static metadata about this engine (name, version, features).
    fn engine_info(&self) -> EngineInfo;

    /// Checks whether the engine is available and healthy.
    ///
    /// For SSH-based engines, this verifies network connectivity and
    /// that the engine binary is accessible on the remote host.
    async fn health_check(&self) -> Result<HealthStatus, CoreError>;

    /// Cancel a running search identified by `run_id`.
    ///
    /// Default implementation is a no-op — the MCP server layer uses
    /// `JoinHandle::abort()` to terminate the local task.
    /// pFind adapter can override this to SSH-kill the remote process.
    async fn cancel(&self, _run_id: Uuid) -> Result<(), CoreError> {
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn engine_info_serde_roundtrip() {
        let info = EngineInfo {
            name: "pFind".to_string(),
            version: "3.1.0".to_string(),
            supported_features: vec![
                "open_search".to_string(),
                "modification_localization".to_string(),
            ],
        };
        let json = serde_json::to_string_pretty(&info).unwrap();
        let back: EngineInfo = serde_json::from_str(&json).unwrap();
        assert_eq!(info, back);
    }

    #[test]
    fn engine_info_empty_features() {
        let info = EngineInfo {
            name: "Comet".to_string(),
            version: "2024.01.0".to_string(),
            supported_features: vec![],
        };
        let json = serde_json::to_string(&info).unwrap();
        let back: EngineInfo = serde_json::from_str(&json).unwrap();
        assert!(back.supported_features.is_empty());
    }

    #[test]
    fn health_status_healthy_serde_roundtrip() {
        let status = HealthStatus::Healthy;
        let json = serde_json::to_string(&status).unwrap();
        let back: HealthStatus = serde_json::from_str(&json).unwrap();
        assert_eq!(status, back);
    }

    #[test]
    fn health_status_degraded_serde_roundtrip() {
        let status = HealthStatus::Degraded {
            reason: "Disk usage above 90%".to_string(),
        };
        let json = serde_json::to_string_pretty(&status).unwrap();
        let back: HealthStatus = serde_json::from_str(&json).unwrap();
        assert_eq!(status, back);
    }

    #[test]
    fn health_status_unavailable_serde_roundtrip() {
        let status = HealthStatus::Unavailable {
            reason: "SSH connection refused".to_string(),
        };
        let json = serde_json::to_string_pretty(&status).unwrap();
        let back: HealthStatus = serde_json::from_str(&json).unwrap();
        assert_eq!(status, back);
    }
}
