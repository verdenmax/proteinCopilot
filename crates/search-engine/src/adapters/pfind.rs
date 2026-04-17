//! pFind search engine adapter (stub — not yet implemented).
//!
//! This module will implement the full pFind integration:
//! - `SearchParams → .cfg` file generation
//! - SSH remote execution via `tokio::process::Command`
//! - pFind output parsing → `SearchResult`
//! - Progress tracking via log file polling
//!
//! # Prerequisites for implementation
//!
//! - pFind `.cfg` file format specification
//! - pFind output file format samples (`.spectra`, `.protein`, etc.)
//! - SSH access to a server with pFind installed
//! - Test fixture files from a real pFind run

use std::path::PathBuf;

use protein_copilot_core::engine::{EngineInfo, HealthStatus, SearchEngineAdapter};
use protein_copilot_core::error::CoreError;
use protein_copilot_core::progress::ProgressCallback;
use protein_copilot_core::search_params::SearchParams;
use protein_copilot_core::search_result::SearchResult;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use uuid::Uuid;

/// SSH connection configuration for remote pFind execution.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SshConfig {
    /// Remote host address (e.g., "compute-01.lab.edu").
    pub host: String,
    /// SSH port (default 22).
    pub port: u16,
    /// SSH username.
    pub user: String,
    /// Path to SSH private key file.
    pub key_path: PathBuf,
    /// Path to pFind executable on the remote server.
    pub pfind_executable: PathBuf,
    /// Remote working directory for search files.
    pub work_dir: PathBuf,
}

impl Default for SshConfig {
    fn default() -> Self {
        Self {
            host: "localhost".to_string(),
            port: 22,
            user: String::new(),
            key_path: PathBuf::new(),
            pfind_executable: PathBuf::from("/usr/local/bin/pfind"),
            work_dir: PathBuf::from("/tmp/pfind_work"),
        }
    }
}

/// pFind search engine adapter.
///
/// Executes pFind on a remote server via SSH. Currently a stub —
/// all methods return `todo!()` errors pending real implementation.
///
/// # Implementation Roadmap
///
/// 1. `PFindConfig`: map `SearchParams` → pFind `.cfg` format
/// 2. SSH file transfer: upload `.cfg`, download results
/// 3. Remote execution: `ssh user@host pfind <cfg_path>`
/// 4. Result parsing: parse `.spectra` output → `SearchResult`
/// 5. Progress polling: read pFind log for completion percentage
pub struct PFindAdapter {
    /// SSH connection configuration.
    pub ssh_config: SshConfig,
}

impl PFindAdapter {
    /// Creates a new pFind adapter with the given SSH configuration.
    pub fn new(ssh_config: SshConfig) -> Self {
        Self { ssh_config }
    }
}

#[async_trait::async_trait]
impl SearchEngineAdapter for PFindAdapter {
    async fn search(
        &self,
        _params: &SearchParams,
        _input_files: &[PathBuf],
        _on_progress: ProgressCallback,
    ) -> Result<SearchResult, CoreError> {
        Err(CoreError::SearchEngineError {
            engine: "pFind".to_string(),
            detail: "pFind adapter not yet implemented".to_string(),
            suggestion: "Use SimpleSearch engine for now, or provide pFind configuration details"
                .to_string(),
        })
    }

    fn engine_info(&self) -> EngineInfo {
        EngineInfo {
            name: "pFind".to_string(),
            version: "3.x (not connected)".to_string(),
            supported_features: vec![
                "open_search".to_string(),
                "modification_localization".to_string(),
            ],
        }
    }

    async fn health_check(&self) -> Result<HealthStatus, CoreError> {
        Ok(HealthStatus::Unavailable {
            reason: "pFind adapter not yet implemented".to_string(),
        })
    }

    async fn cancel(&self, _run_id: Uuid) -> Result<(), CoreError> {
        Err(CoreError::SearchEngineError {
            engine: "pFind".to_string(),
            detail: "cancel not yet implemented".to_string(),
            suggestion: "pFind remote cancellation requires SSH integration".to_string(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::progress::noop_progress;

    #[test]
    fn ssh_config_default() {
        let config = SshConfig::default();
        assert_eq!(config.port, 22);
        assert_eq!(config.host, "localhost");
    }

    #[test]
    fn pfind_engine_info() {
        let adapter = PFindAdapter::new(SshConfig::default());
        let info = adapter.engine_info();
        assert_eq!(info.name, "pFind");
    }

    #[tokio::test]
    async fn pfind_health_check_unavailable() {
        let adapter = PFindAdapter::new(SshConfig::default());
        let status = adapter.health_check().await.unwrap();
        assert!(matches!(status, HealthStatus::Unavailable { .. }));
    }

    #[tokio::test]
    async fn pfind_search_returns_not_implemented() {
        let adapter = PFindAdapter::new(SshConfig::default());
        let params = SearchParams {
            database_path: "/db.fasta".to_string(),
            enzyme: protein_copilot_core::search_params::Enzyme::Trypsin,
            missed_cleavages: 2,
            fixed_modifications: vec![],
            variable_modifications: vec![],
            precursor_tolerance: protein_copilot_core::search_params::MassTolerance {
                value: 10.0,
                unit: protein_copilot_core::search_params::ToleranceUnit::Ppm,
            },
            fragment_tolerance: protein_copilot_core::search_params::MassTolerance {
                value: 0.02,
                unit: protein_copilot_core::search_params::ToleranceUnit::Da,
            },
            decoy_strategy: protein_copilot_core::search_params::DecoyStrategy::Reverse,
            acquisition_mode: None,
            max_variable_modifications: 3,
            min_peptide_length: 7,
            max_peptide_length: 50,
            engine: None,
        };
        let result = adapter.search(&params, &[], noop_progress()).await;
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("not yet implemented"));
    }
}
