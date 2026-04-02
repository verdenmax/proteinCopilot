//! Run metadata for tracking analysis provenance.
//!
//! Every search run is tracked by a [`RunMetadata`] struct that captures:
//! - **What** was searched (input files, parameters, engine)
//! - **When** it was started and how long it took
//! - **Status** of the run (pending, running, completed, failed)
//!
//! This ensures full reproducibility and auditability per §2.5 of
//! copilot-instructions.md.

use std::path::PathBuf;

use chrono::{DateTime, Utc};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use thiserror::Error;
use uuid::Uuid;

use crate::engine::EngineInfo;
use crate::search_params::SearchParams;

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Errors related to run metadata validation.
#[derive(Debug, Error)]
pub enum RunMetadataError {
    /// `duration_sec` contains NaN or Infinity.
    #[error("duration_sec contains non-finite value")]
    NonFiniteDuration,

    /// `duration_sec` is negative.
    #[error("duration_sec must be non-negative, got {value}")]
    NegativeDuration {
        /// The actual value.
        value: f64,
    },

    /// A required field is empty.
    #[error("{field} must not be empty")]
    EmptyField {
        /// Name of the field.
        field: &'static str,
    },

    /// Delegated search params validation failed.
    #[error("params_used validation failed: {0}")]
    InvalidParams(String),
}

// ---------------------------------------------------------------------------
// RunStatus
// ---------------------------------------------------------------------------

/// Status of an analysis run.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub enum RunStatus {
    /// Run is queued but not yet started.
    Pending,
    /// Run is currently executing.
    Running,
    /// Run completed successfully.
    Completed,
    /// Run failed with an error.
    Failed {
        /// Description of why the run failed.
        reason: String,
    },
}

impl std::fmt::Display for RunStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RunStatus::Pending => write!(f, "Pending"),
            RunStatus::Running => write!(f, "Running"),
            RunStatus::Completed => write!(f, "Completed"),
            RunStatus::Failed { reason } => write!(f, "Failed: {reason}"),
        }
    }
}

// ---------------------------------------------------------------------------
// RunMetadata
// ---------------------------------------------------------------------------

/// Metadata for a single analysis run.
///
/// Generated at the start of each search and updated as the run progresses.
/// Stored alongside results to enable full provenance tracking.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct RunMetadata {
    /// Unique identifier for this run (auto-generated UUIDv4).
    pub run_id: Uuid,
    /// Timestamp when the run was created.
    pub created_at: DateTime<Utc>,
    /// Input spectrum file paths.
    pub input_files: Vec<PathBuf>,
    /// Search parameters used for this run.
    pub params_used: SearchParams,
    /// Search engine that executed this run.
    pub engine_info: EngineInfo,
    /// Total duration in seconds (`None` if not yet completed).
    pub duration_sec: Option<f64>,
    /// Current status of the run.
    pub status: RunStatus,
}

impl RunMetadata {
    /// Creates a new `RunMetadata` with auto-generated `run_id` and timestamp.
    ///
    /// The initial status is [`RunStatus::Pending`] and `duration_sec` is `None`.
    pub fn new(params: SearchParams, engine_info: EngineInfo, input_files: Vec<PathBuf>) -> Self {
        Self {
            run_id: Uuid::new_v4(),
            created_at: Utc::now(),
            input_files,
            params_used: params,
            engine_info,
            duration_sec: None,
            status: RunStatus::Pending,
        }
    }

    /// Validates run metadata fields.
    ///
    /// Checks:
    /// - `duration_sec`, if present, is finite and non-negative
    /// - `input_files` is not empty
    /// - `engine_info.name` and `engine_info.version` are not empty
    pub fn validate(&self) -> Result<(), RunMetadataError> {
        if let Some(d) = self.duration_sec {
            if !d.is_finite() {
                return Err(RunMetadataError::NonFiniteDuration);
            }
            if d < 0.0 {
                return Err(RunMetadataError::NegativeDuration { value: d });
            }
        }
        if self.input_files.is_empty() {
            return Err(RunMetadataError::EmptyField {
                field: "input_files",
            });
        }
        if self.engine_info.name.trim().is_empty() {
            return Err(RunMetadataError::EmptyField {
                field: "engine_info.name",
            });
        }
        if self.engine_info.version.trim().is_empty() {
            return Err(RunMetadataError::EmptyField {
                field: "engine_info.version",
            });
        }
        self.params_used
            .validate()
            .map_err(|e| RunMetadataError::InvalidParams(e.to_string()))
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::EngineInfo;
    use crate::search_params::{DecoyStrategy, Enzyme, MassTolerance, SearchParams, ToleranceUnit};

    fn sample_params() -> SearchParams {
        SearchParams {
            database_path: "/data/human.fasta".to_string(),
            enzyme: Enzyme::Trypsin,
            missed_cleavages: 2,
            fixed_modifications: vec![],
            variable_modifications: vec![],
            precursor_tolerance: MassTolerance {
                value: 20.0,
                unit: ToleranceUnit::Ppm,
            },
            fragment_tolerance: MassTolerance {
                value: 0.02,
                unit: ToleranceUnit::Da,
            },
            decoy_strategy: DecoyStrategy::Reverse,
            acquisition_mode: None,
        }
    }

    fn sample_engine_info() -> EngineInfo {
        EngineInfo {
            name: "pFind".to_string(),
            version: "3.1.0".to_string(),
            supported_features: vec!["open_search".to_string()],
        }
    }

    fn sample_run_metadata() -> RunMetadata {
        RunMetadata::new(
            sample_params(),
            sample_engine_info(),
            vec![
                PathBuf::from("/data/sample1.mzML"),
                PathBuf::from("/data/sample2.mzML"),
            ],
        )
    }

    #[test]
    fn new_generates_unique_run_ids() {
        let r1 = sample_run_metadata();
        let r2 = sample_run_metadata();
        assert_ne!(r1.run_id, r2.run_id);
    }

    #[test]
    fn new_sets_pending_status() {
        let r = sample_run_metadata();
        assert_eq!(r.status, RunStatus::Pending);
        assert!(r.duration_sec.is_none());
    }

    #[test]
    fn new_captures_timestamp() {
        let before = Utc::now();
        let r = sample_run_metadata();
        let after = Utc::now();
        assert!(r.created_at >= before && r.created_at <= after);
    }

    #[test]
    fn run_metadata_serde_roundtrip() {
        let r = sample_run_metadata();
        let json = serde_json::to_string_pretty(&r).unwrap();
        let back: RunMetadata = serde_json::from_str(&json).unwrap();
        assert_eq!(r.run_id, back.run_id);
        assert_eq!(r.status, back.status);
        assert_eq!(r.input_files.len(), back.input_files.len());
        assert_eq!(r.params_used.database_path, back.params_used.database_path);
        assert_eq!(r.engine_info.name, back.engine_info.name);
    }

    #[test]
    fn run_status_display() {
        assert_eq!(RunStatus::Pending.to_string(), "Pending");
        assert_eq!(RunStatus::Running.to_string(), "Running");
        assert_eq!(RunStatus::Completed.to_string(), "Completed");
        let failed = RunStatus::Failed {
            reason: "out of memory".to_string(),
        };
        assert!(failed.to_string().contains("out of memory"));
    }

    #[test]
    fn run_status_serde_roundtrip() {
        for status in [
            RunStatus::Pending,
            RunStatus::Running,
            RunStatus::Completed,
            RunStatus::Failed {
                reason: "timeout".to_string(),
            },
        ] {
            let json = serde_json::to_string(&status).unwrap();
            let back: RunStatus = serde_json::from_str(&json).unwrap();
            assert_eq!(status, back);
        }
    }

    #[test]
    fn run_metadata_with_completed_status() {
        let mut r = sample_run_metadata();
        r.status = RunStatus::Completed;
        r.duration_sec = Some(245.5);
        let json = serde_json::to_string(&r).unwrap();
        let back: RunMetadata = serde_json::from_str(&json).unwrap();
        assert_eq!(back.status, RunStatus::Completed);
        assert_eq!(back.duration_sec, Some(245.5));
    }

    #[test]
    fn run_metadata_with_failed_status() {
        let mut r = sample_run_metadata();
        r.status = RunStatus::Failed {
            reason: "SSH disconnected".to_string(),
        };
        r.duration_sec = Some(10.0);
        let json = serde_json::to_string(&r).unwrap();
        let back: RunMetadata = serde_json::from_str(&json).unwrap();
        if let RunStatus::Failed { reason } = &back.status {
            assert!(reason.contains("SSH"));
        } else {
            panic!("expected Failed status");
        }
    }

    // -- Validation -----------------------------------------------------

    #[test]
    fn validate_passes_for_valid_data() {
        assert!(sample_run_metadata().validate().is_ok());
    }

    #[test]
    fn validate_passes_with_valid_duration() {
        let mut r = sample_run_metadata();
        r.duration_sec = Some(245.5);
        assert!(r.validate().is_ok());
    }

    #[test]
    fn validate_rejects_nan_duration() {
        let mut r = sample_run_metadata();
        r.duration_sec = Some(f64::NAN);
        let err = r.validate().unwrap_err();
        assert!(err.to_string().contains("duration_sec"));
    }

    #[test]
    fn validate_rejects_infinity_duration() {
        let mut r = sample_run_metadata();
        r.duration_sec = Some(f64::INFINITY);
        let err = r.validate().unwrap_err();
        assert!(err.to_string().contains("duration_sec"));
    }

    #[test]
    fn validate_rejects_negative_duration() {
        let mut r = sample_run_metadata();
        r.duration_sec = Some(-1.0);
        let err = r.validate().unwrap_err();
        assert!(err.to_string().contains("duration_sec"));
        assert!(err.to_string().contains("-1"));
    }

    #[test]
    fn validate_rejects_empty_input_files() {
        let mut r = sample_run_metadata();
        r.input_files = vec![];
        let err = r.validate().unwrap_err();
        assert!(err.to_string().contains("input_files"));
    }

    #[test]
    fn validate_rejects_empty_engine_name() {
        let mut r = sample_run_metadata();
        r.engine_info.name = "".to_string();
        let err = r.validate().unwrap_err();
        assert!(err.to_string().contains("engine_info.name"));
    }

    #[test]
    fn validate_delegates_to_params_used() {
        let mut r = sample_run_metadata();
        r.params_used.database_path = "".to_string();
        let err = r.validate().unwrap_err();
        assert!(err.to_string().contains("params_used"));
    }
}
