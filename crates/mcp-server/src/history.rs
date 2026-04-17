//! Search history persistence — JSON files in `~/.protein-copilot/history/`.
//!
//! Each completed/failed/cancelled search is persisted as a JSON file
//! containing metadata and summary statistics (not the full PSM list).
//! The MCP server loads history on startup and provides it via `list_searches`.

use std::fs;
use std::path::{Path, PathBuf};

use chrono::{DateTime, Utc};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use uuid::Uuid;

use protein_copilot_core::engine::EngineInfo;
use protein_copilot_core::search_params::SearchParams;

/// Maximum number of history entries before FIFO eviction.
const MAX_HISTORY: usize = 500;

/// Summary metadata persisted for each search run.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SearchHistoryEntry {
    /// Unique run identifier.
    pub run_id: Uuid,
    /// Final status: "Completed", "Failed: ...", "Cancelled".
    pub status: String,
    /// When the search was started.
    pub created_at: DateTime<Utc>,
    /// Total elapsed time in seconds.
    pub elapsed_sec: f64,
    /// Search engine used.
    pub engine_info: EngineInfo,
    /// Input spectrum files.
    pub input_files: Vec<PathBuf>,
    /// Search parameters used.
    pub params_used: SearchParams,
    /// Total PSM count (before FDR filtering).
    pub total_psms: Option<u64>,
    /// PSMs at 1% FDR.
    pub psms_at_1pct_fdr: Option<u64>,
    /// Identification rate (psms_at_fdr / total_spectra).
    pub identification_rate: Option<f64>,
    /// Protein groups identified.
    pub protein_groups: Option<u64>,
}

/// Returns the history directory path, creating it if needed.
pub fn history_dir() -> Option<PathBuf> {
    let home = dirs::home_dir()?;
    let dir = home.join(".protein-copilot").join("history");
    fs::create_dir_all(&dir).ok()?;
    Some(dir)
}

/// Save a history entry to disk.
pub fn save_entry(entry: &SearchHistoryEntry) {
    let Some(dir) = history_dir() else {
        tracing::warn!("cannot determine history directory; skipping persistence");
        return;
    };
    let path = dir.join(format!("{}.json", entry.run_id));
    match serde_json::to_string_pretty(entry) {
        Ok(json) => {
            if let Err(e) = fs::write(&path, json) {
                tracing::warn!("failed to write history {}: {e}", path.display());
            }
        }
        Err(e) => tracing::warn!("failed to serialize history: {e}"),
    }
    evict_oldest(&dir);
}

/// Load all history entries from disk, sorted by created_at descending.
pub fn load_all() -> Vec<SearchHistoryEntry> {
    let Some(dir) = history_dir() else {
        return Vec::new();
    };
    let Ok(read_dir) = fs::read_dir(&dir) else {
        return Vec::new();
    };

    let mut entries = Vec::new();
    for entry in read_dir.flatten() {
        let path = entry.path();
        if path.extension().and_then(|e| e.to_str()) != Some("json") {
            continue;
        }
        match fs::read_to_string(&path) {
            Ok(content) => match serde_json::from_str::<SearchHistoryEntry>(&content) {
                Ok(e) => entries.push(e),
                Err(err) => tracing::warn!("corrupt history file {}: {err}", path.display()),
            },
            Err(e) => tracing::warn!("cannot read {}: {e}", path.display()),
        }
    }
    entries.sort_by(|a, b| b.created_at.cmp(&a.created_at));
    entries
}

/// FIFO eviction: delete oldest entries when over MAX_HISTORY.
fn evict_oldest(dir: &Path) {
    let mut entries = load_all();
    while entries.len() > MAX_HISTORY {
        if let Some(oldest) = entries.pop() {
            let path = dir.join(format!("{}.json", oldest.run_id));
            let _ = fs::remove_file(path);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::engine::EngineInfo;
    use protein_copilot_core::search_params::*;

    fn test_entry(run_id: Uuid) -> SearchHistoryEntry {
        SearchHistoryEntry {
            run_id,
            status: "Completed".to_string(),
            created_at: Utc::now(),
            elapsed_sec: 5.0,
            engine_info: EngineInfo {
                name: "SimpleSearch".to_string(),
                version: "0.1.0".to_string(),
                supported_features: vec![],
            },
            input_files: vec![PathBuf::from("test.mgf")],
            params_used: SearchParams {
                enzyme: Enzyme::Trypsin,
                missed_cleavages: 2,
                fixed_modifications: vec![],
                variable_modifications: vec![],
                precursor_tolerance: MassTolerance {
                    value: 10.0,
                    unit: ToleranceUnit::Ppm,
                },
                fragment_tolerance: MassTolerance {
                    value: 0.02,
                    unit: ToleranceUnit::Da,
                },
                database_path: "/data/test.fasta".to_string(),
                decoy_strategy: DecoyStrategy::Reverse,
                acquisition_mode: None,
                max_variable_modifications: 3,
                min_peptide_length: 7,
                max_peptide_length: 50,
                engine: None,
            },
            total_psms: Some(100),
            psms_at_1pct_fdr: Some(80),
            identification_rate: Some(0.32),
            protein_groups: Some(40),
        }
    }

    #[test]
    fn history_entry_serde_roundtrip() {
        let entry = test_entry(Uuid::new_v4());
        let json = serde_json::to_string_pretty(&entry).unwrap();
        let back: SearchHistoryEntry = serde_json::from_str(&json).unwrap();
        assert_eq!(entry.run_id, back.run_id);
        assert_eq!(entry.status, back.status);
        assert_eq!(entry.total_psms, back.total_psms);
    }

    #[test]
    fn history_dir_is_created() {
        let dir = history_dir();
        assert!(dir.is_some());
        assert!(dir.unwrap().exists());
    }
}
