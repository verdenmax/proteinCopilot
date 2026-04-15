//! Import external proteomics search results and match to mzML scans.
//!
//! Supports:
//! - Custom JSON (hela.json format)
//! - DIA-NN report.parquet
//! - pFind .spectra (skeleton)

pub mod converter;
pub mod custom_json;
pub mod diann;
pub mod error;
pub mod pfind;
pub mod scan_matcher;
pub mod unimod;

use std::collections::HashMap;
use std::path::Path;

use protein_copilot_core::search_params::Modification;
use serde::{Deserialize, Serialize};

pub use error::ResultImportError;

/// A PSM imported from an external search result file.
///
/// RT is always in minutes (converted from source format at parse time).
/// `matched_scan` is `None` until scan matching is performed.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ImportedPsm {
    pub sequence: String,
    pub charge: i32,
    pub precursor_mz: f64,
    /// Retention time in minutes (converted from source format at parse time).
    pub rt_min: f64,
    pub modifications: Vec<Modification>,
    /// Search engine score (e.g. DIA-NN Q.Value). `None` if not available.
    pub score: Option<f64>,
    /// q-value for FDR control. `None` if not available.
    pub q_value: Option<f64>,
    pub protein_accessions: Vec<String>,
    /// Raw file name without extension (used for mzML association).
    pub raw_name: String,
    /// Filled by scan matcher.
    pub matched_scan: Option<u32>,
    /// RT delta in minutes between PSM and matched MS2 scan.
    pub rt_delta_min: Option<f64>,
}

/// Result of the import operation.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
pub struct ImportResult {
    pub run_id: String,
    pub match_report: MatchReport,
    pub imported_psm_count: usize,
    pub unique_peptides: usize,
    pub protein_count: usize,
    pub raw_files: Vec<String>,
}

/// Scan matching quality report.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
pub struct MatchReport {
    pub total_psms: usize,
    pub matched: usize,
    pub unmatched: usize,
    pub median_rt_delta_min: f64,
    pub max_rt_delta_min: f64,
    pub per_file: HashMap<String, FileMatchStats>,
}

/// Per-file matching statistics.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
pub struct FileMatchStats {
    pub total: usize,
    pub matched: usize,
    pub ms2_count: usize,
}

/// Trait for format-specific parsers.
pub trait ResultParser: Send + Sync {
    fn parse(
        &self,
        path: &Path,
        unimod: &unimod::UnimodDb,
    ) -> Result<Vec<ImportedPsm>, ResultImportError>;
}

/// Supported import formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImportFormat {
    CustomJson,
    DiannParquet,
    PFindSpectra,
}

/// Detect format from file extension.
pub fn detect_format(path: &Path) -> Result<ImportFormat, ResultImportError> {
    match path.extension().and_then(|e| e.to_str()) {
        Some("json") => Ok(ImportFormat::CustomJson),
        Some("parquet") => Ok(ImportFormat::DiannParquet),
        Some("spectra") => Ok(ImportFormat::PFindSpectra),
        _ => Err(ResultImportError::FormatDetectionFailed {
            path: path.to_path_buf(),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn detect_format_json() {
        assert_eq!(
            detect_format(Path::new("test.json")).unwrap(),
            ImportFormat::CustomJson
        );
    }

    #[test]
    fn detect_format_parquet() {
        assert_eq!(
            detect_format(Path::new("report.parquet")).unwrap(),
            ImportFormat::DiannParquet
        );
    }

    #[test]
    fn detect_format_spectra() {
        assert_eq!(
            detect_format(Path::new("result.spectra")).unwrap(),
            ImportFormat::PFindSpectra
        );
    }

    #[test]
    fn detect_format_unknown() {
        assert!(detect_format(Path::new("data.csv")).is_err());
    }
}
