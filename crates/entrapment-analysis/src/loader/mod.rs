//! Result loaders for reading search-engine output into [`UnifiedPsm`] records.
//!
//! Supports:
//! - **DIA-NN parquet** — via [`diann_parquet::load_diann_parquet`]
//! - **Generic TSV/TXT** — via [`generic_tsv::load_generic_tsv`] with configurable column mapping

pub mod diann_parquet;
pub mod generic_tsv;

use std::path::Path;

use serde::{Deserialize, Serialize};

pub use generic_tsv::TsvColumnMap;

use crate::error::EntrapmentError;
use crate::types::UnifiedPsm;

/// Supported search-result file formats.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ResultFormat {
    /// DIA-NN `report.parquet` format.
    DiannParquet,
    /// Generic tab-/delimiter-separated text format.
    GenericTsv,
}

impl ResultFormat {
    /// Detect the result format from the file extension.
    ///
    /// - `.parquet` → [`ResultFormat::DiannParquet`]
    /// - `.tsv` or `.txt` → [`ResultFormat::GenericTsv`]
    /// - Other extensions produce a [`EntrapmentError::LoaderError`].
    pub fn from_path(path: &Path) -> Result<Self, EntrapmentError> {
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();

        match ext.as_str() {
            "parquet" => Ok(Self::DiannParquet),
            "tsv" | "txt" => Ok(Self::GenericTsv),
            other => Err(EntrapmentError::LoaderError {
                format: other.to_string(),
                detail: format!(
                    "unsupported result file extension '.{other}' for '{}'",
                    path.display()
                ),
            }),
        }
    }
}

/// Load PSMs from a search-result file.
///
/// Dispatches to the appropriate format-specific loader based on `format`.
/// For [`ResultFormat::GenericTsv`], an optional [`TsvColumnMap`] configures
/// column header mapping; if `None`, defaults are used.
pub fn load_psms(
    path: &Path,
    format: &ResultFormat,
    tsv_config: Option<&TsvColumnMap>,
) -> Result<Vec<UnifiedPsm>, EntrapmentError> {
    match format {
        ResultFormat::DiannParquet => diann_parquet::load_diann_parquet(path),
        ResultFormat::GenericTsv => {
            let default_map = TsvColumnMap::default();
            let map = tsv_config.unwrap_or(&default_map);
            generic_tsv::load_generic_tsv(path, map)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::path::PathBuf;

    #[test]
    fn test_format_detection() {
        // .parquet → DiannParquet
        let fmt = ResultFormat::from_path(Path::new("report.parquet"));
        assert!(fmt.is_ok());
        assert!(matches!(
            fmt.as_ref().ok(),
            Some(ResultFormat::DiannParquet)
        ));

        // .tsv → GenericTsv
        let fmt = ResultFormat::from_path(Path::new("results.tsv"));
        assert!(fmt.is_ok());
        assert!(matches!(fmt.as_ref().ok(), Some(ResultFormat::GenericTsv)));

        // .txt → GenericTsv
        let fmt = ResultFormat::from_path(Path::new("results.txt"));
        assert!(fmt.is_ok());
        assert!(matches!(fmt.as_ref().ok(), Some(ResultFormat::GenericTsv)));

        // .xyz → error
        let fmt = ResultFormat::from_path(Path::new("results.xyz"));
        assert!(fmt.is_err());

        // No extension → error
        let fmt = ResultFormat::from_path(Path::new("results"));
        assert!(fmt.is_err());

        // Case-insensitive
        let fmt = ResultFormat::from_path(Path::new("report.PARQUET"));
        assert!(matches!(
            fmt.as_ref().ok(),
            Some(ResultFormat::DiannParquet)
        ));
    }

    #[test]
    fn test_load_tsv() {
        let dir = tempfile::tempdir().ok();
        let dir = dir.as_ref();
        // Fallback to a fixed temp path if tempdir fails (shouldn't happen)
        let tsv_path: PathBuf = if let Some(d) = dir {
            d.path().join("test_results.tsv")
        } else {
            PathBuf::from("/dev/shm/test_results.tsv")
        };

        // Write a TSV file with default column names
        {
            let mut f = std::fs::File::create(&tsv_path)
                .map_err(|e| format!("failed to create temp file: {e}"))
                .unwrap(); // OK in test code
            writeln!(f, "peptide\tcharge\tprecursor_mz\tretention_time\tscan_number\tspectrum_file\tprotein_ids\tq_value").unwrap();
            writeln!(
                f,
                "PEPTIDEK\t2\t450.123\t12.5\t1001\trun1.raw\tP12345;P67890\t0.01"
            )
            .unwrap();
            writeln!(
                f,
                "ANOTHERPEPTIDE\t3\t600.789\t25.3\t2002\trun2.raw\tQ11111\t0.005"
            )
            .unwrap();
        }

        let map = TsvColumnMap::default();
        let psms = generic_tsv::load_generic_tsv(&tsv_path, &map).unwrap(); // OK in test code

        assert_eq!(psms.len(), 2);

        // Row 1
        assert_eq!(psms[0].peptide, "PEPTIDEK");
        assert_eq!(psms[0].charge, Some(2));
        assert!((psms[0].precursor_mz.unwrap_or(0.0) - 450.123).abs() < 1e-6); // OK in test code
        assert!((psms[0].retention_time.unwrap_or(0.0) - 12.5).abs() < 1e-6); // OK in test code
        assert_eq!(psms[0].scan_number, Some(1001));
        assert_eq!(psms[0].spectrum_file.as_deref(), Some("run1.raw"));
        assert_eq!(psms[0].protein_ids, "P12345;P67890");
        assert!((psms[0].q_value.unwrap_or(1.0) - 0.01).abs() < 1e-6); // OK in test code

        // Row 2
        assert_eq!(psms[1].peptide, "ANOTHERPEPTIDE");
        assert_eq!(psms[1].charge, Some(3));
        assert!((psms[1].precursor_mz.unwrap_or(0.0) - 600.789).abs() < 1e-6); // OK in test code
        assert!((psms[1].retention_time.unwrap_or(0.0) - 25.3).abs() < 1e-6); // OK in test code
        assert_eq!(psms[1].scan_number, Some(2002));
        assert_eq!(psms[1].spectrum_file.as_deref(), Some("run2.raw"));
        assert_eq!(psms[1].protein_ids, "Q11111");
        assert!((psms[1].q_value.unwrap_or(1.0) - 0.005).abs() < 1e-6); // OK in test code
    }

    #[test]
    fn test_load_tsv_missing_optional_columns() {
        let dir = tempfile::tempdir().ok();
        let dir = dir.as_ref();
        let tsv_path: PathBuf = if let Some(d) = dir {
            d.path().join("test_minimal.tsv")
        } else {
            PathBuf::from("/dev/shm/test_minimal.tsv")
        };

        {
            let mut f = std::fs::File::create(&tsv_path).unwrap(); // OK in test code
            writeln!(f, "peptide\tprotein_ids").unwrap();
            writeln!(f, "PEPTIDEK\tP12345").unwrap();
        }

        let map = TsvColumnMap::default();
        let psms = generic_tsv::load_generic_tsv(&tsv_path, &map).unwrap(); // OK in test code

        assert_eq!(psms.len(), 1);
        assert_eq!(psms[0].peptide, "PEPTIDEK");
        assert_eq!(psms[0].protein_ids, "P12345");
        assert_eq!(psms[0].charge, None);
        assert_eq!(psms[0].precursor_mz, None);
        assert_eq!(psms[0].retention_time, None);
        assert_eq!(psms[0].scan_number, None);
        assert_eq!(psms[0].spectrum_file, None);
        assert_eq!(psms[0].q_value, None);
    }

    #[test]
    fn test_load_psms_dispatches_tsv() {
        let dir = tempfile::tempdir().ok();
        let dir = dir.as_ref();
        let tsv_path: PathBuf = if let Some(d) = dir {
            d.path().join("dispatch_test.tsv")
        } else {
            PathBuf::from("/dev/shm/dispatch_test.tsv")
        };

        {
            let mut f = std::fs::File::create(&tsv_path).unwrap(); // OK in test code
            writeln!(f, "peptide\tprotein_ids\tq_value").unwrap();
            writeln!(f, "TESTPEPTIDE\tP99999\t0.02").unwrap();
        }

        let format = ResultFormat::from_path(&tsv_path).unwrap(); // OK in test code
        let psms = load_psms(&tsv_path, &format, None).unwrap(); // OK in test code

        assert_eq!(psms.len(), 1);
        assert_eq!(psms[0].peptide, "TESTPEPTIDE");
    }
}
