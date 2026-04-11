//! Error types for the result-import crate.

use std::path::PathBuf;

/// Errors that can occur during external result import.
#[derive(Debug, thiserror::Error)]
pub enum ResultImportError {
    #[error("file not found: {path}")]
    FileNotFound { path: PathBuf },

    #[error("format detection failed for {path}: expected one of .json, .parquet, .spectra")]
    FormatDetectionFailed { path: PathBuf },

    #[error("unknown Unimod record_id: {0} — check unimod.xml or use builtin database")]
    UnknownUnimodId(u32),

    #[error("mzML file not found for raw name '{raw_name}' in directory {dir} — available files: {available}")]
    MzmlNotFound {
        raw_name: String,
        dir: PathBuf,
        available: String,
    },

    #[error("no MS2 scan found near RT={rt_min:.3} min (±{tolerance_min} min) with precursor_mz={precursor_mz:.4}")]
    NoMatchingScan {
        rt_min: f64,
        tolerance_min: f64,
        precursor_mz: f64,
    },

    #[error("missing required column '{column}' in parquet file — expected columns: {expected}")]
    MissingColumn { column: String, expected: String },

    #[error("JSON parse error: {0}")]
    JsonError(#[from] serde_json::Error),

    #[error("XML parse error: {0}")]
    XmlError(#[from] quick_xml::Error),

    #[error("Parquet read error: {0}")]
    ParquetError(#[from] parquet::errors::ParquetError),

    #[error("Arrow error: {0}")]
    ArrowError(#[from] arrow::error::ArrowError),

    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),

    #[error("invalid modification position {position} for sequence of length {seq_len}")]
    InvalidModPosition { position: usize, seq_len: usize },

    #[error("spectrum-io error: {0}")]
    SpectrumIo(String),

    #[error("{0}")]
    Other(String),
}
