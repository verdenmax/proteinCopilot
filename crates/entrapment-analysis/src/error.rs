//! Error types for the entrapment analysis crate.

use std::path::PathBuf;

use thiserror::Error;

#[derive(Debug, Error)]
#[allow(clippy::enum_variant_names)]
pub enum EntrapmentError {
    #[error("config error: {detail}")]
    ConfigError { detail: String },

    #[error("config file I/O error at {path}: {detail}")]
    ConfigIoError { path: PathBuf, detail: String },

    #[error("FASTA error at {path}: {detail}")]
    FastaError { path: PathBuf, detail: String },

    #[error("loader error for {format}: {detail}")]
    LoaderError { format: String, detail: String },

    #[error("I/O error at {path}: {detail}")]
    IoError { path: PathBuf, detail: String },

    #[error("output error: {detail}")]
    OutputError { detail: String },

    #[error("report error: {detail}")]
    ReportError { detail: String },

    #[error("tagging error: {detail}")]
    TaggingError { detail: String },

    #[error("provenance error: {detail}")]
    ProvenanceError { detail: String },

    #[error("spectrum read error for {path}: {detail}")]
    SpectrumError { path: PathBuf, detail: String },
}
