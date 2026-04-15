//! Error types for the fasta-db crate.

use std::path::PathBuf;

#[derive(Debug, thiserror::Error)]
pub enum FastaDbError {
    #[error("network error downloading {url}: {detail}")]
    NetworkError { url: String, detail: String },

    #[error("I/O error at {path}: {source}")]
    IoError {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("cache registry error: {detail}")]
    RegistryError { detail: String },

    #[error("unknown database '{id}'; available: {}", available.join(", "))]
    UnknownDatabase { id: String, available: Vec<String> },

    #[error("download failed for '{id}': {detail}")]
    DownloadFailed { id: String, detail: String },
}
