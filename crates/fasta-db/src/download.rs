//! HTTPS download with streaming, SHA256 hashing, and protein counting.

use std::path::Path;

use sha2::{Digest, Sha256};
use tokio::io::AsyncWriteExt;

use crate::FastaDbError;

/// Result of a successful download.
#[derive(Debug, Clone)]
pub struct DownloadResult {
    pub file_size_bytes: u64,
    pub protein_count: u64,
    pub sha256: String,
}

/// Downloads a FASTA file from `url` to `dest_path` using streaming.
///
/// Writes to a temporary `.part` file first, then renames atomically.
/// Returns download metadata (size, protein count, SHA256).
pub async fn download_fasta(url: &str, dest_path: &Path) -> Result<DownloadResult, FastaDbError> {
    tracing::info!(url = url, dest = %dest_path.display(), "starting FASTA download");

    let response = reqwest::get(url)
        .await
        .map_err(|e| FastaDbError::NetworkError {
            url: url.to_string(),
            detail: e.to_string(),
        })?;

    if !response.status().is_success() {
        return Err(FastaDbError::NetworkError {
            url: url.to_string(),
            detail: format!("HTTP {}", response.status()),
        });
    }

    let part_path = dest_path.with_extension("fasta.part");

    // Ensure parent directory exists
    if let Some(parent) = dest_path.parent() {
        tokio::fs::create_dir_all(parent)
            .await
            .map_err(|e| FastaDbError::IoError {
                path: parent.to_path_buf(),
                source: e,
            })?;
    }

    let bytes = response
        .bytes()
        .await
        .map_err(|e| FastaDbError::NetworkError {
            url: url.to_string(),
            detail: format!("failed to read response body: {e}"),
        })?;

    if bytes.is_empty() {
        return Err(FastaDbError::DownloadFailed {
            id: url.to_string(),
            detail: "server returned empty response".to_string(),
        });
    }

    // Write to .part file
    let mut file =
        tokio::fs::File::create(&part_path)
            .await
            .map_err(|e| FastaDbError::IoError {
                path: part_path.clone(),
                source: e,
            })?;
    file.write_all(&bytes)
        .await
        .map_err(|e| FastaDbError::IoError {
            path: part_path.clone(),
            source: e,
        })?;
    file.flush().await.map_err(|e| FastaDbError::IoError {
        path: part_path.clone(),
        source: e,
    })?;

    let file_size_bytes = bytes.len() as u64;
    let protein_count = count_fasta_proteins(&bytes);
    let sha256 = compute_sha256(&bytes);

    // Atomic rename
    tokio::fs::rename(&part_path, dest_path)
        .await
        .map_err(|e| FastaDbError::IoError {
            path: dest_path.to_path_buf(),
            source: e,
        })?;

    tracing::info!(
        proteins = protein_count,
        size_bytes = file_size_bytes,
        "FASTA download complete"
    );

    Ok(DownloadResult {
        file_size_bytes,
        protein_count,
        sha256,
    })
}

/// Counts the number of protein entries (lines starting with `>`) in FASTA content.
pub(crate) fn count_fasta_proteins(content: &[u8]) -> u64 {
    content
        .iter()
        .enumerate()
        .filter(|(i, &b)| b == b'>' && (*i == 0 || content[*i - 1] == b'\n'))
        .count() as u64
}

/// Computes the SHA256 hex digest of a byte slice.
pub(crate) fn compute_sha256(data: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(data);
    format!("{:x}", hasher.finalize())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn count_proteins_in_fasta_content() {
        let content = b">sp|P1|TEST\nACDE\n>sp|P2|TEST2\nFGHI\n";
        let count = count_fasta_proteins(content);
        assert_eq!(count, 2);
    }

    #[test]
    fn count_proteins_empty() {
        let count = count_fasta_proteins(b"");
        assert_eq!(count, 0);
    }

    #[test]
    fn compute_sha256_deterministic() {
        let hash1 = compute_sha256(b"hello world");
        let hash2 = compute_sha256(b"hello world");
        assert_eq!(hash1, hash2);
        assert_eq!(hash1.len(), 64);
    }
}
