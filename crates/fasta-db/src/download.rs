//! HTTPS download with streaming, SHA256 hashing, and protein counting.

use std::path::Path;

use sha2::{Digest, Sha256};
use tokio::io::AsyncWriteExt;

use crate::FastaDbError;

/// Maximum FASTA download size (10 GiB).
///
/// A generous bound that guards against unbounded disk usage from a malicious
/// or misconfigured server while comfortably accommodating the largest real
/// proteome databases.
const MAX_FASTA_DOWNLOAD_BYTES: u64 = 10 * 1024 * 1024 * 1024;

/// Maximum time to wait for the TCP/TLS connection to be established.
const CONNECT_TIMEOUT: std::time::Duration = std::time::Duration::from_secs(30);

/// Maximum idle time between received body bytes. A stalled server (no data for
/// this long) fails the read instead of hanging forever; this is per-read, not a
/// total cap, so large legitimate downloads still complete.
const READ_TIMEOUT: std::time::Duration = std::time::Duration::from_secs(60);

/// Builds the HTTP client used for FASTA downloads with bounded connect and
/// read timeouts so a slow or stalled server can never hang the download.
fn build_download_client() -> Result<reqwest::Client, FastaDbError> {
    reqwest::Client::builder()
        .connect_timeout(CONNECT_TIMEOUT)
        .read_timeout(READ_TIMEOUT)
        .build()
        .map_err(|e| FastaDbError::NetworkError {
            url: String::new(),
            detail: format!("failed to build HTTP client: {e}"),
        })
}

/// Result of a successful download.
#[derive(Debug, Clone)]
pub struct DownloadResult {
    pub file_size_bytes: u64,
    pub protein_count: u64,
    pub sha256: String,
}

/// Incrementally counts FASTA protein headers (a `>` at the start of a line)
/// across streamed chunks that may split lines at arbitrary byte boundaries.
struct StreamingProteinCounter {
    count: u64,
    prev_byte: u8,
}

impl StreamingProteinCounter {
    /// `prev_byte` starts as `b'\n'` so a leading `>` at offset 0 counts.
    fn new() -> Self {
        Self {
            count: 0,
            prev_byte: b'\n',
        }
    }

    fn update(&mut self, chunk: &[u8]) {
        for &b in chunk {
            if b == b'>' && self.prev_byte == b'\n' {
                self.count += 1;
            }
            self.prev_byte = b;
        }
    }

    fn count(&self) -> u64 {
        self.count
    }
}

/// Downloads a FASTA file from `url` to `dest_path` using streaming.
///
/// The response body is consumed chunk-by-chunk (never fully buffered in
/// memory): each chunk is written to a temporary `.part` file while the SHA256
/// digest and protein count are updated incrementally. The total size is bounded
/// by [`MAX_FASTA_DOWNLOAD_BYTES`]. On success the `.part` file is renamed
/// atomically onto `dest_path`. Returns download metadata (size, protein count,
/// SHA256).
pub async fn download_fasta(url: &str, dest_path: &Path) -> Result<DownloadResult, FastaDbError> {
    tracing::info!(url = url, dest = %dest_path.display(), "starting FASTA download");

    let client = build_download_client()?;
    let mut response = client
        .get(url)
        .send()
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

    // Write to .part file
    let mut file =
        tokio::fs::File::create(&part_path)
            .await
            .map_err(|e| FastaDbError::IoError {
                path: part_path.clone(),
                source: e,
            })?;

    let mut hasher = Sha256::new();
    let mut counter = StreamingProteinCounter::new();
    let mut total_bytes: u64 = 0;

    loop {
        let chunk = match response.chunk().await {
            Ok(Some(chunk)) => chunk,
            Ok(None) => break,
            Err(e) => {
                let _ = tokio::fs::remove_file(&part_path).await;
                return Err(FastaDbError::NetworkError {
                    url: url.to_string(),
                    detail: format!("failed to read response body: {e}"),
                });
            }
        };

        total_bytes += chunk.len() as u64;
        if total_bytes > MAX_FASTA_DOWNLOAD_BYTES {
            let _ = tokio::fs::remove_file(&part_path).await;
            return Err(FastaDbError::DownloadFailed {
                id: url.to_string(),
                detail: format!(
                    "download exceeds maximum allowed size of {MAX_FASTA_DOWNLOAD_BYTES} bytes"
                ),
            });
        }

        if let Err(e) = file.write_all(&chunk).await {
            let _ = tokio::fs::remove_file(&part_path).await;
            return Err(FastaDbError::IoError {
                path: part_path.clone(),
                source: e,
            });
        }
        hasher.update(&chunk);
        counter.update(&chunk);
    }

    if let Err(e) = file.flush().await {
        let _ = tokio::fs::remove_file(&part_path).await;
        return Err(FastaDbError::IoError {
            path: part_path.clone(),
            source: e,
        });
    }

    if total_bytes == 0 {
        let _ = tokio::fs::remove_file(&part_path).await;
        return Err(FastaDbError::DownloadFailed {
            id: url.to_string(),
            detail: "server returned empty response".to_string(),
        });
    }

    let file_size_bytes = total_bytes;
    let protein_count = counter.count();
    let sha256 = format!("{:x}", hasher.finalize());

    // Atomic rename
    if let Err(e) = tokio::fs::rename(&part_path, dest_path).await {
        let _ = tokio::fs::remove_file(&part_path).await;
        return Err(FastaDbError::IoError {
            path: dest_path.to_path_buf(),
            source: e,
        });
    }

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
///
/// One-shot reference implementation retained for unit tests; the production
/// download path uses [`StreamingProteinCounter`] for chunk-boundary safety.
#[cfg(test)]
pub(crate) fn count_fasta_proteins(content: &[u8]) -> u64 {
    content
        .iter()
        .enumerate()
        .filter(|(i, &b)| b == b'>' && (*i == 0 || content[*i - 1] == b'\n'))
        .count() as u64
}

/// Computes the SHA256 hex digest of a byte slice.
///
/// One-shot reference implementation retained for unit tests; the production
/// download path hashes incrementally over streamed chunks.
#[cfg(test)]
pub(crate) fn compute_sha256(data: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(data);
    format!("{:x}", hasher.finalize())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test-only helper: feed a sequence of chunks through the streaming
    /// protein counter (mirrors the production download loop).
    fn count_proteins_streamed(chunks: &[&[u8]]) -> u64 {
        let mut counter = StreamingProteinCounter::new();
        for chunk in chunks {
            counter.update(chunk);
        }
        counter.count()
    }

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

    #[test]
    fn streamed_protein_count_matches_oneshot_across_boundaries() {
        let full: &[u8] = b">P1\nACDE\n>P2\nFGHI\n";
        let oneshot = count_fasta_proteins(full);
        assert_eq!(oneshot, 2);

        // Split right before the second '>' header.
        let split_before_header: Vec<&[u8]> = vec![b">P1\nACDE\n", b">P2\n", b"FGHI\n"];
        assert_eq!(count_proteins_streamed(&split_before_header), oneshot);

        // Split with a '>' landing at the very END of a chunk (its preceding
        // '\n' is in the same chunk, but the header continues in the next).
        let split_on_header: Vec<&[u8]> = vec![b">P1\nAC", b"DE\n>", b"P2\nFGHI\n"];
        assert_eq!(count_proteins_streamed(&split_on_header), oneshot);

        // Mid-line split and byte-by-byte both agree.
        let split_midline: Vec<&[u8]> = vec![b">P", b"1\nACDE", b"\n>P2\nFG", b"HI\n"];
        assert_eq!(count_proteins_streamed(&split_midline), oneshot);

        let byte_by_byte: Vec<&[u8]> = full.chunks(1).collect();
        assert_eq!(count_proteins_streamed(&byte_by_byte), oneshot);
    }

    #[test]
    fn streamed_sha256_matches_oneshot() {
        let full: &[u8] = b">P1\nACDE\n>P2\nFGHI\nMOREDATA\n";
        let expected = compute_sha256(full);

        // Feed the same bytes in awkward chunks to an incremental hasher.
        let chunks: [&[u8]; 3] = [b">P1\nAC", b"DE\n>P2\nFG", b"HI\nMOREDATA\n"];
        let mut hasher = Sha256::new();
        for chunk in &chunks {
            hasher.update(chunk);
        }
        let streamed = format!("{:x}", hasher.finalize());

        assert_eq!(streamed, expected);
    }

    #[test]
    fn download_client_builds_with_timeouts() {
        // The download client must be constructible and carry bounded
        // connect/read timeouts so a stalled server cannot hang indefinitely.
        let client = build_download_client();
        assert!(client.is_ok(), "client must build: {:?}", client.err());
        assert!(CONNECT_TIMEOUT > std::time::Duration::ZERO);
        assert!(READ_TIMEOUT > std::time::Duration::ZERO);
    }
}
