//! FASTA database registry, download, and cache management.
//!
//! Provides a built-in registry of common proteomics databases (UniProt Swiss-Prot
//! for major model organisms + cRAP contaminants), HTTPS download with streaming,
//! and local file caching with metadata tracking.

pub mod cache;
pub mod download;
pub mod error;
pub mod registry;

pub use cache::{CacheManager, CacheRegistry, CachedDatabase};
pub use error::FastaDbError;
pub use registry::DatabaseEntry;

use serde::{Deserialize, Serialize};
use std::path::Path;

/// Status of a database in the combined view.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
pub struct DatabaseStatus {
    pub id: String,
    pub species: String,
    pub db_type: String,
    pub description: String,
    pub status: DownloadStatus,
}

/// Whether a database is available for download or already cached.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(tag = "state")]
pub enum DownloadStatus {
    Available,
    Downloaded {
        file_name: String,
        file_size_bytes: u64,
        protein_count: u64,
        downloaded_at: String,
    },
}

/// Information about a cached database.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
pub struct DatabaseInfo {
    pub id: String,
    pub species: String,
    pub db_type: String,
    pub path: String,
    pub protein_count: u64,
    pub file_size_bytes: u64,
    pub downloaded_at: String,
    pub sha256: String,
    pub first_accessions: Vec<String>,
}

/// Result returned after downloading a database.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
pub struct DownloadDatabaseResult {
    pub id: String,
    pub path: String,
    pub protein_count: u64,
    pub file_size_bytes: u64,
}

/// Lists all built-in databases with their cache status.
pub fn list_databases(cache_dir: &Path) -> Result<Vec<DatabaseStatus>, FastaDbError> {
    let cache = CacheManager::new(cache_dir.to_path_buf());
    let registry = cache.load_registry()?;

    let statuses = self::registry::all_databases()
        .iter()
        .map(|entry| {
            let status = if let Some(cached) = registry.databases.get(entry.id) {
                if cache.fasta_path(entry.id).exists() {
                    DownloadStatus::Downloaded {
                        file_name: cached.file_name.clone(),
                        file_size_bytes: cached.file_size_bytes,
                        protein_count: cached.protein_count,
                        downloaded_at: cached.downloaded_at.to_rfc3339(),
                    }
                } else {
                    DownloadStatus::Available
                }
            } else {
                DownloadStatus::Available
            };

            DatabaseStatus {
                id: entry.id.to_string(),
                species: entry.species.to_string(),
                db_type: entry.db_type.to_string(),
                description: entry.description.to_string(),
                status,
            }
        })
        .collect();

    Ok(statuses)
}

/// Downloads a database by ID into the cache directory.
pub async fn download_database(
    database_id: &str,
    cache_dir: &Path,
    force: bool,
) -> Result<DownloadDatabaseResult, FastaDbError> {
    let entry = registry::get_database(database_id).ok_or_else(|| FastaDbError::UnknownDatabase {
        id: database_id.to_string(),
        available: registry::all_database_ids()
            .into_iter()
            .map(String::from)
            .collect(),
    })?;

    let cache = CacheManager::new(cache_dir.to_path_buf());

    // Return cached path if already downloaded and not forcing
    if !force {
        if let Some(cached) = cache.get_cached(database_id)? {
            let path = cache.fasta_path(database_id);
            if path.exists() {
                return Ok(DownloadDatabaseResult {
                    id: database_id.to_string(),
                    path: path.to_string_lossy().to_string(),
                    protein_count: cached.protein_count,
                    file_size_bytes: cached.file_size_bytes,
                });
            }
        }
    }

    let dest_path = cache.fasta_path(database_id);
    let result = download::download_fasta(entry.url, &dest_path).await?;

    let cached_entry = CachedDatabase {
        id: database_id.to_string(),
        file_name: format!("{}.fasta", database_id),
        downloaded_at: chrono::Utc::now(),
        file_size_bytes: result.file_size_bytes,
        protein_count: result.protein_count,
        sha256: result.sha256,
    };
    cache.save_entry(&cached_entry)?;

    Ok(DownloadDatabaseResult {
        id: database_id.to_string(),
        path: dest_path.to_string_lossy().to_string(),
        protein_count: result.protein_count,
        file_size_bytes: result.file_size_bytes,
    })
}

/// Returns detailed info about a cached database.
pub fn get_database_info(
    database_id: &str,
    cache_dir: &Path,
) -> Result<DatabaseInfo, FastaDbError> {
    let entry = registry::get_database(database_id).ok_or_else(|| FastaDbError::UnknownDatabase {
        id: database_id.to_string(),
        available: registry::all_database_ids()
            .into_iter()
            .map(String::from)
            .collect(),
    })?;

    let cache = CacheManager::new(cache_dir.to_path_buf());
    let cached = cache
        .get_cached(database_id)?
        .ok_or_else(|| FastaDbError::DownloadFailed {
            id: database_id.to_string(),
            detail: "database not yet downloaded; use download_database first".to_string(),
        })?;

    let fasta_path = cache.fasta_path(database_id);
    let first_accessions = read_first_accessions(&fasta_path, 5)?;

    Ok(DatabaseInfo {
        id: database_id.to_string(),
        species: entry.species.to_string(),
        db_type: entry.db_type.to_string(),
        path: fasta_path.to_string_lossy().to_string(),
        protein_count: cached.protein_count,
        file_size_bytes: cached.file_size_bytes,
        downloaded_at: cached.downloaded_at.to_rfc3339(),
        sha256: cached.sha256,
        first_accessions,
    })
}

/// Reads the first N accessions from a FASTA file.
fn read_first_accessions(path: &Path, n: usize) -> Result<Vec<String>, FastaDbError> {
    use std::io::{BufRead, BufReader};

    let file = std::fs::File::open(path).map_err(|e| FastaDbError::IoError {
        path: path.to_path_buf(),
        source: e,
    })?;
    let reader = BufReader::new(file);
    let mut accessions = Vec::with_capacity(n);

    for line in reader.lines() {
        if accessions.len() >= n {
            break;
        }
        let line = line.map_err(|e| FastaDbError::IoError {
            path: path.to_path_buf(),
            source: e,
        })?;
        if let Some(header) = line.strip_prefix('>') {
            let acc = header.split_whitespace().next().unwrap_or(header);
            accessions.push(acc.to_string());
        }
    }

    Ok(accessions)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn list_databases_returns_all_with_status() {
        let dir = TempDir::new().unwrap();
        let result = list_databases(dir.path());
        assert!(result.is_ok());
        let list = result.unwrap();
        assert!(list.len() >= 6);
        // All should be Available (nothing downloaded)
        assert!(list
            .iter()
            .all(|d| matches!(d.status, DownloadStatus::Available)));
    }

    #[test]
    fn get_database_info_not_cached_returns_error() {
        let dir = TempDir::new().unwrap();
        let result = get_database_info("human_swissprot", dir.path());
        assert!(result.is_err());
    }
}
