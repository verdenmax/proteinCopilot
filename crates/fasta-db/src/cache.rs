//! Cache directory management and registry.json read/write.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};

use crate::FastaDbError;

/// Metadata for a single cached (downloaded) database file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CachedDatabase {
    pub id: String,
    pub file_name: String,
    pub downloaded_at: DateTime<Utc>,
    pub file_size_bytes: u64,
    pub protein_count: u64,
    pub sha256: String,
}

/// On-disk registry tracking all cached databases.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CacheRegistry {
    pub version: u32,
    pub databases: HashMap<String, CachedDatabase>,
}

impl Default for CacheRegistry {
    fn default() -> Self {
        Self {
            version: 1,
            databases: HashMap::new(),
        }
    }
}

/// Manages the local cache directory and registry.json.
pub struct CacheManager {
    cache_dir: PathBuf,
}

impl CacheManager {
    pub fn new(cache_dir: PathBuf) -> Self {
        Self { cache_dir }
    }

    /// Returns the path to `registry.json`.
    fn registry_path(&self) -> PathBuf {
        self.cache_dir.join("registry.json")
    }

    /// Returns the expected FASTA file path for a given database ID.
    pub fn fasta_path(&self, database_id: &str) -> PathBuf {
        self.cache_dir.join(format!("{}.fasta", database_id))
    }

    /// Returns the cache directory path.
    pub fn cache_dir(&self) -> &Path {
        &self.cache_dir
    }

    /// Loads the registry from disk. Returns empty registry if file doesn't exist.
    pub fn load_registry(&self) -> Result<CacheRegistry, FastaDbError> {
        let path = self.registry_path();
        if !path.exists() {
            return Ok(CacheRegistry::default());
        }
        let content = std::fs::read_to_string(&path).map_err(|e| FastaDbError::IoError {
            path: path.clone(),
            source: e,
        })?;
        serde_json::from_str(&content).map_err(|e| FastaDbError::RegistryError {
            detail: format!("failed to parse {}: {}", path.display(), e),
        })
    }

    /// Saves/updates a single entry in the registry (read-modify-write).
    pub fn save_entry(&self, entry: &CachedDatabase) -> Result<(), FastaDbError> {
        std::fs::create_dir_all(&self.cache_dir).map_err(|e| FastaDbError::IoError {
            path: self.cache_dir.clone(),
            source: e,
        })?;

        let mut registry = self.load_registry()?;
        registry.databases.insert(entry.id.clone(), entry.clone());

        let json =
            serde_json::to_string_pretty(&registry).map_err(|e| FastaDbError::RegistryError {
                detail: format!("serialization error: {e}"),
            })?;
        std::fs::write(self.registry_path(), json).map_err(|e| FastaDbError::IoError {
            path: self.registry_path(),
            source: e,
        })
    }

    /// Checks whether a database is cached (entry in registry AND file exists on disk).
    pub fn is_cached(&self, database_id: &str) -> bool {
        let Ok(registry) = self.load_registry() else {
            return false;
        };
        registry.databases.contains_key(database_id) && self.fasta_path(database_id).exists()
    }

    /// Returns the cached metadata for a database, if available.
    pub fn get_cached(&self, database_id: &str) -> Result<Option<CachedDatabase>, FastaDbError> {
        let registry = self.load_registry()?;
        Ok(registry.databases.get(database_id).cloned())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn empty_cache_returns_empty_registry() {
        let dir = TempDir::new().unwrap();
        let cache = CacheManager::new(dir.path().to_path_buf());
        let reg = cache.load_registry().unwrap();
        assert!(reg.databases.is_empty());
    }

    #[test]
    fn save_and_load_round_trips() {
        let dir = TempDir::new().unwrap();
        let cache = CacheManager::new(dir.path().to_path_buf());

        let entry = CachedDatabase {
            id: "human_swissprot".to_string(),
            file_name: "human_swissprot.fasta".to_string(),
            downloaded_at: chrono::Utc::now(),
            file_size_bytes: 12345678,
            protein_count: 20422,
            sha256: "abc123".to_string(),
        };
        cache.save_entry(&entry).unwrap();

        let reg = cache.load_registry().unwrap();
        assert_eq!(reg.databases.len(), 1);
        assert_eq!(reg.databases["human_swissprot"].protein_count, 20422);
    }

    #[test]
    fn is_cached_returns_correct_status() {
        let dir = TempDir::new().unwrap();
        let cache = CacheManager::new(dir.path().to_path_buf());
        assert!(!cache.is_cached("human_swissprot"));

        let entry = CachedDatabase {
            id: "human_swissprot".to_string(),
            file_name: "human_swissprot.fasta".to_string(),
            downloaded_at: chrono::Utc::now(),
            file_size_bytes: 100,
            protein_count: 5,
            sha256: "test".to_string(),
        };
        std::fs::create_dir_all(dir.path()).unwrap();
        std::fs::write(dir.path().join("human_swissprot.fasta"), b">P1\nACDE\n").unwrap();
        cache.save_entry(&entry).unwrap();

        assert!(cache.is_cached("human_swissprot"));
    }

    #[test]
    fn fasta_path_returns_expected_location() {
        let dir = TempDir::new().unwrap();
        let cache = CacheManager::new(dir.path().to_path_buf());
        let path = cache.fasta_path("human_swissprot");
        assert!(path.ends_with("human_swissprot.fasta"));
        assert!(path.starts_with(dir.path()));
    }
}
