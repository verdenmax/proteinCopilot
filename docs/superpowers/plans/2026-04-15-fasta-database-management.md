# FASTA Database Management Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add built-in FASTA database registry with HTTPS download, local caching, and 3 MCP tools (list_databases, download_database, get_database_info).

**Architecture:** New `crates/fasta-db` lib crate handles registry, download, and cache management. MCP server imports it and exposes 3 tools via `#[rmcp::tool]`. Agent prompts updated to auto-suggest databases by species.

**Tech Stack:** Rust, reqwest (HTTPS streaming), serde_json (registry), sha2 (integrity), tokio (async I/O), rmcp (MCP tools)

---

## File Structure

```
crates/fasta-db/                    ← NEW lib crate
├── Cargo.toml
├── src/
│   ├── lib.rs                      ← Public API: list, download, info
│   ├── registry.rs                 ← Built-in database definitions (static)
│   ├── cache.rs                    ← Cache directory management, registry.json R/W
│   ├── download.rs                 ← HTTPS download with streaming + progress
│   └── error.rs                    ← FastaDbError enum
crates/mcp-server/
├── Cargo.toml                      ← Add fasta-db dependency
├── src/tools.rs                    ← Add 3 tool handlers + Input/Response structs
Cargo.toml                          ← Add workspace dep for fasta-db, reqwest, sha2
.github/prompts/basic-search.prompt.md  ← Add database suggestion guidance
```

---

### Task 1: Scaffold `fasta-db` Crate with Error Types

**Files:**
- Create: `crates/fasta-db/Cargo.toml`
- Create: `crates/fasta-db/src/lib.rs`
- Create: `crates/fasta-db/src/error.rs`
- Modify: `Cargo.toml` (workspace deps)

**Context:** Follow the pattern of existing crates (e.g., `crates/xic/Cargo.toml`). The workspace uses `version.workspace = true` and shared deps.

- [ ] **Step 1: Add workspace dependencies**

Add to root `Cargo.toml` under `[workspace.dependencies]`:

```toml
# HTTP client (for FASTA download)
reqwest = { version = "0.12", default-features = false, features = ["rustls-tls", "stream"] }

# Hashing
sha2 = "0.10"

# Internal
protein-copilot-fasta-db = { path = "crates/fasta-db" }
```

- [ ] **Step 2: Create `crates/fasta-db/Cargo.toml`**

```toml
[package]
name = "protein-copilot-fasta-db"
description = "FASTA database registry, download, and cache management"
version.workspace = true
edition.workspace = true
rust-version.workspace = true
license.workspace = true

[dependencies]
serde = { workspace = true }
serde_json = { workspace = true }
thiserror = { workspace = true }
tokio = { workspace = true }
tracing = { workspace = true }
chrono = { workspace = true }
reqwest = { workspace = true }
sha2 = { workspace = true }
```

- [ ] **Step 3: Create `crates/fasta-db/src/error.rs`**

```rust
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
```

- [ ] **Step 4: Create `crates/fasta-db/src/lib.rs`**

```rust
//! FASTA database registry, download, and cache management.
//!
//! Provides a built-in registry of common proteomics databases (UniProt Swiss-Prot
//! for major model organisms + cRAP contaminants), HTTPS download with streaming,
//! and local file caching with metadata tracking.

pub mod cache;
pub mod download;
pub mod error;
pub mod registry;

pub use error::FastaDbError;
```

- [ ] **Step 5: Verify it compiles**

Run: `cargo check -p protein-copilot-fasta-db`
Expected: Compiles (with unused module warnings, acceptable at this stage)

- [ ] **Step 6: Commit**

```bash
git add crates/fasta-db/ Cargo.toml Cargo.lock
git commit -m "feat(fasta-db): scaffold crate with error types"
```

---

### Task 2: Built-in Database Registry

**Files:**
- Create: `crates/fasta-db/src/registry.rs`

**Context:** Static registry of 6 databases. Each entry has id, species, taxonomy_id, db_type, description, and download URL. URL pattern: `https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(reviewed:true)+AND+(model_organism:{taxonomy_id})`

- [ ] **Step 1: Write the test**

Add to `crates/fasta-db/src/registry.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registry_has_expected_entries() {
        let entries = all_databases();
        assert!(entries.len() >= 6, "should have at least 6 built-in databases");

        let human = get_database("human_swissprot");
        assert!(human.is_some(), "human_swissprot must exist");
        let human = human.unwrap();
        assert_eq!(human.taxonomy_id, 9606);
        assert!(human.url.contains("rest.uniprot.org"));

        let crap = get_database("crap");
        assert!(crap.is_some(), "crap must exist");
    }

    #[test]
    fn unknown_database_returns_none() {
        assert!(get_database("nonexistent_db").is_none());
    }

    #[test]
    fn all_ids_returns_complete_list() {
        let ids = all_database_ids();
        assert!(ids.contains(&"human_swissprot"));
        assert!(ids.contains(&"crap"));
        assert_eq!(ids.len(), all_databases().len());
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p protein-copilot-fasta-db -- registry`
Expected: FAIL — functions not defined

- [ ] **Step 3: Implement the registry**

Write `crates/fasta-db/src/registry.rs`:

```rust
//! Built-in database registry — static definitions of common proteomics FASTA databases.

/// A built-in database entry in the registry.
#[derive(Debug, Clone)]
pub struct DatabaseEntry {
    pub id: &'static str,
    pub species: &'static str,
    pub taxonomy_id: u32,
    pub db_type: &'static str,
    pub description: &'static str,
    pub url: &'static str,
}

const BUILTIN_DATABASES: &[DatabaseEntry] = &[
    DatabaseEntry {
        id: "human_swissprot",
        species: "Homo sapiens",
        taxonomy_id: 9606,
        db_type: "Swiss-Prot",
        description: "Human reviewed proteome (UniProt Swiss-Prot)",
        url: "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(reviewed:true)+AND+(model_organism:9606)",
    },
    DatabaseEntry {
        id: "mouse_swissprot",
        species: "Mus musculus",
        taxonomy_id: 10090,
        db_type: "Swiss-Prot",
        description: "Mouse reviewed proteome (UniProt Swiss-Prot)",
        url: "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(reviewed:true)+AND+(model_organism:10090)",
    },
    DatabaseEntry {
        id: "ecoli_swissprot",
        species: "Escherichia coli (K12)",
        taxonomy_id: 83333,
        db_type: "Swiss-Prot",
        description: "E. coli K12 reviewed proteome (UniProt Swiss-Prot)",
        url: "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(reviewed:true)+AND+(model_organism:83333)",
    },
    DatabaseEntry {
        id: "yeast_swissprot",
        species: "Saccharomyces cerevisiae",
        taxonomy_id: 559292,
        db_type: "Swiss-Prot",
        description: "Yeast reviewed proteome (UniProt Swiss-Prot)",
        url: "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(reviewed:true)+AND+(model_organism:559292)",
    },
    DatabaseEntry {
        id: "arabidopsis_swissprot",
        species: "Arabidopsis thaliana",
        taxonomy_id: 3702,
        db_type: "Swiss-Prot",
        description: "Arabidopsis reviewed proteome (UniProt Swiss-Prot)",
        url: "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(reviewed:true)+AND+(model_organism:3702)",
    },
    DatabaseEntry {
        id: "crap",
        species: "Contaminants",
        taxonomy_id: 0,
        db_type: "cRAP",
        description: "Common Repository of Adventitious Proteins (contaminants)",
        url: "https://ftp.thegpm.org/fasta/cRAP/crap.fasta",
    },
];

/// Returns all built-in database entries.
pub fn all_databases() -> &'static [DatabaseEntry] {
    BUILTIN_DATABASES
}

/// Returns all built-in database IDs.
pub fn all_database_ids() -> Vec<&'static str> {
    BUILTIN_DATABASES.iter().map(|e| e.id).collect()
}

/// Looks up a database entry by ID. Returns `None` if not found.
pub fn get_database(id: &str) -> Option<&'static DatabaseEntry> {
    BUILTIN_DATABASES.iter().find(|e| e.id == id)
}

// tests at bottom...
```

- [ ] **Step 4: Run tests**

Run: `cargo test -p protein-copilot-fasta-db -- registry`
Expected: 3 tests PASS

- [ ] **Step 5: Commit**

```bash
git add crates/fasta-db/src/registry.rs
git commit -m "feat(fasta-db): built-in database registry with 6 entries"
```

---

### Task 3: Cache Management (registry.json)

**Files:**
- Create: `crates/fasta-db/src/cache.rs`

**Context:** Manages the `.proteincopilot/databases/` directory and `registry.json` file. Tracks which databases have been downloaded with metadata.

- [ ] **Step 1: Write the tests**

```rust
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
        // Also create the actual file so is_cached checks file existence
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p protein-copilot-fasta-db -- cache`
Expected: FAIL — types not defined

- [ ] **Step 3: Implement cache.rs**

```rust
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

        let json = serde_json::to_string_pretty(&registry).map_err(|e| {
            FastaDbError::RegistryError {
                detail: format!("serialization error: {e}"),
            }
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

// tests at bottom...
```

- [ ] **Step 4: Run tests**

Run: `cargo test -p protein-copilot-fasta-db -- cache`
Expected: 4 tests PASS

- [ ] **Step 5: Commit**

```bash
git add crates/fasta-db/src/cache.rs
git commit -m "feat(fasta-db): cache manager with registry.json persistence"
```

---

### Task 4: HTTPS Download with Streaming

**Files:**
- Create: `crates/fasta-db/src/download.rs`

**Context:** Downloads FASTA from UniProt REST API via reqwest streaming. Writes to temp file first, then atomically renames. Counts proteins (lines starting with `>`) and computes SHA256 during write.

- [ ] **Step 1: Write the tests**

```rust
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
        assert_eq!(hash1.len(), 64); // hex string
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p protein-copilot-fasta-db -- download`
Expected: FAIL — functions not defined

- [ ] **Step 3: Implement download.rs**

```rust
//! HTTPS download with streaming, SHA256 hashing, and protein counting.

use std::path::Path;

use sha2::{Digest, Sha256};
use tokio::io::AsyncWriteExt;
use tracing;

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

    let response = reqwest::get(url).await.map_err(|e| FastaDbError::NetworkError {
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

    let bytes = response.bytes().await.map_err(|e| FastaDbError::NetworkError {
        url: url.to_string(),
        detail: format!("failed to read response body: {e}"),
    })?;

    // Write to .part file
    let mut file = tokio::fs::File::create(&part_path)
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
    content.iter().enumerate().filter(|(i, &b)| {
        b == b'>' && (*i == 0 || content[*i - 1] == b'\n')
    }).count() as u64
}

/// Computes the SHA256 hex digest of a byte slice.
pub(crate) fn compute_sha256(data: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(data);
    format!("{:x}", hasher.finalize())
}

// tests at bottom...
```

- [ ] **Step 4: Run tests**

Run: `cargo test -p protein-copilot-fasta-db -- download`
Expected: 3 tests PASS

- [ ] **Step 5: Commit**

```bash
git add crates/fasta-db/src/download.rs
git commit -m "feat(fasta-db): HTTPS download with SHA256 and protein counting"
```

---

### Task 5: Public API (lib.rs Functions)

**Files:**
- Modify: `crates/fasta-db/src/lib.rs`

**Context:** Wire up the public API that MCP tools will call: `list_databases()`, `download_database()`, `get_database_info()`. These compose registry + cache + download.

- [ ] **Step 1: Write tests**

Add to `crates/fasta-db/src/lib.rs`:

```rust
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
        assert!(list.iter().all(|d| matches!(d.status, DownloadStatus::Available)));
    }

    #[test]
    fn get_database_info_not_cached_returns_error() {
        let dir = TempDir::new().unwrap();
        let result = get_database_info("human_swissprot", dir.path());
        assert!(result.is_err());
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p protein-copilot-fasta-db -- tests`
Expected: FAIL — functions not defined

- [ ] **Step 3: Implement public API in lib.rs**

```rust
//! FASTA database registry, download, and cache management.

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
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatabaseStatus {
    pub id: String,
    pub species: String,
    pub db_type: String,
    pub description: String,
    pub status: DownloadStatus,
}

/// Whether a database is available for download or already cached.
#[derive(Debug, Clone, Serialize, Deserialize)]
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
#[derive(Debug, Clone, Serialize, Deserialize)]
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
#[derive(Debug, Clone, Serialize, Deserialize)]
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
        available: registry::all_database_ids().into_iter().map(String::from).collect(),
    })?;

    let cache = CacheManager::new(cache_dir.to_path_buf());

    // Return cached path if already downloaded and not forcing
    if !force && cache.is_cached(database_id) {
        let path = cache.fasta_path(database_id);
        let cached = cache.get_cached(database_id)?.unwrap();
        return Ok(DownloadDatabaseResult {
            id: database_id.to_string(),
            path: path.to_string_lossy().to_string(),
            protein_count: cached.protein_count,
            file_size_bytes: cached.file_size_bytes,
        });
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
pub fn get_database_info(database_id: &str, cache_dir: &Path) -> Result<DatabaseInfo, FastaDbError> {
    let entry = registry::get_database(database_id).ok_or_else(|| FastaDbError::UnknownDatabase {
        id: database_id.to_string(),
        available: registry::all_database_ids().into_iter().map(String::from).collect(),
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

// tests at bottom...
```

- [ ] **Step 4: Run tests**

Run: `cargo test -p protein-copilot-fasta-db`
Expected: All tests PASS (registry + cache + download + lib tests)

- [ ] **Step 5: Commit**

```bash
git add crates/fasta-db/src/lib.rs
git commit -m "feat(fasta-db): public API — list, download, info"
```

---

### Task 6: MCP Tools Integration

**Files:**
- Modify: `crates/mcp-server/Cargo.toml` (add fasta-db dep)
- Modify: `crates/mcp-server/src/tools.rs` (add 3 tools + Input/Response structs)

**Context:** Follow the pattern of existing tools like `list_presets` (simple, no input) and `run_search` (async with Parameters). Use `#[rmcp::tool]` attribute macro. Input structs derive `Deserialize + JsonSchema`. Response types derive `Serialize + JsonSchema`. Cache dir defaults to `.proteincopilot/databases/` relative to CWD.

- [ ] **Step 1: Add dependency to mcp-server Cargo.toml**

Add to `[dependencies]` in `crates/mcp-server/Cargo.toml`:

```toml
protein-copilot-fasta-db = { workspace = true }
```

- [ ] **Step 2: Add Input/Response structs to tools.rs**

Add near the other Input structs (around line 460):

```rust
// --- FASTA Database Management ---

#[derive(Debug, Deserialize, JsonSchema)]
struct ListDatabasesInput {
    /// Override cache directory. Default: .proteincopilot/databases/
    #[serde(default)]
    cache_dir: Option<String>,
}

#[derive(Debug, Deserialize, JsonSchema)]
struct DownloadDatabaseInput {
    /// Database ID (e.g. "human_swissprot", "mouse_swissprot", "crap")
    database_id: String,
    /// Override cache directory. Default: .proteincopilot/databases/
    #[serde(default)]
    cache_dir: Option<String>,
    /// Force re-download even if already cached
    #[serde(default)]
    force: Option<bool>,
}

#[derive(Debug, Deserialize, JsonSchema)]
struct GetDatabaseInfoInput {
    /// Database ID to get info for
    database_id: String,
    /// Override cache directory. Default: .proteincopilot/databases/
    #[serde(default)]
    cache_dir: Option<String>,
}
```

- [ ] **Step 3: Add helper for default cache dir**

```rust
fn default_cache_dir(override_dir: &Option<String>) -> std::path::PathBuf {
    if let Some(ref dir) = override_dir {
        std::path::PathBuf::from(dir)
    } else {
        std::path::PathBuf::from(".proteincopilot/databases")
    }
}
```

- [ ] **Step 4: Add the 3 tool handlers**

Add inside the `#[rmcp::tool]` impl block alongside existing tools:

```rust
/// List all available FASTA databases and their cache status.
#[rmcp::tool(
    name = "list_databases",
    description = "List all built-in FASTA protein databases (Human, Mouse, E.coli, Yeast, Arabidopsis, cRAP contaminants) with download status. Shows which databases are cached locally and which are available for download."
)]
fn list_databases(
    &self,
    Parameters(input): Parameters<ListDatabasesInput>,
) -> Result<Json<Vec<protein_copilot_fasta_db::DatabaseStatus>>, ErrorData> {
    let cache_dir = default_cache_dir(&input.cache_dir);
    protein_copilot_fasta_db::list_databases(&cache_dir)
        .map(Json)
        .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))
}

/// Download a FASTA database by ID.
#[rmcp::tool(
    name = "download_database",
    description = "Download a FASTA protein database by ID (e.g. 'human_swissprot', 'mouse_swissprot', 'ecoli_swissprot', 'yeast_swissprot', 'arabidopsis_swissprot', 'crap'). Downloads from UniProt via HTTPS and caches locally. Returns the local file path for use as database_path in search parameters. Use list_databases first to see available options."
)]
async fn download_database(
    &self,
    Parameters(input): Parameters<DownloadDatabaseInput>,
) -> Result<Json<protein_copilot_fasta_db::DownloadDatabaseResult>, ErrorData> {
    let cache_dir = default_cache_dir(&input.cache_dir);
    let force = input.force.unwrap_or(false);
    protein_copilot_fasta_db::download_database(&input.database_id, &cache_dir, force)
        .await
        .map(Json)
        .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))
}

/// Get detailed info about a cached FASTA database.
#[rmcp::tool(
    name = "get_database_info",
    description = "Get detailed information about a downloaded FASTA database: protein count, file size, SHA256 hash, download date, and first 5 protein accessions. The database must be downloaded first using download_database."
)]
fn get_database_info(
    &self,
    Parameters(input): Parameters<GetDatabaseInfoInput>,
) -> Result<Json<protein_copilot_fasta_db::DatabaseInfo>, ErrorData> {
    let cache_dir = default_cache_dir(&input.cache_dir);
    protein_copilot_fasta_db::get_database_info(&input.database_id, &cache_dir)
        .map(Json)
        .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e))
}
```

- [ ] **Step 5: Verify compilation**

Run: `cargo check -p protein-copilot-mcp-server`
Expected: Compiles without errors

- [ ] **Step 6: Run full workspace tests**

Run: `cargo test --workspace`
Expected: All tests pass (existing + new fasta-db tests)

- [ ] **Step 7: Commit**

```bash
git add crates/mcp-server/Cargo.toml crates/mcp-server/src/tools.rs
git commit -m "feat: add list_databases, download_database, get_database_info MCP tools"
```

---

### Task 7: Agent Prompt Update

**Files:**
- Modify: `.github/prompts/basic-search.prompt.md`

**Context:** Add guidance for the Agent to auto-suggest databases when users mention species names. The Agent should call `list_databases` to check availability, then suggest `download_database` if needed.

- [ ] **Step 1: Add database suggestion section**

Add a new section to `basic-search.prompt.md` (after the parameter recommendation step, before search execution):

```markdown
## Step 2.5: Database Selection

If the user hasn't specified a FASTA database path:

1. **Detect species intent** from user's message:
   - "人"/"human"/"人类" → `human_swissprot`
   - "小鼠"/"mouse" → `mouse_swissprot`
   - "大肠杆菌"/"E.coli" → `ecoli_swissprot`
   - "酵母"/"yeast" → `yeast_swissprot`
   - "拟南芥"/"Arabidopsis" → `arabidopsis_swissprot`

2. **Check availability**: Call `list_databases` to see if the database is cached locally.

3. **If not cached**: Suggest downloading:
   > "检测到您需要搜索 [物种] 蛋白质组，需要下载 UniProt Swiss-Prot 数据库（约 X MB）。是否下载？"

4. **If cached**: Use the cached path directly as `database_path`.

5. **Contaminants**: Always mention cRAP contaminant database availability. If user wants to include contaminants, download `crap` separately and inform user of both paths.

6. **Set database_path**: After confirming, set `SearchParams.database_path` to the downloaded FASTA path.
```

- [ ] **Step 2: Verify prompt file is valid**

Run: `cat .github/prompts/basic-search.prompt.md | head -5`
Expected: File exists and has content

- [ ] **Step 3: Commit**

```bash
git add .github/prompts/basic-search.prompt.md
git commit -m "docs: add database auto-suggestion to basic-search prompt"
```

---

### Task 8: Integration Test (Real Download)

**Files:**
- Create: `crates/fasta-db/tests/integration.rs`

**Context:** An integration test that actually downloads a small database (E.coli ~2 MB) to verify the full pipeline works. This test is `#[ignore]` by default (requires network).

- [ ] **Step 1: Create the integration test**

```rust
//! Integration tests for fasta-db (requires network access).
//! Run with: cargo test -p protein-copilot-fasta-db -- --ignored

use protein_copilot_fasta_db::{download_database, get_database_info, list_databases};
use tempfile::TempDir;

#[tokio::test]
#[ignore = "requires network access"]
async fn download_ecoli_end_to_end() {
    let dir = TempDir::new().unwrap();
    let cache_dir = dir.path();

    // Step 1: List — ecoli should be Available
    let list = list_databases(cache_dir).unwrap();
    let ecoli = list.iter().find(|d| d.id == "ecoli_swissprot").unwrap();
    assert!(
        matches!(ecoli.status, protein_copilot_fasta_db::DownloadStatus::Available),
        "ecoli should not be cached yet"
    );

    // Step 2: Download
    let result = download_database("ecoli_swissprot", cache_dir, false)
        .await
        .unwrap();
    assert!(result.protein_count > 100, "E.coli should have >100 proteins");
    assert!(result.file_size_bytes > 0);

    // Step 3: List — ecoli should now be Downloaded
    let list = list_databases(cache_dir).unwrap();
    let ecoli = list.iter().find(|d| d.id == "ecoli_swissprot").unwrap();
    assert!(
        matches!(ecoli.status, protein_copilot_fasta_db::DownloadStatus::Downloaded { .. }),
        "ecoli should be cached now"
    );

    // Step 4: Info
    let info = get_database_info("ecoli_swissprot", cache_dir).unwrap();
    assert_eq!(info.protein_count, result.protein_count);
    assert!(!info.first_accessions.is_empty());
    assert!(info.sha256.len() == 64);

    // Step 5: Re-download with force=false — should return cached
    let result2 = download_database("ecoli_swissprot", cache_dir, false)
        .await
        .unwrap();
    assert_eq!(result2.protein_count, result.protein_count);

    // Step 6: Re-download with force=true — should re-download
    let result3 = download_database("ecoli_swissprot", cache_dir, true)
        .await
        .unwrap();
    assert!(result3.protein_count > 0);
}

#[tokio::test]
#[ignore = "requires network access"]
async fn unknown_database_returns_error() {
    let dir = TempDir::new().unwrap();
    let result = download_database("nonexistent_db", dir.path(), false).await;
    assert!(result.is_err());
    let err = result.unwrap_err().to_string();
    assert!(err.contains("unknown database"));
    assert!(err.contains("human_swissprot")); // should list available
}
```

- [ ] **Step 2: Run integration test (network required)**

Run: `cargo test -p protein-copilot-fasta-db -- --ignored --nocapture`
Expected: Tests PASS (downloads ~2 MB E.coli database)

- [ ] **Step 3: Run full workspace tests (unit only)**

Run: `cargo test --workspace`
Expected: All unit tests pass, integration tests skipped (ignored)

- [ ] **Step 4: Final clippy check**

Run: `cargo clippy --all-targets -- -D warnings`
Expected: Clean

- [ ] **Step 5: Commit**

```bash
git add crates/fasta-db/tests/integration.rs
git commit -m "test(fasta-db): add end-to-end integration test with real download"
```
