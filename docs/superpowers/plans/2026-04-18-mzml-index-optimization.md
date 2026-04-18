# mzML Index Optimization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate MCP Client timeouts on large mzML files by adding disk-persisted index caching, SIMD-accelerated fallback scanning, and timeout configuration documentation.

**Architecture:** Three-layer index resolution in `IndexedMzMLReader::open()`: disk cache (.mzml.idx) → native mzML `<indexList>` → accelerated byte-scan fallback. After any non-cached build, persist to disk. Add `memchr` dependency for SIMD byte search. Update `.mcp.json` with timeout field.

**Tech Stack:** Rust, `memchr` crate (SIMD byte search), binary I/O (`std::io::{Read, Write, Seek}`)

---

## File Map

| Action | File | Responsibility |
|--------|------|----------------|
| Create | `crates/spectrum-io/src/disk_cache.rs` | Binary .idx file read/write (PCIX format) |
| Modify | `crates/spectrum-io/src/index.rs` | Add `build_index_by_byte_scan()` alongside existing line-scan |
| Modify | `crates/spectrum-io/src/indexed_mzml.rs` | Three-layer open() with disk cache persist |
| Modify | `crates/spectrum-io/src/lib.rs` | Export `disk_cache` module |
| Modify | `crates/spectrum-io/src/error.rs` | Add `DiskCacheError` variant |
| Modify | `crates/spectrum-io/Cargo.toml` | Add `memchr` dependency |
| Modify | `Cargo.toml` (workspace) | Add `memchr` to workspace dependencies |
| Modify | `.mcp.json` | Add `timeout: 300` |
| Modify | `README.md` | Add timeout config section |

---

### Task 1: Add `memchr` dependency

**Files:**
- Modify: `Cargo.toml` (workspace root, line ~35)
- Modify: `crates/spectrum-io/Cargo.toml` (line ~14)

- [ ] **Step 1: Add `memchr` to workspace dependencies**

In `Cargo.toml` (workspace root), add under `[workspace.dependencies]` after the `flate2` line:

```toml
memchr = "2"
```

- [ ] **Step 2: Add `memchr` to spectrum-io dependencies**

In `crates/spectrum-io/Cargo.toml`, add under `[dependencies]` after `flate2`:

```toml
memchr = { workspace = true }
```

- [ ] **Step 3: Verify it compiles**

Run: `cargo check -p protein-copilot-spectrum-io`
Expected: compiles successfully with `memchr` available

- [ ] **Step 4: Commit**

```bash
git add Cargo.toml crates/spectrum-io/Cargo.toml
git commit -m "chore(spectrum-io): add memchr dependency for SIMD byte scanning"
```

---

### Task 2: Add `DiskCacheError` variant to error.rs

**Files:**
- Modify: `crates/spectrum-io/src/error.rs`

- [ ] **Step 1: Add `DiskCacheError` variant**

In `crates/spectrum-io/src/error.rs`, add a new variant to the `SpectrumIoError` enum after the `IndexParseError` variant (before the closing `}`):

```rust
    /// Disk cache (.idx file) I/O error. Non-fatal: used for logging.
    #[error("disk cache error for {path}: {detail}")]
    DiskCacheError {
        /// The mzML file path (not the .idx path).
        path: PathBuf,
        /// What went wrong.
        detail: String,
    },
```

- [ ] **Step 2: Add the `From<SpectrumIoError> for CoreError` match arm**

In the same file, add a match arm in the `impl From<SpectrumIoError> for CoreError` block, after the `IndexParseError` arm:

```rust
            SpectrumIoError::DiskCacheError { path, detail } => {
                protein_copilot_core::error::CoreError::SpectrumParseError {
                    format: "mzML".to_string(),
                    detail: format!("{}: {detail}", path.display()),
                    suggestion: "Delete the .idx cache file and retry".to_string(),
                }
            }
```

- [ ] **Step 3: Verify it compiles**

Run: `cargo check -p protein-copilot-spectrum-io`
Expected: compiles successfully

- [ ] **Step 4: Commit**

```bash
git add crates/spectrum-io/src/error.rs
git commit -m "feat(spectrum-io): add DiskCacheError variant for .idx file errors"
```

---

### Task 3: Implement `disk_cache.rs` — binary .idx read/write

**Files:**
- Create: `crates/spectrum-io/src/disk_cache.rs`
- Modify: `crates/spectrum-io/src/lib.rs` (add module export)

- [ ] **Step 1: Write the failing tests first**

Create `crates/spectrum-io/src/disk_cache.rs` with the test module and all test functions but empty implementations that return `todo!()`:

```rust
//! Disk-persisted scan index cache using the PCIX binary format.
//!
//! Format (little-endian):
//! ```text
//! [magic: 4 bytes "PCIX"]
//! [version: u8 = 1]
//! [source_file_size: u64]
//! [source_file_mtime: u64]
//! [entry_count: u32]
//! [entries: (scan_number: u32, byte_offset: u64) × entry_count]
//! ```

use std::collections::HashMap;
use std::io::{Read, Write};
use std::path::Path;

use crate::error::SpectrumIoError;
use crate::index::{IndexSource, ScanIndex};

const MAGIC: &[u8; 4] = b"PCIX";
const VERSION: u8 = 1;
const HEADER_SIZE: usize = 4 + 1 + 8 + 8 + 4; // 25 bytes
const ENTRY_SIZE: usize = 4 + 8; // 12 bytes per entry

/// Returns the sidecar .idx path for a given mzML path.
///
/// Example: `/data/experiment.mzML` → `/data/experiment.mzML.idx`
pub fn idx_path(mzml_path: &Path) -> std::path::PathBuf {
    let mut p = mzml_path.as_os_str().to_owned();
    p.push(".idx");
    std::path::PathBuf::from(p)
}

/// Attempts to load a cached ScanIndex from the `.idx` sidecar file.
///
/// Returns `Ok(Some(index))` if the cache exists and matches the source file's
/// size and modification time. Returns `Ok(None)` if the cache is missing or
/// stale. Returns `Err` only on I/O errors reading the cache file.
pub fn load_index(
    mzml_path: &Path,
    expected_size: u64,
    expected_mtime: u64,
) -> Result<Option<ScanIndex>, SpectrumIoError> {
    todo!()
}

/// Persists a ScanIndex to the `.idx` sidecar file.
///
/// If writing fails (e.g. read-only filesystem), logs a warning and returns
/// the error — callers should treat this as non-fatal.
pub fn save_index(
    mzml_path: &Path,
    index: &ScanIndex,
    file_size: u64,
    file_mtime: u64,
) -> Result<(), SpectrumIoError> {
    todo!()
}

/// Reads the source file's size and modification time (Unix epoch seconds).
pub fn file_metadata(path: &Path) -> Result<(u64, u64), SpectrumIoError> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    fn make_index(scans: &[(u32, u64)]) -> ScanIndex {
        let offsets: HashMap<u32, u64> = scans.iter().copied().collect();
        ScanIndex::new(offsets, IndexSource::BuiltFromScan)
    }

    #[test]
    fn round_trip_save_load() {
        let dir = tempfile::tempdir().unwrap();
        let mzml = dir.path().join("test.mzml");
        std::fs::write(&mzml, b"fake mzml content for size").unwrap();
        let (size, mtime) = file_metadata(&mzml).unwrap();

        let original = make_index(&[(1, 100), (5, 500), (10, 1000)]);
        save_index(&mzml, &original, size, mtime).unwrap();

        let loaded = load_index(&mzml, size, mtime).unwrap().expect("should load");
        assert_eq!(loaded.len(), 3);
        assert_eq!(loaded.get_offset(1), Some(100));
        assert_eq!(loaded.get_offset(5), Some(500));
        assert_eq!(loaded.get_offset(10), Some(1000));
    }

    #[test]
    fn stale_size_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let mzml = dir.path().join("test.mzml");
        std::fs::write(&mzml, b"original").unwrap();
        let (size, mtime) = file_metadata(&mzml).unwrap();

        let idx = make_index(&[(1, 0)]);
        save_index(&mzml, &idx, size, mtime).unwrap();

        // File size changed
        let result = load_index(&mzml, size + 999, mtime).unwrap();
        assert!(result.is_none(), "stale size should return None");
    }

    #[test]
    fn stale_mtime_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let mzml = dir.path().join("test.mzml");
        std::fs::write(&mzml, b"content").unwrap();
        let (size, mtime) = file_metadata(&mzml).unwrap();

        let idx = make_index(&[(1, 0)]);
        save_index(&mzml, &idx, size, mtime).unwrap();

        // mtime changed
        let result = load_index(&mzml, size, mtime + 1).unwrap();
        assert!(result.is_none(), "stale mtime should return None");
    }

    #[test]
    fn missing_idx_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let mzml = dir.path().join("nonexistent.mzml");
        let result = load_index(&mzml, 0, 0).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn corrupt_magic_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let mzml = dir.path().join("test.mzml");
        std::fs::write(&mzml, b"content").unwrap();

        let idx_file = idx_path(&mzml);
        std::fs::write(&idx_file, b"BAAD_not_pcix_magic").unwrap();

        let result = load_index(&mzml, 7, 0).unwrap();
        assert!(result.is_none(), "corrupt magic should return None");
    }

    #[test]
    fn truncated_file_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let mzml = dir.path().join("test.mzml");
        std::fs::write(&mzml, b"content").unwrap();

        let idx_file = idx_path(&mzml);
        // Write valid magic but truncated header
        std::fs::write(&idx_file, b"PCIX").unwrap();

        let result = load_index(&mzml, 7, 0).unwrap();
        assert!(result.is_none(), "truncated file should return None");
    }

    #[test]
    fn empty_index_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let mzml = dir.path().join("test.mzml");
        std::fs::write(&mzml, b"x").unwrap();
        let (size, mtime) = file_metadata(&mzml).unwrap();

        let original = make_index(&[]);
        save_index(&mzml, &original, size, mtime).unwrap();

        let loaded = load_index(&mzml, size, mtime).unwrap().expect("should load empty");
        assert_eq!(loaded.len(), 0);
        assert!(loaded.is_empty());
    }

    #[test]
    fn large_offsets_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let mzml = dir.path().join("test.mzml");
        std::fs::write(&mzml, b"x").unwrap();
        let (size, mtime) = file_metadata(&mzml).unwrap();

        // Use offsets >4GB to test u64 handling
        let original = make_index(&[(1, 5_000_000_000), (2, 8_000_000_000)]);
        save_index(&mzml, &original, size, mtime).unwrap();

        let loaded = load_index(&mzml, size, mtime).unwrap().expect("should load");
        assert_eq!(loaded.get_offset(1), Some(5_000_000_000));
        assert_eq!(loaded.get_offset(2), Some(8_000_000_000));
    }

    #[test]
    fn idx_path_appends_extension() {
        let p = idx_path(Path::new("/data/sample.mzML"));
        assert_eq!(p, std::path::PathBuf::from("/data/sample.mzML.idx"));
    }
}
```

- [ ] **Step 2: Add module declaration to lib.rs**

In `crates/spectrum-io/src/lib.rs`, add after `mod util;` (line 30):

```rust
pub mod disk_cache;
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `cargo test -p protein-copilot-spectrum-io disk_cache -- --no-capture 2>&1 | head -30`
Expected: FAIL — all tests hit `todo!()`

- [ ] **Step 4: Implement `file_metadata()`**

Replace the `todo!()` in `file_metadata` with:

```rust
pub fn file_metadata(path: &Path) -> Result<(u64, u64), SpectrumIoError> {
    let meta = std::fs::metadata(path).map_err(|e| SpectrumIoError::DiskCacheError {
        path: path.to_path_buf(),
        detail: format!("cannot read metadata: {e}"),
    })?;
    let size = meta.len();
    let mtime = meta
        .modified()
        .map_err(|e| SpectrumIoError::DiskCacheError {
            path: path.to_path_buf(),
            detail: format!("cannot read mtime: {e}"),
        })?
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    Ok((size, mtime))
}
```

- [ ] **Step 5: Implement `save_index()`**

Replace the `todo!()` in `save_index` with:

```rust
pub fn save_index(
    mzml_path: &Path,
    index: &ScanIndex,
    file_size: u64,
    file_mtime: u64,
) -> Result<(), SpectrumIoError> {
    let out_path = idx_path(mzml_path);
    let scans = index.scan_numbers();
    let entry_count = scans.len() as u32;

    let mut buf = Vec::with_capacity(HEADER_SIZE + scans.len() * ENTRY_SIZE);
    buf.extend_from_slice(MAGIC);
    buf.push(VERSION);
    buf.extend_from_slice(&file_size.to_le_bytes());
    buf.extend_from_slice(&file_mtime.to_le_bytes());
    buf.extend_from_slice(&entry_count.to_le_bytes());

    for &scan in &scans {
        let offset = index.get_offset(scan).unwrap_or(0);
        buf.extend_from_slice(&scan.to_le_bytes());
        buf.extend_from_slice(&offset.to_le_bytes());
    }

    std::fs::write(&out_path, &buf).map_err(|e| SpectrumIoError::DiskCacheError {
        path: mzml_path.to_path_buf(),
        detail: format!("failed to write {}: {e}", out_path.display()),
    })?;

    tracing::info!(
        path = %out_path.display(),
        entries = entry_count,
        "saved disk index cache"
    );
    Ok(())
}
```

- [ ] **Step 6: Implement `load_index()`**

Replace the `todo!()` in `load_index` with:

```rust
pub fn load_index(
    mzml_path: &Path,
    expected_size: u64,
    expected_mtime: u64,
) -> Result<Option<ScanIndex>, SpectrumIoError> {
    let cache_path = idx_path(mzml_path);

    let data = match std::fs::read(&cache_path) {
        Ok(d) => d,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(e) => {
            tracing::warn!(
                path = %cache_path.display(),
                error = %e,
                "failed to read disk index cache"
            );
            return Ok(None);
        }
    };

    // Validate minimum size
    if data.len() < HEADER_SIZE {
        tracing::warn!(path = %cache_path.display(), "disk cache too small, ignoring");
        return Ok(None);
    }

    // Check magic
    if &data[0..4] != MAGIC {
        tracing::warn!(path = %cache_path.display(), "invalid magic in disk cache, ignoring");
        return Ok(None);
    }

    // Check version
    let version = data[4];
    if version != VERSION {
        tracing::warn!(
            path = %cache_path.display(),
            version,
            "unsupported disk cache version, ignoring"
        );
        return Ok(None);
    }

    // Read and check source file metadata
    let stored_size = u64::from_le_bytes(data[5..13].try_into().unwrap());
    let stored_mtime = u64::from_le_bytes(data[13..21].try_into().unwrap());

    if stored_size != expected_size || stored_mtime != expected_mtime {
        tracing::info!(
            path = %cache_path.display(),
            "disk cache stale (size/mtime mismatch), rebuilding"
        );
        return Ok(None);
    }

    let entry_count = u32::from_le_bytes(data[21..25].try_into().unwrap()) as usize;
    let expected_len = HEADER_SIZE + entry_count * ENTRY_SIZE;

    if data.len() < expected_len {
        tracing::warn!(
            path = %cache_path.display(),
            expected = expected_len,
            actual = data.len(),
            "disk cache truncated, ignoring"
        );
        return Ok(None);
    }

    let mut offsets = HashMap::with_capacity(entry_count);
    for i in 0..entry_count {
        let base = HEADER_SIZE + i * ENTRY_SIZE;
        let scan = u32::from_le_bytes(data[base..base + 4].try_into().unwrap());
        let offset = u64::from_le_bytes(data[base + 4..base + 12].try_into().unwrap());
        offsets.insert(scan, offset);
    }

    tracing::info!(
        path = %cache_path.display(),
        entries = entry_count,
        "loaded disk index cache"
    );
    // Disk-cached indexes were originally built by some method; mark as
    // NativeIndex so no further disk-save is triggered if the caller checks.
    Ok(Some(ScanIndex::new(offsets, IndexSource::NativeIndex)))
}
```

- [ ] **Step 7: Run tests to verify they all pass**

Run: `cargo test -p protein-copilot-spectrum-io disk_cache -- --no-capture`
Expected: 9 tests pass

- [ ] **Step 8: Commit**

```bash
git add crates/spectrum-io/src/disk_cache.rs crates/spectrum-io/src/lib.rs
git commit -m "feat(spectrum-io): add disk-persisted scan index cache (PCIX format)

Binary .idx sidecar files enable O(1) index load across server restarts.
Staleness detection via source file size + mtime comparison."
```

---

### Task 4: Implement `build_index_by_byte_scan()` in index.rs

**Files:**
- Modify: `crates/spectrum-io/src/index.rs`

- [ ] **Step 1: Write the failing test**

Add to the `mod tests` block at the bottom of `crates/spectrum-io/src/index.rs`:

```rust
    #[test]
    fn byte_scan_matches_line_scan() {
        let path =
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/small.mzml");
        let line_idx = build_index_by_scanning(&path).unwrap();
        let byte_idx = build_index_by_byte_scan(&path).unwrap();

        assert_eq!(byte_idx.len(), line_idx.len());
        assert_eq!(byte_idx.source(), IndexSource::BuiltFromScan);
        for scan in line_idx.scan_numbers() {
            assert_eq!(
                byte_idx.get_offset(scan),
                line_idx.get_offset(scan),
                "offset mismatch for scan {scan}"
            );
        }
    }

    #[test]
    fn byte_scan_indexed_mzml_matches() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/small_indexed.mzml");
        if !path.exists() {
            generate_indexed_fixture(&path);
        }
        let byte_idx = build_index_by_byte_scan(&path).unwrap();
        assert_eq!(byte_idx.len(), 10);
        for scan in 1..=10 {
            assert!(byte_idx.get_offset(scan).is_some(), "missing scan {scan}");
        }
    }
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p protein-copilot-spectrum-io byte_scan -- --no-capture 2>&1 | head -20`
Expected: FAIL — `build_index_by_byte_scan` not found

- [ ] **Step 3: Implement `extract_scan_from_tag_bytes()` helper**

Add above the existing `extract_id_attr` function in `crates/spectrum-io/src/index.rs`:

```rust
/// Extracts scan number from a `<spectrum ...>` tag in raw bytes.
///
/// Searches for ` id="...scan=N..."` within the byte slice and parses the
/// scan number. Falls back to `fallback_scan` if no id attribute is found.
fn extract_scan_from_tag_bytes(tag_bytes: &[u8], fallback_scan: u32) -> u32 {
    // Search for ` id="` in the tag bytes (must have leading space to avoid "nativeID=")
    let needle = b" id=\"";
    let id_start = match memchr::memmem::find(tag_bytes, needle) {
        Some(pos) => pos + needle.len(),
        None => {
            // Try single quotes
            let needle_sq = b" id='";
            match memchr::memmem::find(tag_bytes, needle_sq) {
                Some(pos) => pos + needle_sq.len(),
                None => return fallback_scan,
            }
        }
    };

    // Find the closing quote
    let remaining = &tag_bytes[id_start..];
    let end = remaining
        .iter()
        .position(|&b| b == b'"' || b == b'\'')
        .unwrap_or(remaining.len().min(256));
    let id_value = &remaining[..end];

    // Parse "scan=N" from the id value
    let scan_needle = b"scan=";
    if let Some(pos) = memchr::memmem::find(id_value, scan_needle) {
        let digits_start = pos + scan_needle.len();
        let digits: Vec<u8> = id_value[digits_start..]
            .iter()
            .take_while(|b| b.is_ascii_digit())
            .copied()
            .collect();
        if let Ok(s) = std::str::from_utf8(&digits) {
            if let Ok(n) = s.parse::<u32>() {
                return n;
            }
        }
    }

    fallback_scan
}
```

- [ ] **Step 4: Implement `build_index_by_byte_scan()`**

Add the function in `crates/spectrum-io/src/index.rs` after `build_index_by_scanning`:

```rust
/// Builds a ScanIndex using SIMD-accelerated byte-level scanning.
///
/// This is significantly faster than [`build_index_by_scanning`] on large files
/// because it avoids per-line UTF-8 conversion and String allocation. Uses
/// `memchr::memmem` for SIMD-accelerated pattern search.
///
/// Handles `<spectrum ` tags that span buffer boundaries by keeping an overlap
/// region between buffer fills.
pub fn build_index_by_byte_scan(path: &Path) -> Result<ScanIndex, SpectrumIoError> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    let file = File::open(path).map_err(|e| {
        if e.kind() == std::io::ErrorKind::NotFound {
            SpectrumIoError::FileNotFound {
                path: path.to_path_buf(),
            }
        } else {
            SpectrumIoError::IoError {
                path: path.to_path_buf(),
                source: e,
            }
        }
    })?;

    let mut reader = BufReader::with_capacity(256 * 1024, file);
    let needle = b"<spectrum ";
    let mut offsets = HashMap::new();
    let mut fallback_scan: u32 = 0;
    let mut global_pos: u64 = 0;

    loop {
        let buf = reader.fill_buf().map_err(|e| SpectrumIoError::IoError {
            path: path.to_path_buf(),
            source: e,
        })?;
        if buf.is_empty() {
            break;
        }

        let buf_len = buf.len();
        let mut search_start = 0;

        while let Some(pos) = memchr::memmem::find(&buf[search_start..], needle) {
            let local_pos = search_start + pos;
            let abs_pos = global_pos + local_pos as u64;
            fallback_scan += 1;

            // Extract scan number from the tag starting at this position.
            // Limit how much of the buffer we inspect (a <spectrum> tag is
            // typically <512 bytes).
            let tag_end = (local_pos + 512).min(buf_len);
            let tag_slice = &buf[local_pos..tag_end];
            let scan = extract_scan_from_tag_bytes(tag_slice, fallback_scan);

            if let Some(prev) = offsets.insert(scan, abs_pos) {
                tracing::warn!(
                    "duplicate scan {} in byte-scan: offset {} replaced by {}",
                    scan,
                    prev,
                    abs_pos
                );
            }
            search_start = local_pos + needle.len();
        }

        // Keep overlap to catch tags that span buffer boundaries.
        // We keep `needle.len() - 1` bytes of overlap so that a partial
        // match at the buffer boundary is found in the next fill.
        let overlap = needle.len() - 1;
        let consumed = if buf_len > overlap {
            buf_len - overlap
        } else {
            buf_len // Tiny buffer: consume all to avoid infinite loop
        };
        global_pos += consumed as u64;
        reader.consume(consumed);
    }

    Ok(ScanIndex::new(offsets, IndexSource::BuiltFromScan))
}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `cargo test -p protein-copilot-spectrum-io byte_scan -- --no-capture`
Expected: 2 tests pass — `byte_scan_matches_line_scan` and `byte_scan_indexed_mzml_matches`

- [ ] **Step 6: Run all spectrum-io tests to avoid regressions**

Run: `cargo test -p protein-copilot-spectrum-io`
Expected: all existing tests continue to pass

- [ ] **Step 7: Commit**

```bash
git add crates/spectrum-io/src/index.rs
git commit -m "feat(spectrum-io): add SIMD-accelerated byte-scan index builder

Uses memchr::memmem for block-level byte scanning instead of per-line
String allocation. Expected 5-10x speedup on large mzML files."
```

---

### Task 5: Wire three-layer resolution into `IndexedMzMLReader::open()`

**Files:**
- Modify: `crates/spectrum-io/src/indexed_mzml.rs` (lines 33-48)

- [ ] **Step 1: Write the failing test**

Add to the `mod tests` block in `crates/spectrum-io/src/indexed_mzml.rs`:

```rust
    #[test]
    fn open_creates_disk_cache() {
        let dir = tempfile::tempdir().unwrap();
        let src = fixture_path();
        let copy = dir.path().join("test.mzml");
        std::fs::copy(&src, &copy).unwrap();

        let idx_file = crate::disk_cache::idx_path(&copy);
        assert!(!idx_file.exists(), "idx should not exist before open");

        let _reader = IndexedMzMLReader::open(&copy).unwrap();
        assert!(idx_file.exists(), "idx should exist after open");
    }

    #[test]
    fn open_uses_disk_cache_on_second_call() {
        let dir = tempfile::tempdir().unwrap();
        let src = fixture_path();
        let copy = dir.path().join("test.mzml");
        std::fs::copy(&src, &copy).unwrap();

        // First open: builds index + saves cache
        let reader1 = IndexedMzMLReader::open(&copy).unwrap();

        // Second open: should load from disk cache
        let reader2 = IndexedMzMLReader::open(&copy).unwrap();
        assert_eq!(reader1.index().len(), reader2.index().len());

        // Verify scans match
        for scan in reader1.index().scan_numbers() {
            assert_eq!(
                reader1.index().get_offset(scan),
                reader2.index().get_offset(scan),
            );
        }
    }
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p protein-copilot-spectrum-io open_creates_disk_cache -- --no-capture 2>&1 | head -20`
Expected: FAIL — no .idx file is created

- [ ] **Step 3: Update `IndexedMzMLReader::open()` to use three-layer resolution**

Replace the current `open` method in `crates/spectrum-io/src/indexed_mzml.rs` (lines 38-48) with:

```rust
    pub fn open(path: &Path) -> Result<Self, SpectrumIoError> {
        // Layer 1: Try disk cache
        if let Ok((file_size, file_mtime)) = crate::disk_cache::file_metadata(path) {
            match crate::disk_cache::load_index(path, file_size, file_mtime) {
                Ok(Some(cached_index)) => {
                    tracing::info!(
                        path = %path.display(),
                        scans = cached_index.len(),
                        "loaded index from disk cache"
                    );
                    return Ok(Self {
                        index: cached_index,
                        path: path.to_path_buf(),
                    });
                }
                Ok(None) => {} // Cache miss — continue to layer 2
                Err(e) => {
                    tracing::warn!(error = %e, "disk cache load error, continuing without cache");
                }
            }
        }

        // Layer 2: Try native <indexList>
        let index = if let Some(native) = build_index_from_native_mzml(path)? {
            native
        } else {
            // Layer 3: Byte-scan fallback (accelerated)
            build_index_by_byte_scan(path)?
        };

        // Persist to disk cache for future opens
        if let Ok((file_size, file_mtime)) = crate::disk_cache::file_metadata(path) {
            if let Err(e) = crate::disk_cache::save_index(path, &index, file_size, file_mtime) {
                tracing::warn!(error = %e, "failed to persist index cache (non-fatal)");
            }
        }

        Ok(Self {
            index,
            path: path.to_path_buf(),
        })
    }
```

- [ ] **Step 4: Add missing import**

In `crates/spectrum-io/src/indexed_mzml.rs`, add to the import block (after line 16):

```rust
use crate::index::build_index_by_byte_scan;
```

- [ ] **Step 5: Run the new tests**

Run: `cargo test -p protein-copilot-spectrum-io open_creates_disk_cache open_uses_disk_cache -- --no-capture`
Expected: 2 tests pass

- [ ] **Step 6: Run all spectrum-io tests**

Run: `cargo test -p protein-copilot-spectrum-io`
Expected: all tests pass

- [ ] **Step 7: Commit**

```bash
git add crates/spectrum-io/src/indexed_mzml.rs
git commit -m "feat(spectrum-io): three-layer index resolution in IndexedMzMLReader::open

Priority: disk cache → native indexList → byte-scan fallback.
Persists index to .idx sidecar after building from layers 2 or 3."
```

---

### Task 6: Update `.mcp.json` with timeout configuration

**Files:**
- Modify: `.mcp.json`

- [ ] **Step 1: Add timeout field**

Update `.mcp.json` to:

```json
{
  "mcpServers": {
    "protein-copilot": {
      "command": "cargo",
      "args": ["run", "--release", "-p", "protein-copilot-mcp-server"],
      "timeout": 300,
      "env": {
        "RUST_LOG": "info"
      }
    }
  }
}
```

- [ ] **Step 2: Commit**

```bash
git add .mcp.json
git commit -m "chore: increase MCP server timeout to 300s for large file operations"
```

---

### Task 7: Update README with timeout documentation

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Add performance/configuration section**

Add a new section after `## 当前进度` in `README.md`:

```markdown
## 大文件性能优化

处理大型 mzML 文件（>1GB）时，ProteinCopilot 使用三层索引加速：

1. **磁盘索引缓存**（`.mzml.idx`）— 首次打开后自动生成，后续毫秒级加载
2. **mzML 原生索引**（`<indexList>`）— 直接读取文件末尾，秒级完成
3. **SIMD 字节扫描**（fallback）— 使用 `memchr` 加速全文件扫描

### MCP 超时配置

对于 8GB+ 的大文件，建议在 `.mcp.json` 中增加超时时间：

```json
{
  "mcpServers": {
    "protein-copilot": {
      "timeout": 300
    }
  }
}
```

默认超时 60 秒可能不足以完成首次索引构建。设置 `timeout: 300`（5 分钟）可避免超时错误。
```

- [ ] **Step 2: Commit**

```bash
git add README.md
git commit -m "docs: add large file performance optimization and timeout config to README"
```

---

### Task 8: Full build + test verification

**Files:** None (verification only)

- [ ] **Step 1: Run full workspace build**

Run: `cargo build --workspace`
Expected: compiles with 0 errors

- [ ] **Step 2: Run clippy**

Run: `cargo clippy --workspace -- -D warnings`
Expected: 0 warnings

- [ ] **Step 3: Run all tests**

Run: `cargo test --workspace`
Expected: all tests pass (existing 672+ plus new ~13 tests)

- [ ] **Step 4: Verify .idx file cleanup in test fixtures**

Run: `ls crates/spectrum-io/tests/fixtures/*.idx 2>/dev/null && echo "WARNING: .idx files in fixtures" || echo "OK: no stale .idx files"`
Expected: "OK: no stale .idx files" (tempdir tests auto-clean)

- [ ] **Step 5: Add `.idx` to .gitignore**

If not already present, add `*.idx` to `.gitignore` to prevent committing generated index files:

```
*.idx
```

- [ ] **Step 6: Final commit**

```bash
git add .gitignore
git commit -m "chore: add *.idx to .gitignore (generated disk index caches)"
```
