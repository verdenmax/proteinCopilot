//! Disk-persisted scan index cache in PCIX binary format.
//!
//! Writes a sidecar `.idx` file next to each mzML file so the scan index
//! only needs to be built once. On subsequent opens the index is loaded
//! from disk in O(n) time instead of re-scanning the full mzML.
//!
//! # PCIX Binary Format v2 (little-endian)
//!
//! ```text
//! [magic: 4 bytes "PCIX"]
//! [version: u8 = 2]
//! [source_file_size: u64 LE]
//! [source_file_mtime: u64 LE (Unix epoch secs)]
//! [entry_count: u32 LE]
//! [entries: (scan_number: u32 LE, byte_offset: u64 LE, rt_seconds: f64 LE,
//!            ms_level: u8, has_isolation: u8,
//!            target_mz: f64 LE, lower_offset: f64 LE, upper_offset: f64 LE)
//!           × entry_count]
//! ```
//!
//! Staleness is detected by comparing the source file's size and mtime
//! against the values stored in the cache header.

use std::collections::HashMap;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

use crate::error::SpectrumIoError;
use crate::index::{IndexSource, ScanIndex, ScanMeta};

/// Magic bytes identifying the PCIX format.
const MAGIC: &[u8; 4] = b"PCIX";

/// Current format version.
const VERSION: u8 = 2;

/// Header size: 4 (magic) + 1 (version) + 8 (size) + 8 (mtime) + 4 (count) = 25 bytes.
const HEADER_SIZE: usize = 4 + 1 + 8 + 8 + 4;

/// Size of a single v2 entry:
/// 4 (scan) + 8 (offset) + 8 (rt_seconds) + 1 (ms_level)
/// + 1 (has_isolation) + 8 (target_mz) + 8 (lower) + 8 (upper) = 46 bytes.
const ENTRY_SIZE: usize = 4 + 8 + 8 + 1 + 1 + 8 + 8 + 8;

/// Computes the expected PCIX file size for `entry_count` entries using checked
/// arithmetic. Returns `None` on `usize` overflow (reachable on 32-bit targets
/// where a corrupt/huge `entry_count` would otherwise panic or wrap), letting
/// the loader treat the cache as a miss instead of crashing.
fn expected_total_size(entry_count: usize) -> Option<usize> {
    entry_count
        .checked_mul(ENTRY_SIZE)
        .and_then(|x| x.checked_add(HEADER_SIZE))
}

/// Returns the sidecar `.idx` cache path for a given mzML file.
///
/// Simply appends `.idx` to the original path, e.g.
/// `data/sample.mzML` → `data/sample.mzML.idx`.
pub fn idx_path(mzml_path: &Path) -> PathBuf {
    let mut p = mzml_path.as_os_str().to_owned();
    p.push(".idx");
    PathBuf::from(p)
}

/// Returns `(file_size, mtime_secs)` for the file at `path`.
///
/// `mtime_secs` is the file's last-modified time as seconds since the
/// Unix epoch. Returns a [`SpectrumIoError::DiskCacheError`] on failure.
pub fn file_metadata(path: &Path) -> Result<(u64, u64), SpectrumIoError> {
    let meta = fs::metadata(path).map_err(|e| SpectrumIoError::DiskCacheError {
        path: path.to_path_buf(),
        detail: format!("failed to read metadata: {e}"),
    })?;

    let size = meta.len();
    let mtime = meta
        .modified()
        .map_err(|e| SpectrumIoError::DiskCacheError {
            path: path.to_path_buf(),
            detail: format!("failed to read mtime: {e}"),
        })?
        .duration_since(std::time::UNIX_EPOCH)
        .map_err(|e| SpectrumIoError::DiskCacheError {
            path: path.to_path_buf(),
            detail: format!("mtime before Unix epoch: {e}"),
        })?
        .as_secs();

    Ok((size, mtime))
}

/// Loads a cached [`ScanIndex`] from the sidecar `.idx` file.
///
/// Returns:
/// - `Ok(Some(index))` — cache hit: file exists, is valid, and matches
///   the expected source file size and mtime.
/// - `Ok(None)` — cache miss: file missing, corrupt, truncated, wrong
///   version, or stale (size/mtime mismatch). This is **not** an error.
/// - `Err(...)` — unexpected I/O error during reading.
///
/// Loaded indexes use [`IndexSource::NativeIndex`] since they represent
/// pre-built data.
pub fn load_index(
    mzml_path: &Path,
    expected_size: u64,
    expected_mtime: u64,
) -> Result<Option<ScanIndex>, SpectrumIoError> {
    let cache_path = idx_path(mzml_path);

    // --- Read the entire file ---
    let data = match fs::read(&cache_path) {
        Ok(d) => d,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            tracing::warn!(
                path = %mzml_path.display(),
                "disk cache miss: .idx file not found"
            );
            return Ok(None);
        }
        Err(e) => {
            tracing::warn!(
                path = %mzml_path.display(),
                error = %e,
                "disk cache miss: failed to read .idx file"
            );
            return Ok(None);
        }
    };

    // --- Validate header ---
    if data.len() < HEADER_SIZE {
        tracing::warn!(
            path = %mzml_path.display(),
            bytes = data.len(),
            "disk cache miss: .idx file truncated (too small for header)"
        );
        return Ok(None);
    }

    // Magic
    if &data[0..4] != MAGIC {
        tracing::warn!(
            path = %mzml_path.display(),
            "disk cache miss: bad magic bytes"
        );
        return Ok(None);
    }

    // Version
    let version = data[4];
    if version != VERSION {
        tracing::warn!(
            path = %mzml_path.display(),
            version,
            "disk cache miss: unsupported version"
        );
        return Ok(None);
    }

    // Source file size
    let stored_size = u64::from_le_bytes(data[5..13].try_into().map_err(|_| {
        SpectrumIoError::DiskCacheError {
            path: mzml_path.to_path_buf(),
            detail: "failed to parse source_file_size".to_string(),
        }
    })?);

    // Source file mtime
    let stored_mtime = u64::from_le_bytes(data[13..21].try_into().map_err(|_| {
        SpectrumIoError::DiskCacheError {
            path: mzml_path.to_path_buf(),
            detail: "failed to parse source_file_mtime".to_string(),
        }
    })?);

    // Staleness check
    if stored_size != expected_size {
        tracing::warn!(
            path = %mzml_path.display(),
            stored_size,
            expected_size,
            "disk cache miss: source file size changed"
        );
        return Ok(None);
    }

    if stored_mtime != expected_mtime {
        tracing::warn!(
            path = %mzml_path.display(),
            stored_mtime,
            expected_mtime,
            "disk cache miss: source file mtime changed"
        );
        return Ok(None);
    }

    // Entry count
    let entry_count = u32::from_le_bytes(data[21..25].try_into().map_err(|_| {
        SpectrumIoError::DiskCacheError {
            path: mzml_path.to_path_buf(),
            detail: "failed to parse entry_count".to_string(),
        }
    })?) as usize;

    // Validate total size (checked to avoid usize overflow on 32-bit targets
    // when entry_count is corrupt/huge).
    let expected_total = match expected_total_size(entry_count) {
        Some(t) => t,
        None => {
            tracing::warn!(
                path = %mzml_path.display(),
                entry_count,
                "disk cache miss: declared entry_count overflows usize"
            );
            return Ok(None);
        }
    };
    if data.len() < expected_total {
        tracing::warn!(
            path = %mzml_path.display(),
            expected_total,
            actual = data.len(),
            "disk cache miss: .idx file truncated (not enough entry data)"
        );
        return Ok(None);
    }

    // --- Parse entries ---
    // Cap pre-allocation by the bytes actually present so a corrupt entry_count
    // cannot trigger a giant allocation.
    let capacity = entry_count.min(data.len() / ENTRY_SIZE);
    let mut entries = HashMap::with_capacity(capacity);
    let mut cursor = &data[HEADER_SIZE..];

    for _ in 0..entry_count {
        let mut scan_buf = [0u8; 4];
        let mut offset_buf = [0u8; 8];
        let mut rt_buf = [0u8; 8];
        let mut ms_level_buf = [0u8; 1];
        let mut has_iso_buf = [0u8; 1];
        let mut target_buf = [0u8; 8];
        let mut lower_buf = [0u8; 8];
        let mut upper_buf = [0u8; 8];

        let read =
            |cursor: &mut &[u8], buf: &mut [u8], field: &str| -> Result<(), SpectrumIoError> {
                cursor
                    .read_exact(buf)
                    .map_err(|e| SpectrumIoError::DiskCacheError {
                        path: mzml_path.to_path_buf(),
                        detail: format!("failed to read {field}: {e}"),
                    })
            };

        read(&mut cursor, &mut scan_buf, "scan_number")?;
        read(&mut cursor, &mut offset_buf, "byte_offset")?;
        read(&mut cursor, &mut rt_buf, "rt_seconds")?;
        read(&mut cursor, &mut ms_level_buf, "ms_level")?;
        read(&mut cursor, &mut has_iso_buf, "has_isolation")?;
        read(&mut cursor, &mut target_buf, "target_mz")?;
        read(&mut cursor, &mut lower_buf, "lower_offset")?;
        read(&mut cursor, &mut upper_buf, "upper_offset")?;

        let scan = u32::from_le_bytes(scan_buf);
        let isolation_window = if has_iso_buf[0] != 0 {
            Some((
                f64::from_le_bytes(target_buf),
                f64::from_le_bytes(lower_buf),
                f64::from_le_bytes(upper_buf),
            ))
        } else {
            None
        };

        entries.insert(
            scan,
            ScanMeta {
                offset: u64::from_le_bytes(offset_buf),
                rt_seconds: f64::from_le_bytes(rt_buf),
                ms_level: ms_level_buf[0],
                isolation_window,
            },
        );
    }

    tracing::info!(
        path = %mzml_path.display(),
        scans = entry_count,
        "disk cache loaded"
    );

    Ok(Some(ScanIndex::from_meta(
        entries,
        IndexSource::NativeIndex,
    )))
}

/// Saves a [`ScanIndex`] to the sidecar `.idx` file.
///
/// The cache records the source file's size and mtime so that stale
/// caches can be detected on the next load.
///
/// Failure is considered non-fatal — callers should log a warning and
/// continue without caching.
pub fn save_index(
    mzml_path: &Path,
    index: &ScanIndex,
    file_size: u64,
    file_mtime: u64,
) -> Result<(), SpectrumIoError> {
    let cache_path = idx_path(mzml_path);

    let count = index.len();
    if count > u32::MAX as usize {
        return Err(SpectrumIoError::DiskCacheError {
            path: mzml_path.to_path_buf(),
            detail: format!(
                "index has {} entries, exceeds u32::MAX for PCIX format",
                count
            ),
        });
    }
    let entry_count = count as u32;

    let total_size = HEADER_SIZE + (entry_count as usize) * ENTRY_SIZE;
    let mut buf = Vec::with_capacity(total_size);

    // Header
    buf.extend_from_slice(MAGIC);
    buf.push(VERSION);
    buf.extend_from_slice(&file_size.to_le_bytes());
    buf.extend_from_slice(&file_mtime.to_le_bytes());
    buf.extend_from_slice(&entry_count.to_le_bytes());

    // Entries (sorted by scan number for deterministic output)
    let mut entries: Vec<(u32, &ScanMeta)> = index.iter_meta().map(|(&s, m)| (s, m)).collect();
    entries.sort_by_key(|&(scan, _)| scan);

    for (scan, meta) in &entries {
        buf.extend_from_slice(&scan.to_le_bytes());
        buf.extend_from_slice(&meta.offset.to_le_bytes());
        buf.extend_from_slice(&meta.rt_seconds.to_le_bytes());
        buf.push(meta.ms_level);
        match meta.isolation_window {
            Some((target, lower, upper)) => {
                buf.push(1u8);
                buf.extend_from_slice(&target.to_le_bytes());
                buf.extend_from_slice(&lower.to_le_bytes());
                buf.extend_from_slice(&upper.to_le_bytes());
            }
            None => {
                buf.push(0u8);
                buf.extend_from_slice(&0.0f64.to_le_bytes());
                buf.extend_from_slice(&0.0f64.to_le_bytes());
                buf.extend_from_slice(&0.0f64.to_le_bytes());
            }
        }
    }

    fs::write(&cache_path, &buf).map_err(|e| SpectrumIoError::DiskCacheError {
        path: mzml_path.to_path_buf(),
        detail: format!("failed to write .idx file: {e}"),
    })?;

    tracing::info!(
        path = %mzml_path.display(),
        scans = entry_count,
        cache_path = %cache_path.display(),
        "disk cache saved"
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    /// Helper: create a fake mzML path inside a tempdir.
    fn fake_mzml(dir: &Path) -> PathBuf {
        let p = dir.join("sample.mzML");
        fs::write(&p, b"fake mzML content").ok();
        p
    }

    /// Helper: build a ScanIndex from a slice of (scan, offset) pairs.
    fn build_index(entries: &[(u32, u64)]) -> ScanIndex {
        let offsets: HashMap<u32, u64> = entries.iter().copied().collect();
        ScanIndex::new(offsets, IndexSource::NativeIndex)
    }

    #[test]
    fn round_trip_save_load() {
        let dir = tempfile::tempdir().unwrap();
        let mzml = fake_mzml(dir.path());

        let index = build_index(&[(1, 100), (5, 5000), (10, 99999)]);
        let size = 123456u64;
        let mtime = 1700000000u64;

        save_index(&mzml, &index, size, mtime).unwrap();
        let loaded = load_index(&mzml, size, mtime).unwrap();

        let loaded = loaded.expect("should load cached index");
        assert_eq!(loaded.len(), 3);
        assert_eq!(loaded.get_offset(1), Some(100));
        assert_eq!(loaded.get_offset(5), Some(5000));
        assert_eq!(loaded.get_offset(10), Some(99999));
        assert_eq!(loaded.source(), IndexSource::NativeIndex);
    }

    #[test]
    fn stale_size_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let mzml = fake_mzml(dir.path());

        let index = build_index(&[(1, 100)]);
        let size = 50000u64;
        let mtime = 1700000000u64;

        save_index(&mzml, &index, size, mtime).unwrap();

        // Load with a different expected size
        let result = load_index(&mzml, size + 999, mtime).unwrap();
        assert!(result.is_none(), "stale size should return None");
    }

    #[test]
    fn stale_mtime_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let mzml = fake_mzml(dir.path());

        let index = build_index(&[(1, 100)]);
        let size = 50000u64;
        let mtime = 1700000000u64;

        save_index(&mzml, &index, size, mtime).unwrap();

        // Load with a different expected mtime
        let result = load_index(&mzml, size, mtime + 1).unwrap();
        assert!(result.is_none(), "stale mtime should return None");
    }

    #[test]
    fn missing_idx_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let mzml = dir.path().join("nonexistent.mzML");

        let result = load_index(&mzml, 1000, 1700000000).unwrap();
        assert!(result.is_none(), "missing .idx file should return None");
    }

    #[test]
    fn corrupt_magic_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let mzml = fake_mzml(dir.path());
        let cache = idx_path(&mzml);

        // Write a file with bad magic bytes but enough size for header
        let mut bad_data = vec![0u8; HEADER_SIZE];
        bad_data[0..4].copy_from_slice(b"BAAD");
        fs::write(&cache, &bad_data).unwrap();

        let result = load_index(&mzml, 0, 0).unwrap();
        assert!(result.is_none(), "corrupt magic should return None");
    }

    #[test]
    fn truncated_file_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let mzml = fake_mzml(dir.path());
        let cache = idx_path(&mzml);

        // Write only the magic bytes (4 bytes, less than HEADER_SIZE)
        fs::write(&cache, MAGIC).unwrap();

        let result = load_index(&mzml, 0, 0).unwrap();
        assert!(result.is_none(), "truncated file should return None");
    }

    #[test]
    fn empty_index_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let mzml = fake_mzml(dir.path());

        let index = build_index(&[]);
        let size = 100u64;
        let mtime = 1700000000u64;

        save_index(&mzml, &index, size, mtime).unwrap();
        let loaded = load_index(&mzml, size, mtime).unwrap();

        let loaded = loaded.expect("should load empty cached index");
        assert_eq!(loaded.len(), 0);
        assert!(loaded.is_empty());
    }

    #[test]
    fn large_offsets_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let mzml = fake_mzml(dir.path());

        // Use offsets > 4GB to verify u64 handling
        let index = build_index(&[(1, 5_000_000_000), (2, 8_000_000_000)]);
        let size = 9_000_000_000u64;
        let mtime = 1700000000u64;

        save_index(&mzml, &index, size, mtime).unwrap();
        let loaded = load_index(&mzml, size, mtime).unwrap();

        let loaded = loaded.expect("should load index with large offsets");
        assert_eq!(loaded.len(), 2);
        assert_eq!(loaded.get_offset(1), Some(5_000_000_000));
        assert_eq!(loaded.get_offset(2), Some(8_000_000_000));
    }

    #[test]
    fn idx_path_appends_extension() {
        let path = Path::new("/data/experiment/sample.mzML");
        let result = idx_path(path);
        assert_eq!(result, PathBuf::from("/data/experiment/sample.mzML.idx"));
    }

    /// Helper: build a ScanIndex from ScanMeta entries.
    fn build_meta_index(entries: &[(u32, ScanMeta)]) -> ScanIndex {
        let map: HashMap<u32, ScanMeta> = entries.iter().cloned().collect();
        ScanIndex::from_meta(map, IndexSource::NativeIndex)
    }

    #[test]
    fn v2_round_trip_with_metadata() {
        let dir = tempfile::tempdir().unwrap();
        let mzml = fake_mzml(dir.path());

        let index = build_meta_index(&[
            (
                1,
                ScanMeta {
                    offset: 100,
                    rt_seconds: 120.5,
                    ms_level: 2,
                    isolation_window: Some((500.0, 1.0, 1.0)),
                },
            ),
            (
                5,
                ScanMeta {
                    offset: 5000,
                    rt_seconds: 300.0,
                    ms_level: 1,
                    isolation_window: None,
                },
            ),
            (
                10,
                ScanMeta {
                    offset: 99999,
                    rt_seconds: 600.0,
                    ms_level: 2,
                    isolation_window: Some((750.5, 12.5, 12.5)),
                },
            ),
        ]);
        let size = 123456u64;
        let mtime = 1700000000u64;

        save_index(&mzml, &index, size, mtime).unwrap();
        let loaded = load_index(&mzml, size, mtime).unwrap();

        let loaded = loaded.expect("should load cached index");
        assert_eq!(loaded.len(), 3);
        assert_eq!(loaded.get_offset(1), Some(100));
        assert_eq!(loaded.get_offset(5), Some(5000));

        let meta1 = loaded.get_meta(1).unwrap();
        assert_eq!(meta1.ms_level, 2);
        assert!((meta1.rt_seconds - 120.5).abs() < 0.001);
        let iw = meta1.isolation_window.unwrap();
        assert!((iw.0 - 500.0).abs() < 0.001);

        let meta5 = loaded.get_meta(5).unwrap();
        assert_eq!(meta5.ms_level, 1);
        assert!(meta5.isolation_window.is_none());
    }

    #[test]
    fn v1_cache_triggers_rebuild() {
        let dir = tempfile::tempdir().unwrap();
        let mzml = fake_mzml(dir.path());
        let cache = idx_path(&mzml);

        // Write a valid v1 format file (version=1)
        let mut buf = Vec::new();
        buf.extend_from_slice(MAGIC);
        buf.push(1); // version 1
        buf.extend_from_slice(&100u64.to_le_bytes());
        buf.extend_from_slice(&200u64.to_le_bytes());
        buf.extend_from_slice(&0u32.to_le_bytes());
        std::fs::write(&cache, &buf).unwrap();

        let result = load_index(&mzml, 100, 200).unwrap();
        assert!(result.is_none(), "v1 cache should be rejected by v2 loader");
    }

    #[test]
    fn expected_total_size_guards_against_overflow() {
        // The historical `HEADER_SIZE + entry_count * ENTRY_SIZE` overflows
        // `usize` for a corrupt/huge entry_count (panics in debug, wraps in
        // release on 32-bit targets where usize == u32). The checked helper
        // returns `None` so the loader treats it as a cache miss.
        assert_eq!(expected_total_size(usize::MAX), None);
        // Proof the unchecked multiply would indeed overflow for this input.
        assert!(usize::MAX.checked_mul(ENTRY_SIZE).is_none());
        // Normal counts still compute exactly.
        assert_eq!(expected_total_size(0), Some(HEADER_SIZE));
        assert_eq!(expected_total_size(3), Some(HEADER_SIZE + 3 * ENTRY_SIZE));
    }

    #[test]
    fn huge_entry_count_short_body_is_cache_miss() {
        // A header declaring a u32::MAX entry_count but with a short body must
        // be a clean cache miss — never a panic or a multi-GB allocation.
        let dir = tempfile::tempdir().unwrap();
        let mzml = dir.path().join("sample.mzML");
        let cache = idx_path(&mzml);
        let size = 1234u64;
        let mtime = 1_700_000_000u64;

        let mut buf = Vec::new();
        buf.extend_from_slice(MAGIC);
        buf.push(VERSION);
        buf.extend_from_slice(&size.to_le_bytes());
        buf.extend_from_slice(&mtime.to_le_bytes());
        // Absurd entry_count; no entry bytes follow, so the body is short.
        buf.extend_from_slice(&u32::MAX.to_le_bytes());
        std::fs::write(&cache, &buf).unwrap();

        let result = load_index(&mzml, size, mtime).unwrap();
        assert!(
            result.is_none(),
            "huge entry_count with short body must be a cache miss"
        );
    }
}
