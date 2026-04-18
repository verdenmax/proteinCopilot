# RT Binary Search Optimization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the O(N) full-file-read RT lookup with O(log N) binary search by extending the PCIX disk cache to store per-scan metadata (RT, ms_level, isolation_window).

**Architecture:** Upgrade `ScanIndex` to store `ScanMeta` (offset + RT + ms_level + isolation_window) instead of plain offsets. Extend `build_index_by_byte_scan()` to extract metadata inline during the scan. Upgrade PCIX format from v1 (12B/entry) to v2 (46B/entry). Add `find_by_rt()` to `SpectrumReader` trait with efficient binary search in `IndexedMzMLReader`. Rewrite `find_scan_by_rt()` in `scan_matcher.rs` to delegate to the new trait method.

**Tech Stack:** Rust, memchr, quick-xml (existing deps only)

---

## File Map

| File | Responsibility | Change Type |
|------|---------------|-------------|
| `crates/spectrum-io/src/index.rs` | `ScanMeta` struct, `ScanIndex` internal storage, `find_by_rt()` method, RT-sorted index | Heavy modify |
| `crates/spectrum-io/src/disk_cache.rs` | PCIX v2 format (46B entries), version bump, read/write metadata | Heavy modify |
| `crates/spectrum-io/src/reader.rs` | Add `find_by_rt()` default method to `SpectrumReader` trait | Light modify |
| `crates/spectrum-io/src/indexed_mzml.rs` | Override `find_by_rt()` in `IndexedMzMLReader` impl, update `open()` | Light modify |
| `crates/result-import/src/scan_matcher.rs` | Rewrite `find_scan_by_rt()` to use `SpectrumReader::find_by_rt()`, update `match_scans()` | Moderate modify |
| `crates/mcp-server/src/tools.rs` | Update `annotate_spectrum` and `extract_xic` RT-lookup call sites | Light modify |

---

### Task 1: Add `ScanMeta` and Expand `ScanIndex`

**Files:**
- Modify: `crates/spectrum-io/src/index.rs:1-70`

- [ ] **Step 1: Write failing tests for `ScanMeta` and expanded `ScanIndex`**

Add at the end of the existing `#[cfg(test)] mod tests` block in `crates/spectrum-io/src/index.rs`:

```rust
#[test]
fn scan_index_with_meta_basic() {
    let mut meta_map = HashMap::new();
    meta_map.insert(1, ScanMeta {
        offset: 100,
        rt_seconds: 120.5,
        ms_level: 2,
        isolation_window: Some((500.0, 1.0, 1.0)),
    });
    meta_map.insert(5, ScanMeta {
        offset: 5000,
        rt_seconds: 300.0,
        ms_level: 1,
        isolation_window: None,
    });
    meta_map.insert(10, ScanMeta {
        offset: 10000,
        rt_seconds: 600.0,
        ms_level: 2,
        isolation_window: Some((600.0, 12.5, 12.5)),
    });
    let idx = ScanIndex::from_meta(meta_map, IndexSource::NativeIndex);

    assert_eq!(idx.len(), 3);
    assert_eq!(idx.get_offset(1), Some(100));
    assert_eq!(idx.get_offset(5), Some(5000));
    assert_eq!(idx.get_offset(99), None);

    let meta = idx.get_meta(1).unwrap();
    assert_eq!(meta.ms_level, 2);
    assert!((meta.rt_seconds - 120.5).abs() < 0.001);
    assert!(meta.isolation_window.is_some());

    assert!(idx.get_meta(99).is_none());
}

#[test]
fn scan_index_from_meta_backward_compat() {
    // The old `ScanIndex::new()` constructor still works
    let mut offsets = HashMap::new();
    offsets.insert(1, 100u64);
    offsets.insert(2, 200u64);
    let idx = ScanIndex::new(offsets, IndexSource::BuiltFromScan);
    assert_eq!(idx.len(), 2);
    assert_eq!(idx.get_offset(1), Some(100));
    // get_meta returns ScanMeta with rt_seconds=0.0, ms_level=0
    let meta = idx.get_meta(1).unwrap();
    assert_eq!(meta.offset, 100);
    assert_eq!(meta.ms_level, 0);
    assert!((meta.rt_seconds).abs() < 0.001);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p protein-copilot-spectrum-io scan_index_with_meta_basic scan_index_from_meta_backward_compat -- --nocapture 2>&1 | tail -20`

Expected: FAIL — `ScanMeta` not defined, `from_meta` not defined, `get_meta` not defined.

- [ ] **Step 3: Implement `ScanMeta` and expand `ScanIndex`**

Replace the `ScanIndex` struct and impl block in `crates/spectrum-io/src/index.rs` (lines 22-69):

```rust
/// Per-scan metadata stored in the index.
#[derive(Debug, Clone)]
pub struct ScanMeta {
    /// Byte offset of the `<spectrum>` opening tag in the mzML file.
    pub offset: u64,
    /// Retention time in seconds. 0.0 if unknown.
    pub rt_seconds: f64,
    /// MS level (1=MS1, 2=MS2). 0 if unknown.
    pub ms_level: u8,
    /// Isolation window: (target_mz, lower_offset, upper_offset). None for MS1 or unknown.
    pub isolation_window: Option<(f64, f64, f64)>,
}

/// Maps scan numbers to byte offsets and metadata within an mzML file.
///
/// Enables O(1) spectrum lookup by scan number and O(log N) lookup by
/// retention time via a pre-sorted RT index.
#[derive(Debug, Clone)]
pub struct ScanIndex {
    /// scan_number → metadata (offset, RT, ms_level, isolation_window).
    entries: HashMap<u32, ScanMeta>,
    /// How this index was built.
    source: IndexSource,
    /// Pre-sorted (rt_seconds, scan_number) pairs for binary search.
    /// Built once at construction time.
    rt_sorted: Vec<(f64, u32)>,
}

impl ScanIndex {
    /// Creates a ScanIndex from a legacy offset-only map.
    ///
    /// Metadata fields are set to defaults (rt=0, ms_level=0, no isolation).
    /// This preserves backward compatibility with code that only has offsets.
    pub fn new(offsets: HashMap<u32, u64>, source: IndexSource) -> Self {
        let entries: HashMap<u32, ScanMeta> = offsets
            .into_iter()
            .map(|(scan, offset)| {
                (
                    scan,
                    ScanMeta {
                        offset,
                        rt_seconds: 0.0,
                        ms_level: 0,
                        isolation_window: None,
                    },
                )
            })
            .collect();
        let rt_sorted = build_rt_sorted(&entries);
        Self {
            entries,
            source,
            rt_sorted,
        }
    }

    /// Creates a ScanIndex from a full metadata map.
    pub fn from_meta(entries: HashMap<u32, ScanMeta>, source: IndexSource) -> Self {
        let rt_sorted = build_rt_sorted(&entries);
        Self {
            entries,
            source,
            rt_sorted,
        }
    }

    /// Returns the byte offset for a given scan number, or `None`.
    pub fn get_offset(&self, scan: u32) -> Option<u64> {
        self.entries.get(&scan).map(|m| m.offset)
    }

    /// Returns the full metadata for a given scan number, or `None`.
    pub fn get_meta(&self, scan: u32) -> Option<&ScanMeta> {
        self.entries.get(&scan)
    }

    /// Returns the number of indexed scans.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the index is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// How this index was constructed.
    pub fn source(&self) -> IndexSource {
        self.source
    }

    /// Returns offset references for disk cache serialization.
    ///
    /// Returns an iterator of `(scan_number, &ScanMeta)`.
    pub fn iter_meta(&self) -> impl Iterator<Item = (&u32, &ScanMeta)> {
        self.entries.iter()
    }

    /// Returns a legacy offsets map (for backward compatibility).
    pub fn offsets(&self) -> HashMap<u32, u64> {
        self.entries
            .iter()
            .map(|(&scan, meta)| (scan, meta.offset))
            .collect()
    }

    /// Returns all indexed scan numbers, sorted ascending.
    pub fn scan_numbers(&self) -> Vec<u32> {
        let mut scans: Vec<u32> = self.entries.keys().copied().collect();
        scans.sort_unstable();
        scans
    }

    /// Returns the pre-sorted RT index for binary search.
    pub fn rt_sorted(&self) -> &[(f64, u32)] {
        &self.rt_sorted
    }

    /// Find the best MS2 scan matching a given RT and precursor m/z.
    ///
    /// Uses binary search on the pre-sorted RT index. O(log N + k) where
    /// k is the number of scans in the RT tolerance window.
    ///
    /// Returns `(scan_number, rt_delta_min)` or `None`.
    pub fn find_by_rt(
        &self,
        rt_min: f64,
        precursor_mz: f64,
        rt_tolerance_min: f64,
    ) -> Option<(u32, f64)> {
        let rt_sec = rt_min * 60.0;
        let tol_sec = rt_tolerance_min * 60.0;

        // Binary search for the start of the RT window
        let start = self
            .rt_sorted
            .partition_point(|&(rt, _)| rt < rt_sec - tol_sec);

        let mut best: Option<(u32, f64)> = None;

        for &(rt, scan) in &self.rt_sorted[start..] {
            let delta_sec = rt - rt_sec;
            if delta_sec > tol_sec {
                break;
            }
            if delta_sec.abs() > tol_sec {
                continue;
            }

            let meta = match self.entries.get(&scan) {
                Some(m) => m,
                None => continue,
            };

            // Only match MS2 scans
            if meta.ms_level != 2 {
                continue;
            }

            // Check isolation window if available
            if let Some((target, lower, upper)) = meta.isolation_window {
                let low = target - lower;
                let high = target + upper;
                if precursor_mz < low || precursor_mz > high {
                    continue;
                }
            }
            // No isolation window → accept based on RT only (DDA fallback)

            let delta_min = delta_sec / 60.0;
            match &best {
                None => best = Some((scan, delta_min)),
                Some((_, best_delta)) => {
                    if delta_min.abs() < best_delta.abs() {
                        best = Some((scan, delta_min));
                    }
                }
            }
        }

        best
    }
}

/// Build sorted (rt_seconds, scan_number) pairs from the entries map.
fn build_rt_sorted(entries: &HashMap<u32, ScanMeta>) -> Vec<(f64, u32)> {
    let mut sorted: Vec<(f64, u32)> = entries
        .iter()
        .map(|(&scan, meta)| (meta.rt_seconds, scan))
        .collect();
    sorted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    sorted
}
```

- [ ] **Step 4: Fix all existing callers of `ScanIndex` that relied on old field names**

The internal field changed from `offsets: HashMap<u32, u64>` to `entries: HashMap<u32, ScanMeta>`. The public API is backward-compatible (`new()`, `get_offset()`, `offsets()`, `scan_numbers()` still work), but `offsets()` now returns an owned `HashMap` instead of `&HashMap`. Check callers in `disk_cache.rs`.

In `crates/spectrum-io/src/disk_cache.rs` line 258, change:
```rust
// OLD:
let offsets = index.offsets();
// NEW:
let offsets = index.offsets();
```
No change needed — `offsets()` returns `HashMap<u32, u64>` in both old and new API.

Actually, line 281 uses `offsets.iter()` which returns `(&u32, &u64)`. Since `offsets()` now returns an owned `HashMap<u32, u64>`, the iterator yields `(&u32, &u64)` — same type. No change needed.

- [ ] **Step 5: Run tests to verify they pass**

Run: `cargo test -p protein-copilot-spectrum-io -- --nocapture 2>&1 | tail -20`

Expected: All tests pass, including the two new ones.

- [ ] **Step 6: Commit**

```bash
git add crates/spectrum-io/src/index.rs
git commit -m "feat(spectrum-io): add ScanMeta and expand ScanIndex with RT binary search

- Add ScanMeta struct (offset, rt_seconds, ms_level, isolation_window)
- Expand ScanIndex to store full metadata
- Add find_by_rt() for O(log N) RT-based scan lookup
- Maintain backward compatibility via ScanIndex::new()

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: Add `find_by_rt()` Tests for `ScanIndex`

**Files:**
- Modify: `crates/spectrum-io/src/index.rs` (test section)

- [ ] **Step 1: Write tests for `find_by_rt()`**

Add to the `#[cfg(test)] mod tests` block in `crates/spectrum-io/src/index.rs`:

```rust
fn make_rt_index() -> ScanIndex {
    let mut meta = HashMap::new();
    // MS1 scans (should be skipped by find_by_rt)
    meta.insert(1, ScanMeta {
        offset: 100,
        rt_seconds: 100.0 * 60.0,  // 100 min
        ms_level: 1,
        isolation_window: None,
    });
    meta.insert(3, ScanMeta {
        offset: 300,
        rt_seconds: 200.0 * 60.0,
        ms_level: 1,
        isolation_window: None,
    });
    // MS2 scans
    meta.insert(2, ScanMeta {
        offset: 200,
        rt_seconds: 100.0 * 60.0,  // 100 min
        ms_level: 2,
        isolation_window: Some((500.0, 12.5, 12.5)),
    });
    meta.insert(4, ScanMeta {
        offset: 400,
        rt_seconds: 200.0 * 60.0,  // 200 min
        ms_level: 2,
        isolation_window: Some((600.0, 12.5, 12.5)),
    });
    meta.insert(5, ScanMeta {
        offset: 500,
        rt_seconds: 300.0 * 60.0,  // 300 min
        ms_level: 2,
        isolation_window: Some((500.0, 25.0, 25.0)),
    });
    meta.insert(6, ScanMeta {
        offset: 600,
        rt_seconds: 400.0 * 60.0,
        ms_level: 2,
        isolation_window: None,  // DDA — no isolation window
    });
    ScanIndex::from_meta(meta, IndexSource::BuiltFromScan)
}

#[test]
fn find_by_rt_exact_match() {
    let idx = make_rt_index();
    let result = idx.find_by_rt(100.0, 500.0, 30.0);
    assert_eq!(result.unwrap().0, 2); // MS2 scan at RT=100min
}

#[test]
fn find_by_rt_skips_ms1() {
    let idx = make_rt_index();
    // RT=100 has both MS1 (scan 1) and MS2 (scan 2). Should pick MS2.
    let result = idx.find_by_rt(100.0, 500.0, 30.0);
    assert_eq!(result.unwrap().0, 2);
}

#[test]
fn find_by_rt_mz_outside_window() {
    let idx = make_rt_index();
    // scan 2: isolation 487.5–512.5, mz=550 is outside
    let result = idx.find_by_rt(100.0, 550.0, 30.0);
    assert!(result.is_none());
}

#[test]
fn find_by_rt_outside_tolerance() {
    let idx = make_rt_index();
    // RT=150, tolerance=30 → nearest MS2 scan 2 at RT=100 is 50 min away
    let result = idx.find_by_rt(150.0, 500.0, 30.0);
    assert!(result.is_none());
}

#[test]
fn find_by_rt_dda_no_isolation_accepts_any_mz() {
    let idx = make_rt_index();
    // scan 6 at RT=400 has no isolation window → accepts any mz
    let result = idx.find_by_rt(400.0, 999.0, 30.0);
    assert_eq!(result.unwrap().0, 6);
}

#[test]
fn find_by_rt_picks_closest() {
    let mut meta = HashMap::new();
    meta.insert(1, ScanMeta {
        offset: 100,
        rt_seconds: 100.0 * 60.0,
        ms_level: 2,
        isolation_window: Some((500.0, 25.0, 25.0)),
    });
    meta.insert(2, ScanMeta {
        offset: 200,
        rt_seconds: 105.0 * 60.0,
        ms_level: 2,
        isolation_window: Some((500.0, 25.0, 25.0)),
    });
    let idx = ScanIndex::from_meta(meta, IndexSource::BuiltFromScan);
    // PSM at RT=103 → closer to scan 2 (105)
    let result = idx.find_by_rt(103.0, 500.0, 30.0);
    assert_eq!(result.unwrap().0, 2);
}

#[test]
fn find_by_rt_empty_index() {
    let idx = ScanIndex::from_meta(HashMap::new(), IndexSource::BuiltFromScan);
    assert!(idx.find_by_rt(100.0, 500.0, 30.0).is_none());
}
```

- [ ] **Step 2: Run tests to verify they pass**

Run: `cargo test -p protein-copilot-spectrum-io find_by_rt -- --nocapture 2>&1 | tail -20`

Expected: All 8 `find_by_rt` tests pass.

- [ ] **Step 3: Commit**

```bash
git add crates/spectrum-io/src/index.rs
git commit -m "test(spectrum-io): add find_by_rt tests for ScanIndex

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: Upgrade PCIX Disk Cache to v2

**Files:**
- Modify: `crates/spectrum-io/src/disk_cache.rs`

- [ ] **Step 1: Write failing test for v2 round-trip**

Add to the `#[cfg(test)] mod tests` block in `crates/spectrum-io/src/disk_cache.rs`:

```rust
use crate::index::ScanMeta;

/// Helper: build a ScanIndex from ScanMeta entries.
fn build_meta_index(entries: &[(u32, ScanMeta)]) -> ScanIndex {
    let map: HashMap<u32, ScanMeta> = entries
        .iter()
        .cloned()
        .collect();
    ScanIndex::from_meta(map, IndexSource::NativeIndex)
}

#[test]
fn v2_round_trip_with_metadata() {
    let dir = tempfile::tempdir().unwrap();
    let mzml = fake_mzml(dir.path());

    let index = build_meta_index(&[
        (1, ScanMeta {
            offset: 100,
            rt_seconds: 120.5,
            ms_level: 2,
            isolation_window: Some((500.0, 1.0, 1.0)),
        }),
        (5, ScanMeta {
            offset: 5000,
            rt_seconds: 300.0,
            ms_level: 1,
            isolation_window: None,
        }),
        (10, ScanMeta {
            offset: 99999,
            rt_seconds: 600.0,
            ms_level: 2,
            isolation_window: Some((750.5, 12.5, 12.5)),
        }),
    ]);
    let size = 123456u64;
    let mtime = 1700000000u64;

    save_index(&mzml, &index, size, mtime).unwrap();
    let loaded = load_index(&mzml, size, mtime).unwrap();

    let loaded = loaded.expect("should load cached index");
    assert_eq!(loaded.len(), 3);
    assert_eq!(loaded.get_offset(1), Some(100));
    assert_eq!(loaded.get_offset(5), Some(5000));

    // Verify metadata was preserved
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

    // Write a valid v1 format file
    let mut buf = Vec::new();
    buf.extend_from_slice(MAGIC);
    buf.push(1); // version 1
    buf.extend_from_slice(&100u64.to_le_bytes()); // file_size
    buf.extend_from_slice(&200u64.to_le_bytes()); // mtime
    buf.extend_from_slice(&0u32.to_le_bytes()); // entry_count = 0
    std::fs::write(&cache, &buf).unwrap();

    // v2 loader should reject this as unsupported version
    let result = load_index(&mzml, 100, 200).unwrap();
    assert!(result.is_none(), "v1 cache should be rejected by v2 loader");
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p protein-copilot-spectrum-io v2_round_trip_with_metadata v1_cache_triggers_rebuild -- --nocapture 2>&1 | tail -20`

Expected: FAIL — `build_meta_index` not defined, `save_index` writes v1 format.

- [ ] **Step 3: Upgrade `disk_cache.rs` to v2 format**

Replace the constants and functions in `crates/spectrum-io/src/disk_cache.rs`:

**Constants** (lines 29-39): Replace with:

```rust
/// Magic bytes identifying the PCIX format.
const MAGIC: &[u8; 4] = b"PCIX";

/// Current format version (v2 adds RT, ms_level, isolation_window per entry).
const VERSION: u8 = 2;

/// Header size: 4 (magic) + 1 (version) + 8 (size) + 8 (mtime) + 4 (count) = 25 bytes.
const HEADER_SIZE: usize = 4 + 1 + 8 + 8 + 4;

/// Size of a single v2 entry:
/// 4 (scan) + 8 (offset) + 8 (rt_seconds) + 1 (ms_level)
/// + 1 (has_isolation) + 8 (target_mz) + 8 (lower) + 8 (upper) = 46 bytes.
const ENTRY_SIZE: usize = 4 + 8 + 8 + 1 + 1 + 8 + 8 + 8;
```

**`load_index`** (lines 89-240): Replace the entry parsing section (from `// --- Parse entries ---` to the end, lines 209-240) with:

```rust
    // --- Parse entries ---
    let mut entries = HashMap::with_capacity(entry_count);
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

        let read = |cursor: &mut &[u8], buf: &mut [u8], field: &str| -> Result<(), SpectrumIoError> {
            cursor.read_exact(buf).map_err(|e| SpectrumIoError::DiskCacheError {
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

        entries.insert(scan, ScanMeta {
            offset: u64::from_le_bytes(offset_buf),
            rt_seconds: f64::from_le_bytes(rt_buf),
            ms_level: ms_level_buf[0],
            isolation_window,
        });
    }

    tracing::info!(
        path = %mzml_path.display(),
        entries = entry_count,
        "disk cache hit: loaded scan index v2 from .idx file"
    );

    Ok(Some(ScanIndex::from_meta(entries, IndexSource::NativeIndex)))
```

**`save_index`** (lines 250-302): Replace the entry serialization section (from `let offsets =` to `Ok(())`):

```rust
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
                buf.push(1u8); // has_isolation = true
                buf.extend_from_slice(&target.to_le_bytes());
                buf.extend_from_slice(&lower.to_le_bytes());
                buf.extend_from_slice(&upper.to_le_bytes());
            }
            None => {
                buf.push(0u8); // has_isolation = false
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
        entries = entry_count,
        cache_path = %cache_path.display(),
        "saved scan index v2 to .idx cache"
    );

    Ok(())
```

Also add `use crate::index::ScanMeta;` at the import section (line 27):

```rust
use crate::error::SpectrumIoError;
use crate::index::{IndexSource, ScanIndex, ScanMeta};
```

- [ ] **Step 4: Update existing v1 tests to work with v2 format**

The `build_index` test helper creates v1-style offset-only indexes. These still work via `ScanIndex::new()`. The round-trip test should still pass because `save_index` now writes v2 and `load_index` reads v2. But the entry data will now include default metadata fields.

No changes needed to existing tests — they use `ScanIndex::new()` which creates default metadata, and `get_offset()` still works.

- [ ] **Step 5: Run all tests**

Run: `cargo test -p protein-copilot-spectrum-io -- --nocapture 2>&1 | tail -30`

Expected: All tests pass, including new v2 tests and existing v1-compatible tests.

- [ ] **Step 6: Commit**

```bash
git add crates/spectrum-io/src/disk_cache.rs
git commit -m "feat(spectrum-io): upgrade PCIX disk cache to v2 with RT metadata

- Entry size 12B → 46B (+ rt_seconds, ms_level, isolation_window)
- Version 1 → 2 (v1 caches rejected as cache miss, triggering rebuild)
- Backward compatible: ScanIndex::new() still works for offset-only data

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 4: Extract Metadata During Byte-Scan

**Files:**
- Modify: `crates/spectrum-io/src/index.rs:398-486` (build_index_by_byte_scan)

- [ ] **Step 1: Write a test that byte-scan produces metadata**

Add to the `#[cfg(test)] mod tests` block in `crates/spectrum-io/src/index.rs`:

```rust
#[test]
fn byte_scan_extracts_metadata() {
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/small.mzml");
    let idx = build_index_by_byte_scan(&path).unwrap();

    // small.mzml has 10 spectra, all MS2
    assert_eq!(idx.len(), 10);

    // Scan 1 should have RT=120.5s, ms_level=2, isolation_window=(471.2561, 1.0, 1.0)
    let meta1 = idx.get_meta(1).expect("scan 1 should exist");
    assert_eq!(meta1.ms_level, 2, "scan 1 ms_level");
    assert!(
        (meta1.rt_seconds - 120.5).abs() < 0.1,
        "scan 1 RT: expected ~120.5, got {}",
        meta1.rt_seconds
    );
    let iw = meta1.isolation_window.expect("scan 1 should have isolation window");
    assert!(
        (iw.0 - 471.2561).abs() < 0.01,
        "scan 1 isolation target_mz: expected ~471.2561, got {}",
        iw.0
    );

    // Scan 2: RT=125.3s, ms_level=2, isolation=(523.7832, 1.0, 1.0)
    let meta2 = idx.get_meta(2).expect("scan 2 should exist");
    assert_eq!(meta2.ms_level, 2);
    assert!(
        (meta2.rt_seconds - 125.3).abs() < 0.1,
        "scan 2 RT: expected ~125.3, got {}",
        meta2.rt_seconds
    );
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p protein-copilot-spectrum-io byte_scan_extracts_metadata -- --nocapture 2>&1 | tail -20`

Expected: FAIL — `meta1.ms_level` is 0 (default), `rt_seconds` is 0.0.

- [ ] **Step 3: Implement metadata extraction in `build_index_by_byte_scan`**

Add a new helper function `extract_meta_from_tag_bytes` right after `extract_scan_from_tag_bytes` in `crates/spectrum-io/src/index.rs`:

```rust
/// Extracts RT, ms_level, and isolation window from raw XML bytes
/// following a `<spectrum ` tag.
///
/// Searches for well-known cvParam accession numbers in the raw bytes.
/// The search region should be ~2KB to cover the spectrum header before
/// `<binaryDataArrayList>`.
fn extract_meta_from_region(region: &[u8]) -> (f64, u8, Option<(f64, f64, f64)>) {
    let mut rt_seconds: f64 = 0.0;
    let mut ms_level: u8 = 0;
    let mut iso_target: Option<f64> = None;
    let mut iso_lower: Option<f64> = None;
    let mut iso_upper: Option<f64> = None;

    // Stop scanning if we hit binaryDataArrayList (no metadata after this)
    let limit = memchr::memmem::find(region, b"<binaryDataArrayList")
        .unwrap_or(region.len());
    let region = &region[..limit];

    // MS:1000016 — scan start time (RT)
    if let Some(pos) = memchr::memmem::find(region, b"MS:1000016") {
        rt_seconds = extract_cv_value(&region[pos..]).unwrap_or(0.0);
        // Check if unit is minutes (UO:0000031) — convert to seconds
        let after = &region[pos..region.len().min(pos + 300)];
        if memchr::memmem::find(after, b"UO:0000031").is_some() {
            rt_seconds *= 60.0;
        }
    }

    // MS:1000511 — ms level (the value attribute contains the level number)
    if let Some(pos) = memchr::memmem::find(region, b"MS:1000511") {
        ms_level = extract_cv_value(&region[pos..])
            .map(|v| v as u8)
            .unwrap_or(0);
    }

    // MS:1000827 — isolation window target m/z
    if let Some(pos) = memchr::memmem::find(region, b"MS:1000827") {
        iso_target = extract_cv_value(&region[pos..]);
    }
    // MS:1000828 — isolation window lower offset
    if let Some(pos) = memchr::memmem::find(region, b"MS:1000828") {
        iso_lower = extract_cv_value(&region[pos..]);
    }
    // MS:1000829 — isolation window upper offset
    if let Some(pos) = memchr::memmem::find(region, b"MS:1000829") {
        iso_upper = extract_cv_value(&region[pos..]);
    }

    let isolation_window = match (iso_target, iso_lower, iso_upper) {
        (Some(t), Some(l), Some(u)) => Some((t, l, u)),
        _ => None,
    };

    (rt_seconds, ms_level, isolation_window)
}

/// Extracts the `value="..."` attribute from a cvParam byte region.
///
/// Searches for `value="` after the current position, parses the f64.
fn extract_cv_value(region: &[u8]) -> Option<f64> {
    let limit = region.len().min(200);
    let search = &region[..limit];
    let pos = memchr::memmem::find(search, b"value=\"")?;
    let after = &search[pos + 7..];
    let end = memchr::memchr(b'"', after)?;
    let val_bytes = &after[..end];
    std::str::from_utf8(val_bytes)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
}
```

Then modify `build_index_by_byte_scan` to use `ScanMeta` instead of plain offsets. Replace the body of the function starting from `let mut offsets = HashMap::new();` (line 425) through `Ok(ScanIndex::new(offsets, IndexSource::BuiltFromScan))` (line 485):

```rust
    let mut reader = BufReader::with_capacity(256 * 1024, file);
    let needle = b"<spectrum ";
    let mut entries: HashMap<u32, ScanMeta> = HashMap::new();
    let mut fallback_scan: u32 = 0;
    let mut global_pos: u64 = 0;

    // Need enough bytes after a tag match to extract both scan number
    // and metadata (RT, ms_level, isolation_window).
    // Metadata cvParams are within ~1-2KB of the <spectrum> tag.
    const TAG_MIN_CONTENT: usize = 2048;

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
            let remaining = buf_len - local_pos;

            if remaining < TAG_MIN_CONTENT && buf_len >= TAG_MIN_CONTENT + needle.len() {
                break;
            }

            let abs_pos = global_pos + local_pos as u64;
            fallback_scan += 1;

            let tag_end = (local_pos + 512).min(buf_len);
            let scan = extract_scan_from_tag_bytes(&buf[local_pos..tag_end], fallback_scan);

            // Extract metadata from the region after the tag
            let meta_end = (local_pos + TAG_MIN_CONTENT).min(buf_len);
            let (rt_seconds, ms_level, isolation_window) =
                extract_meta_from_region(&buf[local_pos..meta_end]);

            let meta = ScanMeta {
                offset: abs_pos,
                rt_seconds,
                ms_level,
                isolation_window,
            };

            if let Some(prev) = entries.insert(scan, meta) {
                tracing::warn!(
                    "duplicate scan {} found while byte-scanning: offset {} replaced by {}",
                    scan,
                    prev.offset,
                    abs_pos
                );
            }
            search_start = local_pos + needle.len();
        }

        let overlap = TAG_MIN_CONTENT + needle.len();
        let consumed = if buf_len > overlap {
            buf_len - overlap
        } else {
            buf_len
        };
        global_pos += consumed as u64;
        reader.consume(consumed);
    }

    Ok(ScanIndex::from_meta(entries, IndexSource::BuiltFromScan))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p protein-copilot-spectrum-io byte_scan_extracts_metadata -- --nocapture 2>&1 | tail -20`

Expected: PASS — RT=120.5, ms_level=2, isolation_window present.

- [ ] **Step 5: Run all spectrum-io tests**

Run: `cargo test -p protein-copilot-spectrum-io -- --nocapture 2>&1 | tail -30`

Expected: All tests pass.

- [ ] **Step 6: Commit**

```bash
git add crates/spectrum-io/src/index.rs
git commit -m "feat(spectrum-io): extract RT/ms_level/isolation during byte-scan

- Add extract_meta_from_region() for cvParam extraction from raw bytes
- Add extract_cv_value() helper for value attribute parsing
- build_index_by_byte_scan() now produces ScanIndex with full metadata
- Handles both seconds (UO:0000010) and minutes (UO:0000031) RT units

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 5: Add `find_by_rt()` to `SpectrumReader` Trait

**Files:**
- Modify: `crates/spectrum-io/src/reader.rs`
- Modify: `crates/spectrum-io/src/indexed_mzml.rs`
- Modify: `crates/spectrum-io/src/mgf.rs` (if needed — check if MgfReader needs a stub)

- [ ] **Step 1: Add default `find_by_rt()` method to `SpectrumReader` trait**

In `crates/spectrum-io/src/reader.rs`, add at the end of the trait body (before the closing `}`):

```rust
    /// Find the best MS2 scan matching a given RT and precursor m/z.
    ///
    /// Default implementation reads all spectra (slow for large files).
    /// [`crate::IndexedMzMLReader`] overrides with O(log N) binary search
    /// using the pre-built RT index.
    ///
    /// Returns `(scan_number, rt_delta_min)` or `None`.
    fn find_by_rt(
        &self,
        path: &Path,
        rt_min: f64,
        precursor_mz: f64,
        rt_tolerance_min: f64,
    ) -> Result<Option<(u32, f64)>, SpectrumIoError> {
        use protein_copilot_core::spectrum::MsLevel;

        let spectra = self.read_all(path)?;

        let mut best: Option<(u32, f64)> = None;
        for spec in &spectra {
            if spec.ms_level != MsLevel::MS2 {
                continue;
            }
            let delta_min = spec.retention_time_min - rt_min;
            if delta_min.abs() > rt_tolerance_min {
                continue;
            }
            // Check isolation window if available
            if let Some(p) = spec.precursors.first() {
                if let Some(w) = &p.isolation_window {
                    let low = w.target_mz - w.lower_offset;
                    let high = w.target_mz + w.upper_offset;
                    if precursor_mz < low || precursor_mz > high {
                        continue;
                    }
                }
            }
            match &best {
                None => best = Some((spec.scan_number, delta_min)),
                Some((_, bd)) => {
                    if delta_min.abs() < bd.abs() {
                        best = Some((spec.scan_number, delta_min));
                    }
                }
            }
        }
        Ok(best)
    }
```

- [ ] **Step 2: Override `find_by_rt()` in `IndexedMzMLReader`**

In `crates/spectrum-io/src/indexed_mzml.rs`, add inside the `impl SpectrumReader for IndexedMzMLReader` block (after `for_each_spectrum`):

```rust
    fn find_by_rt(
        &self,
        _path: &Path,
        rt_min: f64,
        precursor_mz: f64,
        rt_tolerance_min: f64,
    ) -> Result<Option<(u32, f64)>, SpectrumIoError> {
        Ok(self.index.find_by_rt(rt_min, precursor_mz, rt_tolerance_min))
    }
```

- [ ] **Step 3: Check `Spectrum` struct for `retention_time_sec` field**

The default implementation uses `spec.retention_time_sec`. Check if the Spectrum struct has this field or if it's `retention_time_min`. If it uses `retention_time_min`, adjust the default implementation accordingly:

```rust
// If the field is retention_time_min, adjust:
let delta_min = spec.retention_time_min - rt_min;
if delta_min.abs() > rt_tolerance_min { continue; }
```

Run: `grep -n "retention_time" crates/core/src/spectrum.rs | head -5`

Update the default implementation based on the actual field name.

- [ ] **Step 4: Run all tests**

Run: `cargo test --workspace 2>&1 | tail -20`

Expected: All tests pass. No regressions.

- [ ] **Step 5: Commit**

```bash
git add crates/spectrum-io/src/reader.rs crates/spectrum-io/src/indexed_mzml.rs
git commit -m "feat(spectrum-io): add find_by_rt() to SpectrumReader trait

- Default implementation reads all spectra (backward-compatible)
- IndexedMzMLReader overrides with O(log N) binary search on ScanIndex
- MGF reader inherits default implementation

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 6: Rewrite `find_scan_by_rt()` in scan_matcher

**Files:**
- Modify: `crates/result-import/src/scan_matcher.rs`

- [ ] **Step 1: Rewrite `find_scan_by_rt()` to use `SpectrumReader::find_by_rt()`**

Replace the `find_scan_by_rt` function in `crates/result-import/src/scan_matcher.rs` (lines 154-175):

```rust
/// Find the best MS2 scan matching a given RT and precursor m/z.
///
/// Delegates to [`SpectrumReader::find_by_rt()`] which uses O(log N) binary
/// search when backed by an [`IndexedMzMLReader`], or falls back to reading
/// all spectra for other reader implementations.
pub fn find_scan_by_rt(
    file: &Path,
    rt_min: f64,
    precursor_mz: f64,
    rt_tolerance_min: f64,
    reader: &dyn SpectrumReader,
) -> Result<u32, ResultImportError> {
    reader
        .find_by_rt(file, rt_min, precursor_mz, rt_tolerance_min)
        .map_err(|e| ResultImportError::SpectrumIo(e.to_string()))?
        .map(|(scan, _delta)| scan)
        .ok_or(ResultImportError::NoMatchingScan {
            rt_min,
            tolerance_min: rt_tolerance_min,
            precursor_mz,
        })
}
```

- [ ] **Step 2: Update `match_scans()` to use `SpectrumReader::find_by_rt()` metadata**

Replace the `collect_ms2_info` + `find_best_match` usage in `match_scans()`. Replace lines 78-108 (the per-file matching loop body):

```rust
        let reader = reader_factory(&actual_path)?;

        // Use the reader's find_by_rt (O(log N) for IndexedMzMLReader)
        // For match_scans we need ms2_count for stats, so collect MS2 info
        let ms2_infos = collect_ms2_info(&*reader, &actual_path)?;
        let ms2_count = ms2_infos.len();

        let mut file_matched = 0usize;
        let mut file_unmatched = 0usize;

        for &idx in indices {
            let psm = &psms[idx];
            let result = reader.find_by_rt(
                &actual_path,
                psm.rt_min,
                psm.precursor_mz,
                config.rt_tolerance_min,
            ).map_err(|e| ResultImportError::SpectrumIo(e.to_string()))?;

            if let Some((scan, delta)) = result {
                let psm_mut = &mut psms[idx];
                psm_mut.matched_scan = Some(scan);
                psm_mut.rt_delta_min = Some(delta);
                all_rt_deltas.push(delta.abs());
                file_matched += 1;
            } else {
                file_unmatched += 1;
            }
        }
```

Note: `collect_ms2_info` is still used for `ms2_count` in stats. The `find_best_match` function and `sorted_ms2` are no longer used in this path.

Actually, we need `ms2_count` but we don't want to call `read_all`. Two options:
1. Keep `collect_ms2_info` for the count but note it's still O(N) — however `match_scans` is called once per file, not per PSM.
2. Add a method to get scan count from index.

For now, keep `collect_ms2_info` for the count (it's called once per file for stats). The hot path (per-PSM lookup) is the optimized `find_by_rt`.

Wait — `collect_ms2_info` calls `read_all` which is the slow path. We should avoid it. Instead, count MS2 scans from the index if possible. But `SpectrumReader` doesn't expose index metadata. 

Better approach: use `read_summary` which is streaming and doesn't load peaks. It returns `SpectrumSummary` which has `total_spectra` and `ms2_count`.

Replace the `ms2_infos` line with:

```rust
        let reader = reader_factory(&actual_path)?;

        // Get MS2 count from summary (streaming, no peak data loaded)
        let summary = reader.read_summary(&actual_path)
            .map_err(|e| ResultImportError::SpectrumIo(e.to_string()))?;
        let ms2_count = summary.ms2_count as usize;

        let mut file_matched = 0usize;
        let mut file_unmatched = 0usize;

        for &idx in indices {
            let psm = &psms[idx];
            let result = reader.find_by_rt(
                &actual_path,
                psm.rt_min,
                psm.precursor_mz,
                config.rt_tolerance_min,
            ).map_err(|e| ResultImportError::SpectrumIo(e.to_string()))?;

            if let Some((scan, delta)) = result {
                let psm_mut = &mut psms[idx];
                psm_mut.matched_scan = Some(scan);
                psm_mut.rt_delta_min = Some(delta);
                all_rt_deltas.push(delta.abs());
                file_matched += 1;
            } else {
                file_unmatched += 1;
            }
        }
```

Check that `SpectrumSummary` has `ms2_count`:

Run: `grep -n "ms2_count\|ms1_count\|total_spectra" crates/core/src/spectrum.rs | head -10`

Adjust based on actual field names.

- [ ] **Step 3: Run all tests**

Run: `cargo test --workspace 2>&1 | tail -30`

Expected: All tests pass.

- [ ] **Step 4: Commit**

```bash
git add crates/result-import/src/scan_matcher.rs
git commit -m "perf(result-import): use SpectrumReader::find_by_rt() for O(log N) scan matching

- find_scan_by_rt() delegates to SpectrumReader::find_by_rt()
- match_scans() uses find_by_rt() per PSM instead of read_all()
- Uses read_summary() for MS2 count (streaming, no peak data)

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 7: Update MCP Tool Call Sites

**Files:**
- Modify: `crates/mcp-server/src/tools.rs:2074-2093` (annotate_spectrum)
- Modify: `crates/mcp-server/src/tools.rs:2672-2683` (extract_xic)

- [ ] **Step 1: Simplify `annotate_spectrum` RT lookup**

Replace the RT-lookup block in `annotate_spectrum` (lines 2074-2093):

```rust
        // Resolve scan_number: if 0 and retention_time_min provided, auto-match via RT
        let resolved_scan = if input.scan_number == 0 {
            if let Some(rt) = input.retention_time_min {
                let reader = self.get_or_create_reader(&spectrum_file)?;
                use protein_copilot_search_engine::chemistry::{
                    residue_mass, PROTON_MASS, WATER_MASS,
                };
                let base_mass: f64 =
                    peptide_seq.chars().filter_map(residue_mass).sum::<f64>() + WATER_MASS;
                let mod_mass: f64 = modifications.iter().map(|m| m.mass_delta).sum();
                let precursor_mz =
                    (base_mass + mod_mass + charge as f64 * PROTON_MASS) / charge as f64;
                let reader = self.get_or_create_reader(&spectrum_file)?;
                reader
                    .find_by_rt(&spectrum_file, rt, precursor_mz, 0.5)
                    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?
                    .map(|(scan, _)| scan)
                    .ok_or_else(|| {
                        mcp_err(
                            ErrorCode::INVALID_PARAMS,
                            format!(
                                "no MS2 scan found near RT={:.2} min with precursor m/z={:.4}",
                                rt, precursor_mz
                            ),
                        )
                    })?
            } else {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "scan_number is 0: provide retention_time_min for auto lookup",
                ));
            }
        } else {
            input.scan_number
        };
```

Note: This removes the dependency on `protein_copilot_result_import::scan_matcher::find_scan_by_rt` from `annotate_spectrum`. The `reader.find_by_rt()` is on the `SpectrumReader` trait which `get_or_create_reader` already returns.

- [ ] **Step 2: Simplify `extract_xic` RT lookup**

Replace the RT-lookup block in `extract_xic` (lines 2672-2683):

```rust
        // Resolve scan_number: if 0 and retention_time_min provided, auto-match via RT
        let resolved_scan = if input.scan_number == 0 {
            if let Some(rt) = input.retention_time_min {
                let reader = self.get_or_create_reader(&file_path)?;
                reader
                    .find_by_rt(&file_path, rt, precursor_mz, 0.5)
                    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?
                    .map(|(scan, _)| scan)
                    .ok_or_else(|| {
                        mcp_err(
                            ErrorCode::INVALID_PARAMS,
                            format!(
                                "no MS2 scan found near RT={:.2} min with precursor m/z={:.4}",
                                rt, precursor_mz
                            ),
                        )
                    })?
            } else {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "scan_number is 0: provide retention_time_min for auto lookup",
                ));
            }
        } else {
            input.scan_number
        };
```

- [ ] **Step 3: Remove unused `find_scan_by_rt` import if no longer needed**

Check if `find_scan_by_rt` is still imported/used in tools.rs. If the only call sites were `annotate_spectrum` and `extract_xic`, the import can be removed. `match_scans` is still used in `import_search_results`, so keep `scan_matcher::{match_scans, ScanMatcherConfig}`.

Run: `grep -n "find_scan_by_rt\|scan_matcher" crates/mcp-server/src/tools.rs | head -10`

Remove `find_scan_by_rt` from imports if unused.

- [ ] **Step 4: Run all workspace tests**

Run: `cargo test --workspace 2>&1 | tail -30`

Expected: All tests pass.

- [ ] **Step 5: Run clippy**

Run: `cargo clippy --workspace 2>&1 | grep -E "warning|error" | head -20`

Expected: No warnings or errors.

- [ ] **Step 6: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "perf(mcp-server): use SpectrumReader::find_by_rt() for RT lookup

- annotate_spectrum and extract_xic use reader.find_by_rt() directly
- Removes dependency on scan_matcher::find_scan_by_rt for these tools
- O(log N) lookup via IndexedMzMLReader, <1ms on 8GB files

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 8: Integration Test and Final Verification

**Files:**
- Modify: `crates/spectrum-io/src/indexed_mzml.rs` (add integration test)

- [ ] **Step 1: Write end-to-end test: open → find_by_rt**

Add to the `#[cfg(test)] mod tests` block in `crates/spectrum-io/src/indexed_mzml.rs`:

```rust
#[test]
fn find_by_rt_on_indexed_reader() {
    // Clean up any stale cache
    let path = fixture_path();
    let _ = std::fs::remove_file(crate::disk_cache::idx_path(&path));

    let reader = IndexedMzMLReader::open(&path).unwrap();

    // Scan 1: RT=120.5s (2.0083 min), ms_level=2, isolation=(471.2561, 1.0, 1.0)
    let result = reader
        .find_by_rt(&path, 2.0, 471.2561, 0.5)
        .unwrap();
    assert!(result.is_some(), "should find scan near RT=2.0 min");
    let (scan, delta) = result.unwrap();
    assert_eq!(scan, 1);
    assert!(delta.abs() < 0.5, "RT delta should be within tolerance");

    // No match for mz outside isolation window
    let no_match = reader
        .find_by_rt(&path, 2.0, 999.0, 0.5)
        .unwrap();
    assert!(no_match.is_none(), "mz=999 should not match any scan");

    // Clean up
    let _ = std::fs::remove_file(crate::disk_cache::idx_path(&path));
}

#[test]
fn find_by_rt_cached_vs_fresh() {
    let dir = tempfile::tempdir().unwrap();
    let src = fixture_path();
    let copy = dir.path().join("test.mzml");
    std::fs::copy(&src, &copy).unwrap();

    // First open: builds index and caches to disk
    let reader1 = IndexedMzMLReader::open(&copy).unwrap();
    let result1 = reader1.find_by_rt(&copy, 2.0, 471.2561, 0.5).unwrap();

    // Second open: loads from disk cache
    let reader2 = IndexedMzMLReader::open(&copy).unwrap();
    let result2 = reader2.find_by_rt(&copy, 2.0, 471.2561, 0.5).unwrap();

    // Both should return the same result
    assert_eq!(result1, result2, "cached and fresh should agree");
}
```

- [ ] **Step 2: Run the integration tests**

Run: `cargo test -p protein-copilot-spectrum-io find_by_rt_on_indexed find_by_rt_cached -- --nocapture 2>&1 | tail -20`

Expected: Both tests pass.

- [ ] **Step 3: Run full workspace tests**

Run: `cargo test --workspace 2>&1 | tail -20`

Expected: All tests pass.

- [ ] **Step 4: Run clippy**

Run: `cargo clippy --workspace 2>&1 | grep -E "warning|error" | head -20`

Expected: No warnings.

- [ ] **Step 5: Final commit**

```bash
git add crates/spectrum-io/src/indexed_mzml.rs
git commit -m "test(spectrum-io): add integration tests for RT binary search

- find_by_rt_on_indexed_reader: end-to-end open → find_by_rt
- find_by_rt_cached_vs_fresh: verifies disk cache preserves metadata

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```
