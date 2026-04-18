# RT Binary Search Optimization Design

## Problem

`find_scan_by_rt()` calls `reader.read_all()` to load the **entire** mzML file (all spectra including peak data) into memory, just to extract RT/ms_level/isolation_window metadata for scan matching. On an 8GB file this takes ~30 seconds and consumes GB of RAM.

### I/O Benchmark (7.5GB file, SSD)

| Operation | Count | Time | Per-op |
|-----------|-------|------|--------|
| Binary search seeks | 17 | 5 ms | 311 µs |
| RT window reads | 50-100 | 0.6-4 ms | 6-84 µs |
| All-scan header reads | 100K | 7.5 s | 75 µs |
| Sequential full read | 1 | 4.85 s | — |

**Conclusion:** 100K random seeks is slower than reading the entire file sequentially. Binary search (17 seeks + ~50 window reads ≈ 10ms) is 500× faster.

## Design

### Core Idea

Extend the PCIX disk cache to store per-scan metadata (RT, ms_level, isolation_window) alongside byte offsets. On load, build an in-memory sorted RT index for O(log N) binary search. Zero additional I/O for RT lookups.

### PCIX v2 Binary Format

Version 1→2, incompatible (v1 files trigger cache miss → rebuild).

```
Header (25 bytes, unchanged structure):
  [4B magic "PCIX"]
  [1B version = 2]
  [8B source_file_size]
  [8B source_file_mtime]
  [4B entry_count]

Entry (46 bytes each, was 12):
  [4B scan_number: u32]
  [8B byte_offset: u64]
  [8B rt_seconds: f64]
  [1B ms_level: u8]          // 1=MS1, 2=MS2, etc.
  [1B has_isolation: u8]     // 0=no, 1=yes
  [8B target_mz: f64]       // only if has_isolation=1
  [8B lower_offset: f64]    // only if has_isolation=1
  [8B upper_offset: f64]    // only if has_isolation=1
```

When `has_isolation=0`, the three f64 fields are written as 0.0 (fixed layout simplifies parsing).

Size for 100K scans: 25 + 46 × 100,000 = **4.6 MB** (vs 1.2 MB for v1).

### ScanMeta Structure

```rust
/// Per-scan metadata stored in the index.
#[derive(Debug, Clone)]
pub struct ScanMeta {
    pub offset: u64,
    pub rt_seconds: f64,
    pub ms_level: u8,
    pub isolation_window: Option<(f64, f64, f64)>, // (target_mz, lower, upper)
}
```

### ScanIndex Changes

Internal storage changes from `HashMap<u32, u64>` to `HashMap<u32, ScanMeta>`.

New derived structure for binary search:
```rust
/// Pre-sorted (rt_seconds, scan_number) pairs for binary search.
rt_sorted: Vec<(f64, u32)>
```

Built once during `ScanIndex::new()` or via a lazy initializer.

### New API

```rust
impl ScanIndex {
    /// Find the best MS2 scan matching a given RT and precursor m/z.
    ///
    /// Uses binary search on the pre-sorted RT index. O(log N + k) where
    /// k is the number of scans in the RT tolerance window.
    pub fn find_by_rt(
        &self,
        rt_min: f64,          // target RT in minutes
        precursor_mz: f64,    // precursor m/z to match isolation window
        rt_tolerance_min: f64, // RT tolerance in minutes
    ) -> Option<(u32, f64)>;  // (scan_number, rt_delta_min)

    /// Get scan metadata by scan number.
    pub fn get_meta(&self, scan: u32) -> Option<&ScanMeta>;

    /// Get byte offset (backward-compatible convenience).
    pub fn get_offset(&self, scan: u32) -> Option<u64>;
}
```

### Byte-Scan Metadata Extraction

`build_index_by_byte_scan()` currently only extracts scan_number and byte_offset. Extend it to also parse from the same buffer region (typically within 1-2KB after `<spectrum `):

- **RT**: `MS:1000016` (scan start time) — value attribute, convert minutes to seconds if UO:0000031
- **ms_level**: `MS:1000511` (ms level 1) or `MS:1000512` (ms level 2)
- **isolation_window**: `MS:1000827` (target), `MS:1000828` (lower), `MS:1000829` (upper)

These cvParams are always in the spectrum header before `<binaryDataArrayList>`, so the existing 2KB buffer window is sufficient. If a tag's header is truncated at buffer boundary, use defaults (rt=0, ms_level=0, no isolation) and let disk cache preserve whatever was parsed.

### Native Index Metadata

`build_index_from_native_mzml()` only provides byte offsets (the `<indexList>` doesn't contain RT). After native index construction, a single sequential pass reads each scan's header (seek + 2KB read) to fill metadata. For 100K scans this takes ~7.5s on first build, but is persisted to v2 cache for future opens.

**Alternative:** Since byte-scan also fills metadata and takes 3-8s, we could always use byte-scan instead of native index when v2 cache doesn't exist. This avoids the sequential header pass entirely. **Recommended: prefer byte-scan over native index for initial v2 cache build.**

### find_scan_by_rt Rewrite

```rust
// Before (scan_matcher.rs):
pub fn find_scan_by_rt(..., reader: &dyn SpectrumReader) -> Result<u32, ...> {
    let ms2_infos = collect_ms2_info(reader, file)?;  // reads ALL spectra!
    // ...binary search on ms2_infos...
}

// After:
pub fn find_scan_by_rt(..., reader: &dyn SpectrumReader) -> Result<u32, ...> {
    // Cast reader to IndexedMzMLReader to access ScanIndex
    // OR: add find_by_rt to SpectrumReader trait
    // OR: pass ScanIndex directly
}
```

**Decision:** Add `find_by_rt()` to `SpectrumReader` trait with a default implementation that falls back to the current `read_all` approach. `IndexedMzMLReader` overrides with the efficient binary search path.

### match_scans Optimization

`match_scans()` in `scan_matcher.rs` currently calls `collect_ms2_info()` which does `reader.read_all()`. This should be updated to use `ScanIndex` metadata directly when available, eliminating the full file read.

### Data Flow

```
open(file.mzML)
  ├─ Layer 1: Disk cache (.mzml.idx v2) → load ScanMeta + build rt_sorted → done
  ├─ Layer 2: byte-scan (extracts offset + RT + ms_level + isolation in one pass)
  │    └─ save to .idx v2 cache
  └─ ready: ScanIndex with full metadata

find_by_rt(rt=45.3min, mz=500.0, tol=0.5)
  ├─ binary search rt_sorted for [44.8, 45.8] min window
  ├─ filter MS2 only, check isolation_window contains mz=500.0
  ├─ pick closest RT match
  └─ return (scan_number, rt_delta) in <1ms
```

### Performance Target

| Scenario | Before | After |
|----------|--------|-------|
| `find_scan_by_rt()` first call | 30-60s (read_all) | <1ms (cached index) |
| `find_scan_by_rt()` subsequent | 30-60s (read_all each time!) | <1ms |
| `match_scans()` 10K PSMs | 30-60s | <10ms |
| Index build (first open, no cache) | 3-8s (byte-scan) | 3-8s (same, +meta) |
| Index load (cached) | <1ms | <1ms |

### Backward Compatibility

- PCIX v1 files are treated as cache miss (version != 2), triggering rebuild
- `get_offset()` remains available for existing callers
- `SpectrumReader::find_by_rt()` has default impl, so MGF reader doesn't break
- `collect_ms2_info()` / `match_scans()` updated to use new path when reader supports it

### Files Changed

| File | Change |
|------|--------|
| `crates/spectrum-io/src/index.rs` | Add `ScanMeta`, expand `ScanIndex`, add `find_by_rt()`, `rt_sorted` |
| `crates/spectrum-io/src/disk_cache.rs` | PCIX v2 format (read/write), version check |
| `crates/spectrum-io/src/indexed_mzml.rs` | Pass metadata through byte-scan, update `open()` |
| `crates/spectrum-io/src/reader.rs` | Add `find_by_rt()` to `SpectrumReader` trait |
| `crates/result-import/src/scan_matcher.rs` | Rewrite `find_scan_by_rt`, update `match_scans` |
| `crates/mcp-server/src/tools.rs` | Update `annotate_spectrum` / `extract_xic` RT lookup |
