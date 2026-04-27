# XIC Indexed Extraction — Design Spec

## Problem

XIC extraction (`extract_xic` and `extract_xic_with_raw`) scans the **entire** mzML file
sequentially via `for_each_spectrum()`, even though only ~30 scans are needed.
For a 7.5 GB file this takes ~120 seconds per call. `annotate_spectrum` calls
both functions in sequence, resulting in **~240 seconds of redundant I/O** per
annotation.

The `IndexedMzMLReader` already provides:
- Per-scan metadata (RT, ms_level, isolation_window) from a disk-cached index
  → **sub-millisecond** query
- `read_spectrum(scan)` via byte-offset seek → **O(1), milliseconds**

## Solution

Replace the full-file sequential scan with **index-planned targeted reads**,
and merge the two XIC functions into a single call.

## Design

### 1. New trait method: `list_scan_meta()`

**File:** `crates/spectrum-io/src/reader.rs`

```rust
/// Scan metadata for XIC planning and batch read optimization.
pub struct ScanMetaInfo {
    pub scan_number: u32,
    pub ms_level: u8,
    pub rt_min: f64,
    pub isolation_window: Option<(f64, f64, f64)>,
}

// On SpectrumReader trait:
fn list_scan_meta(&self, path: &Path) -> Result<Vec<ScanMetaInfo>, SpectrumIoError> {
    // Default: stream all spectra, extract metadata only
    // (slow for large files — override with index-based implementation)
}
```

**`IndexedMzMLReader` override:** reads from in-memory `ScanIndex::iter_meta()`.
Zero I/O, sub-millisecond even for 100k+ scans.

**`IndexedMgfReader` override:** reads from its in-memory index.

### 2. Merged function: `extract_xic_unified()`

**File:** `crates/xic/src/extract.rs`

```rust
pub struct XicUnifiedResult {
    pub xic_data: XicData,
    pub raw_scans: RawScanData,
    pub ion_metadata: Vec<IonMetadataEntry>,
}

pub fn extract_xic_unified(
    reader: &dyn SpectrumReader,
    file_path: &Path,
    target_scan: u32,
    peptide_sequence: &str,
    charge: i32,
    precursor_mz: f64,
    modifications: &[Modification],
    params: &ExtractionParams,
    ms1_mz_window_da: f64,
) -> Result<XicUnifiedResult, XicError>
```

Replaces both `extract_xic()` and `extract_xic_with_raw()`.

### 3. Indexed extraction algorithm

```
Step 1: reader.read_spectrum(target_scan)
        → target_rt, target_isolation_window

Step 2: reader.list_scan_meta(file_path)
        → Vec<ScanMetaInfo> (sub-ms from index)

Step 3: Plan reads
        ├─ MS2 light scans:
        │   filter: ms_level == 2 AND same_isolation_window(target_window)
        │   sort by scan_number (proxy for time order)
        │   find target position → take [pos - n_cycles, pos + n_cycles]
        │
        ├─ MS2 heavy scans (DIA+SILAC only):
        │   compute heavy_precursor_mz
        │   filter: ms_level == 2 AND window_contains_mz(heavy_mz)
        │   find center closest to target_rt → take ±n_cycles
        │
        └─ MS1 scans:
            compute RT range = union of light/heavy MS2 RT ranges
            filter: ms_level == 1 AND rt_min in [rt_lo, rt_hi]

Step 4: reader.read_spectrum(scan) for each planned scan
        → ~20-30 reads × O(1) seek = <1 second total

Step 5: Process intensities + capture raw peaks (same logic as before)
```

### 4. MCP server changes

**File:** `crates/mcp-server/src/tools.rs`

`annotate_spectrum` handler:
1. Get cached `IndexedMzMLReader` from `reader_cache` (existing LRU cache)
2. Call `extract_xic_unified(reader, ...)` — single call, ~<1s
3. Check if XIC is meaningful: `fragment_xic_traces[0].data_points.len() > 1`
4. Render accordingly (XIC mode or annotation-only mode)

`extract_xic` MCP tool: similarly updated to use `extract_xic_unified()`.

### 5. Old function handling

- `extract_xic()` → deprecated, replaced by `extract_xic_unified()`
- `extract_xic_with_raw()` → removed (logic merged into unified function)
- All call sites in `tools.rs` updated to use the new function
- The old functions can be removed since they are internal crate APIs

### 6. `list_ms2_meta()` remains

`list_ms2_meta()` stays as-is — it's used by DIA extraction and other modules.
`list_scan_meta()` is a superset; `list_ms2_meta()` can be refactored later
to delegate to it, but that's out of scope for this change.

## Performance

| Scenario | Before | After |
|----------|--------|-------|
| Single annotation (7.5 GB mzML) | ~240s (2× full scan) | <1s (30 seeks) |
| Batch 12 annotations (same file) | ~25 min | ~12s |
| XIC probe (DDA, no meaningful XIC) | ~120s (1× full scan) | <1s |

## Testing

1. **Unit test `list_scan_meta()`**: verify indexed reader returns correct metadata
   matching what streaming reader returns for the same file (small test fixture)
2. **Unit test `extract_xic_unified()`**: verify output matches old `extract_xic()`
   for the same inputs on small test fixture
3. **Integration test**: annotate a spectrum from the small.mzml test fixture,
   verify HTML output is generated
4. **Regression**: existing `extract_xic` and `annotate_spectrum` tests continue to pass

## Files Changed

| File | Change |
|------|--------|
| `crates/spectrum-io/src/reader.rs` | Add `ScanMetaInfo` struct, add `list_scan_meta()` to trait with default impl |
| `crates/spectrum-io/src/indexed_mzml.rs` | Override `list_scan_meta()` using index |
| `crates/spectrum-io/src/indexed_mgf.rs` | Override `list_scan_meta()` using index |
| `crates/xic/src/extract.rs` | Add `extract_xic_unified()`, deprecate old functions |
| `crates/xic/src/lib.rs` | Add `XicUnifiedResult` type, re-export |
| `crates/mcp-server/src/tools.rs` | Update `annotate_spectrum` and `extract_xic` tools |
