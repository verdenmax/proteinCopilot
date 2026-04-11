# RT-Based Scan Lookup for annotate_spectrum & extract_xic

## Problem

Currently `annotate_spectrum` and `extract_xic` require a `scan_number`. When users have external results (DIA-NN, custom JSON) containing only peptide + RT + precursor_mz, they must go through a multi-step import workflow just to obtain scan numbers:

1. Prepare custom JSON → `import_search_results` → run_id
2. `export_results` → psm.tsv with scan numbers
3. `annotate_spectrum(run_id, scan_number)` per PSM

This is unnecessary — the mzML file has all the information needed to resolve RT → scan.

## Solution

Add optional `retention_time_min` parameter to both tools. When provided without `scan_number`, the tool automatically finds the best matching MS2 scan in the mzML file using RT proximity + precursor m/z isolation window filtering.

## Scope

- `annotate_spectrum` Mode 2: accept `retention_time_min` as alternative to `scan_number`
- `extract_xic`: accept `retention_time_min` as alternative to `scan_number`
- Extract reusable `find_scan_by_rt()` from existing `scan_matcher.rs` logic

## Design

### New public function in `result-import/scan_matcher.rs`

```rust
/// Find the best MS2 scan matching a given RT and precursor m/z.
///
/// Algorithm:
/// 1. Read all MS2 spectra info (scan_number, rt, isolation_window)
/// 2. Binary search for RT-proximate candidates within tolerance
/// 3. Filter by precursor_mz falling within isolation window
/// 4. Return scan with smallest RT delta
pub fn find_scan_by_rt(
    file: &Path,
    rt_min: f64,
    precursor_mz: f64,
    rt_tolerance_min: f64,
    reader: &dyn SpectrumReader,
) -> Result<u32, ResultImportError>
```

This reuses the existing `Ms2Info` struct and binary search logic already in `scan_matcher.rs`, but exposed as a single-lookup function.

### MCP Tool Changes

**`annotate_spectrum`** — new optional field:
```
retention_time_min: Option<f64>  // RT in minutes for auto scan lookup
```

Resolution logic:
- If `scan_number > 0`: use it directly (existing behavior)
- If `scan_number == 0` and `retention_time_min` is provided: call `find_scan_by_rt()` with the given RT + precursor_mz (from `peptide_sequence` + `charge`, or explicit precursor_mz)
- Requires `file_path` to be an mzML file (MGF has no isolation windows)

**`extract_xic`** — same new field, same resolution logic.

### Match Strategy

- RT tolerance: 0.5 min (default, matches `import_search_results`)
- When multiple scans match: pick the one with smallest |RT_observed - RT_target|
- precursor_mz must fall within the MS2's isolation window (target ± offsets)
- If no match found: return `NoMatchingScan` error with RT value and tolerance

### Error Cases

| Condition | Error |
|-----------|-------|
| RT provided but file is MGF | "RT-based scan lookup requires mzML (isolation window needed)" |
| No MS2 scan within tolerance | "No MS2 scan found near RT={x} min (±{tol} min) with precursor_mz={mz}" |
| Neither scan_number nor RT | Existing error unchanged |

### What Does NOT Change

- Mode 1 (run_id) — unchanged
- `scan_matcher::match_scans()` batch function — unchanged
- `import_search_results` workflow — still available for bulk imports
