# SILAC DIA Dual-Window Extraction Design

## Problem

In DIA + SILAC experiments, heavy-labeled peptides have shifted precursor m/z
(e.g., K+8.014 Da, R+10.008 Da per residue). DIA isolation windows are
typically 2–25 Da wide, so the heavy precursor often falls in a **different**
MS2 isolation window than the light precursor.

**Current behavior (bug):** Both light and heavy fragment ions are extracted
from MS2 scans matching the **light** peptide's isolation window. Heavy
fragment ions are therefore absent or incorrect — their intensities are 0
because the heavy fragments only exist in the MS2 scan covering the heavy
precursor's window.

**Correct behavior:** Light and heavy MS2 fragment ions must be extracted from
their respective isolation windows. This means they come from **different MS2
scans** with slightly different retention times (same DIA cycle, different
window acquisition time).

## Scope

- **Affected tools:** `extract_xic`, `extract_xic_with_raw`, `annotate_spectrum`
- **Affected mode:** DIA + SILAC only
- **Not affected:** DDA mode, non-SILAC mode, MS1 extraction (full scan covers both)

## Design

### 1. XIC Extraction — Two-Pass Strategy

**In `extract_xic()` and `extract_xic_with_raw()`:**

#### Pass 1 (existing, unchanged for light)
Stream all spectra. For MS2 scans matching the **light** target window:
- Extract light fragment ion intensities → `ms2_light_points`
- Collect raw MS2 peaks (for `_with_raw` variant)

#### Pass 2 (new, for heavy when SILAC + DIA)
Only when `label_type.is_some()` AND `target_window.is_some()` (DIA mode):

1. Compute `heavy_precursor_mz` using `compute_heavy_precursor_mz()`
2. Determine the `heavy_target_window` — the DIA isolation window that covers
   `heavy_precursor_mz`. This is found by scanning MS2 spectra and checking
   which isolation window contains `heavy_precursor_mz`.
3. Stream spectra again. For MS2 scans matching the **heavy** target window:
   - Extract heavy fragment ion intensities → `ms2_heavy_points`
   - Collect raw MS2 peaks tagged as heavy

#### Data Separation

Current `ms2_points` tuple:
```rust
Vec<(u32, f64, Vec<(f64, Option<f64>)>, Vec<(f64, Option<f64>)>)>
//   scan  RT   light_intensities       heavy_intensities
```

New separate collections:
```rust
// Light MS2 points (from light window)
ms2_light_points: Vec<(u32, f64, Vec<(f64, Option<f64>)>)>
//                     scan  RT   light_intensities

// Heavy MS2 points (from heavy window) — different scans!
ms2_heavy_points: Vec<(u32, f64, Vec<(f64, Option<f64>)>)>
//                     scan  RT   heavy_intensities
```

Light and heavy MS2 data points have **independent** RT and scan_number
sequences because they come from different isolation windows in the same DIA
cycle.

#### Windowing (n_cycles)
- Light traces: windowed around `target_scan` in `ms2_light_points` (as today)
- Heavy traces: windowed around the nearest heavy scan to `target_scan`'s RT
  in `ms2_heavy_points`

#### RT Range for MS1 Trimming
Use the union of light and heavy MS2 RT ranges to trim MS1 data points.

### 2. Spectrum Annotation — Dual Scan

**In `annotate_spectrum` handler (tools.rs):**

When `label_type.is_some()` AND the data is DIA:

1. Light annotation: use `target_scan` (user-specified) as today
2. Heavy annotation: find the MS2 scan in the same DIA cycle whose isolation
   window covers `heavy_precursor_mz`
   - Use the `find_heavy_dia_scan()` helper (new function in xic crate)
   - This scans nearby spectra (within ±1 DIA cycle of target_scan's RT) and
     finds MS2 scans whose isolation window covers heavy_precursor_mz
3. Generate two spectrum annotation subplots:
   - Top: light scan annotation (target_scan peaks, light b/y ion matches)
   - Bottom: heavy scan annotation (heavy_scan peaks, heavy b/y ion matches)

### 3. Heavy Window Discovery

New helper function:

```rust
/// Find the DIA isolation window and representative scan covering a given m/z.
///
/// Searches MS2 scans near `reference_rt_min` for one whose isolation window
/// contains `target_mz`. Returns (scan_number, isolation_window) if found.
pub fn find_heavy_dia_scan(
    file_path: &Path,
    reference_rt_min: f64,
    target_mz: f64,
    rt_tolerance_min: f64,  // e.g. 1.0 min (one DIA cycle)
) -> Result<Option<(u32, IsolationWindow)>, XicError>
```

This function:
1. Reads spectra near `reference_rt_min` (within ±rt_tolerance_min)
2. For each MS2 scan, checks if its isolation window contains `target_mz`
3. Returns the closest matching scan (by RT) if found

### 4. Fallback Behavior

If `heavy_precursor_mz` is not covered by any DIA isolation window:
- **XIC:** Skip heavy MS2 traces. MS1 heavy trace still works (no window
  filtering). Add a warning field to `XicData`.
- **Annotation:** Skip heavy subplot. Add a note in HTML: "Heavy precursor
  m/z (XXX.XX) is outside all DIA MS2 isolation windows."

### 5. `XicData` Changes

Add optional warning field:
```rust
pub struct XicData {
    // ... existing fields ...
    /// Warning message when heavy extraction is incomplete.
    pub heavy_warning: Option<String>,
}
```

### 6. DDA Behavior (Unchanged)

In DDA mode (`target_window.is_none()`):
- No isolation window filtering → both light and heavy extracted from same scans
- This is actually correct for DDA: the selected precursor fragmented whatever
  was in the narrow window, and heavy fragments (if present) are in the same
  spectrum
- No changes needed

### 7. Unified HTML Template Updates

The unified.html template and xic.html template need minor updates:
- Heavy MS2 traces now carry their own scan numbers (different from light)
- The `buildCustomData()` helper already handles per-point scan numbers
- No structural HTML changes needed, just the data will be correct
- Add heavy_warning display if present

### 8. `RawScanData` Changes (for `extract_xic_with_raw`)

```rust
pub struct RawScanData {
    pub ms1_scans: Vec<RawScan>,
    pub ms2_scans: Vec<RawScan>,         // light window
    pub ms2_heavy_scans: Vec<RawScan>,   // heavy window (new)
}
```

## Testing

1. **Unit test:** `find_heavy_dia_scan` with mock DIA spectra
2. **Unit test:** Two-pass extraction produces different scan numbers for
   light/heavy MS2 traces
3. **Unit test:** Fallback when heavy m/z is outside all windows
4. **Integration:** End-to-end XIC with SILAC + DIA data
5. **Regression:** Non-SILAC and DDA paths unchanged

## Files to Modify

| File | Change |
|------|--------|
| `crates/xic/src/extract.rs` | Two-pass logic, separate light/heavy MS2 collections |
| `crates/xic/src/lib.rs` | `heavy_warning` field, `RawScanData.ms2_heavy_scans` |
| `crates/xic/src/heavy.rs` | `find_heavy_dia_scan()` helper |
| `crates/mcp-server/src/tools.rs` | Annotation handler dual-scan logic |
| `crates/report/src/unified_visualize.rs` | Dual subplot for annotation |
| `crates/report/templates/unified.html` | Heavy warning display |
