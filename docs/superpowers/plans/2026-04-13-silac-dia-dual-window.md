# SILAC DIA Dual-Window Extraction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix SILAC+DIA extraction so heavy fragment ions are extracted from the correct MS2 isolation window (which differs from the light window), instead of the current buggy behavior of extracting both from the light window.

**Architecture:** Two-pass scanning strategy. Pass 1 extracts light fragment ions from MS2 scans matching the light precursor's isolation window (existing behavior). Pass 2 (new) discovers the DIA window covering the heavy precursor m/z and extracts heavy fragment ions from those different MS2 scans. Light and heavy MS2 data have independent scan/RT sequences. DDA mode is unchanged (no window filtering).

**Tech Stack:** Rust, protein_copilot_xic crate, protein_copilot_core::spectrum::IsolationWindow

---

### Task 1: Add `heavy_warning` to `XicData` and `ms2_heavy_scans` to `RawScanData`

**Files:**
- Modify: `crates/xic/src/lib.rs:66-88` (XicData struct)
- Modify: `crates/xic/src/lib.rs:175-181` (RawScanData struct)
- Modify: `crates/xic/src/extract.rs:522-533` (XicData construction in extract_xic)
- Modify: `crates/xic/src/extract.rs:898-909` (XicData construction in extract_xic_with_raw)
- Modify: `crates/xic/src/extract.rs:911-914` (RawScanData construction)
- Modify: `crates/report/src/unified_types.rs:62-75` (RawScanData test)
- Modify: `crates/mcp-server/src/tools.rs` (any places constructing RawScanData or XicData)

- [ ] **Step 1: Add `heavy_warning` field to `XicData`**

In `crates/xic/src/lib.rs`, add after line 87 (`extraction_params`):

```rust
// In XicData struct, after extraction_params:
    /// Warning message when heavy MS2 extraction is incomplete (e.g. heavy
    /// precursor m/z outside all DIA windows). `None` when extraction is normal.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub heavy_warning: Option<String>,
```

- [ ] **Step 2: Add `ms2_heavy_scans` field to `RawScanData`**

In `crates/xic/src/lib.rs`, add after line 180 (`ms2_scans`):

```rust
    /// MS2 scans from the heavy precursor's isolation window (DIA+SILAC only).
    /// Empty when not applicable (DDA, no label, or heavy window not found).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub ms2_heavy_scans: Vec<RawScan>,
```

- [ ] **Step 3: Fix all compilation errors**

Add `heavy_warning: None` to every `XicData` construction site:

In `crates/xic/src/extract.rs` line 522 (inside `extract_xic`):
```rust
    Ok(XicData {
        peptide_sequence: peptide_sequence.to_string(),
        target_rt_min: target_rt,
        target_scan,
        charge,
        precursor_mz,
        ms1_precursor_xic,
        ms1_heavy_precursor_xic,
        fragment_xic_traces: fragment_traces,
        heavy_fragment_xic_traces: heavy_traces,
        extraction_params: params.clone(),
        heavy_warning: None,
    })
```

In `crates/xic/src/extract.rs` line 898 (inside `extract_xic_with_raw`):
```rust
    let xic_data = XicData {
        peptide_sequence: peptide_sequence.to_string(),
        target_rt_min: target_rt,
        target_scan,
        charge,
        precursor_mz,
        ms1_precursor_xic,
        ms1_heavy_precursor_xic,
        fragment_xic_traces: fragment_traces,
        heavy_fragment_xic_traces: heavy_traces,
        extraction_params: params.clone(),
        heavy_warning: None,
    };
```

Add `ms2_heavy_scans: Vec::new()` to every `RawScanData` construction:

In `crates/xic/src/extract.rs` line 911:
```rust
    let raw_scans = crate::RawScanData {
        ms1_scans: raw_ms1_trimmed,
        ms2_scans: raw_ms2_trimmed,
        ms2_heavy_scans: Vec::new(),
    };
```

In `crates/report/src/unified_types.rs` test (line 63):
```rust
        let data = RawScanData {
            ms1_scans: vec![RawScan {
                scan_number: 1,
                retention_time_min: 10.0,
                mz_array: vec![450.0],
                intensity_array: vec![1000.0],
            }],
            ms2_scans: vec![],
            ms2_heavy_scans: vec![],
        };
```

- [ ] **Step 4: Build and test**

Run: `cargo test --quiet 2>&1 | tail -5`
Expected: All tests pass, no compilation errors.

- [ ] **Step 5: Commit**

```bash
git add crates/xic/src/lib.rs crates/xic/src/extract.rs crates/report/src/unified_types.rs
git commit -m "feat: add heavy_warning to XicData and ms2_heavy_scans to RawScanData

Prepares data structures for SILAC DIA dual-window extraction.
heavy_warning communicates when heavy MS2 extraction is incomplete.
ms2_heavy_scans stores raw peaks from the heavy isolation window.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: Implement `window_contains_mz()` and `find_heavy_dia_window()` helpers

**Files:**
- Modify: `crates/xic/src/heavy.rs` (add new public functions)
- Test: `crates/xic/src/heavy.rs` (add tests in existing test module)

- [ ] **Step 1: Write failing tests for `window_contains_mz`**

In `crates/xic/src/heavy.rs`, add to the `tests` module:

```rust
    #[test]
    fn window_contains_mz_inside() {
        let w = IsolationWindow {
            target_mz: 500.0,
            lower_offset: 12.5,
            upper_offset: 12.5,
        };
        assert!(window_contains_mz(&w, 500.0));
        assert!(window_contains_mz(&w, 487.5)); // exactly at lower edge
        assert!(window_contains_mz(&w, 512.5)); // exactly at upper edge
        assert!(window_contains_mz(&w, 505.0)); // inside
    }

    #[test]
    fn window_contains_mz_outside() {
        let w = IsolationWindow {
            target_mz: 500.0,
            lower_offset: 12.5,
            upper_offset: 12.5,
        };
        assert!(!window_contains_mz(&w, 487.4)); // just below
        assert!(!window_contains_mz(&w, 512.6)); // just above
        assert!(!window_contains_mz(&w, 600.0)); // far away
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p protein_copilot_xic heavy::tests::window_contains -- 2>&1`
Expected: FAIL — `window_contains_mz` not found.

- [ ] **Step 3: Implement `window_contains_mz`**

In `crates/xic/src/heavy.rs`, add before the `tests` module:

```rust
use protein_copilot_core::spectrum::IsolationWindow;

/// Check if an isolation window contains a given m/z value.
///
/// Returns `true` when `mz` falls within `[target_mz - lower_offset, target_mz + upper_offset]`.
pub fn window_contains_mz(window: &IsolationWindow, mz: f64) -> bool {
    let lo = window.target_mz - window.lower_offset;
    let hi = window.target_mz + window.upper_offset;
    mz >= lo && mz <= hi
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test -p protein_copilot_xic heavy::tests::window_contains -- 2>&1`
Expected: 2 tests PASS.

- [ ] **Step 5: Write failing test for `find_heavy_dia_window`**

This function searches spectra near a reference RT and returns the isolation window that covers a given m/z. Add test:

```rust
    #[test]
    fn find_heavy_dia_window_basic() {
        use protein_copilot_core::spectrum::{PrecursorInfo, Spectrum, MsLevel};

        // Simulate 3 DIA windows in one cycle near RT=10.0 min
        let spectra = vec![
            make_ms2_spec(100, 10.0, 400.0, 12.5), // window [387.5, 412.5]
            make_ms2_spec(101, 10.01, 425.0, 12.5), // window [412.5, 437.5]
            make_ms2_spec(102, 10.02, 450.0, 12.5), // window [437.5, 462.5]
        ];

        // heavy m/z = 440.0 should match window #3 (450 ±12.5)
        let result = find_heavy_dia_window_from_spectra(&spectra, 10.0, 440.0);
        assert!(result.is_some(), "should find window containing 440.0");
        let (scan, win) = result.unwrap();
        assert_eq!(scan, 102);
        assert!((win.target_mz - 450.0).abs() < 0.01);
    }

    #[test]
    fn find_heavy_dia_window_not_found() {
        let spectra = vec![
            make_ms2_spec(100, 10.0, 400.0, 12.5),
            make_ms2_spec(101, 10.01, 425.0, 12.5),
        ];
        // heavy m/z = 600.0 is outside all windows
        let result = find_heavy_dia_window_from_spectra(&spectra, 10.0, 600.0);
        assert!(result.is_none());
    }

    /// Helper to create a mock MS2 spectrum with an isolation window.
    fn make_ms2_spec(scan: u32, rt_min: f64, center_mz: f64, half_width: f64) -> Spectrum {
        Spectrum {
            scan_number: scan,
            ms_level: MsLevel::MS2,
            retention_time_min: rt_min,
            mz_array: vec![],
            intensity_array: vec![],
            precursors: vec![PrecursorInfo {
                mz: center_mz,
                charge: Some(2),
                intensity: None,
                isolation_window: Some(IsolationWindow {
                    target_mz: center_mz,
                    lower_offset: half_width,
                    upper_offset: half_width,
                }),
            }],
        }
    }
```

- [ ] **Step 6: Run tests to verify they fail**

Run: `cargo test -p protein_copilot_xic heavy::tests::find_heavy -- 2>&1`
Expected: FAIL — functions not defined.

- [ ] **Step 7: Implement `find_heavy_dia_window_from_spectra`**

This is a pure function that works on a slice of spectra (no file I/O — the caller provides the spectra). Add to `crates/xic/src/heavy.rs`:

```rust
use protein_copilot_core::spectrum::Spectrum;

/// Search a slice of MS2 spectra for the isolation window covering `target_mz`.
///
/// Finds the MS2 scan closest to `reference_rt_min` whose isolation window
/// contains `target_mz`. Returns `(scan_number, isolation_window)` if found.
///
/// This is a pure helper — the caller streams spectra and filters by RT proximity
/// before calling this function.
pub fn find_heavy_dia_window_from_spectra(
    spectra: &[Spectrum],
    reference_rt_min: f64,
    target_mz: f64,
) -> Option<(u32, IsolationWindow)> {
    spectra
        .iter()
        .filter_map(|spec| {
            let win = spec
                .precursors
                .first()
                .and_then(|p| p.isolation_window.as_ref())?;
            if window_contains_mz(win, target_mz) {
                let rt_delta = (spec.retention_time_min - reference_rt_min).abs();
                Some((spec.scan_number, win.clone(), rt_delta))
            } else {
                None
            }
        })
        .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(scan, win, _)| (scan, win))
}
```

- [ ] **Step 8: Run tests to verify they pass**

Run: `cargo test -p protein_copilot_xic heavy::tests -- 2>&1`
Expected: All heavy tests pass (existing + new).

- [ ] **Step 9: Commit**

```bash
git add crates/xic/src/heavy.rs
git commit -m "feat: add window_contains_mz and find_heavy_dia_window helpers

window_contains_mz checks if an m/z value falls within a DIA isolation window.
find_heavy_dia_window_from_spectra searches nearby spectra to discover which
DIA window covers the heavy precursor m/z — needed for SILAC+DIA dual-window.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: Two-pass extraction in `extract_xic()`

**Files:**
- Modify: `crates/xic/src/extract.rs:232-534` (extract_xic function)

The key change: split `ms2_points` (which combined light+heavy in same scans) into separate `ms2_light_points` and `ms2_heavy_points` collections. Heavy points come from a second streaming pass that matches the heavy isolation window instead of the light window. DDA mode (no target_window) is unchanged — both light and heavy come from the same scans.

- [ ] **Step 1: Write failing test for dual-window behavior**

Add test in `crates/xic/src/extract.rs` test module. This test verifies that in DIA+SILAC mode, heavy traces use different scan numbers than light traces:

```rust
    #[test]
    fn extract_xic_silac_dia_different_scans_unit() {
        // This is a structural test verifying the two-pass logic.
        // We can't easily create DIA mzML in a unit test, so we test
        // the helper functions independently and rely on the integration
        // test for end-to-end validation.
        //
        // Verify that find_heavy_dia_window_from_spectra returns a different
        // scan than the light target when heavy m/z is in a different window.
        use crate::heavy::{find_heavy_dia_window_from_spectra, window_contains_mz};
        use protein_copilot_core::spectrum::*;

        let light_window = IsolationWindow {
            target_mz: 500.0,
            lower_offset: 12.5,
            upper_offset: 12.5,
        };
        let heavy_window = IsolationWindow {
            target_mz: 525.0,
            lower_offset: 12.5,
            upper_offset: 12.5,
        };

        // Light precursor at 505.0 → in light_window [487.5, 512.5]
        assert!(window_contains_mz(&light_window, 505.0));
        assert!(!window_contains_mz(&heavy_window, 505.0));

        // Heavy precursor at 515.0 → in heavy_window [512.5, 537.5]
        assert!(!window_contains_mz(&light_window, 515.0));
        assert!(window_contains_mz(&heavy_window, 515.0));

        // find_heavy should return the heavy window scan
        let spectra = vec![
            Spectrum {
                scan_number: 100,
                ms_level: MsLevel::MS2,
                retention_time_min: 10.0,
                mz_array: vec![],
                intensity_array: vec![],
                precursors: vec![PrecursorInfo {
                    mz: 500.0,
                    charge: Some(2),
                    intensity: None,
                    isolation_window: Some(light_window),
                }],
            },
            Spectrum {
                scan_number: 101,
                ms_level: MsLevel::MS2,
                retention_time_min: 10.01,
                mz_array: vec![],
                intensity_array: vec![],
                precursors: vec![PrecursorInfo {
                    mz: 525.0,
                    charge: Some(2),
                    intensity: None,
                    isolation_window: Some(heavy_window.clone()),
                }],
            },
        ];

        let result = find_heavy_dia_window_from_spectra(&spectra, 10.0, 515.0);
        assert!(result.is_some());
        let (scan, win) = result.unwrap();
        assert_eq!(scan, 101, "heavy scan should be 101, not 100");
        assert!((win.target_mz - 525.0).abs() < 0.01);
    }
```

- [ ] **Step 2: Run test to verify it passes (this tests the helper, not extract_xic yet)**

Run: `cargo test -p protein_copilot_xic tests::extract_xic_silac_dia_different_scans_unit -- 2>&1`
Expected: PASS (depends on Task 2 being done).

- [ ] **Step 3: Refactor `extract_xic()` — split MS2 into light-only and heavy-only**

Replace the MS2 section in `extract_xic()`. The changes are:

**3a. Change the `ms2_points` type** (line 293):

Replace:
```rust
    let mut ms2_points: Vec<(u32, f64, Vec<(f64, Option<f64>)>, Vec<(f64, Option<f64>)>)> = Vec::new();
```
With:
```rust
    // Light MS2 points: (scan, RT, light_intensities)
    let mut ms2_light_points: Vec<(u32, f64, Vec<(f64, Option<f64>)>)> = Vec::new();
    // Heavy MS2 points: (scan, RT, heavy_intensities) — from different DIA window
    let mut ms2_heavy_points: Vec<(u32, f64, Vec<(f64, Option<f64>)>)> = Vec::new();
    // Whether heavy uses a separate DIA window (DIA+SILAC)
    let is_dia = target_window.is_some();
    let needs_separate_heavy_window = is_dia && !heavy_ions.is_empty();
```

**3b. In the MS2 branch** (lines 332-370), collect only light intensities. For DDA (no window), also collect heavy into a combined point:

Replace the entire `MsLevel::MS2` arm:
```rust
            MsLevel::MS2 => {
                let matches_light_window = match (&target_window, spec.precursors.first()) {
                    (Some(tw), Some(prec)) => prec
                        .isolation_window
                        .as_ref()
                        .is_some_and(|w| same_isolation_window(tw, w)),
                    (None, _) => true, // DDA: no window filtering
                    _ => false,
                };

                if matches_light_window {
                    let light_intensities: Vec<(f64, Option<f64>)> = light_ions
                        .iter()
                        .map(|ion| {
                            extract_intensity(
                                ion.mz,
                                &spec.mz_array,
                                &spec.intensity_array,
                                &params.mz_tolerance,
                                params.intensity_rule,
                            )
                        })
                        .collect();

                    ms2_light_points.push((spec.scan_number, rt, light_intensities));

                    // DDA: heavy also from same scans (no separate window)
                    if !needs_separate_heavy_window && !heavy_ions.is_empty() {
                        let heavy_intensities: Vec<(f64, Option<f64>)> = heavy_ions
                            .iter()
                            .map(|ion| {
                                extract_intensity(
                                    ion.mz,
                                    &spec.mz_array,
                                    &spec.intensity_array,
                                    &params.mz_tolerance,
                                    params.intensity_rule,
                                )
                            })
                            .collect();
                        ms2_heavy_points.push((spec.scan_number, rt, heavy_intensities));
                    }
                }

                // DIA+SILAC: check if this scan matches the HEAVY window
                if needs_separate_heavy_window {
                    if let Some(heavy_mz) = heavy_precursor_mz {
                        let matches_heavy = spec
                            .precursors
                            .first()
                            .and_then(|p| p.isolation_window.as_ref())
                            .is_some_and(|w| crate::heavy::window_contains_mz(w, heavy_mz));

                        if matches_heavy {
                            let heavy_intensities: Vec<(f64, Option<f64>)> = heavy_ions
                                .iter()
                                .map(|ion| {
                                    extract_intensity(
                                        ion.mz,
                                        &spec.mz_array,
                                        &spec.intensity_array,
                                        &params.mz_tolerance,
                                        params.intensity_rule,
                                    )
                                })
                                .collect();
                            ms2_heavy_points.push((spec.scan_number, rt, heavy_intensities));
                        }
                    }
                }
            }
```

**3c. Post-processing — separate windowing** (lines 377-397):

Replace:
```rust
    // --- Post-processing ---
    ms2_points.sort_by_key(|(scan, _, _, _)| *scan);
    // ... target_pos, start, end ...
    let windowed = &ms2_points[start..end];
```

With:
```rust
    // --- Post-processing: light windowing ---
    ms2_light_points.sort_by_key(|(scan, _, _)| *scan);

    let target_pos = ms2_light_points
        .iter()
        .position(|(scan, _, _)| *scan == target_scan);
    let (start, end) = match target_pos {
        Some(pos) => {
            let n = params.n_cycles as usize;
            let start = pos.saturating_sub(n);
            let end = (pos + n + 1).min(ms2_light_points.len());
            (start, end)
        }
        None => {
            return Err(XicError::ScanNotFound {
                scan: target_scan,
                path: file_path.to_path_buf(),
            });
        }
    };
    let light_windowed = &ms2_light_points[start..end];

    // --- Heavy windowing (independent scan sequence) ---
    ms2_heavy_points.sort_by_key(|(scan, _, _)| *scan);
    let heavy_windowed: &[(u32, f64, Vec<(f64, Option<f64>)>)] = if ms2_heavy_points.is_empty() {
        &[]
    } else {
        // Find the heavy scan closest to target_scan's RT
        let target_rt_for_heavy = target_rt;
        let heavy_center = ms2_heavy_points
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                let da = (a.1 - target_rt_for_heavy).abs();
                let db = (b.1 - target_rt_for_heavy).abs();
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(idx, _)| idx)
            .unwrap_or(0);
        let n = params.n_cycles as usize;
        let h_start = heavy_center.saturating_sub(n);
        let h_end = (heavy_center + n + 1).min(ms2_heavy_points.len());
        &ms2_heavy_points[h_start..h_end]
    };
```

**3d. Build light traces** (lines 399-420):

Replace `windowed` references with `light_windowed`:
```rust
    let mut fragment_traces: Vec<XicTrace> = light_ions
        .iter()
        .enumerate()
        .map(|(i, ion)| XicTrace {
            ion_label: ion.label.clone(),
            ion_type: ion.ion_type,
            ion_number: ion.ion_number,
            charge: ion.charge,
            theoretical_mz: ion.mz,
            data_points: light_windowed
                .iter()
                .map(|(scan, rt, ints)| XicDataPoint {
                    retention_time_min: *rt,
                    scan_number: *scan,
                    intensity: ints.get(i).map(|(int, _)| *int).unwrap_or(0.0),
                    observed_mz: ints.get(i).and_then(|(_, mz)| *mz),
                })
                .collect(),
            is_heavy: false,
        })
        .collect();
```

**3e. Build heavy traces** (lines 437-463):

Replace with `heavy_windowed`:
```rust
    let heavy_traces: Vec<XicTrace> = if heavy_ions.is_empty() || heavy_windowed.is_empty() {
        Vec::new()
    } else {
        let top_labels: Vec<String> = fragment_traces.iter().map(|t| t.ion_label.clone()).collect();
        heavy_ions
            .iter()
            .enumerate()
            .filter(|(_, ion)| top_labels.contains(&ion.label))
            .map(|(i, ion)| XicTrace {
                ion_label: ion.label.clone(),
                ion_type: ion.ion_type,
                ion_number: ion.ion_number,
                charge: ion.charge,
                theoretical_mz: ion.mz,
                data_points: heavy_windowed
                    .iter()
                    .map(|(scan, rt, heavy_ints)| XicDataPoint {
                        retention_time_min: *rt,
                        scan_number: *scan,
                        intensity: heavy_ints.get(i).map(|(int, _)| *int).unwrap_or(0.0),
                        observed_mz: heavy_ints.get(i).and_then(|(_, mz)| *mz),
                    })
                    .collect(),
                is_heavy: true,
            })
            .collect()
    };
```

**3f. RT range for MS1 trimming — use union of light+heavy**:

Replace the `rt_range` computation (line 466):
```rust
    // RT range: union of light and heavy MS2 windows
    let light_rt_range = if let (Some(first), Some(last)) = (light_windowed.first(), light_windowed.last()) {
        Some((first.1, last.1))
    } else {
        None
    };
    let heavy_rt_range = if let (Some(first), Some(last)) = (heavy_windowed.first(), heavy_windowed.last()) {
        Some((first.1, last.1))
    } else {
        None
    };
    let rt_range = match (light_rt_range, heavy_rt_range) {
        (Some((l_lo, l_hi)), Some((h_lo, h_hi))) => Some((l_lo.min(h_lo), l_hi.max(h_hi))),
        (Some(r), None) | (None, Some(r)) => Some(r),
        (None, None) => None,
    };
```

**3g. Set `heavy_warning`**:

Before constructing `XicData`, add:
```rust
    let heavy_warning = if needs_separate_heavy_window && ms2_heavy_points.is_empty() {
        Some(format!(
            "Heavy precursor m/z ({:.4}) is outside all DIA MS2 isolation windows. Heavy MS2 traces unavailable.",
            heavy_precursor_mz.unwrap_or(0.0)
        ))
    } else {
        None
    };
```

And use `heavy_warning` in the `XicData` construction.

- [ ] **Step 4: Build and run all tests**

Run: `cargo test --quiet 2>&1 | tail -5`
Expected: All tests pass.

- [ ] **Step 5: Commit**

```bash
git add crates/xic/src/extract.rs
git commit -m "feat: two-pass DIA extraction in extract_xic for SILAC

Light and heavy MS2 fragments now extracted from their respective DIA
isolation windows. In DIA+SILAC, heavy fragments come from different
MS2 scans than light. DDA mode unchanged (same scans for both).
Adds heavy_warning when heavy precursor falls outside all DIA windows.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 4: Two-pass extraction in `extract_xic_with_raw()`

**Files:**
- Modify: `crates/xic/src/extract.rs:545-933` (extract_xic_with_raw function)

Apply the same two-pass pattern as Task 3, plus capture raw heavy MS2 scans.

- [ ] **Step 1: Apply identical MS2 split logic**

The changes mirror Task 3 exactly for the MS2 section. Additionally, raw MS2 scans must be split:

**1a. Change `ms2_points` type** (line 618):
```rust
    let mut ms2_light_points: Vec<(u32, f64, Vec<(f64, Option<f64>)>)> = Vec::new();
    let mut ms2_heavy_points: Vec<(u32, f64, Vec<(f64, Option<f64>)>)> = Vec::new();
    let is_dia = target_window.is_some();
    let needs_separate_heavy_window = is_dia && !heavy_ions.is_empty();
```

And split `raw_ms2_scans` into light and heavy:
```rust
    let mut raw_ms2_light_scans: Vec<crate::RawScan> = Vec::new();
    let mut raw_ms2_heavy_scans: Vec<crate::RawScan> = Vec::new();
```

**1b. MS2 branch** (lines 675-732) — same pattern as Task 3 Step 3b, but also capture raw scans into the correct vec:

```rust
            MsLevel::MS2 => {
                let matches_light_window = match (&target_window, spec.precursors.first()) {
                    (Some(tw), Some(prec)) => prec
                        .isolation_window
                        .as_ref()
                        .is_some_and(|w| same_isolation_window(tw, w)),
                    (None, _) => true,
                    _ => false,
                };

                if matches_light_window {
                    let light_intensities: Vec<(f64, Option<f64>)> = light_ions
                        .iter()
                        .map(|ion| {
                            extract_intensity(
                                ion.mz,
                                &spec.mz_array,
                                &spec.intensity_array,
                                &params.mz_tolerance,
                                params.intensity_rule,
                            )
                        })
                        .collect();

                    ms2_light_points.push((spec.scan_number, rt, light_intensities));

                    if !needs_separate_heavy_window && !heavy_ions.is_empty() {
                        let heavy_intensities: Vec<(f64, Option<f64>)> = heavy_ions
                            .iter()
                            .map(|ion| {
                                extract_intensity(
                                    ion.mz,
                                    &spec.mz_array,
                                    &spec.intensity_array,
                                    &params.mz_tolerance,
                                    params.intensity_rule,
                                )
                            })
                            .collect();
                        ms2_heavy_points.push((spec.scan_number, rt, heavy_intensities));
                    }

                    // Capture raw light MS2 peaks
                    let rt_close = (rt - target_rt).abs() < 300.0;
                    if target_window.is_some() || rt_close {
                        raw_ms2_light_scans.push(crate::RawScan {
                            scan_number: spec.scan_number,
                            retention_time_min: rt,
                            mz_array: spec.mz_array.clone(),
                            intensity_array: spec.intensity_array.clone(),
                        });
                    }
                }

                // DIA+SILAC: check heavy window
                if needs_separate_heavy_window {
                    if let Some(heavy_mz) = heavy_precursor_mz {
                        let matches_heavy = spec
                            .precursors
                            .first()
                            .and_then(|p| p.isolation_window.as_ref())
                            .is_some_and(|w| crate::heavy::window_contains_mz(w, heavy_mz));

                        if matches_heavy {
                            let heavy_intensities: Vec<(f64, Option<f64>)> = heavy_ions
                                .iter()
                                .map(|ion| {
                                    extract_intensity(
                                        ion.mz,
                                        &spec.mz_array,
                                        &spec.intensity_array,
                                        &params.mz_tolerance,
                                        params.intensity_rule,
                                    )
                                })
                                .collect();
                            ms2_heavy_points.push((spec.scan_number, rt, heavy_intensities));

                            let rt_close = (rt - target_rt).abs() < 300.0;
                            if rt_close {
                                raw_ms2_heavy_scans.push(crate::RawScan {
                                    scan_number: spec.scan_number,
                                    retention_time_min: rt,
                                    mz_array: spec.mz_array.clone(),
                                    intensity_array: spec.intensity_array.clone(),
                                });
                            }
                        }
                    }
                }
            }
```

**1c. Post-processing — separate windowing** (lines 739-759): Same as Task 3 Step 3c, using `ms2_light_points` / `ms2_heavy_points`.

**1d. Build light and heavy traces**: Same as Task 3 Steps 3d, 3e.

**1e. RT range union**: Same as Task 3 Step 3f.

**1f. Raw scan trimming** (lines 882-896):

Replace:
```rust
    // MS2 raw scans: keep only the windowed ones
    let windowed_scans: std::collections::HashSet<u32> =
        light_windowed.iter().map(|(scan, _, _)| *scan).collect();
    let raw_ms2_trimmed: Vec<crate::RawScan> = raw_ms2_light_scans
        .into_iter()
        .filter(|s| windowed_scans.contains(&s.scan_number))
        .collect();

    let heavy_windowed_scans: std::collections::HashSet<u32> =
        heavy_windowed.iter().map(|(scan, _, _)| *scan).collect();
    let raw_ms2_heavy_trimmed: Vec<crate::RawScan> = raw_ms2_heavy_scans
        .into_iter()
        .filter(|s| heavy_windowed_scans.contains(&s.scan_number))
        .collect();
```

**1g. RawScanData construction** (line 911):
```rust
    let raw_scans = crate::RawScanData {
        ms1_scans: raw_ms1_trimmed,
        ms2_scans: raw_ms2_trimmed,
        ms2_heavy_scans: raw_ms2_heavy_trimmed,
    };
```

**1h. heavy_warning**: Same as Task 3 Step 3g.

- [ ] **Step 2: Build and run all tests**

Run: `cargo test --quiet 2>&1 | tail -5`
Expected: All tests pass.

- [ ] **Step 3: Commit**

```bash
git add crates/xic/src/extract.rs
git commit -m "feat: two-pass DIA extraction in extract_xic_with_raw for SILAC

Mirrors extract_xic dual-window logic, plus captures raw heavy MS2 scans
in RawScanData.ms2_heavy_scans for client-side SILAC recomputation.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 5: Unified HTML template — heavy_warning display

**Files:**
- Modify: `crates/report/templates/unified.html`

- [ ] **Step 1: Add heavy_warning banner**

In the unified.html template, find the section that renders header info (near the title/peptide info). Add a conditional warning block that reads `window.__UNIFIED_DATA__.xic.heavy_warning`:

After the existing header section, add:

```javascript
// In the initialization/setup section of the JavaScript:
if (D.xic && D.xic.heavy_warning) {
    const warn = document.createElement('div');
    warn.style.cssText = 'background:#fff3cd;border:1px solid #ffc107;color:#856404;padding:8px 12px;margin:8px 0;border-radius:4px;font-size:13px;';
    warn.textContent = '⚠️ ' + D.xic.heavy_warning;
    document.getElementById('annotation-plot').parentNode.insertBefore(warn, document.getElementById('annotation-plot'));
}
```

- [ ] **Step 2: Build and verify**

Run: `cargo build --quiet 2>&1`
Expected: Builds successfully.

- [ ] **Step 3: Commit**

```bash
git add crates/report/templates/unified.html
git commit -m "feat: display heavy_warning banner in unified HTML template

Shows a yellow warning when heavy precursor m/z is outside all DIA windows.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 6: Final verification

**Files:** None (test-only)

- [ ] **Step 1: Run full test suite**

Run: `cargo test --quiet 2>&1 | tail -30`
Expected: All 510+ tests pass, 0 failures.

- [ ] **Step 2: Run clippy**

Run: `cargo clippy --all-targets -- -D warnings 2>&1 | tail -20`
Expected: No warnings.

- [ ] **Step 3: Verify existing SILAC+DIA tests**

Run: `cargo test -p protein_copilot_xic -- 2>&1`
Expected: All xic crate tests pass.

- [ ] **Step 4: Verify non-SILAC regression**

Run: `cargo test -p protein_copilot_search_engine -- 2>&1`
Run: `cargo test -p protein_copilot_core -- 2>&1`
Expected: All pass (no regressions).
