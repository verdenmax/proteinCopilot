# XIC Indexed Extraction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace full-file sequential XIC scans with index-planned targeted reads, reducing annotation time from ~240s to <1s per 7.5 GB mzML file.

**Architecture:** Add `list_scan_meta()` to the `SpectrumReader` trait (indexed override reads from memory), then implement `extract_xic_unified()` that plans reads from metadata and issues ~30 O(1) seek reads instead of streaming 100k+ spectra. Merge the two-pass probe+full pattern into a single call.

**Tech Stack:** Rust, protein-copilot-spectrum-io (indexed readers), protein-copilot-xic (extraction), protein-copilot-mcp-server (tool handlers)

**Spec:** `docs/superpowers/specs/2026-04-27-xic-indexed-extraction-design.md`

---

### Task 1: Add `ScanMetaInfo` struct and `list_scan_meta()` trait method

**Files:**
- Modify: `crates/spectrum-io/src/reader.rs`
- Modify: `crates/spectrum-io/src/indexed_mzml.rs`
- Modify: `crates/spectrum-io/src/indexed_mgf.rs`

- [ ] **Step 1: Add `ScanMetaInfo` struct to reader.rs**

Add a new public struct above the `SpectrumReader` trait definition:

```rust
/// Scan metadata for XIC planning and batch read optimization.
///
/// Contains the minimum metadata needed to plan targeted reads:
/// which scans to read based on ms_level, RT range, and isolation window.
#[derive(Debug, Clone)]
pub struct ScanMetaInfo {
    /// Scan number (1-based).
    pub scan_number: u32,
    /// MS level: 1 = MS1, 2 = MS2.
    pub ms_level: u8,
    /// Retention time in minutes.
    pub rt_min: f64,
    /// Isolation window as (target_mz, lower_offset, upper_offset). None for MS1.
    pub isolation_window: Option<(f64, f64, f64)>,
}
```

- [ ] **Step 2: Add `list_scan_meta()` default method to `SpectrumReader` trait**

Add this method inside the `SpectrumReader` trait block, after `list_ms2_meta()`:

```rust
    /// Returns metadata for all scans (MS1 + MS2).
    ///
    /// Used by XIC extraction to plan targeted reads: identify which scans
    /// match the target isolation window and RT range, then read only those.
    ///
    /// Default implementation streams all spectra (slow for large files).
    /// [`crate::IndexedMzMLReader`] overrides to read from in-memory index
    /// — zero I/O, sub-millisecond even for 100k+ scans.
    fn list_scan_meta(&self, path: &Path) -> Result<Vec<ScanMetaInfo>, SpectrumIoError> {
        use protein_copilot_core::spectrum::MsLevel;

        let mut metas = Vec::new();
        self.for_each_spectrum(path, &mut |spec| {
            let iso = spec.precursors.first().and_then(|p| {
                p.isolation_window
                    .as_ref()
                    .map(|w| (w.target_mz, w.lower_offset, w.upper_offset))
            });
            let ms_level = match spec.ms_level {
                MsLevel::MS1 => 1,
                MsLevel::MS2 => 2,
            };
            metas.push(ScanMetaInfo {
                scan_number: spec.scan_number,
                ms_level,
                rt_min: spec.retention_time_min,
                isolation_window: iso,
            });
            Ok(true)
        })?;
        Ok(metas)
    }
```

- [ ] **Step 3: Override `list_scan_meta()` in `IndexedMzMLReader`**

Add this method inside the `impl SpectrumReader for IndexedMzMLReader` block in `crates/spectrum-io/src/indexed_mzml.rs`, after `list_ms2_meta()`:

```rust
    fn list_scan_meta(
        &self,
        _path: &Path,
    ) -> Result<Vec<crate::reader::ScanMetaInfo>, SpectrumIoError> {
        Ok(self
            .index
            .iter_meta()
            .map(|(&scan, meta)| crate::reader::ScanMetaInfo {
                scan_number: scan,
                ms_level: meta.ms_level,
                rt_min: meta.rt_seconds / 60.0,
                isolation_window: meta.isolation_window,
            })
            .collect())
    }
```

- [ ] **Step 4: Override `list_scan_meta()` in `IndexedMgfReader`**

Add this method inside the `impl SpectrumReader for IndexedMgfReader` block in `crates/spectrum-io/src/indexed_mgf.rs`, after `for_each_spectrum()`:

```rust
    fn list_scan_meta(
        &self,
        _path: &Path,
    ) -> Result<Vec<crate::reader::ScanMetaInfo>, SpectrumIoError> {
        Ok(self
            .index
            .iter_meta()
            .map(|(&scan, meta)| crate::reader::ScanMetaInfo {
                scan_number: scan,
                ms_level: meta.ms_level,
                rt_min: meta.rt_seconds / 60.0,
                isolation_window: meta.isolation_window,
            })
            .collect())
    }
```

- [ ] **Step 5: Build and verify compilation**

Run: `cargo build -p protein-copilot-spectrum-io 2>&1 | tail -5`
Expected: `Finished` with no errors

- [ ] **Step 6: Add test for `list_scan_meta()` in indexed_mzml.rs**

Add this test inside the existing `mod tests` block at the bottom of `crates/spectrum-io/src/indexed_mzml.rs`:

```rust
    #[test]
    fn list_scan_meta_returns_all_scans() {
        let path =
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/small.mzml");
        let reader = IndexedMzMLReader::open(&path).unwrap();
        let metas = reader.list_scan_meta(&path).unwrap();

        // small.mzml has both MS1 and MS2 scans
        assert!(!metas.is_empty(), "should return scan metadata");

        let ms1_count = metas.iter().filter(|m| m.ms_level == 1).count();
        let ms2_count = metas.iter().filter(|m| m.ms_level == 2).count();
        assert!(ms1_count > 0, "should have MS1 scans");
        assert!(ms2_count > 0, "should have MS2 scans");

        // MS2 scans should have isolation windows
        for m in metas.iter().filter(|m| m.ms_level == 2) {
            assert!(
                m.isolation_window.is_some(),
                "MS2 scan {} should have isolation window",
                m.scan_number
            );
        }

        // Total should match index size
        assert_eq!(
            metas.len(),
            reader.index().len(),
            "list_scan_meta count must match index size"
        );
    }
```

- [ ] **Step 7: Run tests**

Run: `cargo test -p protein-copilot-spectrum-io list_scan_meta 2>&1 | tail -10`
Expected: `test ... ok`

- [ ] **Step 8: Commit**

```bash
git add crates/spectrum-io/src/reader.rs crates/spectrum-io/src/indexed_mzml.rs crates/spectrum-io/src/indexed_mgf.rs
git commit -m "feat(spectrum-io): add list_scan_meta() to SpectrumReader trait

Add ScanMetaInfo struct and list_scan_meta() method that returns metadata
for all scans (MS1 + MS2). IndexedMzMLReader and IndexedMgfReader override
with sub-millisecond index reads. Default implementation streams the file.

This enables XIC extraction to plan targeted reads instead of full-file scans."
```

---

### Task 2: Add `XicUnifiedResult` type and re-export

**Files:**
- Modify: `crates/xic/src/lib.rs`

- [ ] **Step 1: Add `XicUnifiedResult` struct to lib.rs**

Add after the `RawScan` struct definition (around line 182):

```rust
/// Combined result from unified XIC extraction.
///
/// Replaces the separate outputs of the old `extract_xic()` and
/// `extract_xic_with_raw()` functions. Contains XIC traces, raw scan
/// data for client-side SILAC recomputation, and ion metadata.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XicUnifiedResult {
    /// XIC trace data (MS1 precursor + MS2 fragment traces, light + heavy).
    pub xic_data: XicData,
    /// Raw peak arrays from scans in the XIC window.
    pub raw_scans: RawScanData,
    /// Ion metadata with K/R counts for client-side SILAC calculation.
    pub ion_metadata: Vec<IonMetadataEntry>,
}
```

- [ ] **Step 2: Build to verify**

Run: `cargo build -p protein-copilot-xic 2>&1 | tail -3`
Expected: `Finished` with no errors

- [ ] **Step 3: Commit**

```bash
git add crates/xic/src/lib.rs
git commit -m "feat(xic): add XicUnifiedResult type

Combined return type for the new extract_xic_unified() function,
replacing separate extract_xic() and extract_xic_with_raw() outputs."
```

---

### Task 3: Implement `extract_xic_unified()`

This is the core optimization. The new function uses index metadata to plan
~30 targeted reads instead of streaming 100k+ spectra.

**Files:**
- Modify: `crates/xic/src/extract.rs`

- [ ] **Step 1: Add the `extract_xic_unified()` function**

Add the following function after `extract_xic_with_raw()` (around line 1163), before `#[cfg(test)]`:

```rust
/// Extract XIC data using index-planned targeted reads.
///
/// Replaces both `extract_xic()` and `extract_xic_with_raw()` with a single
/// call that:
/// 1. Reads target scan metadata via indexed O(1) lookup
/// 2. Queries `list_scan_meta()` to plan which scans to read (sub-ms from index)
/// 3. Reads only the ~30 planned scans via O(1) seek reads
/// 4. Processes intensities and captures raw peaks in one pass
///
/// Performance: ~<1s for a 7.5 GB mzML (vs ~240s with full-file streaming).
#[allow(clippy::too_many_arguments)]
pub fn extract_xic_unified(
    reader: &dyn protein_copilot_spectrum_io::SpectrumReader,
    file_path: &Path,
    target_scan: u32,
    peptide_sequence: &str,
    charge: i32,
    precursor_mz: f64,
    modifications: &[Modification],
    params: &ExtractionParams,
    ms1_mz_window_da: f64,
) -> Result<crate::XicUnifiedResult, XicError> {
    if peptide_sequence.is_empty() {
        return Err(XicError::InvalidPeptide {
            detail: "peptide sequence is empty".to_string(),
        });
    }
    if charge <= 0 {
        return Err(XicError::InvalidPeptide {
            detail: format!("charge must be > 0, got {charge}"),
        });
    }

    // --- Step 1: Read target scan ---
    let target_spectrum = reader.read_spectrum(file_path, target_scan)?;
    if target_spectrum.ms_level != MsLevel::MS2 {
        return Err(XicError::InvalidPeptide {
            detail: format!("scan {} is not MS2 — XIC extraction requires an MS2 scan", target_scan),
        });
    }

    let target_rt = target_spectrum.retention_time_min;
    let target_window = target_spectrum
        .precursors
        .first()
        .and_then(|p| p.isolation_window.as_ref())
        .cloned();

    // --- Build target ion list ---
    let light_ions = build_target_ions(peptide_sequence, modifications, charge);
    let effective_label = params.label_type.as_ref().filter(|label| {
        protein_copilot_core::label::total_heavy_delta(peptide_sequence, label).abs() > 1e-6
    });
    let heavy_ions = match &effective_label {
        Some(label) => crate::heavy::compute_heavy_target_ions(&light_ions, peptide_sequence, label),
        None => Vec::new(),
    };
    let heavy_precursor_mz = effective_label.map(|label| {
        crate::heavy::compute_heavy_precursor_mz(precursor_mz, charge, peptide_sequence, label)
    });

    // DIA detection
    let is_dia = target_window
        .as_ref()
        .map(|w| (w.lower_offset + w.upper_offset) > 1.0)
        .unwrap_or(false);
    let needs_separate_heavy_window = is_dia && !heavy_ions.is_empty();

    // --- Step 2: Query index for all scan metadata ---
    let all_meta = reader.list_scan_meta(file_path)?;

    // --- Step 3: Plan reads ---
    // 3a. MS2 light scans: matching isolation window, sorted by scan number
    let mut ms2_light_meta: Vec<&protein_copilot_spectrum_io::reader::ScanMetaInfo> = all_meta
        .iter()
        .filter(|m| {
            m.ms_level == 2
                && match (&target_window, &m.isolation_window) {
                    (Some(tw), Some(mw)) => {
                        let other = IsolationWindow {
                            target_mz: mw.0,
                            lower_offset: mw.1,
                            upper_offset: mw.2,
                        };
                        same_isolation_window(tw, &other)
                    }
                    (None, _) => true, // DDA: no window filtering
                    _ => false,
                }
        })
        .collect();
    ms2_light_meta.sort_by_key(|m| m.scan_number);

    // Find target position in the light MS2 sequence
    let target_pos = ms2_light_meta
        .iter()
        .position(|m| m.scan_number == target_scan)
        .ok_or_else(|| XicError::ScanNotFound {
            scan: target_scan,
            path: file_path.to_path_buf(),
        })?;

    let n = params.n_cycles as usize;
    let light_start = target_pos.saturating_sub(n);
    let light_end = (target_pos + n + 1).min(ms2_light_meta.len());
    let light_planned: Vec<u32> = ms2_light_meta[light_start..light_end]
        .iter()
        .map(|m| m.scan_number)
        .collect();

    // 3b. MS2 heavy scans (DIA+SILAC): matching heavy isolation window
    let heavy_planned: Vec<u32> = if needs_separate_heavy_window {
        if let Some(heavy_mz) = heavy_precursor_mz {
            let mut ms2_heavy_meta: Vec<&protein_copilot_spectrum_io::reader::ScanMetaInfo> = all_meta
                .iter()
                .filter(|m| {
                    m.ms_level == 2
                        && m.isolation_window.as_ref().is_some_and(|w| {
                            let lo = w.0 - w.1;
                            let hi = w.0 + w.2;
                            heavy_mz >= lo && heavy_mz <= hi
                        })
                })
                .collect();
            ms2_heavy_meta.sort_by_key(|m| m.scan_number);

            if ms2_heavy_meta.is_empty() {
                Vec::new()
            } else {
                // Find center closest to target RT
                let heavy_center = ms2_heavy_meta
                    .iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| {
                        let da = (a.rt_min - target_rt).abs();
                        let db = (b.rt_min - target_rt).abs();
                        da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .map(|(idx, _)| idx)
                    .unwrap_or(0);
                let h_start = heavy_center.saturating_sub(n);
                let h_end = (heavy_center + n + 1).min(ms2_heavy_meta.len());
                ms2_heavy_meta[h_start..h_end]
                    .iter()
                    .map(|m| m.scan_number)
                    .collect()
            }
        } else {
            Vec::new()
        }
    } else {
        Vec::new()
    };

    // 3c. MS1 scans: RT range covering both light and heavy MS2 windows
    let light_rt_range = if light_planned.is_empty() {
        None
    } else {
        let first = ms2_light_meta[light_start].rt_min;
        let last = ms2_light_meta[light_end - 1].rt_min;
        Some((first, last))
    };
    // For heavy, compute from the heavy_planned scan RTs
    let heavy_rt_range = if heavy_planned.is_empty() {
        None
    } else {
        let rts: Vec<f64> = heavy_planned
            .iter()
            .filter_map(|&s| all_meta.iter().find(|m| m.scan_number == s).map(|m| m.rt_min))
            .collect();
        rts.iter()
            .copied()
            .reduce(|a, b| a.min(b))
            .zip(rts.iter().copied().reduce(|a, b| a.max(b)))
    };
    let ms1_rt_range = match (light_rt_range, heavy_rt_range) {
        (Some((l_lo, l_hi)), Some((h_lo, h_hi))) => Some((l_lo.min(h_lo), l_hi.max(h_hi))),
        (Some(r), None) | (None, Some(r)) => Some(r),
        (None, None) => None,
    };

    let ms1_planned: Vec<u32> = match ms1_rt_range {
        Some((rt_lo, rt_hi)) => all_meta
            .iter()
            .filter(|m| m.ms_level == 1 && m.rt_min >= rt_lo && m.rt_min <= rt_hi)
            .map(|m| m.scan_number)
            .collect(),
        None => Vec::new(),
    };

    tracing::info!(
        light_ms2 = light_planned.len(),
        heavy_ms2 = heavy_planned.len(),
        ms1 = ms1_planned.len(),
        total = light_planned.len() + heavy_planned.len() + ms1_planned.len(),
        "XIC indexed: planned targeted reads"
    );

    // Compute dynamic MS1 trim window (same logic as extract_xic_with_raw)
    let k_count = peptide_sequence.chars().filter(|&c| c == 'K' || c == 'k').count() as f64;
    let r_count = peptide_sequence.chars().filter(|&c| c == 'R' || c == 'r').count() as f64;
    let max_heavy_shift =
        (k_count * 8.015 + r_count * 10.009) / (charge.unsigned_abs() as f64).max(1.0);
    let effective_ms1_window = ms1_mz_window_da.max(max_heavy_shift + 5.0);

    // --- Step 4: Read planned scans and extract intensities ---
    let mut ms2_light_points: Vec<Ms2Point> = Vec::new();
    let mut ms2_heavy_points: Vec<Ms2Point> = Vec::new();
    let mut ms1_light_points: Vec<XicDataPoint> = Vec::new();
    let mut ms1_heavy_points: Vec<XicDataPoint> = Vec::new();
    let mut raw_ms1_scans: Vec<crate::RawScan> = Vec::new();
    let mut raw_ms2_light_scans: Vec<crate::RawScan> = Vec::new();
    let mut raw_ms2_heavy_scans: Vec<crate::RawScan> = Vec::new();

    // Read MS1 scans
    for &scan_num in &ms1_planned {
        let spec = match reader.read_spectrum(file_path, scan_num) {
            Ok(s) => s,
            Err(_) => continue,
        };
        let rt = spec.retention_time_min;

        let (light_int, light_obs) = extract_intensity(
            precursor_mz,
            &spec.mz_array,
            &spec.intensity_array,
            &params.mz_tolerance,
            params.intensity_rule,
        );
        ms1_light_points.push(XicDataPoint {
            retention_time_min: rt,
            scan_number: spec.scan_number,
            intensity: light_int,
            observed_mz: light_obs,
        });

        if let Some(heavy_mz) = heavy_precursor_mz {
            let (heavy_int, heavy_obs) = extract_intensity(
                heavy_mz,
                &spec.mz_array,
                &spec.intensity_array,
                &params.mz_tolerance,
                params.intensity_rule,
            );
            ms1_heavy_points.push(XicDataPoint {
                retention_time_min: rt,
                scan_number: spec.scan_number,
                intensity: heavy_int,
                observed_mz: heavy_obs,
            });
        }

        // Capture raw MS1 peaks (trimmed to dynamic window around precursor)
        let (trimmed_mz, trimmed_int) = trim_peaks_to_window(
            &spec.mz_array,
            &spec.intensity_array,
            precursor_mz,
            effective_ms1_window,
        );
        if !trimmed_mz.is_empty() {
            raw_ms1_scans.push(crate::RawScan {
                scan_number: spec.scan_number,
                retention_time_min: rt,
                mz_array: trimmed_mz,
                intensity_array: trimmed_int,
            });
        }
    }

    // Read MS2 light scans
    for &scan_num in &light_planned {
        let spec = match reader.read_spectrum(file_path, scan_num) {
            Ok(s) => s,
            Err(_) => continue,
        };
        let rt = spec.retention_time_min;

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

        // DDA: heavy also from same scans
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

        // Capture raw MS2 light peaks
        raw_ms2_light_scans.push(crate::RawScan {
            scan_number: spec.scan_number,
            retention_time_min: rt,
            mz_array: spec.mz_array.clone(),
            intensity_array: spec.intensity_array.clone(),
        });
    }

    // Read MS2 heavy scans (DIA+SILAC only)
    for &scan_num in &heavy_planned {
        let spec = match reader.read_spectrum(file_path, scan_num) {
            Ok(s) => s,
            Err(_) => continue,
        };
        let rt = spec.retention_time_min;

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

        // Capture raw MS2 heavy peaks
        raw_ms2_heavy_scans.push(crate::RawScan {
            scan_number: spec.scan_number,
            retention_time_min: rt,
            mz_array: spec.mz_array.clone(),
            intensity_array: spec.intensity_array.clone(),
        });
    }

    // --- Step 5: Build XIC traces (same post-processing as before) ---
    ms2_light_points.sort_by_key(|(scan, _, _)| *scan);
    ms2_heavy_points.sort_by_key(|(scan, _, _)| *scan);

    // Build fragment XIC traces from light points
    let mut fragment_traces: Vec<XicTrace> = light_ions
        .iter()
        .enumerate()
        .map(|(i, ion)| XicTrace {
            ion_label: ion.label.clone(),
            ion_type: ion.ion_type,
            ion_number: ion.ion_number,
            charge: ion.charge,
            theoretical_mz: ion.mz,
            data_points: ms2_light_points
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

    // Top-N selection by total intensity
    fragment_traces.sort_by(|a, b| {
        let a_total: f64 = a.data_points.iter().map(|p| p.intensity).sum();
        let b_total: f64 = b.data_points.iter().map(|p| p.intensity).sum();
        b_total.partial_cmp(&a_total).unwrap_or(std::cmp::Ordering::Equal)
    });
    let top_n = params.top_n_ions.min(fragment_traces.len());
    fragment_traces.truncate(top_n);
    fragment_traces.retain(|t| t.data_points.iter().any(|p| p.intensity > 0.0));

    // Build heavy traces matching top-N selection
    let heavy_traces: Vec<XicTrace> = if heavy_ions.is_empty() || ms2_heavy_points.is_empty() {
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
                data_points: ms2_heavy_points
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

    // MS1 precursor XIC traces
    let ms1_precursor_xic = if ms1_light_points.is_empty() {
        None
    } else {
        Some(XicTrace {
            ion_label: "precursor".to_string(),
            ion_type: IonType::Precursor,
            ion_number: 0,
            charge: charge as u32,
            theoretical_mz: precursor_mz,
            data_points: ms1_light_points,
            is_heavy: false,
        })
    };

    let ms1_heavy_precursor_xic = if ms1_heavy_points.is_empty() {
        None
    } else {
        Some(XicTrace {
            ion_label: "precursor (heavy)".to_string(),
            ion_type: IonType::Precursor,
            ion_number: 0,
            charge: charge as u32,
            theoretical_mz: heavy_precursor_mz.unwrap_or(precursor_mz),
            data_points: ms1_heavy_points,
            is_heavy: true,
        })
    };

    let heavy_warning = if needs_separate_heavy_window && ms2_heavy_points.is_empty() {
        Some(format!(
            "Heavy precursor m/z ({:.4}) is outside all DIA MS2 isolation windows. Heavy MS2 traces unavailable.",
            heavy_precursor_mz.unwrap_or(0.0)
        ))
    } else {
        None
    };

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
        heavy_warning,
    };

    let raw_scans = crate::RawScanData {
        ms1_scans: raw_ms1_scans,
        ms2_scans: raw_ms2_light_scans,
        ms2_heavy_scans: raw_ms2_heavy_scans,
    };

    let ion_metadata = compute_ion_metadata(
        &xic_data
            .fragment_xic_traces
            .iter()
            .map(|t| TargetIon {
                label: t.ion_label.clone(),
                ion_type: t.ion_type,
                ion_number: t.ion_number,
                charge: t.charge,
                mz: t.theoretical_mz,
            })
            .collect::<Vec<_>>(),
        peptide_sequence,
    );

    Ok(crate::XicUnifiedResult {
        xic_data,
        raw_scans,
        ion_metadata,
    })
}
```

- [ ] **Step 2: Build to verify compilation**

Run: `cargo build -p protein-copilot-xic 2>&1 | tail -5`
Expected: `Finished` with no errors

- [ ] **Step 3: Commit**

```bash
git add crates/xic/src/extract.rs
git commit -m "feat(xic): implement extract_xic_unified() with index-planned reads

New function uses list_scan_meta() to identify ~30 needed scans, then
reads each via O(1) indexed seek. Merges probe + full extraction into
a single call. Expected speedup: ~240s → <1s per 7.5 GB mzML file."
```

---

### Task 4: Update MCP server `annotate_spectrum` tool

Replace the two-call probe+full pattern with a single `extract_xic_unified()` call,
passing the cached `IndexedMzMLReader`.

**Files:**
- Modify: `crates/mcp-server/src/tools.rs` (lines ~2538-2660)

- [ ] **Step 1: Replace XIC extraction in annotate_spectrum**

Find the block starting at line 2539 (`let mut render_mode = "annotation";`) through
line 2659 (end of the `Err(_) => {` block). Replace the entire `if is_mzml { ... }` block with:

```rust
        let mut render_mode = "annotation";
        if is_mzml {
            let xic_params = protein_copilot_xic::ExtractionParams {
                mz_tolerance: input.extraction_tolerance.unwrap_or(MassTolerance {
                    value: 20.0,
                    unit: ToleranceUnit::Ppm,
                }),
                n_cycles: input.n_cycles.unwrap_or(5),
                top_n_ions: input.top_n_ions.unwrap_or(usize::MAX),
                label_type: input.label_type.clone(),
                intensity_rule: protein_copilot_xic::IntensityRule::MaxInWindow,
            };

            // Get cached indexed reader for O(1) scan lookups
            let cached_reader = self.get_or_create_reader(&spectrum_file)?;

            match protein_copilot_xic::extract::extract_xic_unified(
                cached_reader.as_ref(),
                &spectrum_file,
                resolved_scan,
                &peptide_seq,
                charge,
                annotation.theoretical_mz,
                &modifications,
                &xic_params,
                20.0,
            ) {
                Ok(unified_result) => {
                    // Check if fragment XIC is meaningful (DIA has >1 point, DDA has 1)
                    let xic_meaningful = unified_result
                        .xic_data
                        .fragment_xic_traces
                        .first()
                        .map(|t| t.data_points.len() > 1)
                        .unwrap_or(false);

                    if xic_meaningful {
                        // Refine precursor_mz from MS1 scan
                        let mut annotation = annotation.clone();
                        if let Some(observed) = find_precursor_in_ms1(
                            &unified_result.raw_scans.ms1_scans,
                            annotation.retention_time_min,
                            annotation.theoretical_mz,
                            20.0,
                        ) {
                            annotation.precursor_mz = observed;
                            annotation.delta_mass_ppm = (observed
                                - annotation.theoretical_mz)
                                / annotation.theoretical_mz
                                * 1e6;
                        }

                        // Refine heavy precursor delta from MS1
                        if let Some(ref mut ha) = annotation.heavy_annotation {
                            let heavy_theo = ha.precursor_mz;
                            if let Some(obs_heavy) = find_precursor_in_ms1(
                                &unified_result.raw_scans.ms1_scans,
                                annotation.retention_time_min,
                                heavy_theo,
                                20.0,
                            ) {
                                ha.delta_mass_ppm = Some(
                                    (obs_heavy - heavy_theo) / heavy_theo * 1e6,
                                );
                            } else {
                                ha.delta_mass_ppm = None;
                            }
                        }

                        let unified_data =
                            protein_copilot_report::unified_types::UnifiedViewData {
                                source_file: source_file.clone(),
                                annotation,
                                xic: Some(unified_result.xic_data),
                                raw_scans: Some(unified_result.raw_scans),
                                ion_metadata: unified_result.ion_metadata,
                                peptide_info: make_peptide_info(
                                    &peptide_seq,
                                    charge,
                                    annotation_theo_mz,
                                ),
                            };

                        ReportGenerator::render_unified(
                            &unified_data,
                            &out_path,
                            plotly_mode,
                        )
                        .map_err(|e| {
                            mcp_core_err(protein_copilot_core::error::CoreError::from(e))
                        })?;
                        render_mode = "unified+xic";
                    } else {
                        render_unified_without_xic()?;
                        render_mode = "unified";
                    }
                }
                Err(_) => {
                    render_unified_without_xic()?;
                    render_mode = "unified";
                }
            }
        } else {
            // Non-mzML file — legacy annotation only
            ReportGenerator::render_annotation(&annotation, &out_path)
                .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        }
```

- [ ] **Step 2: Update `extract_xic` tool handler**

Find the `extract_xic` handler (around line 2993). Replace the `extract_xic` call:

```rust
        // Old:
        // let xic_data = protein_copilot_xic::extract::extract_xic(
        //     &file_path, ...
        // )

        // New: use cached indexed reader
        let cached_reader = self.get_or_create_reader(&file_path)?;
        let unified_result = protein_copilot_xic::extract::extract_xic_unified(
            cached_reader.as_ref(),
            &file_path,
            resolved_scan,
            &peptide,
            charge,
            precursor_mz,
            &modifications,
            &params,
            20.0,
        )
        .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?;

        let xic_data = unified_result.xic_data;
```

The rest of the handler (`render_xic_html`, response building) stays the same since
it already uses `xic_data` — the variable name and type are unchanged.

- [ ] **Step 3: Build to verify compilation**

Run: `cargo build -p protein-copilot-mcp-server 2>&1 | tail -10`
Expected: `Finished` with no errors

- [ ] **Step 4: Run all tests**

Run: `cargo test --workspace 2>&1 | tail -20`
Expected: All existing tests pass

- [ ] **Step 5: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "feat(mcp-server): use extract_xic_unified() in annotate/extract tools

Replace two-call probe+full pattern with single extract_xic_unified() call.
Pass cached IndexedMzMLReader for O(1) scan lookups. Single annotation of
7.5 GB mzML drops from ~240s to <1s."
```

---

### Task 5: Clean up old functions and update docs

**Files:**
- Modify: `crates/xic/src/extract.rs`
- Modify: `docs/superpowers/specs/2026-04-27-xic-indexed-extraction-design.md` (mark complete)

- [ ] **Step 1: Mark old functions as deprecated**

Add `#[deprecated]` attributes to both old functions. In `crates/xic/src/extract.rs`:

Before `pub fn extract_xic(`:
```rust
#[deprecated(note = "Use extract_xic_unified() instead — index-planned reads, single call")]
```

Before `pub fn extract_xic_with_raw(`:
```rust
#[deprecated(note = "Use extract_xic_unified() instead — merged into single call")]
```

- [ ] **Step 2: Verify no remaining callers of old functions**

Run: `grep -rn 'extract_xic\b\|extract_xic_with_raw\b' crates/ --include='*.rs' | grep -v 'deprecated\|extract_xic_unified\|#\[test\]\|///\|//!'`
Expected: No output (no non-deprecated, non-unified callers remain)

- [ ] **Step 3: Run full test suite**

Run: `cargo test --workspace 2>&1 | tail -20`
Expected: All tests pass (deprecation warnings are OK)

- [ ] **Step 4: Commit**

```bash
git add crates/xic/src/extract.rs
git commit -m "chore(xic): deprecate extract_xic() and extract_xic_with_raw()

Both functions are superseded by extract_xic_unified() which uses
index-planned targeted reads and merges probe + full extraction."
```

---

### Task 6: Verify end-to-end with real data

This task is manual verification — not automated tests.

- [ ] **Step 1: Build release binary**

Run: `cargo build --release -p protein-copilot-mcp-server 2>&1 | tail -5`
Expected: `Finished` with no errors

- [ ] **Step 2: Test annotate_spectrum on real DIA mzML**

Use the MCP server to annotate one of the previously tested peptides:

```bash
# Start MCP server and send annotate_spectrum request
# Expected: completes in <5s instead of ~120s
```

Verify the output HTML is generated and looks correct compared to the
previously generated annotations in `output/xic_silac/`.

- [ ] **Step 3: Compare output quality**

Open both the new and old annotation HTML files side by side. Verify:
- Same fragment ions matched
- Same XIC trace shapes
- Same SILAC mirror plot structure
- Raw scan data present for client-side recomputation
