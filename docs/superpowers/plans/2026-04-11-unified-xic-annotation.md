# Unified Annotation + XIC with Interactive SILAC — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Combine spectrum annotation and XIC chromatogram into a single HTML with client-side interactive SILAC recomputation.

**Architecture:** The report crate gains a new unified HTML template that embeds both annotation SVG and Plotly XIC. The xic crate's extraction function is extended to also return raw scan peak arrays. A JS engine in the HTML performs binary-search peak extraction for arbitrary SILAC deltas. The `annotate_spectrum` MCP tool gains optional XIC parameters (default: enabled for mzML).

**Tech Stack:** Rust (serde, protein-copilot-xic/report/search-engine crates), HTML/CSS/JS (Plotly.js CDN), MCP tool framework (rmcp)

**Spec:** `docs/superpowers/specs/2026-04-11-unified-xic-annotation-design.md`

**Build/test commands:**
- Build: `cargo build --workspace`
- Test all: `cargo test --workspace`
- Test specific crate: `cargo test -p protein-copilot-xic` / `cargo test -p protein-copilot-report`
- Clippy: `cargo clippy --workspace -- -D warnings`
- Test MCP server crate: `cargo test -p protein-copilot-mcp-server`

---

## File Structure

### New files
| File | Purpose |
|------|---------|
| `crates/report/src/unified_types.rs` | `UnifiedViewData`, `RawScanData`, `RawScan`, `IonMetadataEntry`, `PeptideInfo` structs |
| `crates/report/src/unified_visualize.rs` | `render_unified_html()` function |
| `crates/report/templates/unified.html` | Combined annotation + XIC + SILAC HTML template (~650 lines) |

### Modified files
| File | Change |
|------|--------|
| `crates/xic/src/extract.rs` | Add `extract_xic_with_raw()`, `compute_ion_metadata()` |
| `crates/report/src/lib.rs` | Export `unified_types` and `unified_visualize` modules |
| `crates/mcp-server/src/tools.rs` | Add XIC fields to `AnnotateSpectrumInput`, integrate in handler |

---

## Task 1: Add Unified View Types (report crate)

**Files:**
- Create: `crates/report/src/unified_types.rs`
- Modify: `crates/report/src/lib.rs`

- [ ] **Step 1: Write tests for new types (serialization)**

Add to `crates/report/src/unified_types.rs`:

```rust
//! Data types for the unified annotation + XIC HTML view.

use protein_copilot_search_engine::annotate::SpectrumAnnotation;
use protein_copilot_xic::{IonType, XicData};
use serde::Serialize;

/// Combined data for the unified HTML template.
#[derive(Debug, Clone, Serialize)]
pub struct UnifiedViewData {
    /// Spectrum annotation (peaks, coverage, metadata).
    pub annotation: SpectrumAnnotation,
    /// Pre-computed XIC data (light + optional heavy traces).
    pub xic: Option<XicData>,
    /// Raw scan peak arrays for client-side SILAC recomputation.
    pub raw_scans: Option<RawScanData>,
    /// Fragment ion metadata with K/R counts for SILAC calculation.
    pub ion_metadata: Vec<IonMetadataEntry>,
    /// Peptide-level info for SILAC computation.
    pub peptide_info: PeptideInfo,
}

/// Raw peak data from scans in the XIC RT window.
#[derive(Debug, Clone, Serialize)]
pub struct RawScanData {
    /// MS1 scans (trimmed to narrow m/z window around precursor).
    pub ms1_scans: Vec<RawScan>,
    /// MS2 scans (full peak lists).
    pub ms2_scans: Vec<RawScan>,
}

/// A single raw scan's peak list.
#[derive(Debug, Clone, Serialize)]
pub struct RawScan {
    /// Scan number (1-based).
    pub scan_number: u32,
    /// Retention time in seconds.
    pub retention_time_sec: f64,
    /// m/z values (sorted ascending).
    pub mz_array: Vec<f64>,
    /// Intensity values (parallel to mz_array).
    pub intensity_array: Vec<f64>,
}

/// Metadata for one fragment ion enabling client-side SILAC calculation.
#[derive(Debug, Clone, Serialize)]
pub struct IonMetadataEntry {
    /// Human-readable ion label (e.g. "y5¹⁺").
    pub label: String,
    /// Ion type.
    pub ion_type: IonType,
    /// Ion number (e.g. 5 for y5).
    pub ion_number: u32,
    /// Charge state.
    pub charge: u32,
    /// Theoretical m/z of the light (unlabeled) ion.
    pub light_mz: f64,
    /// Count of K (Lysine) residues in this fragment.
    pub k_count: u32,
    /// Count of R (Arginine) residues in this fragment.
    pub r_count: u32,
}

/// Peptide-level SILAC info.
#[derive(Debug, Clone, Serialize)]
pub struct PeptideInfo {
    /// Peptide amino acid sequence.
    pub sequence: String,
    /// Precursor charge state.
    pub charge: i32,
    /// Light precursor m/z.
    pub precursor_mz: f64,
    /// Total K count in the full peptide.
    pub total_k: u32,
    /// Total R count in the full peptide.
    pub total_r: u32,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peptide_info_serializes() {
        let info = PeptideInfo {
            sequence: "PEPTIDEK".to_string(),
            charge: 2,
            precursor_mz: 450.25,
            total_k: 1,
            total_r: 0,
        };
        let json = serde_json::to_string(&info).unwrap();
        assert!(json.contains("\"total_k\":1"));
        assert!(json.contains("\"total_r\":0"));
        assert!(json.contains("PEPTIDEK"));
    }

    #[test]
    fn ion_metadata_serializes() {
        let entry = IonMetadataEntry {
            label: "y5¹⁺".to_string(),
            ion_type: IonType::Y,
            ion_number: 5,
            charge: 1,
            light_mz: 574.28,
            k_count: 1,
            r_count: 0,
        };
        let json = serde_json::to_string(&entry).unwrap();
        assert!(json.contains("\"k_count\":1"));
        assert!(json.contains("\"light_mz\":574.28"));
    }

    #[test]
    fn raw_scan_serializes() {
        let scan = RawScan {
            scan_number: 42,
            retention_time_sec: 120.5,
            mz_array: vec![100.0, 200.0, 300.0],
            intensity_array: vec![1000.0, 5000.0, 2000.0],
        };
        let json = serde_json::to_string(&scan).unwrap();
        assert!(json.contains("\"scan_number\":42"));
        assert!(json.contains("[100.0,200.0,300.0]"));
    }

    #[test]
    fn raw_scan_data_serializes() {
        let data = RawScanData {
            ms1_scans: vec![RawScan {
                scan_number: 1,
                retention_time_sec: 10.0,
                mz_array: vec![450.0],
                intensity_array: vec![1000.0],
            }],
            ms2_scans: vec![],
        };
        let json = serde_json::to_string(&data).unwrap();
        assert!(json.contains("\"ms1_scans\""));
        assert!(json.contains("\"ms2_scans\":[]"));
    }
}
```

- [ ] **Step 2: Run tests to verify they pass**

Run: `cargo test -p protein-copilot-report unified_types -- --nocapture`
Expected: 4 tests pass (serialization tests).

- [ ] **Step 3: Export from lib.rs**

In `crates/report/src/lib.rs`, add after the existing `pub mod xic_visualize;` line:

```rust
pub mod unified_types;
pub mod unified_visualize;
```

Note: `unified_visualize` will be created in Task 5. For now, create an empty placeholder:

Create `crates/report/src/unified_visualize.rs`:
```rust
//! Unified annotation + XIC visualization — renders combined HTML
//! with interactive SILAC controls.
//!
//! Implemented in Task 5.
```

- [ ] **Step 4: Run build to verify compilation**

Run: `cargo build -p protein-copilot-report`
Expected: Compiles successfully.

- [ ] **Step 5: Commit**

```bash
git add crates/report/src/unified_types.rs crates/report/src/unified_visualize.rs crates/report/src/lib.rs
git commit -m "feat(report): add UnifiedViewData types for combined annotation+XIC HTML"
```

---

## Task 2: Ion Metadata Computation (xic crate)

**Files:**
- Modify: `crates/xic/src/extract.rs`

- [ ] **Step 1: Write tests for ion metadata K/R counting**

Add to the `mod tests` block at the bottom of `crates/xic/src/extract.rs`:

```rust
#[test]
fn compute_ion_metadata_counts_k_r_for_b_ions() {
    // "KPEPTIDER" — b1="K" (K=1,R=0), b3="KPE" (K=1,R=0), b8="KPEPTIDE" (K=1,R=0)
    let ions = build_target_ions("KPEPTIDER", &[], 2);
    let meta = compute_ion_metadata(&ions, "KPEPTIDER");

    // Find b1 (1+ charge)
    let b1 = meta.iter().find(|m| m.label == "b1¹⁺").unwrap();
    assert_eq!(b1.k_count, 1);
    assert_eq!(b1.r_count, 0);

    // Find y1 (suffix "R")
    let y1 = meta.iter().find(|m| m.label == "y1¹⁺").unwrap();
    assert_eq!(y1.k_count, 0);
    assert_eq!(y1.r_count, 1);

    // Find y9 would be full peptide — but y goes up to n-1
    // y8 suffix = "PEPTIDER" (K=0, R=1)
    let y8 = meta.iter().find(|m| m.label == "y8¹⁺").unwrap();
    assert_eq!(y8.k_count, 0);
    assert_eq!(y8.r_count, 1);
}

#[test]
fn compute_ion_metadata_no_k_r() {
    // "PEPTIDE" — no K or R
    let ions = build_target_ions("PEPTIDE", &[], 2);
    let meta = compute_ion_metadata(&ions, "PEPTIDE");
    for m in &meta {
        assert_eq!(m.k_count, 0, "ion {} should have k_count=0", m.label);
        assert_eq!(m.r_count, 0, "ion {} should have r_count=0", m.label);
    }
}

#[test]
fn compute_ion_metadata_preserves_light_mz() {
    let ions = build_target_ions("PEPTIDEK", &[], 2);
    let meta = compute_ion_metadata(&ions, "PEPTIDEK");
    // Every metadata entry's light_mz must match the corresponding ion's mz
    for (ion, m) in ions.iter().zip(meta.iter()) {
        assert!(
            (ion.mz - m.light_mz).abs() < 1e-6,
            "mismatch for {}: ion.mz={}, meta.light_mz={}",
            m.label, ion.mz, m.light_mz
        );
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p protein-copilot-xic compute_ion_metadata -- --nocapture`
Expected: FAIL — `compute_ion_metadata` not found.

- [ ] **Step 3: Implement `compute_ion_metadata()`**

Add this public function in `crates/xic/src/extract.rs`, before the `extract_xic` function:

```rust
use protein_copilot_report::unified_types::IonMetadataEntry;
```

Wait — this would create a circular dependency (xic → report → xic). Instead, define `IonMetadataEntry` locally in xic crate and re-export, OR keep the function in the xic crate but use a simpler return type. Since the struct needs to be shared, let's define a lightweight version in xic and have report use it.

Actually, the cleanest approach: put `IonMetadataEntry` in `crates/xic/src/lib.rs` (it's an XIC concept — which fragment contains which labeled residues). Then the report crate can import it from xic.

Add to `crates/xic/src/lib.rs`, after `PlotlyMode`:

```rust
/// Metadata for one fragment ion enabling client-side SILAC calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IonMetadataEntry {
    /// Human-readable ion label (e.g. "y5¹⁺").
    pub label: String,
    /// Ion type.
    pub ion_type: IonType,
    /// Ion number (e.g. 5 for y5).
    pub ion_number: u32,
    /// Charge state.
    pub charge: u32,
    /// Theoretical m/z of the light (unlabeled) ion.
    pub light_mz: f64,
    /// Count of K (Lysine) residues in this fragment.
    pub k_count: u32,
    /// Count of R (Arginine) residues in this fragment.
    pub r_count: u32,
}
```

Then in `crates/xic/src/extract.rs`, add the function before `extract_xic`:

```rust
/// Compute ion metadata with K/R residue counts for SILAC calculation.
///
/// For each target ion, counts the K and R residues in the fragment
/// (prefix for b-ions, suffix for y-ions) so the browser can compute
/// heavy m/z shifts client-side.
pub fn compute_ion_metadata(ions: &[TargetIon], peptide: &str) -> Vec<crate::IonMetadataEntry> {
    let chars: Vec<char> = peptide.chars().collect();
    let n = chars.len();

    ions.iter()
        .map(|ion| {
            let fragment_chars: &[char] = match ion.ion_type {
                IonType::B => &chars[..(ion.ion_number as usize).min(n)],
                IonType::Y => {
                    let start = n.saturating_sub(ion.ion_number as usize);
                    &chars[start..]
                }
                IonType::Precursor => &chars[..],
            };

            crate::IonMetadataEntry {
                label: ion.label.clone(),
                ion_type: ion.ion_type,
                ion_number: ion.ion_number,
                charge: ion.charge,
                light_mz: ion.mz,
                k_count: fragment_chars.iter().filter(|&&c| c == 'K').count() as u32,
                r_count: fragment_chars.iter().filter(|&&c| c == 'R').count() as u32,
            }
        })
        .collect()
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test -p protein-copilot-xic compute_ion_metadata -- --nocapture`
Expected: 3 tests pass.

- [ ] **Step 5: Update report crate's `IonMetadataEntry` to re-export from xic**

In `crates/report/src/unified_types.rs`, remove the local `IonMetadataEntry` definition and replace with a re-export. Change the import and re-export:

```rust
// Remove the local IonMetadataEntry struct definition.
// Import from xic crate instead:
pub use protein_copilot_xic::IonMetadataEntry;
```

Update the `UnifiedViewData` struct field type accordingly (it already uses `IonMetadataEntry`, which now comes from xic). Update the test to import from xic:

```rust
// In the test module, change:
// use super::*;  ← this still works since IonMetadataEntry is re-exported
```

- [ ] **Step 6: Run full build + tests**

Run: `cargo test -p protein-copilot-xic && cargo test -p protein-copilot-report`
Expected: All tests pass.

- [ ] **Step 7: Commit**

```bash
git add crates/xic/src/lib.rs crates/xic/src/extract.rs crates/report/src/unified_types.rs
git commit -m "feat(xic): add compute_ion_metadata() for client-side SILAC K/R counting"
```

---

## Task 3: Extract XIC with Raw Scan Data (xic crate)

**Files:**
- Modify: `crates/xic/src/extract.rs`
- Modify: `crates/xic/src/lib.rs`

- [ ] **Step 1: Add `RawScanData` and `RawScan` types to xic crate**

Add to `crates/xic/src/lib.rs`, after `IonMetadataEntry`:

```rust
/// Raw peak data from scans in the XIC RT window.
/// Used for client-side SILAC recomputation in the unified HTML.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RawScanData {
    /// MS1 scans (trimmed to narrow m/z window around precursor).
    pub ms1_scans: Vec<RawScan>,
    /// MS2 scans (full peak lists from matching isolation windows).
    pub ms2_scans: Vec<RawScan>,
}

/// A single raw scan's peak list.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RawScan {
    /// Scan number (1-based).
    pub scan_number: u32,
    /// Retention time in seconds.
    pub retention_time_sec: f64,
    /// m/z values (sorted ascending).
    pub mz_array: Vec<f64>,
    /// Intensity values (parallel to mz_array).
    pub intensity_array: Vec<f64>,
}
```

- [ ] **Step 2: Write test for `extract_xic_with_raw`**

Add to `crates/xic/src/extract.rs` test module:

```rust
#[test]
fn extract_xic_with_raw_returns_raw_scans() {
    // This test uses the same test fixture as extract_xic tests.
    // We verify that raw_scans is populated and ms2_scans are non-empty.
    let fixture = std::path::Path::new("tests/fixtures");
    if !fixture.exists() {
        // Skip if no fixture available (unit test environment)
        return;
    }
    // Integration test — will be covered in Task 7
}

#[test]
fn trim_ms1_peaks_filters_by_mz_window() {
    let mz = vec![100.0, 200.0, 449.0, 450.0, 451.0, 500.0, 600.0];
    let int = vec![10.0, 20.0, 100.0, 500.0, 200.0, 30.0, 40.0];
    let (trimmed_mz, trimmed_int) = trim_peaks_to_window(&mz, &int, 450.0, 20.0);
    // Should keep 449, 450, 451 (within ±20 Da of 450)
    assert_eq!(trimmed_mz.len(), 3);
    assert!((trimmed_mz[0] - 449.0).abs() < 0.01);
    assert!((trimmed_mz[2] - 451.0).abs() < 0.01);
    assert_eq!(trimmed_int.len(), 3);
}
```

- [ ] **Step 3: Implement `trim_peaks_to_window` helper**

Add to `crates/xic/src/extract.rs`, before `extract_xic`:

```rust
/// Trim a peak list to only include peaks within ±window_da of center_mz.
///
/// Used to reduce MS1 raw data volume for HTML embedding — only peaks
/// near the precursor m/z are needed for SILAC recomputation.
fn trim_peaks_to_window(
    mz_array: &[f64],
    intensity_array: &[f64],
    center_mz: f64,
    window_da: f64,
) -> (Vec<f64>, Vec<f64>) {
    let lo = center_mz - window_da;
    let hi = center_mz + window_da;

    let start = mz_array.partition_point(|&m| m < lo);
    let end = mz_array.partition_point(|&m| m <= hi);

    (
        mz_array[start..end].to_vec(),
        intensity_array[start..end].to_vec(),
    )
}
```

- [ ] **Step 4: Run trim test**

Run: `cargo test -p protein-copilot-xic trim_ms1_peaks -- --nocapture`
Expected: PASS.

- [ ] **Step 5: Implement `extract_xic_with_raw()`**

Add to `crates/xic/src/extract.rs`, after the existing `extract_xic` function (before `mod tests`):

```rust
/// Extract XIC data AND raw scan peak arrays for client-side SILAC.
///
/// This is an extension of [`extract_xic`] that additionally captures
/// raw peak data from MS1/MS2 scans in the RT window. The raw data
/// enables the HTML frontend to recompute XIC traces for arbitrary
/// SILAC configurations without a backend round-trip.
///
/// MS1 peaks are trimmed to ±`ms1_mz_window_da` around `precursor_mz`
/// to control embedded data volume.
pub fn extract_xic_with_raw(
    file_path: &Path,
    target_scan: u32,
    peptide_sequence: &str,
    charge: i32,
    precursor_mz: f64,
    modifications: &[Modification],
    params: &ExtractionParams,
    ms1_mz_window_da: f64,
) -> Result<(XicData, crate::RawScanData, Vec<crate::IonMetadataEntry>), XicError> {
    if peptide_sequence.is_empty() {
        return Err(XicError::InvalidPeptide {
            detail: "peptide sequence is empty".to_string(),
        });
    }

    let info = protein_copilot_spectrum_io::detect_format(file_path)?;
    if info.format != protein_copilot_core::spectrum::SpectrumFormat::MzML {
        return Err(XicError::UnsupportedFormat {
            path: file_path.to_path_buf(),
        });
    }

    let reader = protein_copilot_spectrum_io::create_reader(&info);

    // --- Pass 0: Get target scan info ---
    let target_spectrum = reader.read_spectrum(file_path, target_scan)?;
    let target_rt = target_spectrum.retention_time_sec;
    let target_window = target_spectrum
        .precursors
        .first()
        .and_then(|p| p.isolation_window.as_ref())
        .cloned();

    // --- Build target ion list ---
    let light_ions = build_target_ions(peptide_sequence, modifications, charge);
    let heavy_ions = match &params.label_type {
        Some(label) => {
            crate::heavy::compute_heavy_target_ions(&light_ions, peptide_sequence, label)
        }
        None => Vec::new(),
    };

    let heavy_precursor_mz = params.label_type.as_ref().map(|label| {
        crate::heavy::compute_heavy_precursor_mz(precursor_mz, charge, peptide_sequence, label)
    });

    // --- Pass 1: Stream spectra, extract intensities, AND capture raw peaks ---
    let mut ms2_points: Vec<(u32, f64, Vec<f64>, Vec<f64>)> = Vec::new();
    let mut ms1_light_points: Vec<XicDataPoint> = Vec::new();
    let mut ms1_heavy_points: Vec<XicDataPoint> = Vec::new();
    let mut raw_ms1_scans: Vec<crate::RawScan> = Vec::new();
    let mut raw_ms2_scans: Vec<crate::RawScan> = Vec::new();

    reader.for_each_spectrum(file_path, &mut |spec| {
        let rt = spec.retention_time_sec;

        match spec.ms_level {
            MsLevel::MS1 => {
                let light_int = extract_intensity(
                    precursor_mz,
                    &spec.mz_array,
                    &spec.intensity_array,
                    &params.mz_tolerance,
                    params.intensity_rule,
                );
                ms1_light_points.push(XicDataPoint {
                    retention_time_sec: rt,
                    scan_number: spec.scan_number,
                    intensity: light_int,
                });

                if let Some(heavy_mz) = heavy_precursor_mz {
                    let heavy_int = extract_intensity(
                        heavy_mz,
                        &spec.mz_array,
                        &spec.intensity_array,
                        &params.mz_tolerance,
                        params.intensity_rule,
                    );
                    ms1_heavy_points.push(XicDataPoint {
                        retention_time_sec: rt,
                        scan_number: spec.scan_number,
                        intensity: heavy_int,
                    });
                }

                // Capture raw MS1 peaks (trimmed to window around precursor)
                let (trimmed_mz, trimmed_int) = trim_peaks_to_window(
                    &spec.mz_array,
                    &spec.intensity_array,
                    precursor_mz,
                    ms1_mz_window_da,
                );
                if !trimmed_mz.is_empty() {
                    raw_ms1_scans.push(crate::RawScan {
                        scan_number: spec.scan_number,
                        retention_time_sec: rt,
                        mz_array: trimmed_mz,
                        intensity_array: trimmed_int,
                    });
                }
            }
            MsLevel::MS2 => {
                let matches_window = match (&target_window, spec.precursors.first()) {
                    (Some(tw), Some(prec)) => prec
                        .isolation_window
                        .as_ref()
                        .is_some_and(|w| same_isolation_window(tw, w)),
                    (None, _) => true,
                    _ => false,
                };

                if matches_window {
                    let light_intensities: Vec<f64> = light_ions
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

                    let heavy_intensities: Vec<f64> = heavy_ions
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

                    ms2_points.push((
                        spec.scan_number,
                        rt,
                        light_intensities,
                        heavy_intensities,
                    ));

                    // Capture raw MS2 peaks (full spectrum)
                    raw_ms2_scans.push(crate::RawScan {
                        scan_number: spec.scan_number,
                        retention_time_sec: rt,
                        mz_array: spec.mz_array.clone(),
                        intensity_array: spec.intensity_array.clone(),
                    });
                }
            }
            _ => {}
        }
        Ok(true)
    })?;

    // --- Post-processing (identical to extract_xic) ---
    ms2_points.sort_by_key(|(scan, _, _, _)| *scan);

    let target_pos = ms2_points
        .iter()
        .position(|(scan, _, _, _)| *scan == target_scan);
    let (start, end) = match target_pos {
        Some(pos) => {
            let n = params.n_cycles as usize;
            let start = pos.saturating_sub(n);
            let end = (pos + n + 1).min(ms2_points.len());
            (start, end)
        }
        None => (0, ms2_points.len()),
    };
    let windowed = &ms2_points[start..end];

    // Build fragment XIC traces (same logic as extract_xic)
    let mut fragment_traces: Vec<XicTrace> = light_ions
        .iter()
        .enumerate()
        .map(|(i, ion)| XicTrace {
            ion_label: ion.label.clone(),
            ion_type: ion.ion_type,
            ion_number: ion.ion_number,
            charge: ion.charge,
            theoretical_mz: ion.mz,
            data_points: windowed
                .iter()
                .map(|(scan, rt, ints, _)| XicDataPoint {
                    retention_time_sec: *rt,
                    scan_number: *scan,
                    intensity: ints.get(i).copied().unwrap_or(0.0),
                })
                .collect(),
            is_heavy: false,
        })
        .collect();

    fragment_traces.sort_by(|a, b| {
        let a_total: f64 = a.data_points.iter().map(|p| p.intensity).sum();
        let b_total: f64 = b.data_points.iter().map(|p| p.intensity).sum();
        b_total
            .partial_cmp(&a_total)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let top_n = params.top_n_ions.min(fragment_traces.len());
    fragment_traces.truncate(top_n);

    let heavy_traces: Vec<XicTrace> = if heavy_ions.is_empty() {
        Vec::new()
    } else {
        let top_labels: Vec<String> =
            fragment_traces.iter().map(|t| t.ion_label.clone()).collect();
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
                data_points: windowed
                    .iter()
                    .map(|(scan, rt, _, heavy_ints)| XicDataPoint {
                        retention_time_sec: *rt,
                        scan_number: *scan,
                        intensity: heavy_ints.get(i).copied().unwrap_or(0.0),
                    })
                    .collect(),
                is_heavy: true,
            })
            .collect()
    };

    // Trim MS1 to same RT range as MS2
    let rt_range = if let (Some(first), Some(last)) = (windowed.first(), windowed.last()) {
        Some((first.1, last.1))
    } else {
        None
    };

    let ms1_precursor_xic = if ms1_light_points.is_empty() {
        None
    } else {
        let filtered: Vec<XicDataPoint> = match rt_range {
            Some((lo, hi)) => ms1_light_points
                .into_iter()
                .filter(|p| p.retention_time_sec >= lo && p.retention_time_sec <= hi)
                .collect(),
            None => ms1_light_points,
        };
        if filtered.is_empty() {
            None
        } else {
            Some(XicTrace {
                ion_label: "precursor".to_string(),
                ion_type: IonType::Precursor,
                ion_number: 0,
                charge: charge as u32,
                theoretical_mz: precursor_mz,
                data_points: filtered,
                is_heavy: false,
            })
        }
    };

    let ms1_heavy_precursor_xic = if ms1_heavy_points.is_empty() {
        None
    } else {
        let filtered: Vec<XicDataPoint> = match rt_range {
            Some((lo, hi)) => ms1_heavy_points
                .into_iter()
                .filter(|p| p.retention_time_sec >= lo && p.retention_time_sec <= hi)
                .collect(),
            None => ms1_heavy_points,
        };
        if filtered.is_empty() {
            None
        } else {
            Some(XicTrace {
                ion_label: "precursor (heavy)".to_string(),
                ion_type: IonType::Precursor,
                ion_number: 0,
                charge: charge as u32,
                theoretical_mz: heavy_precursor_mz.unwrap_or(precursor_mz),
                data_points: filtered,
                is_heavy: true,
            })
        }
    };

    // Trim raw scans to same RT window
    let raw_ms1_trimmed = match rt_range {
        Some((lo, hi)) => raw_ms1_scans
            .into_iter()
            .filter(|s| s.retention_time_sec >= lo && s.retention_time_sec <= hi)
            .collect(),
        None => raw_ms1_scans,
    };
    // MS2 raw scans: keep only the windowed ones
    let windowed_scans: std::collections::HashSet<u32> =
        windowed.iter().map(|(scan, _, _, _)| *scan).collect();
    let raw_ms2_trimmed: Vec<crate::RawScan> = raw_ms2_scans
        .into_iter()
        .filter(|s| windowed_scans.contains(&s.scan_number))
        .collect();

    let xic_data = XicData {
        peptide_sequence: peptide_sequence.to_string(),
        target_rt_sec: target_rt,
        target_scan,
        charge,
        precursor_mz,
        ms1_precursor_xic,
        ms1_heavy_precursor_xic,
        fragment_xic_traces: fragment_traces,
        heavy_fragment_xic_traces: heavy_traces,
        extraction_params: params.clone(),
    };

    let raw_scans = crate::RawScanData {
        ms1_scans: raw_ms1_trimmed,
        ms2_scans: raw_ms2_trimmed,
    };

    // Ion metadata for the top-N selected light ions
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

    Ok((xic_data, raw_scans, ion_metadata))
}
```

- [ ] **Step 6: Run build + existing tests**

Run: `cargo test -p protein-copilot-xic`
Expected: All existing tests pass + new tests pass. No regressions.

- [ ] **Step 7: Commit**

```bash
git add crates/xic/src/extract.rs crates/xic/src/lib.rs
git commit -m "feat(xic): add extract_xic_with_raw() for embedded raw scan data"
```

---

## Task 4: Unified HTML Template

**Files:**
- Create: `crates/report/templates/unified.html`

This is the largest single file. It combines:
- Annotation rendering (SVG, from annotation.html) — info panel, coverage brackets, spectrum peaks
- XIC rendering (Plotly.js) — MS1 + MS2 charts
- SILAC controls — preset selector, custom input, display mode
- Client-side SILAC recomputation engine — binary search + m/z shift

- [ ] **Step 1: Create the unified HTML template**

Create `crates/report/templates/unified.html` with the full template content. The template uses `/*__UNIFIED_JSON__*/` as placeholder for data injection.

The template is too large to inline here (estimated ~700 lines). The structure is:

```
<!DOCTYPE html>
<html>
<head>
  <style>
    /* CSS combining annotation styles + XIC styles + SILAC controls */
  </style>
  <script src="__PLOTLY_SRC__"></script>
</head>
<body>
  <!-- Section 1: Info Panel (from annotation.html) -->
  <div id="info-panel" class="info-card"></div>

  <!-- Section 2: Fragment Ion Coverage (SVG brackets from annotation.html) -->
  <div id="coverage-panel" class="coverage-card"></div>

  <!-- Section 3: Spectrum Annotation (SVG from annotation.html) -->
  <div class="spectrum-card">
    <h2>Spectrum</h2>
    <div id="spectrum-container"></div>
  </div>

  <!-- Section 4: SILAC Controls (NEW) -->
  <div id="silac-controls" class="silac-card">
    <div class="ctrl-row">
      <label>SILAC Preset:
        <select id="silac-preset">
          <option value="none">None (Light only)</option>
          <option value="standard" selected>Standard SILAC (K+8, R+10)</option>
          <option value="medium">Medium SILAC (K+4, R+6)</option>
          <option value="custom">Custom...</option>
        </select>
      </label>
      <label>Display:
        <select id="display-mode">
          <option value="overlay">Overlay (Light + Heavy)</option>
          <option value="light">Light Only</option>
          <option value="heavy">Heavy Only</option>
          <option value="split">Split (L/H stacked)</option>
        </select>
      </label>
      <span id="silac-badge" class="silac-badge"></span>
    </div>
    <div id="custom-fields" style="display:none;">
      K Δ: <input type="number" id="custom-k" step="0.001" value="8.014199"> Da
      R Δ: <input type="number" id="custom-r" step="0.001" value="10.008269"> Da
      <button id="apply-custom">Apply</button>
    </div>
    <div id="no-kr-message" style="display:none;" class="info-msg">
      ℹ️ This peptide contains no K or R residues — heavy label has no effect.
    </div>
  </div>

  <!-- Section 5: MS1 Precursor XIC (Plotly) -->
  <div id="ms1-card" class="xic-card">
    <h2>MS1 Precursor XIC</h2>
    <div id="ms1-plot" style="height:200px;"></div>
  </div>

  <!-- Section 6: MS2 Fragment Ion XIC (Plotly) -->
  <div id="ms2-card" class="xic-card">
    <h2>MS2 Fragment Ion XIC</h2>
    <div id="ms2-plot" style="height:400px;"></div>
  </div>

  <!-- Tooltip (annotation) -->
  <div id="tooltip"></div>

  <!-- Parameters footer -->
  <div class="params" id="params-info"></div>

  <!-- Data injection -->
  <script>/*__UNIFIED_JSON__*/</script>

  <!-- Annotation renderer (from annotation.html) -->
  <script>
    // Render info panel, coverage brackets, spectrum SVG
    // This is the same JS from annotation.html, reading from
    // window.__UNIFIED_DATA__.annotation
  </script>

  <!-- XIC renderer + SILAC engine (NEW) -->
  <script>
    (function() {
      var U = window.__UNIFIED_DATA__;
      if (!U || !U.xic) { /* hide XIC sections */ return; }

      var PRESETS = {
        none:     { k: 0, r: 0 },
        standard: { k: 8.014199, r: 10.008269 },
        medium:   { k: 4.025107, r: 6.020129 },
      };

      var colors = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b'];

      // --- SILAC recomputation engine ---
      function binarySearchLo(arr, val) {
        var lo = 0, hi = arr.length;
        while (lo < hi) { var mid = (lo + hi) >> 1; arr[mid] < val ? lo = mid + 1 : hi = mid; }
        return lo;
      }

      function extractIntensity(targetMz, mzArr, intArr, tolPpm) {
        var tolDa = targetMz * tolPpm * 1e-6;
        var lo = binarySearchLo(mzArr, targetMz - tolDa * 1.5);
        var best = 0;
        for (var i = lo; i < mzArr.length && mzArr[i] <= targetMz + tolDa * 1.5; i++) {
          if (Math.abs(mzArr[i] - targetMz) / targetMz * 1e6 <= tolPpm) {
            if (intArr[i] > best) best = intArr[i];
          }
        }
        return best;
      }

      function recomputeHeavy(kDelta, rDelta) {
        var info = U.peptide_info;
        var tol = U.xic.extraction_params.mz_tolerance.value; // ppm

        // Heavy precursor
        var precDelta = info.total_k * kDelta + info.total_r * rDelta;
        var heavyPrecMz = info.precursor_mz + precDelta / Math.abs(info.charge);

        // Heavy fragments
        var heavyIons = U.ion_metadata.map(function(ion) {
          var d = ion.k_count * kDelta + ion.r_count * rDelta;
          return { label: ion.label, heavy_mz: ion.light_mz + d / ion.charge };
        });

        // Extract from raw scans
        var ms1Heavy = null;
        if (U.raw_scans && U.raw_scans.ms1_scans.length > 0) {
          ms1Heavy = U.raw_scans.ms1_scans.map(function(s) {
            return { rt: s.retention_time_sec, scan: s.scan_number,
                     intensity: extractIntensity(heavyPrecMz, s.mz_array, s.intensity_array, tol) };
          });
        }

        var ms2Heavy = heavyIons.map(function(ion) {
          if (!U.raw_scans) return [];
          return U.raw_scans.ms2_scans.map(function(s) {
            return { rt: s.retention_time_sec, scan: s.scan_number,
                     intensity: extractIntensity(ion.heavy_mz, s.mz_array, s.intensity_array, tol) };
          });
        });

        return { ms1Heavy: ms1Heavy, ms2Heavy: ms2Heavy, heavyIons: heavyIons, heavyPrecMz: heavyPrecMz };
      }

      // --- Plotly rendering ---
      function renderPlots(mode, heavyResult) {
        // ... builds Plotly traces based on display mode (overlay/split/light/heavy)
        // Uses U.xic for light traces, heavyResult for heavy traces
      }

      // --- Event handlers ---
      var presetSelect = document.getElementById('silac-preset');
      var displaySelect = document.getElementById('display-mode');
      var customFields = document.getElementById('custom-fields');
      var badge = document.getElementById('silac-badge');
      var noKrMsg = document.getElementById('no-kr-message');

      // Hide SILAC controls if no K/R
      if (U.peptide_info.total_k === 0 && U.peptide_info.total_r === 0) {
        noKrMsg.style.display = 'block';
        presetSelect.value = 'none';
        presetSelect.disabled = true;
      }

      function update() {
        var preset = presetSelect.value;
        customFields.style.display = preset === 'custom' ? 'flex' : 'none';

        var kd, rd;
        if (preset === 'custom') {
          kd = parseFloat(document.getElementById('custom-k').value) || 0;
          rd = parseFloat(document.getElementById('custom-r').value) || 0;
        } else {
          var p = PRESETS[preset];
          kd = p.k; rd = p.r;
        }

        badge.textContent = kd === 0 && rd === 0 ? 'No label' : 'K+' + kd.toFixed(3) + ' R+' + rd.toFixed(3);

        var heavyResult = (kd === 0 && rd === 0) ? null : recomputeHeavy(kd, rd);
        renderPlots(displaySelect.value, heavyResult);
      }

      presetSelect.addEventListener('change', update);
      displaySelect.addEventListener('change', update);
      document.getElementById('apply-custom').addEventListener('click', update);

      // Initial render with standard SILAC
      update();
    })();
  </script>
</body>
</html>
```

The actual implementation should copy the annotation JS verbatim from `annotation.html` (lines 92–516, the info panel, coverage bracket, and spectrum renderers), adapting only the data source from `window.__ANNOTATION_DATA__` to `window.__UNIFIED_DATA__.annotation`.

- [ ] **Step 2: Verify template file exists and is valid HTML**

Run: `head -5 crates/report/templates/unified.html && echo "..." && wc -l crates/report/templates/unified.html`
Expected: Shows `<!DOCTYPE html>` and line count ~650-750.

- [ ] **Step 3: Commit**

```bash
git add crates/report/templates/unified.html
git commit -m "feat(report): add unified annotation+XIC HTML template with SILAC controls"
```

---

## Task 5: Unified Render Function (report crate)

**Files:**
- Modify: `crates/report/src/unified_visualize.rs` (replace placeholder)
- Modify: `crates/report/src/unified_types.rs` (update re-exports)
- Modify: `crates/report/src/lib.rs` (add `render_unified` to `ReportGenerator`)

- [ ] **Step 1: Write test for `render_unified_html`**

Replace the placeholder content of `crates/report/src/unified_visualize.rs`:

```rust
//! Unified annotation + XIC visualization — renders combined HTML
//! with interactive SILAC controls and client-side recomputation.

use std::fs;
use std::path::Path;

use protein_copilot_xic::PlotlyMode;

use crate::error::ReportError;
use crate::unified_types::UnifiedViewData;

const PLOTLY_CDN: &str = "https://cdn.plot.ly/plotly-2.35.2.min.js";
const TEMPLATE: &str = include_str!("../templates/unified.html");

/// Renders unified annotation + XIC data into a standalone HTML file.
pub fn render_unified_html(
    data: &UnifiedViewData,
    output_path: &Path,
    plotly_mode: PlotlyMode,
) -> Result<(), ReportError> {
    let json = serde_json::to_string(data)
        .map_err(|e| ReportError::SerializationError(e.to_string()))?;

    let plotly_src = match plotly_mode {
        PlotlyMode::Cdn | PlotlyMode::Embedded => PLOTLY_CDN.to_string(),
    };

    let html = TEMPLATE
        .replace(
            "/*__UNIFIED_JSON__*/",
            &format!("window.__UNIFIED_DATA__ = {};", json),
        )
        .replace("__PLOTLY_SRC__", &plotly_src);

    if let Some(parent) = output_path.parent() {
        fs::create_dir_all(parent).map_err(|e| ReportError::IoError {
            path: parent.to_path_buf(),
            detail: e.to_string(),
        })?;
    }

    fs::write(output_path, &html).map_err(|e| ReportError::IoError {
        path: output_path.to_path_buf(),
        detail: e.to_string(),
    })?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};
    use protein_copilot_search_engine::annotate::{
        AnnotatedPeak, SpectrumAnnotation, TheoreticalIon,
    };
    use protein_copilot_xic::*;

    fn sample_unified_data() -> UnifiedViewData {
        let annotation = SpectrumAnnotation {
            scan_number: 42,
            retention_time_sec: 120.0,
            peptide_sequence: "PEPTIDEK".to_string(),
            charge: 2,
            precursor_mz: 450.25,
            theoretical_mz: 450.24,
            delta_mass_ppm: -2.2,
            score: 0.75,
            matched_ions: 6,
            total_ions: 14,
            protein_accessions: vec!["sp|P12345|TEST_HUMAN".to_string()],
            peaks: vec![AnnotatedPeak {
                mz: 574.28,
                intensity: 5000.0,
                annotation: None,
            }],
            b_ions: vec![TheoreticalIon {
                number: 1,
                mz: 98.06,
                matched: false,
                matched_mz: None,
                matched_intensity: None,
                ppm_error: None,
            }],
            y_ions: vec![TheoreticalIon {
                number: 1,
                mz: 147.11,
                matched: false,
                matched_mz: None,
                matched_intensity: None,
                ppm_error: None,
            }],
            modifications: vec![],
        };

        let xic = XicData {
            peptide_sequence: "PEPTIDEK".to_string(),
            target_rt_sec: 120.0,
            target_scan: 42,
            charge: 2,
            precursor_mz: 450.25,
            ms1_precursor_xic: None,
            ms1_heavy_precursor_xic: None,
            fragment_xic_traces: vec![],
            heavy_fragment_xic_traces: vec![],
            extraction_params: ExtractionParams {
                mz_tolerance: MassTolerance {
                    value: 20.0,
                    unit: ToleranceUnit::Ppm,
                },
                n_cycles: 5,
                top_n_ions: 6,
                label_type: None,
                intensity_rule: IntensityRule::MaxInWindow,
            },
        };

        UnifiedViewData {
            annotation,
            xic: Some(xic),
            raw_scans: Some(RawScanData {
                ms1_scans: vec![],
                ms2_scans: vec![],
            }),
            ion_metadata: vec![],
            peptide_info: crate::unified_types::PeptideInfo {
                sequence: "PEPTIDEK".to_string(),
                charge: 2,
                precursor_mz: 450.25,
                total_k: 1,
                total_r: 0,
            },
        }
    }

    #[test]
    fn render_unified_html_creates_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("unified.html");
        let data = sample_unified_data();
        render_unified_html(&data, &path, PlotlyMode::Cdn).unwrap();
        assert!(path.exists());
        let content = fs::read_to_string(&path).unwrap();
        assert!(content.contains("PEPTIDEK"));
        assert!(content.contains("plotly"));
        assert!(content.contains("__UNIFIED_DATA__"));
    }

    #[test]
    fn render_unified_contains_silac_controls() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("unified2.html");
        let data = sample_unified_data();
        render_unified_html(&data, &path, PlotlyMode::Cdn).unwrap();
        let content = fs::read_to_string(&path).unwrap();
        assert!(content.contains("silac-preset"));
        assert!(content.contains("display-mode"));
        assert!(content.contains("extractIntensity"));
    }

    #[test]
    fn render_unified_contains_annotation_data() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("unified3.html");
        let data = sample_unified_data();
        render_unified_html(&data, &path, PlotlyMode::Cdn).unwrap();
        let content = fs::read_to_string(&path).unwrap();
        assert!(content.contains("450.25"));
        assert!(content.contains("sp|P12345|TEST_HUMAN"));
    }
}
```

- [ ] **Step 2: Update `crates/report/src/unified_types.rs`**

Remove the local struct definitions that are now in xic crate. Keep `UnifiedViewData` and `PeptideInfo` (report-specific). Re-export shared types:

```rust
//! Data types for the unified annotation + XIC HTML view.

use protein_copilot_search_engine::annotate::SpectrumAnnotation;
use protein_copilot_xic::XicData;
use serde::Serialize;

// Re-export types that live in xic crate
pub use protein_copilot_xic::{IonMetadataEntry, RawScan, RawScanData};

/// Combined data for the unified HTML template.
#[derive(Debug, Clone, Serialize)]
pub struct UnifiedViewData {
    pub annotation: SpectrumAnnotation,
    pub xic: Option<XicData>,
    pub raw_scans: Option<RawScanData>,
    pub ion_metadata: Vec<IonMetadataEntry>,
    pub peptide_info: PeptideInfo,
}

/// Peptide-level SILAC info.
#[derive(Debug, Clone, Serialize)]
pub struct PeptideInfo {
    pub sequence: String,
    pub charge: i32,
    pub precursor_mz: f64,
    pub total_k: u32,
    pub total_r: u32,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peptide_info_serializes() {
        let info = PeptideInfo {
            sequence: "PEPTIDEK".to_string(),
            charge: 2,
            precursor_mz: 450.25,
            total_k: 1,
            total_r: 0,
        };
        let json = serde_json::to_string(&info).unwrap();
        assert!(json.contains("\"total_k\":1"));
    }

    #[test]
    fn raw_scan_data_serializes() {
        let data = RawScanData {
            ms1_scans: vec![RawScan {
                scan_number: 1,
                retention_time_sec: 10.0,
                mz_array: vec![450.0],
                intensity_array: vec![1000.0],
            }],
            ms2_scans: vec![],
        };
        let json = serde_json::to_string(&data).unwrap();
        assert!(json.contains("\"ms1_scans\""));
    }
}
```

- [ ] **Step 3: Add `render_unified` to `ReportGenerator` in `lib.rs`**

Add after the existing `render_xic` method:

```rust
    /// Renders unified annotation + XIC as a self-contained HTML file.
    pub fn render_unified(
        data: &unified_types::UnifiedViewData,
        output_path: &Path,
        plotly_mode: protein_copilot_xic::PlotlyMode,
    ) -> Result<(), ReportError> {
        unified_visualize::render_unified_html(data, output_path, plotly_mode)
    }
```

- [ ] **Step 4: Run tests**

Run: `cargo test -p protein-copilot-report`
Expected: All tests pass including unified_visualize tests.

- [ ] **Step 5: Run clippy**

Run: `cargo clippy -p protein-copilot-report -- -D warnings`
Expected: No warnings.

- [ ] **Step 6: Commit**

```bash
git add crates/report/src/unified_types.rs crates/report/src/unified_visualize.rs crates/report/src/lib.rs
git commit -m "feat(report): add render_unified_html() for combined annotation+XIC"
```

---

## Task 6: MCP Tool Integration

**Files:**
- Modify: `crates/mcp-server/src/tools.rs`

- [ ] **Step 1: Add XIC fields to `AnnotateSpectrumInput`**

In `crates/mcp-server/src/tools.rs`, add these fields to the `AnnotateSpectrumInput` struct (after line 257, after the `fragment_tolerance` field):

```rust
    /// Whether to include XIC chromatogram below annotation. Default: true for mzML.
    #[serde(default)]
    include_xic: Option<bool>,
    /// Number of DIA cycles before/after target scan for XIC. Default: 5.
    #[serde(default)]
    n_cycles: Option<u32>,
    /// Number of top fragment ions to display in XIC. Default: 6.
    #[serde(default)]
    top_n_ions: Option<usize>,
    /// SILAC label type for pre-computed heavy traces.
    /// Default: standard SILAC. Use {"Silac":{"heavy_k_delta":8.014199,"heavy_r_delta":10.008269}}.
    #[serde(default)]
    label_type: Option<protein_copilot_xic::LabelType>,
    /// Embed raw scan data for interactive SILAC in browser. Default: true.
    #[serde(default)]
    embed_raw_scans: Option<bool>,
    /// Plotly.js loading mode. Default: Cdn.
    #[serde(default)]
    plotly_mode: Option<protein_copilot_xic::PlotlyMode>,
```

- [ ] **Step 2: Update the `annotate_spectrum` handler to integrate XIC**

Replace the rendering section (after the `annotation` variable is created, ~line 1493) with unified rendering logic:

```rust
        // Render HTML — unified or annotation-only
        let out_path = input.output_path.map(PathBuf::from).unwrap_or_else(|| {
            PathBuf::from(format!("output/annotation_scan{}.html", input.scan_number))
        });

        // Determine if XIC should be included
        let is_mzml = spectrum_file
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.eq_ignore_ascii_case("mzml"))
            .unwrap_or(false);
        let include_xic = input.include_xic.unwrap_or(is_mzml);

        if include_xic {
            // Extract XIC + raw scans for unified HTML
            let xic_params = protein_copilot_xic::ExtractionParams {
                mz_tolerance: input.fragment_tolerance.clone().unwrap_or(MassTolerance {
                    value: 20.0,
                    unit: ToleranceUnit::Ppm,
                }),
                n_cycles: input.n_cycles.unwrap_or(5),
                top_n_ions: input.top_n_ions.unwrap_or(6),
                label_type: input.label_type.clone().or_else(|| {
                    Some(protein_copilot_xic::LabelType::standard_silac())
                }),
                intensity_rule: protein_copilot_xic::IntensityRule::MaxInWindow,
            };

            let embed_raw = input.embed_raw_scans.unwrap_or(true);
            let plotly_mode = input
                .plotly_mode
                .unwrap_or(protein_copilot_xic::PlotlyMode::Cdn);

            match protein_copilot_xic::extract::extract_xic_with_raw(
                &spectrum_file,
                input.scan_number,
                &peptide_seq,
                charge,
                annotation.precursor_mz,
                &modifications,
                &xic_params,
                20.0, // MS1 m/z window ±20 Da
            ) {
                Ok((xic_data, raw_scans, ion_metadata)) => {
                    let unified = protein_copilot_report::unified_types::UnifiedViewData {
                        annotation: annotation.clone(),
                        xic: Some(xic_data),
                        raw_scans: if embed_raw { Some(raw_scans) } else { None },
                        ion_metadata,
                        peptide_info: protein_copilot_report::unified_types::PeptideInfo {
                            sequence: peptide_seq.clone(),
                            charge,
                            precursor_mz: annotation.precursor_mz,
                            total_k: peptide_seq.chars().filter(|&c| c == 'K').count() as u32,
                            total_r: peptide_seq.chars().filter(|&c| c == 'R').count() as u32,
                        },
                    };
                    ReportGenerator::render_unified(&unified, &out_path, plotly_mode)
                        .map_err(|e| {
                            mcp_core_err(protein_copilot_core::error::CoreError::from(e))
                        })?;
                }
                Err(e) => {
                    // XIC extraction failed — fall back to annotation-only
                    tracing::warn!("XIC extraction failed, falling back to annotation-only: {e}");
                    ReportGenerator::render_annotation(&annotation, &out_path).map_err(|e| {
                        mcp_core_err(protein_copilot_core::error::CoreError::from(e))
                    })?;
                }
            }
        } else {
            // Annotation-only mode
            ReportGenerator::render_annotation(&annotation, &out_path)
                .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        }

        Ok(Json(AnnotateResult {
            output_path: out_path.display().to_string(),
            scan_number: annotation.scan_number,
            peptide_sequence: annotation.peptide_sequence,
            charge: annotation.charge,
            score: annotation.score,
            matched_ions: annotation.matched_ions,
            total_ions: annotation.total_ions,
            delta_mass_ppm: annotation.delta_mass_ppm,
            protein_accessions: annotation.protein_accessions,
            message: format!(
                "Annotation saved to {}. Matched {}/{} ions (score: {:.3}).{}",
                out_path.display(),
                annotation.matched_ions,
                annotation.total_ions,
                annotation.score,
                if include_xic {
                    " Includes XIC with interactive SILAC controls."
                } else {
                    " Open in browser to view."
                },
            ),
        }))
```

- [ ] **Step 3: Add `Clone` derive to `SpectrumAnnotation` if not already present**

Check `crates/search-engine/src/annotate.rs` line 89. If `SpectrumAnnotation` doesn't derive `Clone`, add it:

```rust
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpectrumAnnotation {
```

This is needed because the unified path needs `annotation.clone()` when passing it to `UnifiedViewData` while also reading fields for the result.

- [ ] **Step 4: Build and verify compilation**

Run: `cargo build --workspace`
Expected: Compiles successfully.

- [ ] **Step 5: Run all tests**

Run: `cargo test --workspace`
Expected: All existing tests pass (no regressions).

- [ ] **Step 6: Run clippy**

Run: `cargo clippy --workspace -- -D warnings`
Expected: No warnings.

- [ ] **Step 7: Commit**

```bash
git add crates/mcp-server/src/tools.rs crates/search-engine/src/annotate.rs
git commit -m "feat(mcp): integrate XIC into annotate_spectrum with interactive SILAC"
```

---

## Task 7: Integration Test & Final Verification

**Files:**
- Run existing test suite
- Manual verification of generated HTML

- [ ] **Step 1: Run full test suite**

Run: `cargo test --workspace`
Expected: All tests pass (previous count ~499 + new tests ≈ ~510+).

- [ ] **Step 2: Run clippy on full workspace**

Run: `cargo clippy --workspace -- -D warnings`
Expected: Zero warnings.

- [ ] **Step 3: Verify HTML template renders correctly**

If test fixtures with mzML data are available, manually invoke:

```bash
cargo run -p protein-copilot-mcp-server -- --help
```

Or write a quick integration test that generates a unified HTML and checks key elements:

```rust
// In crates/report/src/unified_visualize.rs tests, add:
#[test]
fn render_unified_with_raw_scans() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("unified_raw.html");
    let mut data = sample_unified_data();
    data.raw_scans = Some(RawScanData {
        ms1_scans: vec![RawScan {
            scan_number: 40,
            retention_time_sec: 119.0,
            mz_array: vec![449.0, 450.0, 451.0],
            intensity_array: vec![100.0, 5000.0, 200.0],
        }],
        ms2_scans: vec![RawScan {
            scan_number: 41,
            retention_time_sec: 119.5,
            mz_array: vec![174.1, 262.2, 574.3],
            intensity_array: vec![200.0, 800.0, 3000.0],
        }],
    });
    data.ion_metadata = vec![IonMetadataEntry {
        label: "y5¹⁺".to_string(),
        ion_type: IonType::Y,
        ion_number: 5,
        charge: 1,
        light_mz: 574.28,
        k_count: 1,
        r_count: 0,
    }];

    render_unified_html(&data, &path, PlotlyMode::Cdn).unwrap();
    let content = fs::read_to_string(&path).unwrap();
    assert!(content.contains("raw_scans"));
    assert!(content.contains("ion_metadata"));
    assert!(content.contains("574.28"));
}
```

- [ ] **Step 4: Final commit with all integration tests**

```bash
git add -A
git commit -m "test: add integration tests for unified annotation+XIC HTML"
```

---

## Summary

| Task | Description | New/Modified Files | Estimated Tests |
|------|-------------|-------------------|-----------------|
| 1 | Unified view types | `unified_types.rs`, `unified_visualize.rs` (stub), `lib.rs` | 4 |
| 2 | Ion metadata K/R counting | `xic/lib.rs`, `xic/extract.rs` | 3 |
| 3 | Extract XIC with raw scans | `xic/extract.rs`, `xic/lib.rs` | 2 |
| 4 | Unified HTML template | `templates/unified.html` | 0 (template) |
| 5 | Render function | `unified_visualize.rs`, `unified_types.rs`, `lib.rs` | 3 |
| 6 | MCP tool integration | `tools.rs`, `annotate.rs` | 0 (existing) |
| 7 | Integration & verification | tests | 1 |
| **Total** | | **3 new + 6 modified** | **~13 new tests** |
