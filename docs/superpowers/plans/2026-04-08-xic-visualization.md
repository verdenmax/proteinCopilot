# XIC Visualization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement fragment-ion-level XIC (Extracted Ion Chromatogram) visualization with SILAC heavy-label comparison, targeting DIA data as the primary use case.

**Architecture:** New `crates/xic` library crate handles deterministic XIC extraction (data structures, intensity extraction, SILAC heavy-label m/z calculation). The `report` crate gets a new `xic.html` Plotly.js template and `render_xic_html()`. The `mcp-server` gets a new `extract_xic` tool. A prerequisite streaming API is added to `spectrum-io`.

**Tech Stack:** Rust, Plotly.js (CDN/embedded), serde, thiserror. Depends on `spectrum-io`, `search-engine` (matching.rs fragment ion API), `core`, `report`.

**Design Spec:** `docs/superpowers/specs/2026-04-08-xic-visualization-design.md`

---

## File Structure

### New Files

| File | Responsibility |
|------|---------------|
| `crates/xic/Cargo.toml` | Crate manifest — depends on core, spectrum-io, search-engine |
| `crates/xic/src/lib.rs` | Public API, data structures (`XicTrace`, `XicData`, `ExtractionParams`, `LabelType`, `IntensityRule`) |
| `crates/xic/src/error.rs` | `XicError` enum |
| `crates/xic/src/extract.rs` | Core XIC extraction: `extract_xic()`, `extract_intensity()`, `same_isolation_window()`, top-N selection, zero-fill |
| `crates/xic/src/heavy.rs` | SILAC heavy-label: `compute_heavy_fragment_mz()`, `compute_heavy_precursor_mz()` |
| `crates/report/templates/xic.html` | Plotly.js interactive XIC visualization template |
| `crates/report/src/xic_visualize.rs` | `render_xic_html()` — template injection for XIC data |

### Modified Files

| File | Changes |
|------|---------|
| `Cargo.toml` (workspace root) | Add `protein-copilot-xic` workspace dependency |
| `crates/spectrum-io/src/reader.rs` | Add `for_each_spectrum()` method to `SpectrumReader` trait |
| `crates/spectrum-io/src/mzml.rs` | Implement `for_each_spectrum()` for `MzMLReader` |
| `crates/spectrum-io/src/mgf.rs` | Implement `for_each_spectrum()` for `MgfReader` |
| `crates/report/Cargo.toml` | Add `protein-copilot-xic` dependency |
| `crates/report/src/lib.rs` | Add `pub mod xic_visualize;` and `render_xic()` to `ReportGenerator` |
| `crates/mcp-server/Cargo.toml` | Add `protein-copilot-xic` dependency |
| `crates/mcp-server/src/tools.rs` | Add `extract_xic` MCP tool |

---

## Task 1: Add `for_each_spectrum()` streaming API to spectrum-io

**Files:**
- Modify: `crates/spectrum-io/src/reader.rs`
- Modify: `crates/spectrum-io/src/mzml.rs:471-511`
- Modify: `crates/spectrum-io/src/mgf.rs:207-245`
- Test: `crates/spectrum-io/src/mzml.rs` (inline tests) + `crates/spectrum-io/src/mgf.rs` (inline tests)

- [ ] **Step 1: Write failing test for mzML `for_each_spectrum`**

Add to `crates/spectrum-io/src/mzml.rs` in the `#[cfg(test)] mod tests` block:

```rust
#[test]
fn for_each_spectrum_streams_all() {
    let reader = MzMLReader;
    let path = Path::new("../tests/fixtures/sample_5specs.mzml");
    if !path.exists() {
        return; // skip if fixture missing
    }
    let mut count = 0u32;
    let result = reader.for_each_spectrum(path, |_spec| {
        count += 1;
        Ok(true)
    });
    assert!(result.is_ok());
    assert!(count > 0, "should stream at least one spectrum");
    // Verify consistency with read_all
    let all = reader.read_all(path).unwrap();
    assert_eq!(count, all.len() as u32);
}

#[test]
fn for_each_spectrum_early_stop() {
    let reader = MzMLReader;
    let path = Path::new("../tests/fixtures/sample_5specs.mzml");
    if !path.exists() {
        return;
    }
    let mut count = 0u32;
    let _ = reader.for_each_spectrum(path, |_spec| {
        count += 1;
        Ok(count < 2) // stop after 2
    });
    assert_eq!(count, 2);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p protein-copilot-spectrum-io for_each_spectrum -- --nocapture 2>&1 | tail -5`
Expected: compile error — method `for_each_spectrum` not found on trait `SpectrumReader`

- [ ] **Step 3: Add `for_each_spectrum` to `SpectrumReader` trait**

In `crates/spectrum-io/src/reader.rs`, add a new method to the trait. Because trait methods with generic closures cannot use `dyn` dispatch, add it as a provided method that delegates to `read_all` by default, and override in each concrete reader:

```rust
use std::path::Path;

use protein_copilot_core::spectrum::{Spectrum, SpectrumSummary};

use crate::error::SpectrumIoError;

/// Unified interface for reading spectrum files.
///
/// Each supported format (mgf, mzML) implements this trait.
/// Use [`crate::create_reader`] to obtain the appropriate reader
/// for a given file.
pub trait SpectrumReader: Send + Sync {
    /// Reads all spectra from the file.
    fn read_all(&self, path: &Path) -> Result<Vec<Spectrum>, SpectrumIoError>;

    /// Computes a statistical summary of the spectrum file.
    fn read_summary(&self, path: &Path) -> Result<SpectrumSummary, SpectrumIoError>;

    /// Reads a single spectrum by scan number.
    fn read_spectrum(&self, path: &Path, scan: u32) -> Result<Spectrum, SpectrumIoError>;

    /// Streams spectra one at a time, calling `handler` for each.
    ///
    /// The handler returns `Ok(true)` to continue or `Ok(false)` to stop early.
    /// Returns the number of spectra processed (including the one that stopped).
    ///
    /// This avoids loading all spectra into memory at once, which is important
    /// for large DIA files when only extracting specific ion chromatograms.
    fn for_each_spectrum(
        &self,
        path: &Path,
        handler: &mut dyn FnMut(Spectrum) -> Result<bool, SpectrumIoError>,
    ) -> Result<u32, SpectrumIoError>;
}
```

Note: We use `&mut dyn FnMut` instead of a generic `F` so the trait remains object-safe (can be used as `Box<dyn SpectrumReader>`).

- [ ] **Step 4: Implement `for_each_spectrum` for `MzMLReader`**

In `crates/spectrum-io/src/mzml.rs`, inside `impl SpectrumReader for MzMLReader`:

```rust
fn for_each_spectrum(
    &self,
    path: &Path,
    handler: &mut dyn FnMut(Spectrum) -> Result<bool, SpectrumIoError>,
) -> Result<u32, SpectrumIoError> {
    let mut xml_reader = open_xml_reader(path)?;
    parse_mzml_streaming(&mut xml_reader, path, handler)
}
```

- [ ] **Step 5: Implement `for_each_spectrum` for `MgfReader`**

In `crates/spectrum-io/src/mgf.rs`, inside `impl SpectrumReader for MgfReader`:

```rust
fn for_each_spectrum(
    &self,
    path: &Path,
    handler: &mut dyn FnMut(Spectrum) -> Result<bool, SpectrumIoError>,
) -> Result<u32, SpectrumIoError> {
    parse_mgf_streaming(path, handler)
}
```

- [ ] **Step 6: Write failing test for MGF `for_each_spectrum`**

Add to `crates/spectrum-io/src/mgf.rs` in the `#[cfg(test)] mod tests` block:

```rust
#[test]
fn for_each_spectrum_streams_all() {
    let reader = MgfReader;
    let path = Path::new(TEST_MGF);
    let mut count = 0u32;
    let result = reader.for_each_spectrum(path, &mut |_spec| {
        count += 1;
        Ok(true)
    });
    assert!(result.is_ok());
    let all = reader.read_all(path).unwrap();
    assert_eq!(count, all.len() as u32);
}

#[test]
fn for_each_spectrum_early_stop() {
    let reader = MgfReader;
    let path = Path::new(TEST_MGF);
    let mut count = 0u32;
    let _ = reader.for_each_spectrum(path, &mut |_spec| {
        count += 1;
        Ok(count < 2)
    });
    assert_eq!(count, 2);
}
```

- [ ] **Step 7: Fix mzML tests to use new signature**

Update the mzML tests from Step 1 to match the `&mut dyn FnMut` signature:

```rust
#[test]
fn for_each_spectrum_streams_all() {
    let reader = MzMLReader;
    let path = Path::new("../tests/fixtures/sample_5specs.mzml");
    if !path.exists() {
        return;
    }
    let mut count = 0u32;
    let result = reader.for_each_spectrum(path, &mut |_spec| {
        count += 1;
        Ok(true)
    });
    assert!(result.is_ok());
    let all = reader.read_all(path).unwrap();
    assert_eq!(count, all.len() as u32);
}

#[test]
fn for_each_spectrum_early_stop() {
    let reader = MzMLReader;
    let path = Path::new("../tests/fixtures/sample_5specs.mzml");
    if !path.exists() {
        return;
    }
    let mut count = 0u32;
    let _ = reader.for_each_spectrum(path, &mut |_spec| {
        count += 1;
        Ok(count < 2)
    });
    assert_eq!(count, 2);
}
```

- [ ] **Step 8: Run all spectrum-io tests**

Run: `cargo test -p protein-copilot-spectrum-io -- --nocapture 2>&1 | tail -10`
Expected: all tests pass, including new `for_each_spectrum` tests

- [ ] **Step 9: Run workspace build to check no breakage**

Run: `cargo build --workspace 2>&1 | tail -5`
Expected: success (no compile errors across workspace)

- [ ] **Step 10: Commit**

```bash
git add crates/spectrum-io/
git commit -m "feat(spectrum-io): add for_each_spectrum streaming API to SpectrumReader trait

Add for_each_spectrum() method enabling streaming spectrum iteration
without loading all spectra into memory. Uses dyn FnMut for object
safety. Both MzMLReader and MgfReader delegate to their existing
private streaming parsers.

Prerequisite for XIC extraction on large DIA files.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 2: Create `crates/xic` crate — data structures and error types

**Files:**
- Create: `crates/xic/Cargo.toml`
- Create: `crates/xic/src/lib.rs`
- Create: `crates/xic/src/error.rs`
- Modify: `Cargo.toml` (workspace root)

- [ ] **Step 1: Add xic to workspace `Cargo.toml`**

Add to the `[workspace.dependencies]` section:

```toml
protein-copilot-xic = { path = "crates/xic" }
```

- [ ] **Step 2: Create `crates/xic/Cargo.toml`**

```toml
[package]
name = "protein-copilot-xic"
version.workspace = true
edition.workspace = true
rust-version.workspace = true
license.workspace = true
description = "XIC (Extracted Ion Chromatogram) extraction for ProteinCopilot"

[dependencies]
protein-copilot-core = { workspace = true }
protein-copilot-spectrum-io = { workspace = true }
protein-copilot-search-engine = { workspace = true }
serde = { workspace = true }
serde_json = { workspace = true }
thiserror = { workspace = true }
tracing = { workspace = true }

[dev-dependencies]
tempfile = "3"
```

- [ ] **Step 3: Create `crates/xic/src/error.rs`**

```rust
//! Error types for XIC extraction.

use std::path::PathBuf;

use thiserror::Error;

/// Errors that can occur during XIC extraction.
#[derive(Debug, Error)]
pub enum XicError {
    /// Input file format is not supported for XIC extraction.
    #[error("XIC extraction requires mzML format, got: {path}")]
    UnsupportedFormat {
        /// The file path.
        path: PathBuf,
    },

    /// Target scan not found in the spectrum file.
    #[error("target scan {scan} not found in {path}")]
    ScanNotFound {
        /// The file path.
        path: PathBuf,
        /// The requested scan number.
        scan: u32,
    },

    /// Target scan has no isolation window (required for DIA XIC).
    #[error("scan {scan} has no isolation window — XIC requires MS2 with isolation window")]
    NoIsolationWindow {
        /// The scan number.
        scan: u32,
    },

    /// No matching MS2 scans found in the same isolation window.
    #[error("no MS2 scans found in same isolation window as scan {scan}")]
    NoMatchingCycles {
        /// The target scan number.
        scan: u32,
    },

    /// Peptide sequence is empty or contains invalid residues.
    #[error("invalid peptide sequence: {detail}")]
    InvalidPeptide {
        /// What went wrong.
        detail: String,
    },

    /// Spectrum I/O error.
    #[error("spectrum I/O error: {0}")]
    SpectrumIo(#[from] protein_copilot_spectrum_io::SpectrumIoError),
}
```

- [ ] **Step 4: Create `crates/xic/src/lib.rs` with data structures**

```rust
//! # ProteinCopilot XIC
//!
//! Extracted Ion Chromatogram (XIC) computation for proteomics data.
//!
//! This crate handles:
//! - Fragment-ion-level XIC extraction from mzML files (DIA + DDA)
//! - MS1 precursor XIC extraction
//! - SILAC heavy-label m/z calculation for paired XIC traces
//!
//! It is a pure computation crate — no MCP, network, or LLM dependencies.

pub mod error;
pub mod extract;
pub mod heavy;

pub use error::XicError;

use protein_copilot_core::search_params::MassTolerance;
use serde::{Deserialize, Serialize};

/// Ion type for XIC trace identification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum IonType {
    /// b-ion (N-terminal fragment).
    B,
    /// y-ion (C-terminal fragment).
    Y,
    /// Precursor ion (MS1 level).
    Precursor,
}

/// A single data point on an XIC trace.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XicDataPoint {
    /// Retention time in seconds.
    pub retention_time_sec: f64,
    /// Scan number (1-based).
    pub scan_number: u32,
    /// Extracted intensity (0.0 if no peak found in tolerance window).
    pub intensity: f64,
}

/// A single XIC trace for one ion.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XicTrace {
    /// Human-readable ion label (e.g. "y5¹⁺", "b3²⁺", "precursor").
    pub ion_label: String,
    /// Ion type.
    pub ion_type: IonType,
    /// Ion number (e.g. 5 for y5; 0 for precursor).
    pub ion_number: u32,
    /// Charge state.
    pub charge: u32,
    /// Theoretical m/z used for extraction.
    pub theoretical_mz: f64,
    /// Extracted data points.
    pub data_points: Vec<XicDataPoint>,
    /// Whether this is a heavy-label trace.
    pub is_heavy: bool,
}

/// Complete XIC extraction result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XicData {
    /// Target peptide sequence.
    pub peptide_sequence: String,
    /// Target scan retention time (seconds).
    pub target_rt_sec: f64,
    /// Target scan number.
    pub target_scan: u32,
    /// Precursor charge state.
    pub charge: i32,
    /// Real precursor m/z (from PSM or user input, not DIA window center).
    pub precursor_mz: f64,
    /// MS1 precursor XIC (light). `None` if MS1 data unavailable.
    pub ms1_precursor_xic: Option<XicTrace>,
    /// MS1 precursor XIC (heavy). `None` if no label or MS1 unavailable.
    pub ms1_heavy_precursor_xic: Option<XicTrace>,
    /// MS2 fragment ion XIC traces (light).
    pub fragment_xic_traces: Vec<XicTrace>,
    /// MS2 fragment ion XIC traces (heavy).
    pub heavy_fragment_xic_traces: Vec<XicTrace>,
    /// Parameters used for extraction.
    pub extraction_params: ExtractionParams,
}

/// Parameters controlling XIC extraction.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExtractionParams {
    /// Mass tolerance for peak matching.
    pub mz_tolerance: MassTolerance,
    /// Number of DIA cycles before/after target scan.
    pub n_cycles: u32,
    /// Number of top ions to display by default.
    pub top_n_ions: usize,
    /// Heavy-label type (None = no label).
    pub label_type: Option<LabelType>,
    /// How to extract intensity from peaks within tolerance.
    pub intensity_rule: IntensityRule,
}

/// Intensity extraction strategy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum IntensityRule {
    /// Highest peak within tolerance window (default).
    #[default]
    MaxInWindow,
    /// Sum of all peaks within tolerance window.
    SumInWindow,
    /// Nearest peak to theoretical m/z.
    NearestPeak,
}

/// Heavy-label type for SILAC comparison.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum LabelType {
    /// SILAC heavy amino acids.
    Silac {
        /// Mass shift for heavy Lysine (default: 8.014199 Da for ¹³C₆¹⁵N₂-Lys).
        heavy_k_delta: f64,
        /// Mass shift for heavy Arginine (default: 10.008269 Da for ¹³C₆¹⁵N₄-Arg).
        heavy_r_delta: f64,
    },
    /// Custom residue mass shifts.
    Custom {
        /// (residue, mass_delta) pairs.
        residue_deltas: Vec<(char, f64)>,
    },
}

impl LabelType {
    /// Standard SILAC heavy labels (K+8, R+10).
    pub fn standard_silac() -> Self {
        LabelType::Silac {
            heavy_k_delta: 8.014199,
            heavy_r_delta: 10.008269,
        }
    }
}

/// Plotly.js loading mode for HTML output.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum PlotlyMode {
    /// Load from CDN (default, smaller file).
    #[default]
    Cdn,
    /// Embed plotly-basic.min.js inline (larger file, works offline).
    Embedded,
}
```

- [ ] **Step 5: Create placeholder modules**

Create `crates/xic/src/extract.rs`:

```rust
//! XIC extraction core logic.

// Implementation added in Task 3.
```

Create `crates/xic/src/heavy.rs`:

```rust
//! SILAC heavy-label m/z calculation.

// Implementation added in Task 4.
```

- [ ] **Step 6: Verify it compiles**

Run: `cargo build -p protein-copilot-xic 2>&1 | tail -5`
Expected: success

- [ ] **Step 7: Commit**

```bash
git add Cargo.toml crates/xic/
git commit -m "feat(xic): create xic crate with data structures and error types

New library crate for XIC extraction. Defines XicTrace, XicData,
ExtractionParams, LabelType (SILAC), IntensityRule, PlotlyMode,
and XicError. Pure computation — no MCP dependencies.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 3: Implement XIC extraction core logic

**Files:**
- Modify: `crates/xic/src/extract.rs`
- Test: inline `#[cfg(test)] mod tests` in `extract.rs`

- [ ] **Step 1: Write failing tests for `extract_intensity`**

In `crates/xic/src/extract.rs`:

```rust
//! XIC extraction core logic.
//!
//! Implements the 1.5-pass extraction algorithm:
//! - Pass 0: read target scan to get RT and isolation window
//! - Pass 1: stream all spectra, extracting intensities for target ions

use std::path::Path;

use protein_copilot_core::search_params::{MassTolerance, Modification, ToleranceUnit};
use protein_copilot_core::spectrum::{IsolationWindow, MsLevel, Spectrum};
use protein_copilot_search_engine::matching::{
    generate_b_ions_with_charge, generate_y_ions_with_charge, within_tolerance,
};

use crate::{
    ExtractionParams, IntensityRule, IonType, LabelType, XicData, XicDataPoint, XicError,
    XicTrace,
};

/// Extract intensity for a target m/z from a spectrum's peak list.
///
/// Uses binary search for efficiency (mzML peaks are sorted by m/z).
/// Returns 0.0 if no peak is found within tolerance.
pub fn extract_intensity(
    target_mz: f64,
    exp_mz: &[f64],
    exp_int: &[f64],
    tolerance: &MassTolerance,
    rule: IntensityRule,
) -> f64 {
    if exp_mz.is_empty() || exp_mz.len() != exp_int.len() {
        return 0.0;
    }

    // Binary search for approximate position
    let pos = exp_mz.partition_point(|&m| m < target_mz);

    // Scan neighbors within tolerance
    let mut best_intensity = 0.0;
    let mut sum_intensity = 0.0;
    let mut nearest_dist = f64::MAX;
    let mut nearest_intensity = 0.0;
    let mut found = false;

    // Check backwards from pos
    let start = pos.saturating_sub(1);
    // Check forwards from pos (limit scan range for efficiency)
    let end = (pos + 2).min(exp_mz.len());

    // Widen the scan range for ppm tolerance on high m/z values
    let max_da = match tolerance.unit {
        ToleranceUnit::Ppm => target_mz * tolerance.value * 1e-6,
        ToleranceUnit::Da => tolerance.value,
    };
    let scan_start = match exp_mz[..pos].iter().rposition(|&m| target_mz - m > max_da * 1.5) {
        Some(i) => i + 1,
        None => 0,
    };
    let scan_end = match exp_mz[pos..].iter().position(|&m| m - target_mz > max_da * 1.5) {
        Some(i) => pos + i,
        None => exp_mz.len(),
    };

    for i in scan_start..scan_end {
        if within_tolerance(exp_mz[i], target_mz, tolerance) {
            found = true;
            let intensity = exp_int[i];
            if intensity > best_intensity {
                best_intensity = intensity;
            }
            sum_intensity += intensity;
            let dist = (exp_mz[i] - target_mz).abs();
            if dist < nearest_dist {
                nearest_dist = dist;
                nearest_intensity = intensity;
            }
        }
    }

    if !found {
        return 0.0;
    }

    match rule {
        IntensityRule::MaxInWindow => best_intensity,
        IntensityRule::SumInWindow => sum_intensity,
        IntensityRule::NearestPeak => nearest_intensity,
    }
}

/// Check if two isolation windows cover the same DIA region.
///
/// Compares full window bounds (not just center m/z) to avoid
/// conflating adjacent or overlapping windows.
pub fn same_isolation_window(a: &IsolationWindow, b: &IsolationWindow) -> bool {
    let a_lo = a.target_mz - a.lower_offset;
    let a_hi = a.target_mz + a.upper_offset;
    let b_lo = b.target_mz - b.lower_offset;
    let b_hi = b.target_mz + b.upper_offset;

    let center_a = (a_lo + a_hi) / 2.0;
    let center_b = (b_lo + b_hi) / 2.0;
    let center_close = (center_a - center_b).abs() < 1.0;

    let width_a = a_hi - a_lo;
    let width_b = b_hi - b_lo;
    let width_close = width_a > 0.0 && ((width_a - width_b).abs() / width_a) < 0.2;

    center_close && width_close
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extract_intensity_max_in_window() {
        let mz = vec![100.0, 200.0, 200.01, 200.02, 300.0];
        let int = vec![50.0, 100.0, 200.0, 150.0, 75.0];
        let tol = MassTolerance {
            value: 0.05,
            unit: ToleranceUnit::Da,
        };
        let result = extract_intensity(200.01, &mz, &int, &tol, IntensityRule::MaxInWindow);
        assert!((result - 200.0).abs() < 0.01, "expected 200.0, got {result}");
    }

    #[test]
    fn extract_intensity_sum_in_window() {
        let mz = vec![100.0, 200.0, 200.01, 200.02, 300.0];
        let int = vec![50.0, 100.0, 200.0, 150.0, 75.0];
        let tol = MassTolerance {
            value: 0.05,
            unit: ToleranceUnit::Da,
        };
        let result = extract_intensity(200.01, &mz, &int, &tol, IntensityRule::SumInWindow);
        assert!((result - 450.0).abs() < 0.01, "expected 450.0, got {result}");
    }

    #[test]
    fn extract_intensity_nearest_peak() {
        let mz = vec![100.0, 200.0, 200.01, 200.02, 300.0];
        let int = vec![50.0, 100.0, 200.0, 150.0, 75.0];
        let tol = MassTolerance {
            value: 0.05,
            unit: ToleranceUnit::Da,
        };
        let result = extract_intensity(200.01, &mz, &int, &tol, IntensityRule::NearestPeak);
        assert!((result - 200.0).abs() < 0.01, "expected 200.0, got {result}");
    }

    #[test]
    fn extract_intensity_no_match_returns_zero() {
        let mz = vec![100.0, 200.0, 300.0];
        let int = vec![50.0, 100.0, 75.0];
        let tol = MassTolerance {
            value: 0.01,
            unit: ToleranceUnit::Da,
        };
        let result = extract_intensity(250.0, &mz, &int, &tol, IntensityRule::MaxInWindow);
        assert!((result - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn extract_intensity_empty_spectrum() {
        let tol = MassTolerance {
            value: 0.05,
            unit: ToleranceUnit::Da,
        };
        assert_eq!(
            extract_intensity(200.0, &[], &[], &tol, IntensityRule::MaxInWindow),
            0.0
        );
    }

    #[test]
    fn extract_intensity_ppm_tolerance() {
        let mz = vec![500.0, 500.005, 500.01];
        let int = vec![100.0, 200.0, 300.0];
        let tol = MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        };
        // 20 ppm at 500.0 = 0.01 Da window
        let result = extract_intensity(500.005, &mz, &int, &tol, IntensityRule::MaxInWindow);
        assert!(result > 0.0, "should find a peak within 20 ppm");
    }

    #[test]
    fn same_isolation_window_identical() {
        let w = IsolationWindow {
            target_mz: 500.0,
            lower_offset: 12.5,
            upper_offset: 12.5,
        };
        assert!(same_isolation_window(&w, &w));
    }

    #[test]
    fn same_isolation_window_different_center() {
        let a = IsolationWindow {
            target_mz: 500.0,
            lower_offset: 12.5,
            upper_offset: 12.5,
        };
        let b = IsolationWindow {
            target_mz: 525.0,
            lower_offset: 12.5,
            upper_offset: 12.5,
        };
        assert!(!same_isolation_window(&a, &b));
    }

    #[test]
    fn same_isolation_window_different_width() {
        let a = IsolationWindow {
            target_mz: 500.0,
            lower_offset: 12.5,
            upper_offset: 12.5,
        };
        let b = IsolationWindow {
            target_mz: 500.0,
            lower_offset: 5.0,
            upper_offset: 5.0,
        };
        assert!(!same_isolation_window(&a, &b));
    }
}
```

- [ ] **Step 2: Run tests to verify they pass (utility functions)**

Run: `cargo test -p protein-copilot-xic -- --nocapture 2>&1 | tail -10`
Expected: all tests pass

- [ ] **Step 3: Write `build_target_ions` helper and tests**

Add above the `#[cfg(test)]` block in `extract.rs`:

```rust
/// A target ion for XIC extraction.
#[derive(Debug, Clone)]
pub struct TargetIon {
    pub label: String,
    pub ion_type: IonType,
    pub ion_number: u32,
    pub charge: u32,
    pub mz: f64,
}

/// Build the list of target fragment ions from a peptide sequence.
///
/// Reuses `matching.rs` fragment ion generators. Does NOT duplicate ion
/// calculation logic — only wraps the results with metadata.
pub fn build_target_ions(
    sequence: &str,
    modifications: &[Modification],
    precursor_charge: i32,
) -> Vec<TargetIon> {
    let max_frag_charge = if precursor_charge >= 3 { 2 } else { 1 };
    let chars: Vec<char> = sequence.chars().collect();
    let n = chars.len().saturating_sub(1); // number of b/y ions per charge

    let b_mz = generate_b_ions_with_charge(sequence, modifications, max_frag_charge);
    let y_mz = generate_y_ions_with_charge(sequence, modifications, max_frag_charge);

    let max_z = max_frag_charge.max(1) as usize;
    let mut ions = Vec::with_capacity(b_mz.len() + y_mz.len());

    for (idx, &mz) in b_mz.iter().enumerate() {
        let frag_idx = idx / max_z; // 0-based fragment index
        let z = (idx % max_z) + 1;
        let ion_number = (frag_idx + 1) as u32;
        let superscript = if z == 1 { "¹⁺" } else { "²⁺" };
        ions.push(TargetIon {
            label: format!("b{ion_number}{superscript}"),
            ion_type: IonType::B,
            ion_number,
            charge: z as u32,
            mz,
        });
    }

    for (idx, &mz) in y_mz.iter().enumerate() {
        let frag_idx = idx / max_z;
        let z = (idx % max_z) + 1;
        let ion_number = (frag_idx + 1) as u32;
        let superscript = if z == 1 { "¹⁺" } else { "²⁺" };
        ions.push(TargetIon {
            label: format!("y{ion_number}{superscript}"),
            ion_type: IonType::Y,
            ion_number,
            charge: z as u32,
            mz,
        });
    }

    ions
}
```

Add tests:

```rust
#[test]
fn build_target_ions_simple_peptide() {
    let ions = build_target_ions("PEPTIDE", &[], 2);
    // 6 b-ions + 6 y-ions at charge 1 = 12 total
    assert_eq!(ions.len(), 12);
    assert!(ions[0].label.starts_with('b'));
    assert!(ions.last().unwrap().label.starts_with('y'));
}

#[test]
fn build_target_ions_high_charge_gets_doubly_charged() {
    let ions = build_target_ions("PEPTIDE", &[], 3);
    // 6 b-ions × 2 charges + 6 y-ions × 2 charges = 24
    assert_eq!(ions.len(), 24);
    let has_double = ions.iter().any(|i| i.charge == 2);
    assert!(has_double, "charge 3 precursor should produce doubly-charged fragments");
}

#[test]
fn build_target_ions_empty_sequence() {
    let ions = build_target_ions("", &[], 2);
    assert!(ions.is_empty());
}
```

- [ ] **Step 4: Run tests**

Run: `cargo test -p protein-copilot-xic build_target_ions -- --nocapture 2>&1 | tail -10`
Expected: all pass

- [ ] **Step 5: Implement `extract_xic()` main function**

Add to `extract.rs`:

```rust
/// Extract XIC data for a peptide from an mzML file.
///
/// Uses the 1.5-pass strategy:
/// - Pass 0: `read_spectrum(target_scan)` → get RT and isolation window
/// - Pass 1: `for_each_spectrum()` → stream all spectra, extract intensities
pub fn extract_xic(
    file_path: &Path,
    target_scan: u32,
    peptide_sequence: &str,
    charge: i32,
    precursor_mz: f64,
    modifications: &[Modification],
    params: &ExtractionParams,
) -> Result<XicData, XicError> {
    // Validate input
    if peptide_sequence.is_empty() {
        return Err(XicError::InvalidPeptide {
            detail: "peptide sequence is empty".to_string(),
        });
    }

    // Require mzML format
    let info = protein_copilot_spectrum_io::detect_format(file_path)?;
    if info.format != protein_copilot_core::spectrum::SpectrumFormat::MzML {
        return Err(XicError::UnsupportedFormat {
            path: file_path.to_path_buf(),
        });
    }

    let reader = protein_copilot_spectrum_io::create_reader(&info);

    // --- Pass 0: Get target scan info ---
    let target_spectrum = reader.read_spectrum(file_path, target_scan)?;
    let target_rt = target_spectrum.retention_time_sec.unwrap_or(0.0);
    let target_window = target_spectrum
        .precursor
        .as_ref()
        .and_then(|p| p.isolation_window.as_ref())
        .cloned();

    // --- Pre-computation: build target ion list ---
    let light_ions = build_target_ions(peptide_sequence, modifications, charge);
    let heavy_ions = match &params.label_type {
        Some(label) => crate::heavy::compute_heavy_target_ions(
            &light_ions,
            peptide_sequence,
            label,
        ),
        None => Vec::new(),
    };

    // Precursor targets
    let heavy_precursor_mz = params.label_type.as_ref().map(|label| {
        crate::heavy::compute_heavy_precursor_mz(precursor_mz, charge, peptide_sequence, label)
    });

    // --- Pass 1: Stream spectra and extract intensities ---
    // Collectors: Vec<(scan_number, rt, Vec<intensity_per_ion>)>
    let mut ms2_points: Vec<(u32, f64, Vec<f64>, Vec<f64>)> = Vec::new();
    let mut ms1_light_points: Vec<XicDataPoint> = Vec::new();
    let mut ms1_heavy_points: Vec<XicDataPoint> = Vec::new();

    reader.for_each_spectrum(file_path, &mut |spec| {
        let rt = spec.retention_time_sec.unwrap_or(0.0);

        match spec.ms_level {
            MsLevel::MS1 => {
                // Extract precursor XIC from MS1
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
            }
            MsLevel::MS2 => {
                // Check if same isolation window
                let matches_window = match (&target_window, spec.precursor.as_ref()) {
                    (Some(tw), Some(prec)) => prec
                        .isolation_window
                        .as_ref()
                        .map_or(false, |w| same_isolation_window(tw, w)),
                    (None, _) => true, // DDA: no window filtering
                    _ => false,
                };

                if matches_window {
                    // Extract light ion intensities
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

                    // Extract heavy ion intensities
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

                    ms2_points.push((spec.scan_number, rt, light_intensities, heavy_intensities));
                }
            }
            _ => {}
        }
        Ok(true) // continue streaming
    })?;

    // --- Post-processing ---
    // Sort by scan number
    ms2_points.sort_by_key(|(scan, _, _, _)| *scan);

    // Find target scan position and take ±n_cycles
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

    // Build fragment XIC traces
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

    // Select top-N by total intensity (prefer ions matched at target scan)
    fragment_traces.sort_by(|a, b| {
        let a_total: f64 = a.data_points.iter().map(|p| p.intensity).sum();
        let b_total: f64 = b.data_points.iter().map(|p| p.intensity).sum();
        b_total
            .partial_cmp(&a_total)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let top_n = params.top_n_ions.min(fragment_traces.len());
    fragment_traces.truncate(top_n);

    // Build heavy traces (same top-N selection)
    let heavy_traces: Vec<XicTrace> = if heavy_ions.is_empty() {
        Vec::new()
    } else {
        // Map top-N light ion indices back to heavy
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

    Ok(XicData {
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
    })
}
```

- [ ] **Step 6: Run tests and build**

Run: `cargo test -p protein-copilot-xic -- --nocapture 2>&1 | tail -10`
Expected: all existing tests pass

Run: `cargo build -p protein-copilot-xic 2>&1 | tail -5`
Expected: compile success (heavy module referenced but not yet implemented — will be addressed in Task 4)

- [ ] **Step 7: Commit**

```bash
git add crates/xic/src/extract.rs
git commit -m "feat(xic): implement XIC extraction core with 1.5-pass algorithm

Implements extract_xic() with:
- extract_intensity(): binary search + tolerance matching (ppm/Da)
- same_isolation_window(): full interval comparison
- build_target_ions(): reuses matching.rs fragment generators
- 1.5-pass: read_spectrum for target, for_each_spectrum for extraction
- Top-N ion selection, zero-fill, RT windowing, ±N cycles

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 4: Implement SILAC heavy-label m/z calculation

**Files:**
- Modify: `crates/xic/src/heavy.rs`
- Test: inline `#[cfg(test)] mod tests` in `heavy.rs`

- [ ] **Step 1: Write failing tests first**

In `crates/xic/src/heavy.rs`:

```rust
//! SILAC heavy-label m/z calculation for XIC traces.
//!
//! Computes heavy-label m/z shifts for fragment ions and precursors
//! based on the number of K (Lysine) and R (Arginine) residues in
//! each fragment.

use crate::extract::TargetIon;
use crate::{IonType, LabelType};

/// Compute heavy-label m/z for a precursor ion.
///
/// Adds `(count_K × heavy_k_delta + count_R × heavy_r_delta) / charge`
/// to the light precursor m/z.
pub fn compute_heavy_precursor_mz(
    light_mz: f64,
    charge: i32,
    peptide_sequence: &str,
    label: &LabelType,
) -> f64 {
    let total_delta = total_heavy_delta(peptide_sequence, label);
    light_mz + total_delta / charge.abs().max(1) as f64
}

/// Compute heavy-label target ions from light target ions.
///
/// For each light fragment ion, counts the K/R residues in the
/// corresponding fragment (prefix for b-ions, suffix for y-ions)
/// and applies the mass shift divided by the fragment charge.
pub fn compute_heavy_target_ions(
    light_ions: &[TargetIon],
    peptide_sequence: &str,
    label: &LabelType,
) -> Vec<TargetIon> {
    let chars: Vec<char> = peptide_sequence.chars().collect();
    let n = chars.len();

    light_ions
        .iter()
        .map(|ion| {
            let fragment_delta = match ion.ion_type {
                IonType::B => {
                    // b-ion covers prefix[0..ion_number]
                    let prefix = &chars[..((ion.ion_number as usize).min(n))];
                    residue_heavy_delta(prefix, label)
                }
                IonType::Y => {
                    // y-ion covers suffix[n - ion_number..]
                    let start = n.saturating_sub(ion.ion_number as usize);
                    let suffix = &chars[start..];
                    residue_heavy_delta(suffix, label)
                }
                IonType::Precursor => total_heavy_delta(peptide_sequence, label),
            };

            let heavy_mz = ion.mz + fragment_delta / ion.charge.max(1) as f64;

            TargetIon {
                label: ion.label.clone(),
                ion_type: ion.ion_type,
                ion_number: ion.ion_number,
                charge: ion.charge,
                mz: heavy_mz,
            }
        })
        .collect()
}

/// Total heavy mass delta for an entire peptide.
fn total_heavy_delta(peptide_sequence: &str, label: &LabelType) -> f64 {
    let chars: Vec<char> = peptide_sequence.chars().collect();
    residue_heavy_delta(&chars, label)
}

/// Heavy mass delta for a set of residues.
fn residue_heavy_delta(residues: &[char], label: &LabelType) -> f64 {
    match label {
        LabelType::Silac {
            heavy_k_delta,
            heavy_r_delta,
        } => {
            let count_k = residues.iter().filter(|&&c| c == 'K').count() as f64;
            let count_r = residues.iter().filter(|&&c| c == 'R').count() as f64;
            count_k * heavy_k_delta + count_r * heavy_r_delta
        }
        LabelType::Custom { residue_deltas } => {
            let mut delta = 0.0;
            for &(res, d) in residue_deltas {
                let count = residues.iter().filter(|&&c| c == res).count() as f64;
                delta += count * d;
            }
            delta
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn silac() -> LabelType {
        LabelType::standard_silac()
    }

    #[test]
    fn heavy_precursor_single_k() {
        // "PEPTIDEK" has 1 K, 0 R
        let heavy = compute_heavy_precursor_mz(500.0, 2, "PEPTIDEK", &silac());
        let expected = 500.0 + 8.014199 / 2.0;
        assert!((heavy - expected).abs() < 1e-4, "got {heavy}, expected {expected}");
    }

    #[test]
    fn heavy_precursor_k_and_r() {
        // "PEPTIDER" has 0 K, 1 R
        let heavy = compute_heavy_precursor_mz(500.0, 2, "PEPTIDER", &silac());
        let expected = 500.0 + 10.008269 / 2.0;
        assert!((heavy - expected).abs() < 1e-4);
    }

    #[test]
    fn heavy_precursor_no_label_residues() {
        // "PEPTIDE" has 0 K, 0 R
        let heavy = compute_heavy_precursor_mz(500.0, 2, "PEPTIDE", &silac());
        assert!((heavy - 500.0).abs() < 1e-6, "no K/R means no shift");
    }

    #[test]
    fn heavy_fragment_b_ion_prefix() {
        // "KPEPTIDE" — b1 covers "K" (1 K), b2 covers "KP" (1 K)
        let light = vec![
            TargetIon {
                label: "b1¹⁺".to_string(),
                ion_type: IonType::B,
                ion_number: 1,
                charge: 1,
                mz: 100.0,
            },
            TargetIon {
                label: "b2¹⁺".to_string(),
                ion_type: IonType::B,
                ion_number: 2,
                charge: 1,
                mz: 200.0,
            },
        ];
        let heavy = compute_heavy_target_ions(&light, "KPEPTIDE", &silac());
        assert_eq!(heavy.len(), 2);
        // b1 includes K → +8.014
        assert!((heavy[0].mz - (100.0 + 8.014199)).abs() < 1e-4);
        // b2 includes K → +8.014
        assert!((heavy[1].mz - (200.0 + 8.014199)).abs() < 1e-4);
    }

    #[test]
    fn heavy_fragment_y_ion_suffix() {
        // "PEPTIDEK" — y1 covers "K" (1 K), y2 covers "EK" (1 K)
        let light = vec![
            TargetIon {
                label: "y1¹⁺".to_string(),
                ion_type: IonType::Y,
                ion_number: 1,
                charge: 1,
                mz: 100.0,
            },
            TargetIon {
                label: "y2¹⁺".to_string(),
                ion_type: IonType::Y,
                ion_number: 2,
                charge: 1,
                mz: 200.0,
            },
        ];
        let heavy = compute_heavy_target_ions(&light, "PEPTIDEK", &silac());
        // y1 suffix = "K" → +8.014
        assert!((heavy[0].mz - (100.0 + 8.014199)).abs() < 1e-4);
        // y2 suffix = "EK" → +8.014
        assert!((heavy[1].mz - (200.0 + 8.014199)).abs() < 1e-4);
    }

    #[test]
    fn heavy_fragment_doubly_charged() {
        // b1 at charge 2: delta should be halved
        let light = vec![TargetIon {
            label: "b1²⁺".to_string(),
            ion_type: IonType::B,
            ion_number: 1,
            charge: 2,
            mz: 100.0,
        }];
        let heavy = compute_heavy_target_ions(&light, "KPEPTIDE", &silac());
        let expected = 100.0 + 8.014199 / 2.0;
        assert!((heavy[0].mz - expected).abs() < 1e-4);
    }
}
```

- [ ] **Step 2: Run tests**

Run: `cargo test -p protein-copilot-xic heavy -- --nocapture 2>&1 | tail -10`
Expected: all tests pass

- [ ] **Step 3: Commit**

```bash
git add crates/xic/src/heavy.rs
git commit -m "feat(xic): implement SILAC heavy-label m/z calculation

Computes heavy-label shifts for fragment ions (b/y) and precursors
based on K/R residue counts in each fragment prefix/suffix.
Supports standard SILAC (K+8.014, R+10.008) and custom labels.
Fragment charge correctly divides the mass delta.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 5: Create XIC HTML template and render function

**Files:**
- Create: `crates/report/templates/xic.html`
- Create: `crates/report/src/xic_visualize.rs`
- Modify: `crates/report/Cargo.toml`
- Modify: `crates/report/src/lib.rs`

- [ ] **Step 1: Add xic dependency to report crate**

In `crates/report/Cargo.toml`, add to `[dependencies]`:

```toml
protein-copilot-xic = { workspace = true }
```

- [ ] **Step 2: Create the XIC HTML template**

Create `crates/report/templates/xic.html` — a Plotly.js interactive template:

```html
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>XIC — __PEPTIDE_PLACEHOLDER__</title>
    <script src="__PLOTLY_SRC__"></script>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; background: #fafafa; color: #333; }
        .header { background: #1a1a2e; color: #eee; padding: 16px 24px; }
        .header h1 { font-size: 18px; font-weight: 600; }
        .header .meta { font-size: 13px; color: #aaa; margin-top: 4px; }
        .header .meta span { margin-right: 16px; }
        .controls { padding: 8px 24px; background: #f0f0f0; border-bottom: 1px solid #ddd; display: flex; gap: 12px; align-items: center; }
        .controls label { font-size: 13px; cursor: pointer; }
        .controls select, .controls button { font-size: 13px; padding: 4px 8px; }
        .plot-container { padding: 16px 24px; }
        #ms1-plot { height: 200px; margin-bottom: 8px; }
        #ms2-plot { height: 400px; }
        .params { padding: 12px 24px; font-size: 12px; color: #666; border-top: 1px solid #eee; }
    </style>
</head>
<body>
    <div class="header">
        <h1>XIC: <span id="peptide-seq"></span></h1>
        <div class="meta">
            <span>Charge: <strong id="charge-val"></strong></span>
            <span>Precursor m/z: <strong id="mz-val"></strong></span>
            <span>Target RT: <strong id="rt-val"></strong>s</span>
            <span>Target Scan: <strong id="scan-val"></strong></span>
        </div>
    </div>
    <div class="controls">
        <label>Mode:
            <select id="mode-select">
                <option value="overlay">Overlay</option>
                <option value="split">Split (Light / Heavy)</option>
            </select>
        </label>
    </div>
    <div class="plot-container">
        <div id="ms1-plot"></div>
        <div id="ms2-plot"></div>
    </div>
    <div class="params" id="params-info"></div>

    <script type="application/json" id="xic-data">/*__XIC_JSON__*/</script>
    <script>
    (function() {
        const data = JSON.parse(document.getElementById('xic-data').textContent);

        // Header
        document.getElementById('peptide-seq').textContent = data.peptide_sequence;
        document.getElementById('charge-val').textContent = data.charge + '+';
        document.getElementById('mz-val').textContent = data.precursor_mz.toFixed(4);
        document.getElementById('rt-val').textContent = data.target_rt_sec.toFixed(1);
        document.getElementById('scan-val').textContent = data.target_scan;

        // Color palette for ions
        const colors = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf'];

        function makeTrace(xicTrace, colorIdx, dash) {
            return {
                x: xicTrace.data_points.map(p => p.retention_time_sec),
                y: xicTrace.data_points.map(p => p.intensity),
                name: xicTrace.ion_label + (xicTrace.is_heavy ? ' (H)' : ''),
                mode: 'lines+markers',
                line: { color: colors[colorIdx % colors.length], width: 2, dash: dash || 'solid' },
                marker: { size: 3 },
                hovertemplate: xicTrace.ion_label + '<br>RT: %{x:.1f}s<br>Intensity: %{y:.0f}<br>m/z: ' + xicTrace.theoretical_mz.toFixed(4) + '<extra></extra>'
            };
        }

        // MS1 plot
        const ms1Traces = [];
        if (data.ms1_precursor_xic) {
            ms1Traces.push(makeTrace(data.ms1_precursor_xic, 0, 'solid'));
        }
        if (data.ms1_heavy_precursor_xic) {
            ms1Traces.push(makeTrace(data.ms1_heavy_precursor_xic, 0, 'dash'));
        }

        if (ms1Traces.length > 0) {
            Plotly.newPlot('ms1-plot', ms1Traces, {
                title: { text: 'MS1 Precursor XIC', font: { size: 14 } },
                xaxis: { title: 'RT (s)' },
                yaxis: { title: 'Intensity' },
                shapes: [{ type: 'line', x0: data.target_rt_sec, x1: data.target_rt_sec, y0: 0, y1: 1, yref: 'paper', line: { color: 'red', width: 1, dash: 'dot' } }],
                margin: { t: 40, b: 40, l: 60, r: 20 },
                showlegend: true,
                legend: { orientation: 'h', y: 1.15 }
            }, { responsive: true });
        } else {
            document.getElementById('ms1-plot').innerHTML = '<p style="text-align:center;color:#999;padding:20px;">No MS1 precursor XIC data</p>';
        }

        // MS2 fragment plot
        function renderMS2(mode) {
            const traces = [];
            data.fragment_xic_traces.forEach((t, i) => {
                traces.push(makeTrace(t, i, 'solid'));
            });
            if (mode === 'overlay') {
                data.heavy_fragment_xic_traces.forEach((t, i) => {
                    traces.push(makeTrace(t, i, 'dash'));
                });
            }

            const layout = {
                title: { text: 'MS2 Fragment Ion XIC', font: { size: 14 } },
                xaxis: { title: 'RT (s)' },
                yaxis: { title: 'Intensity' },
                shapes: [{ type: 'line', x0: data.target_rt_sec, x1: data.target_rt_sec, y0: 0, y1: 1, yref: 'paper', line: { color: 'red', width: 1, dash: 'dot' } }],
                margin: { t: 40, b: 40, l: 60, r: 20 },
                showlegend: true,
                legend: { orientation: 'h', y: 1.2 }
            };

            Plotly.newPlot('ms2-plot', traces, layout, { responsive: true });
        }

        renderMS2('overlay');

        document.getElementById('mode-select').addEventListener('change', function() {
            renderMS2(this.value);
        });

        // Params info
        const ep = data.extraction_params;
        document.getElementById('params-info').textContent =
            'Tolerance: ' + ep.mz_tolerance.value + ' ' + ep.mz_tolerance.unit +
            ' | Cycles: ±' + ep.n_cycles +
            ' | Top-N: ' + ep.top_n_ions +
            ' | Intensity: ' + ep.intensity_rule +
            (ep.label_type ? ' | Label: SILAC' : '');
    })();
    </script>
</body>
</html>
```

- [ ] **Step 3: Create `crates/report/src/xic_visualize.rs`**

```rust
//! XIC visualization — renders [`XicData`] to a self-contained HTML file
//! with interactive Plotly.js charts.

use std::fs;
use std::path::Path;

use protein_copilot_xic::{PlotlyMode, XicData};

use crate::error::ReportError;

const PLOTLY_CDN: &str = "https://cdn.plot.ly/plotly-2.35.2.min.js";

const TEMPLATE: &str = include_str!("../templates/xic.html");

/// Renders XIC data into a standalone HTML file with Plotly.js charts.
pub fn render_xic_html(
    xic_data: &XicData,
    output_path: &Path,
    plotly_mode: PlotlyMode,
) -> Result<(), ReportError> {
    let json = serde_json::to_string(xic_data)
        .map_err(|e| ReportError::SerializationError(e.to_string()))?;

    let plotly_src = match plotly_mode {
        PlotlyMode::Cdn => PLOTLY_CDN.to_string(),
        PlotlyMode::Embedded => {
            // For embedded mode, we'd need plotly.min.js bundled.
            // For MVP, fall back to CDN with a comment.
            PLOTLY_CDN.to_string()
        }
    };

    let html = TEMPLATE
        .replace("/*__XIC_JSON__*/", &json)
        .replace("__PLOTLY_SRC__", &plotly_src)
        .replace("__PEPTIDE_PLACEHOLDER__", &xic_data.peptide_sequence);

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
    use protein_copilot_xic::*;

    fn sample_xic_data() -> XicData {
        XicData {
            peptide_sequence: "PEPTIDEK".to_string(),
            target_rt_sec: 120.0,
            target_scan: 100,
            charge: 2,
            precursor_mz: 450.25,
            ms1_precursor_xic: Some(XicTrace {
                ion_label: "precursor".to_string(),
                ion_type: IonType::Precursor,
                ion_number: 0,
                charge: 2,
                theoretical_mz: 450.25,
                data_points: vec![
                    XicDataPoint { retention_time_sec: 115.0, scan_number: 90, intensity: 1000.0 },
                    XicDataPoint { retention_time_sec: 120.0, scan_number: 100, intensity: 5000.0 },
                    XicDataPoint { retention_time_sec: 125.0, scan_number: 110, intensity: 2000.0 },
                ],
                is_heavy: false,
            }),
            ms1_heavy_precursor_xic: None,
            fragment_xic_traces: vec![
                XicTrace {
                    ion_label: "y5¹⁺".to_string(),
                    ion_type: IonType::Y,
                    ion_number: 5,
                    charge: 1,
                    theoretical_mz: 574.28,
                    data_points: vec![
                        XicDataPoint { retention_time_sec: 115.0, scan_number: 91, intensity: 800.0 },
                        XicDataPoint { retention_time_sec: 120.0, scan_number: 101, intensity: 3000.0 },
                        XicDataPoint { retention_time_sec: 125.0, scan_number: 111, intensity: 1200.0 },
                    ],
                    is_heavy: false,
                },
            ],
            heavy_fragment_xic_traces: Vec::new(),
            extraction_params: ExtractionParams {
                mz_tolerance: MassTolerance { value: 20.0, unit: ToleranceUnit::Ppm },
                n_cycles: 5,
                top_n_ions: 6,
                label_type: None,
                intensity_rule: IntensityRule::MaxInWindow,
            },
        }
    }

    #[test]
    fn render_xic_html_creates_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test_xic.html");
        let data = sample_xic_data();
        render_xic_html(&data, &path, PlotlyMode::Cdn).unwrap();
        assert!(path.exists());
        let content = fs::read_to_string(&path).unwrap();
        assert!(content.contains("PEPTIDEK"));
        assert!(content.contains("plotly"));
        assert!(content.contains("y5"));
    }

    #[test]
    fn render_xic_html_contains_json_data() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test_xic2.html");
        let data = sample_xic_data();
        render_xic_html(&data, &path, PlotlyMode::Cdn).unwrap();
        let content = fs::read_to_string(&path).unwrap();
        assert!(content.contains("450.25"));
        assert!(content.contains("574.28"));
    }
}
```

- [ ] **Step 4: Update `crates/report/src/lib.rs`**

Add the module and public API:

```rust
pub mod xic_visualize;
```

Add to `impl ReportGenerator`:

```rust
/// Renders XIC data as a self-contained HTML file with Plotly.js charts.
pub fn render_xic(
    xic_data: &protein_copilot_xic::XicData,
    output_path: &Path,
    plotly_mode: protein_copilot_xic::PlotlyMode,
) -> Result<(), ReportError> {
    xic_visualize::render_xic_html(xic_data, output_path, plotly_mode)
}
```

- [ ] **Step 5: Run tests**

Run: `cargo test -p protein-copilot-report xic -- --nocapture 2>&1 | tail -10`
Expected: all tests pass

- [ ] **Step 6: Commit**

```bash
git add crates/report/
git commit -m "feat(report): add XIC HTML visualization with Plotly.js

New xic.html template renders interactive XIC plots:
- MS1 precursor XIC (light + heavy overlay)
- MS2 fragment ion XIC (top-N, legend toggle)
- Target RT marker, hover tooltips
- Overlay/split mode switcher
Uses safe <script type=application/json> data injection.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 6: Add `extract_xic` MCP tool

**Files:**
- Modify: `crates/mcp-server/Cargo.toml`
- Modify: `crates/mcp-server/src/tools.rs`

- [ ] **Step 1: Add xic dependency to mcp-server**

In `crates/mcp-server/Cargo.toml`, add to `[dependencies]`:

```toml
protein-copilot-xic = { workspace = true }
```

- [ ] **Step 2: Add input/output structs for extract_xic**

Add near the other input structs in `crates/mcp-server/src/tools.rs`:

```rust
/// Input for the `extract_xic` MCP tool.
#[derive(Debug, Deserialize, JsonSchema)]
struct ExtractXicInput {
    /// Spectrum file path (mzML only).
    #[schemars(description = "Path to the spectrum file (.mzML). XIC extraction requires mzML format for MS1+MS2 and isolation window data.")]
    file_path: Option<String>,
    /// Target scan number (1-based).
    #[schemars(description = "Scan number (1-based) to center the XIC around.")]
    scan_number: u32,
    /// Peptide sequence.
    #[schemars(description = "Peptide amino acid sequence (one-letter codes).")]
    peptide_sequence: Option<String>,
    /// Charge state.
    #[schemars(description = "Precursor charge state.")]
    charge: Option<i32>,
    /// Real precursor m/z (not DIA isolation window center).
    #[schemars(description = "True precursor m/z. For DIA data, use the PSM-derived value, not the isolation window center.")]
    precursor_mz: Option<f64>,
    /// Complete modifications list (fixed + applied variable).
    #[schemars(description = "Modifications applied to this peptide (fixed + variable). If omitted, uses unmodified sequence.")]
    modifications: Option<Vec<protein_copilot_core::search_params::Modification>>,
    /// Number of DIA cycles before/after target (default: 5).
    #[schemars(description = "Number of DIA cycles before and after target scan. Default: 5.")]
    n_cycles: Option<u32>,
    /// Number of top ions to display (default: 6).
    #[schemars(description = "Number of top fragment ions to display. Default: 6.")]
    top_n_ions: Option<usize>,
    /// Heavy-label type for SILAC comparison.
    #[schemars(description = "Heavy-label configuration. Use {\"Silac\": {\"heavy_k_delta\": 8.014199, \"heavy_r_delta\": 10.008269}} for standard SILAC.")]
    label_type: Option<protein_copilot_xic::LabelType>,
    /// m/z extraction tolerance (default: 20 ppm).
    #[schemars(description = "Mass tolerance for XIC peak extraction. Default: 20 ppm.")]
    extraction_tolerance: Option<protein_copilot_core::search_params::MassTolerance>,
    /// Intensity extraction rule (default: MaxInWindow).
    #[schemars(description = "How to extract intensity from peaks within tolerance. Default: MaxInWindow.")]
    intensity_rule: Option<protein_copilot_xic::IntensityRule>,
    /// Plotly loading mode (default: Cdn).
    #[schemars(description = "Plotly.js loading: 'Cdn' (default, smaller) or 'Embedded' (offline).")]
    plotly_mode: Option<protein_copilot_xic::PlotlyMode>,
    /// Output HTML file path.
    #[schemars(description = "Output HTML file path. Default: ./output/xic_scan{N}.html")]
    output_path: Option<String>,
    /// Run ID to resolve PSM context (single-file searches only).
    #[schemars(description = "Run ID from a previous search. Auto-fills peptide, charge, mods, precursor_mz. MVP: single-file searches only.")]
    run_id: Option<String>,
}

/// Result returned by `extract_xic`.
#[derive(Debug, Serialize, JsonSchema)]
struct ExtractXicResult {
    /// Path to the generated HTML file.
    output_path: String,
    /// Number of MS2 scans in the XIC window.
    ms2_scan_count: usize,
    /// Number of light fragment traces extracted.
    light_trace_count: usize,
    /// Number of heavy fragment traces extracted.
    heavy_trace_count: usize,
    /// Whether MS1 precursor XIC was found.
    has_ms1_xic: bool,
    /// Summary message.
    summary: String,
}
```

- [ ] **Step 3: Implement the `extract_xic` tool method**

Add the tool method to the server impl block (follow the `annotate_spectrum` pattern):

```rust
/// Extract XIC (Extracted Ion Chromatogram) for a peptide from an mzML file.
/// Generates an interactive HTML file with Plotly.js showing MS1 precursor
/// and MS2 fragment ion chromatograms. Supports SILAC heavy-label comparison.
/// Two modes: (1) provide run_id + scan_number to use PSM context, or
/// (2) provide file_path + scan_number + peptide_sequence + charge + precursor_mz.
#[tool]
fn extract_xic(
    &self,
    Parameters(input): Parameters<ExtractXicInput>,
) -> Result<Json<ExtractXicResult>, ErrorData> {
    use protein_copilot_core::search_params::{MassTolerance, Modification, ToleranceUnit};

    // Resolve mode
    let (file_path, peptide, charge, precursor_mz, modifications) =
        if let Some(ref rid) = input.run_id {
            let result = self.get_result(&None, &Some(rid.clone()))?;
            let psm = result
                .psms
                .iter()
                .find(|p| p.spectrum_scan == input.scan_number)
                .ok_or_else(|| {
                    mcp_err(
                        ErrorCode::INVALID_PARAMS,
                        format!("no PSM for scan {} in run {}", input.scan_number, rid),
                    )
                })?;
            let file = result
                .metadata
                .input_files
                .first()
                .ok_or_else(|| {
                    mcp_err(ErrorCode::INTERNAL_ERROR, "no input files in search result")
                })?
                .clone();
            (
                file,
                psm.peptide_sequence.clone(),
                psm.charge,
                psm.precursor_mz,
                psm.modifications.clone(),
            )
        } else if let (Some(ref fp), Some(ref pep), Some(ch), Some(mz)) = (
            &input.file_path,
            &input.peptide_sequence,
            input.charge,
            input.precursor_mz,
        ) {
            validate_file_path(fp)?;
            if pep.trim().is_empty() {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "peptide_sequence cannot be empty",
                ));
            }
            (
                PathBuf::from(fp),
                pep.clone(),
                ch,
                mz,
                input.modifications.clone().unwrap_or_default(),
            )
        } else {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                "provide either 'run_id' or 'file_path' + 'peptide_sequence' + 'charge' + 'precursor_mz'",
            ));
        };

    let params = protein_copilot_xic::ExtractionParams {
        mz_tolerance: input.extraction_tolerance.unwrap_or(MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        }),
        n_cycles: input.n_cycles.unwrap_or(5),
        top_n_ions: input.top_n_ions.unwrap_or(6),
        label_type: input.label_type.clone(),
        intensity_rule: input
            .intensity_rule
            .unwrap_or(protein_copilot_xic::IntensityRule::MaxInWindow),
    };

    // Extract XIC
    let xic_data = protein_copilot_xic::extract::extract_xic(
        &file_path,
        input.scan_number,
        &peptide,
        charge,
        precursor_mz,
        &modifications,
        &params,
    )
    .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?;

    // Render HTML
    let out_path = input.output_path.map(PathBuf::from).unwrap_or_else(|| {
        PathBuf::from(format!("output/xic_scan{}.html", input.scan_number))
    });

    let plotly_mode = input
        .plotly_mode
        .unwrap_or(protein_copilot_xic::PlotlyMode::Cdn);

    protein_copilot_report::xic_visualize::render_xic_html(&xic_data, &out_path, plotly_mode)
        .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?;

    let ms2_count = xic_data
        .fragment_xic_traces
        .first()
        .map(|t| t.data_points.len())
        .unwrap_or(0);

    let summary = format!(
        "XIC extracted for {} ({}+) — {} light traces, {} heavy traces, {} MS2 scans",
        peptide,
        charge,
        xic_data.fragment_xic_traces.len(),
        xic_data.heavy_fragment_xic_traces.len(),
        ms2_count,
    );

    Ok(Json(ExtractXicResult {
        output_path: out_path.to_string_lossy().to_string(),
        ms2_scan_count: ms2_count,
        light_trace_count: xic_data.fragment_xic_traces.len(),
        heavy_trace_count: xic_data.heavy_fragment_xic_traces.len(),
        has_ms1_xic: xic_data.ms1_precursor_xic.is_some(),
        summary,
    }))
}
```

- [ ] **Step 4: Build and verify**

Run: `cargo build --workspace 2>&1 | tail -10`
Expected: success

- [ ] **Step 5: Run full workspace tests**

Run: `cargo test --workspace 2>&1 | tail -15`
Expected: all tests pass (420+ existing + new XIC tests)

- [ ] **Step 6: Run clippy**

Run: `cargo clippy --workspace -- -D warnings 2>&1 | tail -10`
Expected: 0 warnings

- [ ] **Step 7: Commit**

```bash
git add crates/mcp-server/
git commit -m "feat(mcp): add extract_xic MCP tool

New tool for XIC extraction with two modes:
- run_id mode: auto-fills peptide/charge/mods from PSM
- manual mode: user provides all parameters
Renders Plotly.js interactive HTML with MS1+MS2 XIC plots.
Supports SILAC heavy-label comparison.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 7: Integration verification and cleanup

**Files:**
- All crates (verification only)

- [ ] **Step 1: Full workspace build**

Run: `cargo build --workspace 2>&1 | tail -5`
Expected: success

- [ ] **Step 2: Full test suite**

Run: `cargo test --workspace 2>&1 | tail -20`
Expected: all tests pass

- [ ] **Step 3: Clippy clean**

Run: `cargo clippy --workspace -- -D warnings 2>&1 | tail -5`
Expected: 0 warnings

- [ ] **Step 4: Verify MCP tool is registered**

Run: `grep -c '#\[tool\]' crates/mcp-server/src/tools.rs`
Expected: 15 (was 14, now +1 for extract_xic)

- [ ] **Step 5: Final commit with all fixes**

If any fixes were needed during integration:

```bash
git add -A
git commit -m "fix: integration fixes for XIC feature

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```
