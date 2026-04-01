# Spectrum Annotation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let users annotate a single spectrum with b/y fragment ion matching and generate an interactive HTML visualization.

**Architecture:** `search-engine/annotate.rs` computes `SpectrumAnnotation` (reusing ion generation from matching.rs). `report/visualize.rs` renders it to a self-contained HTML file with SVG + JS interactivity. `annotate_spectrum` MCP tool provides two input modes (from existing PSM or manual peptide).

**Tech Stack:** Rust, serde_json (JSON embedding in HTML), SVG + vanilla JS (visualization), rmcp (MCP tool)

**Spec:** `docs/superpowers/specs/2026-03-31-spectrum-annotation-design.md`

---

## File Map

| Action | File | Responsibility |
|--------|------|----------------|
| Modify | `crates/search-engine/src/matching.rs` | Make `generate_b_ions`, `generate_y_ions`, `within_tolerance` public |
| Create | `crates/search-engine/src/annotate.rs` | `SpectrumAnnotation` data structures + `annotate_spectrum()` function |
| Modify | `crates/search-engine/src/lib.rs` | Export `annotate` module |
| Create | `crates/report/templates/annotation.html` | Self-contained HTML template with embedded JS for rendering |
| Create | `crates/report/src/visualize.rs` | `render_annotation_html()` function |
| Modify | `crates/report/src/lib.rs` | Export `visualize` module, add method to `ReportGenerator` |
| Modify | `crates/mcp-server/src/tools.rs` | `annotate_spectrum` MCP tool |
| Modify | `.github/agents/proteomics-search.agent.md` | Add tool + usage scenarios |

---

### Task 1: Make matching.rs ion generation functions public

**Files:**
- Modify: `crates/search-engine/src/matching.rs:46,66,86`

- [ ] **Step 1: Change visibility of three functions**

In `crates/search-engine/src/matching.rs`, change:

Line 46: `fn within_tolerance(` → `pub fn within_tolerance(`
Line 66: `fn generate_b_ions(` → `pub fn generate_b_ions(`
Line 86: `fn generate_y_ions(` → `pub fn generate_y_ions(`

- [ ] **Step 2: Build to verify no breakage**

Run: `cargo build -p protein-copilot-search-engine 2>&1 | tail -3`
Expected: `Finished` with exit 0.

- [ ] **Step 3: Commit**

```bash
git add crates/search-engine/src/matching.rs
git commit -m "refactor(search-engine): make ion generation and tolerance functions public

Needed by annotate module for spectrum annotation feature.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: Create annotate module with SpectrumAnnotation

**Files:**
- Create: `crates/search-engine/src/annotate.rs`
- Modify: `crates/search-engine/src/lib.rs:24-37`

- [ ] **Step 1: Create `crates/search-engine/src/annotate.rs` with data structures**

```rust
//! Spectrum annotation — detailed b/y ion matching for a single spectrum.
//!
//! Given a spectrum and a peptide sequence, generates a [`SpectrumAnnotation`]
//! containing per-peak ion labels, theoretical ion lists with match status,
//! and overall matching statistics.

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use protein_copilot_core::search_params::{MassTolerance, Modification, ToleranceUnit};
use protein_copilot_core::spectrum::Spectrum;

use crate::chemistry::{peptide_mass, peptide_mz, PROTON_MASS};
use crate::error::SearchEngineError;
use crate::matching::{generate_b_ions, generate_y_ions, within_tolerance};

/// Ion type for fragment annotation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub enum IonType {
    B,
    Y,
}

impl std::fmt::Display for IonType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IonType::B => write!(f, "b"),
            IonType::Y => write!(f, "y"),
        }
    }
}

/// Annotation for a single matched ion.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct IonAnnotation {
    pub ion_type: IonType,
    pub ion_number: u32,
    pub theoretical_mz: f64,
    pub delta_mz: f64,
    pub delta_ppm: f64,
}

/// A peak in the experimental spectrum with optional annotation.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct AnnotatedPeak {
    pub mz: f64,
    pub intensity: f64,
    pub annotation: Option<IonAnnotation>,
}

/// A theoretical ion with its match status.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct TheoreticalIon {
    pub ion_type: IonType,
    pub number: u32,
    pub theoretical_mz: f64,
    pub matched: bool,
    pub matched_mz: Option<f64>,
    pub delta_ppm: Option<f64>,
}

/// Complete annotation result for a single spectrum–peptide match.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SpectrumAnnotation {
    pub scan_number: u32,
    pub retention_time_sec: f64,
    pub peptide_sequence: String,
    pub charge: i32,
    pub precursor_mz: f64,
    pub theoretical_mz: f64,
    pub delta_mass_ppm: f64,
    pub score: f64,
    pub matched_ions: u32,
    pub total_ions: u32,
    pub protein_accessions: Vec<String>,
    pub peaks: Vec<AnnotatedPeak>,
    pub b_ions: Vec<TheoreticalIon>,
    pub y_ions: Vec<TheoreticalIon>,
    pub modifications: Vec<Modification>,
}

/// Annotate a single spectrum against a given peptide sequence.
///
/// Generates theoretical b/y ions, matches them against experimental peaks,
/// and produces a detailed [`SpectrumAnnotation`] suitable for visualization.
pub fn annotate_spectrum(
    spectrum: &Spectrum,
    peptide_sequence: &str,
    charge: i32,
    fragment_tolerance: &MassTolerance,
    fixed_modifications: &[Modification],
    protein_accessions: Vec<String>,
) -> Result<SpectrumAnnotation, SearchEngineError> {
    // Validate inputs
    if spectrum.mz_array.is_empty() {
        return Err(SearchEngineError::ExecutionError {
            detail: "spectrum has no peaks".to_string(),
        });
    }
    if spectrum.precursors.is_empty() {
        return Err(SearchEngineError::ExecutionError {
            detail: "spectrum has no precursor information".to_string(),
        });
    }

    let precursor = &spectrum.precursors[0];

    // Calculate theoretical precursor m/z
    let mut neutral_mass = match peptide_mass(peptide_sequence) {
        Some(m) => m,
        None => {
            return Err(SearchEngineError::ExecutionError {
                detail: format!(
                    "peptide sequence '{}' contains non-standard amino acids",
                    peptide_sequence
                ),
            })
        }
    };

    // Apply fixed modifications
    for modification in fixed_modifications {
        for residue in &modification.residues {
            let count = peptide_sequence.chars().filter(|c| c == residue).count();
            neutral_mass += modification.mass_delta * count as f64;
        }
    }

    let theoretical_mz = peptide_mz(neutral_mass, charge);
    let observed_mz = precursor.mz;
    let delta_mass_ppm = (observed_mz - theoretical_mz) / theoretical_mz * 1e6;

    // Generate theoretical b/y ions
    let b_mzs = generate_b_ions(peptide_sequence);
    let y_mzs = generate_y_ions(peptide_sequence);

    // Match each theoretical ion against experimental peaks
    let exp_mz = &spectrum.mz_array;
    let exp_int = &spectrum.intensity_array;

    let mut b_ions = Vec::with_capacity(b_mzs.len());
    let mut y_ions = Vec::with_capacity(y_mzs.len());

    // Track which experimental peaks are annotated (index → annotation)
    let mut peak_annotations: Vec<Option<IonAnnotation>> = vec![None; exp_mz.len()];

    // Match b-ions
    for (i, &theo_mz) in b_mzs.iter().enumerate() {
        let (matched, matched_idx) = find_best_match(theo_mz, exp_mz, fragment_tolerance);
        let ion_number = (i + 1) as u32;

        if let Some(idx) = matched_idx {
            let obs = exp_mz[idx];
            let dppm = (obs - theo_mz) / theo_mz * 1e6;
            b_ions.push(TheoreticalIon {
                ion_type: IonType::B,
                number: ion_number,
                theoretical_mz: theo_mz,
                matched: true,
                matched_mz: Some(obs),
                delta_ppm: Some(dppm),
            });
            // Annotate the experimental peak (only if not already annotated)
            if peak_annotations[idx].is_none() {
                peak_annotations[idx] = Some(IonAnnotation {
                    ion_type: IonType::B,
                    ion_number,
                    theoretical_mz: theo_mz,
                    delta_mz: obs - theo_mz,
                    delta_ppm: dppm,
                });
            }
        } else {
            b_ions.push(TheoreticalIon {
                ion_type: IonType::B,
                number: ion_number,
                theoretical_mz: theo_mz,
                matched,
                matched_mz: None,
                delta_ppm: None,
            });
        }
    }

    // Match y-ions
    for (i, &theo_mz) in y_mzs.iter().enumerate() {
        let (matched, matched_idx) = find_best_match(theo_mz, exp_mz, fragment_tolerance);
        let ion_number = (i + 1) as u32;

        if let Some(idx) = matched_idx {
            let obs = exp_mz[idx];
            let dppm = (obs - theo_mz) / theo_mz * 1e6;
            y_ions.push(TheoreticalIon {
                ion_type: IonType::Y,
                number: ion_number,
                theoretical_mz: theo_mz,
                matched: true,
                matched_mz: Some(obs),
                delta_ppm: Some(dppm),
            });
            if peak_annotations[idx].is_none() {
                peak_annotations[idx] = Some(IonAnnotation {
                    ion_type: IonType::Y,
                    ion_number,
                    theoretical_mz: theo_mz,
                    delta_mz: obs - theo_mz,
                    delta_ppm: dppm,
                });
            }
        } else {
            y_ions.push(TheoreticalIon {
                ion_type: IonType::Y,
                number: ion_number,
                theoretical_mz: theo_mz,
                matched,
                matched_mz: None,
                delta_ppm: None,
            });
        }
    }

    // Build annotated peaks list
    let peaks: Vec<AnnotatedPeak> = exp_mz
        .iter()
        .zip(exp_int.iter())
        .zip(peak_annotations.into_iter())
        .map(|((&mz, &intensity), annotation)| AnnotatedPeak {
            mz,
            intensity,
            annotation,
        })
        .collect();

    let matched_count = b_ions.iter().filter(|i| i.matched).count()
        + y_ions.iter().filter(|i| i.matched).count();
    let total_count = b_ions.len() + y_ions.len();
    let score = if total_count > 0 {
        matched_count as f64 / total_count as f64
    } else {
        0.0
    };

    Ok(SpectrumAnnotation {
        scan_number: spectrum.scan_number,
        retention_time_sec: spectrum.retention_time_sec,
        peptide_sequence: peptide_sequence.to_string(),
        charge,
        precursor_mz: observed_mz,
        theoretical_mz,
        delta_mass_ppm,
        score,
        matched_ions: matched_count as u32,
        total_ions: total_count as u32,
        protein_accessions,
        peaks,
        b_ions,
        y_ions,
        modifications: fixed_modifications.to_vec(),
    })
}

/// Find the best matching experimental peak for a theoretical m/z.
///
/// Returns `(matched, Option<index>)` where index is the position in
/// the sorted experimental m/z array.
fn find_best_match(
    theoretical_mz: f64,
    exp_mz: &[f64],
    tolerance: &MassTolerance,
) -> (bool, Option<usize>) {
    let pos = exp_mz.partition_point(|&x| x < theoretical_mz);

    let candidates = [pos.checked_sub(1), Some(pos), Some(pos + 1)]
        .into_iter()
        .flatten()
        .filter(|&i| i < exp_mz.len());

    let mut best: Option<(usize, f64)> = None;
    for idx in candidates {
        if within_tolerance(exp_mz[idx], theoretical_mz, tolerance) {
            let diff = (exp_mz[idx] - theoretical_mz).abs();
            if best.map_or(true, |(_, d)| diff < d) {
                best = Some((idx, diff));
            }
        }
    }

    match best {
        Some((idx, _)) => (true, Some(idx)),
        None => (false, None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::spectrum::{MsLevel, PrecursorInfo, Spectrum};

    fn test_spectrum() -> Spectrum {
        // A spectrum with peaks at positions matching b/y ions of "PEPTIDE"
        let b_ions = generate_b_ions("PEPTIDE");
        let y_ions = generate_y_ions("PEPTIDE");

        // Create experimental peaks: all b/y ions + some noise peaks
        let mut mz_array: Vec<f64> = Vec::new();
        let mut intensity_array: Vec<f64> = Vec::new();

        // Add b-ions with slight offsets
        for &mz in &b_ions {
            mz_array.push(mz + 0.001); // small offset
            intensity_array.push(1000.0);
        }
        // Add y-ions
        for &mz in &y_ions {
            mz_array.push(mz - 0.001);
            intensity_array.push(800.0);
        }
        // Add noise peaks
        mz_array.push(150.0);
        intensity_array.push(200.0);
        mz_array.push(350.0);
        intensity_array.push(100.0);

        // Sort by m/z (required)
        let mut pairs: Vec<(f64, f64)> = mz_array
            .into_iter()
            .zip(intensity_array.into_iter())
            .collect();
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));

        let mz_array: Vec<f64> = pairs.iter().map(|p| p.0).collect();
        let intensity_array: Vec<f64> = pairs.iter().map(|p| p.1).collect();

        let theoretical_mz = peptide_mz(peptide_mass("PEPTIDE").unwrap(), 2);

        Spectrum::new(
            1,
            MsLevel::MS2,
            30.0,
            vec![PrecursorInfo {
                mz: theoretical_mz,
                charge: Some(2),
                intensity: None,
                isolation_window: None,
            }],
            mz_array,
            intensity_array,
        )
        .unwrap()
    }

    #[test]
    fn annotate_matches_all_ions() {
        let spectrum = test_spectrum();
        let tol = MassTolerance {
            value: 0.05,
            unit: ToleranceUnit::Da,
        };
        let result = annotate_spectrum(&spectrum, "PEPTIDE", 2, &tol, &[], vec![]).unwrap();

        assert_eq!(result.peptide_sequence, "PEPTIDE");
        assert_eq!(result.charge, 2);
        assert!(result.matched_ions > 0);
        assert!(result.score > 0.5, "Should match most ions, got {}", result.score);

        // All b-ions should be matched
        assert!(
            result.b_ions.iter().all(|i| i.matched),
            "All b-ions should match: {:?}",
            result.b_ions.iter().map(|i| i.matched).collect::<Vec<_>>()
        );
        // All y-ions should be matched
        assert!(
            result.y_ions.iter().all(|i| i.matched),
            "All y-ions should match"
        );
    }

    #[test]
    fn annotate_no_match() {
        let spectrum = test_spectrum();
        let tol = MassTolerance {
            value: 0.05,
            unit: ToleranceUnit::Da,
        };
        // Use a completely different peptide
        let result = annotate_spectrum(&spectrum, "GGGGGGGG", 2, &tol, &[], vec![]).unwrap();
        assert!(result.score < 0.5, "Should have low score for wrong peptide");
    }

    #[test]
    fn annotate_nonstandard_residue_errors() {
        let spectrum = test_spectrum();
        let tol = MassTolerance {
            value: 0.05,
            unit: ToleranceUnit::Da,
        };
        let result = annotate_spectrum(&spectrum, "PEPT*DE", 2, &tol, &[], vec![]);
        assert!(result.is_err());
    }

    #[test]
    fn annotated_peaks_include_noise() {
        let spectrum = test_spectrum();
        let tol = MassTolerance {
            value: 0.05,
            unit: ToleranceUnit::Da,
        };
        let result = annotate_spectrum(&spectrum, "PEPTIDE", 2, &tol, &[], vec![]).unwrap();

        // Should have some unannotated peaks (noise)
        let unannotated = result.peaks.iter().filter(|p| p.annotation.is_none()).count();
        assert!(unannotated > 0, "Should have some unannotated noise peaks");

        // Annotated peaks should have correct ion info
        let annotated: Vec<_> = result
            .peaks
            .iter()
            .filter(|p| p.annotation.is_some())
            .collect();
        assert!(!annotated.is_empty());
        for peak in annotated {
            let ann = peak.annotation.as_ref().unwrap();
            assert!(ann.delta_mz.abs() < 0.1);
        }
    }

    #[test]
    fn annotation_serde_roundtrip() {
        let spectrum = test_spectrum();
        let tol = MassTolerance {
            value: 0.05,
            unit: ToleranceUnit::Da,
        };
        let result = annotate_spectrum(&spectrum, "PEPTIDE", 2, &tol, &[], vec![]).unwrap();
        let json = serde_json::to_string(&result).unwrap();
        let back: SpectrumAnnotation = serde_json::from_str(&json).unwrap();
        assert_eq!(result.scan_number, back.scan_number);
        assert_eq!(result.matched_ions, back.matched_ions);
    }
}
```

- [ ] **Step 2: Add `pub mod annotate;` to `crates/search-engine/src/lib.rs`**

After line 24 (`pub mod adapters;`), add:

```rust
pub mod annotate;
```

- [ ] **Step 3: Add `serde_json` to search-engine dependencies (needed for test)**

In `crates/search-engine/Cargo.toml`, move `serde_json` from dev-dependencies to dependencies:

```toml
[dependencies]
# ... existing deps ...
serde_json = { workspace = true }
```

Remove the duplicate from `[dev-dependencies]` if present.

- [ ] **Step 4: Build and run tests**

Run: `cargo test -p protein-copilot-search-engine -- annotate --nocapture 2>&1 | tail -15`
Expected: 5 tests pass (annotate_matches_all_ions, annotate_no_match, annotate_nonstandard_residue_errors, annotated_peaks_include_noise, annotation_serde_roundtrip).

- [ ] **Step 5: Commit**

```bash
git add -A && git commit -m "feat(search-engine): add annotate module for single spectrum annotation

SpectrumAnnotation data structure with per-peak ion labels.
annotate_spectrum() reuses matching.rs ion generation.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: Create HTML template and render function

**Files:**
- Create: `crates/report/templates/annotation.html`
- Create: `crates/report/src/visualize.rs`
- Modify: `crates/report/src/lib.rs:11-15`

- [ ] **Step 1: Create `crates/report/templates/annotation.html`**

This is the self-contained HTML template. It expects a `window.__ANNOTATION_DATA__` JSON object injected by Rust. The template renders:
1. Info panel (scan, peptide, score)
2. Peptide sequence coverage diagram (b/y ion coverage)
3. Interactive spectrum plot (SVG with hover tooltips)

Create `crates/report/templates/` directory and the HTML file with full content. The template should:
- Read `window.__ANNOTATION_DATA__` for the SpectrumAnnotation JSON
- Render an SVG spectrum plot with colored peaks (gray=unmatched, red=b-ion, blue=y-ion)
- Add text labels above matched peaks (e.g. "b3", "y7")
- Render the peptide sequence coverage diagram above the spectrum
- Add mouseover tooltips showing m/z, intensity, delta ppm
- Be entirely self-contained (no external CSS/JS dependencies)
- Work in any modern browser

- [ ] **Step 2: Create `crates/report/src/visualize.rs`**

```rust
//! Spectrum annotation visualization — renders SpectrumAnnotation to HTML.

use std::fs;
use std::path::Path;

use protein_copilot_search_engine::annotate::SpectrumAnnotation;

use crate::error::ReportError;

const TEMPLATE: &str = include_str!("../templates/annotation.html");

/// Render a spectrum annotation as an interactive HTML file.
pub fn render_annotation_html(
    annotation: &SpectrumAnnotation,
    output_path: &Path,
) -> Result<(), ReportError> {
    let json = serde_json::to_string(annotation).map_err(|e| ReportError::SerializationError(e.to_string()))?;

    let html = TEMPLATE.replace(
        "/*__ANNOTATION_JSON__*/",
        &format!("window.__ANNOTATION_DATA__ = {};", json),
    );

    if let Some(parent) = output_path.parent() {
        fs::create_dir_all(parent).map_err(|e| ReportError::IoError {
            path: parent.to_path_buf(),
            detail: e.to_string(),
        })?;
    }

    fs::write(output_path, html).map_err(|e| ReportError::IoError {
        path: output_path.to_path_buf(),
        detail: e.to_string(),
    })?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_search_engine::annotate::{AnnotatedPeak, TheoreticalIon, IonType};
    use protein_copilot_core::search_params::Modification;

    fn test_annotation() -> SpectrumAnnotation {
        SpectrumAnnotation {
            scan_number: 1,
            retention_time_sec: 30.0,
            peptide_sequence: "PEPTIDE".to_string(),
            charge: 2,
            precursor_mz: 400.0,
            theoretical_mz: 400.001,
            delta_mass_ppm: -2.5,
            score: 0.75,
            matched_ions: 9,
            total_ions: 12,
            protein_accessions: vec!["P001".to_string()],
            peaks: vec![
                AnnotatedPeak { mz: 100.0, intensity: 500.0, annotation: None },
                AnnotatedPeak { mz: 200.0, intensity: 1000.0, annotation: None },
            ],
            b_ions: vec![TheoreticalIon {
                ion_type: IonType::B, number: 1, theoretical_mz: 98.06,
                matched: true, matched_mz: Some(98.061), delta_ppm: Some(10.2),
            }],
            y_ions: vec![TheoreticalIon {
                ion_type: IonType::Y, number: 1, theoretical_mz: 148.06,
                matched: false, matched_mz: None, delta_ppm: None,
            }],
            modifications: vec![],
        }
    }

    #[test]
    fn render_creates_html_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.html");
        render_annotation_html(&test_annotation(), &path).unwrap();
        assert!(path.exists());
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("__ANNOTATION_DATA__"));
        assert!(content.contains("PEPTIDE"));
    }
}
```

- [ ] **Step 3: Add `protein-copilot-search-engine` dependency to report Cargo.toml**

In `crates/report/Cargo.toml`, add:

```toml
protein-copilot-search-engine = { workspace = true }
```

- [ ] **Step 4: Update `crates/report/src/lib.rs` to export visualize module**

After `pub mod summary;`, add:

```rust
pub mod visualize;
```

Add a method to `ReportGenerator`:

```rust
pub fn render_annotation(
    annotation: &protein_copilot_search_engine::annotate::SpectrumAnnotation,
    output_path: &std::path::Path,
) -> Result<(), ReportError> {
    visualize::render_annotation_html(annotation, output_path)
}
```

- [ ] **Step 5: Build and run tests**

Run: `cargo test -p protein-copilot-report -- visualize --nocapture 2>&1 | tail -10`
Expected: `render_creates_html_file` test passes.

- [ ] **Step 6: Commit**

```bash
git add -A && git commit -m "feat(report): add spectrum annotation HTML visualization

Self-contained HTML template with SVG spectrum plot,
peptide coverage diagram, and hover tooltips.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 4: Add annotate_spectrum MCP tool

**Files:**
- Modify: `crates/mcp-server/src/tools.rs`
- Modify: `.github/agents/proteomics-search.agent.md`

- [ ] **Step 1: Add AnnotateSpectrumInput struct and AnnotateResult struct**

Add near the other input structs (around line 183 in tools.rs):

```rust
#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct AnnotateSpectrumInput {
    /// Run ID — use to annotate an existing PSM from a search result.
    #[serde(default)]
    run_id: Option<String>,
    /// Spectrum file path — use for manual annotation mode.
    #[serde(default)]
    file_path: Option<String>,
    /// Scan number (1-based) to annotate.
    scan_number: u32,
    /// Peptide sequence — required for manual mode, ignored if run_id provided.
    #[serde(default)]
    peptide_sequence: Option<String>,
    /// Charge state — required for manual mode, ignored if run_id provided.
    #[serde(default)]
    charge: Option<i32>,
    /// Output HTML file path. Default: ./annotation_scan{N}.html
    #[serde(default)]
    output_path: Option<String>,
    /// Fragment mass tolerance. Default: 0.02 Da.
    #[serde(default)]
    fragment_tolerance: Option<protein_copilot_core::search_params::MassTolerance>,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct AnnotateResult {
    output_path: String,
    scan_number: u32,
    peptide_sequence: String,
    charge: i32,
    score: f64,
    matched_ions: u32,
    total_ions: u32,
    delta_mass_ppm: f64,
    protein_accessions: Vec<String>,
    message: String,
}
```

- [ ] **Step 2: Add annotate_spectrum tool implementation**

Add the tool after `list_searches` (after line 847):

```rust
/// Annotate a single spectrum with b/y ion matching and generate interactive HTML.
#[rmcp::tool(
    name = "annotate_spectrum",
    description = "Annotate a single spectrum with peptide fragment ion matching. Generates an interactive HTML file showing matched b/y ions. Two modes: (1) provide run_id + scan_number to annotate an existing PSM, or (2) provide file_path + scan_number + peptide_sequence + charge for manual annotation."
)]
fn annotate_spectrum(
    &self,
    Parameters(input): Parameters<AnnotateSpectrumInput>,
) -> Result<Json<AnnotateResult>, ErrorData> {
    use protein_copilot_search_engine::annotate;

    let default_tol = protein_copilot_core::search_params::MassTolerance {
        value: 0.02,
        unit: protein_copilot_core::search_params::ToleranceUnit::Da,
    };
    let frag_tol = input.fragment_tolerance.unwrap_or(default_tol);

    // Resolve spectrum + peptide info based on input mode
    let (spectrum, peptide_seq, charge, mods, proteins) = if let Some(ref rid) = input.run_id {
        // Mode 1: From existing search result
        let result = self.get_result(&None, &Some(rid.clone()))?;
        let psm = result
            .psms
            .iter()
            .find(|p| p.spectrum_scan == input.scan_number)
            .ok_or_else(|| {
                mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    format!("no PSM found for scan {} in this search", input.scan_number),
                )
            })?;

        // Read the spectrum from the first input file
        let file_path = result
            .metadata
            .input_files
            .first()
            .ok_or_else(|| mcp_err(ErrorCode::INTERNAL_ERROR, "search has no input files"))?;
        let info = protein_copilot_spectrum_io::detect_format(file_path)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        let reader = protein_copilot_spectrum_io::create_reader(&info);
        let spectrum = reader
            .read_spectrum(file_path, input.scan_number)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

        (
            spectrum,
            psm.peptide_sequence.clone(),
            psm.charge,
            psm.modifications.clone(),
            psm.protein_accessions.clone(),
        )
    } else if let (Some(ref fp), Some(ref seq), Some(ch)) =
        (&input.file_path, &input.peptide_sequence, input.charge)
    {
        // Mode 2: Manual annotation
        let path = Path::new(fp);
        let info = protein_copilot_spectrum_io::detect_format(path)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        let reader = protein_copilot_spectrum_io::create_reader(&info);
        let spectrum = reader
            .read_spectrum(path, input.scan_number)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

        (spectrum, seq.clone(), ch, vec![], vec![])
    } else {
        return Err(mcp_err(
            ErrorCode::INVALID_PARAMS,
            "provide either run_id (for existing PSM) or file_path + peptide_sequence + charge (for manual annotation)",
        ));
    };

    // Run annotation
    let annotation = annotate::annotate_spectrum(
        &spectrum, &peptide_seq, charge, &frag_tol, &mods, proteins,
    )
    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

    // Generate HTML
    let output = input
        .output_path
        .unwrap_or_else(|| format!("./annotation_scan{}.html", input.scan_number));
    let output_path = Path::new(&output);

    protein_copilot_report::ReportGenerator::render_annotation(&annotation, output_path)
        .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

    Ok(Json(AnnotateResult {
        output_path: output.clone(),
        scan_number: annotation.scan_number,
        peptide_sequence: annotation.peptide_sequence,
        charge: annotation.charge,
        score: annotation.score,
        matched_ions: annotation.matched_ions,
        total_ions: annotation.total_ions,
        delta_mass_ppm: annotation.delta_mass_ppm,
        protein_accessions: annotation.protein_accessions,
        message: format!(
            "Annotation saved to {}. Open in browser to view interactive spectrum.",
            output
        ),
    }))
}
```

- [ ] **Step 3: Update Agent instructions**

In `.github/agents/proteomics-search.agent.md`, add `annotate_spectrum` to the tools list and add usage scenario:

Tools list:
```yaml
  - annotate_spectrum
```

Add section:
```markdown
## 谱图标注

当用户想查看某一张谱图的匹配详情时：
  - 用户说"看一下 scan 1234 的匹配情况"
    → 调用 annotate_spectrum(run_id=xxx, scan_number=1234)
    → 告知用户"标注文件已生成，请在浏览器中打开 xxx.html 查看"
    → 基于 score/matched_ions 给出简短解读

  - 用户说"用 PEPTIDEK 去匹配 scan 100"
    → 调用 annotate_spectrum(file_path=xxx, scan_number=100, peptide_sequence="PEPTIDEK", charge=2)
    → 展示匹配结果
```

- [ ] **Step 4: Build and test full workspace**

Run: `cargo fmt && cargo build --workspace && cargo clippy --workspace && cargo test --workspace 2>&1 | tail -5`
Expected: All tests pass, 0 clippy warnings, fmt clean.

- [ ] **Step 5: Commit**

```bash
git add -A && git commit -m "feat(mcp-server): add annotate_spectrum tool with HTML visualization

Two modes: from existing PSM (run_id) or manual (file_path + peptide).
Generates interactive HTML with SVG spectrum plot and b/y coverage.
12 MCP tools total.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 5: End-to-end test and verification

**Files:**
- Modify: `crates/search-engine/tests/e2e_integration.rs`

- [ ] **Step 1: Add annotation e2e test**

Add to `crates/search-engine/tests/e2e_integration.rs`:

```rust
#[tokio::test]
async fn scenario_annotate_psm() {
    use protein_copilot_search_engine::annotate;

    // Run a search first
    let file_info = detect_format(&mgf_path()).unwrap();
    let summary = create_reader(&file_info).read_summary(&mgf_path()).unwrap();
    let mut params = ParamRecommender.recommend(&summary, None).unwrap().decision;
    params.database_path = fasta_path().to_string_lossy().to_string();

    let engine = SimpleSearchEngine::new();
    let result = engine
        .search(&params, &[mgf_path()], noop_progress())
        .await
        .unwrap();

    // Take the first PSM and annotate it
    let psm = result.psms.first().expect("should have at least one PSM");
    let spectrum = create_reader(&detect_format(&mgf_path()).unwrap())
        .read_spectrum(&mgf_path(), psm.spectrum_scan)
        .unwrap();

    let tol = protein_copilot_core::search_params::MassTolerance {
        value: 0.02,
        unit: protein_copilot_core::search_params::ToleranceUnit::Da,
    };

    let annotation = annotate::annotate_spectrum(
        &spectrum,
        &psm.peptide_sequence,
        psm.charge,
        &tol,
        &psm.modifications,
        psm.protein_accessions.clone(),
    )
    .unwrap();

    assert_eq!(annotation.scan_number, psm.spectrum_scan);
    assert_eq!(annotation.peptide_sequence, psm.peptide_sequence);
    assert!(annotation.matched_ions > 0);
    assert!(!annotation.peaks.is_empty());
    assert!(!annotation.b_ions.is_empty());
    assert!(!annotation.y_ions.is_empty());

    // Render to HTML
    let dir = tempfile::tempdir().unwrap();
    let html_path = dir.path().join("annotation.html");
    ReportGenerator::render_annotation(&annotation, &html_path).unwrap();
    assert!(html_path.exists());

    let content = std::fs::read_to_string(&html_path).unwrap();
    assert!(content.contains("__ANNOTATION_DATA__"));
    assert!(content.contains(&psm.peptide_sequence));
}
```

- [ ] **Step 2: Run full verification**

```bash
cargo build --workspace
cargo test --workspace
cargo clippy --workspace
cargo fmt --check
```

Expected: All pass, 0 warnings.

- [ ] **Step 3: Run the full pipeline example with annotation**

```bash
cargo run -p protein-copilot-search-engine --example full_search -- \
  crates/search-engine/tests/fixtures_e2e/test_100.mgf \
  crates/search-engine/tests/fixtures_e2e/test_100.fasta \
  /tmp/e2e_output
```

Then verify the output works.

- [ ] **Step 4: Commit**

```bash
git add -A && git commit -m "test: add annotation e2e test, verify full pipeline

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```
