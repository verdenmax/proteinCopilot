# SILAC DIA Annotation Mirror Plot + XIC Data Source Fix

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix DIA+SILAC annotation to show a mirror plot (light scan ▲ / heavy scan ▼) and fix the XIC data source bug where heavy fragment traces are extracted from the wrong MS2 scans.

**Architecture:** Move `LabelType` and heavy m/z computation to `core` crate (shared). Add `HeavyAnnotation` to `SpectrumAnnotation`. The MCP handler finds the heavy DIA scan and calls `annotate_spectrum` twice. The HTML template renders a mirror plot when `heavy_annotation` is present, and uses `ms2_heavy_scans` for XIC heavy traces.

**Tech Stack:** Rust (core, search-engine, xic, report, mcp-server crates), JavaScript (unified.html template)

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `crates/core/src/label.rs` | **Create** | `LabelType`, `IonType` (shared), heavy delta computation |
| `crates/core/src/lib.rs` | Modify | Add `pub mod label;` |
| `crates/xic/src/lib.rs` | Modify | Remove `LabelType`, `IonType`; re-export from core |
| `crates/xic/src/heavy.rs` | Modify | Use `LabelType`/`IonType` from core |
| `crates/xic/src/extract.rs` | Modify | Use `LabelType`/`IonType` from core |
| `crates/search-engine/src/annotate.rs` | Modify | Add `HeavyAnnotation` struct, `annotate_heavy_spectrum()` |
| `crates/mcp-server/src/tools.rs` | Modify | Annotation handler: find heavy scan, dual annotate |
| `crates/report/templates/unified.html` | Modify | Mirror plot rendering + XIC `ms2_heavy_scans` fix |

---

### Task 1: Move LabelType and heavy delta computation to core crate

**Files:**
- Create: `crates/core/src/label.rs`
- Modify: `crates/core/src/lib.rs`

This task extracts `LabelType` and the pure heavy mass delta functions into the `core` crate so both `xic` and `search-engine` can use them without cross-dependency.

- [ ] **Step 1: Create `crates/core/src/label.rs`**

```rust
//! SILAC heavy-label types and mass delta computation.
//!
//! Shared between `xic` (XIC extraction) and `search-engine` (annotation)
//! crates to avoid circular dependencies.

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

/// Heavy-label type for SILAC or custom isotope labeling.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "type")]
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
    /// Standard SILAC: K+8, R+10.
    pub fn standard_silac() -> Self {
        LabelType::Silac {
            heavy_k_delta: 8.014199,
            heavy_r_delta: 10.008269,
        }
    }
}

/// Heavy mass delta for a slice of residues given a label type.
pub fn residue_heavy_delta(residues: &[char], label: &LabelType) -> f64 {
    match label {
        LabelType::Silac {
            heavy_k_delta,
            heavy_r_delta,
        } => {
            let count_k = residues.iter().filter(|&&c| c == 'K' || c == 'k').count() as f64;
            let count_r = residues.iter().filter(|&&c| c == 'R' || c == 'r').count() as f64;
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

/// Total heavy mass delta for an entire peptide sequence.
pub fn total_heavy_delta(peptide_sequence: &str, label: &LabelType) -> f64 {
    let chars: Vec<char> = peptide_sequence.chars().collect();
    residue_heavy_delta(&chars, label)
}

/// Compute heavy precursor m/z from light precursor m/z.
pub fn compute_heavy_precursor_mz(
    light_mz: f64,
    charge: i32,
    peptide_sequence: &str,
    label: &LabelType,
) -> f64 {
    let delta = total_heavy_delta(peptide_sequence, label);
    light_mz + delta / charge.abs().max(1) as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_silac() {
        let label = LabelType::standard_silac();
        match &label {
            LabelType::Silac { heavy_k_delta, heavy_r_delta } => {
                assert!((heavy_k_delta - 8.014199).abs() < 1e-6);
                assert!((heavy_r_delta - 10.008269).abs() < 1e-6);
            }
            _ => panic!("expected Silac variant"),
        }
    }

    #[test]
    fn test_residue_heavy_delta_silac() {
        let label = LabelType::standard_silac();
        // "KR" → K+8 + R+10 = 18.022468
        let chars: Vec<char> = "KR".chars().collect();
        let delta = residue_heavy_delta(&chars, &label);
        assert!((delta - 18.022468).abs() < 1e-4);
    }

    #[test]
    fn test_no_kr_no_delta() {
        let label = LabelType::standard_silac();
        let delta = total_heavy_delta("AGLDEF", &label);
        assert!((delta - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_heavy_precursor_mz() {
        let label = LabelType::standard_silac();
        // DGFLLDGFPR: 1 R → delta = 10.008269, charge 2 → +5.004
        let heavy = compute_heavy_precursor_mz(569.3, 2, "DGFLLDGFPR", &label);
        assert!((heavy - 569.3 - 10.008269 / 2.0).abs() < 1e-4);
    }
}
```

- [ ] **Step 2: Add `pub mod label;` to `crates/core/src/lib.rs`**

Add this line after the existing `pub mod` declarations (after line 21):

```rust
pub mod label;
```

- [ ] **Step 3: Run tests to verify**

Run: `cargo test -p protein-copilot-core --quiet`
Expected: All existing tests pass + 4 new label tests pass.

- [ ] **Step 4: Commit**

```bash
git add crates/core/src/label.rs crates/core/src/lib.rs
git commit -m "feat: add LabelType and heavy delta computation to core crate"
```

---

### Task 2: Re-export LabelType from core in xic crate

**Files:**
- Modify: `crates/xic/src/lib.rs` — remove `LabelType` and `IonType` enum definitions, re-export from core
- Modify: `crates/xic/src/heavy.rs` — use `LabelType` from core, remove duplicated `residue_heavy_delta`/`total_heavy_delta`
- Modify: `crates/xic/src/extract.rs` — update imports if needed

This task makes xic crate use core's `LabelType` instead of its own copy. `IonType` stays in xic since the xic `IonType` has a `Precursor` variant not present in search-engine's `IonType`, and moving it would complicate things. Only `LabelType` and the delta functions move.

- [ ] **Step 1: Modify `crates/xic/src/lib.rs`**

Remove the `LabelType` enum definition (lines ~123-146 including `impl LabelType`) and replace with a re-export:

```rust
// Replace the LabelType enum and impl block with:
pub use protein_copilot_core::label::LabelType;
```

Keep `IonType` as-is in xic (it has `Precursor` variant specific to XIC).

- [ ] **Step 2: Modify `crates/xic/src/heavy.rs`**

Replace the private functions `total_heavy_delta` and `residue_heavy_delta` with imports from core:

```rust
// At the top of heavy.rs, change:
use crate::{IonType, LabelType};
// To:
use crate::IonType;
use protein_copilot_core::label::{LabelType, residue_heavy_delta, total_heavy_delta};
```

Then **delete** the two private functions `total_heavy_delta` (lines 70-73) and `residue_heavy_delta` (lines 76-95) from `heavy.rs`, since they now come from core.

Update `compute_heavy_precursor_mz` to use `total_heavy_delta` from core (it's already imported, just remove the local function).

Update `compute_heavy_target_ions` — the body calls `residue_heavy_delta` which now comes from the import.

- [ ] **Step 3: Run tests to verify nothing broke**

Run: `cargo test -p protein-copilot-xic --quiet`
Expected: All 27+ xic tests pass. The heavy.rs tests use the same function signatures.

- [ ] **Step 4: Run clippy**

Run: `cargo clippy -p protein-copilot-xic --all-targets -- -D warnings`
Expected: Clean.

- [ ] **Step 5: Commit**

```bash
git add crates/xic/src/lib.rs crates/xic/src/heavy.rs crates/xic/src/extract.rs
git commit -m "refactor: use LabelType from core crate in xic"
```

---

### Task 3: Add HeavyAnnotation to search-engine annotate.rs

**Files:**
- Modify: `crates/search-engine/src/annotate.rs`
- Modify: `crates/search-engine/Cargo.toml` (no change needed — core already a dependency)

- [ ] **Step 1: Add `HeavyAnnotation` struct after `SpectrumAnnotation` (after line 123)**

```rust
/// Annotation data from a heavy-label DIA scan.
///
/// When SILAC + DIA is active, the heavy precursor falls in a different
/// isolation window, so its fragments appear in a different MS2 scan.
/// This struct holds the annotation of that separate scan.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct HeavyAnnotation {
    /// Scan number of the heavy MS2 scan.
    pub scan_number: u32,
    /// Retention time of the heavy scan in minutes.
    pub retention_time_min: f64,
    /// Heavy precursor m/z.
    pub precursor_mz: f64,
    /// All experimental peaks in the heavy scan with optional annotations.
    pub peaks: Vec<AnnotatedPeak>,
    /// Theoretical heavy b-ions with match status.
    pub b_ions: Vec<TheoreticalIon>,
    /// Theoretical heavy y-ions with match status.
    pub y_ions: Vec<TheoreticalIon>,
    /// Match score (matched_ions / total_ions).
    pub score: f64,
    /// Number of matched heavy fragment ions.
    pub matched_ions: u32,
    /// Total number of theoretical heavy fragment ions.
    pub total_ions: u32,
}
```

- [ ] **Step 2: Add `heavy_annotation` field to `SpectrumAnnotation`**

Add after the `modifications` field (line 122):

```rust
    /// Heavy-label annotation from a separate DIA scan (DIA+SILAC only).
    /// `None` for DDA or when no heavy label is configured.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub heavy_annotation: Option<HeavyAnnotation>,
```

- [ ] **Step 3: Set `heavy_annotation: None` in existing `annotate_spectrum` function**

In the `Ok(SpectrumAnnotation { ... })` return at the bottom of `annotate_spectrum()`, add:

```rust
        heavy_annotation: None,
```

- [ ] **Step 4: Add `annotate_heavy_spectrum` public function**

Add this function after the existing `annotate_spectrum`:

```rust
/// Annotate a heavy-label spectrum for mirror plot display.
///
/// Given a heavy MS2 spectrum and the peptide sequence + label type,
/// computes heavy theoretical fragment m/z values and matches them
/// against the heavy spectrum's peaks.
///
/// This is called separately from `annotate_spectrum` because in DIA mode,
/// the heavy fragments appear in a different MS2 scan.
pub fn annotate_heavy_spectrum(
    heavy_spectrum: &Spectrum,
    peptide_sequence: &str,
    charge: i32,
    fragment_tolerance: &MassTolerance,
    fixed_modifications: &[Modification],
    label: &protein_copilot_core::label::LabelType,
) -> Result<HeavyAnnotation, SearchEngineError> {
    use protein_copilot_core::label::{compute_heavy_precursor_mz, residue_heavy_delta};

    if charge <= 0 {
        return Err(SearchEngineError::ExecutionError {
            detail: format!("charge must be >= 1, got {charge}"),
        });
    }
    if heavy_spectrum.mz_array.is_empty() {
        return Err(SearchEngineError::ExecutionError {
            detail: "heavy spectrum has no peaks".to_string(),
        });
    }

    let precursor = heavy_spectrum
        .precursors
        .first()
        .ok_or_else(|| SearchEngineError::ExecutionError {
            detail: "heavy spectrum has no precursor information".to_string(),
        })?;

    // Compute heavy precursor m/z
    let light_precursor_mz = {
        let neutral = peptide_mass(peptide_sequence).ok_or_else(|| {
            SearchEngineError::ExecutionError {
                detail: format!("cannot compute mass for '{}'", peptide_sequence),
            }
        })?;
        let mod_delta = apply_fixed_mod_mass(peptide_sequence, fixed_modifications);
        peptide_mz(neutral + mod_delta, charge)
    };
    let heavy_precursor_mz =
        compute_heavy_precursor_mz(light_precursor_mz, charge, peptide_sequence, label);

    // Generate heavy theoretical fragments
    let max_frag_charge: u32 = if charge >= 3 { 2 } else { 1 };
    let chars: Vec<char> = peptide_sequence.chars().collect();
    let n = chars.len();

    // Heavy b-ions
    let light_b = generate_b_entries(peptide_sequence, fixed_modifications, max_frag_charge)
        .ok_or_else(|| SearchEngineError::ExecutionError {
            detail: format!("cannot generate b-ions for '{}'", peptide_sequence),
        })?;
    let heavy_b_entries: Vec<FragmentEntry> = light_b
        .iter()
        .map(|e| {
            let prefix = &chars[..(e.ion_number as usize).min(n)];
            let delta = residue_heavy_delta(prefix, label);
            FragmentEntry {
                ion_type: e.ion_type,
                ion_number: e.ion_number,
                charge: e.charge,
                mz: e.mz + delta / e.charge as f64,
            }
        })
        .collect();

    // Heavy y-ions
    let light_y = generate_y_entries(peptide_sequence, fixed_modifications, max_frag_charge)
        .ok_or_else(|| SearchEngineError::ExecutionError {
            detail: format!("cannot generate y-ions for '{}'", peptide_sequence),
        })?;
    let heavy_y_entries: Vec<FragmentEntry> = light_y
        .iter()
        .map(|e| {
            let start = n.saturating_sub(e.ion_number as usize);
            let suffix = &chars[start..];
            let delta = residue_heavy_delta(suffix, label);
            FragmentEntry {
                ion_type: e.ion_type,
                ion_number: e.ion_number,
                charge: e.charge,
                mz: e.mz + delta / e.charge as f64,
            }
        })
        .collect();

    // Match against heavy spectrum peaks (reuse existing matching logic)
    let exp_mz = &heavy_spectrum.mz_array;
    let exp_int = &heavy_spectrum.intensity_array;
    let mut used_peaks = vec![false; exp_mz.len()];
    let mut matched_count = 0u32;
    let total_count = (heavy_b_entries.len() + heavy_y_entries.len()) as u32;

    // Match heavy b-ions
    let heavy_b_ions: Vec<TheoreticalIon> = heavy_b_entries
        .iter()
        .map(|e| {
            let (matched, best_idx) = find_best_match(e.mz, exp_mz, fragment_tolerance);
            let (matched_mz, delta_ppm) = if let Some(idx) = best_idx {
                if !used_peaks[idx] {
                    used_peaks[idx] = true;
                    matched_count += 1;
                    let obs = exp_mz[idx];
                    let ppm = (obs - e.mz) / e.mz * 1e6;
                    (Some(obs), Some(ppm))
                } else {
                    (None, None)
                }
            } else {
                (None, None)
            };
            TheoreticalIon {
                ion_type: e.ion_type,
                number: e.ion_number,
                charge: e.charge,
                theoretical_mz: e.mz,
                matched: matched && matched_mz.is_some(),
                matched_mz,
                delta_ppm,
            }
        })
        .collect();

    // Match heavy y-ions
    let heavy_y_ions: Vec<TheoreticalIon> = heavy_y_entries
        .iter()
        .map(|e| {
            let (matched, best_idx) = find_best_match(e.mz, exp_mz, fragment_tolerance);
            let (matched_mz, delta_ppm) = if let Some(idx) = best_idx {
                if !used_peaks[idx] {
                    used_peaks[idx] = true;
                    matched_count += 1;
                    let obs = exp_mz[idx];
                    let ppm = (obs - e.mz) / e.mz * 1e6;
                    (Some(obs), Some(ppm))
                } else {
                    (None, None)
                }
            } else {
                (None, None)
            };
            TheoreticalIon {
                ion_type: e.ion_type,
                number: e.ion_number,
                charge: e.charge,
                theoretical_mz: e.mz,
                matched: matched && matched_mz.is_some(),
                matched_mz,
                delta_ppm,
            }
        })
        .collect();

    // Annotate peaks
    let peaks: Vec<AnnotatedPeak> = exp_mz
        .iter()
        .zip(exp_int.iter())
        .enumerate()
        .map(|(i, (&mz, &intensity))| {
            let annotation = if used_peaks[i] {
                // Find which ion matched this peak
                heavy_b_ions
                    .iter()
                    .chain(heavy_y_ions.iter())
                    .find(|ion| ion.matched_mz == Some(mz))
                    .map(|ion| IonAnnotation {
                        ion_type: ion.ion_type,
                        ion_number: ion.number,
                        charge: ion.charge,
                        theoretical_mz: ion.theoretical_mz,
                        delta_mz: (mz - ion.theoretical_mz).abs(),
                        delta_ppm: ion.delta_ppm.unwrap_or(0.0),
                    })
            } else {
                None
            };
            AnnotatedPeak {
                mz,
                intensity,
                annotation,
            }
        })
        .collect();

    let score = if total_count > 0 {
        matched_count as f64 / total_count as f64
    } else {
        0.0
    };

    Ok(HeavyAnnotation {
        scan_number: heavy_spectrum.scan_number,
        retention_time_min: heavy_spectrum.retention_time_min,
        precursor_mz: heavy_precursor_mz,
        peaks,
        b_ions: heavy_b_ions,
        y_ions: heavy_y_ions,
        score,
        matched_ions: matched_count,
        total_ions: total_count,
    })
}
```

- [ ] **Step 5: Add unit test for `annotate_heavy_spectrum`**

Add in the existing `#[cfg(test)] mod tests` block:

```rust
#[test]
fn test_annotate_heavy_spectrum_basic() {
    use protein_copilot_core::label::LabelType;
    use protein_copilot_core::spectrum::{PrecursorInfo, Spectrum};

    // Heavy spectrum with a peak near heavy b2+ m/z for "AK" prefix
    // AK: A=71.03711 + K=128.09496 = 199.13207 neutral b2
    // Heavy K delta = 8.014199 → heavy b2 neutral = 207.146269
    // b2+ m/z = (207.146269 + 1.007276) / 1 = 208.153545
    let heavy_b2_mz = 208.1535;

    let spectrum = Spectrum {
        scan_number: 100,
        ms_level: 2,
        retention_time_min: 25.0,
        mz_array: vec![heavy_b2_mz, 300.0, 500.0],
        intensity_array: vec![1000.0, 200.0, 150.0],
        precursors: vec![PrecursorInfo {
            mz: 300.0,
            charge: Some(2),
            intensity: Some(5000.0),
            isolation_window: None,
            source_scan: None,
        }],
    };

    let label = LabelType::standard_silac();
    let tolerance = MassTolerance {
        value: 20.0,
        unit: ToleranceUnit::Ppm,
    };

    let result = annotate_heavy_spectrum(
        &spectrum, "AKDEF", 2, &tolerance, &[], &label,
    )
    .unwrap();

    assert_eq!(result.scan_number, 100);
    assert!(result.matched_ions >= 1, "should match at least heavy b2+");
    assert!(result.total_ions > 0);
    assert!(result.score > 0.0);
    // Verify a b-ion was matched
    let matched_b = result.b_ions.iter().filter(|i| i.matched).count();
    assert!(matched_b >= 1, "should have at least one matched heavy b-ion");
}
```

- [ ] **Step 6: Fix any compilation issues in other crates**

The new `heavy_annotation` field on `SpectrumAnnotation` will break any code that constructs the struct directly. Search for construction sites:

Run: `cargo build 2>&1 | grep "heavy_annotation" | head -10`

Expected sites to fix: test code in `report` and `mcp-server` crates that construct `SpectrumAnnotation` — add `heavy_annotation: None` to each.

- [ ] **Step 7: Run tests**

Run: `cargo test --quiet`
Expected: All tests pass including the new one.

- [ ] **Step 8: Commit**

```bash
git add crates/search-engine/src/annotate.rs crates/report/ crates/mcp-server/
git commit -m "feat: add HeavyAnnotation and annotate_heavy_spectrum for DIA+SILAC mirror plot"
```

---

### Task 4: MCP handler — find heavy DIA scan and dual-annotate

**Files:**
- Modify: `crates/mcp-server/src/tools.rs` (annotate_spectrum handler, ~line 1440-1730)

The handler currently reads a single spectrum and calls `annotate_spectrum` once. For DIA+SILAC, it needs to:
1. Detect DIA mode (check isolation window width)
2. If `label_type` is set AND DIA: find the heavy scan via `find_heavy_dia_window_from_spectra`
3. Read the heavy spectrum
4. Call `annotate_heavy_spectrum` on it
5. Attach the result as `annotation.heavy_annotation`

- [ ] **Step 1: Add heavy scan lookup after annotation (around line 1567)**

After the existing `annotate_spectrum` call succeeds and before the output_path logic, insert:

```rust
        // ── DIA + SILAC: find and annotate heavy scan ──
        if let Some(ref label) = input.label_type {
            // Check if the target spectrum is DIA (wide isolation window)
            let is_dia = spectrum
                .precursors
                .first()
                .and_then(|p| p.isolation_window.as_ref())
                .map(|w| (w.lower_offset + w.upper_offset) > 3.0)
                .unwrap_or(false);

            if is_dia {
                // Compute heavy precursor m/z
                let core_label = convert_label_type(label);
                let heavy_prec_mz = protein_copilot_core::label::compute_heavy_precursor_mz(
                    annotation.precursor_mz,
                    charge,
                    &peptide_seq,
                    &core_label,
                );

                // Read nearby spectra to find the heavy DIA window
                let nearby = reader
                    .read_spectra_near_rt(
                        &spectrum_file,
                        spectrum.retention_time_min,
                        1.0, // ±1 min window
                    )
                    .unwrap_or_default();

                if let Some((heavy_scan_num, _window)) =
                    protein_copilot_xic::heavy::find_heavy_dia_window_from_spectra(
                        &nearby,
                        spectrum.retention_time_min,
                        heavy_prec_mz,
                    )
                {
                    // Read the heavy spectrum
                    if let Ok(heavy_spectrum) =
                        reader.read_spectrum(&spectrum_file, heavy_scan_num)
                    {
                        match protein_copilot_search_engine::annotate::annotate_heavy_spectrum(
                            &heavy_spectrum,
                            &peptide_seq,
                            charge,
                            &frag_tol,
                            &modifications,
                            &core_label,
                        ) {
                            Ok(heavy_ann) => {
                                annotation.heavy_annotation = Some(heavy_ann);
                            }
                            Err(e) => {
                                tracing::warn!("Heavy annotation failed: {e}");
                            }
                        }
                    }
                }
            }
        }
```

- [ ] **Step 2: Add `convert_label_type` helper function**

The MCP server uses `protein_copilot_xic::LabelType` (which is now re-exported from core). If the types are the same (re-export), no conversion is needed — just use it directly. Verify with:

```rust
// If xic::LabelType IS core::label::LabelType (re-export), this is enough:
let core_label: &protein_copilot_core::label::LabelType = label;
```

If they're the same type via re-export, simplify Step 1 to use `label` directly instead of `convert_label_type(label)`.

- [ ] **Step 3: Add `read_spectra_near_rt` to SpectrumReader trait (if not present)**

Check if the spectrum reader has a method to read multiple spectra near a given RT. If not, use the existing `read_spectrum` approach — iterate nearby scan numbers:

```rust
// Alternative: read spectra in a scan range around the target
let scan_range_start = resolved_scan.saturating_sub(50);
let scan_range_end = resolved_scan + 50;
let mut nearby_spectra = Vec::new();
for s in scan_range_start..=scan_range_end {
    if let Ok(spec) = reader.read_spectrum(&spectrum_file, s) {
        if spec.ms_level == 2 {
            nearby_spectra.push(spec);
        }
    }
}
```

Use `find_heavy_dia_window_from_spectra(&nearby_spectra, ...)` on this collection.

- [ ] **Step 4: Run build to verify compilation**

Run: `cargo build -p protein-copilot-mcp-server 2>&1 | tail -5`
Expected: Compiles successfully.

- [ ] **Step 5: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "feat: MCP annotation handler finds heavy DIA scan and dual-annotates"
```

---

### Task 5: Fix XIC heavy data source bug in unified.html

**Files:**
- Modify: `crates/report/templates/unified.html` (line ~822-829)

This is the one-line bug: heavy fragment XIC traces are extracted from `ms2_scans` (light window) instead of `ms2_heavy_scans` (heavy window).

- [ ] **Step 1: Fix the ms2Heavy computation (line 822-829)**

Replace:

```javascript
    var ms2Heavy = heavyIons.map(function(ion) {
      if (!U.raw_scans) return [];
      return U.raw_scans.ms2_scans.map(function(s) {
        var r = extractIntensity(ion.heavy_mz, s.mz_array, s.intensity_array, tol);
        return { rt: s.retention_time_min, scan: s.scan_number,
                 intensity: r.intensity, observed_mz: r.observed_mz };
      });
    });
```

With:

```javascript
    // Use ms2_heavy_scans for DIA+SILAC (different isolation window),
    // fall back to ms2_scans for DDA (same window).
    var heavyMs2Source = (U.raw_scans.ms2_heavy_scans && U.raw_scans.ms2_heavy_scans.length > 0)
      ? U.raw_scans.ms2_heavy_scans
      : U.raw_scans.ms2_scans;
    var ms2Heavy = heavyIons.map(function(ion) {
      if (!U.raw_scans) return [];
      return heavyMs2Source.map(function(s) {
        var r = extractIntensity(ion.heavy_mz, s.mz_array, s.intensity_array, tol);
        return { rt: s.retention_time_min, scan: s.scan_number,
                 intensity: r.intensity, observed_mz: r.observed_mz };
      });
    });
```

- [ ] **Step 2: Commit**

```bash
git add crates/report/templates/unified.html
git commit -m "fix: use ms2_heavy_scans for DIA+SILAC XIC heavy traces"
```

---

### Task 6: HTML mirror plot rendering for annotation

**Files:**
- Modify: `crates/report/templates/unified.html`

When `heavy_annotation` is present on the annotation data, render a mirror plot instead of a single spectrum plot. Light peaks go upward (existing), heavy peaks go downward (new).

- [ ] **Step 1: Add heavy annotation data variable**

In the annotation section (Section A, around line 239), add after `var D = ...`:

```javascript
  var DH = D.heavy_annotation || null;  // Heavy annotation data (DIA+SILAC)
```

- [ ] **Step 2: Update the info panel for mirror plot mode (line ~297)**

After the existing `add("Scan / RT", ...)` line, add heavy scan info:

```javascript
    if (DH) {
      add("Heavy Scan / RT", "Scan " + DH.scan_number + "  \u2014  " + fmt(DH.retention_time_min, 2) + " min");
      add("Heavy Precursor", fmt(DH.precursor_mz, 4) + " m/z");
      var heavyBadgeCls = "score-badge " + (DH.score > 0.5 ? "good" : DH.score >= 0.2 ? "mid" : "low");
      var heavyBadge = el("span", { "class": heavyBadgeCls }, fmt(DH.score, 3));
      var heavyScoreWrap = document.createElement("span");
      heavyScoreWrap.appendChild(heavyBadge);
      heavyScoreWrap.appendChild(document.createTextNode("  (" + DH.matched_ions + "/" + DH.total_ions + " ions)"));
      add("Heavy Score", heavyScoreWrap);
    }
```

- [ ] **Step 3: Add heavy coverage row to the coverage panel**

After the existing coverage SVG rendering (after `panel.appendChild(covSvg)`, line ~418), add:

```javascript
    // Heavy coverage (below light coverage, if available)
    if (DH) {
      var heavyLabel = el("div", { style: "margin-top:12px; font-size:13px; color:#e67e22; font-weight:600;" },
        "\u25BC Heavy — Scan " + DH.scan_number + " (" + DH.matched_ions + "/" + DH.total_ions + " ions)");
      panel.appendChild(heavyLabel);

      var covSvgH = svgEl("svg", {
        width: "100%", height: "auto",
        viewBox: "0 0 " + totalW + " " + svgH,
        preserveAspectRatio: "xMidYMid meet",
        style: "max-width:" + totalW + "px"
      });

      // Build heavy match lookups
      var bMatchedH = [], yMatchedH = [];
      var bNumbersH = [], yNumbersH = [];
      for (var ch = 0; ch < n - 1; ch++) {
        var bIonH = DH.b_ions[ch];
        bMatchedH[ch] = !!(bIonH && bIonH.matched);
        bNumbersH[ch] = bIonH ? bIonH.number : (ch + 1);
        var yIdxH = n - 2 - ch;
        var yIonH = DH.y_ions[yIdxH];
        yMatchedH[ch] = !!(yIonH && yIonH.matched);
        yNumbersH[ch] = yIonH ? yIonH.number : (n - 1 - ch);
      }

      // Amino acid letters (same sequence)
      for (var sh = 0; sh < n; sh++) {
        var axh = xOff + sh * cellW + cellW / 2;
        var aah = svgEl("text", {
          x: axh, y: seqY, "text-anchor": "middle", "dominant-baseline": "central",
          "font-size": "20", "font-family": "'Courier New', Courier, monospace",
          "font-weight": "700", fill: "#e67e22"
        });
        aah.textContent = seq[sh];
        covSvgH.appendChild(aah);
      }

      // Draw heavy brackets (orange colors)
      for (var ch2 = 0; ch2 < n - 1; ch2++) {
        var cleavXH = xOff + (ch2 + 1) * cellW;
        var rightCXH = cleavXH + cellW / 2;
        var leftCXH = cleavXH - cellW / 2;
        var bColH = bMatchedH[ch2] ? "#d35400" : "#e0e0e0";
        var bStH = bMatchedH[ch2] ? 2.5 : 1.2;
        var yColH = yMatchedH[ch2] ? "#e67e22" : "#e0e0e0";
        var yStH = yMatchedH[ch2] ? 2.5 : 1.2;

        // y-brackets (above)
        covSvgH.appendChild(svgEl("line", { x1: cleavXH, y1: yVertTop, x2: cleavXH, y2: yVertBot, stroke: yColH, "stroke-width": yStH, "stroke-linecap": "round" }));
        covSvgH.appendChild(svgEl("line", { x1: cleavXH, y1: yHorizY, x2: rightCXH, y2: yHorizY, stroke: yColH, "stroke-width": yStH, "stroke-linecap": "round" }));
        var yLblH = svgEl("text", { x: cleavXH + 2, y: yLabelY, "font-size": "10", "font-weight": "600", fill: yColH, "text-anchor": "start" });
        yLblH.textContent = "y" + subscript(yNumbersH[ch2]);
        covSvgH.appendChild(yLblH);

        // b-brackets (below)
        covSvgH.appendChild(svgEl("line", { x1: cleavXH, y1: bVertTop, x2: cleavXH, y2: bVertBot, stroke: bColH, "stroke-width": bStH, "stroke-linecap": "round" }));
        covSvgH.appendChild(svgEl("line", { x1: cleavXH, y1: bHorizY, x2: leftCXH, y2: bHorizY, stroke: bColH, "stroke-width": bStH, "stroke-linecap": "round" }));
        var bLblH = svgEl("text", { x: cleavXH - 2, y: bLabelY, "font-size": "10", "font-weight": "600", fill: bColH, "text-anchor": "end" });
        bLblH.textContent = "b" + subscript(bNumbersH[ch2]);
        covSvgH.appendChild(bLblH);
      }
      panel.appendChild(covSvgH);
    }
```

- [ ] **Step 4: Modify `renderSpectrum()` to draw mirror plot when `DH` present**

This is the largest change. The existing `renderSpectrum()` function (starts line ~422) draws peaks upward. When `DH` is present, it needs to:
1. Split the SVG into upper half (light, peaks up) and lower half (heavy, peaks down)
2. Share the same m/z x-axis
3. Draw a center axis line separating light/heavy

Replace the spectrum rendering section. Key changes:

a) **Double the SVG height** when mirror mode:
```javascript
    var mirrorMode = !!DH;
    var H = mirrorMode ? 800 : 500;
    // Upper half for light, lower half for heavy
    var lightH = mirrorMode ? (H - margin.top - margin.bottom) / 2 : (H - margin.top - margin.bottom);
    var heavyH = mirrorMode ? lightH : 0;
```

b) **Compute max intensity across BOTH spectra** for shared scale:
```javascript
    var heavyPeaks = DH ? DH.peaks : [];
    if (mirrorMode) {
      for (var hp = 0; hp < heavyPeaks.length; hp++) {
        if (heavyPeaks[hp].mz > maxMz) maxMz = heavyPeaks[hp].mz;
        if (heavyPeaks[hp].intensity > maxInt) maxInt = heavyPeaks[hp].intensity;
      }
    }
```

c) **Draw light peaks upward** from center axis (same as now, but from mid-point):
```javascript
    var lightBaseline = margin.top + lightH;  // Center line
    function yScaleLight(v) { return lightBaseline - (v / intMax) * lightH; }
```

d) **Draw heavy peaks downward** from center axis:
```javascript
    var heavyBaseline = lightBaseline;  // Same center line
    function yScaleHeavy(v) { return heavyBaseline + (v / intMax) * heavyH; }
```

e) **Label the halves:**
```javascript
    if (mirrorMode) {
      // Light label (top-left)
      var lightLabel = svgEl("text", { x: margin.left + 5, y: margin.top + 15,
        "font-size": "13", "font-weight": "700", fill: "#0984e3" });
      lightLabel.textContent = "\u25B2 Light — Scan " + D.scan_number;
      svg.appendChild(lightLabel);

      // Heavy label (bottom-left)
      var heavyLabel = svgEl("text", { x: margin.left + 5, y: lightBaseline + 15,
        "font-size": "13", "font-weight": "700", fill: "#e67e22" });
      heavyLabel.textContent = "\u25BC Heavy — Scan " + DH.scan_number;
      svg.appendChild(heavyLabel);

      // Center divider line
      svg.appendChild(svgEl("line", {
        x1: margin.left, y1: lightBaseline, x2: margin.left + pw, y2: lightBaseline,
        stroke: "#666", "stroke-width": "1.5"
      }));
    }
```

f) **Draw heavy peaks** (new loop after the light peaks loop):
```javascript
    if (mirrorMode && heavyPeaks.length > 0) {
      for (var ph2 = 0; ph2 < heavyPeaks.length; ph2++) {
        var pkH = heavyPeaks[ph2];
        var relIntH = pkH.intensity / intMax;
        var xH = xScale(pkH.mz);
        var annH = pkH.annotation;

        if (annH) {
          var peakHH = relIntH * heavyH;
          if (peakHH < minMatchedH) peakHH = minMatchedH;
          var y2mH = heavyBaseline + peakHH;
          var colH = annH.ion_type === "B" ? "#d35400" : "#e67e22";  // Orange tones

          gMatched.appendChild(svgEl("line", {
            x1: xH, y1: heavyBaseline, x2: xH, y2: y2mH,
            stroke: colH, "stroke-width": "3", "stroke-linecap": "round"
          }));

          var ionPrefH = annH.ion_type === "B" ? "b" : "y";
          var labelTxtH = ionPrefH + subscript(annH.ion_number);
          var labelYH = y2mH + 16;
          if (labelYH > heavyBaseline + heavyH - 4) labelYH = heavyBaseline + heavyH - 4;

          var ltH = svgEl("text", {
            x: xH, y: labelYH, "text-anchor": "middle",
            fill: colH, "font-size": "14", "font-weight": "700"
          });
          ltH.textContent = labelTxtH;
          gLabels.appendChild(ltH);

        } else if (relIntH < 0.05) {
          gNoise.appendChild(svgEl("line", {
            x1: xH, y1: heavyBaseline, x2: xH, y2: heavyBaseline + relIntH * heavyH,
            stroke: "#f5e6d0", "stroke-width": "0.5"
          }));
        } else {
          gUnmatched.appendChild(svgEl("line", {
            x1: xH, y1: heavyBaseline, x2: xH, y2: heavyBaseline + relIntH * heavyH,
            stroke: "#ead0b0", "stroke-width": "1"
          }));
        }

        // Hitbox
        (function(peak, px, ann2) {
          var peakTopH = ann2
            ? heavyBaseline + Math.max((peak.intensity / intMax) * heavyH, minMatchedH)
            : heavyBaseline + (peak.intensity / intMax) * heavyH;
          var hitboxH = svgEl("line", {
            x1: px, y1: heavyBaseline, x2: px, y2: peakTopH,
            stroke: "transparent", "stroke-width": "12", fill: "none", cursor: "pointer"
          });
          hitboxH.addEventListener("mouseenter", function() {
            var relI = (peak.intensity / intMax * 100);
            var html = "<div class=\"tip-row\"><span class=\"tip-label\">m/z</span><span class=\"tip-val\">" + fmt(peak.mz, 4) + "</span></div>";
            html += "<div class=\"tip-row\"><span class=\"tip-label\">Intensity</span><span class=\"tip-val\">" + fmt(peak.intensity, 1) + " (" + fmt(relI, 1) + "%)</span></div>";
            html += "<div class=\"tip-row\"><span class=\"tip-label\">Source</span><span class=\"tip-val\" style=\"color:#e67e22\">Heavy Scan " + DH.scan_number + "</span></div>";
            if (ann2) {
              html += "<hr class=\"tip-sep\">";
              var lbl = (ann2.ion_type === "B" ? "b" : "y") + "<sub>" + ann2.ion_number + "</sub>";
              html += "<div class=\"tip-row\"><span class=\"tip-label\">Ion</span><span class=\"tip-val\">" + lbl + " (Heavy)</span></div>";
              html += "<div class=\"tip-row\"><span class=\"tip-label\">\u0394 ppm</span><span class=\"tip-val\">" + fmt(ann2.delta_ppm, 2) + "</span></div>";
            }
            tooltip.innerHTML = html;
            tooltip.style.display = "block";
          });
          hitboxH.addEventListener("mousemove", function(ev) {
            tooltip.style.left = (ev.clientX + 14) + "px";
            tooltip.style.top = (ev.clientY - 28) + "px";
          });
          hitboxH.addEventListener("mouseleave", function() { tooltip.style.display = "none"; });
          gHitboxes.appendChild(hitboxH);
        })(pkH, xH, annH);
      }
    }
```

- [ ] **Step 5: Commit**

```bash
git add crates/report/templates/unified.html
git commit -m "feat: mirror plot rendering for DIA+SILAC annotation (light ▲ / heavy ▼)"
```

---

### Task 7: Final verification

- [ ] **Step 1: Run clippy on all crates**

Run: `cargo clippy --all-targets -- -D warnings 2>&1 | grep -E "^error" | head -10`
Expected: No errors (pre-existing diann.rs warning is acceptable).

- [ ] **Step 2: Run full test suite**

Run: `cargo test --quiet`
Expected: All 510+ tests pass, 0 failures.

- [ ] **Step 3: End-to-end test with real data**

Test with the DGFLLDGFPR peptide that originally demonstrated the bug:

```
annotate_spectrum(
  file_path: "/home/verden/pfind/2025-fall/code/2da/hela-mix-2da.mzML",
  scan_number: 0,
  retention_time_min: <RT from previous export>,
  peptide_sequence: "DGFLLDGFPR",
  charge: 2,
  label_type: { "Silac": { "heavy_k_delta": 8.014199, "heavy_r_delta": 10.008269 } },
  output_path: "output/test_mirror_DGFLLDGFPR.html"
)
```

Open the output HTML and verify:
- Mirror plot shows light peaks upward, heavy peaks downward
- Heavy peaks are from a different scan number
- Coverage panel shows both light and heavy brackets
- XIC traces use correct data sources

- [ ] **Step 4: Commit any fixes**

```bash
git add -A
git commit -m "fix: final adjustments for DIA+SILAC mirror plot"
```

---

## Notes

- **DDA mode is unchanged**: When `heavy_annotation` is `None` (DDA or no SILAC), the HTML renders exactly as before — single upward spectrum, no mirror.
- **`IonType` stays in two places**: `xic::IonType` has `Precursor` variant, `annotate::IonType` has only `B`/`Y`. These serve different purposes and merging them would be complex for no benefit.
- **The MCP handler scan range** (±50 scans) is a heuristic. DIA cycles typically have 20-40 MS2 scans, so ±50 should cover at least one full cycle to find the heavy window.
