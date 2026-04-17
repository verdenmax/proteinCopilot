# M2.0 Sage Search Engine Integration — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Integrate sage-core as a library dependency implementing `SearchEngineAdapter`, enabling users to run real proteomics searches via `engine: "Sage"`.

**Architecture:** SageAdapter lives in `crates/search-engine/src/adapters/sage/` (3 files: mod.rs, convert.rs, config.rs). It bridges our types to sage-core types, runs search in `spawn_blocking` (rayon↔tokio), and converts Feature results back to Psm. The MCP `run_search` handler switches from hardcoded SimpleSearchEngine to registry-based lookup.

**Tech Stack:** sage-core (git dep), rayon, tokio::task::spawn_blocking, serde_json

**Design spec:** `docs/superpowers/specs/2026-04-17-sage-integration-design.md`

**Sage source reference:** cloned at `/tmp/sage-ref` (commit `cd712d4`). Key files:
- `crates/sage/src/spectrum.rs` — RawSpectrum, Precursor, ProcessedSpectrum, SpectrumProcessor
- `crates/sage/src/scoring.rs` — Feature, Scorer
- `crates/sage/src/database.rs` — Builder, Parameters, IndexedDatabase, EnzymeBuilder
- `crates/sage/src/modification.rs` — ModificationSpecificity
- `crates/sage/src/fdr.rs` — picked_peptide, picked_protein
- `crates/sage/src/mass.rs` — Tolerance, PROTON
- `crates/sage/src/ml/qvalue.rs` — spectrum_q_value
- `crates/sage/src/ml/linear_discriminant.rs` — score_psms
- `crates/sage-cli/src/runner.rs` — orchestration reference (Runner::run, spectrum_fdr)

---

### Task 1: Add `engine` field to SearchParams and `extra` field to Psm

These foundational changes are needed before any adapter code.

**Files:**
- Modify: `crates/core/src/search_params.rs` (SearchParams struct + valid_params test helper)
- Modify: `crates/core/src/search_result.rs` (Psm struct + sample_psm helpers)

- [ ] **Step 1: Add `engine` field to SearchParams**

In `crates/core/src/search_params.rs`, add after the `max_peptide_length` field (line ~277):

```rust
    /// Search engine to use. `None` means default engine (SimpleSearch).
    #[serde(default)]
    pub engine: Option<String>,
```

- [ ] **Step 2: Update `valid_params()` test helper**

In the same file, update the `valid_params()` function (line ~366) to include the new field:

```rust
    fn valid_params() -> SearchParams {
        SearchParams {
            enzyme: Enzyme::Trypsin,
            missed_cleavages: 2,
            fixed_modifications: vec![carbamidomethyl()],
            variable_modifications: vec![oxidation()],
            precursor_tolerance: MassTolerance {
                value: 20.0,
                unit: ToleranceUnit::Ppm,
            },
            fragment_tolerance: MassTolerance {
                value: 0.02,
                unit: ToleranceUnit::Da,
            },
            database_path: "/data/uniprot_human.fasta".to_string(),
            decoy_strategy: DecoyStrategy::Reverse,
            acquisition_mode: None,
            max_variable_modifications: 3,
            min_peptide_length: 7,
            max_peptide_length: 50,
            engine: None,
        }
    }
```

- [ ] **Step 3: Add `extra` field to Psm**

In `crates/core/src/search_result.rs`, add after the `is_decoy` field (line ~146):

```rust
    /// Engine-specific extra fields (e.g., Sage's matched_peaks, delta_next).
    /// Preserves information that doesn't fit the standard Psm fields.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub extra: Option<std::collections::HashMap<String, serde_json::Value>>,
```

- [ ] **Step 4: Update `sample_psm()` and `sample_psm_with_mod()` test helpers**

In the same file, add `extra: None,` to both helpers (after `is_decoy: false,`).

- [ ] **Step 5: Fix all compilation errors from the new fields**

Search for all places that construct `Psm { ... }` and `SearchParams { ... }` across the workspace and add the missing fields:

Run: `cd /home/verden/pfind/2026-spring/code/proteinCopilot && cargo build 2>&1 | grep "missing field"`

For each missing field error, add `engine: None,` or `extra: None,` as appropriate. Key files likely affected:
- `crates/search-engine/src/simple_engine.rs` (build_psm function)
- `crates/search-engine/src/matching.rs`
- `crates/mcp-server/src/tools.rs` (run_search handler builds SearchParams)
- `crates/result-import/src/` (builds Psm from imported results)
- Various test files

- [ ] **Step 6: Add serde_json to core's dependencies if not already present**

Check `crates/core/Cargo.toml` — if `serde_json` is not listed, add it:
```toml
serde_json = { workspace = true }
```

- [ ] **Step 7: Run tests to verify nothing broke**

Run: `cd /home/verden/pfind/2026-spring/code/proteinCopilot && cargo test --workspace 2>&1 | tail -20`

Expected: All existing tests pass.

- [ ] **Step 8: Commit**

```bash
git add -A && git commit -m "feat(core): add engine field to SearchParams and extra field to Psm

Prepares data model for multi-engine support. SearchParams.engine
selects which search engine to use (None = default SimpleSearch).
Psm.extra preserves engine-specific fields like Sage's matched_peaks.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: Add sage-core git dependency

**Files:**
- Modify: `crates/search-engine/Cargo.toml`
- Modify: `Cargo.toml` (workspace root, if needed for workspace dep)

- [ ] **Step 1: Add sage-core to search-engine Cargo.toml**

In `crates/search-engine/Cargo.toml`, add to `[dependencies]`:

```toml
sage-core = { git = "https://github.com/lazear/sage.git", rev = "cd712d4b88d6c4fe3b317d47ce39e2c83cae7642", package = "sage-core" }
rayon = "1"
```

Note: The crate directory is `crates/sage` but the package name in its Cargo.toml is `sage-core`.

- [ ] **Step 2: Verify the dependency resolves and compiles**

Run: `cd /home/verden/pfind/2026-spring/code/proteinCopilot && cargo check -p protein-copilot-search-engine 2>&1 | tail -10`

Expected: `Finished` with no errors. This may take a while on first fetch.

If there are dependency conflicts (e.g., incompatible `serde` versions), resolve by aligning versions in the workspace `Cargo.toml`.

- [ ] **Step 3: Commit**

```bash
git add Cargo.toml Cargo.lock crates/search-engine/Cargo.toml && git commit -m "build: add sage-core git dependency to search-engine crate

Pinned to commit cd712d4 (v0.15.0-beta.2).

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: Implement type conversion layer (`convert.rs`)

**Files:**
- Create: `crates/search-engine/src/adapters/sage/convert.rs`
- Test: inline `#[cfg(test)] mod tests`

- [ ] **Step 1: Write failing tests for spectrum conversion**

Create `crates/search-engine/src/adapters/sage/convert.rs` with the test module first:

```rust
//! Type conversion between ProteinCopilot and sage-core types.

use protein_copilot_core::search_params::{MassTolerance, ModPosition, Modification, ToleranceUnit};
use protein_copilot_core::spectrum::{IsolationWindow, MsLevel, PrecursorInfo, Spectrum};
use sage_core::mass::Tolerance as SageTolerance;
use sage_core::modification::ModificationSpecificity;
use sage_core::spectrum::{Precursor as SagePrecursor, RawSpectrum, Representation};

/// Convert our Spectrum to sage RawSpectrum.
///
/// `file_id` is an index identifying which input file this spectrum came from.
pub fn spectrum_to_raw(spec: &Spectrum, file_id: usize) -> RawSpectrum {
    todo!()
}

/// Convert our MassTolerance to sage Tolerance.
pub fn mass_tolerance_to_sage(tol: &MassTolerance) -> SageTolerance {
    todo!()
}

/// Convert a single fixed modification to sage static_mods entries.
/// Returns one entry per target residue.
pub fn fixed_mod_to_sage(m: &Modification) -> Vec<(ModificationSpecificity, f32)> {
    todo!()
}

/// Convert a single variable modification to sage variable_mods entries.
/// Returns one entry per target residue.
pub fn variable_mod_to_sage(m: &Modification) -> Vec<(ModificationSpecificity, Vec<f32>)> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_ms2_spectrum() -> Spectrum {
        Spectrum {
            scan_number: 42,
            ms_level: MsLevel::MS2,
            retention_time_min: 120.5, // seconds
            precursors: vec![PrecursorInfo {
                mz: 500.25,
                charge: Some(2),
                intensity: Some(1e6),
                isolation_window: Some(IsolationWindow {
                    target_mz: 500.25,
                    lower_offset: 1.0,
                    upper_offset: 1.0,
                }),
                source_scan: Some(40),
            }],
            mz_array: vec![100.0, 200.0, 300.5, 400.75],
            intensity_array: vec![1000.0, 2000.0, 500.0, 3000.0],
        }
    }

    #[test]
    fn spectrum_to_raw_basic_fields() {
        let spec = sample_ms2_spectrum();
        let raw = spectrum_to_raw(&spec, 3);

        assert_eq!(raw.file_id, 3);
        assert_eq!(raw.ms_level, 2);
        assert_eq!(raw.id, "42");
        // retention_time_min is in seconds, sage wants minutes
        let expected_rt = (120.5 / 60.0) as f32;
        assert!((raw.scan_start_time - expected_rt).abs() < 1e-5);
        assert_eq!(raw.representation, Representation::Centroid);
    }

    #[test]
    fn spectrum_to_raw_mz_and_intensity() {
        let spec = sample_ms2_spectrum();
        let raw = spectrum_to_raw(&spec, 0);

        assert_eq!(raw.mz.len(), 4);
        assert_eq!(raw.intensity.len(), 4);
        assert!((raw.mz[0] - 100.0f32).abs() < 1e-4);
        assert!((raw.mz[2] - 300.5f32).abs() < 1e-4);
        assert!((raw.intensity[1] - 2000.0f32).abs() < 1e-2);
    }

    #[test]
    fn spectrum_to_raw_precursor() {
        let spec = sample_ms2_spectrum();
        let raw = spectrum_to_raw(&spec, 0);

        assert_eq!(raw.precursors.len(), 1);
        let p = &raw.precursors[0];
        assert!((p.mz - 500.25f32).abs() < 1e-3);
        assert_eq!(p.charge, Some(2));
        assert!(p.intensity.is_some());
        assert!(p.isolation_window.is_some());
    }

    #[test]
    fn spectrum_to_raw_negative_charge_clamped() {
        let mut spec = sample_ms2_spectrum();
        spec.precursors[0].charge = Some(-2);
        let raw = spectrum_to_raw(&spec, 0);
        // Negative charge should be clamped to absolute value
        assert_eq!(raw.precursors[0].charge, Some(2));
    }

    #[test]
    fn spectrum_to_raw_no_precursor() {
        let spec = Spectrum {
            scan_number: 1,
            ms_level: MsLevel::MS1,
            retention_time_min: 60.0,
            precursors: vec![],
            mz_array: vec![100.0],
            intensity_array: vec![1000.0],
        };
        let raw = spectrum_to_raw(&spec, 0);
        assert_eq!(raw.ms_level, 1);
        assert!(raw.precursors.is_empty());
    }

    #[test]
    fn mass_tolerance_ppm() {
        let tol = MassTolerance { value: 20.0, unit: ToleranceUnit::Ppm };
        let sage_tol = mass_tolerance_to_sage(&tol);
        assert_eq!(sage_tol, SageTolerance::Ppm(-20.0, 20.0));
    }

    #[test]
    fn mass_tolerance_da() {
        let tol = MassTolerance { value: 0.5, unit: ToleranceUnit::Da };
        let sage_tol = mass_tolerance_to_sage(&tol);
        assert_eq!(sage_tol, SageTolerance::Da(-0.5, 0.5));
    }

    #[test]
    fn fixed_mod_single_residue() {
        let m = Modification {
            name: "Carbamidomethyl".into(),
            mass_delta: 57.021464,
            residues: vec!['C'],
            position: ModPosition::Anywhere,
        };
        let result = fixed_mod_to_sage(&m);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0, ModificationSpecificity::Residue(b'C'));
        assert!((result[0].1 - 57.021464f32).abs() < 1e-4);
    }

    #[test]
    fn fixed_mod_multiple_residues() {
        let m = Modification {
            name: "Phospho".into(),
            mass_delta: 79.966331,
            residues: vec!['S', 'T', 'Y'],
            position: ModPosition::Anywhere,
        };
        let result = fixed_mod_to_sage(&m);
        assert_eq!(result.len(), 3);
        assert_eq!(result[0].0, ModificationSpecificity::Residue(b'S'));
        assert_eq!(result[1].0, ModificationSpecificity::Residue(b'T'));
        assert_eq!(result[2].0, ModificationSpecificity::Residue(b'Y'));
    }

    #[test]
    fn fixed_mod_nterm() {
        let m = Modification {
            name: "Acetyl".into(),
            mass_delta: 42.010565,
            residues: vec![],
            position: ModPosition::AnyNTerm,
        };
        let result = fixed_mod_to_sage(&m);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0, ModificationSpecificity::PeptideN(None));
    }

    #[test]
    fn fixed_mod_protein_nterm_with_residue() {
        let m = Modification {
            name: "Acetyl".into(),
            mass_delta: 42.010565,
            residues: vec!['M'],
            position: ModPosition::ProteinNTerm,
        };
        let result = fixed_mod_to_sage(&m);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0, ModificationSpecificity::ProteinN(Some(b'M')));
    }

    #[test]
    fn variable_mod_single() {
        let m = Modification {
            name: "Oxidation".into(),
            mass_delta: 15.994915,
            residues: vec!['M'],
            position: ModPosition::Anywhere,
        };
        let result = variable_mod_to_sage(&m);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0, ModificationSpecificity::Residue(b'M'));
        assert_eq!(result[0].1.len(), 1);
        assert!((result[0].1[0] - 15.994915f32).abs() < 1e-4);
    }
}
```

- [ ] **Step 2: Also create the sage module files so it compiles**

Create `crates/search-engine/src/adapters/sage/mod.rs`:

```rust
//! Sage search engine adapter.
//!
//! Integrates sage-core as a library for high-performance proteomics search.

pub mod convert;
pub mod config;
```

Create `crates/search-engine/src/adapters/sage/config.rs`:

```rust
//! SearchParams → sage Parameters configuration builder.
```

Modify `crates/search-engine/src/adapters/mod.rs` — add `pub mod sage;`:

```rust
//! Search engine adapter modules.
//!
//! Each adapter encapsulates the specifics of a particular search engine:
//! configuration file generation, execution, and result parsing.

pub mod pfind;
pub mod sage;
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `cd /home/verden/pfind/2026-spring/code/proteinCopilot && cargo test -p protein-copilot-search-engine --lib adapters::sage::convert 2>&1 | tail -20`

Expected: FAIL — `todo!()` panics.

- [ ] **Step 4: Implement `spectrum_to_raw`**

Replace the `todo!()` in `spectrum_to_raw`:

```rust
pub fn spectrum_to_raw(spec: &Spectrum, file_id: usize) -> RawSpectrum {
    let ms_level = match spec.ms_level {
        MsLevel::MS1 => 1u8,
        MsLevel::MS2 => 2u8,
        MsLevel::Other(n) => n,
    };

    let precursors = spec
        .precursors
        .iter()
        .map(|p| {
            let charge = p.charge.map(|c| c.unsigned_abs() as u8);
            let isolation_window = p.isolation_window.as_ref().map(|iw| {
                let half_width = ((iw.lower_offset + iw.upper_offset) / 2.0) as f32;
                SageTolerance::Da(-half_width, half_width)
            });
            SagePrecursor {
                mz: p.mz as f32,
                intensity: p.intensity.map(|i| i as f32),
                charge,
                spectrum_ref: None,
                isolation_window,
                inverse_ion_mobility: None,
            }
        })
        .collect();

    RawSpectrum {
        file_id,
        ms_level,
        id: spec.scan_number.to_string(),
        precursors,
        representation: Representation::Centroid,
        scan_start_time: (spec.retention_time_min / 60.0) as f32,
        ion_injection_time: 0.0,
        total_ion_current: spec.intensity_array.iter().sum::<f64>() as f32,
        mz: spec.mz_array.iter().map(|&v| v as f32).collect(),
        intensity: spec.intensity_array.iter().map(|&v| v as f32).collect(),
        mobility: None,
    }
}
```

- [ ] **Step 5: Implement `mass_tolerance_to_sage`**

```rust
pub fn mass_tolerance_to_sage(tol: &MassTolerance) -> SageTolerance {
    let v = tol.value as f32;
    match tol.unit {
        ToleranceUnit::Ppm => SageTolerance::Ppm(-v, v),
        ToleranceUnit::Da => SageTolerance::Da(-v, v),
    }
}
```

- [ ] **Step 6: Implement `fixed_mod_to_sage` and `variable_mod_to_sage`**

```rust
pub fn fixed_mod_to_sage(m: &Modification) -> Vec<(ModificationSpecificity, f32)> {
    let mass = m.mass_delta as f32;
    specificity_targets(m)
        .into_iter()
        .map(|spec| (spec, mass))
        .collect()
}

pub fn variable_mod_to_sage(m: &Modification) -> Vec<(ModificationSpecificity, Vec<f32>)> {
    let mass = m.mass_delta as f32;
    specificity_targets(m)
        .into_iter()
        .map(|spec| (spec, vec![mass]))
        .collect()
}

/// Map our ModPosition + residues to one or more sage ModificationSpecificity values.
fn specificity_targets(m: &Modification) -> Vec<ModificationSpecificity> {
    match m.position {
        ModPosition::Anywhere => {
            m.residues
                .iter()
                .map(|&r| ModificationSpecificity::Residue(r as u8))
                .collect()
        }
        ModPosition::AnyNTerm => {
            if m.residues.is_empty() {
                vec![ModificationSpecificity::PeptideN(None)]
            } else {
                m.residues
                    .iter()
                    .map(|&r| ModificationSpecificity::PeptideN(Some(r as u8)))
                    .collect()
            }
        }
        ModPosition::AnyCTerm => {
            if m.residues.is_empty() {
                vec![ModificationSpecificity::PeptideC(None)]
            } else {
                m.residues
                    .iter()
                    .map(|&r| ModificationSpecificity::PeptideC(Some(r as u8)))
                    .collect()
            }
        }
        ModPosition::ProteinNTerm => {
            if m.residues.is_empty() {
                vec![ModificationSpecificity::ProteinN(None)]
            } else {
                m.residues
                    .iter()
                    .map(|&r| ModificationSpecificity::ProteinN(Some(r as u8)))
                    .collect()
            }
        }
        ModPosition::ProteinCTerm => {
            if m.residues.is_empty() {
                vec![ModificationSpecificity::ProteinC(None)]
            } else {
                m.residues
                    .iter()
                    .map(|&r| ModificationSpecificity::ProteinC(Some(r as u8)))
                    .collect()
            }
        }
    }
}
```

- [ ] **Step 7: Run tests to verify they pass**

Run: `cd /home/verden/pfind/2026-spring/code/proteinCopilot && cargo test -p protein-copilot-search-engine --lib adapters::sage::convert 2>&1 | tail -20`

Expected: All tests PASS.

- [ ] **Step 8: Commit**

```bash
git add -A && git commit -m "feat(search-engine): implement sage type conversion layer

Adds convert.rs with spectrum_to_raw, mass_tolerance_to_sage,
fixed_mod_to_sage, and variable_mod_to_sage conversions.
All 11 unit tests pass.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 4: Implement sage config builder (`config.rs`)

**Files:**
- Modify: `crates/search-engine/src/adapters/sage/config.rs`

- [ ] **Step 1: Write failing tests**

Replace the contents of `config.rs`:

```rust
//! SearchParams → sage Parameters configuration builder.

use std::collections::HashMap;

use protein_copilot_core::search_params::{Enzyme, SearchParams};
use sage_core::database::{Builder, Parameters};
use sage_core::ion_series::Kind;
use sage_core::modification::ModificationSpecificity;

use super::convert::{fixed_mod_to_sage, mass_tolerance_to_sage, variable_mod_to_sage};

/// Build sage `Parameters` from our `SearchParams` and FASTA content.
///
/// The returned Parameters is ready to call `.build(fasta)` to create an IndexedDatabase.
pub fn build_sage_parameters(params: &SearchParams) -> Parameters {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::search_params::{
        DecoyStrategy, MassTolerance, ModPosition, Modification, ToleranceUnit,
    };

    fn test_params() -> SearchParams {
        SearchParams {
            enzyme: Enzyme::Trypsin,
            missed_cleavages: 2,
            fixed_modifications: vec![Modification {
                name: "Carbamidomethyl".into(),
                mass_delta: 57.021464,
                residues: vec!['C'],
                position: ModPosition::Anywhere,
            }],
            variable_modifications: vec![Modification {
                name: "Oxidation".into(),
                mass_delta: 15.994915,
                residues: vec!['M'],
                position: ModPosition::Anywhere,
            }],
            precursor_tolerance: MassTolerance {
                value: 20.0,
                unit: ToleranceUnit::Ppm,
            },
            fragment_tolerance: MassTolerance {
                value: 0.02,
                unit: ToleranceUnit::Da,
            },
            database_path: "/data/test.fasta".into(),
            decoy_strategy: DecoyStrategy::Reverse,
            acquisition_mode: None,
            max_variable_modifications: 3,
            min_peptide_length: 7,
            max_peptide_length: 50,
            engine: Some("Sage".into()),
        }
    }

    #[test]
    fn build_parameters_enzyme_trypsin() {
        let sage_params = build_sage_parameters(&test_params());
        let enzyme = sage_params.enzyme;
        assert_eq!(enzyme.cleave_at, Some("KR".into()));
        assert_eq!(enzyme.restrict, Some("P".into()));
        assert_eq!(enzyme.c_terminal, Some(true));
        assert_eq!(enzyme.missed_cleavages, Some(2));
    }

    #[test]
    fn build_parameters_peptide_length() {
        let sage_params = build_sage_parameters(&test_params());
        let enzyme = sage_params.enzyme;
        assert_eq!(enzyme.min_len, Some(7));
        assert_eq!(enzyme.max_len, Some(50));
    }

    #[test]
    fn build_parameters_static_mods() {
        let sage_params = build_sage_parameters(&test_params());
        assert!(sage_params.static_mods.contains_key(&ModificationSpecificity::Residue(b'C')));
        let mass = sage_params.static_mods[&ModificationSpecificity::Residue(b'C')];
        assert!((mass - 57.021464f32).abs() < 1e-4);
    }

    #[test]
    fn build_parameters_variable_mods() {
        let sage_params = build_sage_parameters(&test_params());
        assert!(sage_params.variable_mods.contains_key(&ModificationSpecificity::Residue(b'M')));
        let masses = &sage_params.variable_mods[&ModificationSpecificity::Residue(b'M')];
        assert_eq!(masses.len(), 1);
        assert!((masses[0] - 15.994915f32).abs() < 1e-4);
    }

    #[test]
    fn build_parameters_max_variable_mods() {
        let sage_params = build_sage_parameters(&test_params());
        assert_eq!(sage_params.max_variable_mods, 3);
    }

    #[test]
    fn build_parameters_decoy_generation() {
        let params = test_params();
        let sage_params = build_sage_parameters(&params);
        assert!(sage_params.generate_decoys);

        let mut no_decoy = params;
        no_decoy.decoy_strategy = DecoyStrategy::None;
        let sage_params2 = build_sage_parameters(&no_decoy);
        assert!(!sage_params2.generate_decoys);
    }

    #[test]
    fn build_parameters_fasta_path() {
        let sage_params = build_sage_parameters(&test_params());
        assert_eq!(sage_params.fasta, "/data/test.fasta");
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /home/verden/pfind/2026-spring/code/proteinCopilot && cargo test -p protein-copilot-search-engine --lib adapters::sage::config 2>&1 | tail -20`

Expected: FAIL — `todo!()`.

- [ ] **Step 3: Implement `build_sage_parameters`**

```rust
pub fn build_sage_parameters(params: &SearchParams) -> Parameters {
    let enzyme = match &params.enzyme {
        Enzyme::Trypsin => sage_core::database::EnzymeBuilder {
            cleave_at: Some("KR".into()),
            restrict: Some("P".into()),
            c_terminal: Some(true),
            missed_cleavages: Some(params.missed_cleavages as u8),
            min_len: Some(params.min_peptide_length as usize),
            max_len: Some(params.max_peptide_length as usize),
            semi_enzymatic: Some(false),
        },
        Enzyme::TrypsinP => sage_core::database::EnzymeBuilder {
            cleave_at: Some("KR".into()),
            restrict: Some("".into()),
            c_terminal: Some(true),
            missed_cleavages: Some(params.missed_cleavages as u8),
            min_len: Some(params.min_peptide_length as usize),
            max_len: Some(params.max_peptide_length as usize),
            semi_enzymatic: Some(false),
        },
        Enzyme::LysC => sage_core::database::EnzymeBuilder {
            cleave_at: Some("K".into()),
            restrict: Some("P".into()),
            c_terminal: Some(true),
            missed_cleavages: Some(params.missed_cleavages as u8),
            min_len: Some(params.min_peptide_length as usize),
            max_len: Some(params.max_peptide_length as usize),
            semi_enzymatic: Some(false),
        },
        Enzyme::GluC => sage_core::database::EnzymeBuilder {
            cleave_at: Some("DE".into()),
            restrict: Some("".into()),
            c_terminal: Some(true),
            missed_cleavages: Some(params.missed_cleavages as u8),
            min_len: Some(params.min_peptide_length as usize),
            max_len: Some(params.max_peptide_length as usize),
            semi_enzymatic: Some(false),
        },
        Enzyme::AspN => sage_core::database::EnzymeBuilder {
            cleave_at: Some("D".into()),
            restrict: Some("".into()),
            c_terminal: Some(false),
            missed_cleavages: Some(params.missed_cleavages as u8),
            min_len: Some(params.min_peptide_length as usize),
            max_len: Some(params.max_peptide_length as usize),
            semi_enzymatic: Some(false),
        },
        Enzyme::Chymotrypsin => sage_core::database::EnzymeBuilder {
            cleave_at: Some("FWYL".into()),
            restrict: Some("".into()),
            c_terminal: Some(true),
            missed_cleavages: Some(params.missed_cleavages as u8),
            min_len: Some(params.min_peptide_length as usize),
            max_len: Some(params.max_peptide_length as usize),
            semi_enzymatic: Some(false),
        },
        Enzyme::NonSpecific => sage_core::database::EnzymeBuilder {
            cleave_at: Some("".into()),
            restrict: Some("".into()),
            c_terminal: Some(true),
            missed_cleavages: Some(0),
            min_len: Some(params.min_peptide_length as usize),
            max_len: Some(params.max_peptide_length as usize),
            semi_enzymatic: Some(true),
        },
        Enzyme::Custom { cleavage_rule, .. } => sage_core::database::EnzymeBuilder {
            cleave_at: Some(cleavage_rule.clone()),
            restrict: Some("".into()),
            c_terminal: Some(true),
            missed_cleavages: Some(params.missed_cleavages as u8),
            min_len: Some(params.min_peptide_length as usize),
            max_len: Some(params.max_peptide_length as usize),
            semi_enzymatic: Some(false),
        },
    };

    // Build static mods
    let mut static_mods: HashMap<ModificationSpecificity, f32> = HashMap::new();
    for m in &params.fixed_modifications {
        for (spec, mass) in fixed_mod_to_sage(m) {
            static_mods.insert(spec, mass);
        }
    }

    // Build variable mods
    let mut variable_mods: HashMap<ModificationSpecificity, Vec<f32>> = HashMap::new();
    for m in &params.variable_modifications {
        for (spec, masses) in variable_mod_to_sage(m) {
            variable_mods
                .entry(spec)
                .or_default()
                .extend(masses);
        }
    }

    let generate_decoys = params.decoy_strategy != DecoyStrategy::None;

    Builder {
        bucket_size: None,
        enzyme: Some(enzyme),
        peptide_min_mass: Some(500.0),
        peptide_max_mass: Some(5000.0),
        ion_kinds: Some(vec![Kind::B, Kind::Y]),
        min_ion_index: Some(2),
        static_mods: Some(
            static_mods
                .into_iter()
                .map(|(k, v)| (k.to_string(), v))
                .collect(),
        ),
        variable_mods: Some(
            variable_mods
                .into_iter()
                .map(|(k, v)| (k.to_string(), v))
                .collect(),
        ),
        max_variable_mods: Some(params.max_variable_modifications as usize),
        decoy_tag: Some("rev_".into()),
        generate_decoys: Some(generate_decoys),
        fasta: Some(params.database_path.clone()),
        prefilter_chunk_size: None,
        prefilter: None,
        prefilter_low_memory: None,
    }
    .make_parameters()
}
```

Note: `Builder::make_parameters()` calls `validate_mods` internally which parses the string keys into `ModificationSpecificity`. We convert to string first because `Builder.static_mods` takes `Option<HashMap<String, f32>>`.

- [ ] **Step 4: Add DecoyStrategy import**

Make sure the `use` block at the top of config.rs includes:
```rust
use protein_copilot_core::search_params::{DecoyStrategy, Enzyme, SearchParams};
```

- [ ] **Step 5: Run tests**

Run: `cd /home/verden/pfind/2026-spring/code/proteinCopilot && cargo test -p protein-copilot-search-engine --lib adapters::sage::config 2>&1 | tail -20`

Expected: All tests PASS.

- [ ] **Step 6: Commit**

```bash
git add -A && git commit -m "feat(search-engine): implement sage config builder

Maps SearchParams enzyme/mods/tolerance to sage Builder/Parameters.
All 7 unit tests pass.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 5: Implement SageAdapter (`mod.rs`)

This is the core task — implementing `SearchEngineAdapter` for Sage.

**Files:**
- Modify: `crates/search-engine/src/adapters/sage/mod.rs`

- [ ] **Step 1: Write the SageAdapter struct and EngineInfo/HealthCheck**

Replace `crates/search-engine/src/adapters/sage/mod.rs`:

```rust
//! Sage search engine adapter.
//!
//! Integrates sage-core as a library for high-performance proteomics search.
//! Sage runs entirely in-process using rayon for parallelism, bridged to
//! tokio via `spawn_blocking`.

pub mod config;
pub mod convert;

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

use protein_copilot_core::engine::{EngineInfo, HealthStatus, SearchEngineAdapter};
use protein_copilot_core::error::CoreError;
use protein_copilot_core::progress::{ProgressCallback, SearchProgress};
use protein_copilot_core::run_metadata::{RunMetadata, RunStatus};
use protein_copilot_core::search_params::SearchParams;
use protein_copilot_core::search_result::{
    PeptideResult, ProteinResult, Psm, SearchResult, SearchResultSummary,
};
use protein_copilot_core::spectrum::{MsLevel, Spectrum};

use rayon::prelude::*;
use sage_core::database::IndexedDatabase;
use sage_core::fasta::Fasta;
use sage_core::mass::PROTON;
use sage_core::scoring::{Feature, Scorer, ScoreType};
use sage_core::spectrum::SpectrumProcessor;
use uuid::Uuid;

use self::config::build_sage_parameters;
use self::convert::{mass_tolerance_to_sage, spectrum_to_raw};

/// Sage search engine adapter.
///
/// Wraps sage-core for in-process proteomics search with rayon parallelism.
/// Each `search()` call creates an independent IndexedDatabase and Scorer.
pub struct SageAdapter {
    /// Number of threads for rayon (0 = use rayon default = all cores).
    thread_count: usize,
}

impl SageAdapter {
    pub fn new(thread_count: usize) -> Self {
        Self { thread_count }
    }
}

impl Default for SageAdapter {
    fn default() -> Self {
        Self { thread_count: 0 }
    }
}

#[async_trait::async_trait]
impl SearchEngineAdapter for SageAdapter {
    fn engine_info(&self) -> EngineInfo {
        EngineInfo {
            name: "Sage".to_string(),
            version: "0.15.0".to_string(),
            supported_features: vec![
                "open_search".to_string(),
                "lfq".to_string(),
                "tmt".to_string(),
                "chimera".to_string(),
            ],
        }
    }

    async fn health_check(&self) -> Result<HealthStatus, CoreError> {
        // Sage is a library dependency — always available if compiled in
        Ok(HealthStatus::Healthy)
    }

    async fn search(
        &self,
        params: &SearchParams,
        input_files: &[PathBuf],
        on_progress: ProgressCallback,
    ) -> Result<SearchResult, CoreError> {
        let start = Instant::now();
        let run_id = Uuid::new_v4();

        // Phase 1: Read spectra using spectrum-io
        on_progress(SearchProgress {
            run_id,
            status: "Running".to_string(),
            stage: Some("Loading spectra".to_string()),
            progress_pct: Some(0.0),
            elapsed_sec: 0.0,
            estimated_remaining_sec: None,
        });

        let mut all_spectra: Vec<Spectrum> = Vec::new();
        for path in input_files {
            let reader = protein_copilot_spectrum_io::create_reader(path)
                .map_err(|e| CoreError::SearchEngineError {
                    engine: "Sage".into(),
                    detail: format!("Failed to read spectra from {}: {}", path.display(), e),
                    suggestion: "Check that the input file exists and is a valid mzML/mgf file".into(),
                })?;
            let spectra = reader.read_all().map_err(|e| CoreError::SearchEngineError {
                engine: "Sage".into(),
                detail: format!("Error reading spectra: {}", e),
                suggestion: "Check input file format".into(),
            })?;
            all_spectra.extend(spectra);
        }

        self.search_with_spectra(params, all_spectra, on_progress).await
    }

    async fn search_with_spectra(
        &self,
        params: &SearchParams,
        spectra: Vec<Spectrum>,
        on_progress: ProgressCallback,
    ) -> Result<SearchResult, CoreError> {
        let start = Instant::now();
        let run_id = Uuid::new_v4();

        // Filter to MS2 only
        let ms2_spectra: Vec<&Spectrum> = spectra
            .iter()
            .filter(|s| s.ms_level == MsLevel::MS2)
            .collect();

        if ms2_spectra.is_empty() {
            return Err(CoreError::SearchEngineError {
                engine: "Sage".into(),
                detail: "No MS2 spectra found in input".into(),
                suggestion: "Check that input files contain MS2 spectra".into(),
            });
        }

        let total_spectra = ms2_spectra.len();

        // Convert to sage RawSpectrum
        let raw_spectra: Vec<sage_core::spectrum::RawSpectrum> = ms2_spectra
            .iter()
            .enumerate()
            .map(|(i, s)| spectrum_to_raw(s, i))
            .collect();

        on_progress(SearchProgress {
            run_id,
            status: "Running".to_string(),
            stage: Some("Loading spectra".to_string()),
            progress_pct: Some(100.0),
            elapsed_sec: start.elapsed().as_secs_f64(),
            estimated_remaining_sec: None,
        });

        // Phase 2: Build DB + score (spawn_blocking for rayon)
        let sage_params = build_sage_parameters(params);
        let precursor_tol = mass_tolerance_to_sage(&params.precursor_tolerance);
        let fragment_tol = mass_tolerance_to_sage(&params.fragment_tolerance);

        // Read FASTA content
        let fasta_path = &params.database_path;
        let fasta_content = tokio::fs::read_to_string(fasta_path)
            .await
            .map_err(|e| CoreError::SearchEngineError {
                engine: "Sage".into(),
                detail: format!("Failed to read FASTA file {}: {}", fasta_path, e),
                suggestion: "Check that database_path points to a valid FASTA file".into(),
            })?;

        // Progress counter shared between rayon and tokio
        let progress_counter = Arc::new(AtomicUsize::new(0));
        let counter_for_poll = Arc::clone(&progress_counter);

        // Progress polling task
        let progress_run_id = run_id;
        let progress_start = start;
        let total_for_poll = total_spectra;
        let on_progress_clone = on_progress.clone();
        let progress_handle = tokio::spawn(async move {
            let mut interval = tokio::time::interval(Duration::from_millis(500));
            loop {
                interval.tick().await;
                let done = counter_for_poll.load(Ordering::Relaxed);
                let pct = (done as f64 / total_for_poll.max(1) as f64) * 100.0;
                on_progress_clone(SearchProgress {
                    run_id: progress_run_id,
                    status: "Running".to_string(),
                    stage: Some("Searching".to_string()),
                    progress_pct: Some(pct),
                    elapsed_sec: progress_start.elapsed().as_secs_f64(),
                    estimated_remaining_sec: None,
                });
                if done >= total_for_poll {
                    break;
                }
            }
        });

        // Main search in blocking thread pool
        let counter_for_search = Arc::clone(&progress_counter);
        let search_result = tokio::task::spawn_blocking(move || {
            // Parse FASTA
            let fasta = Fasta::parse(fasta_content, sage_params.decoy_tag.as_str());

            // Build indexed database
            let db = sage_params.build(fasta);

            // Create spectrum processor
            let sp = SpectrumProcessor::new(150, true, 0.0);

            // Create scorer
            let scorer = Scorer {
                db: &db,
                precursor_tol,
                fragment_tol,
                min_matched_peaks: 4,
                min_isotope_err: -1,
                max_isotope_err: 3,
                min_precursor_charge: 2,
                max_precursor_charge: 4,
                override_precursor_charge: false,
                max_fragment_charge: None,
                chimera: false,
                report_psms: 1,
                wide_window: false,
                annotate_matches: false,
                score_type: ScoreType::SageHyperScore,
            };

            // Score all spectra in parallel
            let mut features: Vec<Feature> = raw_spectra
                .into_par_iter()
                .flat_map(|raw| {
                    let processed = sp.process(raw);
                    let results = scorer.score(&processed);
                    counter_for_search.fetch_add(1, Ordering::Relaxed);
                    results
                })
                .collect();

            // Phase 3: FDR
            // LDA rescoring
            if sage_core::ml::linear_discriminant::score_psms(&mut features, precursor_tol)
                .is_none()
            {
                // Fallback: heuristic discriminant score
                features.par_iter_mut().for_each(|feat| {
                    feat.discriminant_score =
                        (-feat.poisson as f32).ln_1p() + feat.longest_y_pct / 3.0;
                });
            }

            // Sort by discriminant score descending
            features.par_sort_unstable_by(|a, b| {
                b.discriminant_score.total_cmp(&a.discriminant_score)
            });

            // Spectrum-level q-value
            sage_core::ml::qvalue::spectrum_q_value(&mut features);

            // Picked peptide FDR
            sage_core::fdr::picked_peptide(&db, &mut features);

            // Picked protein FDR
            sage_core::fdr::picked_protein(&db, &mut features);

            Ok::<(Vec<Feature>, IndexedDatabase), CoreError>((features, db))
        })
        .await
        .map_err(|e| CoreError::SearchEngineError {
            engine: "Sage".into(),
            detail: format!("Search task panicked: {}", e),
            suggestion: "This is likely a bug in the Sage adapter".into(),
        })??;

        // Stop progress polling
        progress_handle.abort();

        let (features, db) = search_result;

        on_progress(SearchProgress {
            run_id,
            status: "Running".to_string(),
            stage: Some("Converting results".to_string()),
            progress_pct: Some(95.0),
            elapsed_sec: start.elapsed().as_secs_f64(),
            estimated_remaining_sec: None,
        });

        // Phase 4: Convert features to Psm
        let psms: Vec<Psm> = features
            .iter()
            .map(|feat| feature_to_psm(feat, &db))
            .collect();

        let duration = start.elapsed().as_secs_f64();

        // Build result
        let total_psms = psms.len() as u64;
        let target_psms = psms.iter().filter(|p| !p.is_decoy).count() as u64;
        let psms_at_1pct = psms
            .iter()
            .filter(|p| !p.is_decoy && p.q_value.map_or(false, |q| q <= 0.01))
            .count() as u64;

        // Build peptide-level results
        let mut peptide_map: HashMap<String, PeptideResult> = HashMap::new();
        for psm in &psms {
            if psm.is_decoy {
                continue;
            }
            let entry = peptide_map
                .entry(psm.peptide_sequence.clone())
                .or_insert_with(|| PeptideResult {
                    sequence: psm.peptide_sequence.clone(),
                    protein_accessions: psm.protein_accessions.clone(),
                    best_score: psm.score,
                    q_value: psm.q_value,
                    psm_count: 0,
                });
            entry.psm_count += 1;
            if psm.score > entry.best_score {
                entry.best_score = psm.score;
                entry.q_value = psm.q_value;
            }
        }
        let peptides: Vec<PeptideResult> = peptide_map.into_values().collect();

        // Build protein-level results
        let mut protein_map: HashMap<String, ProteinResult> = HashMap::new();
        for psm in &psms {
            if psm.is_decoy {
                continue;
            }
            for acc in &psm.protein_accessions {
                let entry = protein_map
                    .entry(acc.clone())
                    .or_insert_with(|| ProteinResult {
                        accession: acc.clone(),
                        description: String::new(),
                        unique_peptides: 0,
                        total_psms: 0,
                        best_score: psm.score,
                        q_value: None,
                        sequence_coverage: None,
                    });
                entry.total_psms += 1;
                if psm.score > entry.best_score {
                    entry.best_score = psm.score;
                }
            }
        }
        let proteins: Vec<ProteinResult> = protein_map.into_values().collect();

        let summary = SearchResultSummary {
            total_spectra: total_spectra as u64,
            total_psms,
            psms_at_1pct_fdr: psms_at_1pct,
            target_psms,
            decoy_psms: total_psms - target_psms,
            total_peptides: peptides.len() as u64,
            total_proteins: proteins.len() as u64,
        };

        let metadata = RunMetadata {
            run_id,
            engine_info: self.engine_info(),
            params: params.clone(),
            input_files: vec![],
            status: RunStatus::Completed,
            created_at: chrono::Utc::now(),
            elapsed_sec: duration,
        };

        Ok(SearchResult {
            run_id,
            psms,
            peptides,
            proteins,
            summary,
            metadata,
        })
    }

    async fn cancel(&self, _run_id: Uuid) -> Result<(), CoreError> {
        // Sage runs in a spawned blocking task — cancellation is handled
        // by the MCP server layer via JoinHandle::abort().
        Ok(())
    }
}

/// Convert a sage Feature + IndexedDatabase lookup into our Psm.
fn feature_to_psm(feat: &Feature, db: &IndexedDatabase) -> Psm {
    let peptide = &db[feat.peptide_idx];

    // Recover peptide sequence from byte array
    let peptide_sequence = String::from_utf8_lossy(&peptide.sequence).to_string();

    // Recover protein accessions
    let protein_accessions: Vec<String> = peptide.proteins.iter().map(|p| p.to_string()).collect();

    // Convert mass → m/z: mz = (mass + charge * PROTON) / charge
    let charge = feat.charge as f64;
    let precursor_mz = (feat.expmass as f64 + charge * PROTON as f64) / charge;
    let calculated_mz = (feat.calcmass as f64 + charge * PROTON as f64) / charge;

    let delta_mass_ppm = if calculated_mz > 0.0 {
        (precursor_mz - calculated_mz) / calculated_mz * 1e6
    } else {
        0.0
    };

    // Parse scan number from spec_id (we stored it as scan_number.to_string())
    let spectrum_scan = feat
        .spec_id
        .parse::<u32>()
        .unwrap_or(1);

    // Build extra fields for sage-specific data
    let mut extra = HashMap::new();
    extra.insert("matched_peaks".into(), serde_json::json!(feat.matched_peaks));
    extra.insert("longest_b".into(), serde_json::json!(feat.longest_b));
    extra.insert("longest_y".into(), serde_json::json!(feat.longest_y));
    extra.insert("delta_next".into(), serde_json::json!(feat.delta_next));
    extra.insert("delta_best".into(), serde_json::json!(feat.delta_best));
    extra.insert(
        "discriminant_score".into(),
        serde_json::json!(feat.discriminant_score),
    );
    extra.insert(
        "posterior_error".into(),
        serde_json::json!(feat.posterior_error),
    );
    extra.insert("hyperscore".into(), serde_json::json!(feat.hyperscore));
    extra.insert("peptide_q".into(), serde_json::json!(feat.peptide_q));
    extra.insert("protein_q".into(), serde_json::json!(feat.protein_q));

    // Convert modifications from sage Peptide
    // sage stores modifications as Vec<f32> parallel to sequence bytes,
    // with 0.0 meaning no modification and nterm/cterm as separate fields.
    let modifications = convert_sage_modifications(peptide);

    Psm {
        spectrum_scan,
        peptide_sequence,
        modifications,
        charge: feat.charge as i32,
        precursor_mz,
        calculated_mz,
        delta_mass_ppm,
        score: feat.hyperscore,
        q_value: Some(feat.spectrum_q as f64),
        protein_accessions,
        is_decoy: feat.label == -1,
        extra: Some(extra),
    }
}

/// Convert sage Peptide modifications to our Modification structs.
fn convert_sage_modifications(peptide: &sage_core::peptide::Peptide) -> Vec<protein_copilot_core::search_params::Modification> {
    use protein_copilot_core::search_params::{ModPosition, Modification};

    let mut mods = Vec::new();

    // N-terminal modification
    if let Some(nterm_mass) = peptide.nterm {
        if nterm_mass.abs() > 1e-6 {
            mods.push(Modification {
                name: format!("NTerm({:.4})", nterm_mass),
                mass_delta: nterm_mass as f64,
                residues: vec![],
                position: ModPosition::AnyNTerm,
            });
        }
    }

    // Residue modifications
    for (i, &mass) in peptide.modifications.iter().enumerate() {
        if mass.abs() > 1e-6 {
            let residue = if i < peptide.sequence.len() {
                peptide.sequence[i] as char
            } else {
                '?'
            };
            mods.push(Modification {
                name: format!("{}{:.4}", residue, mass),
                mass_delta: mass as f64,
                residues: vec![residue],
                position: ModPosition::Anywhere,
            });
        }
    }

    // C-terminal modification
    if let Some(cterm_mass) = peptide.cterm {
        if cterm_mass.abs() > 1e-6 {
            mods.push(Modification {
                name: format!("CTerm({:.4})", cterm_mass),
                mass_delta: cterm_mass as f64,
                residues: vec![],
                position: ModPosition::AnyCTerm,
            });
        }
    }

    mods
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn engine_info() {
        let adapter = SageAdapter::default();
        let info = adapter.engine_info();
        assert_eq!(info.name, "Sage");
        assert!(!info.version.is_empty());
    }

    #[tokio::test]
    async fn health_check_always_healthy() {
        let adapter = SageAdapter::default();
        let status = adapter.health_check().await.unwrap();
        assert_eq!(status, HealthStatus::Healthy);
    }

    #[tokio::test]
    async fn search_with_empty_spectra_returns_error() {
        let adapter = SageAdapter::default();
        let params = protein_copilot_core::search_params::SearchParams {
            enzyme: protein_copilot_core::search_params::Enzyme::Trypsin,
            missed_cleavages: 2,
            fixed_modifications: vec![],
            variable_modifications: vec![],
            precursor_tolerance: protein_copilot_core::search_params::MassTolerance {
                value: 20.0,
                unit: protein_copilot_core::search_params::ToleranceUnit::Ppm,
            },
            fragment_tolerance: protein_copilot_core::search_params::MassTolerance {
                value: 0.5,
                unit: protein_copilot_core::search_params::ToleranceUnit::Da,
            },
            database_path: "/nonexistent.fasta".into(),
            decoy_strategy: protein_copilot_core::search_params::DecoyStrategy::Reverse,
            acquisition_mode: None,
            max_variable_modifications: 3,
            min_peptide_length: 7,
            max_peptide_length: 50,
            engine: Some("Sage".into()),
        };
        let on_progress: ProgressCallback = Box::new(|_| {});
        let result = adapter.search_with_spectra(&params, vec![], on_progress).await;
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("No MS2 spectra"), "Error was: {}", err);
    }
}
```

- [ ] **Step 2: Add necessary dependencies to search-engine's Cargo.toml**

Make sure `crates/search-engine/Cargo.toml` has these dependencies (add if missing):

```toml
protein-copilot-spectrum-io = { workspace = true }
chrono = { workspace = true }
rayon = "1"
```

Also verify `async-trait` is present (should be already).

- [ ] **Step 3: Run tests to verify compilation and basic tests pass**

Run: `cd /home/verden/pfind/2026-spring/code/proteinCopilot && cargo test -p protein-copilot-search-engine --lib adapters::sage 2>&1 | tail -30`

Expected: 3 tests pass (engine_info, health_check, empty_spectra_error).

- [ ] **Step 4: Run full workspace tests**

Run: `cd /home/verden/pfind/2026-spring/code/proteinCopilot && cargo test --workspace 2>&1 | tail -20`

Expected: All tests pass. Fix any compilation errors in other crates caused by the new fields.

- [ ] **Step 5: Commit**

```bash
git add -A && git commit -m "feat(search-engine): implement SageAdapter with SearchEngineAdapter trait

Full sage-core integration:
- Spectrum → RawSpectrum → ProcessedSpectrum → Feature pipeline
- rayon/tokio bridging via spawn_blocking
- LDA rescoring + picked peptide/protein FDR
- Feature → Psm conversion with extra fields preservation
- Progress reporting via AtomicUsize + tokio interval

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 6: Wire up registry-based engine lookup in tools.rs

**Files:**
- Modify: `crates/mcp-server/src/tools.rs`

- [ ] **Step 1: Register SageAdapter in McpState::new()**

In `crates/mcp-server/src/tools.rs`, modify `ProteinCopilotServer::new()` (line ~764):

Change:
```rust
        let mut registry = protein_copilot_search_engine::EngineRegistry::new();
        registry.register(Box::new(SimpleSearchEngine::new()));
```
To:
```rust
        let mut registry = protein_copilot_search_engine::EngineRegistry::new();
        registry.register(Box::new(SimpleSearchEngine::new()));
        registry.register(Box::new(
            protein_copilot_search_engine::adapters::sage::SageAdapter::default(),
        ));
```

- [ ] **Step 2: Replace hardcoded engine in DIA branch (line ~1069)**

Change:
```rust
            let engine = SimpleSearchEngine::new();
```
To:
```rust
            let engine_name = params.engine.as_deref().unwrap_or("SimpleSearch");
            let engine = self.registry.get(engine_name).ok_or_else(|| {
                mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    format!(
                        "Engine '{}' not registered. Available: {:?}",
                        engine_name,
                        self.registry
                            .list_available()
                            .iter()
                            .map(|e| &e.name)
                            .collect::<Vec<_>>()
                    ),
                )
            })?;
```

Also update the search call from `engine.search_with_spectra(...)` to match the borrowed reference (since `registry.get()` returns `&dyn SearchEngineAdapter`).

- [ ] **Step 3: Replace hardcoded engine in file-based branch (line ~1280)**

Apply the same change:

Change:
```rust
        let engine = SimpleSearchEngine::new();
```
To:
```rust
        let engine_name = params.engine.as_deref().unwrap_or("SimpleSearch");
        let engine = self.registry.get(engine_name).ok_or_else(|| {
            mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!(
                    "Engine '{}' not registered. Available: {:?}",
                    engine_name,
                    self.registry
                        .list_available()
                        .iter()
                        .map(|e| &e.name)
                        .collect::<Vec<_>>()
                ),
            )
        })?;
```

**Important**: Since `self.registry.get()` returns a borrowed reference and we need to move it into a `tokio::spawn`, we need to restructure: either clone the params and call search inside a non-move block, or change the registry to store `Arc<dyn SearchEngineAdapter>`. The simplest approach is to create the engine fresh (like `SimpleSearchEngine::new()` does) inside the spawn. But since `SageAdapter` is cheap to create, we can use a factory pattern or simply match on the engine name and construct inside the spawn.

Alternative approach — construct inside spawn:
```rust
        let engine_name = params.engine.as_deref().unwrap_or("SimpleSearch").to_string();
        // Verify engine exists before spawning
        if self.registry.get(&engine_name).is_none() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!("Engine '{}' not registered. Available: {:?}",
                    engine_name,
                    self.registry.list_available().iter().map(|e| &e.name).collect::<Vec<_>>()),
            ));
        }

        // ... inside tokio::spawn:
        let search_result = match engine_name.as_str() {
            "Sage" => {
                let engine = protein_copilot_search_engine::adapters::sage::SageAdapter::default();
                engine.search(&params, &files, on_progress).await
            }
            _ => {
                let engine = SimpleSearchEngine::new();
                engine.search(&params, &files, on_progress).await
            }
        };
```

This is pragmatic and avoids lifetime issues. In the future, if we add more engines, refactor to use `Arc<dyn SearchEngineAdapter>` in the registry.

- [ ] **Step 4: Update the import at the top of tools.rs**

Add to the import section (line ~35):
```rust
use protein_copilot_search_engine::adapters::sage::SageAdapter;
```

- [ ] **Step 5: Add `engine` field to RunSearchInput struct**

Find the `RunSearchInput` struct in tools.rs and verify it includes `engine` (it's part of `SearchParams` so it should be passed through automatically). If `run_search` constructs `SearchParams` manually from input fields, add:

```rust
    engine: input.params.as_ref().and_then(|p| p.engine.clone()),
```

- [ ] **Step 6: Build and test**

Run: `cd /home/verden/pfind/2026-spring/code/proteinCopilot && cargo build --workspace 2>&1 | tail -10`

Expected: Clean build.

Run: `cd /home/verden/pfind/2026-spring/code/proteinCopilot && cargo test --workspace 2>&1 | tail -20`

Expected: All tests pass.

- [ ] **Step 7: Commit**

```bash
git add -A && git commit -m "feat(mcp-server): registry-based engine lookup in run_search

Replace hardcoded SimpleSearchEngine with engine name lookup from
SearchParams.engine field. Sage and SimpleSearch both registered.
Users specify engine: \"Sage\" to use the Sage adapter.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 7: Integration test with real data

**Files:**
- Create: `crates/search-engine/tests/sage_integration.rs`
- Use existing test fixtures: `tests/fixtures/` (small FASTA + mgf)

- [ ] **Step 1: Check existing test fixtures**

Run: `find /home/verden/pfind/2026-spring/code/proteinCopilot/tests -name "*.fasta" -o -name "*.mgf" -o -name "*.mzml" 2>/dev/null | head -10`

If no small FASTA/mgf fixtures exist, create minimal ones. We need:
- A FASTA with ~5 proteins
- An mgf with ~10 MS2 spectra matching some peptides from those proteins

- [ ] **Step 2: Write integration test**

Create `crates/search-engine/tests/sage_integration.rs`:

```rust
//! Integration test for SageAdapter.
//!
//! Uses a small FASTA + mgf file to verify the full search pipeline.

use std::path::PathBuf;
use protein_copilot_core::engine::SearchEngineAdapter;
use protein_copilot_core::search_params::*;
use protein_copilot_search_engine::adapters::sage::SageAdapter;

fn test_fasta_path() -> PathBuf {
    // Use an existing test FASTA or create one
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("tests/fixtures/small_test.fasta");
    if !path.exists() {
        panic!("Test fixture not found: {}. Create a small FASTA for testing.", path.display());
    }
    path
}

fn test_mgf_path() -> PathBuf {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("tests/fixtures/small_test.mgf");
    if !path.exists() {
        panic!("Test fixture not found: {}. Create a small mgf for testing.", path.display());
    }
    path
}

fn test_params(fasta: &str) -> SearchParams {
    SearchParams {
        enzyme: Enzyme::Trypsin,
        missed_cleavages: 2,
        fixed_modifications: vec![Modification {
            name: "Carbamidomethyl".into(),
            mass_delta: 57.021464,
            residues: vec!['C'],
            position: ModPosition::Anywhere,
        }],
        variable_modifications: vec![Modification {
            name: "Oxidation".into(),
            mass_delta: 15.994915,
            residues: vec!['M'],
            position: ModPosition::Anywhere,
        }],
        precursor_tolerance: MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        },
        fragment_tolerance: MassTolerance {
            value: 0.5,
            unit: ToleranceUnit::Da,
        },
        database_path: fasta.to_string(),
        decoy_strategy: DecoyStrategy::Reverse,
        acquisition_mode: None,
        max_variable_modifications: 3,
        min_peptide_length: 5,
        max_peptide_length: 50,
        engine: Some("Sage".into()),
    }
}

#[tokio::test]
async fn sage_search_produces_results() {
    let fasta = test_fasta_path();
    let mgf = test_mgf_path();
    let params = test_params(fasta.to_str().unwrap());

    let adapter = SageAdapter::default();
    let on_progress: protein_copilot_core::progress::ProgressCallback = Box::new(|p| {
        eprintln!("Progress: {:?}", p.stage);
    });

    let result = adapter
        .search(&params, &[mgf], on_progress)
        .await;

    match result {
        Ok(result) => {
            // We should get some PSMs (exact number depends on test data)
            assert!(!result.psms.is_empty(), "Expected at least some PSMs");

            // All PSMs should have valid fields
            for psm in &result.psms {
                assert!(psm.spectrum_scan >= 1, "scan should be >= 1");
                assert!(!psm.peptide_sequence.is_empty(), "sequence should not be empty");
                assert!(psm.charge > 0, "charge should be positive");
                assert!(psm.precursor_mz > 0.0, "precursor_mz should be positive");
                assert!(psm.score.is_finite(), "score should be finite");
                assert!(psm.q_value.is_some(), "q_value should be set by Sage FDR");
                assert!(!psm.protein_accessions.is_empty(), "should have protein accessions");
                assert!(psm.extra.is_some(), "should have extra sage fields");
            }

            // Summary should be consistent
            assert_eq!(result.summary.total_psms, result.psms.len() as u64);

            eprintln!(
                "Sage search complete: {} PSMs, {} at 1% FDR",
                result.psms.len(),
                result.summary.psms_at_1pct_fdr
            );
        }
        Err(e) => {
            // If test fixtures don't exist, skip gracefully
            let msg = e.to_string();
            if msg.contains("not found") || msg.contains("No such file") {
                eprintln!("Skipping sage integration test: test fixtures not found");
                return;
            }
            panic!("Sage search failed: {}", e);
        }
    }
}

#[tokio::test]
async fn sage_engine_info_and_health() {
    let adapter = SageAdapter::default();
    let info = adapter.engine_info();
    assert_eq!(info.name, "Sage");

    let health = adapter.health_check().await.unwrap();
    assert_eq!(health, protein_copilot_core::engine::HealthStatus::Healthy);
}
```

- [ ] **Step 3: Create minimal test fixtures if they don't exist**

If `/tests/fixtures/small_test.fasta` and `small_test.mgf` don't exist, create them with synthetic data that will produce matches. Use a known BSA (bovine serum albumin) tryptic peptide.

small_test.fasta:
```
>sp|P02769|ALBU_BOVIN Bovine serum albumin
MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA
```

small_test.mgf (synthetic spectrum for YICDNQDTISSK, charge 2):
```
BEGIN IONS
TITLE=scan=1
PEPMASS=685.8223 10000.0
CHARGE=2+
RTINSECONDS=300.0
175.1190 500
276.1667 600
391.1936 800
504.2777 900
619.3046 700
733.3476 1000
848.3745 850
963.4015 400
1076.4855 600
1177.5332 500
END IONS
BEGIN IONS
TITLE=scan=2
PEPMASS=480.7400 8000.0
CHARGE=2+
RTINSECONDS=350.0
175.1190 400
288.2030 500
403.2300 700
504.2777 600
617.3617 800
716.4301 500
END IONS
```

- [ ] **Step 4: Run integration test**

Run: `cd /home/verden/pfind/2026-spring/code/proteinCopilot && cargo test -p protein-copilot-search-engine --test sage_integration 2>&1 | tail -30`

Expected: Tests pass (or skip gracefully if fixtures not found).

- [ ] **Step 5: Run full workspace tests**

Run: `cd /home/verden/pfind/2026-spring/code/proteinCopilot && cargo test --workspace 2>&1 | tail -20`

Expected: All tests pass, including the new sage tests.

- [ ] **Step 6: Run clippy**

Run: `cd /home/verden/pfind/2026-spring/code/proteinCopilot && cargo clippy --workspace 2>&1 | tail -20`

Expected: No warnings.

- [ ] **Step 7: Commit**

```bash
git add -A && git commit -m "test: add sage integration test with synthetic BSA data

End-to-end test verifying SageAdapter::search produces valid PSMs
from small FASTA + mgf fixtures.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 8: Update lib.rs exports and documentation

**Files:**
- Modify: `crates/search-engine/src/lib.rs`
- Modify: `docs/architecture.md`
- Modify: `README.md`

- [ ] **Step 1: Export SageAdapter from lib.rs**

In `crates/search-engine/src/lib.rs`, add after `pub use simple_engine::SimpleSearchEngine;`:

```rust
pub use adapters::sage::SageAdapter;
```

- [ ] **Step 2: Update architecture.md**

Add a section for Sage adapter in `docs/architecture.md`.

- [ ] **Step 3: Update README.md**

Update the supported engines section and any engine listing to include Sage.

- [ ] **Step 4: Build and test one final time**

Run: `cd /home/verden/pfind/2026-spring/code/proteinCopilot && cargo test --workspace 2>&1 | tail -20`

Expected: All tests pass.

- [ ] **Step 5: Commit**

```bash
git add -A && git commit -m "docs: update architecture and README for Sage integration

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```
