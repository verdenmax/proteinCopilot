# DIA Phase 2 — Search Engine Integration & MCP Tool Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Integrate dia-extraction into the search pipeline and expose as MCP tool, enabling end-to-end DIA search.

**Architecture:** Modify `match_spectrum()` to iterate all precursors (not just first), add MS level filtering in the search loop, add `extract_dia_precursors` MCP tool that reads spectra → runs DIA extraction → caches enhanced spectra for `run_search`, and add `acquisition_mode` to `SearchParams`.

**Tech Stack:** Rust, rmcp (MCP framework), serde, schemars, tokio

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `crates/core/src/search_params.rs` | Modify | Add `acquisition_mode` field |
| `crates/search-engine/src/matching.rs` | Modify | Multi-precursor matching |
| `crates/search-engine/src/simple_engine.rs` | Modify | MS level filtering |
| `crates/mcp-server/Cargo.toml` | Modify | Add dia-extraction dependency |
| `crates/mcp-server/src/tools.rs` | Modify | Add `extract_dia_precursors` MCP tool |

---

### Task 1: Add `acquisition_mode` to SearchParams

**Files:**
- Modify: `crates/core/src/search_params.rs`

- [ ] **Step 1: Add field to SearchParams**

In `crates/core/src/search_params.rs`, add `acquisition_mode` to the `SearchParams` struct after the `decoy_strategy` field:

```rust
    /// Target-decoy strategy for FDR estimation.
    pub decoy_strategy: DecoyStrategy,
    /// Data acquisition mode. `None` = auto-detect or not applicable.
    #[serde(default)]
    pub acquisition_mode: Option<AcquisitionMode>,
```

Add the import at the top of the file alongside existing imports:

```rust
use crate::spectrum::AcquisitionMode;
```

- [ ] **Step 2: Fix all compilation errors from new field**

Search for all `SearchParams {` constructors in the workspace and add `acquisition_mode: None,` to each. Key locations:

Known locations (9+):
- `crates/core/src/search_params.rs` — `default()` and tests
- `crates/core/src/search_result.rs` — test helper
- `crates/core/src/ai_decision.rs` — test helper
- `crates/core/src/run_metadata.rs` — test helper
- `crates/param-recommend/src/preset.rs` — 5 preset constructors (standard, phospho, tmt, silac, open)
- `crates/param-recommend/src/rules.rs` — `recommend()` method
- `crates/search-engine/src/simple_engine.rs` — test helpers
- `crates/mcp-server/src/tools.rs` — any manual construction

Run: `cargo build 2>&1 | grep "missing field"` to catch any missed locations.

- [ ] **Step 3: Run tests**

```bash
cd /home/verden/pfind/2026-spring/code/proteinCopilot
cargo test --quiet
cargo clippy -- -D warnings
```

Expected: all tests pass, zero warnings.

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "feat(core): add acquisition_mode to SearchParams

Optional field with #[serde(default)] for backward compatibility.
Used by DIA extraction to communicate detected/overridden mode.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: Multi-precursor matching in search engine

**Files:**
- Modify: `crates/search-engine/src/matching.rs`

- [ ] **Step 1: Write tests for multi-precursor matching**

Add to the `#[cfg(test)] mod tests` in `matching.rs`:

```rust
#[test]
fn test_match_spectrum_multiple_precursors() {
    // Spectrum with two precursors at different m/z values
    let spectrum = Spectrum::new(
        1,
        MsLevel::MS2,
        100.0,
        vec![
            PrecursorInfo {
                mz: 500.0,
                charge: Some(2),
                intensity: None,
                isolation_window: None,
                source_scan: None,
            },
            PrecursorInfo {
                mz: 600.0,
                charge: Some(2),
                intensity: None,
                isolation_window: None,
                source_scan: None,
            },
        ],
        vec![175.119, 262.151, 276.134, 363.166],
        vec![100.0, 80.0, 90.0, 70.0],
    )
    .unwrap();

    // Create a peptide that matches the second precursor (m/z 600.0)
    // but not the first (m/z 500.0)
    let peptides = vec![DigestedPeptide {
        sequence: "PEPTIDE".to_string(),
        protein_accession: "test".to_string(),
        neutral_mass: 799.36, // (600.0 - 1.007276) * 2
        missed_cleavages: 0,
        start: 0,
        end: 7,
    }];

    let tol = MassTolerance {
        value: 20.0,
        unit: ToleranceUnit::Ppm,
    };
    let frag_tol = MassTolerance {
        value: 0.02,
        unit: ToleranceUnit::Da,
    };

    let results = match_spectrum_all(
        &spectrum,
        &peptides,
        &tol,
        &frag_tol,
        &[],
    );
    // Should find a match from the second precursor
    assert!(!results.is_empty());
}

#[test]
fn test_match_spectrum_all_empty_precursors() {
    let spectrum = Spectrum::new(
        1,
        MsLevel::MS2,
        100.0,
        vec![],
        vec![100.0, 200.0],
        vec![50.0, 50.0],
    )
    .unwrap();

    let results = match_spectrum_all(
        &spectrum,
        &[],
        &MassTolerance { value: 10.0, unit: ToleranceUnit::Ppm },
        &MassTolerance { value: 0.02, unit: ToleranceUnit::Da },
        &[],
    );
    assert!(results.is_empty());
}
```

- [ ] **Step 2: Implement `match_spectrum_all()`**

Add a new public function in `matching.rs` that iterates all precursors:

```rust
/// Matches a single spectrum against all candidate peptides, trying ALL precursors.
///
/// Unlike `match_spectrum()` which only uses the first precursor, this function
/// iterates over every precursor in `spectrum.precursors` and collects the best
/// match for each. Used for DIA data where a spectrum may have multiple
/// candidate precursor ions extracted from MS1.
///
/// Returns a Vec of matches (one per precursor that produced a hit).
/// Downstream FDR control handles quality filtering.
pub fn match_spectrum_all(
    spectrum: &Spectrum,
    candidates: &[DigestedPeptide],
    precursor_tolerance: &MassTolerance,
    fragment_tolerance: &MassTolerance,
    fixed_mods: &[Modification],
) -> Vec<PeptideMatch> {
    let mut all_matches = Vec::new();

    for precursor in &spectrum.precursors {
        let observed_mz = precursor.mz;

        let charge_states: Vec<i32> = if let Some(c) = precursor.charge {
            vec![c]
        } else {
            vec![2, 3, 1, 4]
        };

        let mut best_match: Option<PeptideMatch> = None;

        for peptide in candidates {
            let modified_mass =
                peptide.neutral_mass + apply_fixed_mods(&peptide.sequence, fixed_mods);

            for &charge in &charge_states {
                if charge == 0 {
                    continue;
                }
                let theoretical_mz = peptide_mz(modified_mass, charge);

                if within_tolerance(observed_mz, theoretical_mz, precursor_tolerance) {
                    let b_ions = generate_b_ions(&peptide.sequence);
                    let y_ions = generate_y_ions(&peptide.sequence);

                    let total_theoretical = (b_ions.len() + y_ions.len()) as u32;
                    if total_theoretical == 0 {
                        continue;
                    }

                    let all_ions: Vec<f64> = b_ions.into_iter().chain(y_ions).collect();
                    let matched =
                        count_matched_ions(&all_ions, &spectrum.mz_array, fragment_tolerance);

                    let score = matched as f64 / total_theoretical as f64;
                    let delta_ppm = calc_delta_ppm(observed_mz, theoretical_mz);

                    if !score.is_finite() || !delta_ppm.is_finite() {
                        continue;
                    }

                    let is_better = match &best_match {
                        None => true,
                        Some(prev) => score > prev.score,
                    };

                    if is_better {
                        best_match = Some(PeptideMatch {
                            peptide: peptide.clone(),
                            charge,
                            observed_mz,
                            theoretical_mz,
                            delta_mass_ppm: delta_ppm,
                            score,
                            matched_ions: matched,
                            total_ions: total_theoretical,
                        });
                    }
                }
            }
        }

        if let Some(m) = best_match {
            all_matches.push(m);
        }
    }

    all_matches
}
```

Note: Keep the original `match_spectrum()` function unchanged — it's used by existing DDA code and the annotation module.

- [ ] **Step 3: Re-export IsotopePatternExtractor**

In `crates/dia-extraction/src/lib.rs`, add to the re-exports section:

```rust
pub use isotope::IsotopePatternExtractor;
```

This is needed by Task 4 (MCP tool) which constructs `IsotopePatternExtractor::default()`.

- [ ] **Step 4: Run tests**

```bash
cd /home/verden/pfind/2026-spring/code/proteinCopilot
cargo test -p protein-copilot-search-engine --quiet
cargo clippy -p protein-copilot-search-engine -- -D warnings
```

Expected: all tests pass including the two new ones.

- [ ] **Step 5: Commit**

```bash
git add crates/search-engine/src/matching.rs crates/dia-extraction/src/lib.rs
git commit -m "feat(search-engine): add match_spectrum_all for multi-precursor DIA

Iterates all precursors in a spectrum, returning best match per precursor.
Original match_spectrum() retained for DDA backward compatibility.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: MS level filtering + DIA-aware search loop

**Files:**
- Modify: `crates/search-engine/src/simple_engine.rs`

- [ ] **Step 1: Add MS level filtering in run_search()**

In `simple_engine.rs`, after reading spectra (Step 3 comment), add filtering to only search MS2 spectra, and use `match_spectrum_all` when a spectrum has multiple precursors:

Replace the matching loop (Step 4 section, approximately lines 120-140) with:

```rust
    // Step 4: Match each spectrum against peptide candidates
    // Filter to MS2 only (MS1 survey scans have no precursors to match)
    let ms2_spectra: Vec<&Spectrum> = all_spectra
        .iter()
        .filter(|s| s.ms_level == MsLevel::MS2)
        .collect();
    let total_spectra = ms2_spectra.len();
    let mut psms: Vec<Psm> = Vec::new();

    for (i, spectrum) in ms2_spectra.iter().enumerate() {
        if i % 50 == 0 || i + 1 == total_spectra {
            let pct = 0.15 + 0.75 * (i as f64 / total_spectra.max(1) as f64);
            report(
                &format!("Matching spectra ({}/{})", i + 1, total_spectra),
                pct,
            );
        }

        if spectrum.precursors.len() > 1 {
            // DIA mode: multiple precursors, collect all matches
            let matches = match_spectrum_all(
                spectrum,
                &all_peptides,
                &params.precursor_tolerance,
                &params.fragment_tolerance,
                &params.fixed_modifications,
            );
            for m in &matches {
                psms.push(build_psm(spectrum, m, &params.fixed_modifications));
            }
        } else {
            // DDA mode: single precursor, use original function
            if let Some(m) = match_spectrum(
                spectrum,
                &all_peptides,
                &params.precursor_tolerance,
                &params.fragment_tolerance,
                &params.fixed_modifications,
            ) {
                psms.push(build_psm(spectrum, &m, &params.fixed_modifications));
            }
        }
    }
```

Add the import at the top of the file:

```rust
use crate::matching::match_spectrum_all;
```

Also update the summary to use `ms2_spectra.len()` instead of `all_spectra.len()`:

```rust
    let summary = build_summary(&psms, ms2_spectra.len() as u64, duration);
```

- [ ] **Step 2: Run tests**

```bash
cd /home/verden/pfind/2026-spring/code/proteinCopilot
cargo test --quiet
cargo clippy -- -D warnings
```

Expected: all existing tests still pass (DDA behavior unchanged — spectra with 0-1 precursors use original path).

- [ ] **Step 3: Commit**

```bash
git add crates/search-engine/src/simple_engine.rs
git commit -m "feat(search-engine): MS level filtering + DIA multi-precursor search

- Skip MS1 spectra in search loop (only match MS2)
- Use match_spectrum_all() for spectra with multiple precursors
- Retain single-precursor path for DDA backward compatibility

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 4: `extract_dia_precursors` MCP tool

**Files:**
- Modify: `crates/mcp-server/Cargo.toml`
- Modify: `crates/mcp-server/src/tools.rs`

- [ ] **Step 1: Add dia-extraction dependency**

In `crates/mcp-server/Cargo.toml`, add under `[dependencies]`:

```toml
protein-copilot-dia-extraction = { workspace = true }
```

- [ ] **Step 2: Add input/output structs**

In `crates/mcp-server/src/tools.rs`, add the input struct near other input structs:

```rust
#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ExtractDiaPrecursorsInput {
    /// Path to the spectrum file (.mzML)
    file_path: String,
    /// Output mode: "multi" (multiple precursors per spectrum) or "pseudo" (one precursor per spectrum). Default: "pseudo"
    #[serde(default = "default_output_mode")]
    output_mode: String,
    /// Minimum charge state to consider (default: 2)
    min_charge: Option<i32>,
    /// Maximum charge state to consider (default: 5)
    max_charge: Option<i32>,
    /// Override acquisition mode detection: "DDA" or "DIA". If not set, auto-detects.
    acquisition_mode: Option<String>,
}

fn default_output_mode() -> String {
    "pseudo".to_string()
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct DiaExtractionOutput {
    detected_mode: String,
    ms1_count: u32,
    ms2_count: u32,
    total_precursors_extracted: u32,
    avg_precursors_per_ms2: f64,
    charge_distribution: std::collections::HashMap<i32, u32>,
    output_spectra_count: u32,
    run_id: String,
    message: String,
}
```

- [ ] **Step 3: Add a DIA spectra cache to ProteinCopilotServer**

Add a field to `ProteinCopilotServer` to cache extracted spectra so `run_search` can use them:

```rust
pub struct ProteinCopilotServer {
    tool_router: ToolRouter<Self>,
    registry: protein_copilot_search_engine::EngineRegistry,
    run_cache: RunCache,
    /// Cache of DIA-extracted spectra, keyed by run_id from extract_dia_precursors.
    dia_cache: Arc<Mutex<HashMap<Uuid, Vec<Spectrum>>>>,
}
```

Update `new()`:

```rust
Self {
    tool_router: Self::tool_router(),
    registry,
    run_cache: Arc::new(Mutex::new(OrderedRunCache::new())),
    dia_cache: Arc::new(Mutex::new(HashMap::new())),
}
```

Add imports:

```rust
use protein_copilot_dia_extraction::{
    extract_dia_precursors, DiaExtractionConfig, IsotopePatternExtractor,
};
use protein_copilot_core::spectrum::AcquisitionMode;
```

Note: `IsotopePatternExtractor` needs to be re-exported from `dia-extraction/src/lib.rs`. Add `pub use isotope::IsotopePatternExtractor;` to `crates/dia-extraction/src/lib.rs`.

- [ ] **Step 4: Implement the MCP tool**

Add the tool method to the `impl ProteinCopilotServer` block:

```rust
#[rmcp::tool(
    name = "extract_dia_precursors",
    description = "Extract candidate precursor ions from DIA mass spectrometry data. Reads mzML file, detects DIA mode from isolation window widths, extracts precursor candidates from MS1 isotope patterns, and caches enhanced spectra for use with run_search. Returns extraction statistics."
)]
fn extract_dia_precursors(
    &self,
    #[tool(params)] input: ExtractDiaPrecursorsInput,
) -> Result<CallToolResult, ErrorData> {
    let path = Path::new(&input.file_path);
    let info = protein_copilot_spectrum_io::detect_format(path)
        .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
    let reader = protein_copilot_spectrum_io::create_reader(&info);
    let spectra = reader
        .read_all(path)
        .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

    // Configure extractor
    let mut extractor = IsotopePatternExtractor::default();
    if let Some(min_c) = input.min_charge {
        extractor.min_charge = min_c;
    }
    if let Some(max_c) = input.max_charge {
        extractor.max_charge = max_c;
    }

    // Configure extraction
    let acq_mode = input.acquisition_mode.as_deref().and_then(|m| match m {
        "DDA" | "dda" => Some(AcquisitionMode::DDA),
        "DIA" | "dia" => Some(AcquisitionMode::DIA),
        _ => None,
    });
    let config = DiaExtractionConfig {
        acquisition_mode: acq_mode,
        ..DiaExtractionConfig::default()
    };

    let result = extract_dia_precursors(&spectra, &extractor, &config)
        .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, &e.to_string()))?;

    let output_spectra = if input.output_mode == "multi" {
        result.enhanced_spectra.clone()
    } else {
        result.expand_to_pseudo_spectra()
    };

    let output_count = output_spectra.len() as u32;
    let run_id = Uuid::new_v4();

    // Cache for run_search
    if let Ok(mut cache) = self.dia_cache.lock() {
        cache.insert(run_id, output_spectra);
    }

    let output = DiaExtractionOutput {
        detected_mode: format!("{}", result.detected_mode),
        ms1_count: result.stats.ms1_count,
        ms2_count: result.stats.ms2_count,
        total_precursors_extracted: result.stats.total_precursors_extracted,
        avg_precursors_per_ms2: result.stats.avg_precursors_per_ms2,
        charge_distribution: result.stats.charge_distribution,
        output_spectra_count: output_count,
        run_id: run_id.to_string(),
        message: format!(
            "DIA extraction complete. {} precursors extracted from {} MS2 spectra. \
             Use run_id '{}' with run_search to search these spectra.",
            result.stats.total_precursors_extracted,
            result.stats.ms2_count,
            run_id
        ),
    };

    Ok(CallToolResult::success(vec![Content::json(output).map_err(
        |e| mcp_err(ErrorCode::INTERNAL_ERROR, &e.to_string()),
    )?]))
}
```

- [ ] **Step 5: Build and verify**

```bash
cd /home/verden/pfind/2026-spring/code/proteinCopilot
cargo build
cargo test --quiet
cargo clippy -- -D warnings
```

Expected: builds and all tests pass.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "feat(mcp-server): add extract_dia_precursors MCP tool

- Reads mzML, auto-detects DIA, extracts precursors from MS1 isotope patterns
- Supports multi-precursor and pseudo-spectra output modes
- Caches enhanced spectra for downstream run_search usage
- Configurable charge range and acquisition mode override

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 5: Update agent/prompt documentation

**Files:**
- Modify: `.github/agents/proteomics-search-assistant.agent.md`
- Modify: `docs/architecture.md`

- [ ] **Step 1: Update agent definition**

In `.github/agents/proteomics-search-assistant.agent.md`, add `extract_dia_precursors` to the tool list and add DIA workflow guidance:

Add to the tools section:
```markdown
- **extract_dia_precursors**: Extract candidate precursor ions from DIA data. 
  Reads mzML, detects DIA mode, extracts precursors from MS1 isotope patterns.
  Use before run_search for DIA data. Returns a run_id for the extracted spectra.
```

Add DIA workflow guidance:
```markdown
### DIA Data Workflow
1. Use `read_spectra` to check if data is DIA (wide isolation windows)
2. Call `extract_dia_precursors` to extract candidate precursors from MS1
3. Use the returned run_id with `run_search` to search the extracted spectra
```

- [ ] **Step 2: Update architecture.md**

Add a section about the dia-extraction crate in the architecture documentation. Include the data flow diagram from the design spec.

- [ ] **Step 3: Commit**

```bash
git add -A
git commit -m "docs: update agent and architecture for DIA extraction

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```
