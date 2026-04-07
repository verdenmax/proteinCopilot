# DIA Precursor Tools & Param-Recommend DIA Detection

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add per-spectrum precursor extraction MCP tool so users can inspect DIA isotope patterns for individual spectra, and enhance param-recommend to detect DIA data and guide users appropriately.

**Architecture:** Two independent features. Feature A adds a new library function in dia-extraction + new MCP tool + prompt skill. Feature B adds DIA detection heuristics to param-recommend by extending SpectrumSummary with isolation window statistics. Both features reuse existing infrastructure (IsotopePatternExtractor, correlation module).

**Tech Stack:** Rust, protein-copilot workspace, MCP (rmcp), schemars/serde

---

### Task 1: Add `median_isolation_window_da` to SpectrumSummary

**Context:** The `SpectrumSummary` struct lacks isolation window information. Adding `median_isolation_window_da: Option<f64>` enables both param-recommend DIA detection and gives users a quick signal of DDA vs DIA when reading spectra. The field is `Option` because MGF files never have isolation windows.

**Files:**
- Modify: `crates/core/src/spectrum.rs` (SpectrumSummary struct)
- Modify: `crates/spectrum-io/src/util.rs` (SummaryAccumulator)
- Test: existing tests in `crates/spectrum-io/src/mzml.rs` and `crates/spectrum-io/src/mgf.rs`

- [ ] **Step 1: Add field to SpectrumSummary**

In `crates/core/src/spectrum.rs`, add after `median_peaks_per_spectrum`:

```rust
    /// Median number of peaks per spectrum.
    pub median_peaks_per_spectrum: u32,
    /// Median isolation window width in Da (`None` if no isolation windows found).
    /// Useful for DIA detection: DDA windows are typically < 3 Da, DIA > 5 Da.
    #[serde(default)]
    pub median_isolation_window_da: Option<f64>,
}
```

- [ ] **Step 2: Update SummaryAccumulator to collect isolation window widths**

In `crates/spectrum-io/src/util.rs`, add `isolation_widths: Vec<f64>` to `SummaryAccumulator`:

```rust
pub(crate) struct SummaryAccumulator {
    total: u64,
    ms1_count: u64,
    ms2_count: u64,
    mz_min: f64,
    mz_max: f64,
    rt_min: f64,
    rt_max: f64,
    charge_dist: HashMap<i32, u64>,
    peak_counts: Vec<u32>,
    isolation_widths: Vec<f64>,
}
```

Initialize in `new()`:

```rust
isolation_widths: Vec::new(),
```

In `observe()`, after the precursor charge loop (after line 125), add:

```rust
        for p in &s.precursors {
            if let Some(iw) = &p.isolation_window {
                let width = iw.lower_offset + iw.upper_offset;
                if width.is_finite() && width > 0.0 {
                    self.isolation_widths.push(width);
                }
            }
        }
```

In `into_summary()`, compute median isolation window width before building the struct (after computing `median_peaks`):

```rust
        self.isolation_widths.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median_iw = if self.isolation_widths.is_empty() {
            None
        } else {
            let mid = self.isolation_widths.len() / 2;
            Some(self.isolation_widths[mid])
        };
```

And add the field to the SpectrumSummary construction:

```rust
            median_peaks_per_spectrum: median_peaks,
            median_isolation_window_da: median_iw,
```

- [ ] **Step 3: Fix all SpectrumSummary constructors in tests**

Search for `SpectrumSummary {` and add `median_isolation_window_da: None,` to each test constructor across:
- `crates/core/src/spectrum.rs` (test functions)
- `crates/param-recommend/src/rules.rs` (test functions)
- Any other files

- [ ] **Step 4: Run full workspace tests**

Run: `cargo test --workspace`
Expected: ALL PASS

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(core): add median_isolation_window_da to SpectrumSummary

Extends SpectrumSummary with median isolation window width in Da,
computed from precursor isolation windows during summary accumulation.
Returns None for MGF files (no isolation windows) or files without
precursor isolation data. Enables DIA detection in param-recommend.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: Add DIA detection to param-recommend

**Context:** The param-recommend engine never detects DIA or sets `acquisition_mode`. With `median_isolation_window_da` now available in SpectrumSummary, we can detect DIA (window > 5 Da) and:
1. Set `acquisition_mode: Some(AcquisitionMode::DIA)` in recommended params
2. Add a DIA-specific note to the explanation, telling the user to run `extract_dia_precursors` first

**Files:**
- Modify: `crates/param-recommend/src/rules.rs` (recommend function + new DIA detection helper)
- Test: `crates/param-recommend/src/rules.rs` (in `mod tests`)

- [ ] **Step 1: Write failing tests**

Add tests in `crates/param-recommend/src/rules.rs` mod tests:

```rust
#[test]
fn recommend_detects_dia_from_wide_isolation_window() {
    let summary = make_summary_with_iw(Some(25.0)); // 25 Da = DIA
    let result = recommend(&summary, None).unwrap();
    assert_eq!(
        result.decision.acquisition_mode,
        Some(AcquisitionMode::DIA)
    );
    assert!(result.explanation.contains("DIA"));
}

#[test]
fn recommend_detects_dda_from_narrow_isolation_window() {
    let summary = make_summary_with_iw(Some(2.0)); // 2 Da = DDA
    let result = recommend(&summary, None).unwrap();
    assert_eq!(result.decision.acquisition_mode, None);
}

#[test]
fn recommend_no_isolation_window_stays_none() {
    let summary = make_summary_with_iw(None); // MGF, no IW
    let result = recommend(&summary, None).unwrap();
    assert_eq!(result.decision.acquisition_mode, None);
}
```

Also add a helper `make_summary_with_iw(iw: Option<f64>) -> SpectrumSummary` that creates a valid summary with the given `median_isolation_window_da`.

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p protein-copilot-param-recommend recommend_detects_dia -- --nocapture`
Expected: FAIL — `acquisition_mode` is always None

- [ ] **Step 3: Add DIA detection logic in `recommend()`**

In `crates/param-recommend/src/rules.rs`, add after Step 3 (enzyme hint override, ~line 56):

```rust
    // Step 3b: Detect DIA from isolation window width or user hints
    let is_dia = detect_dia(summary, hints);
    if is_dia {
        base.acquisition_mode = Some(AcquisitionMode::DIA);
    }
```

Add import at top of file:

```rust
use protein_copilot_core::spectrum::AcquisitionMode;
```

Add the detection function:

```rust
/// Detect DIA acquisition mode from summary data or user hints.
///
/// Detection heuristics:
/// 1. User hint `experiment_type` containing "DIA" → definitive
/// 2. `median_isolation_window_da` > 5.0 Da → strong DIA signal
fn detect_dia(summary: &SpectrumSummary, hints: Option<&UserHints>) -> bool {
    // Check user hint first
    if let Some(exp_type) = hints.and_then(|h| h.experiment_type.as_deref()) {
        if exp_type.to_lowercase().contains("dia") {
            return true;
        }
    }

    // Auto-detect from isolation window width
    if let Some(median_iw) = summary.median_isolation_window_da {
        return median_iw > 5.0;
    }

    false
}
```

- [ ] **Step 4: Append DIA note to explanation**

In `build_explanation()`, or after it returns (in `recommend()`, after line 59), if `is_dia`, append:

```rust
    let dia_note = if is_dia {
        "\n\n🔬 DIA data detected (wide isolation windows). \
         Recommended workflow: call `extract_dia_precursors` first to \
         extract candidate precursors from MS1 isotope patterns, then \
         pass the resulting `dia_run_id` to `run_search`."
    } else {
        ""
    };
```

Append `dia_note` to the `final_explanation`:

```rust
    let final_explanation = if warnings.is_empty() {
        format!("{explanation}{dia_note}")
    } else {
        format!("{explanation}\n\n⚠ Warnings:\n{}{dia_note}", warnings.join("\n"))
    };
```

- [ ] **Step 5: Run all param-recommend tests**

Run: `cargo test -p protein-copilot-param-recommend -- --nocapture`
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add crates/param-recommend/src/rules.rs
git commit -m "feat(param-recommend): detect DIA from isolation window width

Adds DIA detection heuristic: median isolation window > 5 Da triggers
acquisition_mode = DIA in recommended params. Also detects DIA from
user hint experiment_type containing 'DIA'. Appends workflow guidance
to explanation when DIA is detected.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: Add `extract_single_spectrum_precursors()` library function

**Context:** The existing `extract_dia_precursors()` processes an entire file. We need a focused function that extracts precursor candidates from a single MS2 spectrum by correlating it to its MS1 and running isotope pattern detection. This function lives in the dia-extraction library crate.

**Files:**
- Modify: `crates/dia-extraction/src/lib.rs` (add new public function)
- Test: `crates/dia-extraction/src/lib.rs` (in `mod tests`)

- [ ] **Step 1: Define the return type**

Add to `crates/dia-extraction/src/config.rs` (or `lib.rs`):

```rust
/// Result of extracting precursors from a single MS2 spectrum.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SingleSpectrumExtractionResult {
    /// The MS2 scan number that was analyzed.
    pub ms2_scan: u32,
    /// The MS1 scan used for isotope pattern extraction.
    pub ms1_scan_used: u32,
    /// How the MS1 was selected.
    pub correlation_method: String,
    /// Isolation window of the MS2 spectrum (if available).
    pub isolation_window: Option<IsolationWindow>,
    /// Extracted precursor candidates from isotope pattern analysis.
    pub precursors: Vec<PrecursorInfo>,
}
```

Add necessary imports (serde, schemars, IsolationWindow, PrecursorInfo).

- [ ] **Step 2: Write failing test**

In `crates/dia-extraction/src/lib.rs` mod tests:

```rust
#[test]
fn test_extract_single_spectrum_precursors() {
    let ms1 = make_ms1(
        1,
        10.0,
        vec![500.0, 500.502, 501.003, 600.0, 600.335, 600.669],
        vec![10000.0, 5000.0, 2500.0, 8000.0, 4000.0, 2000.0],
    );
    let ms2 = make_ms2_dia(2, 10.1, 550.0, 50.0, 50.0, Some(1));

    let all_spectra = vec![ms1, ms2];
    let extractor = IsotopePatternExtractor::default();

    let result = extract_single_spectrum_precursors(
        &all_spectra,
        2, // ms2 scan_number
        &extractor,
    )
    .unwrap();

    assert_eq!(result.ms2_scan, 2);
    assert_eq!(result.ms1_scan_used, 1);
    assert_eq!(result.correlation_method, "source_scan");
    assert!(!result.precursors.is_empty(), "should find isotope clusters");
}

#[test]
fn test_extract_single_spectrum_scan_not_found() {
    let ms1 = make_ms1(1, 10.0, vec![500.0], vec![1000.0]);
    let all_spectra = vec![ms1];
    let extractor = IsotopePatternExtractor::default();

    let result = extract_single_spectrum_precursors(&all_spectra, 99, &extractor);
    assert!(result.is_err());
}
```

- [ ] **Step 3: Implement `extract_single_spectrum_precursors()`**

In `crates/dia-extraction/src/lib.rs`:

```rust
/// Extract candidate precursors from a single MS2 spectrum.
///
/// Finds the specified MS2 spectrum, correlates it to the nearest MS1,
/// and runs isotope pattern detection within the MS2's isolation window.
pub fn extract_single_spectrum_precursors(
    spectra: &[Spectrum],
    ms2_scan_number: u32,
    extractor: &dyn PrecursorExtractor,
) -> Result<SingleSpectrumExtractionResult, DiaExtractionError> {
    // Find the target MS2 spectrum
    let ms2 = spectra
        .iter()
        .find(|s| s.scan_number == ms2_scan_number && s.ms_level == MsLevel::MS2)
        .ok_or_else(|| DiaExtractionError::ScanNotFound {
            scan: ms2_scan_number,
        })?;

    // Get isolation window from the MS2 spectrum
    let isolation_window = ms2
        .precursors
        .first()
        .and_then(|p| p.isolation_window.clone());

    let iw = isolation_window.clone().ok_or_else(|| {
        DiaExtractionError::NoIsolationWindow {
            scan: ms2_scan_number,
        }
    })?;

    // Collect MS1 spectra references
    let ms1_refs: Vec<&Spectrum> = spectra
        .iter()
        .filter(|s| s.ms_level == MsLevel::MS1)
        .collect();

    if ms1_refs.is_empty() {
        return Err(DiaExtractionError::NoMs1Spectra);
    }

    // Correlate this MS2 to an MS1
    let ms2_refs: Vec<&Spectrum> = vec![ms2];
    let correlations = correlation::correlate_ms1_ms2(&ms1_refs, &ms2_refs);

    let ms1_idx = correlations[0].ok_or(DiaExtractionError::NoMs1Spectra)?;
    let ms1 = ms1_refs[ms1_idx];

    // Determine correlation method
    let method = determine_correlation_method(ms2, ms1);

    // Extract precursors from the correlated MS1
    let precursors = extractor.extract(ms1, &iw);

    Ok(SingleSpectrumExtractionResult {
        ms2_scan: ms2_scan_number,
        ms1_scan_used: ms1.scan_number,
        correlation_method: method,
        isolation_window,
        precursors,
    })
}

/// Determine which correlation method was used (for reporting).
fn determine_correlation_method(ms2: &Spectrum, ms1: &Spectrum) -> String {
    // Check if source_scan matches
    if let Some(source_scan) = ms2.precursors.first().and_then(|p| p.source_scan) {
        if ms1.scan_number == source_scan {
            return "source_scan".to_string();
        }
    }
    // Check if scan order matches (ms1.scan < ms2.scan)
    if ms1.scan_number < ms2.scan_number {
        return "scan_order".to_string();
    }
    "rt_nearest".to_string()
}
```

- [ ] **Step 4: Add new error variants**

In `crates/dia-extraction/src/error.rs`:

```rust
    /// Requested scan number not found as MS2.
    #[error("MS2 scan {scan} not found")]
    ScanNotFound { scan: u32 },

    /// MS2 spectrum has no isolation window (needed for precursor extraction).
    #[error("MS2 scan {scan} has no isolation window")]
    NoIsolationWindow { scan: u32 },
```

- [ ] **Step 5: Re-export the new type from `lib.rs`**

Add to `crates/dia-extraction/src/lib.rs` pub use block:

```rust
pub use config::{DiaExtractionConfig, DiaExtractionResult, ExtractionStats, SingleSpectrumExtractionResult};
```

- [ ] **Step 6: Run dia-extraction tests**

Run: `cargo test -p protein-copilot-dia-extraction -- --nocapture`
Expected: ALL PASS

- [ ] **Step 7: Commit**

```bash
git add crates/dia-extraction/
git commit -m "feat(dia-extraction): add single-spectrum precursor extraction

New public function extract_single_spectrum_precursors() extracts
precursor candidates from a single MS2 spectrum by correlating it
to the nearest MS1 and running isotope pattern detection within
the isolation window. Reports which correlation method was used.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 4: Add `extract_spectrum_precursors` MCP tool

**Context:** Wraps the library function from Task 3 as an MCP tool, allowing users to extract precursors from a single spectrum interactively.

**Files:**
- Modify: `crates/mcp-server/src/tools.rs` (add input struct + tool method)
- Modify: `crates/mcp-server/Cargo.toml` (if new deps needed — unlikely, dia-extraction already a dep)

- [ ] **Step 1: Add input/output structs**

In `crates/mcp-server/src/tools.rs`, after `ExtractDiaPrecursorsInput`:

```rust
#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ExtractSpectrumPrecursorsInput {
    /// Path to the spectrum file (.mzML)
    file_path: String,
    /// Scan number of the MS2 spectrum to analyze (1-based)
    scan_number: u32,
    /// Minimum charge state to consider (default: 2)
    min_charge: Option<i32>,
    /// Maximum charge state to consider (default: 5)
    max_charge: Option<i32>,
}
```

The output will be the `SingleSpectrumExtractionResult` from dia-extraction (already has Serialize + JsonSchema).

- [ ] **Step 2: Add the MCP tool method**

In the `impl ProteinCopilotServer` block, add:

```rust
    /// Extract precursor candidates from a single MS2 spectrum using MS1 isotope analysis.
    #[rmcp::tool(
        name = "extract_spectrum_precursors",
        description = "Extract precursor candidates from a single MS2 spectrum by analyzing its correlated MS1 isotope patterns. Input: file path + MS2 scan number. Output: list of precursor candidates with m/z, charge, intensity, plus which MS1 was used and how it was correlated. Useful for inspecting DIA extraction quality or understanding precursor detection for a specific spectrum."
    )]
    fn extract_spectrum_precursors(
        &self,
        Parameters(input): Parameters<ExtractSpectrumPrecursorsInput>,
    ) -> Result<Json<SingleSpectrumExtractionResult>, ErrorData> {
        use protein_copilot_dia_extraction::{
            extract_single_spectrum_precursors, IsotopePatternExtractor,
        };

        let path = Path::new(&input.file_path);
        let info = protein_copilot_spectrum_io::detect_format(path)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        let reader = protein_copilot_spectrum_io::create_reader(&info);
        let spectra = reader
            .read_all(path)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

        let mut extractor = IsotopePatternExtractor::default();
        if let Some(min_c) = input.min_charge {
            extractor.min_charge = min_c;
        }
        if let Some(max_c) = input.max_charge {
            extractor.max_charge = max_c;
        }

        let result = extract_single_spectrum_precursors(&spectra, input.scan_number, &extractor)
            .map_err(|e| mcp_err(ErrorCode::INVALID_PARAMS, e.to_string()))?;

        Ok(Json(result))
    }
```

- [ ] **Step 3: Register tool in tool_router**

The `#[rmcp::tool]` macro should auto-register. Verify by running clippy.

- [ ] **Step 4: Build and run clippy**

Run: `cargo clippy --workspace --all-targets -- -D warnings`
Expected: PASS

- [ ] **Step 5: Run full test suite**

Run: `cargo test --workspace`
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "feat(mcp-server): add extract_spectrum_precursors MCP tool

New tool extracts precursor candidates from a single MS2 spectrum
by correlating to its MS1 and running isotope pattern detection.
Returns precursor list, MS1 scan used, correlation method, and
isolation window. Enables interactive DIA precursor inspection.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 5: Add DIA precursor extraction prompt skill

**Context:** Create a prompt skill that guides the LLM through DIA precursor inspection and extraction workflows. This prompt teaches the LLM when and how to use `extract_spectrum_precursors` vs `extract_dia_precursors`.

**Files:**
- Create: `.github/prompts/dia-precursor-extraction.prompt.md`
- Modify: `.github/agents/proteomics-search.agent.md` (add new tool reference)

- [ ] **Step 1: Create the prompt skill**

Create `.github/prompts/dia-precursor-extraction.prompt.md`:

```markdown
---
description: "DIA precursor extraction and inspection — guides users through extracting, analyzing, and validating precursor candidates from DIA mass spectrometry data"
---

# DIA Precursor Extraction

You are a DIA data analysis specialist. Help users extract and understand precursor candidates from DIA mass spectrometry data.

## Available Tools

| Tool | Purpose |
|------|---------|
| `read_spectra` | Read file summary (check ms1_count, median_isolation_window_da for DIA detection) |
| `get_spectrum` | View a single spectrum's raw data and precursors |
| `extract_spectrum_precursors` | Extract precursor candidates from a single MS2 scan using MS1 isotope patterns |
| `extract_dia_precursors` | Extract precursors from entire file (for search pipeline) |
| `run_search` | Run search with `dia_run_id` from full-file extraction |

## Workflow: Single Spectrum Inspection

When a user wants to inspect precursor extraction for a specific scan:

1. **Read the spectrum**: `get_spectrum(file_path, scan_number)` → check MS level, isolation window, existing precursors
2. **Extract precursors**: `extract_spectrum_precursors(file_path, scan_number)` → get isotope pattern candidates
3. **Interpret results**:
   - Report each precursor candidate: m/z, charge state, intensity
   - Explain the correlation method used (source_scan > scan_order > rt_nearest)
   - Note the isolation window range
   - If no precursors found, explain possible reasons (no isotope clusters, MS1 too sparse, wrong charge range)
4. **Suggest tuning**: If results are poor, suggest adjusting min_charge/max_charge

## Workflow: Full File DIA Search

When a user wants to search DIA data:

1. **Read file summary**: `read_spectra(file_path)` → check `median_isolation_window_da`
   - If > 5 Da → DIA data, proceed with extraction
   - If < 3 Da → DDA data, skip extraction and use normal search
   - If None → MGF format (no isolation windows), use normal search
2. **Extract all precursors**: `extract_dia_precursors(file_path)` → get `dia_run_id`
3. **Review extraction stats**: Check `avg_precursors_per_ms2`, `charge_distribution`
   - If avg < 1.0 → extraction may have issues, suggest inspecting individual spectra
4. **Search**: `run_search(dia_run_id=..., database_path=..., params=...)`

## Interpreting Isotope Patterns

When explaining results to users:

- **Charge state z**: Isotope peaks are spaced by 1.003/z Da
  - z=2: peaks every ~0.502 Da
  - z=3: peaks every ~0.334 Da
  - z=4: peaks every ~0.251 Da
- **Pattern quality**: More isotope peaks = higher confidence
  - 2 peaks = minimum (possible false positive)
  - 3+ peaks = good confidence
  - 5+ peaks = high confidence
- **Monoisotopic peak**: The first (lightest) peak in the cluster is the monoisotopic m/z
- **Intensity pattern**: Typically decreasing, but for larger peptides the 2nd peak may be more intense

## Decision Tree

```
User has DIA data?
├── Yes: Wants to inspect specific scan?
│   ├── Yes → extract_spectrum_precursors (single scan)
│   └── No: Wants full search?
│       └── Yes → extract_dia_precursors → run_search(dia_run_id)
└── Unsure → read_spectra → check median_isolation_window_da
```
```

- [ ] **Step 2: Update proteomics-search agent**

In `.github/agents/proteomics-search.agent.md`, add `extract_spectrum_precursors` to the tools list and update the DIA workflow section to mention single-spectrum inspection.

- [ ] **Step 3: Commit**

```bash
git add .github/prompts/dia-precursor-extraction.prompt.md .github/agents/proteomics-search.agent.md
git commit -m "docs: add DIA precursor extraction prompt skill

New prompt skill guides LLM through DIA precursor inspection workflows:
single-spectrum extraction, full-file extraction, result interpretation,
and isotope pattern analysis. Updates agent to reference new tool.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 6: Final validation + documentation sync

**Context:** Run the full workspace verification and update documentation.

**Files:**
- Modify: `docs/mcp-tools.md` (add new tool)
- Verify: all workspace tests + clippy

- [ ] **Step 1: Run full workspace validation**

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
```

Expected: All pass, 0 warnings

- [ ] **Step 2: Update docs/mcp-tools.md**

Add `extract_spectrum_precursors` tool documentation with inputs, outputs, and usage examples. Update tool count from 13 to 14.

- [ ] **Step 3: Commit**

```bash
git add -A
git commit -m "docs: update MCP tools documentation for 14 tools

Added extract_spectrum_precursors tool documentation. Updated tool
count from 13 to 14.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---
