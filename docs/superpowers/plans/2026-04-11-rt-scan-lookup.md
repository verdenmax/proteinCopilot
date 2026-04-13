# RT-Based Scan Lookup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Allow `annotate_spectrum` and `extract_xic` to accept `retention_time_min` instead of `scan_number`, auto-matching the best MS2 scan via RT + precursor m/z.

**Architecture:** Extract a public `find_scan_by_rt()` function from existing `scan_matcher.rs` private helpers. Add optional `retention_time_min` field to both MCP tool input structs. When `scan_number` is 0 and RT is provided, resolve the scan automatically before proceeding with existing logic.

**Tech Stack:** Rust, result-import crate (scan_matcher), mcp-server crate (tools.rs)

---

### Task 1: Extract public `find_scan_by_rt()` in scan_matcher.rs

**Files:**
- Modify: `crates/result-import/src/scan_matcher.rs`
- Modify: `crates/result-import/src/error.rs`

- [ ] **Step 1: Add `NoMatchingScan` error variant**

In `crates/result-import/src/error.rs`, add after the `MzmlNotFound` variant:

```rust
    #[error("no MS2 scan found near RT={rt_min:.3f} min (±{tolerance_min} min) with precursor_mz={precursor_mz:.4f}")]
    NoMatchingScan {
        rt_min: f64,
        tolerance_min: f64,
        precursor_mz: f64,
    },
```

- [ ] **Step 2: Make `Ms2Info` and `collect_ms2_info` public**

In `crates/result-import/src/scan_matcher.rs`:

Change `struct Ms2Info` from private to public:
```rust
/// MS2 spectrum info extracted from mzML for scan matching.
#[derive(Debug, Clone)]
pub struct Ms2Info {
    pub scan_number: u32,
    pub rt_min: f64,
    /// (target_mz, lower_offset, upper_offset)
    pub isolation_window: Option<(f64, f64, f64)>,
}
```

Change `fn collect_ms2_info` from private to public:
```rust
/// Collect (scan, rt, isolation_window) for all MS2 spectra from a reader.
pub fn collect_ms2_info(
    reader: &dyn SpectrumReader,
    path: &Path,
) -> Result<Vec<Ms2Info>, ResultImportError> {
```

Change `fn find_best_match` from private to public:
```rust
/// Find the best matching MS2 for a given RT and precursor m/z.
///
/// Returns `(scan_number, rt_delta_min)` or `None` if no match found.
pub fn find_best_match(
    sorted_ms2: &[Ms2Info],
    psm_rt_min: f64,
    psm_precursor_mz: f64,
    rt_tolerance_min: f64,
) -> Option<(u32, f64)> {
```

- [ ] **Step 3: Add `find_scan_by_rt()` public function**

Add this function at the end of the public API section (after `match_scans`), before the private helpers:

```rust
/// Find the best MS2 scan matching a given RT and precursor m/z.
///
/// Single-lookup convenience wrapper around `collect_ms2_info` + `find_best_match`.
/// Used by `annotate_spectrum` and `extract_xic` when the user provides RT instead
/// of scan_number.
///
/// # Arguments
/// * `file` — path to the mzML file
/// * `rt_min` — target retention time in minutes
/// * `precursor_mz` — precursor m/z to match against isolation windows
/// * `rt_tolerance_min` — RT tolerance in minutes (default: 0.5)
/// * `reader` — spectrum reader for the file
pub fn find_scan_by_rt(
    file: &Path,
    rt_min: f64,
    precursor_mz: f64,
    rt_tolerance_min: f64,
    reader: &dyn SpectrumReader,
) -> Result<u32, ResultImportError> {
    let mut ms2_infos = collect_ms2_info(reader, file)?;
    ms2_infos.sort_by(|a, b| {
        a.rt_min
            .partial_cmp(&b.rt_min)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    find_best_match(&ms2_infos, rt_min, precursor_mz, rt_tolerance_min)
        .map(|(scan, _delta)| scan)
        .ok_or(ResultImportError::NoMatchingScan {
            rt_min,
            tolerance_min: rt_tolerance_min,
            precursor_mz,
        })
}
```

- [ ] **Step 4: Write tests for `find_scan_by_rt`**

Add at the bottom of `crates/result-import/src/scan_matcher.rs` (in a `#[cfg(test)] mod tests` block, or add to existing tests):

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn find_best_match_selects_closest_rt_within_window() {
        let ms2s = vec![
            Ms2Info { scan_number: 10, rt_min: 1.0, isolation_window: Some((500.0, 1.0, 1.0)) },
            Ms2Info { scan_number: 20, rt_min: 2.0, isolation_window: Some((500.0, 1.0, 1.0)) },
            Ms2Info { scan_number: 30, rt_min: 3.0, isolation_window: Some((500.0, 1.0, 1.0)) },
        ];
        // Closest to RT=2.1 with mz=500.0 → scan 20
        let result = find_best_match(&ms2s, 2.1, 500.0, 0.5);
        assert_eq!(result, Some((20, -0.09999999999999998))); // delta = 2.0 - 2.1

        // mz outside isolation window → no match
        let result = find_best_match(&ms2s, 2.1, 600.0, 0.5);
        assert_eq!(result, None);
    }

    #[test]
    fn find_best_match_no_candidates_in_tolerance() {
        let ms2s = vec![
            Ms2Info { scan_number: 10, rt_min: 1.0, isolation_window: Some((500.0, 1.0, 1.0)) },
        ];
        // RT=5.0 is far from RT=1.0 with tolerance=0.5
        let result = find_best_match(&ms2s, 5.0, 500.0, 0.5);
        assert_eq!(result, None);
    }

    #[test]
    fn find_best_match_no_isolation_window_accepts_by_rt() {
        let ms2s = vec![
            Ms2Info { scan_number: 10, rt_min: 2.0, isolation_window: None },
        ];
        // No isolation window → accept based on RT only (DDA fallback)
        let result = find_best_match(&ms2s, 2.1, 999.0, 0.5);
        assert_eq!(result, Some((10, -0.09999999999999998)));
    }

    #[test]
    fn find_best_match_empty_list() {
        let result = find_best_match(&[], 2.0, 500.0, 0.5);
        assert_eq!(result, None);
    }
}
```

- [ ] **Step 5: Run tests**

Run: `cargo test -p protein-copilot-result-import --quiet`
Expected: All pass (existing + 4 new)

- [ ] **Step 6: Commit**

```bash
git add crates/result-import/src/scan_matcher.rs crates/result-import/src/error.rs
git commit -m "feat: extract public find_scan_by_rt() from scan_matcher

Public API for single RT→scan lookup, reusing existing binary search +
isolation window matching. Used by annotate_spectrum and extract_xic
when user provides RT instead of scan_number.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: Add `retention_time_min` to `annotate_spectrum`

**Files:**
- Modify: `crates/mcp-server/src/tools.rs` (AnnotateSpectrumInput struct + handler)

- [ ] **Step 1: Add field to AnnotateSpectrumInput**

In `crates/mcp-server/src/tools.rs`, in the `AnnotateSpectrumInput` struct, after the `charge` field (~line 248), add:

```rust
    /// Retention time in minutes — alternative to scan_number for auto scan lookup.
    /// When scan_number is 0 and this is provided, the tool finds the closest MS2
    /// scan matching this RT and precursor_mz in the mzML file.
    #[serde(default)]
    retention_time_min: Option<f64>,
```

- [ ] **Step 2: Add RT-based scan resolution in the handler**

In the `annotate_spectrum` handler function, find the line that reads the spectrum (after Mode 1/Mode 2 resolution):
```rust
        let reader = self.get_or_create_reader(&spectrum_file)?;
        let spectrum = reader
            .read_spectrum(&spectrum_file, input.scan_number)
```

Replace the scan_number resolution block. Insert **before** `let reader = ...`:

```rust
        // Resolve scan_number: if 0 and retention_time_min provided, auto-match
        let resolved_scan = if input.scan_number == 0 {
            if let Some(rt) = input.retention_time_min {
                let reader = self.get_or_create_reader(&spectrum_file)?;
                // Need precursor_mz for isolation window matching
                let precursor_mz_for_lookup = protein_copilot_search_engine::matching::calculate_precursor_mz(
                    &peptide_seq, charge, &modifications,
                );
                protein_copilot_result_import::scan_matcher::find_scan_by_rt(
                    &spectrum_file,
                    rt,
                    precursor_mz_for_lookup,
                    0.5, // default RT tolerance
                    &*reader,
                )
                .map_err(|e| mcp_err(ErrorCode::INVALID_PARAMS, e))?
            } else {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "scan_number is 0: provide either a valid scan_number or retention_time_min for auto lookup",
                ));
            }
        } else {
            input.scan_number
        };
```

Then update the subsequent code to use `resolved_scan` instead of `input.scan_number`:
- `reader.read_spectrum(&spectrum_file, resolved_scan)` 
- `output_path` default: `format!("output/annotation_scan{}.html", resolved_scan)`
- The `AnnotateResult` at the end: `scan_number: annotation.scan_number` (already uses annotation which has correct scan)

- [ ] **Step 3: Check `calculate_precursor_mz` exists**

Verify `protein_copilot_search_engine::matching::calculate_precursor_mz` exists. If not, compute it inline:

```rust
// Inline calculation if no helper exists:
let base_mass: f64 = peptide_seq.chars().map(|c| {
    protein_copilot_search_engine::matching::amino_acid_mass(c)
}).sum::<f64>() + 18.010565; // water
let mod_mass: f64 = modifications.iter().map(|m| m.mass_delta).sum();
let precursor_mz_for_lookup = (base_mass + mod_mass + charge as f64 * 1.007276) / charge as f64;
```

Alternatively, just use the `precursor_mz` if available from Mode 1 (run_id PSM already has it). For Mode 2, we can compute from sequence. The actual value only needs to be approximate — it's checked against a wide isolation window (typically 2-25 Da).

- [ ] **Step 4: Run build**

Run: `cargo build --workspace`
Expected: Compiles successfully

- [ ] **Step 5: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "feat: annotate_spectrum supports retention_time_min auto lookup

When scan_number=0 and retention_time_min is provided, automatically
finds the best MS2 scan via RT + precursor_mz isolation window matching.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: Add `retention_time_min` to `extract_xic`

**Files:**
- Modify: `crates/mcp-server/src/tools.rs` (ExtractXicInput struct + handler)

- [ ] **Step 1: Add field to ExtractXicInput**

In the `ExtractXicInput` struct (~line 336), after `scan_number`, add:

```rust
    /// Retention time in minutes — alternative to scan_number.
    /// When scan_number is 0 and this is provided, auto-finds closest MS2 scan.
    #[serde(default)]
    #[schemars(description = "Retention time in minutes. When scan_number is 0, auto-finds the closest MS2 scan matching this RT and precursor_mz.")]
    retention_time_min: Option<f64>,
```

- [ ] **Step 2: Add RT-based scan resolution in extract_xic handler**

Find the `extract_xic` handler function. Locate where it uses `input.scan_number`. Add the same resolution pattern as Task 2:

```rust
        // Resolve scan_number: if 0 and retention_time_min provided, auto-match
        let resolved_scan = if input.scan_number == 0 {
            if let Some(rt) = input.retention_time_min {
                let precursor_mz = input.precursor_mz.ok_or_else(|| {
                    mcp_err(ErrorCode::INVALID_PARAMS, 
                        "precursor_mz is required when using retention_time_min for scan lookup")
                })?;
                let reader = self.get_or_create_reader(&file_path)?;
                protein_copilot_result_import::scan_matcher::find_scan_by_rt(
                    &file_path,
                    rt,
                    precursor_mz,
                    0.5,
                    &*reader,
                )
                .map_err(|e| mcp_err(ErrorCode::INVALID_PARAMS, e))?
            } else {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    "scan_number is 0: provide either a valid scan_number or retention_time_min for auto lookup",
                ));
            }
        } else {
            input.scan_number
        };
```

Then replace all `input.scan_number` with `resolved_scan` in the handler body.

- [ ] **Step 3: Run build + tests**

Run: `cargo build --workspace && cargo test --workspace --quiet`
Expected: All 509+ tests pass

- [ ] **Step 4: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "feat: extract_xic supports retention_time_min auto lookup

When scan_number=0 and retention_time_min is provided, auto-finds
the closest MS2 scan via RT + precursor_mz matching.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 4: Update MCP tool descriptions + final verification

**Files:**
- Modify: `crates/mcp-server/src/tools.rs` (tool description strings)

- [ ] **Step 1: Update annotate_spectrum description**

Change the tool description (~line 1431) to mention the new capability:

```rust
        description = "Annotate a single spectrum with peptide fragment ion matching. Generates an interactive HTML file showing matched b/y ions. Two modes: (1) provide run_id + scan_number to annotate an existing PSM, or (2) provide file_path + scan_number + peptide_sequence + charge for manual annotation. In mode 2, you can set scan_number=0 and provide retention_time_min to auto-find the matching scan."
```

- [ ] **Step 2: Update extract_xic scan_number description**

Update the `scan_number` field's `schemars(description)` in `ExtractXicInput`:

```rust
    #[schemars(description = "Scan number (1-based) to center the XIC around. Set to 0 with retention_time_min to auto-find scan.")]
    scan_number: u32,
```

And in `AnnotateSpectrumInput`, update the scan_number doc comment:

```rust
    /// Scan number (1-based) to annotate. Set to 0 with retention_time_min for auto lookup.
    scan_number: u32,
```

- [ ] **Step 3: Run full test suite**

Run: `cargo test --workspace`
Expected: All tests pass

- [ ] **Step 4: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "docs: update tool descriptions for RT-based scan lookup

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```
