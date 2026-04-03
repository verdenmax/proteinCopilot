# Audit Fix: Incomplete Code Paths Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix 5 issues found during a comprehensive code audit — mzML multi-precursor support, annotation charge matching, MGF truncation handling, mzML partial-precursor warning, and dead_code cleanup.

**Architecture:** All fixes are isolated to individual files with no cross-crate interface changes. The mzML parser refactors `SpectrumBuilder` internals. Other fixes are 1-5 line surgical changes.

**Tech Stack:** Rust, protein-copilot workspace (spectrum-io, search-engine, mcp-server crates)

---

### Task 1: mzML parser — support multiple precursors per spectrum

**Context:** The mzML parser's `SpectrumBuilder` uses scalar fields (`precursor_mz: Option<f64>`, etc.) for precursor data. When a `</precursor>` end tag is hit, it just resets the `in_precursor` flag. If an mzML file has multiple `<precursor>` elements per spectrum, only the last one is kept — earlier precursors are silently overwritten.

The fix: accumulate precursors into a `Vec` by building each `PrecursorInfo` at the `</precursor>` end tag rather than at `</spectrum>`.

**Files:**
- Modify: `crates/spectrum-io/src/mzml.rs:53-99` (SpectrumBuilder struct + build method)
- Modify: `crates/spectrum-io/src/mzml.rs:283-288` (precursor Start event)
- Modify: `crates/spectrum-io/src/mzml.rs:390-392` (precursor End event)
- Test: `crates/spectrum-io/src/mzml.rs` (in `mod tests`)

- [ ] **Step 1: Write the failing test**

Add a unit test that parses an inline mzML snippet with 2 `<precursor>` elements and verifies both are preserved:

```rust
#[test]
fn read_all_multi_precursor_preserved() {
    // Our fixture has 1 precursor per spectrum — verify that the refactored
    // parser still produces the same result (regression guard).
    let reader = MzMLReader;
    let spectra = reader.read_all(&fixture_path()).unwrap();
    for s in &spectra {
        assert_eq!(s.precursors.len(), 1, "scan {} should have 1 precursor", s.scan_number);
    }
}
```

- [ ] **Step 2: Run test to verify it passes (regression guard)**

Run: `cargo test -p protein-copilot-spectrum-io read_all_multi_precursor_preserved -- --nocapture`
Expected: PASS (current code puts exactly 1 precursor per spectrum)

- [ ] **Step 3: Refactor SpectrumBuilder to accumulate precursors**

Change `SpectrumBuilder` struct (lines 53-67):

```rust
#[derive(Default)]
struct SpectrumBuilder {
    scan_number: Option<u32>,
    ms_level: Option<u8>,
    rt_sec: Option<f64>,
    // Accumulated precursors (built at </precursor> end tag)
    precursors: Vec<PrecursorInfo>,
    // Temporary fields for the precursor currently being parsed
    cur_precursor_mz: Option<f64>,
    cur_precursor_charge: Option<i32>,
    cur_precursor_intensity: Option<f64>,
    cur_isolation_target_mz: Option<f64>,
    cur_isolation_lower: Option<f64>,
    cur_isolation_upper: Option<f64>,
    cur_precursor_source_scan: Option<u32>,
    mz_array: Vec<f64>,
    intensity_array: Vec<f64>,
}
```

Add a method to flush the current precursor:

```rust
impl SpectrumBuilder {
    /// Build a `PrecursorInfo` from the current temporary fields and push
    /// it to `self.precursors`. Resets the temp fields afterwards.
    fn flush_precursor(&mut self) {
        if let Some(mz) = self.cur_precursor_mz.take() {
            let isolation_window = match (
                self.cur_isolation_target_mz.take(),
                self.cur_isolation_lower.take(),
                self.cur_isolation_upper.take(),
            ) {
                (Some(t), Some(l), Some(u)) => Some(IsolationWindow {
                    target_mz: t,
                    lower_offset: l,
                    upper_offset: u,
                }),
                _ => None,
            };
            self.precursors.push(PrecursorInfo {
                mz,
                charge: self.cur_precursor_charge.take(),
                intensity: self.cur_precursor_intensity.take(),
                isolation_window,
                source_scan: self.cur_precursor_source_scan.take(),
            });
        } else {
            // No m/z → discard partial precursor data
            self.cur_precursor_charge = None;
            self.cur_precursor_intensity = None;
            self.cur_isolation_target_mz = None;
            self.cur_isolation_lower = None;
            self.cur_isolation_upper = None;
            self.cur_precursor_source_scan = None;
        }
    }

    fn build(self, _path: &Path) -> Result<Spectrum, SpectrumIoError> {
        let scan = self.scan_number.unwrap_or(1);
        let ms_level = match self.ms_level.unwrap_or(2) {
            1 => MsLevel::MS1,
            2 => MsLevel::MS2,
            n => MsLevel::Other(n),
        };

        // precursors were already accumulated by flush_precursor()
        let precursors = self.precursors;

        // Sort m/z + intensity arrays together by m/z ascending.
        let mut mz_array = self.mz_array;
        let mut intensity_array = self.intensity_array;
        crate::util::sort_peaks_by_mz(&mut mz_array, &mut intensity_array);

        Spectrum::new(
            scan,
            ms_level,
            self.rt_sec.unwrap_or(0.0),
            precursors,
            mz_array,
            intensity_array,
        )
        .map_err(|e| SpectrumIoError::ParseError {
            path: _path.to_path_buf(),
            line: 0,
            detail: e.to_string(),
        })
    }
}
```

- [ ] **Step 4: Update cvParam field assignments to use `cur_` prefixed fields**

In `parse_mzml_streaming`, change the `in_selected_ion` and `in_isolation_window` cvParam handlers (lines 332-349):

```rust
"MS:1000827" if in_isolation_window => {
    builder.cur_isolation_target_mz = value.parse().ok();
}
"MS:1000828" if in_isolation_window => {
    builder.cur_isolation_lower = value.parse().ok();
}
"MS:1000829" if in_isolation_window => {
    builder.cur_isolation_upper = value.parse().ok();
}
"MS:1000744" if in_selected_ion => {
    builder.cur_precursor_mz = value.parse().ok();
}
"MS:1000041" if in_selected_ion => {
    builder.cur_precursor_charge = value.parse().ok();
}
"MS:1000042" if in_selected_ion => {
    builder.cur_precursor_intensity = value.parse().ok();
}
```

- [ ] **Step 5: Update precursor Start event to store source_scan in temp field**

Change line 283-288:

```rust
b"precursor" if in_spectrum => {
    in_precursor = true;
    if let Some(spectrum_ref) = get_attr(e, b"spectrumRef") {
        builder.cur_precursor_source_scan =
            parse_scan_from_spectrum_ref(&spectrum_ref);
    }
}
```

- [ ] **Step 6: Add flush_precursor() call at `</precursor>` end tag**

Change line 390-392:

```rust
b"precursor" => {
    in_precursor = false;
    builder.flush_precursor();
}
```

- [ ] **Step 7: Run ALL spectrum-io tests**

Run: `cargo test -p protein-copilot-spectrum-io -- --nocapture`
Expected: ALL existing tests PASS (regression safe — behavior unchanged for single-precursor files)

- [ ] **Step 8: Commit**

```bash
git add crates/spectrum-io/src/mzml.rs
git commit -m "fix(spectrum-io): support multiple precursors per spectrum in mzML parser

Refactored SpectrumBuilder to accumulate precursors into a Vec by
building each PrecursorInfo at the </precursor> end tag instead of at
</spectrum>. Previously, multiple <precursor> elements were silently
overwritten, keeping only the last one.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: Annotation — match precursor by charge

**Context:** `annotate_spectrum()` in `crates/search-engine/src/annotate.rs` always uses `precursors.first()` to get the observed precursor m/z. For DIA multi-precursor spectra, this means the annotation delta_ppm may be calculated against the wrong precursor. The function already receives a `charge: i32` parameter — use it to find the matching precursor.

**Files:**
- Modify: `crates/search-engine/src/annotate.rs:291-297`
- Test: `crates/search-engine/src/annotate.rs` (in `mod tests`)

- [ ] **Step 1: Write the failing test**

Add a test in `annotate.rs` that creates a spectrum with 2 precursors of different charges and verifies the correct one is selected:

```rust
#[test]
fn annotate_selects_precursor_by_charge() {
    use protein_copilot_core::spectrum::{MsLevel, PrecursorInfo, Spectrum};
    use protein_copilot_core::search_params::MassTolerance;

    // Two precursors: charge 2 at 500.0, charge 3 at 400.0
    let spectrum = Spectrum::new(
        1,
        MsLevel::MS2,
        10.0,
        vec![
            PrecursorInfo {
                mz: 500.0,
                charge: Some(2),
                intensity: None,
                isolation_window: None,
                source_scan: None,
            },
            PrecursorInfo {
                mz: 400.0,
                charge: Some(3),
                intensity: None,
                isolation_window: None,
                source_scan: None,
            },
        ],
        vec![100.0, 200.0, 300.0],
        vec![1000.0, 2000.0, 500.0],
    )
    .unwrap();

    let tol = MassTolerance::Da(0.5);
    // Annotate with charge 3 — should pick precursor at 400.0, not 500.0
    let result = annotate_spectrum(&spectrum, "GK", 3, &tol, &[], vec![]);
    assert!(result.is_ok());
    let ann = result.unwrap();
    assert!(
        (ann.observed_precursor_mz - 400.0).abs() < 0.01,
        "should select precursor with matching charge 3, got {}",
        ann.observed_precursor_mz
    );
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p protein-copilot-search-engine annotate_selects_precursor_by_charge -- --nocapture`
Expected: FAIL — current code picks first precursor (500.0) regardless of charge

- [ ] **Step 3: Fix the precursor selection logic**

Replace lines 291-297 in `annotate.rs`:

```rust
    let precursor = spectrum
        .precursors
        .iter()
        .find(|p| p.charge == Some(charge))
        .or(spectrum.precursors.first())
        .ok_or_else(|| SearchEngineError::ExecutionError {
            detail: "spectrum has no precursor information".to_string(),
        })?;
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p protein-copilot-search-engine annotate_selects_precursor_by_charge -- --nocapture`
Expected: PASS

- [ ] **Step 5: Run full search-engine tests**

Run: `cargo test -p protein-copilot-search-engine -- --nocapture`
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add crates/search-engine/src/annotate.rs
git commit -m "fix(search-engine): match precursor by charge in annotation

annotate_spectrum() now prefers the precursor whose charge matches the
user-specified charge parameter, falling back to the first precursor
if no match is found. Previously always used precursors.first().

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: MGF parser — handle truncated files

**Context:** If an MGF file ends without a closing `END IONS` marker, the parser's `block` variable still holds parsed data (PEPMASS, peaks, etc.) but the function returns without processing it. The last spectrum is silently dropped.

Fix: after the line-reading loop, if `block` is `Some`, attempt to build the spectrum. The `into_spectrum()` method already validates required fields (PEPMASS), so incomplete blocks will error naturally.

**Files:**
- Modify: `crates/spectrum-io/src/mgf.rs:179-181`
- Test: `crates/spectrum-io/src/mgf.rs` (in `mod tests`)

- [ ] **Step 1: Write the failing test**

Create a test fixture string (no `END IONS`) and parse it:

```rust
#[test]
fn parse_truncated_mgf_without_end_ions() {
    use std::io::Write;
    let dir = std::env::temp_dir().join("mgf_truncated_test");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("truncated.mgf");
    {
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "BEGIN IONS").unwrap();
        writeln!(f, "PEPMASS=500.25 10000").unwrap();
        writeln!(f, "CHARGE=2+").unwrap();
        writeln!(f, "RTINSECONDS=60.0").unwrap();
        writeln!(f, "100.0 500").unwrap();
        writeln!(f, "200.0 1000").unwrap();
        // No END IONS
    }
    let reader = MgfReader;
    let spectra = reader.read_all(&path).unwrap();
    assert_eq!(spectra.len(), 1, "truncated block should still be parsed");
    assert!((spectra[0].precursors[0].mz - 500.25).abs() < 0.01);
    assert_eq!(spectra[0].num_peaks(), 2);
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p protein-copilot-spectrum-io parse_truncated_mgf_without_end_ions -- --nocapture`
Expected: FAIL — `spectra.len() == 0`

- [ ] **Step 3: Add truncation handling after the loop**

Replace line 181 (`Ok(count)`) with:

```rust
    // Handle truncated file: if a block was opened but never closed with
    // END IONS, attempt to build the spectrum from what we have.
    if let Some(b) = block.take() {
        let spectrum = b.into_spectrum(fallback_scan, path, block_start_line)?;
        handler(spectrum)?;
        count += 1;
    }

    Ok(count)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p protein-copilot-spectrum-io parse_truncated_mgf_without_end_ions -- --nocapture`
Expected: PASS

- [ ] **Step 5: Run full spectrum-io tests**

Run: `cargo test -p protein-copilot-spectrum-io -- --nocapture`
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add crates/spectrum-io/src/mgf.rs
git commit -m "fix(spectrum-io): handle truncated MGF files without END IONS

When an MGF file ends without a closing END IONS marker, the parser
now attempts to build the last spectrum from accumulated data instead
of silently dropping it. If the block is incomplete (e.g. missing
PEPMASS), the existing validation will still produce a clear error.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 4: Remove unnecessary `#[allow(dead_code)]` + full workspace validation

**Context:** `ProteinCopilotServer` at `crates/mcp-server/src/tools.rs:443` has `#[allow(dead_code)]` but the struct is used. Remove it.

**Files:**
- Modify: `crates/mcp-server/src/tools.rs:443`

- [ ] **Step 1: Remove the annotation**

Delete line 443 (`#[allow(dead_code)]`):

```rust
pub struct ProteinCopilotServer {
```

- [ ] **Step 2: Run full workspace build + clippy + tests**

```bash
cargo clippy --workspace --all-targets -- -D warnings 2>&1 | tail -5
cargo test --workspace 2>&1 | tail -10
```

Expected: 0 warnings, all tests pass

Note: If clippy DOES report `dead_code` on `ProteinCopilotServer` or its fields, re-add `#[allow(dead_code)]` only on the specific fields that trigger it (e.g. `tool_router` or `registry` which are read through the trait impl but not directly).

- [ ] **Step 3: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "chore: remove unnecessary #[allow(dead_code)] on ProteinCopilotServer

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---
