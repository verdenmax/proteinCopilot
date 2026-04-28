# Defensive Programming Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add defensive guards for 8 issues found by audit: division-by-zero, memory limits, data validation, and error reporting improvements.

**Architecture:** Targeted guards at module boundaries. No API changes — only adding validation before dangerous operations and improving log messages. D5 (mz/intensity mismatch) is already handled by `Spectrum::new()` → `validate()` and is excluded.

**Tech Stack:** Rust, tracing crate for logging, existing error types

---

## File Map

| File | Changes |
|------|---------|
| `crates/search-engine/src/adapters/sage/mod.rs` | D1: charge=0 guard in `feature_to_psm` |
| `crates/mcp-server/src/tools.rs` | D2: charge guard after mode resolution; D3: FASTA size limit; D9: DIA spill error upgrade |
| `crates/spectrum-io/src/mzml.rs` | D4: peak array size limit; D8: missing RT warning |
| `crates/spectrum-io/src/util.rs` | D6: empty file warning |
| `crates/entrapment-analysis/src/loader/generic_tsv.rs` | D7: TSV header debug logging |

---

### Task 1: D1 — Sage adapter charge=0 guard

**Files:**
- Modify: `crates/search-engine/src/adapters/sage/mod.rs:543-554` (function signature)
- Modify: `crates/search-engine/src/adapters/sage/mod.rs:346-349` (call site)

- [ ] **Step 1: Change `feature_to_psm` to return `Option<Psm>` and guard charge=0**

In `crates/search-engine/src/adapters/sage/mod.rs`, change the function signature and add guard:

```rust
// Line 543: change return type
fn feature_to_psm(feat: &Feature, db: &IndexedDatabase) -> Option<Psm> {
    let peptide = &db[feat.peptide_idx];

    let peptide_sequence = String::from_utf8_lossy(&peptide.sequence).to_string();
    let protein_accessions: Vec<String> = peptide.proteins.iter().map(|p| p.to_string()).collect();

    // Guard: charge=0 produces Infinity in m/z calculation
    let charge = f64::from(feat.charge);
    if charge == 0.0 {
        tracing::warn!(
            spec_id = %feat.spec_id,
            peptide = %peptide_sequence,
            "sage returned charge=0, skipping PSM"
        );
        return None;
    }
    let proton = f64::from(PROTON);
    let precursor_mz = (f64::from(feat.expmass) + charge * proton) / charge;
    let calculated_mz = (f64::from(feat.calcmass) + charge * proton) / charge;
    // ... rest unchanged, but wrap final return in Some(Psm { ... })
```

At the end of the function, change `Psm { ... }` to `Some(Psm { ... })`.

- [ ] **Step 2: Update call site to use `filter_map`**

In the same file at line 346-349, change `.map()` to `.filter_map()`:

```rust
let psms: Vec<Psm> = features
    .iter()
    .filter_map(|feat| feature_to_psm(feat, &db))
    .collect();
```

- [ ] **Step 3: Run tests**

```bash
cd /home/verden/pfind/2026-spring/code/proteinCopilot
cargo test -p protein-copilot-search-engine --lib 2>&1 | tail -5
```

Expected: All existing tests pass (charge=0 never occurs in test fixtures).

- [ ] **Step 4: Commit**

```bash
git add crates/search-engine/src/adapters/sage/mod.rs
git commit -m "fix(search-engine): guard charge=0 division-by-zero in sage adapter (D1)

feature_to_psm now returns Option<Psm> and skips PSMs with charge=0.
Prevents Infinity in precursor_mz and calculated_mz fields.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: D2 + D3 + D9 — MCP Server defensive guards

**Files:**
- Modify: `crates/mcp-server/src/tools.rs:2270-2322` (D2: charge guard)
- Modify: `crates/mcp-server/src/tools.rs:925-928` (D3: FASTA size limit)
- Modify: `crates/mcp-server/src/tools.rs:1058-1063` (D9: DIA spill error)

- [ ] **Step 1: D2 — Add charge > 0 validation after mode resolution**

In `crates/mcp-server/src/tools.rs`, after the mode 1 / mode 2 resolution block (after line 2309), add a defense-in-depth charge check before the RT auto-lookup code:

```rust
        // Defense-in-depth: validate charge > 0 before any m/z calculation.
        // Mode 2 already validates at line 2291; this catches mode 1 PSMs with charge=0.
        if charge <= 0 {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!("charge must be > 0 (got {charge}); PSM may have invalid charge state"),
            ));
        }
```

Insert this between line 2309 (end of mode resolution) and line 2311 (start of scan resolution).

- [ ] **Step 2: D3 — Add FASTA file size limit**

In `crates/mcp-server/src/tools.rs`, at the beginning of `load_fasta_sequences` function (line 927), add size check before `read_to_string`:

```rust
/// Maximum FASTA file size for in-memory loading (512 MB).
const MAX_FASTA_SIZE: u64 = 512 * 1024 * 1024;

fn load_fasta_sequences(fasta_path: &str) -> Result<HashMap<String, String>, String> {
    // Guard: reject FASTA files larger than 512 MB to prevent OOM
    let metadata = std::fs::metadata(fasta_path)
        .map_err(|e| format!("Failed to read FASTA file metadata {fasta_path}: {e}"))?;
    if metadata.len() > MAX_FASTA_SIZE {
        return Err(format!(
            "FASTA file too large ({:.0} MB, max {} MB): {fasta_path}",
            metadata.len() as f64 / 1_048_576.0,
            MAX_FASTA_SIZE / 1_048_576,
        ));
    }
    let content = std::fs::read_to_string(fasta_path)
        .map_err(|e| format!("Failed to read FASTA file {fasta_path}: {e}"))?;
    // ... rest unchanged
```

- [ ] **Step 3: D9 — Upgrade DIA spill directory error logging**

In `crates/mcp-server/src/tools.rs` at line 1061, change `tracing::warn!` to `tracing::error!` with actionable message:

```rust
fn spill_to_disk(&self, id: Uuid, spectra: &[Spectrum]) -> bool {
    if let Err(e) = std::fs::create_dir_all(&self.spill_dir) {
        tracing::error!(
            dir = %self.spill_dir.display(),
            error = %e,
            "Cannot create DIA spill directory — all spectra kept in memory. \
             Large DIA runs may OOM. Ensure directory is writable."
        );
        return false;
    }
```

- [ ] **Step 4: Run tests**

```bash
cargo test -p protein-copilot-mcp-server --lib 2>&1 | tail -5
```

Expected: All existing tests pass.

- [ ] **Step 5: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "fix(mcp-server): defensive guards — charge validation, FASTA size limit, spill logging (D2/D3/D9)

D2: defense-in-depth charge>0 check before m/z calc in annotate_spectrum
D3: reject FASTA files >512 MB in load_fasta_sequences to prevent OOM
D9: upgrade DIA spill dir failure from warn to error with actionable message

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: D4 + D6 + D8 — spectrum-io defensive guards

**Files:**
- Modify: `crates/spectrum-io/src/mzml.rs:139-167` (D4: peak array limit in `decode_binary_array`)
- Modify: `crates/spectrum-io/src/mzml.rs:107-132` (D8: RT warning in `SpectrumBuilder::build`)
- Modify: `crates/spectrum-io/src/util.rs:144-150` (D6: empty file warning)

- [ ] **Step 1: D4 — Add peak array size limit in `decode_binary_array`**

In `crates/spectrum-io/src/mzml.rs`, add a constant at the top of the file (after imports) and a check at the end of `decode_binary_array` (just before returning):

```rust
/// Maximum peaks per spectrum. Prevents OOM from malformed mzML files
/// declaring extremely large arrays. 500k covers wide DIA windows.
const MAX_PEAKS_PER_SPECTRUM: usize = 500_000;
```

At the end of `decode_binary_array` (before `Ok(values)` for both 32-bit and 64-bit branches), add a check. The cleanest place is after both branches, checking the result. Restructure the function ending:

Replace the two branch endings (lines 181-187 for 64-bit, lines 200-206 for 32-bit) — keep the existing decode logic, but add a check after:

```rust
    // After the if/else block that produces values:
    let values: Vec<f64> = if meta.is_64bit {
        // ... existing 64-bit decode (lines 170-187)
    } else {
        // ... existing 32-bit decode (lines 189-207)
    };

    if values.len() > MAX_PEAKS_PER_SPECTRUM {
        return Err(SpectrumIoError::BinaryDecodeError {
            path: path.to_path_buf(),
            detail: format!(
                "array has {} elements (max {MAX_PEAKS_PER_SPECTRUM}); \
                 file may be corrupt",
                values.len()
            ),
        });
    }

    Ok(values)
```

Note: The existing code already returns `Ok(bytes.chunks_exact(8)...collect())` directly from each branch. You need to capture the result into a `let values = ...;` variable, then check, then return.

- [ ] **Step 2: D8 — Add missing RT warning in SpectrumBuilder::build**

In `crates/spectrum-io/src/mzml.rs`, in `SpectrumBuilder::build()` (line 107), add a warning when RT is None for MS2+ spectra. Add this before the `Spectrum::new()` call:

```rust
    fn build(self, _path: &Path) -> Result<Spectrum, SpectrumIoError> {
        let scan = self.scan_number.unwrap_or(1);
        let ms_level = match self.ms_level.unwrap_or(2) {
            1 => MsLevel::MS1,
            2 => MsLevel::MS2,
            n => MsLevel::Other(n),
        };

        // Warn when MS2+ spectrum has no retention time — RT-based lookups
        // (XIC extraction, scan auto-matching) will use 0.0 as fallback.
        if self.rt_min.is_none() && !matches!(ms_level, MsLevel::MS1) {
            tracing::debug!(
                scan,
                "spectrum missing retention time, defaulting to 0.0"
            );
        }

        // Sort m/z + intensity arrays together by m/z ascending.
        // ... rest unchanged
```

Use `tracing::debug!` (not warn) because some test/synthetic mzML files legitimately lack RT.

- [ ] **Step 3: D6 — Add empty file warning**

In `crates/spectrum-io/src/util.rs`, inside `into_summary()` at line 145, add a warning:

```rust
        if self.total == 0 {
            tracing::warn!(
                path = %path.display(),
                "spectrum file contains 0 spectra — downstream analysis will produce no results"
            );
            self.mz_min = 0.0;
            self.mz_max = 0.0;
            self.rt_min = 0.0;
            self.rt_max = 0.0;
        }
```

- [ ] **Step 4: Run tests**

```bash
cargo test -p protein-copilot-spectrum-io --lib 2>&1 | tail -10
```

Expected: All existing tests pass (6 pre-existing byte-scan failures are unrelated).

- [ ] **Step 5: Commit**

```bash
git add crates/spectrum-io/src/mzml.rs crates/spectrum-io/src/util.rs
git commit -m "fix(spectrum-io): peak array limit, empty file warning, missing RT warning (D4/D6/D8)

D4: reject mzML spectra with >500k peaks to prevent OOM from malformed files
D6: warn when spectrum file contains 0 spectra
D8: debug-log when MS2 spectrum missing retention time (defaults to 0.0)

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 4: D7 — TSV header debug logging

**Files:**
- Modify: `crates/entrapment-analysis/src/loader/generic_tsv.rs:127-133`

- [ ] **Step 1: Add debug logging for unfound optional columns**

In `crates/entrapment-analysis/src/loader/generic_tsv.rs`, after lines 127-133 where optional columns are looked up, add debug logging for each unfound column:

```rust
    let charge_idx = find_column(headers, &column_map.charge);
    let precursor_mz_idx = find_column(headers, &column_map.precursor_mz);
    let retention_time_idx = find_column(headers, &column_map.retention_time);
    let scan_number_idx = find_column(headers, &column_map.scan_number);
    let spectrum_file_idx = find_column(headers, &column_map.spectrum_file);
    let protein_ids_idx = find_column(headers, &column_map.protein_ids);
    let q_value_idx = find_column(headers, &column_map.q_value);

    // Log unfound optional columns to aid debugging when data appears missing
    for (name, idx) in [
        (&column_map.charge, &charge_idx),
        (&column_map.precursor_mz, &precursor_mz_idx),
        (&column_map.retention_time, &retention_time_idx),
        (&column_map.scan_number, &scan_number_idx),
        (&column_map.spectrum_file, &spectrum_file_idx),
        (&column_map.protein_ids, &protein_ids_idx),
        (&column_map.q_value, &q_value_idx),
    ] {
        if idx.is_none() {
            tracing::debug!(
                column = %name,
                file = %path.display(),
                "optional TSV column not found in headers"
            );
        }
    }
```

- [ ] **Step 2: Run tests**

```bash
cargo test -p protein-copilot-entrapment-analysis --lib 2>&1 | tail -5
```

Expected: All existing tests pass.

- [ ] **Step 3: Commit**

```bash
git add crates/entrapment-analysis/src/loader/generic_tsv.rs
git commit -m "fix(entrapment): log unfound optional TSV columns for debugging (D7)

Adds tracing::debug for each optional column not found in TSV headers.
Helps diagnose silent data loss from case-sensitive header mismatches.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 5: Verify all crates compile and test

- [ ] **Step 1: Full workspace build**

```bash
cargo build --workspace 2>&1 | tail -5
```

Expected: No compilation errors.

- [ ] **Step 2: Full workspace test**

```bash
cargo test --workspace 2>&1 | tail -20
```

Expected: All tests pass (except 6 pre-existing spectrum-io byte-scan failures).

- [ ] **Step 3: Clippy check**

```bash
cargo clippy --workspace 2>&1 | tail -10
```

Expected: No new warnings.
