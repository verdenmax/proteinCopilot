# Defensive Programming Fixes — Design Spec

> **Date:** 2026-04-28
> **Scope:** 9 fixes across 6 crates — division-by-zero guards, memory limits, data validation, error reporting

---

## Problem Statement

A comprehensive 4-agent defensive programming audit identified 9 issues in production code:
- 4 critical (🔴): crashes or DoS from user input
- 5 medium (🟡): silent wrong results or poor error messages

All issues are fixable with targeted guards at module boundaries. No architectural changes needed.

---

## Fixes

### D1 🔴 Sage adapter charge=0 division-by-zero

**File:** `crates/search-engine/src/adapters/sage/mod.rs:553`

**Problem:** `(mass + charge * PROTON) / charge` — if sage-core returns `Feature` with `charge=0`, result is `Infinity`. Downstream PSM gets corrupt `calculated_mz` and `delta_mass_ppm`.

**Fix:** Guard before division. If `charge == 0`, log a warning and skip (return `None` from the mapping, filter out).

```rust
let charge = feat.charge;
if charge == 0 {
    tracing::warn!(scan = feat.scan, "sage returned charge=0, skipping PSM");
    return None; // caller uses filter_map
}
```

### D2 🔴 MCP annotate_spectrum charge=0 division-by-zero

**File:** `crates/mcp-server/src/tools.rs:2322`

**Problem:** `(mass + charge * PROTON_MASS) / charge as f64` — user can pass `charge=0` via MCP tool input.

**Fix:** Add charge validation at tool entry point, before any calculation. Return user-friendly error.

```rust
if charge <= 0 {
    return Err("charge must be positive (≥1)".into());
}
```

### D3 🔴 FASTA read_to_string no size limit

**File:** `crates/mcp-server/src/tools.rs:927` (`load_fasta_sequences`)

**Problem:** `std::fs::read_to_string(fasta_path)` with no size check. 2GB+ FASTA = OOM crash.

**Fix:** Check file size before reading. Limit: **512 MB**.

```rust
let metadata = std::fs::metadata(fasta_path)?;
if metadata.len() > 512 * 1024 * 1024 {
    return Err(format!(
        "FASTA file too large ({:.0} MB, max 512 MB): {fasta_path}",
        metadata.len() as f64 / 1_048_576.0
    ));
}
```

### D4 🔴 mzML binary array unbounded allocation

**File:** `crates/spectrum-io/src/mzml.rs` (decode_binary_array / SpectrumBuilder)

**Problem:** Base64-decoded peak arrays have no size limit. Malformed mzML declaring huge `defaultArrayLength` can trigger gigabyte allocations.

**Fix:** After decoding, validate array length ≤ **500,000** peaks. Also validate mz and intensity arrays have equal length.

```rust
const MAX_PEAKS_PER_SPECTRUM: usize = 500_000;
if mz_array.len() > MAX_PEAKS_PER_SPECTRUM {
    return Err(SpectrumIoError::FormatError {
        detail: format!("spectrum has {} peaks (max {})", mz_array.len(), MAX_PEAKS_PER_SPECTRUM),
    });
}
if mz_array.len() != intensity_array.len() {
    return Err(SpectrumIoError::FormatError {
        detail: format!("mz ({}) and intensity ({}) array length mismatch",
            mz_array.len(), intensity_array.len()),
    });
}
```

### D5 🟡 Spectrum mz/intensity length mismatch validation

**File:** `crates/core/src/spectrum.rs`

**Problem:** `Spectrum` struct accepts any `mz_array` and `intensity_array` without length check. Mismatched lengths cause silent data corruption.

**Fix:** Add a `validate()` method (or validate in builder) that checks `mz_array.len() == intensity_array.len()`. mzML parser (D4) calls this. MGF parser also calls this.

### D6 🟡 Empty mzML (0 spectra) not rejected early

**File:** `crates/spectrum-io/src/util.rs` (`SummaryAccumulator::into_summary`)

**Problem:** Empty file produces `SpectrumSummary` with `total_spectra=0`, `mz_range=[0,0]`. Downstream search silently produces no results.

**Fix:** Add warning log when total_spectra == 0. Don't reject (some tools legitimately read empty files for metadata), but warn clearly.

```rust
if self.total == 0 {
    tracing::warn!(path = %path.display(), "spectrum file contains 0 spectra");
}
```

### D7 🟡 Generic TSV header case-sensitive matching without warning

**File:** `crates/entrapment-analysis/src/loader/generic_tsv.rs:128-133`

**Problem:** `find_column(headers, &column_map.charge)` is case-sensitive. If external TSV uses `Charge` instead of `charge`, column is silently missed. All values become `None`.

**Fix:** Add `tracing::debug!` for each optional column not found. This doesn't change behavior but gives users visibility into why data is missing.

```rust
let charge_idx = find_column(headers, &column_map.charge);
if charge_idx.is_none() {
    tracing::debug!(column = %column_map.charge, "optional column not found in TSV headers");
}
```

### D8 🟡 Missing RT silently defaults to 0.0

**File:** `crates/spectrum-io/src/mzml.rs` (SpectrumBuilder)

**Problem:** When `scanStartTime` is missing from mzML, RT defaults to 0.0 without warning. XIC extraction by RT will match wrong scans.

**Fix:** Add warning when RT is None for MS2 spectra (MS1 missing RT is less critical).

```rust
if ms_level >= 2 && rt_sec.is_none() {
    tracing::warn!(scan = scan_number, "MS2 spectrum missing retention time, defaulting to 0.0");
}
```

### D9 🟡 DIA spill directory write failure silent

**File:** `crates/mcp-server/src/tools.rs:1060`

**Problem:** When DIA spill directory can't be created, code silently keeps everything in memory. Large DIA runs OOM without user knowing why.

**Fix:** Upgrade from `tracing::warn!` to `tracing::error!` and include actionable suggestion.

```rust
tracing::error!(
    dir = %self.spill_dir.display(),
    "Cannot create DIA spill directory: {e}. \
     All spectra kept in memory — large DIA runs may OOM. \
     Fix: ensure directory is writable or set a different spill_dir."
);
```

---

## Design Decisions

1. **Guard, don't crash**: All D1-D4 fixes add validation guards. No panics.
2. **Warn, don't reject** for D6-D9: These are data quality issues, not crashes. Warn at appropriate log level to help debugging.
3. **Constants**: `MAX_FASTA_SIZE = 512 MB`, `MAX_PEAKS_PER_SPECTRUM = 500,000` — defined as module-level constants with doc comments.
4. **No behavior changes**: All fixes are additive guards. Existing valid inputs produce identical results.
5. **Test code excluded**: `entrapment_integration.rs` has unsafe indexing but is test-only — not in scope.

---

## Testing Strategy

- **D1**: Unit test with mock Feature(charge=0) → verify PSM skipped
- **D2**: Integration test with charge=0 MCP input → verify error message
- **D3**: Unit test with oversized FASTA metadata → verify rejection
- **D4**: Unit test with array > 500k → verify error
- **D5**: Unit test with mismatched array lengths → verify error
- **D6-D9**: Verify via `cargo test` that existing tests still pass (warn-only changes)
