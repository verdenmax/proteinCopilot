# Integration Tests + Heavy Zero-Offset Validation

**Date**: 2026-04-14
**Status**: Draft
**Scope**: Two focused improvements — workspace-level integration tests and SILAC zero-offset guard

---

## 1. Problem Statement

### 1.1 Missing Integration Tests
All 296 tests are unit-level or crate-internal. There are no workspace-level tests that
exercise cross-crate data flows (e.g., annotation → XIC pipeline, SILAC 4-scenario matrix).
The recently fixed DDA+SILAC bug (scenario ②) had no test to catch it.

### 1.2 Heavy Zero-Offset Bug
When a peptide has no K or R residues (e.g., "PEPTIDE"), SILAC heavy_mz == light_mz.
Current code still attempts heavy scan lookup and XIC extraction, which:
- Wastes computation scanning ±50 nearby spectra
- In DDA mode, may accidentally match the *same* scan as "heavy" (precursor mz matches within tolerance)
- Produces a meaningless mirror plot (identical light and heavy)

---

## 2. Design: Integration Tests

### 2.1 Location
Root directory `tests/` — workspace-level integration tests that import multiple crates.

### 2.2 Test Files

| File | Coverage | Scenarios |
|------|----------|-----------|
| `tests/annotation_scenarios.rs` | 4-scenario annotation matrix | DDA×non-SILAC, DDA×SILAC, DIA×non-SILAC, DIA×SILAC, zero-offset skip |
| `tests/xic_scenarios.rs` | 4-scenario XIC matrix | Same 4 + zero-offset skip + MS1 heavy extraction |
| `tests/search_pipeline.rs` | End-to-end search | read → recommend → search → FDR → summary → annotate |
| `tests/import_pipeline.rs` | External result import | Custom JSON → import → summary → annotate |

### 2.3 Test Data Strategy
**Synthetic Spectrum construction** — no external files for SILAC/DIA tests.

Helper functions build `Spectrum` structs with:
- Controlled m/z arrays (known peaks at theoretical fragment positions)
- Controlled precursor info (m/z, charge, isolation window)
- DDA: narrow window (0.7 Th), DIA: wide window (25 Th)
- SILAC: light scan + heavy scan at shifted precursor m/z

For search pipeline tests, reuse existing `test_100.mgf` + `test_100.fasta` fixtures.

### 2.4 Assertion Strategy
- Annotation: check `score > 0`, `matched_ions > 0`, `heavy_annotation.is_some()` for SILAC, `.is_none()` for non-SILAC
- XIC: check trace count > 0, heavy traces present/absent per scenario
- Zero-offset: `heavy_annotation.is_none()` even when `label_type` provided
- Search pipeline: `psms_at_1pct_fdr > 0`, summary fields non-zero

---

## 3. Design: Heavy Zero-Offset Validation

### 3.1 Check Logic
```rust
let delta = total_heavy_delta(peptide_sequence, label);
if delta.abs() < 1e-6 {
    // Skip: no K/R residues, heavy == light
}
```

K/R count is integer × constant (8.014 or 10.008), so exact zero check with f64 epsilon is safe.

### 3.2 Check Points (3 locations)

1. **`tools.rs` annotation handler** — After computing `heavy_prec_mz`, check if shift is zero.
   If zero: skip heavy scan lookup entirely, log info message, no mirror plot.

2. **`xic/extract.rs`** — After computing `heavy_ions`, check if `heavy_precursor_mz == precursor_mz`
   (within tolerance). If so: clear `heavy_ions` to empty, skip heavy XIC extraction.

3. **`heavy.rs` `compute_heavy_target_ions()`** — No change needed here (pure function,
   caller decides whether to use result).

### 3.3 Behavior When Skipped
- `annotation.heavy_annotation` remains `None`
- XIC: `heavy_fragment_xic_traces` = empty, `ms1_heavy_points` = empty
- HTML: renders as standard (non-mirror) annotation
- Log: `tracing::info!("Skipping heavy annotation: peptide has no K/R, zero SILAC shift")`
- No error, no warning — this is expected behavior for non-labeled peptides in a SILAC experiment

---

## 4. Scope Exclusions
- No new external test fixture files
- No changes to HTML template rendering logic
- No changes to `annotate_heavy_spectrum()` itself (pure function, unchanged)
- No changes to MCP tool schema (label_type is still accepted, just produces standard output)
