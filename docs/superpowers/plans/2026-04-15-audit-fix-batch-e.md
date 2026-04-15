# Audit Fix Batch E — Remaining Confirmed Findings

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the 6 remaining confirmed audit findings (C1/C2, H2-H4, H6) — terminal mod context in annotation, XSS HTML escaping, and duplicate scan handling.

**Architecture:** Three surgical fixes in separate crates: search-engine (annotation terminal mods), report (XSS escaping), spectrum-io (duplicate scans). Each fix is self-contained with tests.

**Tech Stack:** Rust, serde_json, spectrum-io indexed_mgf

---

## Summary of Findings

| ID     | Severity | Crate         | Issue |
|--------|----------|---------------|-------|
| C1/C2  | CRITICAL | search-engine | `apply_fixed_mod_mass()` in annotate.rs skips ProteinNTerm/CTerm mods — no terminal context passed |
| H2-H4  | HIGH     | report        | `escape_json_for_html()` only escapes `</script>` variants; doesn't escape `<`/`>` in JSON values |
| H6     | HIGH     | spectrum-io   | Duplicate scan numbers in MGF silently overwritten; later occurrence wins |

---

### Task 1: Add terminal context to annotation functions (C1/C2)

**Files:**
- Modify: `crates/search-engine/src/annotate.rs:209-231` (apply_fixed_mod_mass)
- Modify: `crates/search-engine/src/annotate.rs:325-332` (annotate_spectrum signature)
- Modify: `crates/search-engine/src/annotate.rs:505-511` (annotate_heavy_spectrum signature)
- Modify: `crates/mcp-server/src/tools.rs` (callers of annotate_spectrum/annotate_heavy_spectrum)
- Test: `crates/search-engine/src/annotate.rs` (existing tests + new)

**Context:** `matching.rs:apply_fixed_mods` correctly accepts `is_protein_nterm: bool, is_protein_cterm: bool` and gates ProteinNTerm/CTerm mods. But `annotate.rs:apply_fixed_mod_mass` always skips them with comment "no terminal context available". The fix adds the same parameters.

- [ ] **Step 1: Add terminal params to `apply_fixed_mod_mass`**

In `crates/search-engine/src/annotate.rs`, change the function signature and body:

```rust
fn apply_fixed_mod_mass(
    sequence: &str,
    fixed_mods: &[Modification],
    is_protein_nterm: bool,
    is_protein_cterm: bool,
) -> f64 {
    use protein_copilot_core::search_params::ModPosition;
    let mut delta = 0.0;
    for m in fixed_mods {
        if m.residues.is_empty() {
            match m.position {
                ModPosition::AnyNTerm | ModPosition::AnyCTerm | ModPosition::Anywhere => {
                    delta += m.mass_delta;
                }
                ModPosition::ProteinNTerm => {
                    if is_protein_nterm {
                        delta += m.mass_delta;
                    }
                }
                ModPosition::ProteinCTerm => {
                    if is_protein_cterm {
                        delta += m.mass_delta;
                    }
                }
            }
        } else {
            for ch in sequence.chars() {
                if m.residues.contains(&ch) {
                    delta += m.mass_delta;
                }
            }
        }
    }
    delta
}
```

- [ ] **Step 2: Update `annotate_spectrum` signature**

Add `is_protein_nterm: bool, is_protein_cterm: bool` parameters after `protein_accessions`:

```rust
pub fn annotate_spectrum(
    spectrum: &Spectrum,
    peptide_sequence: &str,
    charge: i32,
    fragment_tolerance: &MassTolerance,
    fixed_modifications: &[Modification],
    protein_accessions: Vec<String>,
    is_protein_nterm: bool,
    is_protein_cterm: bool,
) -> Result<SpectrumAnnotation, SearchEngineError> {
```

Update the call to `apply_fixed_mod_mass` at line ~361:
```rust
let mod_delta = apply_fixed_mod_mass(peptide_sequence, fixed_modifications, is_protein_nterm, is_protein_cterm);
```

- [ ] **Step 3: Update `annotate_heavy_spectrum` signature**

Same pattern — add the two booleans, update the `apply_fixed_mod_mass` call at line ~532:

```rust
pub fn annotate_heavy_spectrum(
    heavy_spectrum: &Spectrum,
    peptide_sequence: &str,
    charge: i32,
    fragment_tolerance: &MassTolerance,
    fixed_modifications: &[Modification],
    label: &protein_copilot_core::label::LabelType,
    is_protein_nterm: bool,
    is_protein_cterm: bool,
) -> Result<HeavyAnnotation, SearchEngineError> {
```

- [ ] **Step 4: Update all callers in tools.rs**

In `crates/mcp-server/src/tools.rs`, every call to `annotate_spectrum` and `annotate_heavy_spectrum` needs `false, false` appended (conservative default — annotation mode doesn't know terminal position):

Search for `annotate_spectrum(` and `annotate_heavy_spectrum(` and add `, false, false` before the closing `)`.

- [ ] **Step 5: Update existing tests in annotate.rs**

All test calls to `annotate_spectrum` and `annotate_heavy_spectrum` need the new params. Add `false, false` to each existing test call.

- [ ] **Step 6: Add test for terminal mod application**

```rust
#[test]
fn test_terminal_mod_applied_when_context_true() {
    use protein_copilot_core::search_params::{Modification, ModPosition, MassTolerance, ToleranceUnit};
    let nterm_mod = Modification {
        name: "Acetyl".to_string(),
        mass_delta: 42.010565,
        residues: vec![],
        position: ModPosition::ProteinNTerm,
    };
    // With is_protein_nterm=true, mod should be applied
    let delta_applied = apply_fixed_mod_mass("PEPTIDE", &[nterm_mod.clone()], true, false);
    assert!((delta_applied - 42.010565).abs() < 1e-6);

    // With is_protein_nterm=false, mod should be skipped
    let delta_skipped = apply_fixed_mod_mass("PEPTIDE", &[nterm_mod], false, false);
    assert!(delta_skipped.abs() < 1e-6);
}
```

- [ ] **Step 7: Build and test**

Run: `cargo test -p protein-copilot-search-engine && cargo test -p protein-copilot-mcp-server --lib`
Expected: All tests pass

- [ ] **Step 8: Commit**

```bash
git add crates/search-engine/src/annotate.rs crates/mcp-server/src/tools.rs
git commit -m "fix(C1/C2): add terminal mod context to annotation functions

apply_fixed_mod_mass now accepts is_protein_nterm/is_protein_cterm
to correctly gate ProteinNTerm/CTerm fixed modifications.
Callers default to (false, false) — conservative for annotation mode."
```

---

### Task 2: Robust XSS escaping for HTML-embedded JSON (H2-H4)

**Files:**
- Modify: `crates/report/src/lib.rs:22-27` (escape_json_for_html)
- Test: `crates/report/src/lib.rs` (add tests)

**Context:** The current `escape_json_for_html` only replaces `</script>`, `</Script>`, `</SCRIPT>`. This misses case-insensitive variants and doesn't escape raw `<`/`>` that could break out of script context. The standard approach for embedding JSON in HTML `<script>` tags is to escape all `<` as `\u003c` and `>` as `\u003e`.

- [ ] **Step 1: Write failing test**

Add to bottom of `crates/report/src/lib.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn escape_script_close_tag() {
        let input = r#"{"val":"</script><b>"}"#;
        let escaped = escape_json_for_html(input);
        assert!(!escaped.contains("</script>"), "must not contain literal </script>");
        assert!(!escaped.contains('<'), "must not contain literal <");
        assert!(!escaped.contains('>'), "must not contain literal >");
    }

    #[test]
    fn escape_angle_brackets() {
        let input = r#"{"key":"<img onerror=alert(1)>"}"#;
        let escaped = escape_json_for_html(input);
        assert!(!escaped.contains('<'));
        assert!(!escaped.contains('>'));
        // Must still be valid JSON when unescaped
        assert!(escaped.contains(r"\u003c") || escaped.contains(r"\u003C"));
    }

    #[test]
    fn escape_preserves_normal_json() {
        let input = r#"{"score":0.95,"peptide":"ACDK"}"#;
        let escaped = escape_json_for_html(input);
        assert_eq!(input, escaped, "no angle brackets = no changes");
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p protein-copilot-report -- tests::escape`
Expected: `escape_script_close_tag` and `escape_angle_brackets` FAIL

- [ ] **Step 3: Implement robust escaping**

Replace the function in `crates/report/src/lib.rs`:

```rust
/// Escapes JSON for safe embedding inside HTML `<script>` tags.
///
/// Replaces `<` with `\u003c` and `>` with `\u003e` to prevent:
/// - Premature `</script>` tag closure
/// - HTML injection via `<` and `>` in JSON string values
///
/// The escaped string remains valid JSON (parseable by `JSON.parse()`).
pub fn escape_json_for_html(json: &str) -> String {
    json.replace('<', r"\u003c").replace('>', r"\u003e")
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p protein-copilot-report -- tests::escape`
Expected: All 3 tests PASS

- [ ] **Step 5: Run full report tests**

Run: `cargo test -p protein-copilot-report`
Expected: All tests pass (the HTML render tests should still work since `\u003c` is valid JSON)

- [ ] **Step 6: Commit**

```bash
git add crates/report/src/lib.rs
git commit -m "fix(H2-H4): robust XSS escaping for HTML-embedded JSON

Replace ad-hoc </script> escaping with universal < > escaping
using JSON unicode escapes (\u003c, \u003e). Prevents all HTML
injection vectors while preserving valid JSON."
```

---

### Task 3: Warn-and-skip duplicate scan numbers in MGF (H6)

**Files:**
- Modify: `crates/spectrum-io/src/indexed_mgf.rs:146-153`
- Test: `crates/spectrum-io/src/indexed_mgf.rs` (add test)

**Context:** When an MGF file has duplicate SCANS= values, the current code warns but overwrites the earlier entry. This means the first spectrum with that scan number becomes inaccessible. Better behavior: keep the first occurrence and skip subsequent duplicates.

- [ ] **Step 1: Write failing test**

Add to the test module in `crates/spectrum-io/src/indexed_mgf.rs`:

```rust
#[test]
fn duplicate_scan_keeps_first_occurrence() {
    use std::io::Write;
    let dir = tempfile::tempdir().unwrap();
    let mgf_path = dir.path().join("dup.mgf");
    let mut f = std::fs::File::create(&mgf_path).unwrap();
    // Two spectra with same SCANS=5
    write!(f, "BEGIN IONS\nSCANS=5\nPEPMASS=500.0\n100.0 1000\nEND IONS\n").unwrap();
    write!(f, "BEGIN IONS\nSCANS=5\nPEPMASS=600.0\n200.0 2000\nEND IONS\n").unwrap();
    drop(f);

    let index = build_scan_index(&mgf_path).unwrap();
    // Should keep the FIRST occurrence (PEPMASS=500)
    let reader = IndexedMgfReader;
    let spec = reader.read_spectrum(&mgf_path, 5).unwrap();
    assert!(spec.precursors[0].mz > 499.0 && spec.precursors[0].mz < 501.0,
        "should keep first occurrence with PEPMASS=500, got {}", spec.precursors[0].mz);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p protein-copilot-spectrum-io -- duplicate_scan_keeps_first`
Expected: FAIL (currently keeps second occurrence with PEPMASS=600)

- [ ] **Step 3: Fix — keep first, skip subsequent**

In `crates/spectrum-io/src/indexed_mgf.rs`, change lines 146-153:

```rust
        if offsets.contains_key(&scan_num) {
            tracing::warn!(
                "duplicate scan number {} in MGF file {:?}; keeping first occurrence, skipping later",
                scan_num,
                path,
            );
            // Keep the first occurrence — do NOT overwrite
            continue;
        }
        offsets.insert(scan_num, offset);
```

The key change: add `continue;` after the warning instead of falling through to `offsets.insert`.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p protein-copilot-spectrum-io -- duplicate_scan`
Expected: PASS

- [ ] **Step 5: Run full spectrum-io tests**

Run: `cargo test -p protein-copilot-spectrum-io`
Expected: All tests pass

- [ ] **Step 6: Commit**

```bash
git add crates/spectrum-io/src/indexed_mgf.rs
git commit -m "fix(H6): keep first occurrence on duplicate MGF scan numbers

Previously, duplicate SCANS= entries silently overwrote earlier ones.
Now keeps the first occurrence and skips later duplicates with a warning."
```

---

### Task 4: Final verification

- [ ] **Step 1: Full workspace test**

Run: `cargo test --workspace`
Expected: All tests pass (should be ~550+)

- [ ] **Step 2: Clippy clean**

Run: `cargo clippy --all-targets -- -D warnings`
Expected: No warnings

- [ ] **Step 3: Update audit_findings status**

Mark C1, C2, H2, H3, H4, H6 as `fixed` in SQL tracking.

- [ ] **Step 4: Final commit (if any remaining changes)**

```bash
git add -A && git commit -m "chore: audit batch E complete — all confirmed findings fixed"
```
