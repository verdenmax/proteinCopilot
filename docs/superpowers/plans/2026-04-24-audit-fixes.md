# Audit Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix all high/medium priority findings from the 5-agent codebase audit (error handling, HTML templates, MCP tools, API consistency).

**Architecture:** Surgical fixes across multiple crates — no structural changes. Each task is independent and can be committed separately. All fixes are backward-compatible.

**Tech Stack:** Rust, HTML/JS templates, serde

---

## File Map

| File | Action | Responsibility |
|------|--------|----------------|
| `crates/mcp-server/src/tools.rs` | Modify | Fix unwrap L1765, extract magic numbers |
| `crates/result-import/src/unimod.rs` | Modify | Fix unwrap L129 |
| `crates/result-import/src/diann.rs` | Modify | Fix Regex::new().unwrap() L51 |
| `crates/report/src/xic_visualize.rs` | Modify | HTML-escape peptide in title |
| `crates/entrapment-analysis/src/report.rs` | Modify | Strengthen JSON escaping |
| `crates/entrapment-analysis/templates/entrapment_report.html` | Modify | Plotly 2.35.0→2.35.2, add RT column |
| `crates/entrapment-analysis/src/mirror_plot.rs` | Modify | Plotly 2.35.0→2.35.2 |
| `crates/entrapment-analysis/src/multi_report.rs` | Modify | Plotly 2.35.0→2.35.2, add viewport meta |
| `crates/entrapment-analysis/src/types.rs` | Modify | Add doc comment linking UnifiedPsm→Psm |

---

### Task 1: Fix critical unwrap() patterns

**Files:**
- Modify: `crates/mcp-server/src/tools.rs:1765`
- Modify: `crates/result-import/src/unimod.rs:129`
- Modify: `crates/result-import/src/diann.rs:51`

**Context:** Three `unwrap()` calls that can panic at runtime. While two are "guarded" by prior checks, they are fragile — future refactoring could remove the guards. The regex unwrap is an unnecessary risk since regex compilation can fail.

- [ ] **Step 1: Fix tools.rs:1765 — replace unwrap() with safe pattern**

In `crates/mcp-server/src/tools.rs`, around line 1765, replace:

```rust
let first_path = Path::new(input.input_files.first().unwrap());
```

with:

```rust
let first_path = match input.input_files.first() {
    Some(p) => Path::new(p),
    None => {
        return Err(mcp_err(
            ErrorCode::INVALID_PARAMS,
            "input_files is empty (internal error)",
        ));
    }
};
```

- [ ] **Step 2: Fix unimod.rs:129 — replace unwrap() with if-let**

In `crates/result-import/src/unimod.rs`, around line 128-129, replace:

```rust
if key == "site" && val.len() == 1 {
    let ch = val.chars().next().unwrap();
```

with:

```rust
if key == "site" {
    if let Some(ch) = val.chars().next() {
```

Keep the existing body (the `is_ascii_uppercase` check) intact. The `val.len() == 1` check becomes redundant since `chars().next()` handles empty strings safely, but we keep the single-char semantics via the existing check inside.

- [ ] **Step 3: Fix diann.rs:51 — use std::sync::LazyLock for regex**

In `crates/result-import/src/diann.rs`, around line 51, replace:

```rust
let re = Regex::new(r"\(UniMod:(\d+)\)").unwrap();
```

with a `LazyLock` static at module level:

```rust
use std::sync::LazyLock;

static UNIMOD_RE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"\(UniMod:(\d+)\)").expect("UNIMOD_RE is a valid regex literal")
});
```

Then in the function body, use `&*UNIMOD_RE` or `UNIMOD_RE.find(...)` instead of `re`.

- [ ] **Step 4: Build and run tests**

```bash
cargo build -p protein-copilot-mcp-server -p protein-copilot-result-import --quiet 2>&1 | head -20
cargo test -p protein-copilot-result-import --quiet 2>&1 | tail -5
```

Expected: Build succeeds, all result-import tests pass.

- [ ] **Step 5: Commit**

```bash
git add crates/mcp-server/src/tools.rs crates/result-import/src/unimod.rs crates/result-import/src/diann.rs
git commit -m "fix: replace unwrap() with safe patterns in tools.rs, unimod.rs, diann.rs

- tools.rs:1765: unwrap() on input_files.first() → match with error return
- unimod.rs:129: unwrap() on chars().next() → if-let guard
- diann.rs:51: Regex::new().unwrap() → LazyLock static

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: Fix XIC title HTML escaping

**Files:**
- Modify: `crates/report/src/xic_visualize.rs:36`

**Context:** `xic_visualize.rs` line 36 injects `peptide_sequence` directly into HTML `<title>` tag via `__PEPTIDE_PLACEHOLDER__` replacement. Peptide sequences can contain modification notation with `<` and `>` characters (e.g., `PEPTM<ox>IDE`), which could break the title tag. The report crate already has `escape_json_for_html()` but we need an HTML entity escape for this context.

- [ ] **Step 1: Add simple HTML-escape before title injection**

In `crates/report/src/xic_visualize.rs`, add a helper function (or inline the escape):

```rust
fn html_escape_title(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}
```

Then change line 36 from:

```rust
.replace("__PEPTIDE_PLACEHOLDER__", &xic_data.peptide_sequence);
```

to:

```rust
.replace("__PEPTIDE_PLACEHOLDER__", &html_escape_title(&xic_data.peptide_sequence));
```

- [ ] **Step 2: Build and run tests**

```bash
cargo build -p protein-copilot-report --quiet 2>&1 | head -20
cargo test -p protein-copilot-report --quiet 2>&1 | tail -5
```

Expected: Build succeeds, all report tests pass.

- [ ] **Step 3: Commit**

```bash
git add crates/report/src/xic_visualize.rs
git commit -m "fix(report): HTML-escape peptide sequence in XIC title tag

Prevents potential HTML injection when peptide contains modification
notation with < > characters (e.g., PEPTM<ox>IDE).

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: Strengthen entrapment JSON escaping

**Files:**
- Modify: `crates/entrapment-analysis/src/report.rs:159`

**Context:** The entrapment report embeds JSON data in a `<script>` tag. Currently only `</` is escaped to `<\/`, but this doesn't protect against `<img onerror=...>` or other `<` / `>` patterns in JSON values. The report crate already has the correct pattern: replace `<` with `\u003c` and `>` with `\u003e`.

- [ ] **Step 1: Replace weak escaping with full pattern**

In `crates/entrapment-analysis/src/report.rs`, around line 159, replace:

```rust
let safe_json = json.replace("</", "<\\/");
```

with:

```rust
let safe_json = json.replace('<', r"\u003c").replace('>', r"\u003e");
```

- [ ] **Step 2: Build and run tests**

```bash
cargo build -p protein-copilot-entrapment-analysis --quiet 2>&1 | head -20
cargo test -p protein-copilot-entrapment-analysis --quiet 2>&1 | tail -5
```

Expected: Build succeeds, all 180 entrapment tests pass.

- [ ] **Step 3: Commit**

```bash
git add crates/entrapment-analysis/src/report.rs
git commit -m "fix(entrapment): strengthen JSON escaping to match report crate pattern

Replace '</'-only escaping with full < > → \\u003c \\u003e pattern.
Prevents HTML injection via <img onerror=...> or similar vectors
in JSON string values embedded in <script> tags.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 4: Standardize Plotly.js version to 2.35.2

**Files:**
- Modify: `crates/entrapment-analysis/templates/entrapment_report.html:7`
- Modify: `crates/entrapment-analysis/src/mirror_plot.rs:127`
- Modify: `crates/entrapment-analysis/src/multi_report.rs:152`

**Context:** The report crate uses Plotly 2.35.2, but entrapment-analysis uses 2.35.0 in three places. Inconsistent versions can cause subtle rendering differences.

- [ ] **Step 1: Update entrapment_report.html**

In `crates/entrapment-analysis/templates/entrapment_report.html`, line 7, replace:

```html
<script src="https://cdn.plot.ly/plotly-2.35.0.min.js"></script>
```

with:

```html
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
```

- [ ] **Step 2: Update mirror_plot.rs**

In `crates/entrapment-analysis/src/mirror_plot.rs`, line 127, replace:

```html
<script src="https://cdn.plot.ly/plotly-2.35.0.min.js"></script>
```

with:

```html
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
```

- [ ] **Step 3: Update multi_report.rs**

In `crates/entrapment-analysis/src/multi_report.rs`, line 152, replace:

```html
<script src="https://cdn.plot.ly/plotly-2.35.0.min.js"></script>
```

with:

```html
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
```

- [ ] **Step 4: Build and run tests**

```bash
cargo build -p protein-copilot-entrapment-analysis --quiet 2>&1 | head -20
cargo test -p protein-copilot-entrapment-analysis --quiet 2>&1 | tail -5
```

Expected: Build succeeds, all 180 tests pass.

- [ ] **Step 5: Commit**

```bash
git add crates/entrapment-analysis/templates/entrapment_report.html \
        crates/entrapment-analysis/src/mirror_plot.rs \
        crates/entrapment-analysis/src/multi_report.rs
git commit -m "chore(entrapment): standardize Plotly.js version to 2.35.2

Align with report crate; was 2.35.0 in 3 entrapment files.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 5: Add viewport meta to multi_report.rs

**Files:**
- Modify: `crates/entrapment-analysis/src/multi_report.rs:150`

**Context:** The per-PSM provenance HTML generated by `multi_report.rs` is missing the viewport meta tag. All other templates (entrapment_report.html, mirror_plot.rs, annotation.html, xic.html, unified.html) have it.

- [ ] **Step 1: Add viewport meta after charset meta**

In `crates/entrapment-analysis/src/multi_report.rs`, around line 150, change:

```rust
<meta charset="utf-8">
<title>Provenance: {trap} (scan {scan})</title>
```

to:

```rust
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Provenance: {trap} (scan {scan})</title>
```

- [ ] **Step 2: Build and verify**

```bash
cargo build -p protein-copilot-entrapment-analysis --quiet 2>&1 | head -5
```

Expected: Build succeeds.

- [ ] **Step 3: Commit**

```bash
git add crates/entrapment-analysis/src/multi_report.rs
git commit -m "fix(entrapment): add viewport meta to provenance report HTML

Aligns with all other HTML templates for consistent mobile rendering.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 6: Extract magic numbers to named constants

**Files:**
- Modify: `crates/mcp-server/src/tools.rs`

**Context:** Three hardcoded magic numbers in tools.rs hurt readability: DIA detection threshold (5.0 Da), RT auto-lookup tolerance (0.5 min), and FDR 1% filter (0.01). Extract them to named constants near the top of the file alongside existing constants like `MAX_CACHE_SIZE`.

- [ ] **Step 1: Add named constants**

Near line 968 in `crates/mcp-server/src/tools.rs` (where `MAX_CACHE_SIZE` is defined), add:

```rust
/// DIA isolation window detection threshold (Da).
/// Spectra with median isolation window wider than this are classified as DIA.
const DIA_ISOLATION_WINDOW_THRESHOLD_DA: f64 = 5.0;

/// Default RT tolerance (minutes) for auto-scanning MS2 lookup.
const RT_AUTO_LOOKUP_TOLERANCE_MIN: f64 = 0.5;

/// Default FDR threshold (1%) for protein inference filtering.
const FDR_1PCT_THRESHOLD: f64 = 0.01;
```

- [ ] **Step 2: Replace magic numbers with constants**

Line ~1770: replace `if w > 5.0` with `if w > DIA_ISOLATION_WINDOW_THRESHOLD_DA`

Line ~2313: replace `find_by_rt(&spectrum_file, rt, precursor_mz, 0.5)` with `find_by_rt(&spectrum_file, rt, precursor_mz, RT_AUTO_LOOKUP_TOLERANCE_MIN)`

Line ~3302: replace `.is_some_and(|q| q <= 0.01)` with `.is_some_and(|q| q <= FDR_1PCT_THRESHOLD)`

- [ ] **Step 3: Build**

```bash
cargo build -p protein-copilot-mcp-server --quiet 2>&1 | head -5
```

Expected: Build succeeds.

- [ ] **Step 4: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "refactor(mcp): extract magic numbers to named constants

- DIA_ISOLATION_WINDOW_THRESHOLD_DA = 5.0
- RT_AUTO_LOOKUP_TOLERANCE_MIN = 0.5
- FDR_1PCT_THRESHOLD = 0.01

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 7: Add RT column to entrapment table and add UnifiedPsm doc comment

**Files:**
- Modify: `crates/entrapment-analysis/templates/entrapment_report.html`
- Modify: `crates/entrapment-analysis/src/types.rs`

**Context:** `retention_time` is serialized in `PsmRow` but not rendered in the entrapment table. The JS `renderTable()` maps 21 fields by name from JSON but skips `retention_time`. Add the column for completeness. Also add a doc comment to `UnifiedPsm` documenting its relationship to `core::Psm`.

- [ ] **Step 1: Add RT column header to entrapment_report.html**

In the `<thead>` section (around line 79, after the Scan column), add a new header:

```html
<th onclick="sortTable(5)">Scan ⇅</th>
```

becomes:

```html
<th onclick="sortTable(5)">Scan ⇅</th>
<th onclick="sortTable(6)">RT (min) ⇅</th>
```

Then shift all subsequent `sortTable(N)` indices by +1 (Group becomes 7, Level becomes 8, ... through Chimeric becomes 21).

- [ ] **Step 2: Add RT to renderTable() JS**

In the `renderTable()` function, add `p.retention_time` after `p.scan_number`:

```javascript
p.scan_number || '',
p.retention_time || '',    // new
p.group || '',
```

- [ ] **Step 3: Add RT cell to populateTable() JS**

In the `populateTable()` function, add a cell after r[5] (scan):

```javascript
'<td>' + escHtml(r[5]) + '</td>' +
'<td>' + escHtml(r[6]) + '</td>' +  // RT (min) — new
```

Update all subsequent indices: r[6]→r[7] for group, r[7]→r[8] for level, etc., through r[20]→r[21] for chimeric.

- [ ] **Step 4: Update filterTable() JS indices**

In `filterTable()`, update the level index from `r[7]` to `r[8]`, and the haystack search index from `r[4]` to `r[4]` (file stays at index 4, no change needed).

- [ ] **Step 5: Add UnifiedPsm doc comment**

In `crates/entrapment-analysis/src/types.rs`, add a doc comment above `pub struct UnifiedPsm`:

```rust
/// A PSM representation that extends [`protein_copilot_core::search_result::Psm`]
/// with entrapment-specific fields (group, origin file).
///
/// Used as the unified input format across all entrapment analysis stages.
```

- [ ] **Step 6: Build and run tests**

```bash
cargo build -p protein-copilot-entrapment-analysis --quiet 2>&1 | head -20
cargo test -p protein-copilot-entrapment-analysis --quiet 2>&1 | tail -5
```

Expected: Build succeeds, all 180 tests pass.

- [ ] **Step 7: Commit**

```bash
git add crates/entrapment-analysis/templates/entrapment_report.html \
        crates/entrapment-analysis/src/types.rs
git commit -m "feat(entrapment): add RT column to PSM table, document UnifiedPsm

- Add 'RT (min)' column to entrapment report table (22 columns total)
- Update all JS column indices accordingly
- Add doc comment linking UnifiedPsm to core::Psm

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Summary

| Task | Priority | Files | Description |
|------|----------|-------|-------------|
| T1 | 🔴 HIGH | tools.rs, unimod.rs, diann.rs | Fix 3 unwrap() panic risks |
| T2 | 🔴 HIGH | xic_visualize.rs | HTML-escape peptide in XIC title |
| T3 | 🔴 HIGH | entrapment report.rs | Strengthen JSON escaping |
| T4 | 🟡 MEDIUM | 3 entrapment files | Plotly 2.35.0→2.35.2 |
| T5 | 🟡 MEDIUM | multi_report.rs | Add viewport meta |
| T6 | 🟡 MEDIUM | tools.rs | Extract magic numbers |
| T7 | 🟢 LOW | entrapment_report.html, types.rs | Add RT column + UnifiedPsm doc |

**Dependencies:** None — all tasks are independent and can be executed in any order.

**Not included (out of scope):**
- Test coverage gaps (dia-extraction, mcp-server tests) — these are new feature work, not audit fixes
- `detect_format()` rename — would break API, needs separate plan
- entrapment submodule visibility (`pub(crate)`) — optional refactor, low value
