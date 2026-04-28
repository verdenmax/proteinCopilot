# 接口修复 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 修复审计发现的 5 个接口问题：1 个 serde roundtrip bug、2 个 API 可见性泄漏、1 个命名冲突、1 个硬编码脆弱性。

**Architecture:** 纯重构，不新增功能。改动集中在 5 个文件，每个 Task 独立可编译、可测试、可提交。

**Tech Stack:** Rust, serde, schemars

**Spec:** `docs/superpowers/specs/2026-04-28-interface-fixes-design.md`

---

## File Map

| 文件 | 改动 |
|------|------|
| `crates/core/src/progress.rs` | F1: 加 `#[serde(default)]` + roundtrip 测试 |
| `crates/xic/src/extract.rs` | F2+F3: 5 处 `pub` → `pub(crate)` |
| `crates/entrapment-analysis/src/output.rs` | F4: 重命名 `RunMetadata` → `EntrapmentRunMetadata`；F5: 新增列名常量模块 |
| `crates/entrapment-analysis/tests/v3_e2e_provenance.rs` | F4: 跟随重命名 |
| `crates/mcp-server/src/tools.rs` | F4: 跟随重命名；F5: 使用列名常量 |

---

### Task 1: F1 — 修复 `SearchProgress.error_category` serde roundtrip bug

**Files:**
- Modify: `crates/core/src/progress.rs:27`

- [ ] **Step 1: 写 roundtrip 失败测试**

在 `crates/core/src/progress.rs` 文件末尾追加测试模块：

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn search_progress_roundtrip_without_error_category() {
        let progress = SearchProgress {
            run_id: Uuid::new_v4(),
            status: "Running".to_string(),
            stage: Some("Matching spectra".to_string()),
            progress_pct: Some(0.5),
            elapsed_sec: 12.3,
            estimated_remaining_sec: None,
            error_category: None,
            has_diagnostics: false,
        };

        let json = serde_json::to_string(&progress).unwrap();
        // error_category: None should be skipped in JSON
        assert!(!json.contains("error_category"), "None field should be skipped");

        // Roundtrip: deserialize back should NOT fail
        let restored: SearchProgress = serde_json::from_str(&json).unwrap();
        assert_eq!(restored.status, "Running");
        assert!(restored.error_category.is_none());
    }
}
```

- [ ] **Step 2: 运行测试确认失败**

Run: `cargo test -p protein-copilot-core search_progress_roundtrip -- --nocapture 2>&1 | tail -10`

Expected: FAIL — `missing field "error_category"` 反序列化错误

- [ ] **Step 3: 修复 — 加 `#[serde(default)]`**

在 `crates/core/src/progress.rs` 中将第 27 行：

```rust
    #[serde(skip_serializing_if = "Option::is_none")]
```

替换为：

```rust
    #[serde(default, skip_serializing_if = "Option::is_none")]
```

- [ ] **Step 4: 运行测试确认通过**

Run: `cargo test -p protein-copilot-core search_progress_roundtrip -- --nocapture 2>&1 | tail -5`

Expected: `test progress::tests::search_progress_roundtrip_without_error_category ... ok`

- [ ] **Step 5: Commit**

```bash
git add crates/core/src/progress.rs
git commit -m "fix(core): add #[serde(default)] to SearchProgress.error_category

Without this, deserialization fails when error_category is omitted
from JSON (which happens because skip_serializing_if skips None).

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: F2+F3 — XIC 函数可见性降级

**Files:**
- Modify: `crates/xic/src/extract.rs:54,127,259,672`

- [ ] **Step 1: 将 3 个内部辅助函数从 `pub` 改为 `pub(crate)`**

在 `crates/xic/src/extract.rs` 中，替换以下 3 处：

第 54 行：`pub fn extract_intensity(` → `pub(crate) fn extract_intensity(`

第 127 行：`pub fn same_isolation_window(` → `pub(crate) fn same_isolation_window(`

第 198 行：`pub fn compute_ion_metadata(` → `pub(crate) fn compute_ion_metadata(`

注意：`build_target_ions`（第 148 行）保持 `pub`，因为被 `crates/integration-tests/tests/xic_scenarios.rs:6` 引用。

- [ ] **Step 2: 将 2 个废弃函数从 `pub` 改为 `pub(crate)`**

第 261 行：`pub fn extract_xic(` → `pub(crate) fn extract_xic(`

第 673 行：`pub fn extract_xic_with_raw(` → `pub(crate) fn extract_xic_with_raw(`

- [ ] **Step 3: 编译验证**

Run: `cargo build --workspace 2>&1 | tail -5`

Expected: 编译成功，无错误。可能出现新的 dead_code warning（因为 pub(crate) 函数在 crate 内部未被非 deprecated 代码使用），这是预期的。

- [ ] **Step 4: 运行 XIC 测试**

Run: `cargo test -p protein-copilot-xic 2>&1 | tail -5`

Expected: 所有测试通过

- [ ] **Step 5: 运行集成测试**

Run: `cargo test -p integration-tests 2>&1 | tail -10`

Expected: 所有测试通过（`build_target_ions` 仍为 pub，`xic_scenarios.rs` 不受影响）

- [ ] **Step 6: Commit**

```bash
git add crates/xic/src/extract.rs
git commit -m "refactor(xic): restrict internal helpers to pub(crate)

- extract_intensity, same_isolation_window, compute_ion_metadata → pub(crate)
- extract_xic, extract_xic_with_raw (deprecated) → pub(crate)
- build_target_ions stays pub (used by integration tests)

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: F4 — 重命名 Entrapment `RunMetadata` → `EntrapmentRunMetadata`

**Files:**
- Modify: `crates/entrapment-analysis/src/output.rs:25,249,515`
- Modify: `crates/entrapment-analysis/tests/v3_e2e_provenance.rs:8,523`
- Modify: `crates/mcp-server/src/tools.rs:3593,3653`

- [ ] **Step 1: 重命名 struct 定义**

在 `crates/entrapment-analysis/src/output.rs` 第 25 行：

`pub struct RunMetadata {` → `pub struct EntrapmentRunMetadata {`

- [ ] **Step 2: 重命名 `write_run_metadata` 参数类型**

在 `crates/entrapment-analysis/src/output.rs` 第 249 行：

`pub fn write_run_metadata(metadata: &RunMetadata, path: &Path)` → `pub fn write_run_metadata(metadata: &EntrapmentRunMetadata, path: &Path)`

- [ ] **Step 3: 重命名 output.rs 测试中的构造**

在 `crates/entrapment-analysis/src/output.rs` 第 515 行：

`let metadata = RunMetadata {` → `let metadata = EntrapmentRunMetadata {`

- [ ] **Step 4: 更新集成测试 import 和构造**

在 `crates/entrapment-analysis/tests/v3_e2e_provenance.rs` 第 8 行：

`use protein_copilot_entrapment_analysis::output::{write_classified_tsv, write_run_metadata, RunMetadata};`

→

`use protein_copilot_entrapment_analysis::output::{write_classified_tsv, write_run_metadata, EntrapmentRunMetadata};`

第 523 行：

`let metadata = RunMetadata {` → `let metadata = EntrapmentRunMetadata {`

- [ ] **Step 5: 更新 MCP server import 和构造**

在 `crates/mcp-server/src/tools.rs` 第 3593 行：

`output::{self, RunMetadata},` → `output::{self, EntrapmentRunMetadata},`

第 3653 行：

`let metadata = RunMetadata {` → `let metadata = EntrapmentRunMetadata {`

- [ ] **Step 6: 编译验证**

Run: `cargo build --workspace 2>&1 | tail -5`

Expected: 编译成功

- [ ] **Step 7: 运行 entrapment 测试**

Run: `cargo test -p protein-copilot-entrapment-analysis 2>&1 | tail -5`

Expected: 所有测试通过

- [ ] **Step 8: Commit**

```bash
git add crates/entrapment-analysis/src/output.rs \
       crates/entrapment-analysis/tests/v3_e2e_provenance.rs \
       crates/mcp-server/src/tools.rs
git commit -m "refactor(entrapment): rename RunMetadata → EntrapmentRunMetadata

Eliminates naming collision with core::RunMetadata.
No change to serialized JSON format.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 4: F5 — TSV 列名抽取为共享常量

**Files:**
- Modify: `crates/entrapment-analysis/src/output.rs:108-132`
- Modify: `crates/mcp-server/src/tools.rs:3736-3738`

- [ ] **Step 1: 在 output.rs 中添加列名常量模块**

在 `crates/entrapment-analysis/src/output.rs` 的 `write_classified_tsv` 函数前（第 99 行之前），插入：

```rust
/// Column names for the classified entrapment TSV output.
///
/// Used by both `write_classified_tsv` (writer) and MCP tool
/// `analyze_entrapment_stats` (reader) to keep column names in sync.
pub mod columns {
    pub const PEPTIDE: &str = "peptide";
    pub const CHARGE: &str = "charge";
    pub const PRECURSOR_MZ: &str = "precursor_mz";
    pub const RETENTION_TIME: &str = "retention_time";
    pub const SCAN_NUMBER: &str = "scan_number";
    pub const SPECTRUM_FILE: &str = "spectrum_file";
    pub const PROTEIN_IDS: &str = "protein_ids";
    pub const Q_VALUE: &str = "q_value";
    pub const GROUP: &str = "group";
    pub const LEVEL: &str = "level";
    pub const BEST_TARGET_PEPTIDE: &str = "best_target_peptide";
    pub const BEST_TARGET_PROTEIN: &str = "best_target_protein";
    pub const MISMATCHES: &str = "mismatches";
    pub const DELTA_MASS_DA: &str = "delta_mass_da";
    pub const DIFF_POSITIONS: &str = "diff_positions";
    pub const SUBSTITUTION_TYPE: &str = "substitution_type";
    pub const EDIT_DISTANCE: &str = "edit_distance";
    pub const ALIGNMENT_DETAIL: &str = "alignment_detail";
    pub const TRAP_MATCHED: &str = "trap_matched";
    pub const TARGET_MATCHED: &str = "target_matched";
    pub const SHARED_IONS: &str = "shared_ions";
    pub const SHARED_RATIO: &str = "shared_ratio";
    pub const IS_CHIMERIC: &str = "is_chimeric";
}
```

- [ ] **Step 2: 将 `write_classified_tsv` 的 header 改为使用常量**

将 `write_classified_tsv` 中的 `wtr.write_record([ ... ])` 调用（第 108-132 行）替换为：

```rust
    wtr.write_record([
        columns::PEPTIDE,
        columns::CHARGE,
        columns::PRECURSOR_MZ,
        columns::RETENTION_TIME,
        columns::SCAN_NUMBER,
        columns::SPECTRUM_FILE,
        columns::PROTEIN_IDS,
        columns::Q_VALUE,
        columns::GROUP,
        columns::LEVEL,
        columns::BEST_TARGET_PEPTIDE,
        columns::BEST_TARGET_PROTEIN,
        columns::MISMATCHES,
        columns::DELTA_MASS_DA,
        columns::DIFF_POSITIONS,
        columns::SUBSTITUTION_TYPE,
        columns::EDIT_DISTANCE,
        columns::ALIGNMENT_DETAIL,
        columns::TRAP_MATCHED,
        columns::TARGET_MATCHED,
        columns::SHARED_IONS,
        columns::SHARED_RATIO,
        columns::IS_CHIMERIC,
    ])
```

- [ ] **Step 3: 更新 MCP tool 中的列名查找**

在 `crates/mcp-server/src/tools.rs` 中，先在 `analyze_entrapment_stats` 函数内部
（第 3590 行区域）添加 columns import：

将现有的：
```rust
        use protein_copilot_entrapment_analysis::{
            config::EntrapmentConfig,
            loader::{self, ResultFormat},
            output::{self, EntrapmentRunMetadata},
            EntrapmentAnalyzer,
        };
```

注意：这个 use block 在 `classify_entrapment_hits` 函数中，不在 `analyze_entrapment_stats` 中。
需要在 `analyze_entrapment_stats` 函数体内（约第 3710 行后）添加：

```rust
        use protein_copilot_entrapment_analysis::output::columns;
```

然后将第 3736-3738 行的硬编码字符串：

```rust
        let level_idx = headers.iter().position(|h| h == "level");
        let delta_idx = headers.iter().position(|h| h == "delta_mass_da");
        let target_protein_idx = headers.iter().position(|h| h == "best_target_protein");
```

替换为：

```rust
        let level_idx = headers.iter().position(|h| h == columns::LEVEL);
        let delta_idx = headers.iter().position(|h| h == columns::DELTA_MASS_DA);
        let target_protein_idx = headers.iter().position(|h| h == columns::BEST_TARGET_PROTEIN);
```

- [ ] **Step 4: 编译验证**

Run: `cargo build --workspace 2>&1 | tail -5`

Expected: 编译成功

- [ ] **Step 5: 运行 entrapment 测试**

Run: `cargo test -p protein-copilot-entrapment-analysis 2>&1 | tail -5`

Expected: 所有测试通过（TSV header 内容不变，只是改用常量引用）

- [ ] **Step 6: Commit**

```bash
git add crates/entrapment-analysis/src/output.rs \
       crates/mcp-server/src/tools.rs
git commit -m "refactor(entrapment): extract TSV column names as shared constants

Single source of truth for column names between writer
(write_classified_tsv) and reader (analyze_entrapment_stats).

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 5: 全量验证

- [ ] **Step 1: Workspace 全量编译**

Run: `cargo build --workspace 2>&1 | tail -5`

Expected: 编译成功，无 error

- [ ] **Step 2: Workspace 全量测试**

Run: `cargo test --workspace 2>&1 | tail -20`

Expected: 所有测试通过

- [ ] **Step 3: Clippy 检查**

Run: `cargo clippy --workspace 2>&1 | tail -10`

Expected: 无 error（warning 可接受）
