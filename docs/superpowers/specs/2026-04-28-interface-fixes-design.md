# 接口修复设计：Serde 缺陷 + API 可见性 + 命名冲突 + TSV 常量

**来源**: 2026-04-27 全项目审计（4 个并行 agent 审计）

**目标**: 修复 1 个序列化 bug、2 个 API 可见性问题、1 个命名冲突、1 个硬编码脆弱性。

**范围**: 仅修改已有代码，不新增功能。所有改动都是内部重构，不影响 MCP 协议接口。

---

## F1 🔴 `SearchProgress.error_category` 缺少 `#[serde(default)]`

**问题**: `crates/core/src/progress.rs:27` — `error_category` 字段有 `skip_serializing_if`
但没有 `#[serde(default)]`。序列化时 `None` 值被跳过，反序列化时因字段缺失而失败。

**修复**:
```rust
// 修复前
#[serde(skip_serializing_if = "Option::is_none")]
pub error_category: Option<crate::diagnostics::ErrorCategory>,

// 修复后
#[serde(default, skip_serializing_if = "Option::is_none")]
pub error_category: Option<crate::diagnostics::ErrorCategory>,
```

**验证**: 添加 roundtrip 单元测试（序列化后反序列化，确认 `error_category: None` 正确处理）。

---

## F2 🟡 废弃 XIC 函数降为 `pub(crate)`

**问题**: `extract_xic()` (L259) 和 `extract_xic_with_raw()` (L670) 已标记 `#[deprecated]`
但仍为 `pub`，外部 crate 仍可调用。

**修复**: 将两个函数从 `pub` 改为 `pub(crate)`。

**影响**: 仅影响 crate 内部可见性，无外部调用者。

---

## F3 🟡 XIC 内部辅助函数降为 `pub(crate)`

**问题**: 4 个函数为实现细节但暴露为 `pub`：
- `extract_intensity()` (L54) — 仅 extract.rs 内部使用
- `same_isolation_window()` (L127) — 仅 extract.rs 内部使用
- `build_target_ions()` (L148) — extract.rs 内部 + 1 处集成测试
- `compute_ion_metadata()` (L198) — 仅 extract.rs 内部使用

**修复**:
- `extract_intensity`, `same_isolation_window`, `compute_ion_metadata` → `pub(crate)`
- `build_target_ions` → 保持 `pub`（被 `crates/integration-tests/tests/xic_scenarios.rs:6` 引用）

**决策理由**: `build_target_ions` 是构建目标离子列表的核心函数，集成测试使用它来验证 XIC
提取的正确性。保持 pub 是合理的 — 它是一个有意义的 API，不像其他 3 个只是实现细节。

---

## F4 🟡 Entrapment `RunMetadata` 重命名为 `EntrapmentRunMetadata`

**问题**: `crates/entrapment-analysis/src/output.rs:25` 定义了 `RunMetadata`，
与 `crates/core/src/run_metadata.rs:92` 的 `RunMetadata` 同名。
虽然目前不交叉引用，但在 `tools.rs:3593` 中 `use output::RunMetadata` 与
`use protein_copilot_core::RunMetadata` 容易混淆。

**修复**: 重命名为 `EntrapmentRunMetadata`。

**影响范围**（6 处引用）:
1. `entrapment-analysis/src/output.rs:25` — struct 定义
2. `entrapment-analysis/src/output.rs:249` — `write_run_metadata` 参数类型
3. `entrapment-analysis/src/output.rs:515` — struct 构造
4. `entrapment-analysis/tests/v3_e2e_provenance.rs:8` — use 导入
5. `entrapment-analysis/tests/v3_e2e_provenance.rs:523` — struct 构造
6. `mcp-server/src/tools.rs:3593,3653` — use 导入 + struct 构造

**JSON 输出不变**: TSV/JSON 文件格式不受影响（serde 序列化键名不变）。

---

## F5 🟡 TSV 列名抽取为共享常量

**问题**: `mcp-server/src/tools.rs:3736-3738` 硬编码字符串 `"level"`, `"delta_mass_da"`,
`"best_target_protein"` 来查找 TSV 列。`output.rs:108-132` 也硬编码相同字符串写 header。
如果列名修改，两处必须同步变更，否则静默丢数据。

**修复**: 在 `entrapment-analysis/src/output.rs` 中定义列名常量模块，
`write_classified_tsv` 和 MCP tool 都引用同一常量：

```rust
/// Column names for the classified entrapment TSV output.
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

MCP tool 中改为:
```rust
let level_idx = headers.iter().position(|h| h == output::columns::LEVEL);
let delta_idx = headers.iter().position(|h| h == output::columns::DELTA_MASS_DA);
let target_protein_idx = headers.iter().position(|h| h == output::columns::BEST_TARGET_PROTEIN);
```

---

## 改动汇总

| 文件 | 改动类型 | 行数估计 |
|------|----------|----------|
| `core/src/progress.rs` | 加 `#[serde(default)]` + roundtrip 测试 | ~15 行 |
| `xic/src/extract.rs` | 5 处 `pub` → `pub(crate)` | 5 行 |
| `entrapment-analysis/src/output.rs` | 重命名 + 常量模块 | ~35 行 |
| `entrapment-analysis/tests/v3_e2e_provenance.rs` | 跟随重命名 | 2 行 |
| `mcp-server/src/tools.rs` | 跟随重命名 + 用常量 | 5 行 |

**总改动**: ~60 行，5 个文件
