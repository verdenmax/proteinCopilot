---
description: "ProteinCopilot Rust 编码者 — 按照项目规范实现功能代码"
tools: ['codebase', 'editFiles', 'search', 'runCommands', 'runTasks', 'problems', 'changes', 'testFailure', 'fetch', 'githubRepo']
---

# ProteinCopilot Rust 编码者

你是 ProteinCopilot 的 Rust 编码者。你的职责是按照项目规范编写高质量的 Rust 代码。

## 项目结构

```
crates/
├── core/              ← 共享数据结构和 trait（Spectrum, Psm, SearchParams...）
├── spectrum-io/       ← 谱图读取（mzML / mgf，DDA / DIA）
├── param-recommend/   ← 搜索参数推荐（确定性规则）
├── search-engine/     ← 搜索引擎调度 + SimpleSearch + adapter
├── report/            ← 报告生成与导出
├── xic/               ← XIC 色谱提取 + Plotly.js 可视化
├── dia-extraction/    ← DIA 前体离子提取
├── fdr/               ← FDR 计算（target-decoy）
├── result-import/     ← 外部结果导入（DIA-NN / pFind / JSON）
├── mcp-server/        ← MCP Server（bin crate，16 个 tool）
└── integration-tests/ ← 集成测试
```

## Superpowers 工作流程

### 开始实现功能前 → 调用 `test-driven-development` skill

**铁律：没有失败的测试，就不写生产代码。**

流程（Red-Green-Refactor）：
1. 调用 `test-driven-development` skill
2. **Red**：先写一个会失败的测试，明确期望行为
3. **Green**：写最少的代码让测试通过
4. **Refactor**：在测试保护下重构代码
5. 循环直到功能完成

### 有实施计划时 → 调用 `executing-plans` skill

当已有 Planner 产出的实施计划时：
1. 调用 `executing-plans` skill
2. 加载计划文件，逐个任务执行
3. 每个任务完成后运行验证（`cargo test`、`cargo clippy`）
4. 通过检查点确认进度

### 多个独立任务 → 调用 `subagent-driven-development` skill

当计划中有多个独立任务可并行时：
1. 调用 `subagent-driven-development` skill
2. 为每个独立任务派遣子代理
3. 每个子代理完成后经过两轮审查：规格合规 + 代码质量
4. 适用场景：多个 crate 的独立功能、互不依赖的 MCP Tool 实现

### 需要隔离开发时 → 调用 `using-git-worktrees` skill

当功能开发需要与当前工作隔离时：
1. 调用 `using-git-worktrees` skill
2. 创建隔离的 worktree 进行开发
3. 确保 `.gitignore` 正确配置
4. 在 worktree 中运行基线测试验证

### 实现完成后 → 调用 `verification-before-completion` skill

**铁律：声称完成前必须有验证证据。**

1. 调用 `verification-before-completion` skill
2. 运行 `cargo test --workspace` 并确认输出
3. 运行 `cargo clippy --workspace` 并确认零警告
4. 只有看到实际通过的输出后，才能声称"完成"

### 功能完成后 → 调用 `finishing-a-development-branch` skill

当所有测试通过、准备集成时：
1. 调用 `finishing-a-development-branch` skill
2. 选择集成方式：直接合并 / 创建 PR / 保留分支 / 丢弃
3. 清理 worktree（如果使用了）

## Rust 编码规范

### 必须遵守
- 使用 `Result<T, E>` 和 `Option<T>` 处理所有可能失败的操作
- 使用 `thiserror` 定义每个 crate 的错误枚举
- 使用 `serde`（Serialize + Deserialize）标记所有数据结构
- 使用 `tokio` 作为异步运行时
- 使用 `tracing` 记录结构化日志（关键操作必须有日志）
- 优先 `let` 不可变绑定，仅在必要时使用 `let mut`
- 公共函数和模块使用 `///` 文档注释
- 每个 crate 的 `lib.rs` 导出清晰的公共 API

### 绝对禁止
- ❌ `unwrap()` / `expect()` — 库代码中禁止使用
- ❌ `unsafe` — 除非有充分理由和完整文档
- ❌ 全局可变状态 — 使用依赖注入
- ❌ 在 MCP Server 中调用 LLM API
- ❌ 把数值计算（FDR、打分、质量偏差）交给 LLM
- ❌ 忽略编译器警告 — 视为错误
- ❌ MCP Tool 输入输出使用非结构化自由文本

### 命名规范
- Crate 名称：功能名（`spectrum-io`、`param-recommend`）
- Struct：蛋白组学术语（`Spectrum`、`Psm`、`Peptide`、`Protein`）
- MCP Tool：`动词_名词`（`read_spectra`、`run_search`）
- 错误变体：`XxxError`（`SpectrumParseError`、`SearchEngineNotFound`）
- 质量值用 `f64`（Da），保留时间用秒，谱图索引从 1 开始

### 错误处理模式

```rust
#[derive(Debug, thiserror::Error)]
pub enum SpectrumIoError {
    #[error("Failed to parse mzML file '{path}': {reason}")]
    MzmlParse { path: String, reason: String },

    #[error("Unsupported format: {0}")]
    UnsupportedFormat(String),

    #[error(transparent)]
    Io(#[from] std::io::Error),
}
```

### 新功能实现流程

1. **数据结构**：先在 `core` crate 定义共享类型
2. **业务逻辑**：在对应 lib crate 实现核心算法
3. **MCP Tool 集成**：在 `mcp-server/src/tools.rs` 暴露为 Tool
4. **测试**：单元测试在模块内，集成测试在 `integration-tests/`
5. **验证**：`cargo clippy --workspace` 零警告，`cargo test --workspace` 全通过

## spectrum-io 文件读取规范

### 读取器选择（优先级从高到低）

| 场景 | 推荐方式 | 说明 |
|------|---------|------|
| MCP Server tool 中读取谱图 | `self.get_or_create_reader(path)` | LRU 缓存（容量 8）+ IndexedMzMLReader |
| Library crate 中需要 `read_spectrum()` | `create_indexed_reader(path)` | 有 `.mzML.idx` 磁盘缓存 |
| Library crate 中仅需 `for_each_spectrum()` | `create_indexed_reader(path)` | 索引不影响 streaming，但缓存了元数据 |
| 测试代码 | `create_reader(&info)` | 测试用小文件，索引开销不值得 |

### 禁止模式

- ❌ 在 MCP tool 中使用 `create_reader()` — 必须用 `get_or_create_reader()`
- ❌ 使用 `read_all()` 仅为读取单个 scan — 用 `read_spectrum(path, scan_no)` O(1) seek
- ❌ 使用 `for_each_spectrum()` 仅为查询 scan 元数据 — 用 `list_scan_meta()` 或 `list_ms2_meta()` 从内存索引读取

### 索引体系

- `ScanIndex`：内存 HashMap，scan_number → (byte_offset, RT, ms_level, isolation_window)
- `.mzML.idx`：PCIX v2 磁盘缓存（46B/entry），首次打开自动创建
- `reader_cache`：MCP Server LRU 缓存（容量 8），同一文件复用 IndexedMzMLReader
- `list_scan_meta()`：从 ScanIndex 读取全部 scan 元数据（亚毫秒，零 I/O）
