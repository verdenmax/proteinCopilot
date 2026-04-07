# 开发指南

## 项目结构

```text
crates/
├── core/              共享类型（Spectrum, SearchParams, SearchResult 等）
├── spectrum-io/       谱图文件解析（mgf/mzML）
├── param-recommend/   参数推荐规则引擎
├── search-engine/     搜索引擎（SimpleSearch + pFind 预留）
├── dia-extraction/    DIA 前体离子提取（同位素模式检测 + MS1↔MS2 关联）
├── report/            报告生成（摘要 + TSV/JSON 导出）
└── mcp-server/        MCP Server 二进制（组装所有 tool）
```

**依赖方向**：library crate → core。mcp-server 依赖所有 library。Library crate 之间无循环依赖。

## 添加新搜索引擎 Adapter

1. 在 `crates/search-engine/src/adapters/` 下创建新文件（如 `msfragger.rs`）
2. 实现 `SearchEngineAdapter` trait：

```rust
#[async_trait::async_trait]
impl SearchEngineAdapter for MSFraggerAdapter {
    async fn search(&self, params: &SearchParams, input_files: &[PathBuf])
        -> Result<SearchResult, CoreError> { ... }
    fn engine_info(&self) -> EngineInfo { ... }
    async fn health_check(&self) -> Result<HealthStatus, CoreError> { ... }
}
```

3. 在 `mcp-server/src/tools.rs` 的 `ProteinCopilotServer::new()` 中注册：

```rust
registry.register(Box::new(MSFraggerAdapter::new(config)));
```

4. `run_search` tool 通过 `EngineRegistry` 自动发现新引擎。

## 添加新 MCP Tool

1. 在 `crates/mcp-server/src/tools.rs` 的 `#[rmcp::tool_router] impl` 块中添加：

```rust
#[rmcp::tool(
    name = "my_tool",
    description = "Tool description for LLM"
)]
fn my_tool(&self, Parameters(input): Parameters<MyInput>) -> Result<Json<MyOutput>, ErrorData> {
    // 调用 library crate 函数
    let result = my_library::do_something(&input.param)?;
    Ok(Json(result))
}
```

2. 定义输入结构体（derive `Deserialize` + `schemars::JsonSchema`）
3. 错误统一使用 `mcp_core_err()` 转换

## 运行测试

```bash
# 全部测试
cargo test

# 单个 crate
cargo test -p protein-copilot-core
cargo test -p protein-copilot-spectrum-io

# 端到端测试
cargo test -p protein-copilot-search-engine --test e2e_integration

# 带输出
cargo test -- --nocapture
```

## 代码质量

```bash
# Clippy (0 warnings 要求)
cargo clippy --all-targets

# 格式化
cargo fmt --check

# 构建 MCP Server
cargo build --release -p protein-copilot-mcp-server
```

## 关键设计原则

1. **确定性/LLM 分层**：Rust 做计算，LLM 做理解和解释
2. **所有 f64 必须验证**：`is_finite()` + 物理约束
3. **不允许 library 代码 unwrap()**：用 `Result` + `?`
4. **MCP Tool 只做胶水**：参数解析 → 调 library → 返回 JSON
5. **错误必须含 suggestion**：通过 `CoreError::suggestion()` 提供修复建议

## Agent 与 Skill（Prompt）

```text
.github/
├── agents/proteomics-search.agent.md    ← Agent 定义
├── prompts/basic-search.prompt.md       ← Skill: 基础搜索流程
└── prompts/result-interpretation.prompt.md  ← Skill: 结果解读
```

### 关系

- **Agent**（`.agent.md`）：定义 LLM 的角色、可用 tools、工作流程、决策边界。Agent 是"长期身份"，持续整个对话。
- **Skill/Prompt**（`.prompt.md`）：可复用的任务模板，用户通过 `/` 命令触发。Skill 是"短期任务"，执行特定流程。

### Agent 调用 Skill

用户对话中，Agent 可以参考 Skill 的步骤执行操作：
1. 用户说"帮我搜索这个数据" → Agent 按 `basic-search.prompt.md` 的流程执行
2. 搜索完成后说"解读一下结果" → Agent 按 `result-interpretation.prompt.md` 分析

### 编写规范

- Agent 必须列出所有可用 `tools`（frontmatter）
- Agent 必须定义决策边界（什么可以自动执行，什么需要用户确认）
- Skill 必须说明输入要求、预期输出、适用场景
- Agent 调用 MCP Tool 前不能凭空推断数据（§4.3）

### 搜索引擎说明

当前 MVP 使用 **SimpleSearchEngine**（内置简化搜索引擎），不需要 SSH 或外部依赖。
后续接入 pFind 后，Agent 工作流程不变，`run_search` tool 会通过 `EngineRegistry` 自动分发到 pFind adapter。

## 质谱学生物学约定

以下常数和公式已经过审计验证（2026-04-07），与 NIST/UniMod 标准一致。

### 质量常数 (`crates/search-engine/src/chemistry.rs`)

| 常数 | 值 (Da) | 来源 |
|------|---------|------|
| PROTON_MASS | 1.007276 | NIST |
| WATER_MASS | 18.010565 | H₂O 单同位素 |
| C13_C12_MASS_DIFF | 1.003355 | ¹³C - ¹²C（`crates/dia-extraction/src/isotope.rs`） |

### 碎片离子公式

- **b 离子**: `b_n = Σ(residue_1..n)` — 不含水
- **y 离子**: `y_n = Σ(residue_{n+1}..end) + H₂O` — 含水（C 端保留 OH，N 端保留 H）
- **m/z 转换**: `ion_mz = (ion_mass + charge × PROTON_MASS) / charge`
- **当前限制**: 仅生成单电荷碎片（b¹⁺, y¹⁺）
- **⚠️ 已知问题**: `matching.rs` 的碎片离子生成不应用固定修饰（`annotate.rs` 版本正确），待修复
- **未实现**: 中性丢失离子（b°/y° -18 Da, b\*/y\* -17 Da）、磷酸化丢失（-98 Da H₃PO₄）

### PPM 计算

```
delta_ppm = (observed - theoretical) / theoretical × 1e6
```

分母始终使用**理论值**（不是观测值）。

### 修饰应用规则

- **固定修饰**: 自动应用到所有目标残基
- **可变修饰**: 组合枚举，受 `max_variable_modifications` 限制（默认 3）
- **N 端修饰**: `AnyNTerm` 应用于所有肽段 N 端；`ProteinNTerm` 仅应用于蛋白质第一条肽段（需要 `DigestedPeptide.is_protein_nterm` 标志）
- **C 端修饰**: 同理，使用 `is_protein_cterm` 标志

### 酶切规则

| 酶 | 规则 | 异常 |
|----|------|------|
| Trypsin | K/R 后切 | P 前不切 |
| Trypsin/P | K/R 后切 | 无异常 |
| Lys-C | K 后切 | — |
| Glu-C | D/E 后切 | — |
| Asp-N | D **前**切 | — |
| Chymotrypsin | F/W/Y/L 后切 | — |

### DIA 检测

- **隔离窗口宽度** = `lower_offset + upper_offset`（总宽度，非半宽）
- **DIA 判定阈值**: 中位窗口宽度 > 5 Da
- **局限**: 5 Da 窄窗口 DIA 会被误判为 DDA（可通过 `acquisition_mode` 手动指定）

---

## 代码审计验证记录

以下问题在代码审计中被标记，经过分析确认为安全/非问题：

### cancel_search 竞态条件 — 安全

`cancel_search()` 持有 Mutex 锁的同时检查并修改状态为 "Cancelled"。后台搜索任务在写入 "Completed" 前也需要获取同一把锁，并会检查 `status == "Cancelled"` 后跳过覆盖。`PanicGuard` 也检查 `status == "Running"` 才覆盖。因此不存在竞态条件。

### 空谱图/空文件处理 — 安全

`Spectrum::validate()` 会拒绝空 peaks 数据（`NoPeaks` 错误）。mzML/MGF reader 在解析后调用 `validate()`，因此不会产生无效谱图进入下游流程。空文件会在 reader 层返回空 Vec，搜索引擎对空输入返回空结果而非 panic。

### Target-Decoy 未实现 — MVP 已知限制

当前 `DecoyStrategy::Reverse` 被接受但实际由搜索引擎内部处理。Rust 侧未实现独立的 decoy 数据库生成。这是 MVP 的已知限制（参见 `tasks/001` 中 FW-6）。

### FDR 硬编码 1% — MVP 已知限制

`report` crate 中 FDR 阈值硬编码为 1%（PSM/肽段/蛋白质）。可配置 FDR 阈值列为未来工作。

### DIA 5 Da 检测阈值 — 已文档化

窄窗口 DIA（<5 Da）可能被误判为 DDA。用户可通过 `acquisition_mode` 参数手动覆盖。
