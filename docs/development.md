# 开发指南

## 项目结构

```text
crates/
├── core/              共享类型（Spectrum, SearchParams, SearchResult 等）
├── spectrum-io/       谱图文件解析（mgf/mzML）
├── param-recommend/   参数推荐规则引擎
├── search-engine/     搜索引擎（SimpleSearch + pFind 预留）
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
