# ProteinCopilot — Project Copilot Instructions

> **一句话定位**：这是一个 **Rust workspace + 多 MCP Server + Agent/Skill 驱动** 的蛋白质组学质谱智能搜索平台。
> Rust 负责所有确定性计算，LLM（通过 Copilot CLI / Claude Desktop 等 MCP Client）负责意图理解、流程编排、参数推荐和结果解释。

---

## 1. 项目架构总览

```text
用户 <─> MCP Client (Copilot CLI / Claude Desktop)
              │
              ├── .github/agents/*.agent.md    ← 领域 Agent 定义（LLM 读取）
              ├── .github/prompts/*.prompt.md  ← 可复用 Skill / Prompt
              │
              └── MCP Server + Library Crates (Rust)  ← 确定性计算能力
                    ├── spectrum-io             ← 谱图读取与解析（lib crate）
                    ├── param-recommend         ← 搜索参数推荐（lib crate）
                    ├── search-engine           ← 搜索引擎调度（lib crate）
                    ├── report                  ← 报告生成（lib crate）
                    ├── core                    ← 共享数据结构与领域模型
                    └── mcp-server              ← MCP Server 组装（bin crate）
```

### 职责分层

| 层 | 职责 | 实现方式 |
|---|---|---|
| **用户交互层** | 接收自然语言指令、展示结果 | MCP Client（Copilot CLI 等） |
| **AI 编排层** | 意图理解、流程规划、参数推荐、结果解释、失败诊断 | Agent.md + Skill Prompt + LLM |
| **MCP Tool 层** | 确定性计算：谱图解析、搜索执行、FDR 计算、报告生成 | Rust MCP Server |
| **搜索引擎层** | 实际质谱搜索执行 | pFind（主）/ MSFragger / Comet（扩展） |

---

## 2. 核心原则

### 2.1 确定性逻辑与 LLM 逻辑必须严格分层

- **Rust MCP Server 只做确定性计算**：谱图解析、数值计算、搜索引擎调用、FDR 统计、报告模板渲染。
- **LLM 只做推理与解释**：用户意图理解、参数推荐理由、结果解释、失败诊断、下一步建议。
- **绝对禁止**：把核心数值计算（FDR、打分、质量偏差计算等）交给 LLM 完成。
- **绝对禁止**：在 Rust MCP Server 中硬编码 LLM 调用逻辑——所有 AI 决策由 Agent 层发起。

### 2.2 Workspace 与 Crate 结构

- 项目是一个 **Rust workspace**，根目录有 `Cargo.toml` workspace 定义。
- 每个 MCP Server 是一个独立 crate（如 `crates/mcp-spectrum-io`）。
- 共享数据结构放在 `crates/core` crate 中。
- 搜索引擎 adapter 放在 `crates/mcp-search-engine/src/adapters/` 下。
- 所有 crate 必须能独立编译和测试。

### 2.3 MCP 协议规范

- 遵循 **Anthropic MCP 标准协议**（JSON-RPC 2.0 over stdio/SSE）。
- 每个 MCP Server 通过 `tools`、`resources`、`prompts` 暴露能力。
- Tool 的输入输出必须是 **JSON Schema** 可描述的结构化数据。
- 每个 Tool 必须有清晰的 `name`、`description`、`inputSchema`，方便 LLM 理解和调用。
- Tool 的 description 必须包含：功能说明、输入要求、输出格式、典型使用场景。

### 2.4 所有 AI 决策输出必须结构化

当 Agent/Skill 产出涉及 AI 推理的结果时（参数推荐、结果解释、失败诊断），输出必须包含：

```json
{
  "decision": "推荐使用 Trypsin 作为消化酶",
  "confidence": 0.92,
  "explanation": "输入数据的末端碎裂模式符合 Trypsin 消化特征...",
  "input_summary": "检测到 12,345 张谱图，平均母离子质量 1,200 Da...",
  "alternatives": ["Lys-C", "Chymotrypsin"],
  "evidence": ["末端碎裂模式分析", "母离子质量分布"]
}
```

这些字段在 `crates/core` 中定义为统一的 Rust 结构体，以保证各模块一致。

### 2.5 可序列化、可审计、可复现

- 所有结果对象必须实现 `Serialize` + `Deserialize`（serde）。
- 每次分析运行必须生成唯一 `run_id`，关联所有中间和最终结果。
- 搜索参数、输入数据摘要、搜索引擎版本、运行时间等元数据必须记录。
- 输出目录结构必须可自描述（包含 `manifest.json` 或 `run_metadata.json`）。

### 2.6 搜索引擎调用必须通过 Adapter 层

```rust
/// 所有搜索引擎必须实现此 trait
pub trait SearchEngineAdapter: Send + Sync {
    /// 执行搜索，返回标准化结果
    async fn search(&self, params: &SearchParams, input_files: &[PathBuf])
        -> Result<SearchResult, CoreError>;
    /// 返回引擎名称和版本
    fn engine_info(&self) -> EngineInfo;
    /// 检查引擎是否可用
    async fn health_check(&self) -> Result<HealthStatus, CoreError>;
}
```

- pFind 是首要 adapter，后续扩展 MSFragger、Comet。
- 搜索结果必须标准化为统一的 `SearchResult` 结构体。
- 不同引擎的原始输出解析逻辑隔离在各自 adapter 内部。

### 2.7 所有外部依赖可替换

- 搜索引擎、LLM Provider、文件格式解析器等外部依赖必须通过 trait 抽象。
- 使用依赖注入模式，不使用全局可变状态。
- 配置文件指定具体实现（哪个搜索引擎、哪个 LLM）。

### 2.8 优先可测试性和可观测性

- 每个 MCP Tool 必须有单元测试和集成测试。
- 使用 trait-based mocking 测试 adapter 层。
- 关键操作（搜索开始/结束、FDR 计算、参数推荐）必须有结构化日志（`tracing` crate）。
- 每个 MCP Tool 调用必须记录：输入摘要、执行时间、输出摘要、错误信息。

---

## 3. 蛋白质组学领域约束

### 3.1 数据格式

- MVP 优先支持 **mzML** 和 **mgf** 格式。
- 同时支持 **DDA** 和 **DIA** 数据采集模式。
- DIA 谱图通过 `IsolationWindow`（target_mz + lower/upper offset）表示宽隔离窗口。
- 谱图数据结构必须在 `crates/core` 中统一定义。
- 质量值（mass_delta、tolerance）使用 `f64` 类型，单位为 Da（道尔顿）。m/z 值为质荷比，无量纲。
- 保留时间统一使用秒（`retention_time_sec`），mzML 中分钟单位（UO:0000031）自动转换。
- 强度值单位为 detector counts。
- 谱图索引从 1 开始（与质谱学惯例一致）。

### 3.2 搜索参数

搜索参数必须结构化表示，至少包含：

- 消化酶（enzyme）
- 固定修饰（fixed modifications）
- 可变修饰（variable modifications）
- 前体离子质量偏差（precursor mass tolerance）
- 碎片离子质量偏差（fragment mass tolerance）
- 漏切位点数（missed cleavages）
- 数据库路径（FASTA database path）

### 3.3 结果标准化

所有搜索引擎的结果必须标准化为统一结构：

- PSM（Peptide-Spectrum Match）级别结果
- 肽段（Peptide）级别结果
- 蛋白质（Protein）级别结果
- 每级结果包含：score、q-value、spectrum reference、modifications

### 3.4 FDR 控制

- FDR 计算是确定性逻辑，必须在 Rust 中实现。
- 支持 target-decoy 策略。
- 默认 FDR 阈值：PSM 1%、肽段 1%、蛋白质 1%。
- FDR 阈值必须可配置。

---

## 4. Agent 与 Skill 设计规范

### 4.1 Agent 定义

- 所有蛋白质组学领域 Agent 放在 `.github/agents/` 目录。
- Agent 的 `description` 必须明确其在蛋白质搜索流程中的角色。
- Agent 的 `tools` 字段列出它可以调用的 MCP Tool。
- Agent 的指令中必须包含：领域知识、决策边界、何时请求用户确认。

### 4.2 Skill / Prompt 定义

- 可复用的分析流程模板放在 `.github/prompts/` 目录。
- Skill 应该覆盖典型场景：基础搜索、磷酸化搜索、开放搜索、多引擎比较等。
- 每个 Skill 必须说明：输入要求、预期输出、适用场景。

### 4.3 AI 编排原则

- LLM 在推荐参数前，必须先通过 MCP Tool 获取数据特征（谱图数量、质量范围等）。
- LLM 在解释结果前，必须先通过 MCP Tool 获取统计摘要（匹配率、FDR 分布等）。
- LLM 不能"凭空"推荐参数或解释结果——必须基于 MCP Tool 返回的数据。
- 用户可以在任何环节要求 LLM 解释"为什么"。

---

## 5. Rust 编码规范

### 5.1 通用 Rust 规范

- 使用 idiomatic Rust，遵循 [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/)。
- 使用 `Result` 和 `Option`，禁止在库代码中使用 `unwrap()` / `expect()`。
- 使用 `thiserror` 定义领域错误类型，`anyhow` 用于顶层应用。
- 使用 `serde` 序列化所有数据结构。
- 使用 `tokio` 作为异步运行时。
- 使用 `tracing` 进行结构化日志记录。
- 使用 `clap` 处理 CLI 参数（如果需要独立运行模式）。

### 5.2 错误处理

- 每个 crate 定义自己的 `Error` 枚举（使用 `thiserror`）。
- MCP Tool 返回的错误必须包含：错误码、人类可读描述、建议操作。
- 错误消息面向蛋白质组学用户：说明发生了什么以及如何修复。

### 5.3 文件结构

```text
proteinCopilot/
├── Cargo.toml                        ← workspace 定义
├── crates/
│   ├── core/                         ← 共享数据结构和 trait
│   │   ├── src/
│   │   │   ├── lib.rs
│   │   │   ├── spectrum.rs           ← 谱图数据结构
│   │   │   ├── search_params.rs      ← 搜索参数
│   │   │   ├── search_result.rs      ← 标准化搜索结果
│   │   │   ├── ai_decision.rs        ← AI 决策输出结构
│   │   │   ├── error.rs              ← 领域错误类型
│   │   │   ├── engine.rs             ← 搜索引擎 Adapter trait
│   │   │   └── run_metadata.rs       ← 运行元数据
│   │   └── Cargo.toml
│   ├── spectrum-io/                  ← 谱图读取（lib crate，非 MCP）
│   │   ├── src/
│   │   │   ├── lib.rs                ← detect_format + create_reader
│   │   │   ├── reader.rs             ← SpectrumReader trait
│   │   │   ├── mgf.rs               ← MGF 解析器
│   │   │   ├── mzml.rs              ← mzML 解析器
│   │   │   └── error.rs             ← SpectrumIoError
│   │   └── Cargo.toml
│   ├── param-recommend/              ← 参数推荐（lib crate）
│   ├── search-engine/                ← 搜索引擎调度（lib crate）
│   │   └── src/adapters/
│   │       ├── mod.rs
│   │       ├── pfind.rs              ← pFind adapter
│   │       ├── msfragger.rs          ← MSFragger adapter（预留）
│   │       └── comet.rs              ← Comet adapter（预留）
│   ├── report/                       ← 报告生成（lib crate）
│   └── mcp-server/                   ← MCP Server 组装（bin crate）
├── .github/
│   ├── agents/                       ← 领域 Agent 定义
│   ├── prompts/                      ← Skill / Prompt 模板
│   └── copilot-instructions.md       ← 本文件
├── tests/
│   ├── integration/                  ← 集成测试
│   └── fixtures/                     ← 测试用谱图数据
└── docs/                             ← 项目文档
```

### 5.4 命名规范

- Crate 名称：功能名（如 `spectrum-io`、`param-recommend`、`search-engine`）。Library crate 无 `mcp-` 前缀；MCP Server 组装 crate 为 `mcp-server`。
- Struct 名称：使用蛋白质组学领域术语（`Spectrum`、`Psm`、`Peptide`、`Protein`）。
- MCP Tool 名称：`动词_名词` 格式（如 `read_spectra`、`run_search`、`calculate_fdr`）。
- 错误变体名称：`XxxError` 格式，包含上下文（如 `SpectrumParseError`、`SearchEngineNotFound`）。

---

## 6. 禁止事项

- ❌ 不允许在 MCP Server 中直接调用 LLM API。
- ❌ 不允许把核心数值计算（FDR、打分、质量计算）交给 LLM。
- ❌ 不允许使用全局可变状态。
- ❌ 不允许硬编码搜索引擎路径或参数——必须通过配置。
- ❌ 不允许忽略编译器警告——CI 中视为错误。
- ❌ 不允许使用 `unsafe` 除非有充分理由和完整文档。
- ❌ 不允许在 MCP Tool 的输入/输出中使用非结构化的自由文本（必须 JSON Schema 可描述）。
- ❌ 不允许跳过错误处理——所有 `Result` 必须被处理。

---

## 7. MVP 阶段范围

### 7.1 必须实现

1. **`core` crate**：谱图、搜索参数、搜索结果、AI 决策等共享数据结构。✅
2. **`spectrum-io`**：读取 mzML / mgf 文件，返回谱图摘要。支持 DDA/DIA。✅
3. **`param-recommend`**：基于谱图特征生成默认参数建议（确定性规则）。
4. **`search-engine`**：通过 pFind adapter 执行搜索。
5. **`report`**：生成结构化搜索结果摘要（供 LLM 解释）。
6. **Agent 定义**：蛋白质搜索助手 Agent（`.github/agents/`）。
7. **Skill 定义**：基础搜索流程 Prompt（`.github/prompts/`）。

### 7.2 预留但不实现

- `mcp-qc`：质控模块
- `mcp-fdr`：FDR 控制（可先由搜索引擎自带 FDR）
- `mcp-protein-inference`：蛋白推断
- MSFragger / Comet adapter
- 多轮分析对话
- 失败诊断

---

## 8. 开发流程规范

### 8.1 新功能开发流程

每次生成新的开发计划或新功能目标时，**必须**按照以下流程执行：

1. **头脑风暴（Brainstorming）**：调用 `brainstorming` skill，使用 Visual Companion 在浏览器中展示架构图、数据流和设计方案，与用户充分讨论并确认设计。
2. **编写设计文档**：将确认后的设计写入 `docs/superpowers/specs/YYYY-MM-DD-<topic>-design.md`。
3. **编写实施计划**：调用 `writing-plans` skill，基于设计文档生成详细的实施计划。
4. **执行实施**：按计划逐步实施，每个任务完成后进行代码审查。

---

## 9. 参考资料

- [MCP 协议规范](https://modelcontextprotocol.io/)
- [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/)
- [mzML 格式规范](https://www.psidev.info/mzML)
- [pFind 搜索引擎](http://pfind.org/)
- [蛋白质组学数据标准 (PSI)](https://www.psidev.info/)
