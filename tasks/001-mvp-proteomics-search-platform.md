# PRD: ProteinCopilot MVP — AI 驱动的蛋白质质谱搜索平台

> **文件名**：`prd-mvp-proteomics-search.md`
> **版本**：1.0
> **创建日期**：2026-03-27
> **状态**：In Progress — M1.1 ✅ M1.2 ✅ M1.3 ✅ M1.4 ✅ M1.5 ✅ M1.6 ✅ M1.7 ✅（357 tests, 0 warnings）— MVP 完成 + Post-MVP 功能（异步搜索优化、搜索历史持久化、谱图注释可视化）+ DIA 数据支持（前体提取 + 搜索集成 + 端到端工作流）

---

## 1. Introduction / Overview

### 1.1 项目简介

ProteinCopilot 是一个 **AI 驱动的蛋白质组学质谱搜索与结果解释平台**。用户通过 MCP Client（如 Copilot CLI、Claude Desktop）用自然语言描述搜索需求，系统自动完成「谱图导入 → 数据质控 → 参数推荐 → 搜索引擎执行 → 结果汇总 → AI 报告」的全流程。

### 1.2 为什么需要这个项目

传统蛋白质组学质谱搜索流程存在以下痛点：

| 痛点 | 说明 |
|---|---|
| **参数选择门槛高** | 新手不知道该用什么消化酶、质量偏差、修饰组合；资深分析员也经常凭经验试错 |
| **流程割裂** | 谱图转换、搜索引擎配置、FDR 控制、报告生成分散在不同工具中，手动串联 |
| **结果解释困难** | 搜索完成后，"为什么鉴定率低？为什么没找到目标蛋白？" 缺乏即时解答 |
| **复现性差** | 参数、版本、中间数据散落各处，难以回溯和复现 |

ProteinCopilot 的价值在于：**让 LLM 做"脑"（意图理解、参数推荐、结果解释），让 Rust MCP Server 做"手"（确定性计算），用户只需用自然语言描述目标，即可获得端到端的搜索结果和专业解读。**

### 1.3 技术架构概要

```text
用户 ─→ MCP Client (Copilot CLI) ─→ LLM (Agent.md 编排)
                                         │
                    ┌────────────────────┼────────────────────┐
                    ▼                    ▼                    ▼
             spectrum-io          param-recommend       search-engine
             (谱图解析)           (参数推荐)            (pFind adapter)
                    │                    │                    │
                    └────────────────────┼────────────────────┘
                                         ▼
                                    report (结果摘要)
                                         │
                                         ▼
                              LLM 生成自然语言报告 → 用户
```

- **Rust Workspace**：多个 library crate + 一个 MCP Server 二进制
- **MCP 协议**：Anthropic 标准 MCP（JSON-RPC 2.0 over stdio）
- **搜索引擎**：pFind 为主（远程服务器/集群，通过 SSH 调用），后续扩展 MSFragger/Comet
- **数据格式**：mzML + mgf

---

## 2. Goals

### 2.1 主要目标

| # | 目标 | 说明 |
|---|---|---|
| G1 | **端到端自动化** | 用户一句自然语言指令即可完成从谱图导入到 AI 报告的完整流程 |
| G2 | **智能参数推荐** | 基于谱图特征自动推荐搜索参数，并用自然语言解释推荐理由 |
| G3 | **结果可解释** | LLM 基于结构化搜索结果生成人类可读的分析报告，解释"为什么" |
| G4 | **覆盖多层次用户** | 新手获得全程引导，资深分析员获得效率提升，根据用户表达自动调整交互深度 |
| G5 | **可复现** | 每次分析运行记录完整元数据（参数、引擎版本、输入摘要），结果可审计可回溯 |

### 2.2 商业价值

- **降低门槛**：湿实验研究员无需学习搜索引擎命令行即可完成蛋白搜索
- **提升效率**：资深分析员减少手动配置和结果解读时间
- **差异化**：市场上没有将 LLM 深度集成到蛋白质组学搜索流程的产品

---

## 3. User Stories

### US-1：新手研究员的基础搜索

> **作为** 一名湿实验研究员，  
> **我希望** 告诉系统"帮我搜一下这批 HeLa 细胞的蛋白数据"，  
> **以便** 不需要了解搜索引擎参数细节就能获得搜索结果。

**验收标准**：
- AC-1.1：用户提供谱图文件路径和 FASTA 数据库路径后，系统自动读取谱图并返回数据摘要
- AC-1.2：系统基于谱图特征推荐搜索参数，并用自然语言解释推荐理由
- AC-1.3：用户确认参数后，系统调用 pFind 执行搜索
- AC-1.4：搜索完成后，系统生成包含鉴定数量、FDR 分布、修饰统计的结构化报告
- AC-1.5：LLM 用自然语言解读报告，说明结果质量是否符合预期

### US-2：资深分析员的精确搜索

> **作为** 一名资深生信分析员，  
> **我希望** 直接指定"用 Trypsin，10ppm precursor tolerance，磷酸化搜索"，  
> **以便** 快速启动搜索而不需要经过推荐流程。

**验收标准**：
- AC-2.1：用户直接指定参数时，系统跳过推荐流程，验证参数合法性后直接执行搜索
- AC-2.2：搜索参数验证不通过时，系统返回具体的错误信息和修复建议
- AC-2.3：搜索结果与手动运行 pFind 结果一致

### US-3：谱图数据探查

> **作为** 研究员，  
> **我希望** 在搜索前先了解谱图数据的基本特征，  
> **以便** 判断数据质量并决定搜索策略。

**验收标准**：
- AC-3.1：系统能读取 mgf 和 mzML 格式文件
- AC-3.2：返回的摘要包含：谱图总数、MS1/MS2 数量、质量范围、保留时间范围、电荷分布
- AC-3.3：单个 mgf 文件 10 万张谱图的摘要生成不超过 30 秒
- AC-3.4：能按 scan number 提取单张谱图的完整数据

### US-4：搜索参数推荐与解释

> **作为** 研究员，  
> **我希望** 系统基于数据特征推荐搜索参数并解释"为什么推荐这些参数"，  
> **以便** 理解推荐逻辑并做出知情决策。

**验收标准**：
- AC-4.1：推荐结果包含 confidence（置信度）、explanation（解释）、alternatives（替代方案）
- AC-4.2：用户通过自然语言补充约束（如"这是磷酸化实验"），推荐结果相应调整
- AC-4.3：提供预设方案列表（标准搜索、磷酸化搜索、TMT 标记等），用户可直接选用
- AC-4.4：LLM 能将规则引擎的模板化解释润色为用户友好的自然语言

### US-5：搜索执行与进度跟踪

> **作为** 研究员，  
> **我希望** 启动搜索后能了解执行进度，  
> **以便** 知道还要等多久。

**验收标准**：
- AC-5.1：搜索引擎通过 SSH 在远程服务器执行，本地 MCP Server 管理连接
- AC-5.2：支持查询搜索进度（进度百分比、已耗时、预估剩余时间）
- AC-5.3：搜索完成后返回标准化 `SearchResult`，不论使用哪个搜索引擎
- AC-5.4：搜索失败时返回错误码、描述和建议操作

### US-6：结果报告与导出

> **作为** 研究员，  
> **我希望** 搜索完成后获得结构化摘要和可下载的结果文件，  
> **以便** 用于后续分析和发表。

**验收标准**：
- AC-6.1：生成的摘要包含：PSM/肽段/蛋白质鉴定数量（1% FDR）、鉴定率、修饰分布、质量偏差分布
- AC-6.2：结果可导出为 TSV 和 JSON 格式
- AC-6.3：导出目录包含 `run_metadata.json`，记录完整的运行元数据（参数、引擎版本、时间等）
- AC-6.4：LLM 基于摘要数据生成自然语言分析报告

### US-7：FASTA 数据库管理

> **作为** 研究员，  
> **我希望** 系统能自动下载常用蛋白数据库或使用我指定的数据库，  
> **以便** 不需要手动管理 FASTA 文件。

**验收标准**：
- AC-7.1：支持用户手动指定 FASTA 文件路径
- AC-7.2：内置常用数据库列表（UniProt Human、Mouse 等），支持按名称选择
- AC-7.3：支持从 UniProt 自动下载并缓存到本地

---

## 4. Functional Requirements

### FR-1：谱图读取模块（spectrum-io）

| ID | 需求 | 优先级 |
|---|---|---|
| FR-1.1 | 解析 mgf 格式文件，提取谱图（scan number, precursor, m/z array, intensity array, retention time） | P0 |
| FR-1.2 | 解析 mzML 格式文件，支持 base64 + zlib 压缩的 binary data array | P0 |
| FR-1.3 | 自动检测文件格式（通过扩展名 + 文件头 magic bytes） | P0 |
| FR-1.4 | 生成 `SpectrumSummary`（谱图总数、MS 级别分布、质量范围、RT 范围、电荷分布、中位峰数） | P0 |
| FR-1.5 | 按 scan number 提取单张谱图完整数据 | P1 |
| FR-1.6 | 使用 streaming 解析，支持 10 万张谱图级别的文件而不会内存溢出 | P0 |

**MCP Tools**：
- `read_spectra`：输入 file_path → 输出 SpectrumSummary
- `get_spectrum`：输入 file_path + scan_number → 输出 Spectrum

### FR-2：参数推荐模块（param-recommend）

| ID | 需求 | 优先级 |
|---|---|---|
| FR-2.1 | 基于 SpectrumSummary 的质量范围和碎裂特征推断仪器类型 | P0 |
| FR-2.2 | 根据仪器类型推荐 precursor/fragment mass tolerance | P0 |
| FR-2.3 | 提供默认消化酶推荐（Trypsin 为默认） | P0 |
| FR-2.4 | 支持常见修饰组合预设（标准搜索、磷酸化、TMT、SILAC） | P0 |
| FR-2.5 | 接受 `UserHints`（experiment_type, instrument_type, custom_notes）调整推荐 | P0 |
| FR-2.6 | 输出封装为 `AiDecision<SearchParams>`（decision + confidence + explanation + alternatives + evidence） | P0 |
| FR-2.7 | 所有推荐逻辑为确定性规则（相同输入 → 相同输出），不调用 LLM | P0 |

**MCP Tools**：
- `recommend_params`：输入 SpectrumSummary + UserHints? → 输出 AiDecision\<SearchParams\>
- `list_presets`：无输入 → 输出 Vec\<SearchPreset\>

### FR-3：搜索引擎调度模块（search-engine）

| ID | 需求 | 优先级 |
|---|---|---|
| FR-3.1 | 实现 `SearchEngineAdapter` trait，定义统一的搜索接口 | P0 |
| FR-3.2 | 实现 pFind adapter：将 SearchParams 转换为 pFind .cfg 配置格式 | P0 |
| FR-3.3 | 通过 SSH 连接远程服务器执行 pFind（`tokio::process::Command` + SSH） | P0 |
| FR-3.4 | 解析 pFind 输出结果文件，转换为标准化 `SearchResult` | P0 |
| FR-3.5 | 实现搜索进度跟踪（轮询 pFind 日志/进度文件） | P1 |
| FR-3.6 | 实现 health_check：检查远程 pFind 是否可达、版本信息 | P0 |
| FR-3.7 | 实现 `EngineRegistry`：管理多个搜索引擎 adapter | P1 |
| FR-3.8 | SearchParams 在执行前必须通过 `validate()` 校验 | P0 |

**MCP Tools**：
- `run_search`：输入 SearchParams + input_files + engine? → 输出 SearchResult
- `check_engine`：无输入 → 输出 Vec\<(EngineInfo, HealthStatus)\>
- `get_search_status`：输入 run_id → 输出 进度信息

### FR-4：报告生成模块（report）

| ID | 需求 | 优先级 |
|---|---|---|
| FR-4.1 | 从 SearchResult 生成 SearchResultSummary（统计聚合） | P0 |
| FR-4.2 | 摘要包含：PSM/肽段/蛋白数量（1% FDR）、鉴定率、中位 score、中位质量偏差、修饰分布、电荷分布 | P0 |
| FR-4.3 | 导出 SearchResult 为 TSV 格式（PSM 级别、肽段级别、蛋白质级别各一个文件） | P0 |
| FR-4.4 | 导出 SearchResult 为 JSON 格式 | P1 |
| FR-4.5 | 每次导出同时生成 `run_metadata.json`（run_id, params, engine_info, duration, timestamps） | P0 |
| FR-4.6 | 比较两次搜索结果的差异（鉴定数变化、overlap） | P2 |

**MCP Tools**：
- `generate_summary`：输入 SearchResult → 输出 SearchResultSummary
- `export_results`：输入 SearchResult + format + output_path → 输出 ExportResult

### FR-5：Agent 与 Skill 定义

| ID | 需求 | 优先级 |
|---|---|---|
| FR-5.1 | 定义 `proteomics-search.agent.md`：蛋白质搜索助手，编排完整搜索流程 | P0 |
| FR-5.2 | 定义 `basic-search.prompt.md` Skill：step-by-step 基础搜索流程 | P0 |
| FR-5.3 | 定义 `result-interpretation.prompt.md` Skill：引导 LLM 解读搜索结果 | P0 |
| FR-5.4 | Agent 必须在执行搜索前请求用户确认参数 | P0 |
| FR-5.5 | Agent 必须先调用 MCP Tool 获取数据再做推荐/解释，不得"凭空"推断 | P0 |

### FR-6：FASTA 数据库管理

| ID | 需求 | 优先级 |
|---|---|---|
| FR-6.1 | 支持用户指定本地 FASTA 路径 | P0 |
| FR-6.2 | 内置常用数据库列表（UniProt Human、Mouse、E.coli 等） | P1 |
| FR-6.3 | 支持从 UniProt FTP 自动下载并缓存 | P2 |

### FR-7：共享数据结构（core crate）

| ID | 需求 | 优先级 |
|---|---|---|
| FR-7.1 | 定义 Spectrum, SpectrumSummary, PrecursorInfo, MsLevel | P0 |
| FR-7.2 | 定义 SearchParams, Enzyme, Modification, MassTolerance | P0 |
| FR-7.3 | 定义 PSM, PeptideResult, ProteinResult, SearchResult, SearchResultSummary | P0 |
| FR-7.4 | 定义 AiDecision\<T\> 泛型包装器 | P0 |
| FR-7.5 | 定义 SearchEngineAdapter trait (async) | P0 |
| FR-7.6 | 定义 RunMetadata, CoreError | P0 |
| FR-7.7 | 所有类型 derive Serialize + Deserialize + Debug + Clone + JsonSchema | P0 |

---

## 5. Non-Goals（MVP 明确排除）

以下功能在 MVP 中**不实现**，但架构设计预留扩展空间：

| # | 排除项 | 未来阶段 |
|---|---|---|
| NG-1 | 独立的 FDR 控制模块（`mcp-fdr`）——MVP 依赖 pFind 自带的 FDR 计算 | Phase 2 |
| NG-2 | 数据质控模块（`mcp-qc`）——谱图质量评估和异常检测 | Phase 2 |
| NG-3 | 蛋白推断模块（`mcp-protein-inference`）——Parsimony 算法 | Phase 2 |
| NG-4 | MSFragger / Comet 等其他搜索引擎 adapter | Phase 2 |
| NG-5 | 搜索失败的自动诊断 Agent——诊断参数不合理、数据质量差等原因 | Phase 2 |
| NG-6 | 多轮分析对话——在一次会话中迭代调参、重搜、比较 | Phase 3 |
| NG-7 | 开放搜索 / 未知修饰发现 | Phase 3 |
| NG-8 | 定量分析支持（Label-free / TMT / iTRAQ） | Phase 3 |
| NG-9 | 谱图可视化（SVG/HTML）和交互式报告 | Phase 3 |
| NG-10 | Web UI / GUI——MVP 仅通过 MCP Client（Copilot CLI 等）交互 | 未规划 |

---

## 6. Design Considerations

### 6.1 交互设计（通过 MCP Client + Agent.md）

- **自适应交互深度**：Agent 根据用户表达的专业程度自动调整，新手获得详细解释和确认步骤，资深用户获得简洁高效的响应
- **关键决策必须确认**：搜索参数推荐后，Agent 必须向用户展示参数并等待确认，不得自动执行搜索
- **渐进式披露**：默认展示摘要信息，用户可追问"为什么"获取详细解释
- **错误信息面向领域用户**：错误消息说明发生了什么、可能的原因、建议的修复操作

### 6.2 典型交互流程

```text
用户：帮我搜一下这批 HeLa 的磷酸化数据
  │
  ├─ Agent：好的，请提供谱图文件路径
  │  用户：/data/hela_phospho.mgf
  │
  ├─ Agent 调用 read_spectra → 返回 SpectrumSummary
  │  Agent：检测到 98,432 张 MS2 谱图，质量范围 350-2000 Da，保留时间 0-120 分钟
  │
  ├─ Agent 调用 recommend_params(summary, {experiment_type: "phosphorylation"})
  │  Agent：推荐参数：Trypsin, 10ppm precursor, 0.02Da fragment,
  │         固定修饰 Carbamidomethyl(C), 可变修饰 Phospho(STY) + Oxidation(M)
  │         置信度 0.88，理由：质量范围和碎裂模式符合 Orbitrap 高分辨仪器特征...
  │  用户：用这个参数，数据库用 UniProt Human
  │
  ├─ Agent 调用 run_search(params, files)
  │  Agent：搜索已提交到远程服务器，预计 15-30 分钟
  │  ... (用户可查询进度) ...
  │  Agent：搜索完成！
  │
  ├─ Agent 调用 generate_summary(search_result)
  │  Agent：鉴定 12,456 个 PSM (1% FDR)，8,234 个肽段，1,567 个蛋白质
  │         鉴定率 12.7%，磷酸化位点 2,345 个（pS:pT:pY = 85:12:3）
  │         这个结果质量较好，磷酸化位点分布符合预期...
  │
  └─ 用户可继续追问或调参重搜
```

### 6.3 数据流约束

- LLM 向 MCP Tool 传递的必须是**结构化数据**（JSON），不是自由文本
- MCP Tool 返回的必须是**结构化结果**（JSON Schema 可描述），不是自由文本
- LLM 做解释时，必须**先调用 Tool 获取数据**，不得"凭空"推断

---

## 7. Technical Considerations

### 7.1 架构（详见 `docs/architecture.md`）

- **Rust Workspace**：单 MCP Server 二进制（`mcp-server`） + 多 library crate（`core`, `spectrum-io`, `param-recommend`, `search-engine`, `report`）
- **MCP SDK**：`rmcp` crate（v0.16+），使用 `#[tool]` / `#[tool_router]` 宏
- **async runtime**：`tokio`
- **确定性/LLM 分层**：library crate 不依赖 rmcp 和 LLM，mcp-server 只做薄包装

### 7.2 性能

| 场景 | 要求 |
|---|---|
| mgf 摘要生成（10 万张谱图） | ≤ 30 秒 |
| mzML 摘要生成（10 万张谱图） | ≤ 60 秒 |
| 参数推荐（规则引擎） | ≤ 1 秒 |
| pFind 搜索 | 取决于远程服务器，MCP Server 不阻塞 |
| 结果摘要生成 | ≤ 5 秒 |
| TSV 导出（10 万 PSM） | ≤ 10 秒 |

- 谱图文件使用 **streaming 解析**（逐条读取），不一次性加载到内存
- 远程搜索使用 **异步执行**，支持进度查询

### 7.3 远程搜索引擎连接

- pFind 部署在远程 Linux 服务器/集群
- MCP Server 通过 **SSH**（`tokio::process::Command` + ssh）连接远程执行
- SSH 配置（host, user, key, pFind 路径等）通过配置文件或环境变量指定
- 文件传输：谱图和数据库文件假设已在远程服务器上，MCP 传递远程路径
- 结果回传：通过 SSH + scp/rsync 拉取结果文件到本地解析

### 7.4 错误处理

- 每个 crate 定义领域 Error 枚举（`thiserror`）
- MCP Tool 错误响应包含：error_code（机器可读）、message（人类可读）、suggestion（建议操作）
- 错误消息面向蛋白质组学用户，例如：
  ```json
  {
    "error_code": "SEARCH_ENGINE_UNREACHABLE",
    "message": "无法连接到 pFind 服务器 (gpu-server-01)",
    "suggestion": "请检查 SSH 连接配置和服务器状态，运行 check_engine 查看详情"
  }
  ```

### 7.5 可复现性

- 每次搜索运行生成唯一 `run_id`（UUID v4）
- `RunMetadata` 记录：run_id, input_files, params_used, engine_info (name + version), start_time, end_time, duration, status
- 导出结果目录必须包含 `run_metadata.json`，实现自描述

### 7.6 可测试性

| 层 | 测试策略 |
|---|---|
| core | 单元测试：serde round-trip、参数验证、类型构造 |
| spectrum-io | 单元测试 + fixture 文件：mgf/mzML 解析正确性 |
| param-recommend | 单元测试：不同谱图特征 → 推荐结果验证 |
| search-engine | Mock adapter 测试 tool 逻辑；集成测试需 pFind 环境 |
| report | 单元测试：构造 SearchResult → 验证摘要/导出正确性 |
| mcp-server | 集成测试：JSON-RPC 请求 → 响应验证 |

### 7.7 依赖管理

核心 Rust 依赖：

| Crate | 用途 | 范围 |
|---|---|---|
| `serde` + `serde_json` | 序列化 | 全局 |
| `schemars` | JSON Schema 生成 | core |
| `thiserror` | 错误定义 | 全局 |
| `anyhow` | 顶层错误处理 | mcp-server |
| `tokio` | async runtime | mcp-server, search-engine |
| `tracing` | 结构化日志 | mcp-server |
| `rmcp` | MCP SDK | mcp-server |
| `uuid` | run_id 生成 | core |
| `chrono` | 时间戳 | core |
| `quick-xml` | mzML 解析 | spectrum-io |
| `base64` | mzML binary decode | spectrum-io |
| `flate2` | zlib 解压 | spectrum-io |
| `csv` | TSV 导出 | report |

---

## 8. Success Metrics

### 8.1 MVP 验收指标

| # | 指标 | 目标 | 验证方式 |
|---|---|---|---|
| SM-1 | 端到端流程可跑通 | 从谱图导入到 AI 报告的完整流程至少成功运行 1 次 | 手动测试 |
| SM-2 | mgf 解析正确性 | 解析结果与第三方工具（如 pyteomics）一致 | 对比测试 |
| SM-3 | mzML 解析正确性 | 解析结果与第三方工具一致 | 对比测试 |
| SM-4 | 参数推荐覆盖度 | 标准搜索和磷酸化搜索两种场景推荐结果合理 | 领域专家评审 |
| SM-5 | pFind 搜索结果一致性 | MCP 调用 pFind 的结果与手动运行 pFind 结果一致 | 对比测试 |
| SM-6 | 结果摘要准确性 | 摘要中的统计数字（PSM数、蛋白数等）与原始结果一致 | 单元测试 |
| SM-7 | MCP Server 可被 Copilot CLI 连接 | 通过 stdio transport 成功连接并调用 tool | 手动测试 |
| SM-8 | 所有 crate 编译通过 | `cargo build --workspace` 无 error/warning | CI |
| SM-9 | 单元测试覆盖率 | core + spectrum-io + param-recommend 测试通过 | `cargo test` |

### 8.2 质量指标

- `cargo clippy --workspace` 零 warning
- `cargo fmt --check` 通过
- 所有公开 API 有 `///` doc 注释

---

## 9. Open Questions

| # | 问题 | 影响 | 当前假设 |
|---|---|---|---|
| OQ-1 | pFind 的 SSH 远程调用具体协议？是否需要作业调度系统（SLURM/PBS）？ | search-engine adapter 设计 | 假设直接 SSH + 命令行调用 |
| OQ-2 | pFind 的输出文件格式细节（哪些文件、字段定义）？ | 结果解析 parser 实现 | 需要获取 pFind 输出样例 |
| OQ-3 | pFind 配置文件 (.cfg) 的完整字段规范？ | SearchParams → .cfg 转换逻辑 | 需要获取 pFind 文档 |
| OQ-4 | 谱图文件和 FASTA 是否已在远程服务器上？还是需要 MCP Server 上传？ | 文件传输设计 | 假设已在远程 |
| OQ-5 | UniProt 自动下载功能的缓存策略和更新频率？ | FASTA 管理模块 | MVP 先做手动指定，自动下载 P2 |
| OQ-6 | 多文件分析时，是合并后搜索还是分别搜索？ | run_search 的输入设计 | 假设多文件合并为一次搜索 |
| OQ-7 | `rmcp` crate 的 `#[tool]` 宏是否支持复杂嵌套类型作为参数？ | MCP Tool 输入设计 | 需要 POC 验证 |

---

## 10. Implementation Milestones

### 关键路径

```text
M1.1 (core) ──→ M1.2 (spectrum-io) ──→ M1.3 (param-recommend) ──→ M1.4 (search-engine) ──→ M1.5 (report) ──→ M1.6 (agent/skill) ──→ M1.7 (集成验证)
```

### 依赖关系

```text
M1.1 (core) ──┬──→ M1.2 (spectrum-io)
              ├──→ M1.3 (param-recommend)  ← 需要 M1.2 的 SpectrumSummary
              ├──→ M1.4 (search-engine)    ← 需要 M1.2 + M1.3
              └──→ M1.5 (report)           ← 需要 M1.4
                                                │
M1.6 (agent/skill) ← 需要 M1.2 ~ M1.5 的 Tool 定义
                                                │
M1.7 (集成验证)    ← 需要所有 MVP Milestone
```

---

### Milestone 1.1：项目骨架与 `core` crate

> **状态**：✅ 已完成（160 tests, 0 clippy warnings）

> 建立 Rust Workspace 结构，定义所有共享数据类型和 trait，为后续 MCP crate 提供基础。
> 关联 FR：FR-7.1 ~ FR-7.7

#### Task 1.1.1：初始化 Rust Workspace ✅

- **Sub-task 1.1.1.1**：创建根 `Cargo.toml`（workspace 定义），配置 `crates/*` 成员
- **Sub-task 1.1.1.2**：创建 `crates/core/Cargo.toml` 和 `crates/core/src/lib.rs`
- **Sub-task 1.1.1.3**：添加基础依赖：`serde`, `serde_json`, `schemars`, `thiserror`, `uuid`, `chrono`
- **Sub-task 1.1.1.4**：配置 `rustfmt.toml` 和 `.clippy.toml`，确保代码风格一致
- **Sub-task 1.1.1.5**：验证 `cargo build` 和 `cargo test` 通过

#### Task 1.1.2：定义谱图数据结构（`spectrum.rs`）✅

- **Sub-task 1.1.2.1**：定义 `MsLevel` 枚举（MS1, MS2, Other(u8)）
- **Sub-task 1.1.2.2**：定义 `PrecursorInfo` 结构体（mz, charge, intensity, isolation_window: Option\<IsolationWindow\>）
- **Sub-task 1.1.2.3**：定义 `Spectrum` 结构体（scan_number, ms_level, retention_time_sec, precursors: Vec\<PrecursorInfo\>, mz_array, intensity_array）
  - DDA: 1 precursor, DIA: 0~1 with wide isolation window, MS1: empty
  - `IsolationWindow`：target_mz, lower_offset, upper_offset（对齐 mzML）
- **Sub-task 1.1.2.4**：定义 `SpectrumSummary` 结构体（file_path, format, total_spectra, ms1_count, ms2_count, mz_range, rt_range_sec, precursor_charge_distribution, median_peaks_per_spectrum）
- **Sub-task 1.1.2.5**：定义 `SpectrumFileInfo` 结构体（path, format, file_size_bytes），`SpectrumFormat` 枚举（MzML, Mgf）
- **Sub-task 1.1.2.6**：所有结构体 derive `Serialize`, `Deserialize`, `Debug`, `Clone`, `JsonSchema`
- **Sub-task 1.1.2.7**：编写单元测试：序列化/反序列化 round-trip 测试
- **实现补充**：`Spectrum::new()` 构造器 + `validate()` 方法：mz/intensity 数组长度一致、mz 升序、m/z > 0、intensity ≥ 0、scan_number ≥ 1、precursor 字段校验。`SpectrumSummary::validate()` 校验范围一致性。`SpectrumError` 枚举（8 变体）。

#### Task 1.1.3：定义搜索参数结构（`search_params.rs`）✅

- **Sub-task 1.1.3.1**：定义 `Enzyme` 枚举（Trypsin, LysC, GluC, AspN, Chymotrypsin, TrypsinP, NonSpecific, Custom { name, cleavage_rule }）
- **Sub-task 1.1.3.2**：定义 `ModPosition` 枚举（Anywhere, AnyNTerm, AnyCTerm, ProteinNTerm, ProteinCTerm）
- **Sub-task 1.1.3.3**：定义 `Modification` 结构体（name, mass_delta, residues: Vec\<char\>, position: ModPosition）
- **Sub-task 1.1.3.4**：定义 `ToleranceUnit` 枚举（Ppm, Da），`MassTolerance` 结构体（value, unit）
- **Sub-task 1.1.3.5**：定义 `DecoyStrategy` 枚举（Reverse, Shuffle, None）
- **Sub-task 1.1.3.6**：定义 `SearchParams` 结构体（enzyme, missed_cleavages, fixed_modifications, variable_modifications, precursor_tolerance, fragment_tolerance, database_path: String, decoy_strategy: DecoyStrategy）
- **Sub-task 1.1.3.7**：实现 `SearchParams::validate()` → `Result<(), SearchParamsError>`：检查 database_path 非空、tolerance > 0 且 finite、missed_cleavages ≤ 5、modification mass_delta finite
- **Sub-task 1.1.3.8**：编写单元测试：合法参数通过、非法参数返回明确错误
- **实现补充**：`SearchParams` 额外 derive `PartialEq`（支持一致性校验）。`SearchParamsError` 枚举（4 变体）。

#### Task 1.1.4：定义搜索结果结构（`search_result.rs`）✅

- **Sub-task 1.1.4.1**：定义 `Psm` 结构体（spectrum_scan, peptide_sequence, modifications, charge, precursor_mz, calculated_mz, delta_mass_ppm, score, q_value, protein_accessions, is_decoy）
- **Sub-task 1.1.4.2**：定义 `PeptideResult` 结构体（sequence, protein_accessions, best_score, q_value, psm_count）
- **Sub-task 1.1.4.3**：定义 `ProteinResult` 结构体（accession, description, coverage, peptide_count, unique_peptide_count）
- **Sub-task 1.1.4.4**：定义 `SearchResultSummary` 结构体（total_spectra_searched, total_psms, psms_at_1pct_fdr, unique_peptides_at_1pct_fdr, protein_groups_at_1pct_fdr, median_score, median_delta_mass_ppm, identification_rate, modification_distribution, charge_distribution, search_duration_sec）
- **Sub-task 1.1.4.5**：定义 `SearchResult` 结构体（run_id: Uuid, engine_info: EngineInfo, params_used: SearchParams, psms, peptides, proteins, summary, metadata: RunMetadata）
- **Sub-task 1.1.4.6**：编写单元测试：构造完整 SearchResult → JSON round-trip
- **实现补充**：所有类型均有 `validate()` 方法。`SearchResultError` 枚举（10 变体）。`SearchResult::validate()` 校验 run_id/engine_info/params_used 与 metadata 一致性，并委托到 summary.validate() + metadata.validate()。

#### Task 1.1.5：定义 AI 决策输出结构（`ai_decision.rs`）✅

- **Sub-task 1.1.5.1**：定义 `AiDecision<T>` 泛型结构体（decision: T, confidence: f64, explanation: String, input_summary: String, alternatives: Vec\<String\>, evidence: Vec\<String\>）
- **Sub-task 1.1.5.2**：编写单元测试：`AiDecision<String>` 和 `AiDecision<SearchParams>` 的 JSON 输出格式符合 copilot-instructions.md §2.4 规范
- **实现补充**：`AiDecision::validate()` 校验 confidence∈[0,1]、explanation/input_summary 非空。`AiDecisionError` 枚举（2 变体）。

#### Task 1.1.6：定义领域错误类型（`error.rs`）✅

- **Sub-task 1.1.6.1**：定义 `CoreError` 枚举，使用 `#[derive(thiserror::Error)]`：
  - `SpectrumParseError { format: String, detail: String, suggestion: String }`
  - `InvalidSearchParams { field: String, reason: String, suggestion: String }`
  - `SearchEngineError { engine: String, detail: String, suggestion: String }`
  - `FileNotFound { path: PathBuf }`
  - `UnsupportedFormat { format: String, supported: Vec<String> }`
  - `SshConnectionError { host: String, detail: String }`
  - `ResultParseError { engine: String, file: PathBuf, detail: String }`
  - `ValidationError { context: String, detail: String, suggestion: String }`
- **Sub-task 1.1.6.2**：实现 `CoreError::suggestion(&self) -> &str` 方法，返回建议操作
- **实现补充**：`ErrorReport` 可序列化错误摘要（category, message, suggestion）用于 MCP 响应。`From` impl 覆盖所有 5 个模块错误类型 → CoreError。

#### Task 1.1.7：定义搜索引擎 Adapter trait（`engine.rs`）✅

- **Sub-task 1.1.7.1**：定义 `SearchEngineAdapter` async trait：
  - `async fn search(&self, params: &SearchParams, input_files: &[PathBuf]) -> Result<SearchResult, CoreError>`
  - `fn engine_info(&self) -> EngineInfo`
  - `async fn health_check(&self) -> Result<HealthStatus, CoreError>`
- **Sub-task 1.1.7.2**：定义 `EngineInfo` 结构体（name: String, version: String, supported_features: Vec\<String\>）
- **Sub-task 1.1.7.3**：定义 `HealthStatus` 枚举（Healthy, Degraded { reason }, Unavailable { reason }）

#### Task 1.1.8：定义运行元数据结构（`run_metadata.rs`）✅

- **Sub-task 1.1.8.1**：定义 `RunStatus` 枚举（Pending, Running, Completed, Failed { reason }）
- **Sub-task 1.1.8.2**：定义 `RunMetadata` 结构体（run_id: Uuid, created_at: DateTime, input_files: Vec\<PathBuf\>, params_used: SearchParams, engine_info: EngineInfo, duration_sec: Option\<f64\>, status: RunStatus）
- **Sub-task 1.1.8.3**：实现 `RunMetadata::new(params, engine_info, input_files)` 自动生成 run_id + timestamp
- **实现补充**：`RunMetadata::validate()` 校验 duration_sec finite+非负、input_files 非空、engine_info 非空，并委托到 `params_used.validate()`。`RunMetadataError` 枚举（4 变体）。

---

### Milestone 1.2：`spectrum-io` — 谱图读取 Library Crate

> **状态**：✅ 已完成（46 tests, 0 clippy warnings，已通过 165MB 真实 DIA mzML 验证）
>
> 实现谱图文件解析库，支持 mgf 和 mzML。作为 library crate，不依赖 MCP。
> 关联 FR：FR-1.1 ~ FR-1.6 | 关联 US：US-3

#### Task 1.2.1：创建 crate 并定义 Reader trait ✅

- **Sub-task 1.2.1.1**：创建 `crates/spectrum-io/Cargo.toml`，依赖 `core`, `quick-xml`, `base64`, `flate2`
- **Sub-task 1.2.1.2**：定义 `SpectrumReader` trait（`read_all`, `read_summary`, `read_spectrum`）
- **Sub-task 1.2.1.3**：实现 `detect_format(path) -> Result<SpectrumFileInfo>` 和 `create_reader(info) -> Box<dyn SpectrumReader>`
- **实现补充**：定义 `SpectrumIoError` 枚举（9 变体）+ `From<SpectrumIoError> for CoreError` 转换。

#### Task 1.2.2：实现 mgf 格式解析 ✅

- **Sub-task 1.2.2.1**：实现 `MgfReader` 结构体，实现 `SpectrumReader` trait
- **Sub-task 1.2.2.2**：实现 streaming 解析器：逐行读取，遇到 `BEGIN IONS` 开始收集，`END IONS` 结束
- **Sub-task 1.2.2.3**：解析 header 字段：TITLE, PEPMASS（mz + intensity）, CHARGE（支持 `2+`/`3-`/`2` 格式）, RTINSECONDS, SCANS
- **Sub-task 1.2.2.4**：解析 m/z + intensity 行对（空格或 tab 分隔）
- **Sub-task 1.2.2.5**：实现 `read_summary`：streaming 遍历所有谱图，统计不加载全部到内存
- **Sub-task 1.2.2.6**：准备 fixture 文件：`tests/fixtures/small.mgf`（10 张谱图）
- **Sub-task 1.2.2.7**：编写单元测试：解析 fixture → 验证 scan count、precursor mz、peak count
- **实现补充**：解析后自动按 m/z 升序排序（真实数据可能无序）。`read_spectrum` 支持 streaming + 找到即停。`read_summary` 返回前调用 `validate()`。

#### Task 1.2.3：实现 mzML 格式解析 ✅

- **Sub-task 1.2.3.1**：实现 `MzMLReader` 结构体，实现 `SpectrumReader` trait
- **Sub-task 1.2.3.2**：使用 `quick-xml` 的 event-based reader，streaming 解析 `<spectrum>` 元素
- **Sub-task 1.2.3.3**：从 `<cvParam>` 提取：ms level (MS:1000511), scan number, retention time (MS:1000016)
- **Sub-task 1.2.3.4**：从 `<precursor>/<selectedIon>` 提取：precursor mz, charge, intensity；从 `<isolationWindow>` 提取 DIA 隔离窗口
- **Sub-task 1.2.3.5**：解码 `<binaryDataArray>`：读取 base64 → 解码 → 如有 zlib 压缩则解压 → 解析为 f64 数组
- **Sub-task 1.2.3.6**：通过 `<cvParam>` 判断 precision（32-bit or 64-bit float）和压缩方式
- **Sub-task 1.2.3.7**：准备 fixture 文件：`tests/fixtures/small.mzml`（10 张谱图，scan 9 使用 zlib 压缩）
- **Sub-task 1.2.3.8**：编写单元测试：解析 fixture → 验证谱图数量和关键字段
- **实现补充**：RT 单位自动转换（UO:0000031 分钟 → 秒）。自动 m/z 排序。DIA isolation window 完整提取。

#### Task 1.2.4：测试与基准 ✅

- **Sub-task 1.2.4.1**：跨格式一致性测试：mgf 和 mzml 解析同一数据集，验证结果一致（scan number、peak count、m/z、intensity、charge）
- **Sub-task 1.2.4.2**：错误场景测试：不存在的文件、损坏的 mgf、不完整的 mzML、无 binary array 的 mzML、无效 XML
- **Sub-task 1.2.4.3**：真实数据验证：165MB HeLa DIA mzML（2000 spectra, 31K peaks/spectrum）解析成功
- **实现补充**：添加 `examples/read_spectra.rs` CLI 工具用于快速验证任意文件。

---

### Milestone 1.3：`param-recommend` — 参数推荐 Library Crate

> 确定性规则引擎，基于谱图特征推荐搜索参数。不依赖 MCP 和 LLM。
> 关联 FR：FR-2.1 ~ FR-2.7 | 关联 US：US-4

#### Task 1.3.1：创建 crate

- **Sub-task 1.3.1.1**：创建 `crates/param-recommend/Cargo.toml`，仅依赖 `core`
- **Sub-task 1.3.1.2**：创建 `src/lib.rs`，定义模块结构

#### Task 1.3.2：定义 UserHints 和 SearchPreset

- **Sub-task 1.3.2.1**：定义 `UserHints` 结构体（experiment_type: Option\<String\>, instrument_type: Option\<String\>, custom_notes: Option\<String\>），derive JsonSchema
- **Sub-task 1.3.2.2**：定义 `SearchPreset` 结构体（name, description, params: SearchParams, applicable_scenarios: Vec\<String\>）
- **Sub-task 1.3.2.3**：实现内置预设列表：
  - `standard`：标准蛋白搜索（Trypsin, Carbamidomethyl(C), Oxidation(M), 10ppm/0.02Da）
  - `phospho`：磷酸化搜索（+ Phospho(STY)）
  - `tmt`：TMT 标记搜索（+ TMT6plex(K, N-term)）
  - `open`：开放搜索（500Da precursor tolerance）

#### Task 1.3.3：实现推荐规则引擎

- **Sub-task 1.3.3.1**：实现 `ParamRecommender::recommend(summary, hints) -> AiDecision<SearchParams>`
  - 输入 `SpectrumSummary` 为空时（`is_empty()` 返回 true），返回明确错误或默认参数 + 低 confidence 说明
- **Sub-task 1.3.3.2**：仪器类型推断规则：
  - mz_range 上限 > 1500 + median_peaks > 200 → Orbitrap/高分辨 → 10ppm / 0.02Da
  - mz_range 上限 < 1500 + median_peaks < 100 → TOF/低分辨 → 20ppm / 0.1Da
  - 否则 → 通用默认 → 15ppm / 0.05Da
- **Sub-task 1.3.3.3**：消化酶推荐：默认 Trypsin，hints 中指定则覆盖
- **Sub-task 1.3.3.4**：修饰推荐：根据 experiment_type 选择预设修饰组合
- **Sub-task 1.3.3.5**：生成 explanation 模板文字，例如：「根据质量范围 [350-2000 Da] 和中位峰数 [256 peaks/spectrum]，推断为高分辨 Orbitrap 仪器，推荐 precursor tolerance 10 ppm」
- **Sub-task 1.3.3.6**：生成 alternatives 列表和 evidence 列表
- **Sub-task 1.3.3.7**：confidence 计算：有 user_hints 命中 → 0.90+，纯自动推断 → 0.70-0.85
- **Sub-task 1.3.3.8**：SearchParams 语义冲突校验：检测参数组合矛盾（如 NonSpecific 酶 + missed_cleavages > 0、开放搜索 + 窄 fragment tolerance 等），在 explanation 中给出警告

#### Task 1.3.4：实现 list_presets

- **Sub-task 1.3.4.1**：实现 `ParamRecommender::list_presets() -> Vec<SearchPreset>`

#### Task 1.3.5：测试

- **Sub-task 1.3.5.1**：测试高分辨数据输入 → 推荐 10ppm
- **Sub-task 1.3.5.2**：测试低分辨数据输入 → 推荐 20ppm
- **Sub-task 1.3.5.3**：测试 user_hints "phosphorylation" → 推荐包含 Phospho(STY) 修饰
- **Sub-task 1.3.5.4**：测试 confidence 范围合理性
- **Sub-task 1.3.5.5**：验证相同输入 → 完全相同输出（确定性）

#### M1.3 已知局限性（MVP 阶段有意简化，供后续迭代参考）

1. **仅基于统计特征推断**：只使用 m/z 范围和中位 peak 数两个特征推断仪器类型，未利用 charge 分布、RT 分布等其他维度
2. **无碎裂模式分析**：无法从 b/y 离子分布推断消化酶或碎裂方式（HCD/CID/ETD），酶推荐依赖用户 hints 或默认 Trypsin
3. **无 DDA/DIA 自动区分**：不分析 isolation window 宽度来判断采集模式，DIA 数据（如宽窗口 25 Da）与 DDA 数据使用相同推荐逻辑
4. **硬编码阈值规则**：仪器推断使用固定评分阈值，非机器学习模型，边界情况可能误分类
5. **database_path 占位符**：推荐结果中 FASTA 数据库路径为占位值，需调用方（MCP tool 层或 LLM）填充
6. **无跨搜索引擎参数适配**：当前参数格式通用，未针对 pFind/MSFragger/Comet 的特有参数做差异化推荐

**后续改进方向**：
- 利用 charge 分布判断标记类型（SILAC 数据 charge 分布偏向高电荷态）
- 从 isolation window 宽度推断 DDA vs DIA，DIA 可推荐更宽容差
- 引入简单统计模型替代硬编码阈值
- LLM 层可在规则引擎输出基础上做进一步润色和用户交互

---

### Milestone 1.4：`search-engine` — 搜索引擎调度 Library Crate

> **状态**：✅ 已完成（37 tests, 0 warnings）— 采用 SimpleSearchEngine + pFind 预留结构
>
> 实现方案调整：先实现简化内置搜索引擎验证端到端数据流，pFind adapter 保留桩待后续对接。
> 关联 FR：FR-3.1 ~ FR-3.8 | 关联 US：US-2, US-5

#### Task 1.4.1：创建 crate 并定义模块结构 ✅

- **Sub-task 1.4.1.1**：创建 `crates/search-engine/Cargo.toml`，依赖 `core`, `spectrum-io`, `tokio`, `async-trait`
- **Sub-task 1.4.1.2**：EngineRegistry（register/get/list_available/health_check_all）
- **Sub-task 1.4.1.3**：SearchProgress 结构体（run_id, status, progress_pct, elapsed_sec）
- **Sub-task 1.4.1.4**：SearchEngineError（6 变体）+ From for CoreError
- **实现补充**：将 core 的 SearchEngineAdapter trait 从 `async fn in trait` 改为 `#[async_trait]` 以支持 dyn 兼容

#### Task 1.4.2：FASTA 解析 + 酶切消化 ✅

- **Sub-task 1.4.2.1**：`fasta.rs`：解析标准 FASTA 文件（>header + 多行 sequence）
- **Sub-task 1.4.2.2**：`digest.rs`：支持全部 7 种 Enzyme 变体的酶切消化规则
  - Trypsin（K/R 后切，P 除外）、TrypsinP、LysC、GluC、AspN、Chymotrypsin、NonSpecific
- **Sub-task 1.4.2.3**：missed cleavages 支持（0-N）、肽段长度过滤（6-50 aa）
- **Sub-task 1.4.2.4**：`chemistry.rs`：20 种标准氨基酸单同位素质量表 + peptide_mass() + peptide_mz()

#### Task 1.4.3：谱图匹配 + 打分 ✅

- **Sub-task 1.4.3.1**：`matching.rs`：precursor m/z 匹配（ppm/Da 容差）
- **Sub-task 1.4.3.2**：理论 b/y 离子生成（单电荷）
- **Sub-task 1.4.3.3**：碎片离子匹配（binary search 在排序 peak list 中）
- **Sub-task 1.4.3.4**：打分 = matched_ions / total_ions
- **Sub-task 1.4.3.5**：固定修饰质量调整、多电荷态尝试（未知 charge 时试 2→3→1→4）

#### Task 1.4.4：SimpleSearchEngine 组装 ✅

- **Sub-task 1.4.4.1**：`simple_engine.rs` 实现 SearchEngineAdapter trait
- **Sub-task 1.4.4.2**：完整流程：validate → read FASTA → digest → read spectra → match → score → aggregate → SearchResult
- **Sub-task 1.4.4.3**：肽段级聚合（best score, PSM count）、蛋白级聚合（位置追踪 coverage）
- **Sub-task 1.4.4.4**：SearchResultSummary 统计（中位 score、delta ppm、charge/mod 分布、鉴定率）
- **Sub-task 1.4.4.5**：RunMetadata 自动追踪（run_id、duration、Completed 状态）

#### Task 1.4.5：pFind adapter 预留 ✅

- **Sub-task 1.4.5.1**：`adapters/pfind.rs`：PFindAdapter + SshConfig 结构体
- **Sub-task 1.4.5.2**：search() 返回 "not yet implemented" 错误，health_check() 返回 Unavailable
- **Sub-task 1.4.5.3**：实现路线图文档化在模块 doc 中

#### M1.4 已知局限性

1. **SimpleSearchEngine 是验证引擎**：搜索质量不如 pFind/MSFragger，仅用于验证架构和数据流
2. **O(N×M) 全量匹配**：每张谱图对每个候选肽段逐一比较，无索引优化。1 spectrum × 20420 proteins ≈ 0.83 sec，2000 spectra 会需要 ~30 min
3. **无统计学打分**：用 matched_b_y_ions / total_theoretical_ions 比例打分，非 hyperscore/binomial 统计模型
4. **无 FDR 计算**：所有 PSM 默认通过（q_value = None），无 target-decoy 策略执行
5. **电荷范围**：已知 charge 用实际值；未知 charge 尝试 2→3→1→4 四个电荷态
6. **仅支持单电荷碎片离子**：b/y 离子只生成 z=1，不生成 z=2 或中性丢失离子
7. **pFind 待对接**：需提供 .cfg 格式样例和输出文件样例

**真实数据验证**：使用 pBuild MGF (1 spectrum, charge 5+) + UniProt Human reviewed (20420 proteins) 完成端到端搜索，鉴定到 `DQSSWQNSDASQEVGGHQER` (TB182_HUMAN, Δ5.6 ppm, score 0.079)。

#### 原 M1.4 任务 1.4.3-1.4.7（pFind 相关）→ 推迟

原计划的 pFind 配置生成、SSH 远程执行、结果解析、进度跟踪、集成测试推迟到获取 pFind 样例文件后实施。对应的 Task 编号保留，实现时参考原 Sub-task 描述。

---

### Milestone 1.5：`report` — 报告生成 Library Crate

> **状态**：✅ 已完成（10 tests, 0 warnings）
>
> 搜索结果的统计聚合与导出。纯计算，不依赖 MCP。
> 关联 FR：FR-4.1 ~ FR-4.6 | 关联 US：US-6

#### Task 1.5.1：创建 crate ✅

- **Sub-task 1.5.1.1**：创建 `crates/report/Cargo.toml`，依赖 `core`, `serde_json`, `thiserror`
- **实现补充**：ReportError（3 变体）+ From for CoreError

#### Task 1.5.2：实现摘要生成 ✅

- **Sub-task 1.5.2.1**：实现 `ReportGenerator::generate_summary(result) -> SearchResultSummary`
- **Sub-task 1.5.2.2**：统计计算：total_psms, 按 q_value ≤ 0.01 过滤得 psms_at_1pct_fdr, 去重得 unique_peptides, 蛋白分组
- **Sub-task 1.5.2.3**：计算 identification_rate = psms_at_1pct_fdr / total_spectra_searched
- **Sub-task 1.5.2.4**：统计 modification_distribution（每种修饰出现次数）和 charge_distribution
- **Sub-task 1.5.2.5**：使用 `core::util::compute_median` 计算 median_score 和 median_delta_mass_ppm（跨 crate 统一算法）
- **实现补充**：无 q-value 时保留所有 PSM（兼容 SimpleSearchEngine）。`SearchResult.summary` 是引擎侧初步统计，`generate_summary()` 是带 FDR 过滤的最终版本。

#### Task 1.5.3：实现结果导出 ✅

- **Sub-task 1.5.3.1**：实现 `ReportGenerator::export_tsv(result, output_dir)` → 生成 3 个 TSV：psm.tsv, peptide.tsv, protein.tsv
- **Sub-task 1.5.3.2**：实现 `ReportGenerator::export_json(result, output_path)` → 完整 JSON 导出
- **Sub-task 1.5.3.3**：实现 `ReportGenerator::export_metadata(metadata, output_path)` → run_metadata.json
- **实现补充**：TSV 输出使用 `sanitize_tsv()` 转义 tab/换行字符。

#### Task 1.5.4：实现结果比较（P2，推迟）

- 推迟到后续迭代。ComparisonSummary 结构和 compare() 方法未实现。

#### Task 1.5.5：测试 ✅

- **Sub-task 1.5.5.1**：构造 SearchResult fixture → 验证 summary 各字段计算正确（FDR 过滤、median、charge 分布）
- **Sub-task 1.5.5.2**：验证导出的 TSV 文件创建和字段内容正确
- **Sub-task 1.5.5.3**：验证 JSON roundtrip 和 run_metadata.json 字段完整

---

### Milestone 1.6：`mcp-server` — MCP Server 二进制 + Agent/Skill

> **状态**：✅ 已完成（8 MCP tools 注册，Agent + 2 Skill 编写完成）
>
> 组装所有 library crate 为 MCP Server，编写蛋白搜索 Agent 和 Skill。
> 关联 FR：FR-5.1 ~ FR-5.5 | 关联 US：US-1

#### Task 1.6.0：升级 schemars + rust-version ✅

- schemars 0.8 → 1.x（rmcp v1.3 依赖）
- rust-version 1.80 → 1.85

#### Task 1.6.1：搭建 MCP Server 框架 ✅

- **Sub-task 1.6.1.1**：创建 `crates/mcp-server/Cargo.toml`，依赖所有 library + `rmcp` v1.3, `tokio`, `tracing`
- **Sub-task 1.6.1.2**：实现 `main.rs`：初始化 tracing → 构建 ProteinCopilotServer → `serve(stdio())`
- **Sub-task 1.6.1.3**：定义 `AppConfig` 结构体（ssh_config, pfind_config, data_dirs），从配置文件或环境变量加载 → **推迟到 pFind 接入时实现**（当前 SimpleSearchEngine 无需外部配置）
- **Sub-task 1.6.1.4**：验证 MCP Server 启动 → `Copilot CLI` 能发现并列出 tools

#### Task 1.6.2：注册 spectrum-io Tools ✅

- **Sub-task 1.6.2.1**：实现 `#[tool] read_spectra(file_path, format?)` → 调用 `spectrum_io::create_reader` → `reader.read_summary()` → 返回 JSON
- **Sub-task 1.6.2.2**：实现 `#[tool] get_spectrum(file_path, scan_number)` → `reader.read_spectrum()` → 返回 JSON
- **Sub-task 1.6.2.3**：错误转换：`CoreError` → MCP error response（error_code + message + suggestion）

#### Task 1.6.3：注册 param-recommend Tools ✅

- **Sub-task 1.6.3.1**：实现 `#[tool] recommend_params(spectrum_summary, user_hints?)` → `ParamRecommender::recommend()` → 返回 JSON
- **Sub-task 1.6.3.2**：实现 `#[tool] list_presets()` → `ParamRecommender::list_presets()` → 返回 JSON

#### Task 1.6.4：注册 search-engine Tools ✅

- **Sub-task 1.6.4.1**：实现 `#[tool] run_search(params, input_files, engine?)` → `EngineRegistry::get(engine)` → `adapter.search()` → 返回 JSON
- **Sub-task 1.6.4.2**：实现 `#[tool] check_engine()` → 遍历所有 adapter → 返回 EngineInfo + HealthStatus 列表
- **Sub-task 1.6.4.3**：实现 `#[tool] get_search_status(run_id)` → 返回 SearchProgress → **推迟到 pFind 接入时实现**（SimpleSearchEngine 同步执行，无中间进度可查询）

#### Task 1.6.5：注册 report Tools ✅

- **Sub-task 1.6.5.1**：实现 `#[tool] generate_summary(search_result)` → `ReportGenerator::generate_summary()` → 返回 JSON
- **Sub-task 1.6.5.2**：实现 `#[tool] export_results(search_result, format, output_path)` → 调用导出方法 → 返回文件路径

#### Task 1.6.6：编写 proteomics-search.agent.md ✅

- **Sub-task 1.6.6.1**：frontmatter：description, tools 列表
- **Sub-task 1.6.6.2**：角色定义：蛋白质质谱搜索助手，覆盖新手和资深用户
- **Sub-task 1.6.6.3**：标准工作流程指令：
  1. 询问用户谱图文件路径和 FASTA 路径
  2. 调用 `read_spectra` 获取数据摘要，向用户展示
  3. 调用 `recommend_params` 获取推荐参数
  4. 用自然语言向用户解释推荐理由（基于 AiDecision.explanation）
  5. **等待用户确认或修改参数**
  6. 调用 `run_search` 执行搜索
  7. 搜索完成后调用 `generate_summary` 生成摘要
  8. 基于 SearchResultSummary 生成自然语言分析报告
- **Sub-task 1.6.6.4**：决策边界规则：
  - 必须确认：搜索参数、搜索引擎选择
  - 可自动执行：读取谱图、生成摘要
  - 禁止：自行估算数值、跳过参数确认
- **Sub-task 1.6.6.5**：嵌入领域知识片段：常见消化酶说明、FDR 含义、常见修饰类型

#### Task 1.6.7：编写 basic-search.prompt.md ✅

- **Sub-task 1.6.7.1**：mode: agent, description
- **Sub-task 1.6.7.2**：step-by-step 流程（与 agent 的标准工作流程一致，但更紧凑）
- **Sub-task 1.6.7.3**：输入要求：谱图文件路径 + FASTA 路径
- **Sub-task 1.6.7.4**：输出格式：搜索结果摘要表 + 自然语言解读段落

#### Task 1.6.8：编写 result-interpretation.prompt.md ✅

- **Sub-task 1.6.8.1**：mode: agent, description
- **Sub-task 1.6.8.2**：引导 LLM 从 SearchResultSummary 中提取关键洞察的模板：
  - 鉴定率是否在正常范围（标准搜索 20-40%，磷酸化 5-15%）
  - FDR 分布是否健康
  - 修饰分布是否符合实验预期
  - 质量偏差中位数是否在 tolerance 范围内
- **Sub-task 1.6.8.3**：常见问题的解释模式库（鉴定率低的可能原因、异常修饰分布的解释等）

#### Task 1.6.9：配置 MCP Server 注册 ✅

- **Sub-task 1.6.9.1**：创建项目根目录 `.mcp.json`：注册 mcp-server 二进制（stdio transport）
- **Sub-task 1.6.9.2**：验证 Copilot CLI 打开项目后能自动发现 server 并列出所有 tools

---

### Milestone 1.7：端到端集成与验证

> 验证完整流程，编写文档。
> 关联 SM：SM-1, SM-7

#### Task 1.7.1：端到端测试场景

- **Sub-task 1.7.1.1**：准备测试数据：小型 mgf（100 张谱图）+ 小型 FASTA（100 蛋白）
- **Sub-task 1.7.1.2**：场景 A（自动化流程）：用户说"搜索这个数据" → Agent 执行完整流程 → 验证输出报告
- **Sub-task 1.7.1.3**：场景 B（指定参数）：用户直接给参数 → 跳过推荐 → 验证搜索执行
- **Sub-task 1.7.1.4**：场景 C（探查数据）：用户只想看谱图摘要 → 验证 read_spectra 返回正确
- **Sub-task 1.7.1.5**：记录所有问题和优化点到 GitHub Issues

#### Task 1.7.2：文档

- **Sub-task 1.7.2.1**：更新 `README.md`：项目简介、快速开始（安装 → 配置 SSH → 运行）、架构图
- **Sub-task 1.7.2.2**：编写 `docs/mcp-tools.md`：每个 MCP Tool 的名称、输入 schema、输出 schema、使用场景、示例
- **Sub-task 1.7.2.3**：编写 `docs/development.md`：如何添加新搜索引擎 adapter、如何添加新 MCP Tool、测试指南
- **Sub-task 1.7.2.4**：编写 `docs/configuration.md`：SSH 配置、pFind 配置、数据目录配置

---

## Phase 2：增强能力（大纲，后续细化）

### Milestone 2.1：`mcp-qc` — 质控模块
- 谱图质量评估（信噪比、前体离子纯度、保留时间分布异常检测）
- QC 报告生成
- MCP Tool：`run_qc`, `get_qc_report`

### Milestone 2.2：`mcp-fdr` — 独立 FDR 控制
- 实现 target-decoy FDR 计算（不再依赖 pFind 自带 FDR）
- 支持 Percolator / mokapot 重评分 adapter
- PSM / 肽段 / 蛋白质多级 FDR
- MCP Tool：`calculate_fdr`, `rescore`

### Milestone 2.3：`mcp-protein-inference` — 蛋白推断
- Parsimony 蛋白推断算法
- 蛋白分组与排序
- 推断过程可解释性输出
- MCP Tool：`infer_proteins`

### Milestone 2.4：多搜索引擎支持
- MSFragger adapter
- Comet adapter
- 多引擎结果合并与一致性分析
- MCP Tool：`compare_engines`

### Milestone 2.5：失败诊断 Agent
- `failure-diagnosis.agent.md`
- 搜索失败原因诊断（参数不合理、数据质量差、数据库不匹配等）
- 自动建议修复方案

---

## Phase 3：高级功能（大纲，后续细化）

### Milestone 3.1：多轮分析对话
- 支持在一次会话中多次迭代搜索（调参 → 重搜 → 比较）
- 分析历史上下文管理

### Milestone 3.2：开放搜索与未知修饰发现
- 支持 open search 模式
- 未知质量偏移的自动注释与解释

### Milestone 3.3：定量分析支持
- Label-free 定量
- TMT/iTRAQ 定量
- 定量结果的 AI 解读

### Milestone 3.4：可视化与交互式报告
- 谱图可视化输出（SVG/HTML）
- 交互式 Venn 图、火山图等

---

## 变更记录

| 日期 | 变更内容 | 原因 |
|---|---|---|
| 2026-03-26 | 初始任务规划版本 | 项目启动 |
| 2026-03-27 | 重写为 PRD 格式，补充 User Stories、Functional Requirements、Technical Considerations | 按 prd-creation.prompt.md 规范完善 |
| 2026-03-27 | 补充详细的 Milestone/Task/Sub-task 实施计划，与 architecture.md 对齐 | 实施计划细化 |

---

### Post-MVP：异步搜索模式

> **问题**：run_search 同步阻塞，长时间搜索（10+ min）导致 LLM 认为 MCP 超时。
> **方案**：改为异步模式，run_search 立即返回 run_id，后台执行搜索。

#### Task P1：实现异步搜索

- **Sub-task P1.1**：run_search 改为 `tokio::spawn` 后台执行，立即返回 `{run_id, status: "Running"}`
- **Sub-task P1.2**：添加 `progress_cache: HashMap<Uuid, SearchProgress>` 到 server
- **Sub-task P1.3**：搜索完成后将结果存入 `result_cache`，状态更新为 `Completed`
- **Sub-task P1.4**：实现 `get_search_status(run_id)` tool，返回 `SearchProgress`
- **Sub-task P1.5**：更新 Agent 指令：搜索后轮询 `get_search_status` 直到完成
- **Sub-task P1.6**：测试：验证异步搜索 + 状态轮询 + 结果获取的完整流程

#### 工作流变化

```text
之前（同步）：
  run_search() ── 10min 阻塞 ──→ SearchResult（LLM 可能超时）

之后（异步）：
  run_search()          → {run_id, status: "Running"}     （< 1秒）
  get_search_status()   → {progress: 30%, elapsed: 3min}  （轮询）
  get_search_status()   → {status: "Completed"}           （完成）
  generate_summary()    → SearchResultSummary              （结果）
```

---

## DIA 数据支持

### Phase 1: DIA 前体离子提取 ✅

**状态**: 完成（357 tests, 0 warnings）

新增 `dia-extraction` crate:
- DDA/DIA 自动检测（基于隔离窗口宽度中位数）
- MS1↔MS2 三级关联（spectrumRef → 扫描顺序 → RT 近邻）
- 同位素峰模式检测（¹³C-¹²C 间距，z=2-5）
- 伪谱图展开模式（兼容外部搜索引擎）

### Phase 2: 搜索引擎集成 ✅

**状态**: 完成

- `SearchParams` 增加 `acquisition_mode: Option<AcquisitionMode>`
- `match_spectrum_all()` 多母离子匹配
- 搜索引擎 MS2 过滤 + DIA/DDA 自动分流
- `extract_dia_precursors` MCP Tool（13号工具）
- `OrderedDiaCache` FIFO 缓存管理

### Phase 3: run_search 端到端集成 ✅

**状态**: 完成

- `SearchEngineAdapter` 增加 `search_with_spectra()` trait 方法
- `SimpleSearchEngine` 重构：`run_search_on_spectra()` 提取核心逻辑
- `RunSearchInput` 增加 `dia_run_id` 参数
- 完整 DIA 工作流：`extract_dia_precursors` → `run_search(dia_run_id=...)`
