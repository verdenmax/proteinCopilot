# ProteinCopilot — 项目总体任务规划

> **项目定位**：AI 驱动的蛋白质组学质谱搜索与结果解释平台
> **技术栈**：Rust Workspace + 多 MCP Server + Agent.md / Skill Prompt
> **当前阶段**：Phase 1 — MVP
> **维护方式**：本文件随开发进展动态更新

---

## Phase 1：MVP — 端到端基础搜索流程

**目标**：用户通过 Copilot CLI 打开项目，用自然语言描述搜索需求，系统完成"谱图导入 → 参数推荐 → pFind 搜索 → 结果汇总 → AI 报告"的完整流程。

---

### Milestone 1.1：项目骨架与 `core` crate

> 建立 Rust Workspace 结构，定义所有共享数据类型和 trait，为后续 MCP crate 提供基础。

#### Task 1.1.1：初始化 Rust Workspace

- **Sub-task 1.1.1.1**：创建根 `Cargo.toml`（workspace 定义），配置 `crates/*` 成员
- **Sub-task 1.1.1.2**：创建 `crates/core/Cargo.toml` 和 `crates/core/src/lib.rs`
- **Sub-task 1.1.1.3**：添加基础依赖：`serde`, `serde_json`, `thiserror`, `anyhow`, `tokio`, `tracing`, `uuid`
- **Sub-task 1.1.1.4**：配置 `rustfmt.toml` 和 `clippy.toml`，确保代码风格一致
- **Sub-task 1.1.1.5**：验证 `cargo build` 和 `cargo test` 通过

#### Task 1.1.2：定义谱图数据结构（`spectrum.rs`）

- **Sub-task 1.1.2.1**：定义 `Spectrum` 结构体（scan_number, precursor_mz, precursor_charge, mz_array, intensity_array, ms_level, retention_time 等）
- **Sub-task 1.1.2.2**：定义 `SpectrumSummary` 结构体（用于 MCP Tool 返回给 LLM 的摘要：谱图总数、质量范围、MS1/MS2 比例等）
- **Sub-task 1.1.2.3**：定义 `SpectrumFile` 枚举（MzML / Mgf），关联文件路径与元数据
- **Sub-task 1.1.2.4**：所有结构体 derive `Serialize`, `Deserialize`, `Debug`, `Clone`
- **Sub-task 1.1.2.5**：编写单元测试：序列化/反序列化 round-trip 测试

#### Task 1.1.3：定义搜索参数结构（`search_params.rs`）

- **Sub-task 1.1.3.1**：定义 `SearchParams` 结构体（enzyme, fixed_mods, variable_mods, precursor_tolerance, fragment_tolerance, missed_cleavages, database_path 等）
- **Sub-task 1.1.3.2**：定义 `Enzyme` 枚举（Trypsin, LysC, GluC, Chymotrypsin, 等）
- **Sub-task 1.1.3.3**：定义 `Modification` 结构体（name, mass_delta, residues, position）
- **Sub-task 1.1.3.4**：定义 `MassTolerance` 结构体（value, unit: Ppm / Da）
- **Sub-task 1.1.3.5**：实现 `SearchParams::validate()` 方法，检查参数合法性
- **Sub-task 1.1.3.6**：编写单元测试：合法/非法参数的验证测试

#### Task 1.1.4：定义搜索结果结构（`search_result.rs`）

- **Sub-task 1.1.4.1**：定义 `PSM`（Peptide-Spectrum Match）结构体（spectrum_ref, peptide_sequence, score, q_value, modifications, charge, precursor_mz, delta_mass 等）
- **Sub-task 1.1.4.2**：定义 `PeptideResult` 结构体（sequence, protein_accessions, best_score, q_value, psm_count）
- **Sub-task 1.1.4.3**：定义 `ProteinResult` 结构体（accession, description, coverage, peptide_count, unique_peptide_count）
- **Sub-task 1.1.4.4**：定义 `SearchResult` 结构体（run_id, engine_info, params_used, psms, peptides, proteins, summary_stats）
- **Sub-task 1.1.4.5**：定义 `SearchResultSummary` 结构体（供 LLM 解读：总 PSM 数、FDR 分布、匹配率等）
- **Sub-task 1.1.4.6**：编写单元测试

#### Task 1.1.5：定义 AI 决策输出结构（`ai_decision.rs`）

- **Sub-task 1.1.5.1**：定义 `AiDecision<T>` 泛型结构体（decision: T, confidence: f64, explanation: String, input_summary: String, alternatives: Vec<String>, evidence: Vec<String>）
- **Sub-task 1.1.5.2**：编写单元测试：确保 JSON 序列化格式符合 copilot-instructions 中的规范

#### Task 1.1.6：定义领域错误类型（`error.rs`）

- **Sub-task 1.1.6.1**：定义 `CoreError` 枚举（SpectrumParseError, InvalidSearchParams, SearchEngineError, FileNotFound 等），使用 `thiserror`
- **Sub-task 1.1.6.2**：每个错误变体包含上下文信息和建议操作

#### Task 1.1.7：定义搜索引擎 Adapter trait

- **Sub-task 1.1.7.1**：定义 `SearchEngineAdapter` trait（search, engine_info, health_check）
- **Sub-task 1.1.7.2**：定义 `EngineInfo` 结构体（name, version, supported_features）
- **Sub-task 1.1.7.3**：定义 `HealthStatus` 枚举（Healthy, Degraded, Unavailable）

#### Task 1.1.8：定义运行元数据结构（`run_metadata.rs`）

- **Sub-task 1.1.8.1**：定义 `RunMetadata` 结构体（run_id: Uuid, created_at, input_files, params_used, engine_info, duration, status）
- **Sub-task 1.1.8.2**：实现 `RunMetadata::new()` 自动生成 run_id 和 timestamp

---

### Milestone 1.2：`mcp-spectrum-io` — 谱图读取 MCP Server

> 实现第一个 MCP Server，能读取 mzML/mgf 文件并返回谱图数据和摘要。

#### Task 1.2.1：搭建 MCP Server 框架

- **Sub-task 1.2.1.1**：创建 `crates/mcp-spectrum-io/Cargo.toml`，引入 MCP SDK 依赖（调研选择 `rmcp` 或手动实现 JSON-RPC）
- **Sub-task 1.2.1.2**：实现 MCP Server 入口（stdio transport），注册 tools
- **Sub-task 1.2.1.3**：验证 MCP Server 能被 Copilot CLI / Claude Desktop 发现和连接

#### Task 1.2.2：实现 mgf 格式解析

- **Sub-task 1.2.2.1**：实现 `MgfParser`：逐条读取 BEGIN IONS...END IONS 块
- **Sub-task 1.2.2.2**：解析 TITLE, PEPMASS, CHARGE, RTINSECONDS 等字段
- **Sub-task 1.2.2.3**：解析 m/z + intensity 数组
- **Sub-task 1.2.2.4**：返回 `Vec<Spectrum>` 结果
- **Sub-task 1.2.2.5**：编写单元测试：使用小型 fixture mgf 文件测试

#### Task 1.2.3：实现 mzML 格式解析

- **Sub-task 1.2.3.1**：选择 XML 解析库（`quick-xml` 推荐）
- **Sub-task 1.2.3.2**：实现 `MzMLParser`：解析 `<spectrum>` 元素
- **Sub-task 1.2.3.3**：解析 binary data array（base64 decoded, zlib 解压）
- **Sub-task 1.2.3.4**：提取 scan number, precursor info, CV params
- **Sub-task 1.2.3.5**：编写单元测试：使用小型 fixture mzML 文件测试

#### Task 1.2.4：实现 MCP Tool `read_spectra`

- **Sub-task 1.2.4.1**：Tool 输入：`{ file_path: string, format?: "mzml" | "mgf" }`
- **Sub-task 1.2.4.2**：Tool 输出：`SpectrumSummary`（谱图数量、MS 级别分布、质量范围、保留时间范围等）
- **Sub-task 1.2.4.3**：自动检测文件格式（通过扩展名或文件头）
- **Sub-task 1.2.4.4**：错误处理：文件不存在、格式不支持、解析失败

#### Task 1.2.5：实现 MCP Tool `get_spectrum`

- **Sub-task 1.2.5.1**：Tool 输入：`{ file_path: string, scan_number: int }`
- **Sub-task 1.2.5.2**：Tool 输出：单张谱图的完整数据（`Spectrum` 序列化）
- **Sub-task 1.2.5.3**：支持按 scan number 精确提取

#### Task 1.2.6：集成测试

- **Sub-task 1.2.6.1**：准备 `tests/fixtures/` 下的测试 mgf 和 mzML 文件
- **Sub-task 1.2.6.2**：端到端测试：启动 MCP Server → 发送 JSON-RPC 请求 → 验证响应

---

### Milestone 1.3：`mcp-param-recommend` — 参数推荐 MCP Server

> 基于谱图特征生成搜索参数建议，供 LLM 向用户解释和调整。

#### Task 1.3.1：搭建 MCP Server 框架

- **Sub-task 1.3.1.1**：创建 crate，复用 Milestone 1.2 中建立的 MCP Server 框架模式
- **Sub-task 1.3.1.2**：注册 tools

#### Task 1.3.2：实现参数推荐规则引擎

- **Sub-task 1.3.2.1**：基于谱图摘要（MS 级别分布、质量范围、碎裂模式）推断仪器类型
- **Sub-task 1.3.2.2**：根据仪器类型推荐 precursor/fragment tolerance
- **Sub-task 1.3.2.3**：根据质量分布推荐消化酶
- **Sub-task 1.3.2.4**：提供常见修饰组合建议（标准、磷酸化、TMT 等预设方案）
- **Sub-task 1.3.2.5**：所有推荐结果包装为 `AiDecision<SearchParams>` 格式，包含 confidence 和 explanation

#### Task 1.3.3：实现 MCP Tool `recommend_params`

- **Sub-task 1.3.3.1**：Tool 输入：`{ spectrum_summary: SpectrumSummary, user_hints?: object }`
- **Sub-task 1.3.3.2**：Tool 输出：`AiDecision<SearchParams>`（推荐参数 + 推荐理由 + 置信度）
- **Sub-task 1.3.3.3**：`user_hints` 允许用户通过 LLM 传入模糊约束（如"磷酸化实验"、"高精度仪器"）

#### Task 1.3.4：实现 MCP Tool `list_presets`

- **Sub-task 1.3.4.1**：Tool 输出：可用的预设参数方案列表（标准蛋白搜索、磷酸化搜索、TMT 标记等）
- **Sub-task 1.3.4.2**：每个预设包含名称、描述、适用场景

#### Task 1.3.5：单元测试与集成测试

- **Sub-task 1.3.5.1**：测试不同谱图特征输入下的参数推荐正确性
- **Sub-task 1.3.5.2**：测试 user_hints 对推荐结果的影响

---

### Milestone 1.4：`mcp-search-engine` — 搜索引擎调度 MCP Server

> 通过 adapter 层调用 pFind 执行搜索，返回标准化结果。

#### Task 1.4.1：搭建 MCP Server 框架

- **Sub-task 1.4.1.1**：创建 crate，配置 adapter 模块结构

#### Task 1.4.2：实现 pFind Adapter

- **Sub-task 1.4.2.1**：实现 `PFindAdapter` 结构体，实现 `SearchEngineAdapter` trait
- **Sub-task 1.4.2.2**：将 `SearchParams` 转换为 pFind 配置文件格式（.cfg）
- **Sub-task 1.4.2.3**：调用 pFind 可执行文件（`Command::new()` 子进程）
- **Sub-task 1.4.2.4**：监控搜索进度（读取 pFind 的日志/进度输出）
- **Sub-task 1.4.2.5**：解析 pFind 输出结果文件，转换为 `SearchResult` 标准结构
- **Sub-task 1.4.2.6**：实现 `health_check`：检查 pFind 可执行文件是否存在且可运行

#### Task 1.4.3：实现 MCP Tool `run_search`

- **Sub-task 1.4.3.1**：Tool 输入：`{ params: SearchParams, input_files: [string], engine?: string }`
- **Sub-task 1.4.3.2**：Tool 输出：`SearchResult`（标准化搜索结果）
- **Sub-task 1.4.3.3**：支持异步执行（搜索可能耗时数分钟到数小时）
- **Sub-task 1.4.3.4**：生成 `RunMetadata`，记录本次搜索的完整上下文

#### Task 1.4.4：实现 MCP Tool `check_engine`

- **Sub-task 1.4.4.1**：Tool 输出：可用搜索引擎列表及其健康状态
- **Sub-task 1.4.4.2**：返回 `EngineInfo` + `HealthStatus`

#### Task 1.4.5：实现 MCP Tool `get_search_status`

- **Sub-task 1.4.5.1**：对于长时间运行的搜索，返回当前进度
- **Sub-task 1.4.5.2**：Tool 输出：`{ run_id, status, progress_pct, elapsed_time, estimated_remaining }`

#### Task 1.4.6：测试

- **Sub-task 1.4.6.1**：使用 mock adapter 测试 MCP Tool 逻辑
- **Sub-task 1.4.6.2**：pFind adapter 的集成测试（需要 pFind 安装环境）

---

### Milestone 1.5：`mcp-report` — 报告生成 MCP Server

> 将搜索结果转换为结构化摘要，供 LLM 解读和生成用户报告。

#### Task 1.5.1：搭建 MCP Server 框架

- **Sub-task 1.5.1.1**：创建 crate

#### Task 1.5.2：实现 MCP Tool `generate_summary`

- **Sub-task 1.5.2.1**：Tool 输入：`{ search_result: SearchResult }`
- **Sub-task 1.5.2.2**：Tool 输出：`SearchResultSummary`（总 PSM 数、1% FDR 下的鉴定数、蛋白数、修饰分布、质量偏差分布等）
- **Sub-task 1.5.2.3**：输出格式专为 LLM 解读优化（结构化 JSON + 人类可读字段名）

#### Task 1.5.3：实现 MCP Tool `export_results`

- **Sub-task 1.5.3.1**：Tool 输入：`{ search_result: SearchResult, format: "tsv" | "json", output_path: string }`
- **Sub-task 1.5.3.2**：导出标准化结果文件
- **Sub-task 1.5.3.3**：同时导出 `run_metadata.json`

#### Task 1.5.4：实现 MCP Tool `compare_results`（可选，MVP 简化版）

- **Sub-task 1.5.4.1**：比较两次搜索结果的差异（鉴定数变化、overlap 分析）
- **Sub-task 1.5.4.2**：输出比较摘要供 LLM 解释

#### Task 1.5.5：测试

- **Sub-task 1.5.5.1**：使用构造的 `SearchResult` 测试摘要生成正确性
- **Sub-task 1.5.5.2**：测试导出文件格式合规性

---

### Milestone 1.6：Agent 与 Skill 定义

> 编写蛋白质搜索领域的 Agent.md 和 Skill Prompt，使 LLM 能编排完整流程。

#### Task 1.6.1：编写 `proteomics-search.agent.md`

- **Sub-task 1.6.1.1**：定义 Agent 角色：蛋白质质谱搜索助手
- **Sub-task 1.6.1.2**：列出可用 MCP Tools（read_spectra, recommend_params, run_search, generate_summary 等）
- **Sub-task 1.6.1.3**：定义工作流程指令：先读谱图 → 推荐参数 → 用户确认 → 执行搜索 → 生成报告
- **Sub-task 1.6.1.4**：定义决策边界：何时自动执行 vs 何时请求用户确认
- **Sub-task 1.6.1.5**：嵌入领域知识：常见消化酶选择逻辑、修饰类型说明、FDR 阈值含义等

#### Task 1.6.2：编写 `basic-search.prompt.md` Skill

- **Sub-task 1.6.2.1**：定义基础搜索流程的 step-by-step prompt
- **Sub-task 1.6.2.2**：输入要求：谱图文件路径 + FASTA 数据库路径
- **Sub-task 1.6.2.3**：输出：搜索结果摘要 + AI 解读

#### Task 1.6.3：编写 `result-interpretation.prompt.md` Skill

- **Sub-task 1.6.3.1**：定义结果解释的 prompt 模板
- **Sub-task 1.6.3.2**：引导 LLM 从 SearchResultSummary 中提取关键洞察
- **Sub-task 1.6.3.3**：包含常见问题的解释模式（为什么鉴定率低？为什么修饰异常？）

#### Task 1.6.4：配置 MCP Server 注册

- **Sub-task 1.6.4.1**：在项目根目录创建 MCP 配置文件（如 `.mcp.json` 或 `mcp-config.json`），注册所有 MCP Server
- **Sub-task 1.6.4.2**：确保 Copilot CLI 能自动发现并启动各 MCP Server

---

### Milestone 1.7：端到端集成与验证

> 验证完整流程：用户通过 Copilot CLI → Agent 编排 → MCP 执行 → AI 报告。

#### Task 1.7.1：端到端测试场景

- **Sub-task 1.7.1.1**：准备测试数据集（小型 mgf + FASTA）
- **Sub-task 1.7.1.2**：手动测试：通过 Copilot CLI 发起搜索 → 验证完整流程
- **Sub-task 1.7.1.3**：记录流程中的问题和优化点

#### Task 1.7.2：文档

- **Sub-task 1.7.2.1**：编写 `README.md`：项目简介、快速开始、架构说明
- **Sub-task 1.7.2.2**：编写 `docs/mcp-tools.md`：所有 MCP Tool 的使用说明
- **Sub-task 1.7.2.3**：编写 `docs/development.md`：开发者指南（如何添加新 MCP Server / Adapter）

---

## Phase 2：增强能力（大纲，后续细化）

### Milestone 2.1：`mcp-qc` — 质控模块
- 谱图质量评估（信噪比、前体离子纯度、保留时间分布异常检测）
- QC 报告生成

### Milestone 2.2：`mcp-fdr` — 独立 FDR 控制
- 实现 target-decoy FDR 计算
- 支持 Percolator / mokapot 重评分 adapter
- PSM / 肽段 / 蛋白质多级 FDR

### Milestone 2.3：`mcp-protein-inference` — 蛋白推断
- Parsimony 蛋白推断算法
- 蛋白分组与排序
- 推断过程可解释性输出

### Milestone 2.4：多搜索引擎支持
- MSFragger adapter
- Comet adapter
- 多引擎结果合并与一致性分析

### Milestone 2.5：失败诊断 Agent
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

## 任务依赖关系

```text
M1.1 (core) ──┬──> M1.2 (spectrum-io)
              ├──> M1.3 (param-recommend)  ── 依赖 M1.2 的 SpectrumSummary
              ├──> M1.4 (search-engine)    ── 依赖 M1.2 + M1.3
              └──> M1.5 (report)           ── 依赖 M1.4
                                                │
M1.6 (agent/skill) ── 依赖 M1.2 ~ M1.5 的 Tool 定义
                                                │
M1.7 (集成验证)    ── 依赖所有 MVP Milestone
```

---

## 变更记录

| 日期 | 变更内容 | 原因 |
|---|---|---|
| 2026-03-26 | 初始版本创建 | 项目启动 |
