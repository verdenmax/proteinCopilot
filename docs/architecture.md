# ProteinCopilot — Rust Workspace 架构设计 v1

> 本文档从 Architect 视角，定义项目第一版 Rust workspace 的 crate 划分、职责边界、
> 共享 schema 设计、确定性逻辑与 LLM 编排的分层方案，以及关键架构决策。

---

## 1. 系统上下文

```text
┌──────────────────────────────────────────────────────────────────┐
│                        用户 (研究员)                              │
│  "帮我搜一下这批 HeLa 磷酸化数据"                                │
└──────────────────┬───────────────────────────────────────────────┘
                   │ 自然语言
                   ▼
┌──────────────────────────────────────────────────────────────────┐
│              MCP Client  (Copilot CLI / Claude Desktop)          │
│  ┌────────────────────────────────────────────────────────────┐  │
│  │ LLM (GPT-4o / Claude)                                     │  │
│  │  ├─ 读取 .github/agents/proteomics-search.agent.md        │  │
│  │  ├─ 读取 .github/prompts/basic-search.prompt.md           │  │
│  │  └─ 理解意图 → 规划步骤 → 调用 MCP Tools → 解释结果       │  │
│  └────────────────────────┬───────────────────────────────────┘  │
└───────────────────────────┼──────────────────────────────────────┘
                            │ MCP Protocol (JSON-RPC 2.0 / stdio)
                            ▼
┌──────────────────────────────────────────────────────────────────┐
│              ProteinCopilot MCP Server  (Rust)                   │
│                                                                  │
│  ┌─────────┐ ┌──────────────┐ ┌──────────────┐ ┌────────────┐  │
│  │spectrum │ │param-        │ │search-       │ │  report    │  │
│  │   -io   │ │ recommend    │ │  engine      │ │            │  │
│  │         │ │              │ │  ┌────────┐  │ │            │  │
│  │read_    │ │recommend_    │ │  │ pFind  │  │ │generate_   │  │
│  │spectra  │ │params        │ │  │adapter │  │ │summary     │  │
│  │         │ │              │ │  └────────┘  │ │            │  │
│  │get_     │ │list_         │ │run_search    │ │export_     │  │
│  │spectrum │ │presets       │ │check_engine  │ │results     │  │
│  └────┬────┘ └──────┬───────┘ └──────┬───────┘ └─────┬──────┘  │
│       │             │               │               │          │
│  ┌────┴─────────────┴───────────────┴───────────────┴──────┐   │
│  │                     core  (lib crate)                    │   │
│  │  Spectrum · SearchParams · SearchResult · AiDecision     │   │
│  │  SearchEngineAdapter trait · RunMetadata · CoreError      │   │
│  └──────────────────────────────────────────────────────────┘   │
└──────────────────────────────────────────────────────────────────┘
                            │
                            ▼ 子进程调用
                ┌───────────────────────┐
                │   搜索引擎 (pFind)     │
                └───────────────────────┘
```

---

## 2. 架构决策记录 (ADR)

### ADR-001: 单 MCP Server 二进制 + 多领域 Library Crate

**背景**：copilot-instructions.md 中原始设计为「每个 MCP Server 是独立 crate/二进制」，
即 spectrum-io、param-recommend 等各自是独立进程。

**决策**：MVP 阶段采用 **单 MCP Server 二进制 + 多 Library Crate** 架构。

**理由**：

| 方案 | 优点 | 缺点 |
|---|---|---|
| **A: 多 MCP Server 二进制** | 独立部署、独立扩展、故障隔离 | 多进程管理复杂；MCP 配置需注册多个 server；跨 server 数据传递需序列化；开发/调试成本高 |
| **B: 单二进制 + 多 Library（✅ 采用）** | 一个进程、一份 MCP 配置；crate 间可零成本共享内存中的数据；开发调试简单 | 无故障隔离；单进程资源限制 |

**后续演进路径**：由于领域逻辑在独立 library crate 中，未来任何 library 都可以被抽取为独立
MCP Server 二进制，只需新建一个 `main.rs` 注册其 tools。架构上预留了拆分能力。

**后果**：
- `crates/mcp-server/` 是唯一的二进制 crate，组合所有 library
- 每个领域模块是纯 library crate（无 MCP 依赖），可独立编译和测试
- `.mcp.json` 中只注册一个 server

---

### ADR-002: core crate 不依赖 MCP SDK

**决策**：`crates/core` 只包含领域数据结构和 trait，不引入 `rmcp` 依赖。

**理由**：
- core 是所有 crate 的依赖根，保持最小依赖面
- 领域类型应该是纯粹的数据定义，不绑定特定传输协议
- MCP 的 JSON Schema 通过 `schemars` derive 满足，不需要 rmcp 本身

**依赖**：core 仅依赖 `serde`, `serde_json`, `schemars`, `thiserror`, `uuid`, `chrono`

---

### ADR-003: 使用 rmcp 官方 SDK

**决策**：MCP 协议实现使用 `rmcp` crate（modelcontextprotocol/rust-sdk）。

**理由**：
- Anthropic 官方维护，协议合规性有保证
- 提供 `#[tool]` / `#[tool_router]` 宏，自动生成 JSON Schema 和路由
- 支持 stdio transport（Copilot CLI 使用此模式）
- 类型安全，与 serde + schemars 集成良好
- 活跃维护（当前 v0.16.0）

---

### ADR-004: 搜索引擎 Adapter 使用 async trait

**决策**：`SearchEngineAdapter` 是 `async trait`（Rust 1.75+ 原生支持 AFIT）。

**理由**：
- 搜索引擎调用涉及子进程管理、文件 I/O，天然异步
- 使用 `async fn in trait` 而非 `#[async_trait]` 宏，减少间接开销
- `Send + Sync` bound 确保可在 tokio 多线程运行时安全使用

---

## 3. Crate 划分与职责边界

```text
proteinCopilot/
├── Cargo.toml                         ← workspace
├── crates/
│   ├── core/                          ← [lib] 共享领域模型
│   ├── spectrum-io/                   ← [lib] 谱图解析（mgf, mzML）
│   ├── param-recommend/               ← [lib] 参数推荐规则引擎
│   ├── search-engine/                 ← [lib] 搜索引擎 adapter 层
│   ├── report/                        ← [lib] 结果摘要与导出
│   └── mcp-server/                    ← [bin] MCP Server（组装所有 tool）
```

### 3.1 `core`（lib crate）

**职责**：定义所有共享数据结构、领域 trait、错误类型。零业务逻辑。

**对外暴露**：
```text
core::
├── spectrum        ← Spectrum, SpectrumSummary, SpectrumFileInfo, MsLevel
├── search_params   ← SearchParams, Enzyme, Modification, MassTolerance, ToleranceUnit
├── search_result   ← PSM, PeptideResult, ProteinResult, SearchResult, SearchResultSummary
├── ai_decision     ← AiDecision<T>
├── engine          ← SearchEngineAdapter trait, EngineInfo, HealthStatus
├── run_metadata    ← RunMetadata, RunStatus
└── error           ← CoreError
```

**依赖**：`serde`, `serde_json`, `schemars`, `thiserror`, `uuid`, `chrono`
**不依赖**：`rmcp`, `tokio`, `tracing`（保持纯数据层）

**设计原则**：
- 所有 struct derive `Serialize`, `Deserialize`, `Debug`, `Clone`, `JsonSchema`
- `JsonSchema`（schemars）使得 MCP Tool 的 inputSchema 可从 Rust 类型自动生成
- 值对象优先使用 newtype pattern（如 `ScanNumber(u32)`）

---

### 3.2 `spectrum-io`（lib crate）

**职责**：谱图文件的解析和读取。纯确定性计算，无 MCP/网络依赖。

**对外暴露**：
```rust
pub trait SpectrumReader: Send + Sync {
    fn read_all(&self, path: &Path) -> Result<Vec<Spectrum>>;
    fn read_summary(&self, path: &Path) -> Result<SpectrumSummary>;
    fn read_spectrum(&self, path: &Path, scan: ScanNumber) -> Result<Spectrum>;
}

pub struct MgfReader;   // impl SpectrumReader
pub struct MzMLReader;  // impl SpectrumReader

pub fn detect_format(path: &Path) -> Result<SpectrumFileInfo>;
pub fn create_reader(info: &SpectrumFileInfo) -> Box<dyn SpectrumReader>;
```

**依赖**：`core`, `quick-xml`（mzML 解析）, `base64`（mzML binary data）, `flate2`（zlib 解压）

**设计原则**：
- Reader trait 使得未来增加新格式（mzXML, .raw）只需新增实现
- 大文件使用 streaming 解析（逐条读取，不一次性加载全部谱图到内存）
- 格式检测通过扩展名 + magic bytes 双重验证

---

### 3.3 `param-recommend`（lib crate）

**职责**：基于谱图特征推荐搜索参数。**确定性规则引擎**，不调用 LLM。

**对外暴露**：
```rust
pub struct ParamRecommender;

impl ParamRecommender {
    pub fn recommend(&self, summary: &SpectrumSummary, hints: Option<&UserHints>)
        -> AiDecision<SearchParams>;
    pub fn list_presets() -> Vec<SearchPreset>;
}

pub struct UserHints {
    pub experiment_type: Option<String>,   // "phosphorylation", "TMT"
    pub instrument_type: Option<String>,   // "Orbitrap", "TOF"
    pub custom_notes: Option<String>,
}

pub struct SearchPreset {
    pub name: String,
    pub description: String,
    pub params: SearchParams,
    pub applicable_scenarios: Vec<String>,
}
```

**依赖**：`core`
**不依赖**：`rmcp`, `tokio`

**设计原则**：
- 推荐逻辑是纯函数：相同输入 → 相同输出
- `AiDecision<SearchParams>` 的 `explanation` 字段由规则引擎生成模板化文字
  （如「根据质量范围 [400-2000 Da]，推荐 precursor tolerance 10 ppm」）
- LLM 可以对这些模板化解释做进一步润色和扩展
- `UserHints` 是 LLM 将用户自然语言意图转译后传入的结构化提示

---

### 3.4 `search-engine`（lib crate）

**职责**：搜索引擎的调度、调用和结果解析。通过 Adapter 模式封装外部搜索引擎。

**对外暴露**：
```rust
// core 中定义的 trait
pub trait SearchEngineAdapter: Send + Sync {
    async fn search(&self, params: &SearchParams, input_files: &[PathBuf])
        -> Result<SearchResult>;
    fn engine_info(&self) -> EngineInfo;
    async fn health_check(&self) -> Result<HealthStatus>;
}

// search-engine crate 中的实现
pub struct PFindAdapter { config: PFindConfig }
impl SearchEngineAdapter for PFindAdapter { ... }

// 引擎注册与发现
pub struct EngineRegistry { engines: HashMap<String, Box<dyn SearchEngineAdapter>> }
impl EngineRegistry {
    pub fn register(&mut self, adapter: Box<dyn SearchEngineAdapter>);
    pub fn get(&self, name: &str) -> Option<&dyn SearchEngineAdapter>;
    pub fn list_available(&self) -> Vec<EngineInfo>;
}
```

**内部结构**：
```text
search-engine/src/
├── lib.rs
├── registry.rs              ← EngineRegistry
├── adapters/
│   ├── mod.rs
│   ├── pfind.rs             ← PFindAdapter
│   ├── pfind_config.rs      ← pFind 配置文件生成
│   └── pfind_parser.rs      ← pFind 输出结果解析
└── progress.rs              ← 搜索进度追踪
```

**依赖**：`core`, `tokio`（子进程管理）, `tracing`

**设计原则**：
- Adapter 内部逻辑完全隔离：pFind 的 .cfg 格式、输出格式解析等不泄露到外部
- `SearchResult` 是标准化输出——不管哪个引擎，返回相同结构
- 搜索执行是 async：通过 `tokio::process::Command` 管理子进程
- 进度追踪通过轮询日志文件或 stdout 实现

---

### 3.5 `report`（lib crate）

**职责**：将 `SearchResult` 转换为各种输出格式（摘要、TSV、JSON）。

**对外暴露**：
```rust
pub struct ReportGenerator;

impl ReportGenerator {
    pub fn generate_summary(result: &SearchResult) -> SearchResultSummary;
    pub fn export_tsv(result: &SearchResult, output: &Path) -> Result<()>;
    pub fn export_json(result: &SearchResult, output: &Path) -> Result<()>;
    pub fn compare(a: &SearchResult, b: &SearchResult) -> ComparisonSummary;
}
```

**依赖**：`core`, `csv`（TSV 导出）

---

### 3.6 `mcp-server`（bin crate）— 组装层

**职责**：唯一的二进制入口。组装所有 library，注册为 MCP Tools，启动 stdio server。

**内部结构**：
```text
mcp-server/src/
├── main.rs                  ← 入口：初始化 + 启动 MCP Server
├── tools/
│   ├── mod.rs
│   ├── spectrum.rs          ← read_spectra, get_spectrum tool 实现
│   ├── params.rs            ← recommend_params, list_presets tool 实现
│   ├── search.rs            ← run_search, check_engine, get_search_status tool 实现
│   └── report.rs            ← generate_summary, export_results tool 实现
└── config.rs                ← 服务器配置（搜索引擎路径等）
```

**模式**：每个 tool 模块是一个薄包装层——

```rust
// tools/spectrum.rs 示意
#[tool(description = "Read spectra from a mass spectrometry file (mzML/mgf). \
    Returns a summary including spectrum count, mass range, MS level distribution, \
    and retention time range. Use this as the first step to understand input data.")]
async fn read_spectra(&self, #[tool(param)] input: ReadSpectraInput) -> String {
    // 1. 调用 spectrum_io::create_reader()
    // 2. 调用 reader.read_summary()
    // 3. 序列化为 JSON 返回
}
```

**依赖**：`core`, `spectrum-io`, `param-recommend`, `search-engine`, `report`, `rmcp`, `tokio`, `tracing`

**设计原则**：
- Tool 实现只做「胶水」：参数解析 → 调用 library → 序列化响应
- 不在 tool 实现中写业务逻辑
- 错误统一转换为 MCP error response（包含 error_code + message + suggestion）

---

## 4. 共享 Schema 设计

### 4.1 核心类型关系图

```text
SpectrumFileInfo ──▶ SpectrumReader ──▶ Vec<Spectrum>
                                    └──▶ SpectrumSummary
                                              │
                                              ▼
UserHints ─────────▶ ParamRecommender ──▶ AiDecision<SearchParams>
                                              │
                                              ▼
SearchParams ──────▶ SearchEngineAdapter ──▶ SearchResult ──▶ ReportGenerator
  + input_files                │                                    │
                               ▼                                    ▼
                          RunMetadata                     SearchResultSummary
                                                          ComparisonSummary
```

### 4.2 关键类型定义

```rust
// ===== spectrum.rs =====

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub enum MsLevel { MS1, MS2, Other(u8) }

/// 隔离窗口（DDA 窄窗口, DIA 宽窗口），对齐 mzML <isolationWindow>
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct IsolationWindow {
    pub target_mz: f64,
    pub lower_offset: f64,  // Da
    pub upper_offset: f64,  // Da
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct PrecursorInfo {
    pub mz: f64,
    pub charge: Option<i32>,
    pub intensity: Option<f64>,
    pub isolation_window: Option<IsolationWindow>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct Spectrum {
    pub scan_number: u32,
    pub ms_level: MsLevel,
    pub retention_time_sec: f64,
    pub precursors: Vec<PrecursorInfo>,  // DDA: 1, DIA: 0~1(宽窗口), MS1: empty
    pub mz_array: Vec<f64>,
    pub intensity_array: Vec<f64>,
}

/// LLM 可读的谱图数据摘要——这是 AI 编排层了解数据特征的入口
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SpectrumSummary {
    pub file_path: String,
    pub format: SpectrumFormat,  // 改为枚举类型
    pub total_spectra: u64,
    pub ms1_count: u64,
    pub ms2_count: u64,
    pub mz_range: (f64, f64),
    pub rt_range_sec: (f64, f64),
    pub precursor_charge_distribution: HashMap<i32, u64>,
    pub median_peaks_per_spectrum: u32,
}

// ===== search_params.rs =====

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct SearchParams {
    pub enzyme: Enzyme,
    pub missed_cleavages: u32,
    pub fixed_modifications: Vec<Modification>,
    pub variable_modifications: Vec<Modification>,
    pub precursor_tolerance: MassTolerance,
    pub fragment_tolerance: MassTolerance,
    pub database_path: String,
    pub decoy_strategy: DecoyStrategy,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum Enzyme {
    Trypsin,
    LysC,
    GluC,
    AspN,
    Chymotrypsin,
    TrypsinP,            // Trypsin/P (不在 P 前切)
    NonSpecific,
    Custom { name: String, cleavage_rule: String },
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct Modification {
    pub name: String,           // e.g. "Carbamidomethyl", "Oxidation"
    pub mass_delta: f64,        // Da
    pub residues: Vec<char>,    // e.g. ['C'], ['M']
    pub position: ModPosition,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum ModPosition { Anywhere, AnyNTerm, AnyCTerm, ProteinNTerm, ProteinCTerm }

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct MassTolerance {
    pub value: f64,
    pub unit: ToleranceUnit,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum ToleranceUnit { Ppm, Da }

// ===== search_result.rs =====

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SearchResult {
    pub run_id: Uuid,
    pub engine_info: EngineInfo,
    pub params_used: SearchParams,
    pub psms: Vec<Psm>,
    pub peptides: Vec<PeptideResult>,
    pub proteins: Vec<ProteinResult>,
    pub summary: SearchResultSummary,
    pub metadata: RunMetadata,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct Psm {
    pub spectrum_scan: u32,
    pub peptide_sequence: String,
    pub modifications: Vec<Modification>,
    pub charge: i32,
    pub precursor_mz: f64,
    pub calculated_mz: f64,
    pub delta_mass_ppm: f64,
    pub score: f64,
    pub q_value: Option<f64>,
    pub protein_accessions: Vec<String>,
    pub is_decoy: bool,
}

/// LLM 可读的搜索结果统计摘要
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SearchResultSummary {
    pub total_spectra_searched: u64,
    pub total_psms: u64,
    pub psms_at_1pct_fdr: u64,
    pub unique_peptides_at_1pct_fdr: u64,
    pub protein_groups_at_1pct_fdr: u64,
    pub median_score: f64,
    pub median_delta_mass_ppm: f64,
    pub identification_rate: f64,       // psms_at_1pct_fdr / total_spectra
    pub modification_distribution: HashMap<String, u64>,
    pub charge_distribution: HashMap<i32, u64>,
    pub search_duration_sec: f64,
}

// ===== ai_decision.rs =====

/// 所有 AI 辅助决策的标准包装器
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct AiDecision<T> {
    pub decision: T,
    pub confidence: f64,             // 0.0 ~ 1.0
    pub explanation: String,         // 确定性规则引擎生成的解释模板
    pub input_summary: String,       // 决策依据的输入数据摘要
    pub alternatives: Vec<String>,   // 其他可选方案
    pub evidence: Vec<String>,       // 支持该决策的证据列表
}

// ===== engine.rs =====

pub trait SearchEngineAdapter: Send + Sync {
    async fn search(&self, params: &SearchParams, input_files: &[PathBuf])
        -> Result<SearchResult, CoreError>;
    fn engine_info(&self) -> EngineInfo;
    async fn health_check(&self) -> Result<HealthStatus, CoreError>;
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct EngineInfo {
    pub name: String,
    pub version: String,
    pub supported_features: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum HealthStatus {
    Healthy,
    Degraded { reason: String },
    Unavailable { reason: String },
}
```

---

## 5. 确定性逻辑与 LLM 编排的分层方案

### 5.1 职责边界矩阵

| 操作 | 由谁执行 | 位置 | 理由 |
|---|---|---|---|
| 谱图文件解析（mzML/mgf → Spectrum） | Rust | spectrum-io | 确定性 I/O |
| 谱图统计摘要（SpectrumSummary） | Rust | spectrum-io | 确定性计算 |
| 推断仪器类型 | Rust | param-recommend | 基于规则的特征匹配 |
| 推荐搜索参数 | Rust | param-recommend | 确定性规则引擎 |
| **解释推荐理由**（面向用户） | **LLM** | Agent.md | 自然语言润色 |
| **理解用户模糊指示** | **LLM** | Agent.md | 意图理解 |
| **将用户意图转化为 UserHints** | **LLM** | Agent.md | 语义映射 |
| 搜索参数验证 | Rust | core | 确定性校验 |
| 搜索引擎调用 | Rust | search-engine | 子进程管理 |
| pFind 配置文件生成 | Rust | search-engine/pfind | 模板化转换 |
| pFind 结果解析 → 标准化 | Rust | search-engine/pfind | 确定性解析 |
| FDR 计算（target-decoy） | Rust | core / mcp-fdr（Phase 2） | 数值计算 |
| 搜索结果统计摘要 | Rust | report | 确定性聚合 |
| **搜索结果解释**（面向用户） | **LLM** | Agent.md | 自然语言解释 |
| **搜索失败诊断** | **LLM** | Agent.md（Phase 2） | 推理诊断 |
| **规划多步分析流程** | **LLM** | Agent.md | 工作流编排 |
| 导出 TSV/JSON 结果文件 | Rust | report | 确定性 I/O |

### 5.2 数据流分层图

```text
Layer 0: 用户
  │  "搜一下这批磷酸化数据，用 10ppm 精度"
  ▼
Layer 1: LLM 编排层 (Agent.md)
  │  意图理解 → 提取关键信息：
  │    - 实验类型: "phosphorylation"
  │    - 精度要求: "10 ppm"
  │    - 操作: "搜索"
  │
  │  ① 调用 read_spectra(file_path)
  ▼
Layer 2: MCP Tool 层 (Rust)
  │  → spectrum-io 解析文件
  │  ← 返回 SpectrumSummary (JSON)
  │
Layer 1: LLM 编排层
  │  ② 构造 UserHints { experiment_type: "phosphorylation", ... }
  │     调用 recommend_params(spectrum_summary, user_hints)
  ▼
Layer 2: MCP Tool 层
  │  → param-recommend 规则引擎推荐参数
  │  ← 返回 AiDecision<SearchParams> (JSON)
  │
Layer 1: LLM 编排层
  │  ③ LLM 读取 AiDecision，用自然语言向用户解释推荐理由
  │     用户确认或调整参数
  │     调用 run_search(params, input_files)
  ▼
Layer 2: MCP Tool 层
  │  → search-engine 调用 pFind 子进程
  │  ← 返回 SearchResult (JSON)
  │
Layer 1: LLM 编排层
  │  ④ 调用 generate_summary(search_result)
  ▼
Layer 2: MCP Tool 层
  │  → report 生成 SearchResultSummary
  │  ← 返回 SearchResultSummary (JSON)
  │
Layer 1: LLM 编排层
  │  ⑤ LLM 基于 SearchResultSummary 生成用户报告：
  │     - 鉴定了 X 个肽段，Y 个蛋白
  │     - 磷酸化修饰主要出现在 S/T 位点
  │     - FDR 控制在 1%，鉴定率 35%，符合预期
  │     - 建议：可尝试放宽 fragment tolerance 以提高鉴定率
  ▼
Layer 0: 用户
  │  收到结构化报告 + 自然语言解读
```

### 5.3 LLM ↔ Rust 的交互契约

**规则 1：LLM 不计算，只传达**
```
✅ LLM: "根据 MCP Tool 返回的数据，1% FDR 下鉴定了 1,234 个 PSM"
❌ LLM: "我估计大约有 1,200 个 PSM"（自行估算）
```

**规则 2：Rust 不解释，只输出**
```
✅ Rust: SearchResultSummary { identification_rate: 0.35, ... }
❌ Rust: "鉴定率 35%，这个结果还不错"（不做主观判断）
```

**规则 3：LLM 向 Rust 传递结构化指令**
```
✅ LLM → recommend_params({ spectrum_summary: {...}, user_hints: { experiment_type: "phospho" } })
❌ LLM → recommend_params({ instructions: "用户想做磷酸化搜索，请推荐参数" })
```

---

## 6. 依赖图

```text
                    ┌──────────┐
                    │   core   │  (serde, schemars, thiserror, uuid, chrono)
                    └────┬─────┘
              ┌──────────┼──────────┬──────────┐
              ▼          ▼          ▼          ▼
        ┌───────────┐ ┌──────────┐ ┌────────┐ ┌────────┐
        │spectrum-io│ │  param-  │ │search- │ │ report │
        │           │ │recommend │ │ engine │ │        │
        │(quick-xml │ │          │ │(tokio) │ │ (csv)  │
        │ base64    │ │          │ │        │ │        │
        │ flate2)   │ │          │ │        │ │        │
        └─────┬─────┘ └────┬─────┘ └───┬────┘ └───┬────┘
              │            │           │           │
              └────────────┴───────┬───┴───────────┘
                                   ▼
                            ┌────────────┐
                            │ mcp-server │  (rmcp, tokio, tracing, clap)
                            │   [bin]    │
                            └────────────┘
```

依赖方向始终向上（library → core），mcp-server 在最底层聚合所有 library。
**禁止**：library 之间相互依赖（spectrum-io 不依赖 param-recommend）。

---

## 7. MCP Tool 注册表（MVP）

| Tool Name | 所属模块 | 输入 | 输出 | 场景 |
|---|---|---|---|---|
| `read_spectra` | spectrum-io | file_path, format? | SpectrumSummary | 了解数据特征 |
| `get_spectrum` | spectrum-io | file_path, scan_number | Spectrum | 查看单张谱图 |
| `recommend_params` | param-recommend | spectrum_summary, user_hints? | AiDecision\<SearchParams\> | 推荐搜索参数 |
| `list_presets` | param-recommend | (无) | Vec\<SearchPreset\> | 列出预设方案 |
| `run_search` | search-engine | params, input_files, engine? | SearchResult | 执行搜索 |
| `check_engine` | search-engine | (无) | Vec\<(EngineInfo, HealthStatus)\> | 检查引擎可用性 |
| `generate_summary` | report | search_result | SearchResultSummary | 生成结果摘要 |
| `export_results` | report | search_result, format, output_path | ExportResult | 导出结果文件 |

---

## 8. 可演进性设计

### 8.1 新增搜索引擎

```text
1. 在 search-engine/src/adapters/ 下新增 msfragger.rs
2. 实现 SearchEngineAdapter trait
3. 在 EngineRegistry 中注册
4. 无需修改 mcp-server 的 tool 定义（run_search 通过 engine? 参数分发）
```

### 8.2 新增 MCP 模块（如 mcp-qc）

```text
1. 创建 crates/qc/ (lib crate)，依赖 core
2. 在 mcp-server/src/tools/ 下新增 qc.rs，注册新 tool
3. 或者：创建独立 crates/mcp-qc/ (bin crate)，作为独立 MCP Server
```

### 8.3 拆分为多 MCP Server

```text
如果某个模块需要独立部署（如 search-engine 需要在 GPU 服务器上运行）：
1. 创建新 bin crate（如 crates/mcp-search/），仅引入 search-engine + core + rmcp
2. 在 .mcp.json 中注册为独立 server
3. 原有 mcp-server 移除该模块的 tool
```

---

## 9. 风险与缓解

| 风险 | 影响 | 缓解策略 |
|---|---|---|
| pFind 仅支持 Windows/特定 Linux | 限制可用平台 | health_check 提前检测；Adapter 层屏蔽平台差异 |
| 大谱图文件内存溢出 | 崩溃 | streaming 解析；read_spectra 只返回 summary 不加载全部 |
| 搜索耗时过长（小时级） | MCP 超时 | 异步执行 + get_search_status 轮询进度 |
| rmcp SDK 版本不稳定 | API break | Cargo.lock 锁版本；tool 层做薄包装 |
| LLM 产生幻觉参数 | 无效搜索 | SearchParams::validate() 在 Rust 侧强制校验 |

---

## 10. 下一步

1. 初始化 Workspace + core crate（Milestone 1.1 Task 1.1.1 ~ 1.1.8）
2. 验证 rmcp `#[tool_router]` 宏在最小示例中的可用性
3. 实现 spectrum-io 的 mgf 解析（最小可验证路径）
