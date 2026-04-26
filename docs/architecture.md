# ProteinCopilot — Rust Workspace 架构设计 v1

> 本文档从 Architect 视角，定义项目第一版 Rust workspace 的 crate 划分、职责边界、
> 共享 schema 设计、确定性逻辑与 LLM 编排的分层方案，以及关键架构决策。

> 演示版可视化页面见：`docs/architecture.html`
>
> `docs/architecture.md` 保留为完整的架构说明和决策记录。

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

### ADR-005: RT 单位约定 — `retention_time_min` 以分钟为单位

**决策**：`Spectrum.retention_time_min` 字段以 **分钟** 为单位（字段名 `_min` 后缀即分钟）。

**背景**：mzML 文件中 RT 可能以秒（`UO:0000010`）或分钟（`UO:0000031`）为单位。MGF 的 `RTINSECONDS` 以秒为单位。Sage 的 `scan_start_time` 也以分钟为单位。

**规则**：
- mzML parser：秒 → ÷60 转为分钟，分钟 → 直接使用
- MGF parser：RTINSECONDS → ÷60 转为分钟
- Sage adapter：直接传入 `retention_time_min`，不做二次转换
- spectrum-io 内部索引 `rt_seconds`：以秒存储，在边界处正确转换

**后果**：此前 Sage adapter 错误地将分钟再次除以 60（BUG-3），已修复。

---

### ADR-006: DIA 原始数据安全守卫

**决策**：DIA 数据在未经 `extract_dia_precursors` 提取前，不允许直接搜索。

**背景**：DIA 谱图的 `precursors[0].mz` 是隔离窗口中心值（如 500.0 Da），不是真实的前体离子 m/z。直接用于搜索会导致错误匹配。

**守卫层**：
1. `simple_engine.rs`：检测 isolation window 宽度 > 5 Da 的 DIA 谱图，自动跳过
2. `mcp-server/tools.rs`：`run_search` 检查输入文件的 `median_isolation_window_da`，DIA 文件必须提供 `dia_run_id`
3. `report/templates`：显示 `theoretical_mz` 而非 DIA 窗口中心值
4. MCP tool 层：通过 `find_precursor_in_ms1()` 从 MS1 谱图提取真实前体 m/z

---

### ADR-007: Decoy 检测统一化

**决策**：所有 decoy 蛋白检测通过 `core::util::is_decoy_accession()` 统一函数。

**支持前缀**：`REV_`、`SHUF_`、`DECOY_`、`REVERSED_`（覆盖主流搜索引擎的 decoy 命名约定）。

**使用点**：`simple_engine.rs`（PSM 标记）、`mapper.rs`（蛋白分组）、`parsimony.rs`（parsimony 推断）。

**后果**：此前硬编码 `"REV_"` 会漏检其他引擎的 decoy，已修复。

### ADR-008: HTML 安全与 Plotly.js 版本统一

**决策**：所有 HTML 模板遵循统一安全标准，Plotly.js 锁定 2.35.2。

**HTML 安全规范**：
- JSON 嵌入 `<script>` 时，`<` → `\u003c`、`>` → `\u003e`（gold standard：`escape_json_for_html()`）
- 用户数据写入 HTML 标签属性/内容时，`<>&"` → HTML 实体
- 所有 HTML 模板包含 `<meta name="viewport">` 以支持移动端

**Plotly.js 版本**：全项目统一 2.35.2 CDN（`mirror_plot.rs`、`multi_report.rs`、`entrapment_report.html`、`xic_visualize.rs`）

**MCP Server 命名常量**：DIA 阈值、RT 容差、FDR 阈值等魔法数字均提取为 `const`

**后果**：消除 XSS 注入风险，避免 Plotly 版本不一致导致渲染差异。

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
│   ├── dia-extraction/                ← [lib] DIA 前体离子提取
│   ├── fdr/                           ← [lib] FDR 计算（PSM/肽段/蛋白 三级）
│   ├── protein-inference/             ← [lib] 蛋白推断（parsimony + razor + 覆盖率）
│   ├── xic/                           ← [lib] XIC 提取与可视化（Plotly.js HTML）
│   ├── result-import/                 ← [lib] 外部搜索结果导入（DIA-NN / custom JSON）
│   ├── fasta-db/                      ← [lib] FASTA 数据库管理（注册表 + 下载 + 缓存）
│   ├── integration-tests/             ← [lib] 集成测试
│   ├── entrapment-analysis/           ← [lib] 陷阱库分析（L0-L4 分级 + edit distance + fragment provenance + 报告）
│   ├── entrapment-cli/                ← [bin] 陷阱库分析 CLI 工具
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
├── util            ← DECOY_PREFIXES, is_decoy_accession() 统一 decoy 检测
└── error           ← CoreError
```

**依赖**：`serde`, `serde_json`, `schemars`, `thiserror`, `uuid`, `chrono`
**不依赖**：`rmcp`, `tokio`, `tracing`（保持纯数据层）

**设计原则**：
- 所有 struct derive `Serialize`, `Deserialize`, `Debug`, `Clone`, `JsonSchema`
- `JsonSchema`（schemars）使得 MCP Tool 的 inputSchema 可从 Rust 类型自动生成

---

### 3.2 `spectrum-io`（lib crate）

**职责**：谱图文件的解析和读取。纯确定性计算，无 MCP/网络依赖。

**对外暴露**：
```rust
pub trait SpectrumReader: Send + Sync {
    fn read_all(&self, path: &Path) -> Result<Vec<Spectrum>, SpectrumIoError>;
    fn read_summary(&self, path: &Path) -> Result<SpectrumSummary, SpectrumIoError>;
    fn read_spectrum(&self, path: &Path, scan: u32) -> Result<Spectrum, SpectrumIoError>;
    fn for_each_spectrum(&self, path: &Path, handler: &mut dyn FnMut(Spectrum) -> Result<bool, SpectrumIoError>) -> Result<u32, SpectrumIoError>;
    fn list_ms2_meta(&self, path: &Path) -> Result<Vec<Ms2ScanMeta>, SpectrumIoError>;
    fn find_by_rt(&self, path: &Path, rt_min: f64, precursor_mz: f64, rt_tolerance_min: f64) -> Result<Option<(u32, f64)>, SpectrumIoError>;
}

pub struct Ms2ScanMeta { pub scan_number: u32, pub rt_min: f64, pub isolation_window: Option<(f64, f64, f64)> }
pub struct MgfReader;            // impl SpectrumReader
pub struct MzMLReader;           // impl SpectrumReader
pub struct IndexedMzMLReader;    // impl SpectrumReader（O(1) 随机访问 + O(log N) RT 查找）

pub fn detect_format(path: &Path) -> Result<SpectrumFileInfo, SpectrumIoError>;
pub fn create_reader(info: &SpectrumFileInfo) -> Box<dyn SpectrumReader>;
```

**子模块**：

| 模块 | 功能 |
|------|------|
| `index.rs` | ScanMeta + ScanIndex（RT 排序索引 + find_by_rt 二分查找 + 字节扫描元数据提取） |
| `disk_cache.rs` | PCIX v2 磁盘缓存（46B/entry 二进制格式，含 RT + ms_level + 隔离窗口） |
| `indexed_mzml.rs` | IndexedMzMLReader（两层策略：PCIX v2 缓存 → 字节扫描构建） |
| `mgf.rs` | MGF 格式解析器 |
| `mzml.rs` | mzML 流式解析器 |

**依赖**：`core`, `quick-xml`（mzML 解析）, `base64`（mzML binary data）, `flate2`（zlib 解压）, `memchr`（SIMD 字节扫描）

**设计原则**：
- Reader trait 使得未来增加新格式（mzXML, .raw）只需新增实现
- 大文件使用 streaming 解析（逐条读取，不一次性加载全部谱图到内存）
- IndexedMzMLReader 两层索引策略：PCIX v2 缓存 → 字节扫描构建（跳过 native index）
- `find_by_rt()` O(log N) 二分查找：RT 容差窗口 + 隔离窗口/前体 m/z 匹配
- `list_ms2_meta()` 从内存 ScanIndex 直接迭代，零 I/O
- 格式检测通过文件扩展名（大小写不敏感）
- 解析后自动按 m/z 升序排序（真实数据可能无序）
- mzML 保留时间自动单位转换（分钟→秒）
- 支持 DIA 隔离窗口（从 mzML `<isolationWindow>` 提取）

---

### 3.3 `param-recommend`（lib crate）

**职责**：基于谱图特征推荐搜索参数。**确定性规则引擎**，不调用 LLM。

**对外暴露**：
```rust
pub struct ParamRecommender;

impl ParamRecommender {
    pub fn recommend(&self, summary: &SpectrumSummary, hints: Option<&UserHints>)
        -> Result<AiDecision<SearchParams>, ParamRecommendError>;
    pub fn list_presets() -> Vec<SearchPreset>;
}

pub struct UserHints {
    pub experiment_type: Option<String>,   // "phosphorylation", "TMT"
    pub instrument_type: Option<String>,   // "Orbitrap", "TOF"
    pub enzyme: Option<Enzyme>,            // override default enzyme
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
  （如「根据 m/z 范围 [400-2000]，推荐 precursor tolerance 10 ppm」）
- LLM 可以对这些模板化解释做进一步润色和扩展
- `UserHints` 是 LLM 将用户自然语言意图转译后传入的结构化提示

---

### 3.4 `search-engine`（lib crate）

**职责**：搜索引擎的调度、调用和结果解析。包含一个简化的内置搜索引擎（MVP 验证用）、一个 Sage adapter（生产级搜索）和 pFind adapter 预留结构。

**对外暴露**：
```rust
// core 中定义的 trait（使用 #[async_trait]）
#[async_trait]
pub trait SearchEngineAdapter: Send + Sync {
    async fn search(&self, params: &SearchParams, input_files: &[PathBuf])
        -> Result<SearchResult, CoreError>;
    /// DIA 模式：直接接收已提取前体的谱图，跳过文件读取
    async fn search_with_spectra(&self, params: &SearchParams, spectra: Vec<Spectrum>)
        -> Result<SearchResult, CoreError>;
    fn engine_info(&self) -> EngineInfo;
    async fn health_check(&self) -> Result<HealthStatus, CoreError>;
}

// 简化内置搜索引擎（MVP）
pub struct SimpleSearchEngine;
impl SearchEngineAdapter for SimpleSearchEngine { ... }
// 完整流程: FASTA→酶切→precursor匹配→b/y离子打分→SearchResult
// 内部重构: run_search_on_spectra() 提取核心搜索逻辑，供 search() 和 search_with_spectra() 复用

// Sage adapter（生产级，集成 sage-core v0.15.0）
pub struct SageAdapter { thread_count: usize }
impl SearchEngineAdapter for SageAdapter { ... }
// 完整流程: FASTA→sage IndexedDatabase→rayon 并行打分→LDA rescoring→SearchResult
// 用户指定 engine: "Sage" 即可使用

// pFind adapter（预留桩，待对接真实 pFind）
pub struct PFindAdapter { ssh_config: SshConfig }
impl SearchEngineAdapter for PFindAdapter { ... } // 当前返回 not-implemented 错误

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
├── lib.rs                   ← 模块声明 + 公开 API
├── error.rs                 ← SearchEngineError
├── registry.rs              ← EngineRegistry
├── progress.rs              ← SearchProgress
├── chemistry.rs             ← 氨基酸质量表 + 质量计算（共享）
├── fasta.rs                 ← FASTA 数据库解析
├── digest.rs                ← 酶切消化（支持 7 种酶）
├── matching.rs              ← precursor m/z 匹配 + b/y 离子打分（已应用固定+可变修饰）
├── varmod.rs                ← 可变修饰站点发现 + 组合枚举（DFS）
├── annotate.rs              ← 谱图碎片离子注释（SpectrumAnnotation + HTML 可视化）
├── simple_engine.rs         ← SimpleSearchEngine 实现（含 decoy 生成 + FDR 集成）
└── adapters/
    ├── mod.rs
    ├── pfind.rs             ← PFindAdapter + SshConfig（预留桩）
    └── sage/                ← SageAdapter（sage-core 集成）
        ├── mod.rs           ← SageAdapter struct + SearchEngineAdapter impl
        ├── config.rs        ← SearchParams → sage_core::SageParameters 转换
        └── convert.rs       ← Spectrum/MassTolerance 类型转换
```

**依赖**：`core`, `spectrum-io`, `fdr`, `tokio`, `async-trait`, `sage-core`, `rayon`

**设计原则**：
- SimpleSearchEngine 是 MVP 验证引擎，用于测试端到端数据流正确性
- **SageAdapter** 集成 sage-core v0.15.0 作为库依赖（非 CLI 子进程），提供生产级搜索能力
- pFind adapter 预留完整结构（SshConfig、cfg 生成、结果解析），待提供 pFind 样例后对接
- Adapter 内部逻辑完全隔离：各引擎的配置格式、输出格式解析不泄露到外部
- `SearchResult` 是标准化输出——不管哪个引擎，返回相同结构
- 搜索执行是 async，通过 `#[async_trait]` 支持 `Box<dyn SearchEngineAdapter>`
- 氨基酸质量表集中在 `chemistry.rs`，避免重复

**SageAdapter 架构详情**：

Sage adapter 将 sage-core 作为 Rust 库依赖直接调用，避免 CLI 子进程的开销和复杂性。

```text
SearchParams + FASTA path
       │
       ▼
  config::build_sage_parameters()    ← 转换酶/修饰/tolerance 到 sage 类型
       │
       ▼
  Fasta::digest(&sage_params)        ← 酶切 + 构建 IndexedDatabase
       │
       ▼
  convert::spectrum_to_raw()         ← Spectrum → sage RawSpectrum（批量转换）
       │
       ▼
  SpectrumProcessor::process()       ← 谱图预处理（去噪/归一化）
       │
       ▼
  Scorer::score() [rayon parallel]   ← b/y 离子打分（spawn_blocking 桥接 tokio）
       │
       ▼
  LDA rescoring (discriminant_score) ← 线性判别分析重打分
       │
       ▼
  convert → Psm / SearchResult      ← 标准化输出
```

**关键特性**：
- **并行计算**：sage-core 使用 rayon 进行谱图级并行打分，通过 `tokio::task::spawn_blocking` 桥接异步运行时
- **FDR pipeline**：LDA rescoring（带 NaN 保护，poisson 值 clamp 到 -0.999）→ spectrum q-value → picked peptide FDR → picked protein FDR
- **支持特性**：b/y 离子打分、LDA rescoring、open search、LFQ、TMT、chimera scoring
- **酶支持**：全部 7 种标准酶（Trypsin, LysC, GluC, AspN, Chymotrypsin, NonSpecific, Custom）
- **修饰映射**：固定/可变修饰自动映射到 sage `ModificationSpecificity`（5 种位置：Residue, PeptideN/C, ProteinN/C）
- **结果扩展字段**：hyperscore, matched_peaks, delta_next, delta_best, discriminant_score, posterior_error, peptide_q, protein_q 等 sage 特有指标通过 `Psm.extra` 传递
- **蛋白 FDR**：使用 sage picked_protein 算法，protein_groups_at_1pct_fdr 基于 protein_q <= 0.01 过滤
- **引擎选择**：`run_search` 通过 `engine` 参数指定（大小写不敏感），默认 SimpleSearch

**SimpleSearchEngine 搜索算法详情**：

```text
输入: SearchParams + spectrum_files + FASTA

Step 1: 参数校验 (validate)
Step 2: 读取 FASTA → 蛋白质列表
Step 3: 酶切消化 → 候选肽段列表（target）
         ├── 按 Enzyme 规则切割 (Trypsin: K/R 后切, P 除外)
         ├── missed cleavages 0~N
         └── 过滤: 6 ≤ 肽段长度 ≤ 50 aa
Step 3b: Decoy 生成（如 DecoyStrategy != None）
         ├── Reverse: 反转蛋白序列（保留末尾 AA）
         ├── Shuffle: 随机打乱（种子 42）
         └── 酶切 decoy 蛋白 → decoy 肽段（accession 加 REV_/SHUF_ 前缀）
Step 4: 读取谱图文件 (spectrum-io)
Step 5: 逐谱图匹配（target + decoy 肽段一起参与竞争）
         对于每张 MS2 谱图:
           observed_mz = precursor m/z
           对于每个候选肽段:
             modified_mass = peptide_mass + Σ fixed_mod_delta
             对于每种可变修饰组合（DFS 枚举，受 max_variable_mods 限制）:
               total_mass = modified_mass + Σ varmod_delta
               对于每个 charge (已知用实际值, 未知试 2→3→1→4):
                 theoretical_mz = (total_mass + charge × 1.007276) / charge
                 if |observed - theoretical| / theoretical × 1e6 < tolerance_ppm:
                   合并固定+可变修饰 → 生成理论 b/y 离子
                   matched = 在实验 peak list 中 binary search 匹配数
                   score = matched / total_theoretical_ions
                   保留最高分匹配（记录修饰组合）
Step 6: FDR 计算（如有 decoy）
         ├── 竞争式 TDA: FDR(t) = decoys / targets 在每个分数阈值
         ├── q-value 单调化（反向扫描）
         └── 移除 decoy PSM
Step 7: 聚合 → PSM → Peptide → Protein (位置追踪 coverage)
Step 8: 统计 → SearchResultSummary + RunMetadata
```

**复杂度**: O(spectra × peptides × charges)，全量遍历无索引。
**性能实测**: 1 spectrum × 20420 proteins (UniProt Human) = 0.83 sec (release)。
**电荷范围**: 已知 charge 用实际值；未知 charge 尝试 2, 3, 1, 4（按常见频率排序）。
**打分模型**: 匹配碎片离子数 / 总理论碎片数（简化版，非统计学打分）。

---

### 3.5 `report`（lib crate）

**职责**：将 `SearchResult` 转换为统计摘要和各种输出格式（TSV、JSON）。

**对外暴露**：
```rust
pub struct ReportGenerator;

impl ReportGenerator {
    /// 生成带 FDR 过滤的统计摘要（q_value ≤ 0.01）。
    /// 与 SearchResult.summary（引擎侧初步统计）互补。
    pub fn generate_summary(result: &SearchResult) -> SearchResultSummary;

    /// 导出 3 个 TSV 文件：psm.tsv, peptide.tsv, protein.tsv
    pub fn export_tsv(result: &SearchResult, output_dir: &Path) -> Result<(), ReportError>;

    /// 导出完整 SearchResult 为 JSON
    pub fn export_json(result: &SearchResult, output_path: &Path) -> Result<(), ReportError>;

    /// 导出运行元数据为 JSON
    pub fn export_metadata(metadata: &RunMetadata, output_path: &Path) -> Result<(), ReportError>;
}
```

**内部结构**：
```text
report/src/
├── lib.rs                ← ReportGenerator 门面
├── error.rs              ← ReportError（3 变体）
├── summary.rs            ← FDR 过滤 + 统计聚合
├── export.rs             ← TSV/JSON 导出（含 sanitize_tsv 转义）
├── visualize.rs          ← 谱图注释 HTML 渲染（自包含 HTML）
├── unified_types.rs      ← UnifiedViewData / PeptideInfo 类型
└── unified_visualize.rs  ← 统一标注+XIC HTML 渲染
templates/
├── annotation.html       ← 纯标注 HTML 模板（SVG 谱图 + 覆盖图）
└── unified.html          ← 统一标注+XIC HTML 模板（含 SILAC 交互控件 + Plotly.js）
```

**对外暴露**（补充 `render_annotation` + `render_unified`）：
```rust
impl ReportGenerator {
    /// 渲染谱图注释为自包含 HTML 文件。
    pub fn render_annotation(
        annotation: &SpectrumAnnotation,
        output_path: &Path,
    ) -> Result<(), ReportError>;

    /// 渲染统一标注+XIC 视图为自包含 HTML 文件。
    /// 包含：文件名 + Scan/RT、覆盖图（SVG bracket）、谱图、
    /// SILAC 控件、逐离子 L/H 开关网格、MS1/MS2 XIC（Plotly.js）。
    /// 当 xic=None 时自动隐藏 XIC 相关区域（DDA 模式）。
    pub fn render_unified(
        data: &UnifiedViewData,
        output_path: &Path,
        plotly_mode: PlotlyMode,
    ) -> Result<(), ReportError>;
}
```

**依赖**：`core`（共享类型 + compute_median）, `search-engine`（SpectrumAnnotation 类型）, `xic`（IonMetadataEntry / RawScanData 类型）, `serde_json`

---

### 3.6 `dia-extraction`（lib crate）

**职责**：从 DIA（Data-Independent Acquisition）MS1 谱图中提取候选前体离子。通过同位素模式检测识别肽段前体，将 DIA 宽隔离窗口数据转换为类 DDA 格式供搜索引擎使用。

**核心组件**：

- **IsotopePatternExtractor**：检测 MS1 谱图中的同位素包络模式（isotope envelope），基于 averagine 模型估计电荷态和单同位素质量
- **MS1↔MS2 关联**：通过保留时间窗口将 MS1 提取的前体与对应 MS2 谱图关联
- **DDA/DIA 自动检测**：根据隔离窗口宽度自动判断采集模式（窄窗口 → DDA，宽窗口 → DIA）

**数据流**：
```text
mzML 输入
  → spectrum-io 读取所有谱图
  → 分离 MS1 / MS2
  → MS1 同位素模式提取 → 候选前体列表 (mz, charge, intensity)
  → 按保留时间窗口关联 MS1 前体 → MS2 谱图
  → 输出增强谱图：每张 MS2 携带多个候选前体
```

**集成方式**：
- 通过 `extract_dia_precursors` MCP Tool 暴露
- 提取结果写入 `OrderedDiaCache`（FIFO 缓存），返回 `run_id` 供后续 `run_search` 使用
- 搜索引擎的 `match_spectrum_all()` 支持多前体匹配，与 DIA 提取结果无缝衔接

**DIA 端到端工作流**：
```text
① extract_dia_precursors(file_path)
   → spectrum-io 读取 → DDA/DIA 自动检测
   → MS1 同位素模式提取 → MS1↔MS2 关联
   → 缓存增强谱图 → 返回 {dia_run_id, summary}

② run_search(input_files, database_path, dia_run_id=...)
   → 从 OrderedDiaCache 取出已提取的谱图
   → search_with_spectra(params, spectra)
   → run_search_on_spectra() 核心逻辑
   → SearchResult
```

**依赖**：`core`

---

### 3.7 `fdr`（lib crate）

**职责**：独立的 FDR 计算模块。生成 decoy 蛋白序列，实现竞争式 target-decoy 分析（TDA），计算 q-value 并强制单调化。

**核心模块**：
```text
fdr::
├── decoy       ← reverse_sequence()（保留末尾 AA）, shuffle_sequence()（种子 42）
├── calculation ← calculate_fdr()：竞争式 TDA, q-value 反向单调化
└── error       ← FdrError → CoreError 转换
```

**关键设计**：
- `reverse_sequence()` 保留最后一个氨基酸不动（tryptic decoy 兼容性）
- Shuffle 使用确定性种子（42）确保可复现
- FDR 公式：`FDR(t) = decoys_at_t / targets_at_t`（竞争式 TDA）
- q-value 单调化：从最差分数向最佳分数反向扫描，`q[i] = min(raw_fdr[i], q[i+1])`

**依赖**：`core`, `rand 0.8`（确定性 shuffle）

---

### 3.8 `mcp-server`（bin crate）— 组装层

**职责**：唯一的二进制入口。组装所有 library，注册为 MCP Tools，启动 stdio server。

**注册的 19 个 MCP Tools**：

| Tool | 功能 | 对应 Library |
|------|------|-------------|
| `read_spectra` | 读取谱图文件 → SpectrumSummary | spectrum-io |
| `get_spectrum` | 按 scan 读取单张谱图 | spectrum-io |
| `recommend_params` | 推荐搜索参数 → AiDecision | param-recommend |
| `list_presets` | 列出内置预设 | param-recommend |
| `run_search` | 异步执行搜索（立即返回 run_id，支持 `dia_run_id` 参数接收 DIA 提取结果） | search-engine |
| `get_search_status` | 查询搜索进度 | mcp-server (cache) |
| `cancel_search` | 取消正在运行的搜索 | mcp-server (cache) |
| `check_engine` | 检查引擎状态 | search-engine |
| `generate_summary` | FDR 过滤统计摘要 | report |
| `export_results` | 导出 TSV/JSON 文件 | report |
| `list_searches` | 列出历史搜索记录 | mcp-server (cache) |
| `annotate_spectrum` | 谱图碎片离子注释 | search-engine |
| `extract_dia_precursors` | DIA MS1 前体提取 | dia-extraction |
| `extract_spectrum_precursors` | 单张 MS2 母离子提取（调试用） | dia-extraction |
| `extract_xic` | XIC 碎片离子色谱图（支持 SILAC 轻重标记） | xic |
| `import_search_results` | 导入外部搜索结果（DIA-NN / custom JSON） | result-import |

**内部结构**：
```text
mcp-server/src/
├── main.rs       ← 入口：tracing 初始化 + ProteinCopilotServer + serve(stdio)
├── tools.rs      ← 19 个 tool 定义 + EngineRegistry 初始化
│                    使用 #[rmcp::tool_router] + #[rmcp::tool_handler] 宏
└── history.rs    ← 搜索历史持久化（磁盘 JSON）
```

**关键实现细节**：
- `#[rmcp::tool_router]` 自动生成 tool 注册和 JSON Schema
- `#[rmcp::tool_handler]` 自动实现 `list_tools` 和 `call_tool`
- `EngineRegistry` 在启动时注册 SimpleSearchEngine
- 错误通过 `mcp_core_err()` 统一转换，包含 `CoreError::suggestion()`
- `run_search` 入口显式调用 `params.validate()` 提前拦截无效参数
- `reader_cache: Arc<Mutex<HashMap<PathBuf, Arc<dyn IndexedSpectrumReader>>>>` 缓存已创建的索引读取器
- `RunCache = Arc<Mutex<OrderedRunCache>>` — FIFO 驱逐，最多 100 条

### 3.9 `xic`（lib crate）

**职责**：从 mzML 原始数据提取碎片离子的 XIC（Extracted Ion Chromatogram），支持 SILAC 重标记轻重离子对。同时提供客户端 SILAC 所需的 raw scan 数据和离子元数据。

**对外暴露**：
```text
xic::
├── extract::extract_xic()           ← 核心函数：肽段 + mzML → XicData
├── extract::extract_xic_with_raw()  ← 增强版：+ RawScanData + IonMetadataEntry
├── XicData / XicTrace / XicDataPoint
├── IonMetadataEntry                 ← 离子 K/R 计数（供客户端 SILAC 重计算）
├── RawScanData / RawScan            ← MS1/MS2 原始峰数组（嵌入 HTML 供 JS 使用）
├── ExtractionParams                 ← 提取参数（tolerance, n_cycles, top_n_ions 等）
├── LabelType / SilacLabel           ← SILAC 重标记定义
└── PlotlyMode                       ← Plotly.js CDN/Embedded 加载模式
```

**依赖**：`core`, `spectrum-io`
**不依赖**：`rmcp`（纯 library）

**设计要点**：
- MS2 fragment XIC 通过 `same_isolation_window()` 匹配同窗 MS2 扫描（DIA 多点、DDA 单点）
- MS1 precursor XIC 从全扫描中按 m/z 窗口提取（DDA/DIA 均可用）
- `extract_xic_with_raw()` 额外输出修剪后的 raw peaks 和 `compute_ion_metadata()` 的 K/R 计数，供 unified HTML 模板的客户端 SILAC 引擎使用
- MS1 修剪窗口根据肽段 K/R 数量和电荷态动态计算，确保 SILAC 重标前体不被截断
- DDA 无隔离窗口时，对 raw MS2 克隆使用 ±300s RT 近邻预过滤防止内存暴涨
- 入口处验证目标扫描必须为 MS2 级别

### 3.10 `result-import`（lib crate）

**职责**：将外部搜索引擎的结果（DIA-NN、自定义 JSON 等）导入为 ProteinCopilot 标准 `SearchResult`，并匹配 mzML 扫描号。

**对外暴露**：
```text
result-import::
├── detect_format()      ← 自动检测结果格式（parquet / JSON / pFind）
├── ResultParser trait    ← 统一解析接口
├── DiannParser          ← DIA-NN report.parquet 解析器
├── CustomJsonParser     ← hela.json 自定义格式解析器
├── PFindParser          ← pFind 结果解析（预留骨架）
├── ScanMatcher          ← RT + isolation window 匹配 mzML scan number（委托 reader.find_by_rt()）
├── Converter            ← ImportedPsm → core::Psm → SearchResult
├── UnimodDb             ← 22 内置修饰 + Unimod XML 解析器
├── ImportedPsm          ← 中间表示（format-agnostic）
├── ImportResult         ← 解析后结果（PSMs + metadata）
└── MatchReport          ← 匹配统计（per-file matched/unmatched/total MS2）
```

**依赖**：`core`, `spectrum-io`, `arrow`, `parquet`, `regex`, `serde_json`, `quick-xml`
**不依赖**：`rmcp`（纯 library）

**设计要点**：
- **格式检测**：`detect_format()` 基于文件扩展名 + magic bytes 自动识别
- **RT 单位转换**：外部数据 RT（分钟）在解析阶段自动 ×60 转为秒（项目内部统一秒）
- **扫描匹配**：`ScanMatcher` 委托 `reader.find_by_rt()` 进行 O(log N) RT 二分查找，DIA 模式额外检查 isolation window 包含性
- **score 方向**：DIA-NN Q.Value (lower=better) 存为 `1.0 - qvalue` 以满足 "higher=better" 约定
- **UnimodDb**：支持从名称/record_id 查找修饰质量偏移，解析 DIA-NN 带修饰序列（如 `M(UniMod:35)`）
- **下游兼容**：转换后的 SearchResult 与 run_search 产生的结果结构完全一致，可直接用于所有下游 tool
- 返回类型统一使用 `Result<Json<T>, ErrorData>`

**依赖**：`core`, `spectrum-io`, `param-recommend`, `search-engine`, `dia-extraction`, `report`, `protein-inference`, `fdr`, `fasta-db`, `rmcp` v1.3, `tokio`, `tracing`

### 3.11 `protein-inference`（lib crate）

**职责**：从 PSM 列表执行蛋白推断，产出最小蛋白组列表。

**子模块**：

| 模块 | 功能 |
|------|------|
| `mapper` | 构建肽段↔蛋白质双向映射图（I/L 等价、q-value 过滤、decoy 分类） |
| `parsimony` | 贪心集合覆盖：合并不可区分蛋白 → 移除子集 → 最小蛋白集 |
| `razor` | 共享肽段分配给证据最强蛋白组（按 unique 肽段数 > 分数 > 字母序） |
| `coverage` | 序列覆盖率计算（肽段定位到 FASTA 序列，I/L 归一化） |

**核心数据结构**（定义在 `core`）：
- `ProteinGroup`：leader + members + peptides (unique/shared/razor) + score + q_value + coverage
- `InferenceResult`：groups + razor_map + 统计信息

**设计决策**：
- 肽段的 I/L 等价处理：所有肽段序列 I→L 归一化后再比较
- Decoy 前缀固定为 `REV_`，与整个项目一致
- 蛋白组 leader 选择：字母序最先的 accession（确定性输出）
- 蛋白打分：组内最佳肽段分数

**依赖**：`core`, `thiserror`, `tracing`

### 3.12 `fasta-db`（lib crate）

**职责**：FASTA 蛋白数据库管理——内置物种库注册表、HTTPS 下载、本地文件缓存。

**功能**：
- 内置 UniProt 常用物种数据库注册表（Human, Mouse, E. coli, Yeast 等）
- HTTPS 下载 + 自动缓存到 `~/.proteincopilot/databases/`
- 缓存查询（已下载、大小、路径）

**依赖**：`reqwest`, `serde`, `tracing`

### 3.13 `entrapment-analysis`（lib crate）

**职责**：陷阱库命中分类与同源性分析。对搜索结果中的 trap PSM 进行 L0-L4 同源性分级，识别 razor 归属错误，生成 HTML 交互报告。

**核心模块**：

| 模块 | 功能 |
|------|------|
| `config.rs` | YAML 配置解析（SimilarityConfig + ProvenanceConfig + target/trap 规则） |
| `loader/` | 搜索结果加载（DIA-NN parquet + 通用 TSV） |
| `tagger.rs` | target/trap 标记（accession 规则匹配） |
| `digest.rs` | tryptic in-silico digest + k-mer 倒排索引（`find_similar()` 跨长搜索） |
| `similarity.rs` | L0-L4 分级（Phase A 等长 Hamming + Phase B 跨长 Levenshtein） |
| `levenshtein.rs` | Levenshtein edit distance + alignment（v2 新增） |
| `provenance.rs` | 碎片离子溯源引擎（b/y ion 匹配 + TrapOnly/TargetOnly/Shared 分类）（v3 新增） |
| `mod_parser.rs` | UniMod 修饰解析（`(UniMod:X)` → position + delta mass）（v3 新增） |
| `mirror_plot.rs` | trap vs target 镜像图渲染（Plotly.js HTML）（v3 新增） |
| `coelution.rs` | 共洗脱索引（RT 窗口交叉 + DIA 隔离窗口匹配）（v4 新增） |
| `multi_provenance.rs` | 多目标碎片匹配引擎（返回 MirrorData + shift_ions_heavy）（v4 新增） |
| `multi_report.rs` | 双扫描镜像 HTML 报告（轻/重标分离渲染 + 噪声过滤 + RT/q-value/修饰展示 + Heavy 计数 + 嵌合检测改进）（v4 新增） |
| `types.rs` | `ClassifiedPsm` + `SubstitutionType` + `FragmentProvenance` + `MirrorData` + `MultiTargetProvenance`（含 `trap_retention_time_min` + `trap_q_value`）+ `CoElutingCandidate` |
| `output.rs` | TSV 输出（含 substitution_type / edit_distance / provenance 列） |
| `report.rs` | HTML 交互报告（Plotly.js + mDa 显示 + 溯源统计） |

**v2 关键特性**：
- **Levenshtein edit distance**：替代 Hamming-only，支持 indel 跨长比较
- **k-mer 倒排索引**：pigeonhole 原理预筛，过滤 >99% 候选，加速跨长搜索
- **SubstitutionType 注释**：QKSubstitution、IsobaricDipeptideSingle/Dipeptide、NearIsobaric、Distinguishable、LengthMismatch
- **delta_mass 有符号化**：修复 v1 Hamming 路径中使用绝对值的问题
- **BestMatch tiebreaker**：使用 `<` 确保 delta_mass 比较可达

**v3 关键特性**：
- **碎片离子溯源**：trace_provenance() 对每个观测峰分类为 TrapOnly/TargetOnly/Shared/Unassigned
- **UniMod 修饰解析**：parse_modified_sequence() 解析 DIA-NN Modified.Sequence
- **RT-based scan lookup**：find_by_rt() 二分查找，DIA-NN 数据无 scan_number 时的回退方案
- **嵌合谱检测**：shared_ratio > chimera_threshold → is_chimeric 标记
- **容错 mzML 加载**：缺失文件跳过并 warn，不中断批量溯源

**v4 关键特性**：
- **共洗脱索引（CoElutionIndex）**：基于 RT 窗口交叉 + DIA 隔离窗口匹配，O(log N + k) 查询共洗脱 target
- **轻重标 SILAC 搜索**：同时查找 light 和 SILAC heavy 形式的候选（LabelForm::Light / Heavy）
- **多目标碎片匹配**：trace_multi_target() 返回 MirrorData，每个观测峰归属到 0..N 个 target
- **DIA 双扫描镜像**：轻/重标前体落在不同隔离窗口，通过 find_by_rt 分别定位各自的 MS2 扫描
- **MirrorData 结构**：`light: MirrorData` + `heavy: Option<MirrorData>` 独立存储各镜像数据
- **SILAC 偏移 trap 离子**：shift_ions_heavy() 为重标镜像生成偏移后的 trap 理论碎片离子
- **增强 HTML 报告**：双扫描号、候选前体 m/z、噪声过滤（<5%）、离子标注、无边框柱形
- **边界安全**：heavy_scan ≠ light_scan 验证、precursor_mz=None 守卫、parse_ion_label 空串保护
- **溯源报告增强**：Per-PSM header 展示 RT(min) + q-value；Summary 表含 Heavy 计数（H:TrapOnly/H:Shared/H:TargetOnly）；Candidate 表展示修饰；嵌合检测使用轻+重合并 shared fraction

**依赖**：`core`, `spectrum-io`（v3 新增）, `serde`, `serde_yaml`, `arrow`, `parquet`, `tracing`
**不依赖**：`rmcp`

### 3.14 `entrapment-cli`（bin crate）

**职责**：陷阱库分析的独立 CLI 工具（薄壳）。

**子命令**：`analyze` / `report` / `inspect`

**依赖**：`entrapment-analysis`, `clap`, `tracing`

### 3.13 `fdr` 扩展（三级 FDR）

在 MVP 的 PSM 级 FDR 基础上，新增：

| 模块 | 功能 |
|------|------|
| `peptide_fdr` | 肽段级 FDR：提取最优肽段分数 → TDC 竞争 → q-value |
| `protein_fdr` | 蛋白级 FDR：picked-protein 方法（target/decoy 配对竞争 → FDR） |

**Picked-protein 算法**：
1. Target 蛋白组与 `REV_` 前缀的 decoy 配对
2. 每对中高分者胜出（同分 → target 胜）
3. 对胜出者执行标准 TDC FDR 计算
4. 无配对 target 自动胜出（q=0）

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
    pub lower_offset: f64,  // m/z units
    pub upper_offset: f64,  // m/z units
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
    #[serde(default)]
    pub acquisition_mode: Option<AcquisitionMode>,  // DDA/DIA/Unknown
    #[serde(default = "default_max_variable_modifications")]  // default: 3, max: 10
    pub max_variable_modifications: u32,
    #[serde(default = "default_min_peptide_length")]  // default: 7
    pub min_peptide_length: u32,
    #[serde(default = "default_max_peptide_length")]  // default: 50
    pub max_peptide_length: u32,
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
    async fn search_with_spectra(&self, params: &SearchParams, spectra: Vec<Spectrum>)
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
| FDR 计算（target-decoy） | Rust | fdr crate（✅ 已实现） | 数值计算 |
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
              ┌──────────┬─────┼──────────┬──────────┬──────────────┐
              ▼          ▼     ▼          ▼          ▼              ▼
        ┌───────────┐ ┌──────────┐ ┌────────┐ ┌────────┐ ┌───────────┐ ┌─────────────────┐
        │spectrum-io│ │  param-  │ │search- │ │ report │ │   dia-    │ │  entrapment-    │
        │           │ │recommend │ │ engine │ │        │ │extraction │ │   analysis      │
        │(quick-xml │ │          │ │(tokio) │ │ (csv)  │ │           │ │(levenshtein,    │
        │ base64    │ │          │ │        │ │        │ │           │ │ k-mer index)    │
        │ flate2)   │ │          │ │        │ │        │ │           │ │                 │
        └─────┬─────┘ └────┬─────┘ └───┬────┘ └───┬────┘ └─────┬─────┘ └────────┬────────┘
              │            │           │           │            │                │
              └────────────┴───────┬───┴───────────┴────────────┘                │
                                   ▼                                             │
                            ┌────────────┐    ┌────────────────┐                 │
                            │ mcp-server │◀───│ entrapment-cli │◀────────────────┘
                            │   [bin]    │    │     [bin]       │
                            └────────────┘    └────────────────┘
```

依赖方向始终向上（library → core），mcp-server 在最底层聚合所有 library。
entrapment-cli 是独立二进制，仅依赖 entrapment-analysis。
**禁止**：library 之间相互依赖（spectrum-io 不依赖 param-recommend）。

---

## 7. MCP Tool 注册表（28 tools）

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
| `annotate_spectrum` | report + xic | run_id/file+peptide+charge, scan | AnnotateResult + HTML | 谱图注释（DIA 自动含 XIC+SILAC） |
| `extract_dia_precursors` | dia-extraction | file_path, params? | RunId + ExtractionSummary | DIA 前体提取 |
| `extract_spectrum_precursors` | dia-extraction | file_path, scan_number | SingleSpectrumExtractionResult | 单谱图母离子提取 |
| `extract_xic` | xic | run_id?, file_path?, scan_number, peptide?, charge? | HTML file path | XIC 碎片离子色谱图（支持 SILAC 轻重标记） |
| `import_search_results` | result-import | result_path, mzml_files, format? | RunId + ImportSummary | 导入外部搜索结果（DIA-NN / custom JSON） |
| `classify_entrapment_hits` | entrapment-analysis | results_file, config_file, target_fasta, output_dir? | EntrapmentSummary | 运行陷阱库分类流程（L0-L4 + HTML 报告） |
| `analyze_entrapment_stats` | entrapment-analysis | classified_file | DetailedStats | 从已分类 TSV 生成统计分析 |
| `find_similar_targets` | entrapment-analysis | peptide, target_fasta, max_mismatches? | Vec\<SimilarityHit\> | 查找肽段在 target 库中最相似序列 |
| `annotate_provenance` | entrapment-analysis | peptide, target_peptide, file_path, scan_number, fragment_tolerance?, chimera_threshold? | FragmentProvenance + is_chimeric | 对单个 trap PSM 进行碎片离子溯源 |
| `list_searches` | mcp-server | status_filter?, limit? | Vec\<SearchHistoryEntry\> | 搜索历史 |
| `get_search_status` | mcp-server | run_id | SearchProgress | 查询搜索进度 |
| `cancel_search` | mcp-server | run_id | SearchProgress | 取消搜索 |
| `prepare_search` | mcp-server | file_paths, database_path/organism?, user_hints? | SearchParams | 一键搜索准备（read_spectra + recommend_params + download_database） |
| `infer_proteins` | protein-inference | run_id/search_result | ProteinGroups | 蛋白推断（parsimony + razor + 蛋白级 FDR） |
| `diagnose_search` | mcp-server | run_id | DiagnosticReport | 搜索诊断（失败分析 / 质量评估） |
| `get_dia_cache_status` | dia-extraction | run_id | DiaCacheLocation | 检查 DIA 提取缓存状态（memory/disk/not_found） |
| `list_databases` | fasta-db | (无) | Vec\<DatabaseInfo\> | 列出内置 FASTA 数据库及下载状态 |
| `download_database` | fasta-db | database_id | DatabasePath | 下载 FASTA 数据库（UniProt HTTPS + 本地缓存） |
| `get_database_info` | fasta-db | database_id | DatabaseDetail | 查询已下载数据库详情（蛋白数、SHA256 等） |

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

1. ~~初始化 Workspace + core crate~~ ✅
2. ~~验证 rmcp `#[tool_router]` 宏在最小示例中的可用性~~ ✅
3. ~~实现 spectrum-io 的 mgf 解析~~ ✅

> MVP 已完成，BUG-1（碎片离子固定修饰）已修复。
> Post-MVP 功能已完成：异步搜索、索引访问、DIA 支持、XIC、外部结果导入、Biology Audit。
> 陷阱库分析 v1+v2 已完成：L0-L4 分级 + Levenshtein edit distance + k-mer 预筛 + SubstitutionType 注释 + HTML 报告 + CLI + 3 MCP tools。
> 统一标注+XIC 视图已完成：
> - 文件名 + Scan/RT 显示
> - DDA 自动跳过 XIC（基于 fragment trace 数据点判断，非窗口宽度阈值）
> - 客户端 SILAC 重计算引擎（raw peaks 嵌入 HTML）
> - 逐离子 L/H 开关网格 + 批量切换按钮
> - 详见 `tasks/001-mvp-proteomics-search-platform.md`
>
> **跨模块流程审计已完成（2026-04-21）：**
> - BUG-2~9 已修复（DIA precursor m/z 全链路 + RT 单位 + PPM 零值保护 + DIA 搜索安全 + score NaN + Decoy 统一）
> - SILAC 轻重标 MS2 分离验证通过（XIC 提取 / 标注 / MCP tool 三层均正确分离轻重窗口）
> - 740+ tests，0 warnings
>
> **代码审计修复已完成（2026-04-24）：**
> - BUG-10~16 已修复：entrapment Heavy 离子归属 + unwrap 安全化（tools.rs/unimod.rs/diann.rs）+ HTML 注入防护（xic_visualize.rs/report.rs）+ 魔法数字常量化
> - Plotly.js 全项目统一为 2.35.2（mirror_plot.rs / multi_report.rs / entrapment_report.html）
> - 详见 `docs/superpowers/plans/2026-04-24-audit-fixes.md`
>
> **溯源报告增强已完成（2026-04-24）：**
> - MultiTargetProvenance 新增 `trap_retention_time_min` + `trap_q_value` 字段
> - Per-PSM Header 展示 RT(min) 和 q-value
> - Summary 表新增 Heavy 计数列（H:TrapOnly / H:Shared / H:TargetOnly）
> - Candidate 表替换冗余 Spectrum File 列为 Modifications 列
> - 嵌合检测改用轻+重合并 shared fraction
> - 838 tests，0 warnings
