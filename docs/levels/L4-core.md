# L4 — core crate 核心路径（共享数据结构与领域 trait）

承接 [L2 架构](L2-architecture.md)。本篇逐文件精讲 `crates/core`：它定义全平台共享的领域类型与 `SearchEngineAdapter` trait，是依赖图的根。所有结构体、枚举、trait、签名、字段均核对自源码（标注行号），不臆造。

## 1. 用途、位置与依赖

一句话：**core 是 ProteinCopilot 的"领域词典"——把谱图、参数、结果、AI 决策、引擎抽象、运行元数据固化为带校验的 serde 类型，供各上层 crate 复用。**

位置：依赖图最底层，**无任何内部 crate 依赖**（与 `fasta-db` 同为根），被各上层 crate 广泛依赖（见 L2 第 2 节）。core 自身不读文件、不打分、不调 LLM——只提供类型、`validate()` 校验与少量纯函数；一切数值计算留给下游确定性完成。

外部依赖（均来自 workspace）：

```text
serde + serde_json  ->  全部类型 Serialize/Deserialize
schemars            ->  JsonSchema，供 MCP 工具暴露 schema
thiserror           ->  每模块的 *Error 枚举
chrono              ->  RunMetadata.created_at : DateTime<Utc>
uuid                ->  run_id : Uuid (v4)
async-trait         ->  SearchEngineAdapter 的 async 方法
```

`lib.rs` 仅做模块导出（12 个 `pub mod`，无 re-export），下游按 `protein_copilot_core::spectrum::Spectrum` 全路径引用。

贯穿约定（读这些类型前需知道的几条规矩）：

- **校验前置**：几乎每个数据类型都带 `validate()`，把"非有限值、范围越界、计数倒挂、空必填"在构造或反序列化后就拦下，下游拿到的恒为合法数据。
- **serde 默认值**：新增字段一律配 `#[serde(default = ...)]` 保证旧 JSON 仍可反序列化，例如 `SearchParams` 的 `max_variable_modifications=3`、`min_peptide_length=7`、`max_peptide_length=50`、`engine=None`（search_params.rs:228-282）。
- **领域单位固定**：质量偏差用 `f64`（Da/ppm 由 `ToleranceUnit` 标注），保留时间统一"分钟"（`retention_time_min`），`scan_number` 从 1 起（0 触发 `ZeroScanNumber`，spectrum.rs:231）。
- **双 derive**：所有对外类型同时 `Serialize + Deserialize + JsonSchema`，MCP 工具据此自动产出结构化 I/O，不走自由文本。

## 2. 关键类型清单（按模块）

下表覆盖 12 个模块的核心导出，括注字段为最常被下游读写的项；完整字段与约束以源码为准。

| 模块 | 类型 | 一行职责 |
|------|------|---------|
| `spectrum` | `Spectrum` | 单张谱图：scan/ms_level/RT/precursors + mz/intensity 双数组 |
| | `SpectrumSummary` | 谱图文件统计摘要，推参的主输入 |
| | `IsolationWindow` / `PrecursorInfo` | 隔离窗口；母离子（mz/charge/window/source_scan） |
| | `MsLevel` / `AcquisitionMode` / `SpectrumFormat` | MS1/MS2/Other；DDA/DIA/Unknown；MzML/Mgf/Pfb |
| | `SpectrumFileInfo` / `SpectrumError` | 磁盘文件元数据；谱图校验错误 |
| `search_params` | `SearchParams` | 完整搜索配置（酶/修饰/容差/库/decoy/长度/引擎） |
| | `Enzyme` / `ModPosition` / `Modification` | 酶；修饰位点；修饰（名/mass_delta/残基/位点） |
| | `MassTolerance` / `ToleranceUnit` / `DecoyStrategy` | 容差值+单位；Ppm/Da；Reverse/Shuffle/None |
| `search_result` | `SearchResult` | 一次运行的总输出：PSM/肽/蛋白/summary/metadata |
| | `Psm` | 肽谱匹配（scan/序列/charge/score/q_value/is_decoy） |
| | `PeptideResult` / `ProteinResult` | 肽级、蛋白级聚合 |
| | `SearchResultSummary` | 1% FDR 统计摘要，供 AI 解释 |
| `ai_decision` | `AiDecision<T>` | AI 决策统一包装（decision/confidence/explanation...） |
| `engine` | `SearchEngineAdapter` (trait) | 全部搜索引擎的统一异步接口 |
| | `EngineInfo` / `HealthStatus` | 引擎元数据；Healthy/Degraded/Unavailable |
| `label` | `LabelType` | SILAC/自定义重标，附 mass delta 计算 |
| `error` | `CoreError` / `ErrorReport` | MCP 边界统一错误 + 可序列化错误体 |
| `run_metadata` | `RunMetadata` / `RunStatus` | 运行溯源（run_id/时间/参数/引擎/状态） |
| `protein_group` | `ProteinGroup` / `InferenceResult` | 蛋白分组；蛋白推断完整结果 |
| `progress` | `SearchProgress` / `ProgressCallback` | 进度数据；进度回调类型别名 |
| `diagnostics` | `SearchDiagnostics` | 分阶段计时 + 异常检测 + 修复建议 |
| `util` | `is_decoy_accession` / `DECOY_PREFIXES` / `compute_median` | decoy 前缀判定；中位数纯函数 |

## 3. 核心 trait 与函数签名（核对自源码）

SearchEngineAdapter（engine.rs:79-135，`#[async_trait]`，约束 `Send + Sync`）:

```rust
async fn search(&self, params: &SearchParams, input_files: &[PathBuf],
    on_progress: ProgressCallback, diagnostics: &mut SearchDiagnostics)
    -> Result<SearchResult, CoreError>;                                  // :88
fn engine_info(&self) -> EngineInfo;                                     // :97
async fn health_check(&self) -> Result<HealthStatus, CoreError>;        // :103
// 默认实现：返回 SearchEngineError（引擎需 opt-in）
async fn search_with_spectra(&self, params: &SearchParams, spectra: Vec<Spectrum>,
    on_progress: ProgressCallback, diagnostics: &mut SearchDiagnostics)
    -> Result<SearchResult, CoreError>;                                 // :108
// 默认实现：no-op Ok(())
async fn cancel(&self, _run_id: Uuid) -> Result<(), CoreError>;         // :132
```

label.rs 纯函数（无副作用、不调 LLM）:

```rust
pub fn total_heavy_delta(peptide_sequence: &str, label: &LabelType) -> f64;   // label.rs:59
pub fn residue_heavy_delta(residues: &[char], label: &LabelType) -> f64;      // :37
pub fn compute_heavy_precursor_mz(light_mz: f64, charge: i32,
    peptide_sequence: &str, label: &LabelType) -> f64;                        // :65
```

util.rs / progress.rs:

```rust
pub const DECOY_PREFIXES: &[&str] = &["REV_", "SHUF_", "DECOY_", "REVERSED_"]; // util.rs:38
pub fn is_decoy_accession(accession: &str) -> bool;                           // :41
pub fn compute_median(sorted: &[f64]) -> f64;                                 // :10
pub fn compute_median_u32(sorted: &[u32]) -> u32;                             // :26
pub type ProgressCallback = Box<dyn Fn(SearchProgress) + Send + Sync>;        // progress.rs:40
pub fn noop_progress() -> ProgressCallback;                                   // :43
```

## 4. 简化源码片段

(1) AiDecision<T>（ai_decision.rs:70-110）——泛型包装，所有 AI 决策必走此结构:

```rust
pub struct AiDecision<T> {
    pub decision: T,             // 具体决策值（如 SearchParams 或 String）
    pub confidence: f64,         // [0.0, 1.0]
    pub explanation: String,     // 为何这样决策（validate 要求非空）
    pub input_summary: String,   // 基于何种输入（validate 要求非空）
    pub alternatives: Vec<String>,
    pub evidence: Vec<String>,
}
// validate(): confidence 有限且落 [0,1]，explanation/input_summary 去空白后非空
```

(2) LabelType 的 delta 计算（label.rs:37-56）——纯计数累加:

```rust
pub enum LabelType {
    Silac { heavy_k_delta: f64, heavy_r_delta: f64 },  // standard: K+8.014199 / R+10.008269
    Custom { residue_deltas: Vec<(char, f64)> },
}
// residue_heavy_delta: 数 K/R（或自定义残基）个数 * 对应 delta 再累加
let count_k = residues.iter().filter(|&&c| c == 'K' || c == 'k').count() as f64;
let count_r = residues.iter().filter(|&&c| c == 'R' || c == 'r').count() as f64;
count_k * heavy_k_delta + count_r * heavy_r_delta
```

(3) CoreError 在 MCP 边界统一错误（error.rs:24-105 + From 实现）:

```rust
pub enum CoreError {
    SpectrumParseError { format, detail, suggestion },
    InvalidSearchParams { field, reason, suggestion },
    SearchEngineError  { engine, detail, suggestion },
    FileNotFound { path },  UnsupportedFormat { format, supported },
    SshConnectionError { host, detail },
    ResultParseError   { engine, file, detail },
    ValidationError    { context, detail, suggestion },
}
// 各模块 *Error 经 From 收敛到 CoreError；suggestion() 给可执行建议；
// ErrorReport::from(&err) 产出 { category, message, suggestion } 供 MCP JSON 返回。
```

## 5. 典型调用链（谁用 core 的什么）

```text
spectrum-io       -> 产出 core::spectrum::{Spectrum, SpectrumSummary}
param-recommend   -> 吃 SpectrumSummary，产出 AiDecision<SearchParams>
search-engine     -> impl SearchEngineAdapter；search() 返回 core::SearchResult
   |                 内部用 SearchDiagnostics 计时，ProgressCallback 报进度
fdr               -> 读写 Psm.q_value / is_decoy（util::is_decoy_accession）
protein-inference -> 产出 ProteinGroup / InferenceResult
report            -> 读 SearchResult，按 1% FDR 重算 SearchResultSummary
xic               -> 用 label::total_heavy_delta 算重标母离子 m/z
mcp-server        -> 全部类型经 JsonSchema 暴露为工具 I/O；CoreError -> ErrorReport
```

要点：core 只被依赖、不依赖任何内部 crate；FDR、打分、中位数、质量偏差等数值全在下游确定性算出，LLM 仅消费 schema 化摘要。

## 6. 测试入口

测试全部为模块内 `#[cfg(test)] mod tests`（无独立 `tests/` 目录）。覆盖 11 个文件，`mod tests` 起始行：spectrum.rs:464、search_params.rs:350、search_result.rs:507、ai_decision.rs:118、error.rs:230、label.rs:76、util.rs:46、engine.rs:142、progress.rs:48、run_metadata.rs:166、diagnostics.rs:409。`protein_group.rs` 与 `lib.rs` 无测试模块（前者校验逻辑由下游 crate 覆盖）。

跑法（必须 `--offline`，联网会挂）:

```text
cargo test -p protein-copilot-core --offline
# 实测：191 passed; 0 failed（含 serde 往返、validate 边界、From 收敛、label delta）
```

回到 [README](README.md)。
