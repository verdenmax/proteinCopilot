# L4 — param-recommend（搜索参数推荐）

承接 [L2](L2-architecture.md)。本篇逐源码讲清 `crates/param-recommend`：它如何用**确定性规则**把一份 `SpectrumSummary` 翻译成带置信度与理由的 `AiDecision<SearchParams>`。所有签名、阈值、默认值均核对自源码（标注行号），无臆造。

## 1. 用途、位置与依赖

一句话：**纯库 crate，输入谱图统计 + 可选用户提示，输出可解释的搜索参数推荐——同输入必同输出，不调 LLM、不联网、无随机。**

```text
  上游                   param-recommend                        下游
SpectrumSummary --> ParamRecommender::recommend --> AiDecision<SearchParams> --> search-engine
UserHints(可选) -->   (确定性规则: DIA/仪器/容差/修饰)        decision/confidence/        (run_search)
                                                          explanation/evidence/...
```

在 L2 的 15 crate 依赖图里，它夹在 spectrum-io 与 search-engine 之间：spectrum-io 产出 `SpectrumSummary` 喂进来，本 crate 的 `AiDecision<SearchParams>` 经 mcp-server 解包后交给 search-engine 跑搜索。LLM 只负责把用户自然语言意图翻成 `UserHints`，所有数值决策留在这里——这正是"Rust 不调 LLM、LLM 不算参数"分工的落点。

依赖（Cargo.toml）：`core`（`SpectrumSummary` / `SearchParams` / `AiDecision`）+ serde + schemars + thiserror + tracing；dev-dep 加 `spectrum-io`（集成测试读真实 fixture）。模块职责清晰：`lib.rs` 暴露入口、`preset.rs` 定义 5 个预设、`rules.rs` 装全部规则、`hints.rs` 定义提示结构、`error.rs` 定义错误。

## 2. 对外 API

`ParamRecommender` 是无字段单元结构（lib.rs:44），两个方法：

```rust
// lib.rs:56  规则全在 rules::recommend
pub fn recommend(&self, summary: &SpectrumSummary, hints: Option<&UserHints>)
    -> Result<AiDecision<SearchParams>, ParamRecommendError>;
// lib.rs:65
pub fn list_presets() -> Vec<SearchPreset>;
```

- `UserHints`（hints.rs:16）：`experiment_type` / `instrument_type` / `enzyme` / `custom_notes`，全为 `Option`；`custom_notes` 仅存不被规则消费。
- `SearchPreset`（preset.rs:11）：`name` / `description` / `params` / `applicable_scenarios`；`with_database(&str) -> SearchParams`（preset.rs:27）替换占位库路径。
- `ParamRecommendError`（error.rs:7）：`EmptySummary` 与 `InvalidSummary { field, detail }`，并 `From` 成 `CoreError::ValidationError`。

内置 5 个预设（preset.rs:268 `all_presets`，顺序固定）：

| name | 关键修饰（核对自源码） | 容差 |
|------|------|------|
| `standard` | fixed Carbamidomethyl(C)+57.021464；var Oxidation(M)+15.994915 | 10/20 ppm |
| `phospho` | + var Phospho(STY)+79.966331 | 10/20 ppm |
| `tmt` | fixed + TMT6plex(K) 与 N-term 各 +229.162932 | 10/20 ppm |
| `silac` | + var Label 13C(6)15N(2)-K +8.014199、13C(6)15N(4)-R +10.008269 | 10/20 ppm |
| `open` | 仅 Carbamidomethyl(C)+Oxidation(M) | 500 Da / 20 ppm |

全预设共享默认：`Enzyme::Trypsin`、`missed_cleavages 2`、`DecoyStrategy::Reverse`、`max_variable_modifications 3`、`min_peptide_length 7`、`max_peptide_length 50`、`engine None`、`database_path "<database_path>"` 占位（调用方必须替换）。

## 3. 确定性规则（逐条对齐源码）

- **DIA 检测**（rules.rs:113）：`median_isolation_window_da > 5.0` 即 `acquisition_mode = Some(DIA)`，否则保持 `None`。
- **仪器推断**（rules.rs:131）：hint 优先——含 `orbitrap`/`exactive`/`hires` -> 高分辨，含 `tof`/`trap`/`lowres` -> 低分辨；否则按下表打分。

| 特征 | 高分加分 | 低分加分 |
|------|------|------|
| `mz_range[1] > 1800` / `> 1500` / `< 1200` | +2 / +1 | — / — / +1 |
| `median_peaks_per_spectrum > 300` / `> 200` / `< 100` | +2 / +1 | — / — / +1 |

  判定（rules.rs:165）：`hi >= 2` -> 高分辨；否则 `lo >= 2` -> 低分辨；否则 General。

- **容差**（apply_tolerance，rules.rs:224；open 搜索跳过以保 500 Da）：高分辨 10 ppm / 20 ppm；低分辨 20 ppm / 0.1 Da；General 15 ppm / 20 ppm。
- **酶**：默认 Trypsin；`hints.enzyme` 直接覆盖（rules.rs:54）。
- **预设选择与修饰合并**（select_preset，rules.rs:178）：experiment_type 小写后子串匹配 `phospho`/`tmt`/`silac` 选预设，命中 `open` 时以 open 预设（500 Da）为基底，再把对应实验类型的修饰用 `merge_modifications`（rules.rs:207）叠加上去；去重按"全身份"（name+mass_delta+residues+position）而非仅名字，故 TMT-K 与 TMT-N端虽同名也都保留。
- **置信度**（compute_confidence，rules.rs:366）：基线高/低分辨 0.80、General 0.70；`experiment_type` +0.10、`instrument_type` +0.10、`enzyme` +0.05；末尾 `.min(0.95)` 封顶。
- **解释与证据**：build_explanation（rules.rs:265）把仪器判定、容差、实验类型、DIA 提示连成自然语言；build_evidence（rules.rs:301）列出 m/z 区间、中位峰数、谱图计数、按电荷排序的分布；build_alternatives（rules.rs:339）反向列出未选实验类型供 LLM 备选。
- **语义冲突**（detect_conflicts，rules.rs:391）：NonSpecific 却 `missed_cleavages > 0`；非 open 但 `precursor > 100` 配过窄碎片（Da<0.05 或 ppm<10）；`precursor > 100` 且 `missed_cleavages > 3`。命中即把 "Warnings:" 段追加进 explanation。

## 4. 简化源码片段

推荐主流程骨架（rules.rs:20 起）：

```rust
pub(crate) fn recommend(summary, hints) -> Result<AiDecision<SearchParams>, _> {
    if summary.is_empty() { return Err(EmptySummary); }          // 空文件直接拒
    summary.validate().map_err(|e| InvalidSummary { field: "summary", detail: e.to_string() })?;

    let experiment_type = hints.and_then(|h| h.experiment_type.as_deref()).unwrap_or("standard");
    let lower = experiment_type.to_lowercase();
    let is_open = lower.contains("open");

    let mut base = select_preset(&lower, is_open);              // open 基底再叠加 phospho/tmt/silac 修饰
    let instrument = infer_instrument(summary, hints);
    if !is_open { apply_tolerance(&mut base, &instrument); }    // open 保留 500 Da
    if let Some(e) = hints.and_then(|h| h.enzyme.clone()) { base.enzyme = e; }
    if detect_dia(summary) { base.acquisition_mode = Some(AcquisitionMode::DIA); }
    /* ... 构造 explanation / evidence / alternatives / confidence ... */
}
```

DIA 判定 + 仪器打分（rules.rs:113 / :131）：

```rust
fn detect_dia(s) -> bool { s.median_isolation_window_da.map(|w| w > 5.0).unwrap_or(false) }

fn infer_instrument(s, hints) -> InstrumentClass {
    // hint 命中即早退（略）
    let (mut hi, mut lo) = (0, 0);
    if s.mz_range[1] > 1800.0 { hi += 2 } else if s.mz_range[1] > 1500.0 { hi += 1 }
        else if s.mz_range[1] < 1200.0 { lo += 1 }
    if s.median_peaks_per_spectrum > 300 { hi += 2 } else if s.median_peaks_per_spectrum > 200 { hi += 1 }
        else if s.median_peaks_per_spectrum < 100 { lo += 1 }
    if hi >= 2 { HighResolution } else if lo >= 2 { LowResolution } else { General }
}
```

收尾构造 `AiDecision<SearchParams>`（rules.rs:96）：

```rust
Ok(AiDecision {
    decision: base,                 // 推荐出的 SearchParams
    confidence,                     // compute_confidence(...).min(0.95)
    explanation: final_explanation, // 命中冲突则带 "Warnings:" 段
    input_summary,                  // "N spectra, m/z [..], median K peaks/spectrum, RT [..] sec"
    alternatives,                   // 反向列出未选实验类型（open/phospho/tmt/silac/standard）
    evidence,                       // m/z 区间、峰数、谱图计数、电荷分布、仪器判定...
})
```

## 5. 调用链（mcp-server）

三个 tool 共用同一引擎，均经 `From<ParamRecommendError> for CoreError` 转结构化错误：

- `recommend_params`（tools.rs:1468）：summary 来自入参或 `get_or_create_reader(path).read_summary()`，调 `ParamRecommender.recommend(&summary, hints.as_ref())`，再按需注入 `database_path`，返回 `Json<AiDecision<SearchParams>>`。
- `list_presets`（tools.rs:1516）：`ParamRecommender::list_presets()` 包成 `PresetsResponse`。
- `prepare_search`（tools.rs:3622）：读首个文件 summary -> `recommend` -> 设 `engine` -> 按 `database_path` 或 `organism` 解析 FASTA，产出可直喂 `run_search` 的参数。`run_search`（tools.rs:1832）在缺参时也走同一 `recommend` 自动补全。

## 6. 测试入口

```bash
cargo test -p protein-copilot-param-recommend --offline
```

35 单测（lib / preset / rules 模块内）+ 4 集成（tests/integration.rs，经 spectrum-io 读 `small.mgf` / `small.mzml` 全链路）全绿。

—— 返回 [README](README.md)。
