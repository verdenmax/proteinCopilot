# 全链路结构化追踪设计文档

**日期**: 2026-04-28
**状态**: Approved
**范围**: 15 crates · ~90 span · 23 热循环进度 · subscriber 增强

---

## 1. 问题陈述

ProteinCopilot 已有 `tracing` 基础设施（workspace 依赖 + stderr subscriber），但利用率极低：

- **89 条日志** 全是 `warn!/debug!` 用于错误/边界情况
- **零个 `#[instrument]`** — 无自动耗时追踪
- **零个 span 嵌套** — 无法看到调用链父子关系
- **3 个 crate 完全无日志**：`fdr`、`report`、`param-recommend`
- **手动 `Instant::now()` 散落 5 处**，不统一
- 搜索运行时用户看到"一片空白"，无法知道执行到哪、哪步慢、还要多久

## 2. 设计目标

1. 用户在 info 级别即可看到完整的执行流程、决策上下文和耗时
2. 异常时可精确定位到阶段和函数
3. 热循环提供进度百分比和 ETA
4. 不影响 MCP JSON-RPC 协议（所有日志走 stderr）
5. 支持 `PROTEIN_LOG_JSON=1` 切换 JSON 格式输出
6. 零运行时开销（tracing 在编译期过滤不活跃 level）

## 3. 三层追踪架构

### Layer 1: MCP Tool Span（28 个）

每个 `#[rmcp::tool]` 方法入口创建 span，记录：
- `tool_name`: 工具名称
- `run_id`: 运行 ID（如适用）
- 输入摘要（文件数、数据库路径、引擎名等）
- 完成时记录耗时和结果摘要

实现方式：在每个 tool 方法开头手动创建 `info_span!`，方法结束时 span 自动 drop 记录耗时。不用 `#[instrument]` 是因为 rmcp 宏与 tracing instrument 宏存在兼容性问题，且 tool 参数包含大量 JSON 不适合自动 skip。

```rust
pub async fn run_search(&self, params: Parameters<RunSearchInput>) -> Result<...> {
    let span = tracing::info_span!("mcp_tool", name = "run_search",
        run_id = tracing::field::Empty);
    let _enter = span.enter();
    // ... 拿到 run_id 后 ...
    span.record("run_id", run_id.to_string().as_str());
    tracing::info!(engine = %engine, files = input_files.len(), "started");
    // ... 执行 ...
    tracing::info!(psms_1pct = result.psms_at_1pct, total_sec = elapsed, "completed");
}
```

### Layer 2: Library Function Span（~35 个）

核心库函数入口使用 `#[instrument]`，skip 大型数据参数：

```rust
#[tracing::instrument(skip(spectra, on_progress), fields(spectrum_count = spectra.len()))]
pub fn match_all_spectra(spectra: &[Spectrum], ...) -> Vec<Psm> { ... }
```

关键函数列表（按 crate）：

**spectrum-io**（8 span）：
- `IndexedMzMLReader::open` — fields: file_path, file_size, index_source, scan_count
- `build_index_by_byte_scan` — fields: file_size, scan_count
- `disk_cache::load` / `disk_cache::save` — fields: path, scan_count
- `for_each_spectrum` — fields: file_path, total_spectra
- `read_summary` — fields: ms1_count, ms2_count, mz_range, rt_range
- `read_spectrum` — fields: scan_number
- `build_index_by_scanning` — fields: file_path

**search-engine**（15 span）：
- `SimpleSearchEngine::search` / `SageSearchEngine::search` — fields: engine, spectrum_count
- `parse_fasta` — fields: path, protein_count
- `digest` — fields: protein_count, peptide_count, enzyme
- `match_spectrum` — 不加 instrument（调用次数太多），由外层批处理 span 覆盖
- `match_all_spectra` — fields: spectrum_count, matched_count
- `score_psms` — fields: candidate_count
- `generate_b_ions` / `generate_y_ions` — fields: sequence_len, charge
- `match_fragments` — fields: theoretical_count, matched_count
- `annotate_spectrum_impl` — fields: scan, peptide, charge

**xic**（6 span）：
- `extract_xic_unified` — fields: peptide, precursor_mz, scans_planned
- `plan_scans` — fields: rt_center, rt_window, n_cycles
- `extract_ms1_xic` — fields: precursor_mz, points_collected
- `extract_ms2_xic` — fields: ion_count, points_collected
- `render_xic_html` — fields: output_path
- `build_xic_data` — fields: light_scans, heavy_scans

**dia-extraction**（5 span）：
- `extract_dia_precursors` — fields: file_path, ms1_count, ms2_count
- `detect_acquisition_mode` — fields: detected_mode
- `correlate_ms1_ms2` — fields: ms1_count, ms2_count
- `extract_isotope_patterns` — fields: spectrum_count, candidates_found
- `build_pseudo_spectra` — fields: count

**fdr**（3 span）：
- `calculate_fdr` — fields: psm_count, target_count, decoy_count, psms_at_1pct
- `calculate_peptide_fdr` — fields: peptide_count
- `calculate_protein_fdr` — fields: protein_group_count

**report**（4 span）：
- `generate_summary` — fields: psm_count, id_rate, median_ppm
- `export_tsv` — fields: output_dir, psm_count, peptide_count, protein_count
- `export_json` — fields: output_path
- `render_entrapment_report` — fields: output_path

**result-import**（5 span）：
- `detect_format` — fields: file_path, detected_format
- `DiannParser::parse` / `PfindParser::parse` / `CustomJsonParser::parse` — fields: file_path, row_count, valid_psms
- `match_scans` — fields: psm_count, matched_count, unmatched_count

**protein-inference**（6 span）：
- `infer_proteins` — fields: psm_count, q_value_threshold
- `build_peptide_protein_map` — fields: psm_count, unique_peptides
- `run_parsimony` — fields: peptide_count, protein_count, groups_formed
- `assign_razor_peptides` — fields: peptide_count
- `calculate_coverage` — fields: protein_count, mean_coverage
- `score_protein_groups` — fields: group_count

**entrapment-analysis**（8 span）：
- `classify_entrapment_hits` — fields: results_file, target_fasta
- `load_search_results` — fields: file_path, psm_count
- `TargetDigestIndex::from_fasta` — fields: fasta_path, protein_count, peptide_count
- `classify_all` — fields: psm_count, trap_count, target_count
- `find_similar_targets` — fields: query_peptide, max_mismatches, matches_found
- `trace_provenance_batch` — fields: psm_count, mzml_dir
- `generate_report` — fields: output_dir
- `analyze_stats` — fields: classified_file

**param-recommend**（2 span）：
- `recommend_params` — fields: file_path, confidence
- `apply_preset` — fields: preset_name

### Layer 3: Hot Loop Progress（23 处）

在处理大量数据的循环内，每 N 次迭代输出一条 info 日志：

```rust
for (i, spectrum) in spectra.iter().enumerate() {
    // ... process ...
    if (i + 1) % progress_interval == 0 || i + 1 == total {
        let elapsed = start.elapsed().as_secs_f64();
        let rate = (i + 1) as f64 / elapsed;
        let remaining = (total - i - 1) as f64 / rate;
        tracing::info!(
            progress = i + 1,
            total = total,
            pct = format!("{:.1}%", (i + 1) as f64 / total as f64 * 100.0),
            rate_per_sec = format!("{:.0}", rate),
            eta_sec = format!("{:.0}", remaining),
            "matching spectra"
        );
    }
}
```

热循环位置和进度间隔：

| Crate | 循环位置 | 迭代内容 | 进度间隔 |
|-------|---------|---------|---------|
| search-engine | `match_all_spectra` | 谱图匹配 | 每 500 条 |
| search-engine | `digest` | 蛋白消化 | 每 1000 个蛋白 |
| search-engine | `score_psms` | PSM 打分 | 每 5000 条 |
| search-engine | `parse_fasta` | FASTA 读取 | 每 5000 个蛋白 |
| spectrum-io | `for_each_spectrum` | 谱图流式读取 | 每 1000 条 |
| spectrum-io | `build_index_by_byte_scan` | 字节扫描建索引 | 每 5000 条 |
| xic | `extract_ms2_xic` 内循环 | XIC 离子提取 | 每扫描完成 |
| dia-extraction | `extract_isotope_patterns` | 同位素模式提取 | 每 500 条 MS2 |
| dia-extraction | `correlate_ms1_ms2` | MS1-MS2 关联 | 每 1000 条 |
| dia-extraction | `build_pseudo_spectra` | 伪谱构建 | 每 500 条 |
| result-import | `match_scans` | 扫描匹配 | 每 1000 条 PSM |
| result-import | parse 循环 | 结果解析 | 每 5000 行 |
| entrapment | `classify_all` | PSM 分类 | 每 500 条 |
| entrapment | `TargetDigestIndex::from_fasta` | 目标消化 | 每 1000 个蛋白 |
| entrapment | `trace_provenance_batch` | 溯源分析 | 每 50 条 |
| entrapment | `find_similar` 外循环 | 相似性搜索 | 每 100 条 |
| protein-inference | `build_peptide_protein_map` | 肽段-蛋白映射 | 每 5000 条 PSM |
| protein-inference | `run_parsimony` | 简约法迭代 | 每轮迭代 |
| protein-inference | `assign_razor_peptides` | Razor 分配 | 每 1000 条 |
| protein-inference | `calculate_coverage` | 覆盖率计算 | 每 500 个蛋白 |
| fdr | `calculate_fdr` 排序+扫描 | FDR 计算 | 每 5000 条 |
| report | `export_tsv` 写入循环 | TSV 导出 | 每 5000 行 |
| fasta-db | `download_database` | 下载进度 | 每 1MB |

## 4. Subscriber 增强

### 当前 subscriber（main.rs）

```rust
tracing_subscriber::fmt()
    .with_env_filter(EnvFilter::from_default_env())
    .with_writer(std::io::stderr)
    .init();
```

### 增强后

```rust
use tracing_subscriber::{fmt, EnvFilter, layer::SubscriberExt, util::SubscriberInitExt};

let env_filter = EnvFilter::try_from_default_env()
    .unwrap_or_else(|_| EnvFilter::new("info"));

let use_json = std::env::var("PROTEIN_LOG_JSON").map_or(false, |v| v == "1");

if use_json {
    tracing_subscriber::registry()
        .with(env_filter)
        .with(fmt::layer()
            .json()
            .with_writer(std::io::stderr)
            .with_span_events(fmt::format::FmtSpan::CLOSE)
            .with_target(true)
            .with_timer(fmt::time::uptime()))
        .init();
} else {
    tracing_subscriber::registry()
        .with(env_filter)
        .with(fmt::layer()
            .with_writer(std::io::stderr)
            .with_span_events(fmt::format::FmtSpan::CLOSE)
            .with_target(true)
            .with_timer(fmt::time::uptime()))
        .init();
}
```

关键变化：
- `FmtSpan::CLOSE`：span 结束时自动输出耗时
- `with_target(true)`：显示 crate::module 路径
- `with_timer(uptime())`：显示进程启动以来的时间
- JSON 模式：`PROTEIN_LOG_JSON=1` 启用结构化输出
- 默认级别：`info`（而非之前的无默认）

### 依赖补全

需要为以下 crate 添加 `tracing` workspace 依赖：
- `fdr` — Cargo.toml 添加 `tracing.workspace = true`
- `report` — 同上
- `param-recommend` — 同上
- `core` — 同上（用于 diagnostics span）

`mcp-server` 需要额外添加 `tracing-subscriber` 的 `json` feature：
```toml
tracing-subscriber = { workspace = true, features = ["env-filter", "json"] }
```

## 5. Instant::now() 迁移

现有 5 处 `Instant::now()` 手动计时需要迁移为 tracing span：

| 文件 | 当前用法 | 迁移方案 |
|------|---------|---------|
| `search-engine/simple_engine.rs:185` | 搜索总耗时 | 用 `#[instrument]` + `FmtSpan::CLOSE` 自动记录 |
| `search-engine/simple_engine.rs:477` | 阶段耗时 | 改为 `info_span!` 包裹每个阶段 |
| `search-engine/adapters/sage/mod.rs` | Sage 搜索耗时 | 用 `#[instrument]` 替代 |
| `mcp-server/tools.rs:1585` | run_search 任务耗时 | 用 tool span 替代 |
| `mcp-server/tools.rs:1598` | 同上 | 合并到 tool span |

保留 `Instant::now()` 用于 `SearchProgress` 回调中的 `elapsed_sec` 计算（这是给 MCP client 的结构化进度数据，不是日志）。

## 6. 决策上下文 info 日志

以下决策点在 info 级别输出：

| 决策点 | 输出内容 |
|--------|---------|
| 索引来源 | `index source=disk_cache` / `byte_scan` / `native_index` |
| 缓存命中 | `cache hit path="xxx.idx"` / `cache miss, rebuilding` |
| 缓存保存 | `saved disk cache path="xxx.idx" scans=N` |
| 引擎选择 | `engine="Sage"` / `"SimpleSearch"` |
| 参数推荐 | `enzyme=Trypsin tolerance=20ppm confidence=0.92` |
| 数据库 | `database="human_swissprot" proteins=20380` |
| 采集模式 | `acquisition_mode=DIA isolation_window=25Da` |
| DIA 提取方式 | `output_mode=pseudo candidates=45000` |
| 格式检测 | `detected format=diann_parquet` |
| FDR 结果 | `psms@1%FDR=28500 peptides@1%=18000 proteins@1%=2850` |

## 7. 异常/慢速告警

在 info 之外，自动检测异常情况并输出 warn：

| 条件 | 级别 | 消息 |
|------|------|------|
| 热循环速率下降 > 50% | WARN | `slow batch: rate=120/s (normal=530/s)` |
| 单个谱图匹配耗时 > 100ms | WARN | `slow spectrum matching: scan=N elapsed=150ms` |
| FDR 1% 下匹配率 < 10% | WARN | `low identification rate: 8.5%` |
| 索引重建耗时 > 10s | WARN | `slow index rebuild: 15.2s for 7.5GB file` |
| 内存中谱图数 > 100K | WARN | `high memory usage: 120K spectra in memory` |

## 8. 不变的约束

- 所有日志输出到 **stderr**，绝不污染 stdout（MCP JSON-RPC 通道）
- `tracing` 在编译期过滤不活跃 level，零运行时开销
- 不引入 OpenTelemetry 或外部 collector（YAGNI，预留接口即可）
- 不改变现有 `SearchProgress` 回调机制（它服务于 MCP `get_search_status` tool）
- 现有 `warn!/debug!` 日志保留不动，只新增不删除
