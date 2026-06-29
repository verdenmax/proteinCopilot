# L4 — report crate（结果摘要 + 导出）

承接 [L2](L2-architecture.md)。本篇逐 crate 讲清 `crates/report`：把搜索引擎产出的 `SearchResult`
折算为可供 LLM 解释的统计摘要 `SearchResultSummary`，并落盘为 TSV / JSON / 元数据文件；附带把
谱图注释与 XIC 渲染为自包含 HTML。

## 1. 用途 + 位置 + 依赖

纯计算 crate — 不碰 MCP、网络、LLM，只做"格式化与导出"。

- 上游：`core::search_result::SearchResult`（PSM/肽/蛋白三级 + `summary` + `metadata`），由
  `search-engine` 或 `result-import` 产出并缓存于 mcp-server。
- 下游产物：内存里的 `SearchResultSummary`（交给 LLM 解读）+ 五个落盘文件（交给研究者/程序消费）
  + HTML 可视化。
- 依赖：`core`（数据结构、`util::compute_median`、`RunMetadata`）、`search-engine`
  （`annotate::SpectrumAnnotation`）、`xic`（`XicData`、`PlotlyMode`）。
- 位置：`crates/report/src/`；唯一使用方 `crates/mcp-server/src/tools.rs`。

模块边界：`summary.rs`（摘要）、`export.rs`（TSV/JSON/元数据）、`visualize.rs` + `xic_visualize.rs`
+ `unified_visualize.rs` + `xic3d_*`（HTML）、`error.rs`（`ReportError`）。

要点：`core` 侧 `SearchResult.summary` 是引擎计算、未做 FDR 过滤的粗摘要；本 crate 的
`generate_summary` 才按 1% FDR 口径重算，是结果落盘与 LLM 解读前的最后一环，也是工具层对外暴露的
唯一权威统计。整个 crate 无可变全局状态，函数纯输入到输出，便于在 mcp-server 里随取随用。

## 2. 对外 API 与导出文件清单

入口是无状态 struct `ReportGenerator`（方法皆静态），错误统一 `ReportError`：

```rust
impl ReportGenerator {
    pub fn generate_summary(result: &SearchResult) -> SearchResultSummary;
    pub fn export_tsv(result: &SearchResult, output_dir: &Path) -> Result<(), ReportError>;
    pub fn export_json(result: &SearchResult, output_path: &Path) -> Result<(), ReportError>;
    pub fn export_metadata(metadata: &RunMetadata, output_path: &Path) -> Result<(), ReportError>;
    // 以下渲染 self-contained HTML（实参类型见 search-engine / xic）
    pub fn render_annotation(annotation: &SpectrumAnnotation, output_path: &Path)
        -> Result<(), ReportError>;
    pub fn render_xic(xic_data: &XicData, output_path: &Path, plotly_mode: PlotlyMode)
        -> Result<(), ReportError>;
    pub fn render_unified(data: &UnifiedViewData, output_path: &Path, plotly_mode: PlotlyMode)
        -> Result<(), ReportError>;
    pub fn render_xic_3d(data: &Xic3dData, output_path: &Path, plotly_mode: PlotlyMode,
                         max_peaks_per_scan_3d: Option<usize>) -> Result<(), ReportError>;
}
```

导出文件清单（核对 `export.rs` 与 mcp-server 调用）：

- `export_tsv` 在 `output_dir` 下固定写 3 个文件：`psm.tsv`、`peptide.tsv`、`protein.tsv`。
- `export_json` 把整个 `SearchResult` 经 `serde_json::to_string_pretty` 写到调用方给定路径；
  mcp-server 取名 `result.json`。
- `export_metadata` 把 `RunMetadata` 写到给定路径；mcp-server 取名 `run_metadata.json`。

故 `export_results` 工具实际产出五件套：`psm.tsv + peptide.tsv + protein.tsv + result.json
+ run_metadata.json`。`ReportError` 含 `IoError { path, detail }`、`SerializationError(String)`、
`EmptyResult`、`EmptyMs2Window { scan }`、`AnnotationError { scan, detail }`，并
`impl From<ReportError> for CoreError`（映射为 `ValidationError`）。

四个 `render_*` 共享同一套手法：把数据用 `serde_json` 序列化后经 `escape_json_for_html`（把
`<`/`>` 转成 `\u003c`/`\u003e`，既避免 `</script>` 提前闭合与 HTML 注入，又保持结果仍是合法 JSON）
注入由 `include_str!` 内联进二进制的 HTML 模板占位符，从而产出零外链、可直接双击打开的单文件。
Plotly 图统一引 CDN `plotly-2.35.2`；`render_unified` 还把原始谱峰 `raw_scans` 与离子元数据一并下发，
供前端按 SILAC 在浏览器侧重算轻/重通道，`render_xic_3d` 则用 `max_peaks_per_scan_3d` 限制每张谱
的绘制峰数以控体积。

## 3. 摘要内容与 1% FDR 过滤

`generate_summary` 先做 FDR 过滤再统计；过滤规则：任一 PSM 带 `q_value` 则保留 `q_value <= 0.01`
的，否则全保留（无 FDR 时不丢数据）。基于过滤后的 PSM 计算：

- `total_spectra_searched`、`search_duration_sec`：透传自 `result.summary`（引擎侧）。
- `total_psms`：过滤前 `result.psms.len()`。
- `psms_at_1pct_fdr`：过滤后条数。
- `unique_peptides_at_1pct_fdr` / `protein_groups_at_1pct_fdr`：对过滤后 PSM 的
  `peptide_sequence` / `protein_accessions` 取 `HashSet` 去重计数。
- `identification_rate = psms_at_1pct_fdr / total_spectra_searched`（分母为 0 取 0.0）。
- `modification_distribution`（修饰名 -> 计数）、`charge_distribution`（电荷 -> 计数，仅 `charge > 0`）。
- `median_score` / `median_delta_mass_ppm`：仅取有限值，排序后 `compute_median`（偶数取中间两数均值）。

需强调：分布、唯一计数与两个中位数全部只基于过滤后的 PSM 集合，与 `total_psms`（过滤前规模）口径
不同；这样摘要既只反映可信鉴定，又保留了"提交了多少谱图、引擎给了多少候选"的原始基数，便于诊断
低鉴定率究竟出在采集、参数还是数据库。

日志：`tracing::info` 记 `id_rate`/`psms_1pct`；当 `identification_rate < 0.10` 额外 `warn`
提示参数或库可疑。摘要结构（`core::search_result::SearchResultSummary`）：

```rust
pub struct SearchResultSummary {
    total_spectra_searched: u64,
    total_psms: u64,
    psms_at_1pct_fdr: u64,
    unique_peptides_at_1pct_fdr: u64,
    protein_groups_at_1pct_fdr: u64,
    median_score: f64,
    median_delta_mass_ppm: f64,
    identification_rate: f64,
    modification_distribution: HashMap<String, u64>,
    charge_distribution: HashMap<i32, u64>,
    search_duration_sec: f64,
}
```

## 4. 简化源码

FDR 过滤骨架（`summary.rs`）：

```rust
let has_qvalues = result.psms.iter().any(|p| p.q_value.is_some());
let filtered_psms: Vec<_> = if has_qvalues {
    result.psms.iter()
        .filter(|p| p.q_value.is_some_and(|q| q <= 0.01))
        .collect()
} else {
    result.psms.iter().collect() // 无 q 值则全保留
};
```

写文件骨架与 TSV 表头（`export.rs`）：

```rust
fs::create_dir_all(output_dir)?;
let mut psm_file = create_file(&output_dir.join("psm.tsv"))?;
writeln!(psm_file,
    "scan\tsequence\tcharge\tprecursor_mz\tcalculated_mz\t\
     delta_ppm\tscore\tq_value\tproteins\tis_decoy\tmodifications")?;
for psm in &result.psms { writeln!(psm_file, "{}\t{}\t...", ...)?; }
// peptide.tsv 表头: sequence | proteins | best_score | q_value | psm_count
// protein.tsv 表头: accession | description | coverage | peptide_count | unique_peptide_count
```

字段细节：`q_value` 为 `None` 写 `NA`；多值列（`proteins`、`modifications`）用 `;` 连接；
`sanitize_tsv` 把 `\t \n \r ;` 替换为空格，避免破坏分隔。

数值精度固定写死：`precursor_mz`/`calculated_mz` 保留 8 位、`delta_ppm` 2 位、`score` 与肽级
`best_score` 6 位、蛋白 `coverage` 4 位；列顺序与表头亦固定，确保同一结果跨平台导出的 TSV 可逐字节
复现，方便做回归对比与下游脚本解析。

## 5. 调用链（mcp-server）

```text
generate_summary(run_id) -> get_result -> ReportGenerator::generate_summary
                          -> Json<SearchResultSummary>

export_results(run_id, output_dir)
   -> ReportGenerator::export_tsv(result, output_dir)               // psm/peptide/protein.tsv
   -> ReportGenerator::export_json(result, output_dir/result.json)
   -> ReportGenerator::export_metadata(result.metadata, output_dir/run_metadata.json)
   -> Json<ExportResultsOutput { output_dir, files: [5 个文件名] }>
```

两工具都经 `get_result(input.result, input.run_id)` 从直传参数或 `run_cache` 取回 `SearchResult`；
HTML 渲染则由 `annotate_spectrum` 等工具走 `render_unified` / `render_annotation`。

`ExportResultsOutput` 回带 `output_dir` 与 `files` 文件名清单，供上层确认产物落点；中途任一步的
`IoError` / `SerializationError` 都经 `From<ReportError>` 升为 `CoreError::ValidationError`，再由
mcp-server 转成带码、描述与修复建议的结构化错误，不会以裸字符串或 panic 形式外泄。

## 6. 测试入口

```text
cargo test -p protein-copilot-report --offline
```

单元测试覆盖 FDR 过滤、唯一肽计数、鉴定率、中位数、电荷分布、三 TSV 落盘、JSON 往返、元数据
`run_id` 一致、HTML 注入与转义（实测 31 passed）。

回到 [README](README.md)。
