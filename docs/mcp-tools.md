# MCP Tools 参考

ProteinCopilot MCP Server 提供 16 个工具，通过 JSON-RPC over stdio 暴露给 LLM。

> **当前搜索引擎**：MVP 使用内置的 **SimpleSearchEngine**（基于 b/y 离子匹配的简化搜索），
> 后续将接入 pFind 作为生产级搜索引擎。SimpleSearch 足以验证完整流程，但搜索质量和性能不如专业引擎。

## 启动方式

```bash
cargo run --release -p protein-copilot-mcp-server
```

或通过 `.mcp.json` 配置自动发现。

---

## read_spectra

读取质谱文件，返回统计摘要。

**输入**：
```json
{
  "file_path": "/data/sample.mgf"
}
```

**输出**：`SpectrumSummary`

**使用场景**：分析数据特征，为参数推荐提供输入。

---

## get_spectrum

按 scan 号读取单张谱图。

**输入**：
```json
{
  "file_path": "/data/sample.mgf",
  "scan_number": 42
}
```

**输出**：`Spectrum`（含 mz_array、intensity_array、precursors）

---

## recommend_params

基于谱图特征推荐搜索参数。可直接传文件路径 + 数据库路径。

**输入**：
```json
{
  "file_path": "/data/sample.mgf",
  "database_path": "/data/human.fasta",
  "hints": {
    "experiment_type": "phosphorylation"
  }
}
```

**输出**：`AiDecision<SearchParams>`（含 decision、confidence、explanation）

> **注意**：`database_path` 会自动注入到推荐结果的 `decision` 中，LLM 无需手动修改。

---

## list_presets

列出所有内置搜索参数预设。

**输入**：无

**输出**：`{ "presets": [...] }`（5 个预设）

---

## run_search

执行蛋白质数据库搜索。支持两种调用方式：

**简单模式（LLM 推荐）**— 自动推荐参数：
```json
{
  "input_files": ["/data/sample.mgf"],
  "database_path": "/data/human.fasta"
}
```

**高级模式** — 传入 recommend_params 返回的参数：
```json
{
  "params": { ... },
  "input_files": ["/data/sample.mgf"]
}
```

**DIA 模式** — 使用 `extract_dia_precursors` 提取的结果：
```json
{
  "input_files": ["/data/sample.mzML"],
  "database_path": "/data/human.fasta",
  "dia_run_id": "a1b2c3d4-..."
}
```

> 当提供 `dia_run_id` 时，搜索引擎从 DIA 缓存中获取已提取前体的谱图，通过 `search_with_spectra()` 执行搜索。

结果自动缓存，后续用 `run_id` 引用。

**输出**：`{run_id, status: "Running", message: "..."}`（立即返回，搜索在后台执行）

> **异步模式**：run_search 不阻塞。搜索完成后用 `get_search_status(run_id)` 确认，再用 `generate_summary(run_id)` 获取结果。

---

## get_search_status

查询搜索进度。搜索完成后 status 变为 "Completed"。

**输入**：
```json
{
  "run_id": "9f71e493-..."
}
```

**输出**：`SearchProgress`
```json
{
  "run_id": "...",
  "status": "Completed",
  "progress_pct": 1.0,
  "elapsed_sec": 5.2
}
```

status 值：`Running` / `Completed` / `Failed: <reason>` / `Cancelled`

---

## cancel_search

取消正在运行的搜索任务。

**输入**：
```json
{
  "run_id": "9f71e493-..."
}
```

**输出**：`SearchProgress`
```json
{
  "run_id": "...",
  "status": "Cancelled",
  "progress_pct": 0.45,
  "elapsed_sec": 3.1
}
```

> 内部通过 `JoinHandle::abort()` 终止搜索任务。已完成或已失败的搜索无法取消。

---

## check_engine

检查搜索引擎状态。

**输入**：无

**输出**：`{ "engine": {...}, "status": "Healthy" }`

---

## generate_summary

从搜索结果生成 FDR 过滤后的统计摘要。支持通过 `run_id` 引用缓存结果。

**输入**：
```json
{
  "run_id": "7ab6d7d4-df4d-4aa0-..."
}
```

或直接传 `{"result": {...}}`

**输出**：`SearchResultSummary`

---

## export_results

将搜索结果导出为文件。支持通过 `run_id` 引用缓存结果。

**输入**：
```json
{
  "run_id": "7ab6d7d4-df4d-4aa0-...",
  "output_dir": "./output"
}
```

**输出**：导出文件列表

---

## list_searches

列出搜索历史，包括当前活跃的搜索和持久化到磁盘的历史记录（`~/.protein-copilot/history/`）。

**输入**：
```json
{
  "status_filter": "Completed",
  "limit": 10
}
```

所有参数可选。`status_filter` 按状态过滤（`Running` / `Completed` / `Failed` / `Cancelled`），`limit` 限制返回数量。

**输出**：`Vec<SearchHistoryEntry>`
```json
[
  {
    "run_id": "...",
    "status": "Completed",
    "input_files": ["/data/sample.mgf"],
    "started_at": "2026-04-15T10:30:00Z",
    "elapsed_sec": 12.5
  }
]
```

---

## annotate_spectrum

对单张谱图进行 b/y 离子匹配注释，生成交互式 HTML 可视化。

**渲染模式自动选择**：
- **mzML + DIA（fragment XIC >1 点）**：统一视图（标注 + XIC + SILAC 交互控件）
- **mzML + DDA（fragment XIC ≤1 点）**：统一视图（标注 only，跳过无意义的 XIC）
- **mgf 等非 mzML**：纯标注视图（SVG 谱图 + 覆盖图）

> DDA/DIA 判断基于 fragment XIC 数据点数量，不依赖隔离窗口宽度阈值，
> 因此窄窗 DIA（如 Scanning SWATH 2 Da）也能正确显示 XIC。

**统一视图包含**：
- 📄 源文件名 + Scan/RT 合并显示
- Fragment Ion Coverage（SVG bracket 标注）
- 谱图（SVG，b/y 离子着色 + hover tooltip）
- SILAC 预设切换（None / Standard / Medium / Custom）
- 逐离子 L/H 开关网格（Precursor, y₁, b₂... 竖向排列）
- MS1 Precursor XIC + MS2 Fragment Ion XIC（Plotly.js）

**模式一：基于搜索结果**
```json
{
  "run_id": "9f71e493-...",
  "scan_number": 42
}
```

**模式二：直接指定肽段**
```json
{
  "file_path": "/data/sample.mzML",
  "peptide_sequence": "PEPTIDEK",
  "charge": 2,
  "scan_number": 42
}
```

**可选 XIC 参数**（仅 mzML 有效）：
```json
{
  "n_cycles": 5,
  "top_n_ions": 6,
  "label_type": { "Silac": { "k_delta": 8.014199, "r_delta": 10.008269 } },
  "extraction_tolerance": { "value": 20, "unit": "Ppm" },
  "plotly_mode": "Cdn"
}
```

**输出**：`AnnotateResult`
```json
{
  "score": 0.85,
  "matched_ions": 12,
  "total_ions": 18,
  "output_path": "./output/annotation_scan42.html",
  "message": "...Includes XIC + SILAC controls (DIA)."
}
```

---

## extract_dia_precursors

从 DIA mzML 文件中提取候选前体离子。自动检测 DDA/DIA 采集模式，对 DIA 数据执行 MS1 同位素模式检测和 MS1↔MS2 关联。

**输入**：
```json
{
  "file_path": "/data/dia_sample.mzML"
}
```

**输出**：`ExtractionResult`
```json
{
  "dia_run_id": "a1b2c3d4-...",
  "acquisition_mode": "DIA",
  "total_ms1": 500,
  "total_ms2": 5000,
  "precursors_extracted": 12345,
  "charge_distribution": {"2": 6000, "3": 4500, "4": 1500, "5": 345}
}
```

> 提取结果缓存在 `OrderedDiaCache` 中，通过 `dia_run_id` 传给 `run_search` 使用。

---

## extract_spectrum_precursors

对单张 MS2 谱图提取候选母离子。读取 mzML 文件，找到目标 scan，关联最近的 MS1 谱图，在隔离窗口内运行同位素模式分析。

**输入**：
```json
{
  "file_path": "/data/dia_sample.mzML",
  "scan_number": 42
}
```

**输出**：`SingleSpectrumExtractionResult`
```json
{
  "ms2_scan": 42,
  "ms1_scan_used": 41,
  "correlation_method": "scan_order",
  "isolation_window": {"target_mz": 500.0, "lower_offset": 12.5, "upper_offset": 12.5},
  "precursors": [
    {"mz": 499.75, "charge": 2, "intensity": 1234.5},
    {"mz": 503.28, "charge": 3, "intensity": 890.0}
  ]
}
```

> 用于调试和检查单张谱图的母离子提取结果，了解 MS1 关联方法和同位素模式匹配详情。

---

## extract_xic

从 mzML 文件提取肽段的 XIC（Extracted Ion Chromatogram）。生成交互式 HTML 文件，展示 MS1 母离子和 MS2 碎片离子色谱图。支持 SILAC 重标比较。

两种模式：(1) 提供 `run_id` + `scan_number` 使用已有 PSM 上下文；(2) 提供 `file_path` + `scan_number` + `peptide_sequence` + `charge` + `precursor_mz` 手动指定。

**输入**：
```json
{
  "run_id": "abc-123",
  "scan_number": 42,
  "label_type": { "Silac": { "heavy_k_delta": 8.014199, "heavy_r_delta": 10.008269 } },
  "top_n_ions": 6
}
```

**输出**：生成 HTML 文件路径 + XIC 数据摘要。

> 用于验证 PSM 的色谱峰形，检查轻标/重标 SILAC 对的共洗脱行为。

---

## import_search_results

导入外部搜索结果（DIA-NN、pFind、自定义 JSON），匹配到 mzML 谱图。返回 `run_id` 供后续 `annotate_spectrum`、`extract_xic`、`generate_summary` 使用。

**输入**：
```json
{
  "result_file": "/data/diann_report.parquet",
  "mzml_dir": "/data/mzml/",
  "format": "auto",
  "filter_qvalue": 0.01
}
```

**输出**：`ImportResult`
```json
{
  "run_id": "def-456",
  "total_psms_imported": 12345,
  "runs_found": ["sample_01", "sample_02"],
  "scan_match_rate": 0.95
}
```

> 用于将 DIA-NN 或 pFind 的搜索结果导入 ProteinCopilot 进行可视化和进一步分析。

---

## LLM 完整工作流

**最简模式（一步搜索）：**
```
① run_search(input_files, database_path)       →  {run_id, "Running"}
② get_search_status(run_id)                    →  轮询直到 "Completed"
③ generate_summary(run_id)                     →  统计摘要
④ export_results(run_id, output_dir)           →  导出文件
```

**标准模式（分步控制）：**
```
① read_spectra(file_path)                      →  数据摘要
② recommend_params(file_path, database_path)   →  推荐参数
③ run_search(decision, input_files)            →  {run_id, "Running"}
④ get_search_status(run_id)                    →  轮询直到 "Completed"
⑤ generate_summary(run_id)                     →  统计摘要
⑥ export_results(run_id, output_dir)           →  导出文件
```

LLM 全程只需传简单参数（路径、run_id），无需构造复杂 JSON。
搜索不会阻塞 — LLM 可以在等待期间与用户交互。

**DIA 模式（两步搜索）：**
```
① extract_dia_precursors(file_path)                →  {dia_run_id, summary}
② run_search(input_files, database_path, dia_run_id) →  {run_id, "Running"}
③ get_search_status(run_id)                        →  轮询直到 "Completed"
④ generate_summary(run_id)                         →  统计摘要
⑤ export_results(run_id, output_dir)               →  导出文件
```

**单谱图母离子检查：**
```
① extract_spectrum_precursors(file_path, scan_number) →  {precursors, ms1_used, method}
② get_spectrum(file_path, scan_number)                →  查看原始谱图数据
```
