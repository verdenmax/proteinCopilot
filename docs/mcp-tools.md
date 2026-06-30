# MCP 工具参考（ProteinCopilot）

> 本文件由 `scripts/gen_mcp_tools_doc.py` 从 release 二进制的 `tools/list`（JSON Schema）**自动生成**，是面向调用方的权威接口契约——无需阅读源码。工具签名变更后请重跑该脚本。

> 工具总数：**27**。传输：JSON-RPC 2.0 over stdio。生成时间：2026-06-30。


## 启动与接入

```bash
# 直接运行已编译的二进制（推荐发布形态）
./protein-copilot-mcp

# 或从源码运行
cargo run --release -p protein-copilot-mcp-server
```

命令行自检（无需客户端，直接在终端查看工具契约）：

```bash
./protein-copilot-mcp --list-tools          # 文本目录：参数/类型/范围/默认/输出
./protein-copilot-mcp --list-tools --json   # 完整 JSON Schema（机器可读）
./protein-copilot-mcp --help                # 用法
```

在 MCP 客户端（Copilot CLI / Claude Desktop 等）中登记：

```json
{
  "mcpServers": {
    "protein-copilot": {
      "command": "/path/to/protein-copilot-mcp",
      "env": { "RUST_LOG": "info" }
    }
  }
}
```

- 所有工具的输入/输出均为结构化 JSON，类型见下方「参数」表与「共享数据类型」。
- 描述文本为二进制 `#[schemars]` 原文（即客户端实际收到的内容），故为英文。
- 生成的 HTML/TSV 默认写入 `./output/`，可用 `PROTEIN_OUTPUT_DIR` 改基目录；返回路径为绝对路径。
- 搜索为异步：`run_search` 立即返回 `run_id`，用 `get_search_status` 轮询，完成后 `generate_summary` / `export_results` / `infer_proteins`。


## 工具索引

- **读取谱图**：[`read_spectra`](#read_spectra), [`get_spectrum`](#get_spectrum)
- **参数推荐**：[`recommend_params`](#recommend_params), [`list_presets`](#list_presets), [`prepare_search`](#prepare_search)
- **搜索生命周期**：[`run_search`](#run_search), [`get_search_status`](#get_search_status), [`cancel_search`](#cancel_search), [`check_engine`](#check_engine), [`diagnose_search`](#diagnose_search)
- **结果摘要与导出**：[`generate_summary`](#generate_summary), [`export_results`](#export_results), [`list_searches`](#list_searches)
- **蛋白推断**：[`infer_proteins`](#infer_proteins)
- **谱图注释与可视化**：[`annotate_spectrum`](#annotate_spectrum), [`extract_xic`](#extract_xic)
- **DIA 数据提取**：[`extract_dia_precursors`](#extract_dia_precursors), [`extract_spectrum_precursors`](#extract_spectrum_precursors), [`get_dia_cache_status`](#get_dia_cache_status)
- **外部结果导入**：[`import_search_results`](#import_search_results)
- **FASTA 数据库**：[`list_databases`](#list_databases), [`download_database`](#download_database), [`get_database_info`](#get_database_info)
- **entrapment 分析**：[`classify_entrapment_hits`](#classify_entrapment_hits), [`analyze_entrapment_stats`](#analyze_entrapment_stats), [`find_similar_targets`](#find_similar_targets), [`annotate_provenance`](#annotate_provenance)

---

## 工具详情


### 读取谱图


#### `read_spectra`

Read a mass spectrometry file (mgf/mzML) and return a statistical summary including spectrum count, m/z range, RT range, charge distribution, and median peaks per spectrum. Use this as the first step to understand input data.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `file_path` | 是 | string | — | Path to the spectrum file (.mgf or .mzML) |

**输出**：`SpectrumSummary` — Statistical summary of a spectrum file.

This is the primary input for AI-driven parameter recommendation
and data quality assessment. The LLM reads this summary (via MCP tool)
to understand data characteristics before making recommendations.


#### `get_spectrum`

Read a single spectrum from a file by scan number (1-based). Returns the spectrum with m/z array, intensity array, precursor info, and MS level.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `file_path` | 是 | string | — | Path to the spectrum file (.mgf or .mzML) |
| `scan_number` | 是 | integer | — | Scan number to retrieve (1-based) |

**输出**：`Spectrum` — A single mass spectrum with peak data.

`mz_array` and `intensity_array` must always have the same length.
`mz_array` is expected to be sorted in ascending order.

The `precursors` field supports both DDA (typically 1 precursor)
and DIA (0 precursors, or 1 with a wide isolation window).


### 参数推荐


#### `recommend_params`

Recommend search parameters based on spectrum file characteristics. Input: SpectrumSummary from read_spectra + optional UserHints (experiment_type, instrument_type, enzyme). Output: recommended SearchParams with confidence score and explanation. Note: set database_path in params to the FASTA file path.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `database_path` | 否 | string *(可空)* | — | FASTA database path. If provided, sets database_path in the recommended params. |
| `file_path` | 否 | string *(可空)* | — | Path to spectrum file. Used to generate summary if summary is not provided. |
| `hints` | 否 | [`UserHints`](#类型-userhints) *(可空)* | — | Optional user hints (experiment_type, instrument_type, enzyme) |
| `summary` | 否 | [`SpectrumSummary`](#类型-spectrumsummary) *(可空)* | — | Spectrum summary (from read_spectra). If provided, uses this directly. |

**输出**：`AiDecision` — Structured wrapper for all AI-assisted decisions.

Every time the LLM (via Agent/Skill) makes a recommendation or
interpretation, the result is wrapped in this struct. This ensures:

- **Auditability**: `explanation` + `evidence` record *why* a decision
  was made.
- **Calibration**: `confidence` gives a quantitative self-assessment.
- **Exploration**: `alternatives` lists other options considered.
- **Context**: `input_summary` captures the data the decision was based on.

# Example JSON output (per copilot-instructions.md §2.4)

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


#### `list_presets`

List all built-in search parameter presets (standard, phospho, TMT, SILAC, open search). Each preset includes name, description, parameters, and applicable scenarios.


*无参数。*

**输出**：`PresetsResponse` — Presets list response


#### `prepare_search`

One-shot search preparation: reads spectrum files, recommends search parameters, and resolves a FASTA database. Combines read_spectra + recommend_params + download_database into a single call. Provide either 'database_path' (direct FASTA path) or 'organism' (e.g. 'human', 'mouse', 'E.coli', '小鼠') for auto-resolution. Returns ready-to-use SearchParams that can be passed directly to run_search.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `input_files` | 是 | array&lt;string&gt; | — | Paths to spectrum files (.mgf or .mzML) |
| `cache_dir` | 否 | string *(可空)* | — | Override cache directory for database downloads. |
| `database_path` | 否 | string *(可空)* | — | Direct FASTA database path. Takes priority over organism auto-resolution. |
| `engine` | 否 | string *(可空)* | — | Search engine: "Sage" or "SimpleSearch". Default: "SimpleSearch". |
| `hints` | 否 | [`UserHints`](#类型-userhints) *(可空)* | — | Optional user hints (experiment_type, instrument_type, enzyme) |
| `organism` | 否 | string *(可空)* | — | Target organism for auto database resolution (e.g. "human", "mouse", "E.coli", "小鼠"). |

**输出**：`PrepareSearchOutput`


### 搜索生命周期


#### `run_search`

Run a proteomics database search. Returns immediately with a run_id. The search runs in the background. Call get_search_status(run_id) to check progress. When status is Completed, use generate_summary(run_id) and export_results(run_id).


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `input_files` | 是 | array&lt;string&gt; | — | Paths to spectrum files |
| `database_path` | 否 | string *(可空)* | — | FASTA database path (used if params is not provided or params.database_path is placeholder) |
| `dia_run_id` | 否 | string *(可空)* | — | Optional run_id from extract_dia_precursors. When provided, uses cached DIA-extracted spectra instead of reading from input_files. |
| `hints` | 否 | [`UserHints`](#类型-userhints) *(可空)* | — | Optional user hints for auto-recommendation (used when params not provided) |
| `params` | 否 | [`SearchParams`](#类型-searchparams) *(可空)* | — | Search parameters (from recommend_params decision). If not provided, auto-recommends. |

**输出**：`SearchStarted` — Response when search is started asynchronously


#### `get_search_status`

Check the status of a search started by run_search. Returns progress percentage and elapsed time. When status is Completed, use generate_summary(run_id) to get results.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `run_id` | 是 | string | — | Run ID from run_search |

**输出**：`SearchProgress` — Progress information for a running search.


#### `cancel_search`

Cancel a running search. The search task is immediately terminated and status is set to Cancelled.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `run_id` | 是 | string | — | Run ID of the search to cancel. |

**输出**：`SearchProgress` — Progress information for a running search.


#### `check_engine`

Check available search engines and their health status. Returns engine name, version, supported features, and availability.


*无参数。*

**输出**：`EngineStatus` — Engine status response


#### `diagnose_search`

Get diagnostic report for a search run. Works for both failed searches (error analysis) and completed searches (quality assessment). Returns stage metrics, detected anomalies, and repair suggestions. Call after get_search_status shows the search has finished (status is Completed, Failed, or Cancelled). Use has_diagnostics=true from get_search_status to confirm diagnostics are available.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `run_id` | 是 | string | — | The run_id to diagnose (from run_search or get_search_status) |

**输出**：`DiagnoseSearchOutput`


### 结果摘要与导出


#### `generate_summary`

Generate a statistical summary from search results with 1% FDR filtering. Includes identification rate, median score, median delta ppm, modification and charge distributions. Use this after run_search to interpret results.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `result` | 否 | [`SearchResult`](#类型-searchresult) *(可空)* | — | Search result to summarize (provide either this or run_id) |
| `run_id` | 否 | string *(可空)* | — | Run ID from a previous run_search call (server retrieves cached result) |

**输出**：`SearchResultSummary` — Statistical summary of search results for LLM-driven interpretation.

This is the primary input for the AI layer to understand and explain
search quality. All fields are deterministically computed by Rust.


#### `export_results`

Export search results to files. Creates psm.tsv, peptide.tsv, protein.tsv, result.json, and run_metadata.json in the specified output directory.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `output_dir` | 是 | string | — | Output directory path |
| `result` | 否 | [`SearchResult`](#类型-searchresult) *(可空)* | — | Search result to export (provide either this or run_id) |
| `run_id` | 否 | string *(可空)* | — | Run ID from a previous run_search call (server retrieves cached result) |

**输出**：`ExportResultsOutput`


#### `list_searches`

List recent search runs with their status, duration, and key metrics. Includes both active searches and completed history.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `limit` | 否 | integer *(可空)* | — | Maximum results to return. Default 20. |
| `status_filter` | 否 | string *(可空)* | — | Filter by status prefix (e.g. "Completed", "Failed"). Optional. |

**输出**：`ListSearchesResponse`


### 蛋白推断


#### `infer_proteins`

Run protein inference on search results. Performs parsimony analysis, razor peptide assignment, protein-level FDR, and optional sequence coverage. Input: run_id from a previous search or direct SearchResult. Returns protein groups with scores, q-values, and peptide assignments.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `fasta_path` | 否 | string *(可空)* | — | Path to FASTA database for sequence coverage calculation. If not provided, coverage is not calculated. |
| `q_value_threshold` | 否 | number *(可空)* | `0.01` | Q-value threshold for filtering PSMs before inference (default: 0.01). |
| `result` | 否 | [`SearchResult`](#类型-searchresult) *(可空)* | — | Direct SearchResult (alternative to run_id). |
| `run_id` | 否 | string *(可空)* | — | Run ID from a previous search. Uses cached PSMs for inference. |

**输出**：`InferenceResult` — Complete result of protein inference.

Produced by running the parsimony algorithm, razor peptide assignment,
and protein-level FDR on a set of PSMs.


### 谱图注释与可视化


#### `annotate_spectrum`

Annotate a single spectrum with peptide fragment ion matching. Generates an interactive HTML file showing matched b/y ions. Two modes: (1) provide run_id + scan_number to annotate an existing PSM, or (2) provide file_path + scan_number + peptide_sequence + charge for manual annotation. In mode 2, you can set scan_number=0 and provide retention_time_min to auto-find the matching scan.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `scan_number` | 是 | integer | — | Scan number (1-based) to annotate. |
| `charge` | 否 | integer *(可空)* | — | Charge state — required for manual mode. |
| `extraction_tolerance` | 否 | [`MassTolerance`](#类型-masstolerance) *(可空)* | — | m/z extraction tolerance for XIC (default: 20 ppm). |
| `file_path` | 否 | string *(可空)* | — | Spectrum file path — use for manual annotation mode. |
| `fragment_tolerance` | 否 | [`MassTolerance`](#类型-masstolerance) *(可空)* | — | Fragment mass tolerance. Default: 20 ppm. |
| `label_type` | 否 | [`LabelType`](#类型-labeltype) *(可空)* | — | Heavy-label type for SILAC comparison. |
| `n_cycles` | 否 | integer *(可空)* | — | Number of DIA cycles before/after target for XIC (default: 5). |
| `output_path` | 否 | string *(可空)* | — | Output HTML file path. Default: ./annotation_scan{N}.html |
| `peptide_sequence` | 否 | string *(可空)* | — | Peptide sequence — required for manual mode. |
| `plotly_mode` | 否 | [`PlotlyMode`](#类型-plotlymode) *(可空)* | — | Plotly loading mode (default: Cdn). |
| `protein_accessions` | 否 | array *(可空)* | — | Protein accession(s) — optional for manual mode (e.g. ["sp\|P00001\|TEST_HUMAN"]). |
| `retention_time_min` | 否 | number *(可空)* | — | Retention time in minutes — alternative to scan_number for auto scan lookup. |
| `run_id` | 否 | string *(可空)* | — | Run ID — use to annotate an existing PSM from a search result. |
| `top_n_ions` | 否 | integer *(可空)* | — | Number of top fragment ions for XIC (default: all, zero-intensity traces excluded). |

**输出**：`AnnotateResult`


#### `extract_xic`

Extract XIC (Extracted Ion Chromatogram) for a peptide from an mzML file. Generates an interactive HTML file with Plotly.js showing MS1 precursor and MS2 fragment ion chromatograms. Supports SILAC heavy-label comparison. Two modes: (1) provide run_id + scan_number to use PSM context, or (2) provide file_path + scan_number + peptide_sequence + charge + precursor_mz. In mode 2, set scan_number=0 with retention_time_min to auto-find scan. Set view='3d' for a 3D MS2 overview (RT x m/z x intensity sticks) plus per-scan b/y annotated spectra with total peak counts (output: xic3d_scan{N}.html).


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `scan_number` | 是 | integer | — | Scan number (1-based) to center the XIC around. |
| `charge` | 否 | integer *(可空)* | — | Precursor charge state. |
| `extraction_tolerance` | 否 | [`MassTolerance`](#类型-masstolerance) *(可空)* | — | Mass tolerance for XIC peak extraction. Default: 20 ppm. |
| `file_path` | 否 | string *(可空)* | — | Path to the spectrum file (.mzML). XIC extraction requires mzML format for MS1+MS2 and isolation window data. |
| `intensity_rule` | 否 | [`IntensityRule`](#类型-intensityrule) *(可空)* | — | How to extract intensity from peaks within tolerance. Default: MaxInWindow. |
| `label_type` | 否 | [`LabelType`](#类型-labeltype) *(可空)* | — | Heavy-label configuration. Use {"Silac": {"heavy_k_delta": 8.014199, "heavy_r_delta": 10.008269}} for standard SILAC. |
| `max_peaks_per_scan_3d` | 否 | integer *(可空)* | — | Only for view=3d: max non-matched peaks per scan drawn in the 3D overview (display declutter; matched b/y always kept). Default 200. |
| `modifications` | 否 | array *(可空)* | — | Modifications applied to this peptide (fixed + variable). If omitted, uses unmodified sequence. |
| `n_cycles` | 否 | integer *(可空)* | — | Number of DIA cycles before and after target scan. Default: 5. |
| `output_path` | 否 | string *(可空)* | — | Output HTML file path. Default: ./output/xic_scan{N}.html |
| `peptide_sequence` | 否 | string *(可空)* | — | Peptide amino acid sequence (one-letter codes). |
| `plotly_mode` | 否 | [`PlotlyMode`](#类型-plotlymode) *(可空)* | — | Plotly.js loading: 'Cdn' (default, smaller) or 'Embedded' (offline). |
| `precursor_mz` | 否 | number *(可空)* | — | True precursor m/z. For DIA data, use the PSM-derived value, not the isolation window center. |
| `retention_time_min` | 否 | number *(可空)* | — | Retention time in minutes. When scan_number is 0, auto-finds the closest MS2 scan matching this RT and precursor_mz. |
| `run_id` | 否 | string *(可空)* | — | Run ID from a previous search. Auto-fills peptide, charge, mods, precursor_mz. MVP: single-file searches only. |
| `top_n_ions` | 否 | integer *(可空)* | — | Number of top fragment ions to display. Default: all (zero-intensity excluded). |
| `view` | 否 | [`XicView`](#类型-xicview) *(可空)* | — | View mode: 'standard' (default, 2D XIC line chart) or '3d' (3D MS2 overview RT x m/z x intensity + per-scan b/y annotated spectra with total peak counts). |

**输出**：`ExtractXicResult` — Result returned by `extract_xic`.


### DIA 数据提取


#### `extract_dia_precursors`

Extract candidate precursor ions from DIA mass spectrometry data. Reads mzML file, detects DIA mode from isolation window widths, extracts precursor candidates from MS1 isotope patterns, and caches enhanced spectra for use with run_search. Returns extraction statistics.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `file_path` | 是 | string | — | Path to the spectrum file (.mzML) |
| `acquisition_mode` | 否 | string *(可空)* | — | Override acquisition mode detection: "DDA" or "DIA". If not set, auto-detects. |
| `max_charge` | 否 | integer *(可空)* | — | Maximum charge state to consider (default: 5) |
| `min_charge` | 否 | integer *(可空)* | — | Minimum charge state to consider (default: 2) |
| `output_mode` | 否 | string | `"pseudo"` | Output mode: "multi" (multiple precursors per spectrum) or "pseudo" (one precursor per spectrum). Default: "pseudo" |

**输出**：`DiaExtractionOutput`


#### `extract_spectrum_precursors`

Extract candidate precursor ions from a single MS2 spectrum. Reads the mzML file, finds the target MS2 by scan number, correlates it to the nearest MS1, and runs isotope pattern analysis within the isolation window. Returns extracted precursor candidates with charge states and the correlation method used.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `file_path` | 是 | string | — | Path to the spectrum file (.mzML). The file is read to find both the target MS2 scan and nearby MS1 spectra for isotope pattern extraction. |
| `scan_number` | 是 | integer | — | Scan number (1-based) of the MS2 spectrum to extract precursors for. |
| `max_charge` | 否 | integer *(可空)* | — | Maximum charge state to consider (default: 5) |
| `min_charge` | 否 | integer *(可空)* | — | Minimum charge state to consider (default: 2) |

**输出**：`SingleSpectrumExtractionResult` — Result of extracting precursors from a single MS2 spectrum.


#### `get_dia_cache_status`

Check if a DIA extraction result is still cached and available for use with run_search. Returns cache location (memory/disk/not_found) and spectrum count. Call this before run_search(dia_run_id=...) to verify the extraction hasn't been evicted.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `dia_run_id` | 是 | string | — | The dia_run_id returned by extract_dia_precursors |

**输出**：`DiaCacheStatusOutput`


### 外部结果导入


#### `import_search_results`

Import external search results (DIA-NN, custom JSON, pFind) and match to mzML scans. Returns a run_id for use with annotate_spectrum, extract_xic, and generate_summary.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `mzml_dir` | 是 | string | — | Directory containing mzML files. File association: raw_name + '.mzML'. |
| `result_file` | 是 | string | — | Path to external search result file (.json, .parquet, .spectra, .tsv). |
| `filter_qvalue` | 否 | number | `0.01` | Q-value threshold for filtering (DIA-NN). Default: 0.01. |
| `format` | 否 | string | `"auto"` | Result file format. 'auto' detects from extension. Options: auto, custom_json, diann_parquet, pfind_spectra, pfind_tsv. |
| `rt_tolerance_min` | 否 | number | `0.5` | RT tolerance in minutes for scan matching. Default: 0.5. |
| `run_filter` | 否 | string *(可空)* | — | Optional: only import PSMs from this specific run/raw_title. |
| `unimod_path` | 否 | string *(可空)* | — | Path to unimod.xml. If not provided, uses builtin modification database (~22 common mods). |

**输出**：`ImportResult` — Result of the import operation.


### FASTA 数据库


#### `list_databases`

List all built-in FASTA protein databases (Human, Mouse, E.coli, Yeast, Arabidopsis, cRAP contaminants) with download status. Shows which databases are cached locally and which are available for download.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `cache_dir` | 否 | string *(可空)* | — | Override cache directory. Default: .proteincopilot/databases/ |

**输出**：`ListDatabasesOutput`


#### `download_database`

Download a FASTA protein database by ID (e.g. 'human_swissprot', 'mouse_swissprot', 'ecoli_swissprot', 'yeast_swissprot', 'arabidopsis_swissprot', 'crap'). Downloads from UniProt via HTTPS and caches locally. Returns the local file path for use as database_path in search parameters. Use list_databases first to see available options.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `database_id` | 是 | string | — | Database ID (e.g. "human_swissprot", "mouse_swissprot", "crap") |
| `cache_dir` | 否 | string *(可空)* | — | Override cache directory. Default: .proteincopilot/databases/ |
| `force` | 否 | boolean *(可空)* | — | Force re-download even if already cached |

**输出**：`DownloadDatabaseResult` — Result returned after downloading a database.


#### `get_database_info`

Get detailed information about a downloaded FASTA database: protein count, file size, SHA256 hash, download date, and first 5 protein accessions. The database must be downloaded first using download_database.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `database_id` | 是 | string | — | Database ID to get info for |
| `cache_dir` | 否 | string *(可空)* | — | Override cache directory. Default: .proteincopilot/databases/ |

**输出**：`DatabaseInfo` — Information about a cached database.


### entrapment 分析


#### `classify_entrapment_hits`

Classify trap-database PSM hits by homology to target proteome. Reads search results, applies target/trap rules from YAML config, digests target FASTA, and classifies each trap PSM as L0-L4. Optionally traces fragment ion provenance when mzml_dir is provided. Outputs classified.tsv, razor_errors.tsv, run_metadata.json, and entrapment_report.html. Returns summary statistics.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `config_file` | 是 | string | — | Path to YAML config file defining target/trap rules |
| `results_file` | 是 | string | — | Path to search results file (.parquet for DIA-NN or .tsv) |
| `target_fasta` | 是 | string | — | Path to target FASTA database |
| `format` | 否 | string *(可空)* | — | Result format override. Auto-detects from extension if omitted. |
| `mzml_dir` | 否 | string *(可空)* | — | Directory containing mzML files for provenance tracing (optional) |
| `output_dir` | 否 | string *(可空)* | — | Output directory (default: ./output/entrapment/) |

**输出**：`ClassifyEntrapmentOutput`


#### `analyze_entrapment_stats`

Get detailed statistics from a classified entrapment TSV file. Returns level distribution, protein family clusters, and delta-mass analysis. Use after classify_entrapment_hits to interpret results.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `classified_file` | 是 | string | — | Path to classified TSV file (output from classify_entrapment_hits) |

**输出**：`AnalyzeEntrapmentStatsOutput`


#### `find_similar_targets`

Find similar target peptides for a given sequence. Digests the target FASTA, compares the query peptide against target peptides using edit distance (Hamming for same-length, Levenshtein for cross-length). Returns closest matches with mass differences and substitution type annotations. Useful for investigating individual trap PSMs.


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `peptide` | 是 | string | — | Peptide sequence to look up |
| `target_fasta` | 是 | string | — | Path to target FASTA database |
| `max_mismatches` | 否 | integer *(可空)* | — | Maximum mismatches to consider (default: 2) |

**输出**：`FindSimilarTargetsOutput`


#### `annotate_provenance`

Annotate a single spectrum with fragment ion provenance analysis. Generates a mirror plot HTML file showing which peaks come from the trap peptide, target peptide, both (shared), or neither (unassigned).


| 参数 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `file_path` | 是 | string | — | Path to the mzML spectrum file |
| `scan_number` | 是 | integer | — | Scan number (1-based) |
| `trap_sequence` | 是 | string | — | Trap peptide sequence (stripped, no modifications) |
| `chimera_threshold` | 否 | number | `0.3` | Chimera threshold for shared_ratio (default: 0.3) |
| `fragment_tolerance_ppm` | 否 | number | `20.0` | Fragment mass tolerance in ppm (default: 20.0) |
| `max_fragment_charge` | 否 | integer | `2` | Maximum fragment charge state (default: 2) |
| `modifications` | 否 | array&lt;array&lt;object&gt;&gt; | `[]` | Modifications as (position, delta_mass) pairs |
| `output_path` | 否 | string *(可空)* | — | Output HTML file path (default: ./provenance_scan{N}.html) |
| `target_sequence` | 否 | string | `""` | Target peptide sequence (stripped). Empty string if L4. |

**输出**：`AnnotateProvenanceOutput`


---

## 共享数据类型


参数与输出中引用的复合类型定义如下（枚举列出全部取值，结构体列出字段）。


### 类型 AcquisitionMode

Data acquisition mode.


变体：

- `DDA` — Data-Dependent Acquisition — narrow isolation window, single precursor.
- `DIA` — Data-Independent Acquisition — wide isolation window, multiple co-fragmented precursors.
- `Unknown` — Acquisition mode could not be determined.

### 类型 DecoyStrategy

Target-decoy strategy for FDR estimation.


变体：

- `Reverse` — Reverse protein sequences.
- `Shuffle` — Shuffle protein sequences.
- `None` — No decoy database.

### 类型 EngineInfo

Static metadata about a search engine.


| 字段 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `name` | 是 | string | — | Engine name (e.g. "pFind", "MSFragger", "Comet"). |
| `supported_features` | 是 | array&lt;string&gt; | — | Features supported by this engine (e.g. "open_search", "glyco"). |
| `version` | 是 | string | — | Engine version string (e.g. "3.1.0"). |

### 类型 Enzyme

Digestion enzyme used for protein cleavage.


变体：

- `Trypsin` — Trypsin — cleaves after K/R (not before P).
- `LysC` — Lys-C — cleaves after K.
- `GluC` — Glu-C — cleaves after D/E.
- `AspN` — Asp-N — cleaves before D.
- `Chymotrypsin` — Chymotrypsin — cleaves after F/W/Y/L.
- `TrypsinP` — Trypsin/P — cleaves after K/R (including before P).
- `NonSpecific` — No specific cleavage rule.
- `Custom` — User-defined enzyme with custom cleavage rule.：{ `cleavage_rule`: string; `name`: string }

### 类型 IntensityRule

Intensity extraction strategy.


变体：

- `MaxInWindow` — Highest peak within tolerance window (default).
- `SumInWindow` — Sum of all peaks within tolerance window.
- `NearestPeak` — Nearest peak to theoretical m/z.

### 类型 LabelType

Heavy-label type for SILAC or custom isotope labeling.


变体：

- `Silac` — SILAC heavy amino acids.：{ `heavy_k_delta`: number; `heavy_r_delta`: number }
- `Custom` — Custom residue mass shifts.：{ `residue_deltas`: array&lt;array&lt;object&gt;&gt; }

### 类型 MassTolerance

Mass tolerance specification (value + unit).


| 字段 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `unit` | 是 | [`ToleranceUnit`](#类型-toleranceunit) | — | Tolerance unit. |
| `value` | 是 | number | — | Tolerance value (must be positive). |

### 类型 ModPosition

Position where a modification can occur.


变体：

- `Anywhere` — Anywhere on the peptide.
- `AnyNTerm` — Any peptide N-terminus.
- `AnyCTerm` — Any peptide C-terminus.
- `ProteinNTerm` — Protein N-terminus only.
- `ProteinCTerm` — Protein C-terminus only.

### 类型 Modification

A chemical modification (fixed or variable) applied during search.


| 字段 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `mass_delta` | 是 | number | — | Mass shift in Daltons. Positive for mass-increasing modifications (e.g., +57.021 for Carbamidomethyl), negative for mass-decreasing (rare; e.g., -18.011 for dehydration). |
| `name` | 是 | string | — | Modification name (e.g. "Carbamidomethyl", "Oxidation"). |
| `position` | 是 | [`ModPosition`](#类型-modposition) | — | Where on the peptide/protein this modification can occur. |
| `residues` | 是 | array&lt;string&gt; | — | Target residues (e.g. `['C']` for Carbamidomethyl, `['M']` for Oxidation). |

### 类型 PeptideResult

Peptide-level search result, aggregated from PSMs.


| 字段 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `best_score` | 是 | number | — | Best score among all PSMs for this peptide. |
| `protein_accessions` | 是 | array&lt;string&gt; | — | Protein accessions containing this peptide. |
| `psm_count` | 是 | integer | — | Number of PSMs supporting this peptide. |
| `q_value` | 否 | number *(可空)* | — | q-value at peptide level (`None` if not calculated). |
| `sequence` | 是 | string | — | Peptide amino acid sequence. |

### 类型 PlotlyMode

Plotly.js loading mode for HTML output.


变体：

- `Cdn` — Load from CDN (default, smaller file).
- `Embedded` — Embed plotly-basic.min.js inline (larger file, works offline).

### 类型 ProteinResult

Protein-level search result, aggregated from peptides.


| 字段 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `accession` | 是 | string | — | Protein accession (e.g. UniProt ID like "P12345"). |
| `coverage` | 是 | number | — | Sequence coverage (0.0–1.0). |
| `description` | 是 | string | — | Protein description / name. |
| `peptide_count` | 是 | integer | — | Total number of peptides mapped to this protein. |
| `unique_peptide_count` | 是 | integer | — | Number of unique (non-shared) peptides. |

### 类型 Psm

A single Peptide-Spectrum Match — the fundamental unit of a database search.

Each PSM represents the assignment of a peptide sequence to a specific
MS2 spectrum, along with scoring and quality metrics.


| 字段 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `calculated_mz` | 是 | number | — | Theoretical precursor m/z calculated from the peptide. |
| `charge` | 是 | integer | — | Precursor charge state. |
| `delta_mass_ppm` | 是 | number | — | Mass deviation between observed and calculated (ppm). Formula: `(precursor_mz - calculated_mz) / calculated_mz × 1e6` Provided by the search engine adapter. |
| `extra` | 否 | object *(可空)* | — | Engine-specific extra fields (e.g., Sage's matched_peaks, delta_next). Preserves information that doesn't fit the standard Psm fields. |
| `is_decoy` | 是 | boolean | — | Whether this PSM is from the decoy database. |
| `modifications` | 是 | array&lt;[`Modification`](#类型-modification)&gt; | — | Modifications identified on this peptide. |
| `peptide_sequence` | 是 | string | — | Identified peptide sequence (one-letter amino acid codes). |
| `precursor_mz` | 是 | number | — | Observed precursor m/z (mass-to-charge ratio). |
| `protein_accessions` | 是 | array&lt;string&gt; | — | Protein accessions this peptide maps to. |
| `q_value` | 否 | number *(可空)* | — | False discovery rate q-value (`None` if FDR not yet calculated). |
| `score` | 是 | number | — | Search engine score (higher = better match, engine-dependent). |
| `spectrum_scan` | 是 | integer | — | Scan number of the matched spectrum (1-based). |

### 类型 RunMetadata

Metadata for a single analysis run.

Generated at the start of each search and updated as the run progresses.
Stored alongside results to enable full provenance tracking.


| 字段 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `created_at` | 是 | string | — | Timestamp when the run was created. |
| `duration_sec` | 否 | number *(可空)* | — | Total duration in seconds (`None` if not yet completed). |
| `engine_info` | 是 | [`EngineInfo`](#类型-engineinfo) | — | Search engine that executed this run. |
| `input_files` | 是 | array&lt;string&gt; | — | Input spectrum file paths. |
| `params_used` | 是 | [`SearchParams`](#类型-searchparams) | — | Search parameters used for this run. |
| `run_id` | 是 | string | — | Unique identifier for this run (auto-generated UUIDv4). |
| `status` | 是 | [`RunStatus`](#类型-runstatus) | — | Current status of the run. |

### 类型 RunStatus

Status of an analysis run.


变体：

- `Pending` — Run is queued but not yet started.
- `Running` — Run is currently executing.
- `Completed` — Run completed successfully.
- `Failed` — Run failed with an error.：{ `reason`: string }

### 类型 SearchParams

Complete search configuration for a proteomics database search.

Use [`SearchParams::validate()`] after construction or deserialization
to ensure all values are within acceptable ranges.


| 字段 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `acquisition_mode` | 否 | [`AcquisitionMode`](#类型-acquisitionmode) *(可空)* | — | Data acquisition mode. `None` = auto-detect or not applicable. |
| `database_path` | 是 | string | — | Path to the FASTA protein database file. |
| `decoy_strategy` | 是 | [`DecoyStrategy`](#类型-decoystrategy) | — | Target-decoy strategy for FDR estimation. |
| `engine` | 否 | string *(可空)* | — | Search engine to use. Valid values: `"Sage"` (production, sage-core library), `"SimpleSearch"` (built-in MVP engine). Case-insensitive. Default: `"SimpleSearch"`. |
| `enzyme` | 是 | [`Enzyme`](#类型-enzyme) | — | Digestion enzyme. |
| `fixed_modifications` | 是 | array&lt;[`Modification`](#类型-modification)&gt; | — | Modifications always present (e.g. Carbamidomethyl on C). |
| `fragment_tolerance` | 是 | [`MassTolerance`](#类型-masstolerance) | — | Fragment ion mass tolerance. |
| `max_peptide_length` | 否 | integer | `50` | Maximum peptide length in residues (default: 50). Peptides longer than this are excluded from search results. |
| `max_variable_modifications` | 否 | integer | `3` | Maximum number of variable modifications per peptide (default: 3). Limits combinatorial explosion during variable modification enumeration. |
| `min_peptide_length` | 否 | integer | `7` | Minimum peptide length in residues (default: 7). Peptides shorter than this are excluded from search results. |
| `missed_cleavages` | 是 | integer | — | Maximum number of missed cleavage sites (0–5). |
| `precursor_tolerance` | 是 | [`MassTolerance`](#类型-masstolerance) | — | Precursor ion mass tolerance. |
| `variable_modifications` | 是 | array&lt;[`Modification`](#类型-modification)&gt; | — | Modifications that may or may not be present (e.g. Oxidation on M). |

### 类型 SearchResult

Complete search result from one engine run.

This is the top-level output of a search engine invocation, combining
all PSMs, peptide/protein rollups, statistical summary, and metadata.
The `run_id` links this result to the originating run context.


| 字段 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `engine_info` | 是 | [`EngineInfo`](#类型-engineinfo) | — | Search engine metadata (name, version, supported features). |
| `metadata` | 是 | [`RunMetadata`](#类型-runmetadata) | — | Full run metadata for provenance tracking. |
| `params_used` | 是 | [`SearchParams`](#类型-searchparams) | — | Search parameters used for this run. |
| `peptides` | 是 | array&lt;[`PeptideResult`](#类型-peptideresult)&gt; | — | Peptide-level aggregation. |
| `proteins` | 是 | array&lt;[`ProteinResult`](#类型-proteinresult)&gt; | — | Protein-level aggregation. |
| `psms` | 是 | array&lt;[`Psm`](#类型-psm)&gt; | — | All PSMs returned by the search. |
| `run_id` | 是 | string | — | Unique identifier for this analysis run. |
| `summary` | 是 | [`SearchResultSummary`](#类型-searchresultsummary) | — | Statistical summary (engine-side, before FDR filtering).  For a summary with proper FDR filtering, use `ReportGenerator::generate_summary()` from the `report` crate. |

### 类型 SearchResultSummary

Statistical summary of search results for LLM-driven interpretation.

This is the primary input for the AI layer to understand and explain
search quality. All fields are deterministically computed by Rust.


| 字段 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `charge_distribution` | 是 | object | — | Charge state distribution (charge → count). |
| `identification_rate` | 是 | number | — | Identification rate: `psms_at_1pct_fdr / total_spectra_searched`. |
| `median_delta_mass_ppm` | 是 | number | — | Median mass deviation (ppm) across all PSMs. |
| `median_score` | 是 | number | — | Median search engine score across all PSMs. |
| `modification_distribution` | 是 | object | — | Modification frequency distribution (modification name → count). |
| `protein_groups_at_1pct_fdr` | 是 | integer | — | Protein groups at 1% FDR. |
| `psms_at_1pct_fdr` | 是 | integer | — | PSMs passing 1% FDR threshold. |
| `search_duration_sec` | 是 | number | — | Total search duration in seconds. |
| `total_psms` | 是 | integer | — | Total number of PSMs returned by the engine (before FDR filtering). |
| `total_spectra_searched` | 是 | integer | — | Total number of spectra submitted to the search. |
| `unique_peptides_at_1pct_fdr` | 是 | integer | — | Unique peptide sequences at 1% FDR. |

### 类型 SpectrumFormat

Supported spectrum file formats.


变体：

- `MzML` — mzML format (PSI standard).
- `Mgf` — Mascot Generic Format.
- `Pfb` — pXtract3 / pParse2+ binary PFB format.

### 类型 SpectrumSummary

Statistical summary of a spectrum file.

This is the primary input for AI-driven parameter recommendation
and data quality assessment. The LLM reads this summary (via MCP tool)
to understand data characteristics before making recommendations.


| 字段 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `file_path` | 是 | string | — | Source file path. |
| `format` | 是 | [`SpectrumFormat`](#类型-spectrumformat) | — | File format. |
| `median_isolation_window_da` | 否 | number *(可空)* | — | Median isolation window width in Da (`None` if no isolation windows found). Useful for DIA detection: DDA windows are typically < 3 Da, DIA > 5 Da. |
| `median_peaks_per_spectrum` | 是 | integer | — | Median number of peaks per spectrum. |
| `ms1_count` | 是 | integer | — | Number of MS1 spectra. |
| `ms2_count` | 是 | integer | — | Number of MS2 spectra. |
| `mz_range` | 是 | array&lt;number&gt; | — | m/z range: \[min, max\]. |
| `precursor_charge_distribution` | 是 | object | — | Distribution of precursor charge states (charge → count). |
| `rt_range_min` | 是 | array&lt;number&gt; | — | Retention time range: \[min, max\] in minutes. |
| `total_spectra` | 是 | integer | — | Total number of spectra in the file. |

### 类型 ToleranceUnit

Unit for mass tolerance.


变体：

- `Ppm` — Parts per million.
- `Da` — Daltons (absolute mass).

### 类型 UserHints

Structured hints from the user (via LLM translation).

All fields are optional — the rule engine has sensible defaults.
When provided, hints override or adjust the automatic recommendation.


| 字段 | 必填 | 类型 | 默认 | 说明 |
|------|------|------|------|------|
| `custom_notes` | 否 | string *(可空)* | — | Free-form notes from the user (e.g., "use 5ppm tolerance"). |
| `enzyme` | 否 | [`Enzyme`](#类型-enzyme) *(可空)* | — | Digestion enzyme override. If provided, replaces the preset default. |
| `experiment_type` | 否 | string *(可空)* | — | Experiment type (e.g., "phosphorylation", "TMT", "SILAC", "standard"). |
| `instrument_type` | 否 | string *(可空)* | — | Instrument type (e.g., "Orbitrap", "TOF", "QExactive"). |

### 类型 XicView

View mode for the `extract_xic` tool.


变体：

- `standard` — Standard 2D XIC line chart (default).
- `3d` — 3D MS2 overview + per-scan b/y annotated spectra.
