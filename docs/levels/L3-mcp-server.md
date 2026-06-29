# L3 — MCP Server / 工具层

承接 [L2](L2-architecture.md)。本篇讲清最上层 `mcp-server`（workspace 的 MCP 服务 bin crate）如何把全部确定性库 crate 组装成 **27 个 MCP 工具**，经 stdio 暴露给 LLM 客户端（Copilot CLI / Claude Desktop），以及异步搜索模型、三类缓存、结构化 I/O 与历史持久化。

## 1. 职责与位置

`mcp-server` 是四层架构里的「MCP Tool 层」：向上对接 AI 编排层，向下调度 12 个库 crate。它**只做编排与桥接**——解析参数、委托库 crate、返回结构化 JSON 或错误。一切数值计算（打分、FDR、Δppm、蛋白推断）都在库里，工具本身既不算也**不调 LLM**。

源码只有三个文件：

- `src/main.rs`（62 行）：`#[tokio::main]` 入口；`tracing_subscriber` 初始化（读 `RUST_LOG`，默认 `info`；`PROTEIN_LOG_JSON=1` 切 JSON；日志写 **stderr**，stdout 让给 JSON-RPC）；`ProteinCopilotServer::new().serve(stdio()).await`，再 `service.waiting().await`。
- `src/tools.rs`（4551 行）：全部工具定义、输入/输出结构体、缓存与服务器结构体。
- `src/history.rs`（174 行）：搜索历史磁盘持久化。

bin 名 `protein-copilot-mcp`，传输 `rmcp::transport::stdio`。服务器状态全部内聚在一个结构体（依赖注入，无全局可变状态）：

```rust
pub struct ProteinCopilotServer {
    tool_router: ToolRouter<Self>,                 // rmcp 收集的 27 个工具
    registry: EngineRegistry,                      // 注册 SimpleSearch + Sage
    run_cache: Arc<Mutex<OrderedRunCache>>,        // 搜索 run 缓存（上限 100）
    dia_cache: Arc<Mutex<OrderedDiaCache>>,        // DIA 提取缓存（上限 10，可溢出磁盘）
    reader_cache: Arc<Mutex<LruCache<PathBuf, Arc<dyn SpectrumReader>>>>, // 索引读取器 LRU（8）
}
```

`new()` 里 `registry.register(SimpleSearchEngine)` 与 `SageAdapter`，并由 `Self::tool_router()` 收集工具；`#[rmcp::tool_router]` 挂在 impl、`#[rmcp::tool_handler]` 挂在 `ServerHandler` 上。

## 2. 工具清单（27 个，按类别分组）

每个工具基本都是 `Result<Json<T>, ErrorData>`（少数无失败路径的返回裸 `Json<T>`）。下面「名字 — 一句话 + 关键输入/输出」逐一列出（名字与 `description` 抽自源码 `#[rmcp::tool(...)]`）：

**读谱（spectrum-io）**
- `read_spectra` — 读 mgf/mzML 返回统计摘要，分析的第一步。入 `{file_path}` -> `SpectrumSummary`（谱图数、m/z 范围、RT 范围、电荷分布、中位峰数）。
- `get_spectrum` — 按 scan_number（1 起）读单张谱图。入 `{file_path, scan_number}` -> `Spectrum`（m/z 数组、强度数组、前体、MS level）。

**推参（param-recommend）**
- `recommend_params` — 据谱图特征推荐搜索参数。入 `{summary?|file_path?, hints?, database_path?}` -> `AiDecision<SearchParams>`（带 confidence + explanation）。
- `list_presets` — 列内置预设（standard/phospho/TMT/SILAC/open）。`()` -> `PresetsResponse{presets}`。
- `prepare_search` — 一步备搜：read_spectra + recommend_params + 解析 FASTA 三合一。入 `{input_files, hints?, organism?|database_path?, engine?}` -> `PrepareSearchOutput{params, reasoning, confidence, alternatives, evidence, spectra_summary, database_info?}`。

**搜索生命周期（search-engine）**
- `run_search` — 启动数据库搜索，**立即返回 run_id**（后台跑）。入 `{params?, input_files, database_path?, hints?, dia_run_id?}` -> `SearchStarted{run_id, status, message}`。
- `get_search_status` — 轮询进度。入 `{run_id}` -> `SearchProgress{status, stage, progress_pct, elapsed_sec, estimated_remaining_sec, error_category?, has_diagnostics}`。
- `cancel_search` — 取消运行中的搜索（abort 任务）。入 `{run_id}` -> `SearchProgress`（status=Cancelled）。
- `check_engine` — 列可用引擎与健康状态（裸 `Json`）。`()` -> `EngineStatus{engine, status, all_engines}`。
- `diagnose_search` — 搜索诊断，失败做错因分析、完成做质量评估。入 `{run_id}` -> `DiagnoseSearchOutput{overall_status, error_category?, stages, anomalies, suggestions, ...}`。

**结果（report）**
- `generate_summary` — 1% FDR 过滤后的统计摘要。入 `{result?|run_id?}` -> `SearchResultSummary`。
- `export_results` — 导出 psm/peptide/protein.tsv + result.json + run_metadata.json。入 `{result?|run_id?, output_dir}` -> `ExportResultsOutput{output_dir, files}`。
- `list_searches` — 列近期搜索（活动 + 历史，裸 `Json`）。入 `{status_filter?, limit?}` -> `ListSearchesResponse{searches}`。

**推断（protein-inference）**
- `infer_proteins` — parsimony + razor 肽分配 + 蛋白 FDR + 可选覆盖率。入 `{run_id?|result?, q_value_threshold=0.01, fasta_path?}` -> `InferenceResult`（蛋白组 + score + q-value + 肽分配）。

**注释/可视化（xic）**
- `annotate_spectrum` — 单谱 b/y 碎片离子标注，出交互 HTML。入 `{run_id?|file_path?, scan_number, peptide_sequence?, charge?, retention_time_min?, ...}` -> `AnnotateResult{output_path, matched_ions, total_ions, delta_mass_ppm, ...}`。两模式：PSM 上下文 / 手动。
- `extract_xic` — 提取 XIC，Plotly.js HTML（MS1 前体 + MS2 碎片色谱，支持 SILAC 重标，`view='3d'` 出 3D 总览 + 逐 scan b/y 标注）。入 `{run_id?|file_path?, scan_number, ...}` -> `ExtractXicResult{output_path, ms2_scan_count, light_trace_count, heavy_trace_count, has_ms1_xic, ...}`。

**DIA（dia-extraction）**
- `extract_dia_precursors` — 从 DIA 数据提取候选前体并缓存增强谱图供 run_search。入 `{file_path, output_mode=pseudo, min_charge?, max_charge?, acquisition_mode?}` -> `DiaExtractionOutput{detected_mode, ms1_count, ms2_count, total_precursors_extracted, avg_precursors_per_ms2, run_id, ...}`。
- `extract_spectrum_precursors` — 对单张 MS2 提取候选前体（关联最近 MS1 做同位素分析）。入 `{file_path, scan_number, min_charge?, max_charge?}` -> `SingleSpectrumExtractionResult`。
- `get_dia_cache_status` — 查 DIA 提取结果是否仍在缓存（run_search 前自检）。入 `{dia_run_id}` -> `DiaCacheStatusOutput{exists, location(memory/disk/not_found), spectrum_count?, extracted_at?}`。

**导入（result-import）**
- `import_search_results` — 导入 DIA-NN / 自定义 JSON / pFind 并匹配 mzML scan。入 `{result_file, format=auto, mzml_dir, unimod_path?, rt_tolerance_min=0.5, filter_qvalue=0.01, run_filter?}` -> `ImportResult`（返回 run_id，供后续注释/XIC/汇总）。

**数据库（fasta-db）**
- `list_databases` — 列内置 FASTA（Human/Mouse/E.coli/Yeast/Arabidopsis/cRAP）+ 下载状态。入 `{cache_dir?}` -> `ListDatabasesOutput{databases}`。
- `download_database` — 按 ID 从 UniProt HTTPS 下载并缓存。入 `{database_id, cache_dir?, force?}` -> `DownloadDatabaseResult`。
- `get_database_info` — 已下载库详情（蛋白数、大小、SHA256、下载日期、前 5 个 accession）。入 `{database_id, cache_dir?}` -> `DatabaseInfo`。

**entrapment（entrapment-analysis）**
- `classify_entrapment_hits` — 把 trap 命中按与靶蛋白组同源性分 L0-L4，`mzml_dir` 提供时追踪碎片来源。入 `{results_file, format?, config_file, target_fasta, output_dir?, mzml_dir?}` -> `ClassifyEntrapmentOutput{total/target/trap/ambiguous_psms, level_counts(l0..l4), top_razor_families}`；落 classified.tsv / razor_errors.tsv / run_metadata.json / entrapment_report.html。
- `analyze_entrapment_stats` — 从 classified.tsv 出分级分布 / 家族簇 / Δ质量统计。入 `{classified_file}` -> `AnalyzeEntrapmentStatsOutput{total_classified, level_distribution, delta_mass_stats, top_protein_families}`。
- `find_similar_targets` — 用编辑距离（同长 Hamming / 跨长 Levenshtein）找相似靶肽。入 `{peptide, target_fasta, max_mismatches=2}` -> `FindSimilarTargetsOutput{level, best_target_peptide?, mismatches?, delta_mass_da?, substitution_type?, edit_distance?, ...}`。
- `annotate_provenance` — 单谱碎片来源标注（trap / target / shared / unassigned），出 mirror plot HTML。入 `{file_path, scan_number, trap_sequence, target_sequence, modifications, fragment_tolerance_ppm=20, max_fragment_charge=2, chimera_threshold=0.3, output_path?}` -> `AnnotateProvenanceOutput{output_file, trap_matched_count, target_matched_count, shared_count, ...}`。

## 3. 异步搜索模型

`run_search` 是唯一长耗时工具，采用「立即返回 + 后台执行 + 轮询」：

```text
run_search(params, files)                        立即返回（不阻塞 client）
  | run_id = Uuid::new_v4()
  | run_cache.insert(run_id, RunState{ Running, .. })
  | tokio::spawn(async move { ... }) ----------------+ 后台任务
  '-> SearchStarted{ run_id, "Running", message }    |
                                                      |  engine.search(...).await
get_search_status(run_id)  --轮询-->  SearchProgress  |  on_progress 回调写 run_cache
  | status == "Running"    -> progress_pct / stage    |
  | status == "Completed"  -> generate_summary/export |
  | status == "Failed: .." / "Cancelled"              v
cancel_search(run_id) -> handle.abort(); 置 "Cancelled"   原子更新 result + 落历史
```

要点：

- **后台任务**：`tokio::spawn` 内先装 `PanicGuard`——任务异常退出（panic/abort）时若状态仍是 `Running`，其 `Drop` 改写为 `"Failed: task panicked"`，杜绝永久卡 Running。
- **进度回调**：`ProgressCallback` 在每阶段写 `run_cache`，仅当状态仍 `Running` 才覆盖 stage / 进度。
- **原子收尾**：搜索结束在**单次加锁**内同时写 progress + result + diagnostics（异常检测 finalize），再在锁外 `save_entry` 落历史；若已被 `Cancelled` 则不覆盖。
- **DIA 旁路**：带 `dia_run_id` 时，先校验 params / engine / 数据库存在，**再** `dia_cache.remove()` 取走缓存谱图（避免校验失败白白消耗缓存），走 `engine.search_with_spectra(...)`。
- **取消**：`cancel_search` 调 `handle.abort()` 并置 `Cancelled`，且仅当当前为 `Running`。

三类缓存全部 `Arc<Mutex<...>>`，并用**锁中毒自愈** `.lock().unwrap_or_else(|e| e.into_inner())`（全文 15 处），保证一次 panic 持锁不会让 server 整体卡死：

```text
run_cache    OrderedRunCache  上限 100，FIFO 只淘汰非 Running 的 run
dia_cache    OrderedDiaCache  上限 10，溢出 spill 到 .proteincopilot/dia_cache/<id>.bin（bincode）
reader_cache LruCache         上限 8，缓存 IndexedMzMLReader 供 O(1) scan 查询
```

`get_result(direct, run_id)` 统一从「直接传入的 `SearchResult`」或「`run_id` 取缓存」二选一解析；两者都给报歧义，run 未完成报当前状态，被淘汰报「max 100 recent runs」。

## 4. 工具调用伪代码

每个工具薄而同构：解析（schemars 生成 JSON Schema）-> 校验 -> 委托库 -> 结构化结果/错误。

```rust
#[rmcp::tool(
    name = "read_spectra",
    description = "Read a mass spectrometry file (mgf/mzML) and return a statistical summary ..."
)]
fn read_spectra(
    &self,
    Parameters(input): Parameters<ReadSpectraInput>,
) -> Result<Json<SpectrumSummary>, ErrorData> {
    let _span = tracing::info_span!("mcp_tool", name = "read_spectra").entered();
    tracing::info!(file = %input.file_path, "started");
    validate_file_path(&input.file_path)?;            // 空 / 不存在 -> INVALID_PARAMS
    let reader = self.get_or_create_reader(path)?;    // 委托 spectrum-io（带缓存）
    // ... 聚合统计 ...
    Ok(Json(summary))
}
```

输入结构体一律 `#[derive(Deserialize, schemars::JsonSchema)]`，字段用 `///` 注释——它直接成为工具参数说明：

```rust
#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct GetSpectrumInput {
    /// Path to the spectrum file (.mgf or .mzML)
    file_path: String,
    /// Scan number to retrieve (1-based)
    scan_number: u32,
}
```

整体调用链：

```text
MCP Client --JSON-RPC--> rmcp tool_router（派发到具名工具）
  -> 反序列化 Parameters<XxxInput>（部分字段容错：既收 JSON 对象也收 "JSON 字符串"）
  -> validate_file_path / validate_scan_number（前置校验，失败即 INVALID_PARAMS）
  -> 委托库 crate（spectrum-io / param-recommend / search-engine / xic / ...）
  -> Ok(Json<XxxOutput>)  或  Err(mcp_err(ErrorCode::INVALID_PARAMS, msg))
```

错误统一经两个 helper 构造结构化 `ErrorData`（带 `ErrorCode`），库代码不 `unwrap`/`expect`、工具层不 `panic`：

```rust
fn mcp_err(code: ErrorCode, err: impl std::fmt::Display) -> ErrorData { /* ... */ }
fn mcp_core_err(err: CoreError) -> ErrorData { /* 库错误 -> MCP 错误 */ }
```

## 5. 设计约束

- **不在 server 调 LLM**：工具只委托确定性库；意图理解、推参理由、结果解释留给上层 agent/prompt。
- **数值不交给 LLM**：FDR、打分、Δppm、推断全部在 Rust 库里完成。
- **无全局可变状态**：状态内聚在 `ProteinCopilotServer`（`Arc<Mutex>` 注入），无 `static mut`。
- **结构化 I/O**：每个输入/输出皆 `serde` + `schemars::JsonSchema`；因 rmcp 要求 outputSchema 根类型为 `object`，输出一律用具名结构体（如 `AnalyzeEntrapmentStatsOutput`）而非 `Json<serde_json::Value>`。
- **结构化错误**：`ErrorData{ code, message }`，给出码 + 描述 + 修复建议（如「Use list_databases ... or download_database」）。
- **可观测**：每个工具进入即开 `tracing::info_span!("mcp_tool", name = ...)` 并 `info!("started"/"completed")`（共 27 处 span），关键字段（文件、run_id、引擎、容差、命中数）结构化记录。

## 6. 搜索历史持久化（history.rs）

每次搜索结束（Completed / Failed / Cancelled）落一份 JSON 摘要到 `~/.protein-copilot/history/<run_id>.json`，**只存摘要统计、不存完整 PSM 列表**：

```rust
pub struct SearchHistoryEntry {
    run_id: Uuid, status: String, created_at: DateTime<Utc>, elapsed_sec: f64,
    engine_info: EngineInfo, input_files: Vec<PathBuf>, params_used: SearchParams,
    total_psms: Option<u64>, psms_at_1pct_fdr: Option<u64>,
    identification_rate: Option<f64>, protein_groups: Option<u64>,
}
```

- `save_entry()`：写 `<run_id>.json`（pretty JSON），随后 `evict_oldest()` 按 FIFO 把超过 `MAX_HISTORY = 500` 的最旧文件删掉。
- `load_all()`：读目录全部 `.json`，按 `created_at` 倒序返回；坏文件只 `warn` 跳过。
- `list_searches` 把内存里的活动 run 与磁盘历史合并返回，所以重启 server 后历史仍在。

落盘失败、目录不可写、文件损坏均只 `tracing::warn!` 而不致命——历史是辅助信息，绝不阻断主搜索流程。

---

往上看整体分层与数据流见 [L2](L2-architecture.md)；返回文档目录见 [README](README.md)。
