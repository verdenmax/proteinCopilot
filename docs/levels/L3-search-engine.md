# L3 - 搜索引擎调度与匹配

承接 [L2](L2-architecture.md)。本篇深入 `crates/search-engine`（package `protein-copilot-search-engine`）：它如何把 `SearchParams` 与谱图文件，经"酶切 -> 候选 -> 理论谱 -> 容差匹配打分 -> target-decoy -> 聚合"变成统一的 `SearchResult`，以及如何用 `SearchEngineAdapter` 抽象多引擎。

## 1. 职责与流程位置

一句话：**把搜索参数 + 谱图确定性地变成 PSM / 肽 / 蛋白三级结果**。上游是 `param-recommend` 给出的 `SearchParams`；下游是 `fdr`（q-value）、`report`（汇总导出）、`protein-inference`（parsimony 推断）。本 crate 内置 `SimpleSearch`（自研 MVP 引擎），并以库方式集成 `Sage`（sage-core），`pFind` 留有预留 adapter。

```text
SearchParams + input_files
  -> EngineRegistry.get("SimpleSearch" | "Sage")
  -> SearchEngineAdapter::search()
       SimpleSearch: FASTA -> digest -> b/y 打分 -> FDR -> 聚合
       Sage:         sage-core -> rayon 并行打分 -> LDA 重打分
  -> SearchResult (PSM/肽/蛋白 + summary + metadata)
```

## 2. 模块边界

| 文件 | 职责 |
|------|------|
| `lib.rs` | 公共 API 入口；re-export `EngineRegistry / SimpleSearchEngine / SageAdapter / SearchEngineError` |
| `chemistry.rs` | 残基单同位素质量 `residue_mass`、`peptide_mass`、`peptide_mz`，常量 `PROTON_MASS` / `WATER_MASS` |
| `digest.rs` | 体外酶切 `digest` / `digest_with_length`，私有 `find_cleavage_sites` / `split_at_sites`；类型 `DigestedPeptide` |
| `varmod.rs` | 可变修饰位点发现与组合枚举 `find_applicable_sites` / `enumerate_combinations`；末端哨兵 `NTERM_POS` / `CTERM_POS` |
| `matching.rs` | 前体容差匹配 + b/y 理论谱 + 计数打分 `match_spectrum` / `match_spectrum_all` |
| `fasta.rs` | FASTA 解析 `parse_fasta` -> `Vec<FastaEntry>` |
| `simple_engine.rs` | `SimpleSearchEngine`：编排全流程，实现 `SearchEngineAdapter` |
| `annotate.rs` | 单谱注释（可视化 / 质检）`annotate_spectrum` -> `SpectrumAnnotation` |
| `registry.rs` | `EngineRegistry`：按名注册 / 查找 adapter |
| `adapters/sage/` | Sage adapter：`mod.rs` 编排、`config.rs` 参数映射、`convert.rs` 谱图 / 容差转换 |
| `adapters/pfind.rs` | pFind adapter（stub，`SshConfig` 预留远程执行） |
| `error.rs` | `SearchEngineError`（thiserror），`From` 到 `CoreError` |

## 3. 关键数据结构

**SearchParams**（`core::search_params`）-- 一次搜索的完整配置：

```rust
pub struct SearchParams {
    pub enzyme: Enzyme,                       // Trypsin/LysC/GluC/AspN/Chymotrypsin/TrypsinP/NonSpecific/Custom
    pub missed_cleavages: u32,                // <= 5
    pub fixed_modifications: Vec<Modification>,
    pub variable_modifications: Vec<Modification>,
    pub precursor_tolerance: MassTolerance,   // value: f64 + unit: Ppm|Da
    pub fragment_tolerance: MassTolerance,
    pub database_path: String,                // FASTA 路径
    pub decoy_strategy: DecoyStrategy,        // Reverse|Shuffle|None
    pub acquisition_mode: Option<AcquisitionMode>,
    pub max_variable_modifications: u32,      // 默认 3
    pub min_peptide_length: u32,              // 默认 7
    pub max_peptide_length: u32,              // 默认 50
    pub engine: Option<String>,               // "Sage" | "SimpleSearch"
}
```

`Modification { name, mass_delta: f64, residues: Vec<char>, position: ModPosition }`；`ModPosition` 含 `Anywhere / AnyNTerm / AnyCTerm / ProteinNTerm / ProteinCTerm`。`validate()` 检查路径非空、两个容差有限且 > 0、漏切 <= 5、所有 `mass_delta` 有限。

**DigestedPeptide**（`digest.rs`）-- 候选肽：

```rust
pub struct DigestedPeptide {
    pub sequence: String,
    pub protein_accession: String,
    pub neutral_mass: f64,        // 单同位素中性质量 (Da)
    pub is_protein_nterm: bool,   // 是否位于蛋白 N 端
    pub is_protein_cterm: bool,   // 是否位于蛋白 C 端
}
```

**SearchResult**（`core::search_result`）-- 三级结果 + 摘要 + 元数据：

```rust
pub struct SearchResult {
    pub run_id: Uuid,
    pub engine_info: EngineInfo,
    pub params_used: SearchParams,
    pub psms: Vec<Psm>,               // 谱-肽匹配
    pub peptides: Vec<PeptideResult>, // 肽级聚合 (best_score/psm_count)
    pub proteins: Vec<ProteinResult>, // 蛋白级聚合 (含 coverage)
    pub summary: SearchResultSummary,
    pub metadata: RunMetadata,
}
```

`Psm` 携 `spectrum_scan / peptide_sequence / modifications / charge / precursor_mz / calculated_mz / delta_mass_ppm / score / q_value: Option<f64> / protein_accessions / is_decoy`。

## 4. 主流程伪代码（SimpleSearch）

`SimpleSearchEngine::run_search`（`simple_engine.rs`）：

```text
params.validate(); 若 input_files 为空 -> NoInputSpectra
proteins = parse_fasta(database_path)
for protein in proteins:
    all_peptides += digest_with_length(seq, acc, enzyme, mc, min_len, max_len)
若 decoy_strategy != None:
    decoys = fdr::generate_decoys(targets, strategy)   // 注 REV_/SHUF_ 前缀
    for d in decoys: all_peptides += digest_with_length(...)
预扫 input_files 统计 ms2 总数 (create_indexed_reader.read_summary)
for file in input_files:
    reader.for_each_spectrum(file, |spec| {              // 流式, 不全量载入
        if spec.ms_level == MS2:
            collect_psms_for_spectrum(spec, params, all_peptides, &mut psms)
    })
finalize_search_result: FDR -> 去 decoy -> 聚合 -> summary
```

酶切规则（`find_cleavage_sites`，位置为字符索引，UTF-8 安全）：

```text
Trypsin      : K/R 之后切, 后随 P 不切
TrypsinP     : K/R 之后切 (含 P)
LysC         : K 之后切
GluC         : D/E 之后切
AspN         : D 之前切
Chymotrypsin : F/W/Y/L 之后切
NonSpecific  : 枚举 [min_len, max_len] 内所有子串
```

漏切窗口：把相邻片段按 `windows(mc + 1)` 合并，长度落在 `[min, max]` 且全为标准残基者入候选，并记录是否触达蛋白 N/C 端（`is_protein_nterm/cterm`）。

单谱匹配（`match_spectrum`，DDA 取首个 precursor；`collect_psms_for_spectrum` 检测到多 precursor 或 isolation window > 5.0 Da 的 DIA 时改走 `match_spectrum_all` 遍历全部 precursor）：

```text
observed_mz = precursor.mz
charges = precursor.charge.map([c]).unwrap_or([2,3,1,4])
for pep in candidates:
    fixed_delta = apply_fixed_mods(seq, fixed_mods, nterm, cterm)
    sites  = find_applicable_sites(seq, var_mods, nterm, cterm)
    combos = enumerate_combinations(var_mods, sites, max_var_mods)   // 含空组合
    for combo in combos, for z in charges:                           // z<=0 跳过
        mz = peptide_mz(neutral + fixed_delta + combo.mass_delta, z)
        if within_tolerance(observed_mz, mz, precursor_tol):
            (nonpos, pos_deltas) = resolve_combined_mods(...)
            b = generate_b_ions_positional(...)   // z>=3 时额外出 2+ 碎片
            y = generate_y_ions_positional(...)
            matched = count_matched_ions(b ++ y, spec.mz_array, fragment_tol) // 二分查找
            score   = matched / (b.len + y.len)
            按 score 严格变大更新 best_match
```

`within_tolerance`：ppm 用 `|obs-theo|/theo*1e6 <= value`，Da 用 `|obs-theo| <= value`。打分确定性、可复现：`best_match` 仅在 `score` 严格更大时更新；聚合阶段 `aggregate_peptides` 按 `sequence`、`aggregate_proteins` 按 `accession` 排序，规避 HashMap 迭代序。

target-decoy：`is_decoy` 由 `core::util::is_decoy_accession`（前缀 `REV_/SHUF_/DECOY_/REVERSED_`）判定；`finalize_search_result` 把 `ScoredPsm{index, score, is_decoy}` 交 `fdr::calculate_fdr` 算 q-value，再 `retain(|p| !p.is_decoy)` 剔除诱饵。`build_summary` 在有 q-value 时按 `q <= 0.01` 统计 1% FDR 下的 PSM / 肽 / 蛋白数。

每条命中经 `build_psm` 落成 `Psm`：收集命中的固定修饰与 `applied_variable_mods`，记录 `spectrum_scan / charge / precursor_mz`（观测）/ `calculated_mz`（理论）/ `delta_mass_ppm`，并据诱饵前缀置 `is_decoy`。蛋白级 `aggregate_proteins` 还会在蛋白序列上滑动定位每条肽的全部出现区间，累计被覆盖位点算出 `coverage`（序列覆盖率）。

除批量搜索外，`annotate.rs::annotate_spectrum` 复用同一套 b/y 生成与 `within_tolerance` 逻辑，对单张谱图产出 `SpectrumAnnotation`：每个实验峰的可选 `IonAnnotation`、b/y 理论离子（`TheoreticalIon`）的命中状态、以及 `matched_ions / total_ions`，供 `report` / `xic` 做可视化与质检；DIA + SILAC 场景另有 `annotate_heavy_spectrum` 处理落在另一隔离窗口的重标碎片。

全流程是同步核心（`run_search` / `run_search_on_spectra`）由 async `search` 薄包装，谱图读取借 `for_each_spectrum` 流式处理、边读边匹配以压低内存峰值；各阶段经 `on_progress` 回调上报 `SearchProgress`（`stage` + `progress_pct` + `elapsed_sec`），典型刻度为读库 0.02、酶切 0.08、生成诱饵 0.10、匹配 0.15~0.90、算 FDR 0.88、聚合 0.92，供 MCP 客户端轮询展示进度与耗时。每次运行生成独立 `run_id` 并落 `RunMetadata`（参数、输入文件、引擎信息、状态、耗时），以保证可复现、可审计。

## 5. 引擎抽象

`core::engine::SearchEngineAdapter`（`#[async_trait]`，`Send + Sync`）：

```rust
async fn search(&self, params: &SearchParams, input_files: &[PathBuf],
                on_progress: ProgressCallback, diagnostics: &mut SearchDiagnostics)
    -> Result<SearchResult, CoreError>;
fn engine_info(&self) -> EngineInfo;                           // name/version/supported_features
async fn health_check(&self) -> Result<HealthStatus, CoreError>;
async fn search_with_spectra(...) -> ...;   // 默认报错; 引擎选择性实现 (DIA 缓存预载)
async fn cancel(&self, run_id: Uuid) -> Result<(), CoreError>; // 默认 no-op
```

`EngineRegistry`（`registry.rs`）以 `HashMap<String, Box<dyn SearchEngineAdapter>>` 按 `engine_info().name` 注册，提供 `register / get / list_available / health_check_all`。MCP 工具层据 `params.engine`（大小写不敏感）选 `Sage` 或 `SimpleSearch`。两引擎都实现 `search_with_spectra`，使 DIA 提取后的谱图缓存可直接复用、免去二次落盘读取。

- **SimpleSearch**：`engine_info` = `{ "SimpleSearch", "0.1.0", ["basic_search","b_y_scoring"] }`，纯进程内、无外部依赖，算法见第 4 节。
- **Sage**：`{ "Sage", "0.15.0", ["open_search","lfq","tmt","chimera"] }`。`mod.rs` 用 `tokio::task::spawn_blocking` 把 sage-core 桥进 async：`Fasta::parse -> IndexedDatabase -> SpectrumProcessor::new(150,true,0.0) -> Scorer{ ScoreType::SageHyperScore } -> rayon 并行 score`，随后 LDA `score_psms` 重打分、`spectrum_q_value`、`picked_peptide` / `picked_protein` FDR；另起 500ms 轮询任务汇报进度。`config.rs::build_sage_parameters` 把 `SearchParams` 映射为 sage `Parameters`。
- **pFind**：`{ "pFind", "3.x (not connected)", ["open_search","modification_localization"] }`，`adapters/pfind.rs` 为 stub，`search` / `health_check` 返回未实现 / `Unavailable`，`SshConfig` 预留远程 .cfg 生成与执行。

## 6. 修饰处理要点

- **固定修饰** 同时进前体质量（`apply_fixed_mods`）与碎片（`mod_delta_fragment`）。残基特异按命中残基逐一累加；末端 / 全局修饰（`residues` 为空）按 `ModPosition` 分流：`AnyNTerm` 计入 b 离子、`AnyCTerm` 计入 y 离子。
- **可变修饰** 先 `find_applicable_sites` 找位点，再 `enumerate_combinations` 在 `max_variable_modifications` 上限内枚举（含空组合，限制组合爆炸）。同一残基只改一次；N / C 端是两个不同化学位点（`NTERM_POS = usize::MAX`，`CTERM_POS = usize::MAX - 1`），同端互斥、异端可共存。
- **蛋白末端门控**：`ProteinNTerm` / `ProteinCTerm` 修饰仅当 `is_protein_nterm` / `is_protein_cterm` 为真才生效，前体侧（`apply_fixed_mods`）与碎片侧（`mod_delta_fragment`）逻辑镜像，保证母离子与碎片一致。
- **碎片与母离子一致**：`resolve_combined_mods` 把可变修饰拆成"非定域"（末端 / 全局，走 `mod_delta_fragment`）与"逐位点 `position_deltas`"；b/y 离子按各自覆盖到的残基沿链累加 `position_deltas`，避免同种修饰在多位点上被重复计数。

## 7. 错误处理

`SearchEngineError`（`error.rs`，thiserror）覆盖 `InvalidParams / FastaError{path,detail} / IoError / EngineNotFound / ExecutionError / NoInputSpectra`，并 `impl From<SearchEngineError> for CoreError`，转成结构化 `CoreError::SearchEngineError { engine, detail, suggestion }`（每个变体给出修复建议）。库代码不 `unwrap` / `expect`：`residue_mass` / `peptide_mass` 返回 `Option`，含非标准残基的肽被静默跳过；匹配中 `score` / `delta_ppm` 非有限值的候选被丢弃。各阶段经 `SearchDiagnostics`（`begin_stage / end_stage / fail_stage / set_error`）记录 `ErrorCategory`（Parameters / InputData / Database / Engine），并以 `tracing` span 记录蛋白数、肽数、谱图处理速率与 ETA。

往上返回 [README](README.md)。

