# L4 - entrapment-analysis（模块级 API 与代码骨架）

承接 [L3-entrapment](L3-entrapment.md) 与 [L2](L2-architecture.md)。L3 已讲清 L0-L4 分级的概念与判定语义，本篇不复述，只聚焦 `crates/entrapment-analysis`（package `protein-copilot-entrapment-analysis`）的**模块级 API 表面**与**代码骨架**，所有签名、常量、默认值逐一核对源码。

## 1. 用途 + 位置 + 依赖

- 位置：`crates/entrapment-analysis`，package 名 `protein-copilot-entrapment-analysis`，全 crate ~11k LOC（src 10332 + tests 913），为最大的库 crate。
- 用途：把搜索结果中命中陷阱库（trap / entrapment）的 PSM，按其与靶库（target）蛋白的同源性分为 L0-L4 五级，并产出 HTML 报告、碎片离子溯源（provenance）与多靶共洗脱（co-elution）分析。

| 类别 | crate |
|------|-------|
| workspace 运行依赖 | `protein-copilot-core`、`protein-copilot-search-engine`、`protein-copilot-spectrum-io` |
| 外部运行依赖 | `parquet` / `arrow`（v54，读 DIA-NN report.parquet）、`csv`、`regex`、`serde_yaml`、`sha2`、`chrono` |
| 基础设施 | `serde` / `serde_json`、`thiserror`、`tracing` |
| dev-dependencies | `tokio`、`tempfile = "3"` |

依赖方向：`core` 提供 `Psm` / `SearchParams` / `MassTolerance`；`search-engine` 复用其 `fasta::parse_fasta`、`digest::{digest_with_length, is_standard_sequence, residue_mass}`；`spectrum-io` 提供 `create_indexed_reader` 做 O(1) 取谱与 `find_by_rt` 反查 scan。

## 2. 模块级 API 表

| 文件 | 关键导出（签名要点） | 职责一句话 |
|------|------|------|
| `lib.rs` | `EntrapmentAnalyzer::{new, classify, classify_all, summary}`；`trace_provenance_batch(&mut [ClassifiedPsm], &Path, &EntrapmentConfig) -> Result<u32, _>`；`extract_dia_windows`；`trace_multi_target_provenance(..) -> Result<(u32, Vec<MultiTargetProvenance>), _>` | crate 编排入口 |
| `types.rs` | `enum DiscriminabilityLevel{L0..L4}`、`PsmGroup`、`SubstitutionType`；`UnifiedPsm`、`ClassifiedPsm`、`LevelCounts`、`EntrapmentSummary`、`RazorFamily` + v4 多靶类型 | 共享数据结构（全 serde） |
| `config.rs` | `EntrapmentConfig::{from_yaml, from_yaml_str, validate}`；`SimilarityConfig`、`ProvenanceConfig`、`SilacConfig`、`Rule`、`ConflictResolution`、`UnmatchedPolicy` | YAML 配置解析 + 校验 |
| `error.rs` | `enum EntrapmentError`（11 变体：`ConfigError` `FastaError` `LoaderError` `ProvenanceError` `SpectrumError` ...） | 错误枚举 |
| `tagger.rs` | `Tagger::{new(&EntrapmentConfig), tag(&str) -> Result<PsmGroup, _>}` | 蛋白号 -> target/trap/ambiguous |
| `digest.rs` | `TargetDigestIndex::{from_fasta(&Path, u32, u16), find_similar(&str, u16, usize, &SimilarityConfig) -> Vec<SimilarityMatch>, has_exact, has_normalized, peptides_of_length}`；`normalize_li`、`TargetPeptide` | 靶库胰酶切索引 + k-mer 相似搜索 |
| `levenshtein.rs` | `edit_distance(&str, &str) -> u32`；`align(&str, &str) -> AlignmentResult{edit_distance, delta_mass_da, alignment_detail}` | 编辑距离 + 回溯对齐 + 质量差 |
| `similarity.rs` | `classify_single(&UnifiedPsm, PsmGroup, &TargetDigestIndex, &SimilarityConfig) -> ClassifiedPsm`；`hamming_diff` | L0-L4 判级内核 |
| `provenance.rs` | `trace_provenance(..) -> FragmentProvenance`；`enum IonOrigin{TrapOnly, TargetOnly, Shared, Unassigned}`、`AnnotatedPeak` | 单谱 trap/target b/y 离子溯源 |
| `coelution.rs` | `CoElutionIndex::{build(.., max_candidates: usize), find_co_eluting(&UnifiedPsm, &str) -> Vec<CoElutingCandidate>}`；`DiaWindow` | 同 DIA 窗 + RT 内共洗脱靶肽 |
| `multi_provenance.rs` | `trace_multi_target(.., &[CoElutingCandidate], ..) -> MirrorData` | 一个 trap 对多靶候选的镜像注释 |
| `mod_parser.rs` | `parse_modified_sequence(&str) -> (String, Vec<ParsedModification>)`；`unimod_delta_mass(u32) -> Option<f64>`（13 条 UniMod） | DIA-NN Modified.Sequence 解析 |
| `output.rs` | `write_classified_tsv`、`write_razor_errors_tsv`、`write_run_metadata`、`file_sha256(&Path) -> Result<String, _>`、`pub mod columns`（24 列名常量） | TSV / 元数据落盘 |
| `report.rs` | `render_report(&EntrapmentSummary, &[ClassifiedPsm], &Path) -> Result<(), _>` | 汇总 HTML 报告 |
| `mirror_plot.rs` | `render_mirror_plot(&FragmentProvenance, &Path)`、`generate_mirror_html(&FragmentProvenance) -> String` | 单 PSM Plotly 镜像图 |
| `multi_report.rs` | `render_multi_provenance_report`、`render_provenance_summary(&[MultiTargetProvenance], &Path)` + `generate_*_html` | 多靶逐 PSM / 汇总 HTML |
| `loader/{mod,pfind_tsv,diann_parquet,generic_tsv}` | `ResultFormat::from_path`、`load_psms(&Path, &ResultFormat, Option<&TsvColumnMap>) -> Result<Vec<UnifiedPsm>, _>`；`load_pfind_tsv`、`load_diann_parquet`、`load_generic_tsv`、`pfind_tsv::detect(&Path) -> bool` | 多格式结果导入 -> UnifiedPsm |

模块自下而上：`config`/`error`/`types` 是基座；`tagger`/`digest`/`levenshtein`/`similarity` 构成分级内核；`provenance`/`coelution`/`multi_provenance` 做谱级溯源；`output`/`report`/`mirror_plot`/`multi_report`/`loader` 处理 I/O；`lib.rs` 把它们编排成批处理流水线。

`loader` 把三种异构输入统一成 `UnifiedPsm`（剥离修饰的裸序列、分号拼接的蛋白号、分钟单位 RT、1 基 scan 号、`(位置, delta_mass)` 修饰对）：`ResultFormat::from_path` 按扩展名分发，`.parquet` 走 DIA-NN，`.tsv`/`.txt` 先用 `pfind_tsv::detect` 探测表头（同时含 `PeptideSequence` / `ScanNo` / `FileName` 才判为 pFind），否则按 `TsvColumnMap` 配置的列名读通用 TSV。DIA-NN 的 `Modified.Sequence` 经 `mod_parser` 把 `(UniMod:N)` 标注换算成 delta 质量。

## 3. 关键结构 / 默认值 / k-mer 公式

- `DiscriminabilityLevel`：`L0`(razor 错误，精确命中靶库) -> `L1`(I/L 同分异构) -> `L2`(近等重单残基替换) -> `L3`(可区分同源) -> `L4`(真陷阱)。`as_str()` 返回 `"L0".."L4"`，语义详见 L3。
- `SimilarityConfig` 默认（`config.rs` 的 `default_*` 函数）：`max_mismatches = 2`、`delta_mass_threshold_da = 1.0`（L2/L3 分界，别名 `delta_mz_threshold_da`）、`require_tryptic_ends = true`、`max_missed_cleavages = 2`、`len_tolerance = 2`、`enable_dipeptide_check / enable_qk_detection = true`。
- `ProvenanceConfig` 默认：`fragment_tolerance_ppm = 20.0`、`max_fragment_charge = 2`、`chimera_threshold = 0.3`、`min_peaks_for_analysis = 6`、`levels_to_trace = ["L2","L3","L4"]`、`rt_tolerance_min = 0.5`。SILAC 重标默认 `heavy_k_delta = 8.014199`、`heavy_r_delta = 10.008269`。
- k-mer 长度（鸽巢原理，`digest.rs` 构建期）：`kmer_k = (6 / (max_edit_distance + 1)).max(1)`，其中 `6` 为最小酶切肽长；默认 `max_edit_distance = 2` 时 `kmer_k = 2`。靶库胰酶切固定 `Enzyme::Trypsin`、长度窗 `6..=50`。
- `SubstitutionType`（`LIIsomer` / `QKSubstitution` / `IsobaricDipeptide` / `NearIsobaric` / `Distinguishable`）只是 L2/L3 上的信息性标注，**不改变** L0-L4 等级。`summary()` 另对 L0 命中按 `extract_family_name` 聚类，按计数降序取前 10 个 `RazorFamily` 暴露常见 razor 归属错误来源。

## 4. 简化源码片段

为突出主干，下面三段省略了日志、进度上报与字段细节。

**(a) `EntrapmentAnalyzer` 编排（`lib.rs`）**

```rust
pub struct EntrapmentAnalyzer { config: EntrapmentConfig, tagger: Tagger, index: TargetDigestIndex }

pub fn new(config: EntrapmentConfig, fasta_path: &Path) -> Result<Self, EntrapmentError> {
    let tagger = Tagger::new(&config)?;                         // 编译 target/trap 规则 + 载入号
    let index = TargetDigestIndex::from_fasta(fasta_path,
        config.similarity.max_missed_cleavages,
        config.similarity.max_mismatches)?;                    // 胰酶切 6..=50 + k-mer 倒排
    Ok(Self { config, tagger, index })
}
pub fn classify(&self, psm: &UnifiedPsm) -> Result<ClassifiedPsm, EntrapmentError> {
    let group = self.tagger.tag(&psm.protein_ids)?;            // target / trap / ambiguous
    Ok(classify_single(psm, group, &self.index, &self.config.similarity))
}
```

**(b) `find_similar` 的 k-mer 预筛（`digest.rs`）**

```rust
// 鸽巢: edit_distance <= d 的两序列必共享一个长 k 的精确 k-mer
let kmer_k = (6 / (max_edit_distance as usize + 1)).max(1);    // min_len = 6

pub fn find_similar(&self, query, max_edit_dist, len_tolerance, _cfg) -> Vec<SimilarityMatch> {
    let candidate_ids: Vec<u32> = if query.len() < self.kmer_k {
        (0..self.all_peptides.len() as u32).collect()          // 短肽 -> 回退全扫
    } else {
        let mut set = HashSet::new();                           // 倒排表取并集
        for kh in extract_kmers(query.as_bytes(), self.kmer_k) {
            if let Some(ids) = self.kmer_index.get(&kh) { set.extend(ids); }
        }
        set.into_iter().collect()
    };
    // 长度窗 [len-tol, len+tol] + edit_distance() 复核 -> 幸存者 align() -> 按 edit/|dm| 排序
}
```

**(c) `classify_single` 判级骨架（`similarity.rs`）**

```rust
pub fn classify_single(psm, group, index, cfg) -> ClassifiedPsm {
    if group != PsmGroup::Trap          { return /* L4, 无匹配信息 */; }   // 非 trap 直接 L4
    if index.has_exact(&psm.peptide)    { return /* L0 */; }              // 精确命中靶库
    if index.has_normalized(&psm.peptide) { return /* L1 (I->L 归一) */; }
    // Phase A: 同长 Hamming 扫 index.peptides_of_length(len)
    // Phase B: 跨长 index.find_similar() -> best_cross；二者取 edit/|dm| 更优者
    let adj = signed_dm - mod_mass_adjustment(&psm.modifications, &diff_pos); // 扣修饰质量
    let level = if adj.abs() < cfg.delta_mass_threshold_da { L2 } else { L3 };// 默认 1.0 Da
    // 全程无候选 -> BestMatch::None -> L4
}
```

## 5. 调用链

**entrapment-cli**（3 子命令 `Analyze` / `Report`(stub) / `Inspect`）：

```text
entrapment-cli analyze --results --config --target-fasta [--mzml-dir]
  -> EntrapmentConfig::from_yaml()
  -> loader::load_psms(results, format) -> Vec<UnifiedPsm>
  -> EntrapmentAnalyzer::new(config, target_fasta)
  -> .classify_all(&psms) -> Vec<ClassifiedPsm>
  -> [--mzml-dir] trace_provenance_batch() | trace_multi_target_provenance()
  -> .summary() -> EntrapmentSummary
  -> output::{write_classified_tsv, write_razor_errors_tsv, write_run_metadata}
  -> report::render_report() -> entrapment_report.html
```

`Inspect` 子命令绕过 PSM 加载，对单条肽直接 `TargetDigestIndex::from_fasta` + `classify_single`。

**MCP Server（4 个 tool）**：

```text
classify_entrapment_hits -> load_psms + classify_all + [trace_provenance_batch] + render_report + summary
analyze_entrapment_stats -> 读回 classified.tsv（用 output::columns 列名）-> 分级 / 家族 / dm 统计
find_similar_targets     -> from_fasta + classify_single(单 UnifiedPsm, group=Trap)
annotate_provenance      -> get_or_create_reader 取 scan -> trace_provenance -> render_mirror_plot
```

四个 tool 均为确定性计算，HTML / 统计不经 LLM。

`trace_provenance_batch` 的内部约束：只追溯 `group == Trap` 且等级落在 `levels_to_trace` 的 PSM；按 `spectrum_file` 分组复用 `create_indexed_reader`；scan 号缺失时用 `find_by_rt`（容差取 `(rt_stop - rt_start)/2` 或配置 `rt_tolerance_min`）反查；峰数少于 `min_peaks_for_analysis` 则跳过；最后按 `shared_ratio > chimera_threshold` 置 `is_chimeric`。DIA-NN parquet 无 scan 号时会发 `warn` 并整体跳过溯源。

**v4 多靶 SILAC 链**（`trace_multi_target_provenance`，用于 trap 落在多个共洗脱靶肽之间的判别）：

```text
trace_multi_target_provenance(classified, all_psms, all_groups, mzml_dir, config, out_dir)
  -> 逐 run: create_indexed_reader + extract_dia_windows()
  -> CoElutionIndex::build(all_psms, groups, windows, silac, max_co_eluting_candidates)
  -> 逐 eligible trap PSM: 解析 scan(scan_number | find_by_rt)
       -> index.find_co_eluting(trap, run) -> Vec<CoElutingCandidate>(light + heavy)
       -> 读 light/heavy scan -> trace_multi_target() -> MirrorData
       -> [generate_per_psm_reports] render_multi_provenance_report() -> provenance/*.html
  -> render_provenance_summary() -> provenance_summary.html
  返回 (traced_count, Vec<MultiTargetProvenance>)
```

## 6. 测试入口

```bash
cargo test -p protein-copilot-entrapment-analysis --offline
```

| 测试二进制 | 通过数 |
|-----------|--------|
| `unittests src/lib.rs` | 180 |
| `tests/v2_edit_distance.rs` | 4 |
| `tests/v3_e2e_provenance.rs` | 8 |
| `tests/v4_multi_target.rs` | 1 |
| doc-tests | 1 |
| **合计** | **194** |

## 核对签名（带行号）

| 文件:行 | 签名 |
|---------|------|
| `lib.rs:53` | `pub fn new(config: EntrapmentConfig, fasta_path: &Path) -> Result<Self, EntrapmentError>` |
| `lib.rs:68` | `pub fn classify(&self, psm: &UnifiedPsm) -> Result<ClassifiedPsm, EntrapmentError>` |
| `lib.rs:79` | `pub fn classify_all(&self, psms: &[UnifiedPsm]) -> Result<Vec<ClassifiedPsm>, EntrapmentError>` |
| `lib.rs:116` | `pub fn summary(&self, classified: &[ClassifiedPsm]) -> EntrapmentSummary` |
| `lib.rs:192` | `pub fn trace_provenance_batch(classified: &mut [ClassifiedPsm], mzml_dir: &Path, config: &EntrapmentConfig) -> Result<u32, EntrapmentError>` |
| `lib.rs:454` | `pub fn trace_multi_target_provenance(..) -> Result<(u32, Vec<types::MultiTargetProvenance>), EntrapmentError>` |
| `types.rs:19` | `pub enum DiscriminabilityLevel { L0, L1, L2, L3, L4 }` |
| `config.rs:121` | `pub struct SimilarityConfig`（默认 `max_mismatches=2`,`delta_mass_threshold_da=1.0`,`len_tolerance=2`，见 `config.rs:272-286`） |
| `digest.rs:106` | `pub fn from_fasta(path: &Path, max_missed_cleavages: u32, max_edit_distance: u16) -> Result<Self, EntrapmentError>` |
| `digest.rs:207` | `let kmer_k = (min_peptide_len / (max_edit_distance as usize + 1)).max(1);`（`min_peptide_len = 6`） |
| `digest.rs:290` | `pub fn find_similar(&self, query: &str, max_edit_dist: u16, len_tolerance: usize, _config: &SimilarityConfig) -> Vec<SimilarityMatch>` |
| `similarity.rs:248` | `pub fn classify_single(psm: &UnifiedPsm, group: PsmGroup, index: &TargetDigestIndex, config: &SimilarityConfig) -> ClassifiedPsm` |
| `levenshtein.rs:24` | `pub fn edit_distance(a: &str, b: &str) -> u32` |
| `levenshtein.rs:63` | `pub fn align(a: &str, b: &str) -> AlignmentResult` |
| `tagger.rs:152` | `pub fn tag(&self, protein_ids: &str) -> Result<PsmGroup, EntrapmentError>` |
| `coelution.rs:143` | `pub fn find_co_eluting(&self, trap: &UnifiedPsm, run: &str) -> Vec<CoElutingCandidate>` |
| `provenance.rs:102` | `pub fn trace_provenance(observed_mz, observed_intensity, trap_sequence, target_sequence, trap_modifications, fragment_tolerance, max_fragment_charge) -> FragmentProvenance` |
| `report.rs:130` | `pub fn render_report(summary: &EntrapmentSummary, classified: &[ClassifiedPsm], output_path: &Path) -> Result<(), EntrapmentError>` |
| `mod_parser.rs:59` | `pub fn parse_modified_sequence(modified_seq: &str) -> (String, Vec<ParsedModification>)` |
| `loader/mod.rs:69` | `pub fn load_psms(path: &Path, format: &ResultFormat, tsv_config: Option<&TsvColumnMap>) -> Result<Vec<UnifiedPsm>, EntrapmentError>` |

返回目录 [README](README.md)。
