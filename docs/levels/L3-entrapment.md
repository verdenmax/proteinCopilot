# L3 — entrapment 同源性分级

承接 [L2](L2-architecture.md)。本篇讲清 `entrapment-analysis` crate 如何把"外部搜索命中"按其肽段与目标蛋白库的同源关系分到 L0-L4 五级，给出每级的源码判定边界、主流程伪代码、以及消除 HashMap 顺序的确定性手法。代码以 `crates/entrapment-analysis/src/` 为准（包名 `protein-copilot-entrapment-analysis`），CLI 见 `crates/entrapment-cli`（`protein-copilot-entrapment-cli`）。

## 1. 职责与位置

陷阱库（entrapment / trap）分析回答一个问题：命中陷阱库的 PSM 究竟是"真陷阱命中"，还是只是目标蛋白同源肽被错误归属？

- 输入：一份外部搜索结果（DIA-NN parquet / pFind TSV / 通用 TSV）+ 一个目标 FASTA 库 + 一份 YAML 配置（target/trap 规则、相似度阈值）。
- 处理：把每条 PSM 先打 group（target/trap/ambiguous），再对 trap PSM 计算它与目标库消化肽的同源等级 L0-L4。
- 输出：`classified.tsv`（逐 PSM 分级）、`razor_errors.tsv`（L0 razor 错误）、`run_metadata.json`（可复现快照）、`entrapment_report.html`（自包含报告）；给了 `--mzml_dir` 还会做碎片离子溯源与共洗脱多靶报告。

它处于 L2 的"计算/引擎层"，依赖 `core` + `search-engine`（借其 `digest`/`fasta` 与 `residue_mass`）+ `spectrum-io`（溯源读谱）。本子系统经三轮演进：v2 把同长 Hamming 升级为跨长编辑距离 + k-mer 预筛，v3 增碎片离子溯源，v4 增多靶共洗脱镜像报告——故 `types.rs` 同时保留 Hamming（`mismatches`）与 alignment（`edit_distance`/`alignment_detail`）两套字段。

## 2. 模块边界

三条主线——加载与打标、消化与索引、判级与溯源——各文件单一职责、互不越界：

```
loader/{diann_parquet,pfind_tsv,generic_tsv}.rs  外部结果 -> UnifiedPsm
config.rs        YAML 配置 + 阈值默认值
tagger.rs        protein_ids -> PsmGroup（target/trap/ambiguous）
digest.rs        目标库消化 + 索引（exact/normalized/by_length/k-mer）
levenshtein.rs   编辑距离 + 对齐回溯 + delta-mass
similarity.rs    classify_single：L0-L4 判级核心
mod_parser.rs    DIA-NN Modified.Sequence -> (pos, delta)
provenance.rs        单靶碎片溯源（trap/target/shared/unassigned）
multi_provenance.rs  多靶溯源（N 个共洗脱目标）
coelution.rs     同一 DIA 窗 + RT 共洗脱目标查找
mirror_plot.rs / report.rs / multi_report.rs   HTML 报告
output.rs        classified.tsv / razor_errors.tsv / 元数据
lib.rs           门面 EntrapmentAnalyzer + 批量溯源调度
```

`lib.rs` 暴露门面 `EntrapmentAnalyzer`（`new`/`classify`/`classify_all`/`summary`，lib.rs:53/68/79/116）。

## 3. L0-L4 分级定义

枚举见 `types.rs:19`（`DiscriminabilityLevel`），判定逻辑全在 `similarity.rs` 的 `classify_single`（similarity.rs:248）。阈值取自 `SimilarityConfig`（默认：`max_mismatches=2`、`delta_mass_threshold_da=1.0` Da、`len_tolerance=2`，config.rs:272-286）。判级是**优先级短路链**，命中即返回。

分级前先由 `Tagger::tag`（tagger.rs:151）按 `protein_ids`（分号分隔）定 group：任一 accession 命中 target 规则即 `is_target`，命中 trap 规则即 `is_trap`；二者皆中按 `conflict_resolution`（PreferTarget/PreferTrap/MarkAmbiguous）裁决，皆不中按 `unmatched`（Ignore/Target/Trap/Error）兜底。只有 group==Trap 进入下面的 L0-L4 链；target/ambiguous 直接 L4 占位、不计入报告级别统计。

- **L0 — razor 归属错误**（similarity.rs:273）：`index.has_exact(peptide)` 为真，即陷阱肽与某目标肽**序列完全相同**。这条肽不能区分两库，命中陷阱只是 razor 归属把它分给了陷阱蛋白。`mismatches=0, delta_mass=0`。
- **L1 — L/I 同分异构**（similarity.rs:292）：非精确，但 L/I 归一后相等（`has_normalized`，把所有 `I` 换成 `L`，digest.rs:37）。Leu/Ile 单同位素质量相同（113.084064 Da），质谱不可分；`substitution_type=LIIsomer`。
- **L2 — 近等重同源**（similarity.rs:462 / 498）：进入编辑距离扫描后找到最佳目标，且修饰校正后 `|delta_mass| < delta_mass_threshold_da`（默认 < 1.0 Da）。即存在 1-2 处差异但质量上近乎不可分（如 Q<->K ~36.4 mDa、N<->GG、Q<->AG）。
- **L3 — 可区分同源**（similarity.rs:465 / 501）：同样找到了编辑距离 <= `max_mismatches` 的目标，但 `|delta_mass| >= 阈值`，质量上可区分。
- **L4 — 真陷阱命中**（similarity.rs:430）：编辑距离扫描 `BestMatch::None`，目标库内无近邻，无任何匹配信息。

> 边界要点：L2 与 L3 的唯一分水岭是**修饰校正后的 |delta_mass| 与 1.0 Da 阈值**（不是 mismatch 数）；mismatch/长度差只决定是否进入候选（`mm>max_mismatches` 或 L/I-only 会被跳过，similarity.rs:349-354）。另外，非 trap 的 PSM（target/ambiguous）一律直接返回 L4 且无匹配信息（similarity.rs:255），但 `summary` 只把 trap PSM 计入 `level_counts`（lib.rs:129），故报告里的 L4 = 真陷阱命中。此外 `substitution_type`（QKSubstitution / IsobaricDipeptide{N:GG、Q:AG} / NearIsobaric / Distinguishable，见 `categorize_substitution`，similarity.rs:77）仅作注释，不改变 L 级。

判级核心（简化）：

```rust
// similarity.rs:248
pub fn classify_single(psm, group, index, config) -> ClassifiedPsm {
    if group != Trap { return L4_no_match; }          // 255
    if index.has_exact(&psm.peptide)      { return L0; }   // 273
    if index.has_normalized(&psm.peptide) { return L1(LIIsomer); } // 292
    // Phase A：同长 Hamming 扫描；Phase B：跨长 Levenshtein（k-mer 预筛）
    match overall_best {
        None            => L4,                              // 430
        Hamming|Cross   => if adj_abs_dm < config.delta_mass_threshold_da
                               { L2 } else { L3 },          // 462 / 498
    }
}
```

## 4. 主流程伪代码

`entrapment analyze` 子命令（main.rs:156）把六步串起来：

```
load config (YAML)                                  config.rs
psms   = load_psms(results, format)                 loader/*  -> Vec<UnifiedPsm>
index  = TargetDigestIndex::from_fasta(fasta,        digest.rs:106
              max_missed_cleavages=2, max_edit=2)
   |- Trypsin 消化，肽长 6..50，去非标准残基        digest.rs:142
   |- exact_set / normalized_set / by_length         首条蛋白胜出，保 FASTA 序
   |- kmer_k = (6 / (max_edit+1)).max(1)，建倒排索引  digest.rs:207
for psm in psms:                                    classify_all  lib.rs:79
   group = tagger.tag(psm.protein_ids)              tagger.rs:151
   classify_single(psm, group, index, config):
      L0: has_exact         -> 返回
      L1: has_normalized    -> 返回
      else 候选预筛 + 判级:
         Phase A: by_length[len] 同长 Hamming        similarity.rs:341
         Phase B: index.find_similar(...)            digest.rs:290
            query 切 k-mer -> 命中 kmer_index 的肽为候选（鸽巢原理）
            长度窗 [len-tol, len+tol] 过滤 -> Levenshtein 校验 -> align()
         取 Hamming/Cross 更优者，按 |delta_mass| vs 1.0 -> L2 / L3 / L4
(可选) trace_provenance_batch / 多靶共洗脱           provenance / coelution
   coelution: 按 run 分组、rt_start 排序，partition_point 取 RT 重叠候选
summary -> top_razor_families                       lib.rs:116
write classified.tsv / razor_errors.tsv / metadata  output.rs
render_report(summary, classified) -> HTML          report.rs:130
```

k-mer 鸽巢预筛的正确性来自 `kmer_k = min_len/(max_edit+1)`：编辑距离 <= e 时至少有一个长 k 子串完整保留，故真匹配必出现在某条 k-mer 倒排表里，预筛不漏（query 短于 k 时回退全扫，digest.rs:310）。

编辑距离与 delta-mass 同源由 `levenshtein.rs` 提供：`edit_distance`（单行 DP，levenshtein.rs:24）做快筛，`align`（回溯，levenshtein.rs:63）给出替换/插入/删除明细与按残基质量累加的 `delta_mass_da`。明细串含三类 op（替换、`ins:`、`del:`）；跨长匹配的替换位点由对齐串解析而非按位 zip，避免 indel 错位导致质量校正算错（`extract_substitution_positions_from_alignment`，similarity.rs:204）。

## 5. 确定性要点

候选来自 `HashSet`/`HashMap`，迭代序不稳定；本子系统靠"显式排序键"把不确定性消干净：

- **相似度候选**：`find_similar` 末尾按 `(edit_distance, |delta_mass|, target_peptide, target_protein)` 全序排序（digest.rs:374），保证同分候选下"最佳目标"唯一。Hamming 与 Cross 的取舍也用严格不等式 + `|delta_mass|` 兜底（similarity.rs:401-405）。
- **同长 Phase A**：`by_length[len]` 是 `Vec`，按 FASTA 顺序 `push`（digest.rs:172），择优用 `mm < best || (mm==best && abs_dm < best_dm)`（similarity.rs:357），不依赖 HashMap 序；`exact_to_protein`/`normalized_to_original` 亦"首条蛋白胜出"（digest.rs:154-163）。
- **razor 家族**：`top_razor_families` 按 `(count desc, family asc)` 排序后 `truncate(10)`（lib.rs:160），同票稳定。
- **共洗脱**：每个 run 的目标先按 `rt_start` 排序（coelution.rs:115），才能用 `partition_point` 取 `rt_start <= trap_rt_stop` 的上界、再线性过滤 `rt_stop >= trap_rt_start`（coelution.rs:174）。

归根结底，全部确定性都立在两条之上——"目标库按 FASTA 顺序构建索引（含重复肽首条胜出）"与"查询/聚合结果显式排序键"——故同一输入多次运行逐字节一致，报告可审计。

## 6. 跨 crate 交互与错误处理

- 复用 `search-engine`：`digest::digest_with_length`/`residue_mass`、`fasta::parse_fasta`（消化、残基质量、FASTA 解析），不自造算法。
- 复用 `core`：`search_params::Enzyme`、`MassTolerance/ToleranceUnit`（溯源容差）。
- `spectrum-io`：溯源时按 `spectrum_file`/scan 读 MS2。
- 错误：crate 内 `EntrapmentError`（thiserror，error.rs）覆盖 Config/Fasta/Loader/Io/Output/Report/Tagging/Provenance/Spectrum，全部带 `path`/`detail`；库代码不 `unwrap`（CLI 顶层统一 `eprintln!` + `exit(1)`）。可观测性用 `tracing`：消化、分类批次都打 progress/rate/eta span。
- 溯源（v3/v4，可选 `--mzml_dir`）：`ProvenanceConfig` 默认只追 `levels_to_trace=[L2,L3,L4]`、碎片容差 20 ppm、`chimera_threshold=0.3`、谱图至少 6 峰（config.rs:248-262）；`trace_provenance`（provenance.rs:102）把每个观测峰判为 trap/target/shared/unassigned，`shared_ratio` 超阈即标 `is_chimeric`，多靶版再按 DIA 窗 + RT 重叠拉入 N 个共洗脱目标（含 SILAC 重标）出镜像图。
- 可复现：`run_metadata.json` 落输入/FASTA 的 SHA-256、config 快照、版本、各级计数（output.rs / main.rs:259）。

延伸：报告渲染 `render_report` 只取 trap PSM 入明细表，序列化 JSON 注入 HTML 模板并转义 `<`/`>` 防注入（report.rs:136-161）。

回到 [README](README.md)。
