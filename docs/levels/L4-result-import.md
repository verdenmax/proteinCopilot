# L4 — result-import crate

承接 [L3-result-import](L3-result-import.md)，回溯 [L2](L2-architecture.md)。本篇只聚焦 `crates/result-import` 一个 crate 的模块级 API 与代码骨架；三类格式的语义叙述见 L3，此处签名、常量、守卫均按源码核验，行号随附，完整逻辑以源码为准。

## 1. 用途 + 位置 + 依赖

`protein-copilot-result-import` 把 DIA-NN parquet / pFind TSV / 自定义 JSON 三类外部搜索结果解析为中间体 `ImportedPsm`，按 scan 关联回 mzML，再归一为 `core::SearchResult` + `ImportResult`。纯确定性，无 LLM、无全局可变状态。

位置：`crates/result-import/src/`，唯一使用方 `mcp-server/src/tools.rs::import_search_results`。

依赖（Cargo.toml:10-25）：`protein-copilot-core`、`protein-copilot-search-engine`（`chemistry::{peptide_mass, peptide_mz}`）、`protein-copilot-spectrum-io`（`SpectrumReader`）、`arrow 54` / `parquet 54`、`csv 1`、`regex 1`、`quick-xml`、`serde`、`thiserror`、`tracing`、`uuid`。模块导出见 lib.rs:8-15：`converter` / `custom_json` / `diann` / `error` / `pfind` / `pfind_tsv` / `scan_matcher` / `unimod`。

数据流：

```
result file -> detect_format -> ResultParser::parse -> Vec<ImportedPsm>
                                                            |
   matched_scan 缺失 -> scan_matcher::match_scans (RT + 隔离窗)
   matched_scan 已带 -> 直接构造 MatchReport (pFind TSV)
                                                            |
                       converter::build_search_result -> (SearchResult, ImportResult)
```

设计上解析与匹配解耦：四个 parser 只负责把各自格式映射成 `ImportedPsm`，scan 关联（`scan_matcher`）与归一聚合（`converter`）各自独立。新增格式时通常只需实现 `ResultParser` 并在 `detect_format` 注册一个扩展名分支，无须改动下游。

## 2. 模块级 API

| 文件:行 | 项（源码原样） | 职责 |
|---|---|---|
| lib.rs:30 | `struct ImportedPsm` | 解析中间体；RT 恒分钟，11 字段，matched_scan 后填 |
| lib.rs:52/63/74 | `ImportResult` / `MatchReport` / `FileMatchStats` | 导入结果 + 匹配质量 + 每文件统计（三者均 `JsonSchema`） |
| lib.rs:81 | `trait ResultParser` | `parse(&self, &Path, &UnimodDb) -> Result<Vec<ImportedPsm>>`，`Send + Sync` |
| lib.rs:91 | `enum ImportFormat` | `CustomJson` / `DiannParquet` / `PFindSpectra` / `PFindTsv` |
| lib.rs:99 | `fn detect_format` | 扩展名分派；.tsv/.txt 再查表头 |
| diann.rs:24 | `struct DiannParser` | `filter_qvalue` / `run_filter`；arrow 读列 + `Modified.Sequence` 正则解析 |
| pfind_tsv.rs:39 / :225 | `struct PFindTsvParser` + `pub fn detect` | 16 列 TSV，自带 ScanNo；detect 看表头三列 |
| custom_json.rs:36 | `struct CustomJsonParser` | serde 反序列化 hela.json 数组 |
| pfind.rs:12 / :40 | `struct PFindParser` + `pub fn detect` | .spectra 骨架，parse 返回 `Other("not yet implemented")` |
| scan_matcher.rs:43 | `fn match_scans` | 按 raw_name 分组 + 定位 mzML + RT/隔离窗匹配 -> MatchReport |
| scan_matcher.rs:233 | `fn find_best_match` | RT 升序 `partition_point` 二分 + 隔离窗，取最近 |
| scan_matcher.rs:213 / :192 | `fn collect_ms2_info` / `fn find_scan_by_rt` | 读全部 MS2 元数据 / 单点 RT -> scan |
| scan_matcher.rs:288 | `fn validate_raw_name` | 路径穿越守卫（见 §3） |
| converter.rs:23 | `fn build_search_result` | ImportedPsm -> (SearchResult, ImportResult) |
| error.rs:7 | `enum ResultImportError` | thiserror，15 变体 |

辅助公共类型：`ScanMatcherConfig { rt_tolerance_min, mzml_dir }`（scan_matcher.rs:27）、`Ms2Info { scan_number, rt_min, isolation_window }`（:19）、`type ReaderFactory`（:41）。私有助手 `mh_plus_to_precursor_mz`（pfind_tsv.rs:146）、`parse_modified_sequence`（diann.rs:55）见 §3/§4。

## 3. 关键结构 / 常量 / 守卫

`ImportedPsm`（lib.rs:30-48，11 字段）：

```rust
pub struct ImportedPsm {
    pub sequence: String,
    pub charge: i32,
    pub precursor_mz: f64,
    pub rt_min: f64,                 // 各 parser 解析时统一为分钟
    pub modifications: Vec<Modification>,
    pub score: Option<f64>,          // DIA-NN=1-Q.Value, pFind=FinalScore, JSON=None
    pub q_value: Option<f64>,
    pub protein_accessions: Vec<String>,
    pub raw_name: String,            // 不含扩展名，关联 {raw}.mzML
    pub matched_scan: Option<u32>,   // 由 scan_matcher 填
    pub rt_delta_min: Option<f64>,
}
```

常量：`PROTON_MASS: f64 = 1.007276`（pfind_tsv.rs:36）；DIA-NN 默认 `filter_qvalue = Some(0.01)`（diann.rs:34）；UniMod 正则 `\(UniMod:(\d+)\)`（diann.rs:48）。

charge<=0 守卫（双层）：pfind_tsv.rs:111 解析期 `if charge <= 0 { warn; continue }` 直接跳过；diann.rs:187 用 `Some(c) if c > 0` 守卫；converter.rs:196 兜底 `if imported.charge > 0` 才算理论 m/z，否则回退 `precursor_mz`，避免把非法电荷传入 `peptide_mz`。

validate_raw_name 规则（scan_matcher.rs:288-302）：空串、`.`、`..`、含 `/`、`\`、平台分隔符 `MAIN_SEPARATOR`，或 `Path::file_name() != Some(raw_name)` 一律判 `InvalidRawName`，阻止 `../` 逃出 mzml_dir。

## 4. 简化源码片段

detect_format 分派（lib.rs:99-117，简化）：

```rust
match path.extension().and_then(|e| e.to_str()) {
    Some("json")    => Ok(ImportFormat::CustomJson),
    Some("parquet") => Ok(ImportFormat::DiannParquet),
    Some("spectra") => Ok(ImportFormat::PFindSpectra),
    Some("tsv") | Some("txt") =>
        if pfind_tsv::detect(path) { Ok(ImportFormat::PFindTsv) }
        else { Err(FormatDetectionFailed { path: path.into() }) },
    _ => Err(FormatDetectionFailed { path: path.into() }),
}
```

pFind TSV 列映射 + MH+ 换算（pfind_tsv.rs:100-147）：

```rust
let mh_plus: f64 = record.get(i_mh_plus).unwrap_or("0").parse().unwrap_or(0.0);
let charge: i32  = record.get(i_charge).unwrap_or("0").parse().unwrap_or(0);
if charge <= 0 { tracing::warn!(scan = scan_no, "skipping ..."); continue; } // 守卫
let precursor_mz = mh_plus_to_precursor_mz(mh_plus, charge);
// matched_scan = Some(scan_no) — pFind 自带 scan，跳过 RT 匹配

fn mh_plus_to_precursor_mz(mh_plus: f64, charge: i32) -> f64 {
    (mh_plus + (charge as f64 - 1.0) * PROTON_MASS) / charge as f64
}
```

注：`record.get(..).unwrap_or("0")` 是 `Option::unwrap_or`、`parse().unwrap_or(0.0)` 是 `Result::unwrap_or`，均不 panic，不违反禁 `unwrap()/expect()` 规约。

build_search_result 骨架（converter.rs:23-188）：

```rust
let core_psms: Vec<Psm> = psms.iter()
    .filter(|p| p.matched_scan.is_some())          // 仅纳入已匹配
    .map(to_core_psm).collect();
let peptides = aggregate_peptides(&core_psms);
let proteins = aggregate_proteins(&core_psms);
let has_qvalues = core_psms.iter().any(|p| p.q_value.is_some());
let psms_at_fdr = if has_qvalues {                 // 有 q 值按 1% 过滤，否则全计
    core_psms.iter().filter(|p| p.q_value.is_some_and(|q| q <= 0.01)).count()
} else { core_psms.len() };
let engine_info = EngineInfo { name: "imported", version: format_name, .. };
// params_used 为占位 (database_path: "imported")，外部结果不携带搜索参数
(SearchResult { run_id, psms: core_psms, peptides, proteins, summary, .. },
 ImportResult { run_id, match_report, imported_psm_count, raw_files, .. })
```

`to_core_psm` 用序列 + 修饰质量算理论 m/z，再填质量偏差 `delta_mass_ppm = (precursor_mz - calculated_mz) / calculated_mz * 1e6`（converter.rs:204）。私有助手 `aggregate_peptides`（:233）与 `aggregate_proteins`（:273）用 HashMap 按序列/accession 聚合并取 best_score 与最小 q 值，`median`（:294）取排序后中点；summary 的电荷与修饰分布按出现次数累加，均为确定性输出，便于跨次导入复现。

## 5. 调用链（mcp-server import_search_results）

`mcp-server/tools.rs::import_search_results`（tools.rs:3251）顺序，行号原样：

```
3296 unimod = unimod_path ? UnimodDb::from_xml : UnimodDb::builtin
3304 format = (input.format == "auto") ? detect_format : 显式枚举映射
3322 psms = parser.parse(result_path, unimod)         # 各 parser 内完成列映射
3347 if run_filter: psms.retain(raw_name == filter)   # pFind TSV 外层补过滤
3360 if psms.is_empty(): INVALID_PARAMS
3368 all_scans_present = psms.all(matched_scan.is_some())
       true  -> 直接构造 MatchReport(matched = total)            # pFind TSV
       false -> match_scans(&mut psms, config, create_indexed_reader)  # 3402
3413 if match_report.matched == 0: INVALID_PARAMS
3435 raw_names.sort_unstable(); for raw: validate_raw_name(raw)  # 3438
3464 build_search_result(psms, report, format_name, mzml_files)
3469 metadata.status = Completed; duration_sec
3473 run_cache.insert(run_id, sr)   # 落 history，返回 ImportResult(含 run_id)
```

两条分支的依据是 `all_scans_present`：pFind TSV 已携带 `ScanNo`，无须打开 mzML，直接按 raw_name 汇总 `MatchReport`；DIA-NN / JSON 缺 scan，才走 `match_scans` 的 RT + 隔离窗匹配。reader 工厂用 `spectrum-io::create_indexed_reader`（带 `.mzML.idx` 磁盘缓存），`collect_ms2_info` 经 `reader.list_ms2_meta` 从内存索引读 MS2 元数据，避免全量 streaming。`run_filter` 在 3347 处由 MCP 层对 pFind TSV 等不内建过滤的格式补一次 `retain`；若导入 PSM 跨多个原始文件，3448 处告警提示下游注释/XIC 暂只取首个文件。`PFindSpectra` 在 3334 处显式拒绝（未实现）。

## 6. 测试入口

```
cargo test -p protein-copilot-result-import --offline
```

51 个单元测试全通过（converter 5 + custom_json 3 + diann 6 + lib 6 + pfind 2 + pfind_tsv 14 + scan_matcher 10 + unimod 5），覆盖：detect_format 四格式 + 未知扩展、`Modified.Sequence` 各形态（无修饰 / 单修饰 / 多修饰 / N 端 / 裸序列）、MH+ 换算两电荷、修饰串与蛋白串解析、charge-0 跳过、find_best_match 隔离窗 / RT 容差 / 宽 DIA 窗 / DDA 回退 / 空表、match_scans 拒绝 `../evil` 穿越且不触工厂、build_search_result 排除未匹配 + 聚合肽 + raw_files 确定性排序。

---

回到 [README](README.md) 选择其它层级或子系统。
