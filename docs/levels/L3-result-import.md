# L3 — 外部结果导入子系统

承接 [L2](L2-architecture.md)。本篇讲清 `crates/result-import` 如何把 DIA-NN、pFind、自定义
JSON 三类外部搜索结果，归一为 `core::search_result::SearchResult`，并经 MCP 工具
`import_search_results` 产出 `run_id`，复用注释 / XIC / 汇总链路（等价于 L2 §3 的原生搜索旁路）。

## 1. 职责与位置

该子系统让用户跳过 ProteinCopilot 原生搜索，把已有的第三方搜索结果直接接入同一套下游分析；
核心工作是"解析外部格式、按 scan 关联回 mzML、归一为统一结构"三步，使后续注释、XIC、汇总、
蛋白推断无须区分结果来源。

- 输入：外部工具结果文件（`.parquet` / `.tsv` / `.json` / `.spectra`）+ 一个 mzML 目录。
- 输出：统一 `SearchResult`（PSM/肽/蛋白三级 + summary）与 `ImportResult`（含 `run_id`、`MatchReport`）。
- 下游：`run_id` 落 MCP `run_cache`，供 `annotate_spectrum` / `extract_xic` / `generate_summary` /
  `infer_proteins` 复用，使外部结果与原生搜索结果在工具层无差别。
- 依赖：`core`、`search-engine`（`chemistry::{peptide_mass, peptide_mz}`）、`spectrum-io`（scan 匹配）。

位置：`crates/result-import/src/`；唯一使用方 `crates/mcp-server/src/tools.rs::import_search_results`。

## 2. 支持格式：识别与列映射

`detect_format(path)`（`lib.rs`）按扩展名分派为 `ImportFormat`：

```text
.json     -> CustomJson
.parquet  -> DiannParquet
.spectra  -> PFindSpectra
.tsv/.txt -> pfind_tsv::detect(header) ? PFindTsv : FormatDetectionFailed
```

四个输入解析器都实现 `ResultParser`，统一产出 `Vec<ImportedPsm>`，随后汇聚成 **一个**统一出口
`SearchResult`（任务所称 "unified" 即此目标，并非第五种输入解析器）。

- **DIA-NN parquet**（`diann.rs`）。列 `Modified.Sequence` / `Precursor.Charge` / `Precursor.Mz` /
  `RT` / `Q.Value` / `Run`，可选 `Protein.Names`。`Modified.Sequence` 形如
  `_AAAC(UniMod:4)DM(UniMod:35)K_`，用正则 `\(UniMod:(\d+)\)` 解析，剥首尾 `_`，N 端修饰记 position 0。
  `score = 1.0 - Q.Value`（反转为越大越好），`q_value = Q.Value`，`RT` 原样按分钟存；`matched_scan`
  留空，需 RT 匹配。解析时对缺失或非法字段（`Q.Value`、`Run`、`Modified.Sequence` 为空，或
  `charge <= 0`、`precursor_mz <= 0`、`RT` 缺失）逐行跳过并累加告警计数；`filter_qvalue`（默认 0.01）
  滤除高 q 值行，`run_filter` 可只导入指定 `Run`。
- **pFind TSV**（`pfind_tsv.rs`，16 列）。直接带 `ScanNo`，故 `matched_scan = Some(scan_no)`，跳过 RT 匹配。
  `MH+` 经 `mh_plus_to_precursor_mz` 换算 m/z；`Modifications` 形如 `10,Carbamidomethyl[C];` 由内置质量表
  解析；`Proteins` 按 `/` 切分。`detect` 看表头同时含 `PeptideSequence`、`ScanNo`、`FileName`。
  质子质量常数 `PROTON_MASS = 1.007276`；修饰位点串中的 `ProteinN-term`、`AnyN-term`、`C-term` 等
  映射为对应 `ModPosition`，普通残基记 `Anywhere`，未知修饰名质量记 0。
- **custom JSON**（`custom_json.rs`，hela.json 风格）。字段
  `sequence/charge/modify([[pos,id]])/rt/precursor_mz/raw_title/protein_names`；`modify` 经
  `UnimodDb::to_modification` 解析；`score`、`q_value` 均为 `None`；`matched_scan` 留空，需 RT 匹配。
  `modify` 中长度非 2 的项或解析失败的修饰会被跳过并告警（超过 5 条则折叠）。
- **pFind .spectra**（`pfind.rs`）。仅骨架，`parse` 返回 `Other("not yet implemented")`，MCP 层显式拒绝。

## 3. 模块边界

解析与匹配解耦：四个 parser 只产出 `ImportedPsm`，scan 关联与归一各自独立；新增格式时通常只需实现
`ResultParser` 并在 `detect_format` 注册。

| 文件 | 职责 |
|------|------|
| `lib.rs` | `ImportedPsm`/`ImportResult`/`MatchReport`/`FileMatchStats`、`ResultParser` trait、`ImportFormat` + `detect_format` |
| `diann.rs` | `DiannParser`（`filter_qvalue`/`run_filter`）+ arrow 读列 + `Modified.Sequence` 解析 |
| `pfind_tsv.rs` | `PFindTsvParser`：TSV 列映射 + MH+ -> m/z + 修饰/蛋白解析 + `detect` |
| `custom_json.rs` | `CustomJsonParser`：serde 反序列化 |
| `pfind.rs` | `PFindParser` 骨架 + `.spectra` `detect` |
| `scan_matcher.rs` | `match_scans` / `find_best_match` / `find_scan_by_rt` / `collect_ms2_info` / `validate_raw_name` |
| `converter.rs` | `build_search_result`：`ImportedPsm` -> `SearchResult` + `ImportResult` |
| `unimod.rs` | `UnimodDb`（`builtin` / `from_xml` / `to_modification`） |
| `error.rs` | `ResultImportError` |

## 4. 关键数据结构

解析中间体 `ImportedPsm`（RT 恒为分钟，`matched_scan` 在匹配后填）：

```rust
pub struct ImportedPsm {
    pub sequence: String,
    pub charge: i32,
    pub precursor_mz: f64,
    pub rt_min: f64,                  // 分钟（各 parser 解析时统一）
    pub modifications: Vec<Modification>,
    pub score: Option<f64>,          // DIA-NN=1-Q.Value, pFind=FinalScore, JSON=None
    pub q_value: Option<f64>,
    pub protein_accessions: Vec<String>,
    pub raw_name: String,            // 不含扩展名，关联 {raw}.mzML
    pub matched_scan: Option<u32>,   // 由 scan 匹配填充
    pub rt_delta_min: Option<f64>,
}
```

`build_search_result` 只纳入 `matched_scan.is_some()` 的 PSM，转 `core::Psm`，聚合肽/蛋白并算 summary：
- `engine_info = { name: "imported", version: format_name }`；
- `params_used` 为占位（`database_path: "imported"`，外部结果不携带搜索参数）；
- `RunMetadata::new(params, engine_info, input_files)`，`input_files` 是 MCP 层排序后的 mzML 路径；
- `ImportResult.raw_files` 经 `HashSet` 去重后 `sort()`，确定性输出（见 `converter.rs` 的排序测试）。

summary 同样在 `converter` 内确定性算出：若任一 PSM 带 `q_value`，则 `psms_at_1pct_fdr` 及肽/蛋白同名
指标按 `q <= 0.01` 过滤计数，否则（如 custom JSON 无 q 值）全部计入；总谱图数取各文件 MS2 计数之和，
`identification_rate = psms_at_1pct_fdr / total_spectra`；`median_score` 与 `median_delta_mass_ppm`
用中位数函数算，电荷与修饰分布按出现次数聚合。`to_core_psm` 用序列与修饰质量算理论 m/z，再填
`delta_mass_ppm = (precursor_mz - calculated_mz) / calculated_mz * 1e6`，作为质量偏差留给下游展示。

## 5. 主流程伪代码

```text
import_search_results(input):                  # mcp-server/tools.rs
  校验: result_file 存在 / mzml_dir 是目录 / rt_tol>=0 / filter_qvalue in [0,1]
  unimod = unimod_path ? UnimodDb::from_xml : UnimodDb::builtin
  format = (format == "auto") ? detect_format(path) : 显式枚举映射
  psms = parser.parse(path, unimod)            # 列映射在各 parser 内完成
  if run_filter: psms.retain(raw_name == filter)   # pFind TSV 等在外层补过滤
  if psms.is_empty(): 报错

  if psms 全部 matched_scan.is_some():          # pFind TSV：自带 scan
      直接构造 MatchReport(matched = total)
  else:                                         # DIA-NN / JSON：RT + 隔离窗匹配
      match_scans(&mut psms, {rt_tol, mzml_dir}, create_indexed_reader)
  if match_report.matched == 0: 报错

  raw_names = 去重 + sort_unstable
  for raw in raw_names: validate_raw_name(raw); 拼 {raw}.mzML 或 .mzml
  (sr, import) = build_search_result(psms, report, format_name, mzml_files)
  sr.metadata.status = Completed; 记 duration_sec
  run_cache.insert(run_id, sr); 落 history
  return import                                 # 含 run_id
```

列映射的两个确定性换算（parser / converter 内）：

```rust
// charge<=0 守卫：跳过或回退，避免 peptide_mz 出现除零/崩溃
let calculated_mz = if imported.charge > 0 {
    peptide_mass(&seq)
        .map(|m| peptide_mz(m + mod_mass, imported.charge))
        .unwrap_or(imported.precursor_mz)
} else {
    imported.precursor_mz
};

// MH+ -> precursor m/z（pFind TSV）
fn mh_plus_to_precursor_mz(mh: f64, z: i32) -> f64 {
    (mh + (z as f64 - 1.0) * PROTON_MASS) / z as f64
}
```

`match_scans` 先按 `raw_name` 分组，每组定位 `{raw}.mzML`（找不到再试 `.mzml`，仍无则 `MzmlNotFound`
并列出目录内可选文件），经 `collect_ms2_info` 从内存索引读取该文件全部 MS2 的 `(scan, RT, 隔离窗)`。
MCP 层用 `create_indexed_reader`（带 `.mzML.idx` 磁盘缓存）作为 reader 工厂；当下游只给 RT 不给 scan 时，
`find_scan_by_rt` 走 `reader.find_by_rt` 的 O(log N) 二分。

scan 匹配核心 `find_best_match`：MS2 按 RT 升序后 `partition_point` 二分定位窗口起点，沿途比对隔离窗
`[target - lower, target + upper]`，取 RT 最近者；若该 MS2 无隔离窗信息则仅按 RT 接受（DDA 回退）。

## 6. 安全要点

- **raw_name 路径校验**：`raw_name` 来自外部结果文件，在 `mzml_dir.join("{raw}.mzML")` 之前必须经
  `validate_raw_name`——空串、`.`、`..`、含 `/`、`\`、平台分隔符，或 `Path::file_name() != raw_name`
  一律判为 `InvalidRawName`，阻止 `../` 穿越逃出 `mzml_dir`。`match_scans` 内每组先校验；pFind TSV 走
  跳过分支不进 `match_scans`，故 MCP 层在拼 mzML 路径前再补一次 `validate_raw_name`。
- **charge 守卫**：DIA-NN 与 pFind TSV 在解析期对 `charge <= 0` 的行直接跳过（并 `tracing::warn`）；
  `converter` 再设一层兜底——仅 `charge > 0` 才计算理论 m/z，否则回退用 `precursor_mz`，避免把非法电荷
  传入 `peptide_mz`。

## 7. 错误处理

`ResultImportError`（`error.rs`，`thiserror`）覆盖：`FileNotFound`、`FormatDetectionFailed`、
`MissingColumn { column, expected }`、`UnknownUnimodId`、`InvalidModPosition`、
`MzmlNotFound { raw_name, dir, available }`、`InvalidRawName`、`NoMatchingScan`，以及 `#[from]` 透传的
`JsonError` / `XmlError` / `ParquetError` / `ArrowError` / `IoError` 和兜底 `SpectrumIo(String)`、
`Other(String)`。库代码零 `unwrap`/`expect`；MCP 层把它映射为结构化 `ErrorData`
（`INVALID_PARAMS` / `INTERNAL_ERROR`），并对"解析 0 条 PSM""匹配 0 个 scan""`.spectra` 未实现"
给出明确可操作提示。

导入成功后，MCP 层把 `SearchResult` 写入 `run_cache`（满则先驱逐），落一条 history 记录，并将运行
元数据状态置 `Completed`、记录耗时；若导入的 PSM 跨多个原始文件，会额外告警，提示下游注释/XIC
目前只取第一个文件。

---

返回 [README](README.md)。
