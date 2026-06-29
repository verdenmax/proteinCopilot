# L4 — integration-tests crate

回溯 [L2](L2-architecture.md)。本篇聚焦 `crates/integration-tests` 一个 crate：它处在依赖图最顶端，把 core / spectrum-io / param-recommend / search-engine / fdr / report / xic / result-import / entrapment-analysis 串成可运行的端到端断言。所有文件名、测试名、fixtures 均按源码核验，完整逻辑以源码为准。

## 1. 用途 + 位置

`protein-copilot-integration-tests`（`publish = false`，纯测试不发布）有两重身份：

- `src/lib.rs` 编译为名为 `test_helpers` 的库（`[lib] name = "test_helpers"`），提供合成谱图 / 碎片构造器；
- `tests/*.rs` 是跨 crate 集成测试，`use test_helpers::*` 复用这些构造器。

`Cargo.toml` 把 9 个工作区 crate 加 `serde_json / serde_yaml 0.9 / tokio / tempfile 3` 全列为普通 `dependencies`（lib 与 tests 共用，无 `dev-dependencies`）。各 crate 自己的单元测试只覆盖本 crate 内部逻辑；本 crate 专门补上"接缝"：谱图读取产出的 `Spectrum` 能否喂进搜索引擎、搜索产出的 `Psm` 能否被注释与报告消费、外部结果能否还原成统一 PSM。任一 crate 改了公共契约却忘了同步下游，集成测试会第一时间红灯。它不被任何 crate 依赖，只验证别人协作是否正确：

```
test_helpers (synthetic Spectrum / fragments)
   |
   +-> annotation_scenarios   -> search-engine::annotate
   +-> xic_scenarios          -> xic::heavy / xic::extract
   +-> import_pipeline        -> result-import::custom_json
   +-> entrapment_integration -> entrapment-analysis
   +-> search_pipeline        -> spectrum-io + param-recommend + search-engine + report + fdr
```

## 2. 测试文件清单

| 文件 | 测试数 | 覆盖场景（一行） |
|---|---|---|
| tests/annotation_scenarios.rs | 5 | DDA/DIA x 非SILAC/SILAC 的轻链/重链注释 + 无 K/R 零偏移边界 |
| tests/xic_scenarios.rs | 6 | DDA 按前体 m/z、DIA 按窗口包含找重链 scan + 零/非零偏移重链离子 |
| tests/import_pipeline.rs | 3 | custom JSON 导入解析 PSM、空数组、缺文件报错 |
| tests/search_pipeline.rs | 1 | 读谱->推参->搜索->FDR 汇总->注释最佳 PSM->导出 TSV 全链 |
| tests/entrapment_integration.rs | 7 | L0-L4 分级 + V2 字段 + V3 provenance / 23 列 TSV / mod 解析 / 配置兼容 |
| src/lib.rs（test_helpers 单测） | 5 | helper 自检：理论碎片数、合成峰非空、SILAC 值、窗口宽窄 |

测试函数名（与源码一致）：annotation 为 `scenario_1_dda_no_silac_annotates_light_only` ... `scenario_4_dia_silac_heavy_annotation_succeeds` 加 `zero_offset_peptide_no_kr_skips_heavy`；xic 为 `dda_find_heavy_scan_by_precursor_mz` / `dda_no_match_returns_none` / `dia_find_heavy_scan_by_window_containment` / `dia_no_window_contains_heavy_mz` / `zero_offset_no_heavy_target_ions_shift` / `nonzero_offset_heavy_ions_shifted`；import 为 `import_custom_json_parses_psms` / `import_empty_json_array` / `import_nonexistent_file_errors`；search 为 `full_pipeline_read_search_annotate`；entrapment 为 `test_known_peptide_classifications` / `test_v2_fields_end_to_end` / `test_v3_provenance_columns_in_tsv` / `test_v3_mod_parser_integration` / `test_v3_provenance_trace_known_peptides` / `test_v3_classified_psm_provenance_roundtrip` / `test_v3_config_backward_compatible`。

注释与 XIC 两组刻意用 2x2 矩阵（DDA/DIA 各配非 SILAC 与 SILAC）外加无 K/R 的零偏移边界，覆盖隔离窗宽窄与重链质量位移两条正交维度；entrapment 一组按 V1 分级、V2 替换字段、V3 provenance 三代演进逐层叠加，并各自核对 TSV 列数与 YAML 向后兼容。

## 3. fixtures 与缺失行为

本 crate 自身无 fixtures 数据（`tests/helpers/` 目录存在但为空，无共享模块文件）；端到端测试复用两处外部样本，两者缺失行为截然不同：

- `crates/search-engine/tests/fixtures_e2e/`：`test_100.mgf`（约 72 KB，100 张 MS2）加 `test_100.fasta`（约 33 KB），由 `search_pipeline.rs` 经 `search_engine_fixtures().parent()/fixtures_e2e` 定位。这两个文件已入 git（fresh checkout 必有），契约是硬断言：

```rust
assert!(mgf.exists() && fasta.exists(),
        "required fixtures missing: {mgf:?} / {fasta:?}");   // 缺则 test FAIL, 非跳过
```

- 工作区根 `.proteincopilot/databases/human_swissprot.fasta`（约 13 MB，`CARGO_MANIFEST_DIR` 上溯两级；未入 git、被 gitignore）。7 个 entrapment 测试中 3 个依赖它（`test_known_peptide_classifications` / `test_v2_fields_end_to_end` / `test_v3_provenance_columns_in_tsv`），缺失时软跳过（非 `#[ignore]`，提前 `return` 仍记 passed）：

```rust
if !fasta.exists() {
    eprintln!("SKIP: {} not found ...", fasta.display());
    return;                 // 软跳过, 仍 passed
}
```

其余 4 个 entrapment 测试与 import / xic / annotation / 单测全部用合成数据，零外部 I/O。这种"小样本入仓、大数据库下载"的切分让离线 CI 始终能跑：核心全链有确定样本护航，而依赖大库的同源性分级在无库时自动让路，绝不把缺数据误报成失败。

## 4. test_helpers 共享构造器（src/lib.rs）

| 函数 | 作用 |
|---|---|
| `silac_label()` | `LabelType::Silac { heavy_k_delta: 8.014199, heavy_r_delta: 10.008269 }` |
| `default_frag_tolerance()` | 20 ppm `MassTolerance` |
| `dda_window(mz)` / `dia_window(mz)` | 窄 0.35 Th / 宽 12.5 Th 隔离窗 |
| `make_ms2(scan, rt_min, mz, charge, iso, peaks)` | 合成 MS2（峰按 m/z 排序，前体强度 1e6） |
| `make_ms1(scan, rt_min, peaks)` | 合成 MS1 |
| `theoretical_fragments(peptide)` | 1+ b/y 离子 `IonEntry=(u32,f64)`，调 search-engine::chemistry |
| `synthetic_peaks_for_peptide(pep, base)` | 约半数 b/y 加 3 个噪声峰 |
| `heavy_shifted_peaks / search_engine_fixtures` | 近似重链移位 / fixtures 路径基准 |

临时目录用 `tempfile::tempdir()` 直接构造（import / search 导出 / entrapment TSV）。注意：annotation 的重链场景不用近似的 `heavy_shifted_peaks`，而在测试内自定义 `exact_heavy_peaks()` 产精确可匹配的重链峰，以保证 `matched_ions > 0`。

## 5. 端到端骨架（简化）

`search_pipeline.rs` 的六步全链：

```rust
// full_pipeline_read_search_annotate (#[tokio::test])
let info   = detect_format(&mgf).unwrap();              // 1 spectrum-io
let reader = create_reader(&info);                      // 测试用 create_reader (非缓存)
let summary = reader.read_summary(&mgf).unwrap();
let mut params = ParamRecommender.recommend(&summary, None).unwrap().decision;  // 2 推参
params.database_path = fasta.to_string_lossy().to_string();
let result = SimpleSearchEngine::new()                  // 3 搜索
    .search(&params, std::slice::from_ref(&mgf), noop_progress(),
            &mut SearchDiagnostics::new()).await.unwrap();
let fdr = ReportGenerator::generate_summary(&result);   // 4 fdr/report, id_rate in [0,1]
let best = result.psms.iter()
    .max_by(|a, b| a.score.partial_cmp(&b.score).unwrap()).unwrap();
let spec = reader.read_spectrum(&mgf, best.spectrum_scan).unwrap();   // O(1) seek
let ann  = annotate_spectrum(&spec, &best.peptide_sequence, best.charge,
                             &tol, &best.modifications, vec![], false, false).unwrap();  // 5
assert!(ann.matched_ions > 0);
ReportGenerator::export_tsv(&result, dir.path()).unwrap();           // 6 psm/peptide/protein.tsv
```

仅靠 test_helpers 即可独立断言注释（无需 fixtures）：

```rust
// scenario_1_dda_no_silac_annotates_light_only
let peaks = synthetic_peaks_for_peptide("PEPTIDEK", 1000.0);
let spec  = make_ms2(1, 30.0, 458.24, 2, Some(dda_window(458.24)), peaks);
let ann = annotate_spectrum(&spec, "PEPTIDEK", 2, &default_frag_tolerance(),
                            &[], vec![], false, false).unwrap();
assert!(ann.matched_ions > 0);
assert!(ann.heavy_annotation.is_none());   // 非 SILAC -> 无重链注释
```

六步串起五个 crate：前两步把读谱摘要交给确定性推参，中段搜索产出 PSM 后立即过 FDR 汇总，末段回头按 `spectrum_scan` 做 O(1) 取谱、注释最佳命中并落盘三张 TSV，正是真实工具链 `run_search` 之后的典型路径。

## 6. 运行

```
cargo test -p protein-copilot-integration-tests --offline
```

实测 27 个全过：test_helpers 单测 5 + annotation 5 + entrapment 7 + import 3 + search 1 + xic 6（doctests 0）。依赖关系：

- `search_pipeline`（1 个）依赖在仓的 `fixtures_e2e/test_100.mgf|fasta`，缺则硬断言失败；
- 3 个 entrapment 测试依赖未入仓的 `human_swissprot.fasta`，缺则软跳过仍 passed；
- 其余 23 个全用合成数据，离线零外部依赖。

运行时间集中在 entrapment（建 FASTA 索引约 95 s）与 search_pipeline（约 8 s），其余近乎瞬时。

---

回到 [README](README.md) 选择其它层级或子系统。
