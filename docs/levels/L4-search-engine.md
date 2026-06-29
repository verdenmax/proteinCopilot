# L4 - search-engine（模块级 API 与代码骨架）

承接 [L3-search-engine](L3-search-engine.md) 与 [L2](L2-architecture.md)。L3 已讲清"酶切 -> 候选 -> 理论谱 -> 容差打分 -> target-decoy -> 聚合"的整体流程，本篇不复述，只聚焦 `crates/search-engine`（package `protein-copilot-search-engine`）的**模块级 API 表面**与**代码骨架**，所有签名与常量逐一核对源码。

## 1. 用途 + 位置 + 依赖

- 位置：`crates/search-engine`，package 名 `protein-copilot-search-engine`。
- 用途：酶切 / 匹配 / 打分 / 引擎调度 —— 把 `SearchParams` + 谱图文件确定性地变成统一的 `SearchResult`，并用 `SearchEngineAdapter` trait 抽象 SimpleSearch / Sage / pFind 多引擎。

| 类别 | crate |
|------|-------|
| workspace 运行依赖 | `protein-copilot-core`、`protein-copilot-spectrum-io`、`protein-copilot-fdr` |
| 外部运行依赖 | `sage-core`（git pin `rev cd712d4`，v0.15.0-beta.2）、`rayon = "1"` |
| 基础设施 | `tokio`、`async-trait`、`thiserror`、`serde`、`schemars`、`uuid`、`chrono`、`serde_json`、`tracing` |
| dev-dependencies | `protein-copilot-param-recommend`、`protein-copilot-report`、`tempfile = "3"` |

注（核对 Cargo.toml）：`param-recommend` 与 `report` 只出现在 `[dev-dependencies]`，运行期不链接；`spectrum-io` 同时被列入运行与 dev 两处。

## 2. 模块级 API 表

| 文件 | 关键导出（签名要点） | 职责一句话 |
|------|------|------|
| `lib.rs` | re-export `EngineRegistry` `SimpleSearchEngine` `SageAdapter` `SearchProgress` `SearchEngineError` | crate 公共入口 |
| `chemistry.rs` | `residue_mass(char) -> Option<f64>`；`peptide_mass(&str) -> Option<f64>`；`peptide_mz(f64, i32) -> f64`；`is_standard_sequence(&str) -> bool`；`const WATER_MASS = 18.010565`、`PROTON_MASS = 1.007276` | 残基单同位素质量与 m/z 计算 |
| `digest.rs` | `digest(seq, acc, &Enzyme, missed: u32)`（默认长度 6..=50）；`digest_with_length(.., min: u32, max: u32) -> Vec<DigestedPeptide>`；私有 `find_cleavage_sites(&[char], &Enzyme) -> Vec<usize>`、`split_at_sites` | 体外酶切（char 索引，UTF-8 安全） |
| `matching.rs` | `within_tolerance(obs, theo, &MassTolerance) -> bool`；`match_spectrum(..) -> Option<PeptideMatch>`；`match_spectrum_all(..) -> Vec<PeptideMatch>`；`generate_b_ions` `generate_y_ions(&str, &[Modification]) -> Vec<f64>`；`pub(crate) mod_delta_fragment(..) -> f64` | 前体容差匹配 + b/y 理论谱 + 计数打分 |
| `varmod.rs` | `find_applicable_sites(seq, &[Modification], nterm, cterm) -> Vec<(usize, Vec<usize>)>`；`enumerate_combinations(&[Modification], sites, max: u32) -> Vec<ModCombination>`；`is_terminal_pos(usize) -> bool`；`const NTERM_POS = usize::MAX`、`CTERM_POS = usize::MAX - 1` | 可变修饰位点发现与组合枚举 |
| `annotate.rs` | `annotate_spectrum(..) -> Result<SpectrumAnnotation, SearchEngineError>`；`annotate_heavy_spectrum(..) -> Result<HeavyAnnotation, _>`；类型 `IonType{B,Y}` `IonAnnotation` `AnnotatedPeak` `TheoreticalIon` `SpectrumAnnotation`（均 serde + JsonSchema） | 单谱注释（可视化 / 质检） |
| `registry.rs` | `EngineRegistry::{new, register(Box<dyn SearchEngineAdapter>), get(&str) -> Option<&dyn ..>, list_available() -> Vec<EngineInfo>, health_check_all().await, len, is_empty}` | 按名注册 / 查找 adapter |
| `simple_engine.rs` | `SimpleSearchEngine`（ZST）+ `new()`；impl `SearchEngineAdapter`；私有 `run_search` `run_search_on_spectra` `collect_psms_for_spectrum` `finalize_search_result` | 自研 MVP 引擎，编排全流程 |
| `fasta.rs` | `parse_fasta(&Path) -> Result<Vec<FastaEntry>, SearchEngineError>`；`FastaEntry{accession, description, sequence}` | FASTA 解析 |
| `adapters/sage/mod.rs` | `SageAdapter::new(thread_count: usize)`；impl `SearchEngineAdapter`（含 `search_with_spectra` `cancel`）；私有 `feature_to_psm` | Sage 库集成编排 |
| `adapters/sage/config.rs` | `build_sage_parameters(&SearchParams) -> sage_core::database::Parameters` | SearchParams -> Sage Parameters 映射 |
| `adapters/sage/convert.rs` | `spectrum_to_raw`、`mass_tolerance_to_sage`、`fixed_mod_to_sage`、`variable_mod_to_sage` | 谱图 / 容差 / 修饰类型转换 |
| `adapters/pfind.rs` | `PFindAdapter::new(SshConfig)`、`SshConfig`（stub，方法返回未实现错误） | pFind 远程执行预留 |
| `error.rs` | `enum SearchEngineError`：`InvalidParams` `FastaError` `IoError` `EngineNotFound` `ExecutionError` `NoInputSpectra` + `From<_> for CoreError` | 错误枚举与转换 |

数据结构 `DigestedPeptide` `PeptideMatch` `ModSite` `ModCombination` 的字段定义见对应文件，在下方骨架中出现。

模块依赖自下而上分层：`chemistry` 提供质量常量与计算基元；`digest` / `varmod` / `matching` / `annotate` 在其上构建候选肽、修饰组合、理论谱与打分；`fasta` 负责输入解析；`simple_engine` 把这些编排成一次完整搜索；`registry` 与 `adapters/*` 则把不同引擎统一到 `SearchEngineAdapter` 抽象之后。错误统一经 `error::SearchEngineError` 向 `core::CoreError` 收敛。

## 3. 引擎实现对比（engine_info 真实值）

| Adapter | name | version | supported_features | health_check |
|---------|------|---------|--------------------|--------------|
| `SimpleSearchEngine` | `"SimpleSearch"` | `"0.1.0"` | `["basic_search", "b_y_scoring"]` | `Healthy` |
| `SageAdapter` | `"Sage"` | `"0.15.0"` | `["open_search", "lfq", "tmt", "chimera"]` | `Healthy` |
| `PFindAdapter` | `"pFind"` | `"3.x (not connected)"` | `["open_search", "modification_localization"]` | `Unavailable { reason }` |

三者都实现同一 `core::engine::SearchEngineAdapter`；只有 SimpleSearch 与 Sage 真正覆写 `search_with_spectra`，pFind 各方法均返回 `CoreError::SearchEngineError`（其 `search_with_spectra` 走 trait 默认实现报错）。

## 4. 简化源码片段

以下三段分别展示引擎接缝、确定性酶切与打分内核；为突出主干，省略了日志、诊断与进度上报。

**(a) `SearchEngineAdapter` 实现骨架（`simple_engine.rs`）**

```rust
#[async_trait::async_trait]
impl SearchEngineAdapter for SimpleSearchEngine {
    async fn search(&self, params: &SearchParams, input_files: &[PathBuf],
                    on_progress: ProgressCallback, diag: &mut SearchDiagnostics)
        -> Result<SearchResult, CoreError> {
        self.run_search(params, input_files, &*on_progress, diag)
            .map_err(CoreError::from)            // SearchEngineError -> CoreError
    }
    fn engine_info(&self) -> EngineInfo {
        EngineInfo { name: "SimpleSearch".to_string(), version: "0.1.0".to_string(),
            supported_features: vec!["basic_search".to_string(), "b_y_scoring".to_string()] }
    }
    async fn health_check(&self) -> Result<HealthStatus, CoreError> { Ok(HealthStatus::Healthy) }
}
```

**(b) 酶切骨架（`digest.rs`，非 NonSpecific 路径）**

```rust
pub fn digest_with_length(seq, acc, enzyme, missed, min, max) -> Vec<DigestedPeptide> {
    let chars: Vec<char> = seq.chars().collect();      // char 索引 = UTF-8 安全
    let sites = find_cleavage_sites(&chars, enzyme);   // 例: Trypsin 在 K 或 R 后, 但不在 P 前
    let frags = split_at_sites(&chars, &sites);
    let mut out = Vec::new();
    for mc in 0..=missed as usize {                    // 0..=漏切数
        for (i, w) in frags.windows(mc + 1).enumerate() {
            let pep: String = w.concat();
            let len = pep.chars().count() as u32;
            if len >= min && len <= max {
                if let Some(mass) = peptide_mass(&pep) {   // 含非标准残基 -> None, 丢弃
                    out.push(DigestedPeptide { sequence: pep, neutral_mass: mass,
                        protein_accession: acc.to_string(),
                        is_protein_nterm: i == 0,
                        is_protein_cterm: i + mc + 1 == frags.len() });
                }
            }
        }
    }
    out
}
```

**(c) 容差判定 + b 离子生成（`matching.rs`）**

```rust
pub fn within_tolerance(obs: f64, theo: f64, tol: &MassTolerance) -> bool {
    match tol.unit {
        ToleranceUnit::Ppm => theo > 0.0 && ((obs - theo) / theo).abs() * 1e6 <= tol.value,
        ToleranceUnit::Da  => (obs - theo).abs() <= tol.value,
    }
}
// b 离子: 沿前缀累加残基质量, 叠加修饰偏移, 逐电荷 z 输出 m/z
let mut cumulative = 0.0;
for &aa in &chars[..n - 1] {
    cumulative += residue_mass(aa)?;                       // 非标准残基 -> 返回空 Vec
    let neutral = cumulative + mod_delta + cumulative_var; // 固定 + 位点化可变修饰偏移
    for z in 1..=max_z { ions.push((neutral + z as f64 * PROTON_MASS) / z as f64); }
}
```

匹配细节：前体若未带电荷，`match_spectrum` 依次尝试 `z = 2, 3, 1, 4`；当 `z >= 3` 时片段离子额外生成 2+（`max_frag_charge`）。`digest` 对 `NonSpecific` 走特例——直接枚举落在 `[min, max]` 内的全部子串，而非按漏切窗口截断，避免长度被 `mc + 1` 限死。

## 5. 调用链

SimpleSearch（`run_search` 与 `search_with_spectra` 共用 `run_search_on_spectra`）：

```text
SearchParams + input_files
  -> params.validate()
  -> parse_fasta() -> Vec<FastaEntry>
  -> digest_with_length() 逐蛋白 -> Vec<DigestedPeptide>
  -> 若 decoy 开启: protein_copilot_fdr::generate_decoys()
  -> 逐 MS2: collect_psms_for_spectrum()
       -> match_spectrum() | match_spectrum_all()
            -> find_applicable_sites() -> enumerate_combinations()
            -> within_tolerance() -> generate_b_ions_positional + generate_y_ions_positional
            -> count_matched_ions()
  -> protein_copilot_fdr::calculate_fdr() -> psm.q_value
  -> finalize_search_result() -> SearchResult (PSM + 肽 + 蛋白 + summary + metadata)
```

`collect_psms_for_spectrum` 按数据类型分支：多前体、或宽分离窗（`lower_offset + upper_offset` 之和 > 5 Da）判为 DIA，走 `match_spectrum_all` 收集每个前体的最佳匹配；否则按 DDA 走 `match_spectrum` 取单条最佳。`finalize_search_result` 在有 q-value 时以 1% FDR 过滤再汇总 PSM / 肽 / 蛋白计数。

Sage（`adapters/sage/mod.rs`，4 个 phase）：

```text
search() -> create_indexed_reader().read_all() -> search_with_spectra()
  -> Phase1 过滤 MS2 + spectrum_to_raw()
  -> Phase2 build_sage_parameters() -> IndexedDatabase + Scorer
            rayon into_par_iter().map(scorer.score)
  -> Phase3 LDA score_psms -> spectrum_q_value -> picked_peptide + picked_protein
  -> Phase4 feature_to_psm() -> SearchResult
```

## 6. 测试入口

```bash
cargo test -p protein-copilot-search-engine --offline
```

| 测试二进制 | 通过数 |
|-----------|--------|
| `unittests src/lib.rs` | 118 |
| `tests/e2e_integration.rs` | 8 |
| `tests/integration.rs` | 6 |
| `tests/sage_integration.rs` | 3 |
| doc-tests | 0 |
| **合计** | **135** |

返回目录 [README](README.md)。
