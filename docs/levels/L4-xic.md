# L4 -- xic crate（碎片离子 XIC 提取 + Plotly 可视化 + SILAC）

承接 [L3-xic-dia](L3-xic-dia.md) 的宏观视角，上溯 [L2](L2-architecture.md)。本篇只看 `crates/xic` 一个 crate——它的对外函数签名、数据结构字段与核心代码骨架。所有签名/常量以 `crates/xic/src/{lib,extract,heavy}.rs` 源码为准。

## 1. 用途与位置

输入：一份 PSM 上下文（肽序列 + 电荷 + 真实前体 m/z + 目标 MS2 scan + 修饰）加一份 mzML 谱图（经 `SpectrumReader` 定点读）。输出：`XicUnifiedResult`——随保留时间变化的前体（MS1）与 b/y 碎片（MS2）强度曲线，外加原始峰阵列与离子元数据，供 `report` 渲染成自包含 Plotly.js HTML，让人眼核验共洗脱碎片是否齐步起落、轻重峰比是否合理。

纯计算 crate，无 MCP/网络/LLM，只读不改谱图。依赖三层：`core`（`Spectrum`/`MassTolerance`/`LabelType` 等共享类型与 SILAC delta）、`spectrum-io`（`SpectrumReader::{list_scan_meta, read_spectrum}`）、`search-engine`（复用 `generate_b/y_ions_with_charge`、`within_tolerance`，不重复实现碎片计算）。算法是 1.5 趟：Pass 0 读目标 scan 取 RT 与隔离窗，Pass 1 只读规划好的那批 scan。

## 2. 对外 API

唯一入口（旧 `extract_xic` / `extract_xic_with_raw` 已删除，只保留它）：

```rust
#[allow(clippy::too_many_arguments)]
pub fn extract_xic_unified(
    reader: &dyn protein_copilot_spectrum_io::SpectrumReader,
    file_path: &Path,
    target_scan: u32,
    peptide_sequence: &str,
    charge: i32,
    precursor_mz: f64,
    modifications: &[Modification],
    params: &ExtractionParams,
    ms1_mz_window_da: f64,
) -> Result<crate::XicUnifiedResult, XicError>
```

返回 `XicUnifiedResult { xic_data, raw_scans, ion_metadata }`：

- `xic_data: XicData` -- 曲线本体；
- `raw_scans: RawScanData { ms1_scans, ms2_scans, ms2_heavy_scans }` -- 窗口内原始峰阵列，供浏览器端 SILAC 重算（`ms2_heavy_scans` 仅 DIA+SILAC 非空）；
- `ion_metadata: Vec<IonMetadataEntry>` -- 每条 top-N 离子的 K/R 计数 + `light_mz`。

`XicData` 字段：`peptide_sequence`、`target_rt_min`、`target_scan`、`charge: i32`、`precursor_mz`、`ms1_precursor_xic: Option<XicTrace>`、`ms1_heavy_precursor_xic`、`fragment_xic_traces: Vec<XicTrace>`、`heavy_fragment_xic_traces`、`extraction_params`、`heavy_warning: Option<String>`。

`XicTrace`：`ion_label`、`ion_type: IonType`、`ion_number`、`charge: u32`、`theoretical_mz`、`data_points: Vec<XicDataPoint>`、`is_heavy: bool`。`XicDataPoint`：`retention_time_min`、`scan_number`、`intensity`、`observed_mz: Option<f64>`（无峰即 `0.0` + `None`）。

`ExtractionParams`：`mz_tolerance: MassTolerance`、`n_cycles: u32`、`top_n_ions: usize`、`label_type: Option<LabelType>`、`intensity_rule: IntensityRule`。`IntensityRule` 三档（`#[default] MaxInWindow`、`SumInWindow`、`NearestPeak`）决定容差窗内如何取强度。

后处理：碎片轨迹按各自数据点强度之和降序排序，截到 `top_n_ions` 条，再剔除全零轨迹；heavy 轨迹只对留下的 top label 重建。`heavy_warning` 在带 K/R 标签却没抽到任何 heavy MS2 时给出——DIA 提示 heavy 前体落在所有窗口之外，DDA 提示目标 RT 附近找不到选中 heavy 前体的 scan。

## 3. 关键逻辑

- 定点读规划：先 `reader.list_scan_meta(file_path)` 从内存索引拿到全部 scan 元数据（亚毫秒、零额外 I/O），据此规划三组目标 scan——同隔离窗的 light MS2（以 `target_scan` 为中心取前后各 `n_cycles` 个循环）、heavy MS2（DIA+SILAC 时按 heavy 前体 m/z 命中宽窗）、覆盖上述 RT 跨度的 MS1——再逐个 `read_spectrum` 做 O(1) seek。heavy MS2 在宽窗命中后还会挑 RT 最接近目标 scan 的那张作中心，再向两侧扩 `n_cycles`；MS1 则取落在全部已规划 MS2 的 RT 跨度内的所有 MS1 scan。把全文件流式遍历降为约 30 次定点读，单条注释从约 240s 降到 <1s。
- 碎片离子构建：`build_target_ions` 复用 search-engine 的生成器，碎片电荷上限由前体电荷决定——`max_frag_charge = if precursor_charge >= 3 { 2 } else { 1 }`，即 3+ 及以上前体才额外生成 2+ 碎片，离子标签按电荷附 1+/2+ 上标。
- `same_isolation_window` 判据：两窗中心差 < 1.0 Da 且宽度相对差 < 0.2（20%），且宽度 > 0，二者皆成立才算同一 DIA 窗。DIA 检测看 `lower_offset + upper_offset > 1.0`。
- SILAC heavy mirror：轻肽与重肽共洗脱但前体 m/z 不同，重碎片不能从轻 scan 复用。`compute_heavy_target_ions` 按离子类型取前缀（b）/后缀（y）数 K/R，乘 `residue_heavy_delta` 再除以碎片电荷得 heavy m/z。DIA 走 heavy 前体所在宽窗；DDA 走 `find_dda_heavy_scan`——在候选 MS2 里挑前体 m/z 落在 ppm 容差内、且 RT 最近目标的那张 scan。重标定位失败时填 `heavy_warning` 软告警而非报错。
- 三种取强度规则共用一次窗口扫描：`MaxInWindow` 取窗内最高峰（观测 m/z 为该峰）、`SumInWindow` 求和（观测 m/z 取最近峰）、`NearestPeak` 取最近峰；无峰一律返回 `(0.0, None)`。
- top-N 与轻重配对：light 碎片先排序截断，heavy 轨迹按相同 `ion_label` 集合筛出，保证轻重逐离子对应；MS1 前体曲线分别在 `precursor_mz` 与 heavy 前体 m/z 上抽强度，得 `ms1_precursor_xic` 与 `ms1_heavy_precursor_xic`。
- `label_type` 经 `total_heavy_delta(...).abs() > 1e-6` 过滤：无 K/R 的肽即便带标签也跳过 heavy。MS1 trim 半窗动态放大到 `max(ms1_mz_window_da, max_heavy_shift + 5.0)`。

## 4. 简化源码片段

extract_intensity——`partition_point` 二分定位 + 1.5 倍容差扩边：

```rust
let pos = exp_mz.partition_point(|&m| m < target_mz);
let max_da = match tolerance.unit {
    ToleranceUnit::Ppm => target_mz * tolerance.value * 1e-6,
    ToleranceUnit::Da  => tolerance.value,
};
// 自 pos 向两侧各扩 1.5 倍容差，圈定 [scan_start, scan_end)
for i in scan_start..scan_end {
    if within_tolerance(exp_mz[i], target_mz, tolerance) { /* max / sum / nearest */ }
}
match rule {
    IntensityRule::MaxInWindow => (best_intensity, Some(best_mz)),
    IntensityRule::SumInWindow => (sum_intensity, Some(nearest_mz)),
    IntensityRule::NearestPeak => (nearest_intensity, Some(nearest_mz)),
}
```

extract_xic_unified 骨架（规划 + 定点读）：

```rust
let target = reader.read_spectrum(file_path, target_scan)?;   // Pass 0: RT + 隔离窗
let light_ions = build_target_ions(peptide_sequence, modifications, charge);
let all_meta = reader.list_scan_meta(file_path)?;             // 内存索引，零额外 I/O
let light_ms2 = all_meta.iter()
    .filter(|m| m.ms_level == 2 && same_isolation_window(tw, &meta_window))
    .collect();                                               // 再以 target 为中心切 +/- n_cycles
for (scan, _rt) in planned_light_ms2 {
    let spec = reader.read_spectrum(file_path, scan)?;        // O(1) seek
    // 对每条 light_ion 调 extract_intensity ...
}
// heavy MS2 + MS1 同理；最后按总强度降序排序、截 top_n、建 XicTrace
```

heavy 计算（按 b 前缀 / y 后缀数 K/R，再除以碎片电荷）：

```rust
let fragment_delta = match ion.ion_type {
    IonType::B => residue_heavy_delta(&chars[..(ion.ion_number as usize).min(n)], label),
    IonType::Y => residue_heavy_delta(&chars[n.saturating_sub(ion.ion_number as usize)..], label),
    IonType::Precursor => total_heavy_delta(peptide_sequence, label),
};
let heavy_mz = ion.mz + fragment_delta / ion.charge.max(1) as f64;
```

标准 SILAC delta：K +8.014199 Da、R +10.008269 Da（`LabelType::standard_silac`）。

## 5. 调用链与 report

```
mcp-server tool extract_xic / annotate_spectrum
  -> get_or_create_reader(path)            // LRU 缓存 IndexedMzMLReader
  -> xic::extract::extract_xic_unified(..) // 本 crate
  -> report::render_xic / render_unified / render_xic_3d
       -> XicData / RawScan + PlotlyMode -> 自包含 Plotly.js HTML
```

- `extract_xic` 工具：默认传 `mz_tolerance` 20 ppm、`n_cycles` 5、`ms1_mz_window_da` 20.0；2D 走 `render_xic`，3D 走 `render_xic_3d`（用 `RawScanData.ms2_scans` 重建 3D 峰面）。
- `annotate_spectrum` 工具：mzML 输入时顺带跑 `extract_xic_unified`，把 XIC 并入统一注释视图 `render_unified`，并据 MS1 观测峰回修前体 m/z。
- report 侧只消费 `XicData` / `RawScan` / `IonMetadataEntry` / `PlotlyMode`，不反向依赖 xic 的提取逻辑——单向依赖，便于各自演进。同一文件被反复注释时，索引与元数据在 MCP 层 LRU 缓存（容量 8）里复用，省去重复建 `.mzML.idx`。

## 6. 测试入口

```
cargo test -p protein-copilot-xic --offline
```

覆盖 33 个单元测试：`extract.rs` 测 `extract_intensity`（max/sum/nearest/ppm/空谱）、`same_isolation_window`、`build_target_ions`、`compute_ion_metadata`、`trim_peaks_to_window`，以及合成 `.pfb` 的端到端 `pfb_end_to_end_xic_via_indexed_reader` 与 `dda_silac_missing_heavy_scan_emits_warning`；`heavy.rs` 测 `compute_heavy_precursor_mz`、`compute_heavy_target_ions`、`find_dda_heavy_scan`、`find_heavy_dia_window_from_spectra`、`window_contains_mz`。

---

继续阅读：[README](README.md) 总目录 -- 上溯 [L2 架构](L2-architecture.md)、[L3 XIC + DIA](L3-xic-dia.md)。
