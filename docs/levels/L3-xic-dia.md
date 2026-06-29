# L3 — XIC 提取 + DIA 子系统

承接 [L2](L2-architecture.md)。本篇讲清两个互补的库 crate：`xic`（色谱提取，供注释与可视化）与 `dia-extraction`（DIA 前体提取，供搜索）。两者都是纯计算 crate——无 MCP、无网络、无 LLM。源码以 `crates/xic/src/{lib,extract,heavy}.rs`、`crates/dia-extraction/src/{lib,correlation,detection,isotope}.rs`、`crates/core/src/{label,spectrum}.rs` 为准。

## 1. 职责与位置

- **`xic`**：给定一个 PSM（肽序列、电荷、真实前体 m/z、目标 scan），从 mzML 抽出**前体（MS1）+ b/y 碎片（MS2）**的随保留时间变化的强度曲线。可选 SILAC 重标，输出轻/重两套对照色谱。它被注释类 MCP 工具调用，结果交给 `report` 渲染 Plotly.js HTML，供人眼核验鉴定质量——共洗脱的碎片是否齐步起落、轻重峰强度比是否合理，是判断一条鉴定真伪的直观证据。它只读不改谱图，依赖 `search-engine` 复用碎片离子计算，依赖 `spectrum-io` 读盘。
- **`dia-extraction`**：DIA 谱图的宽隔离窗里混着多个共碎裂前体，搜索引擎无法直接用「一谱一前体」的假设去搜。本 crate 从 MS1 同位素模式反推出**候选前体（m/z + 电荷）**，写回每张 MS2 的 `precursors`，再缓存。MCP 工具 `extract_dia_precursors` 返回 `run_id`，随后 `run_search(dia_run_id=...)` 从缓存取增强谱再搜——这是 DIA 进入打分流程的旁路（见 L2 第 3 节）。另有 `extract_single_spectrum_precursors` 供单张 MS2 的探查式提取，附带它用了哪张 MS1、用何种关联方式。

## 2. 模块边界

| crate | 文件 | 职责 |
|-------|------|------|
| xic | `lib.rs` | 公共数据结构（`XicUnifiedResult`/`XicData`/`XicTrace`）+ `IonType`/`IntensityRule`/`PlotlyMode` |
| xic | `extract.rs` | `extract_xic_unified` 主流程、`extract_intensity`、`same_isolation_window`、`build_target_ions`、`trim_peaks_to_window` |
| xic | `heavy.rs` | SILAC 重标镜像：`compute_heavy_target_ions`、`find_dda_heavy_scan`、`window_contains_mz` |
| dia | `lib.rs` | `extract_dia_precursors` / `extract_single_spectrum_precursors` 编排 |
| dia | `detection.rs` | `detect_acquisition_mode`（阈值判 DDA/DIA）、`separate_by_ms_level` |
| dia | `correlation.rs` | `correlate_ms1_ms2`：每张 MS2 关联到一张 MS1 |
| dia | `isotope.rs` | `IsotopePatternExtractor`：同位素簇 -> 候选前体 |

> 注意：旧的 `extract_xic` / `extract_xic_with_raw` 已删除，统一为 `extract_xic_unified` 一次产出色谱 + 原始峰 + 离子元数据。

## 3. 关键数据结构

色谱以「结果 -> 多条 trace -> 多个点」三层组织（`xic/lib.rs`，简化）：

```rust
pub struct XicUnifiedResult {
    pub xic_data: XicData,               // 色谱
    pub raw_scans: RawScanData,          // 窗内原始峰（供浏览器端重算 SILAC）
    pub ion_metadata: Vec<IonMetadataEntry>, // 每个碎片的 K/R 计数
}
pub struct XicData {
    pub precursor_mz: f64,               // 真实前体 m/z（来自 PSM，非 DIA 窗中心）
    pub ms1_precursor_xic: Option<XicTrace>,
    pub ms1_heavy_precursor_xic: Option<XicTrace>,
    pub fragment_xic_traces: Vec<XicTrace>,       // 轻
    pub heavy_fragment_xic_traces: Vec<XicTrace>, // 重
    pub heavy_warning: Option<String>,   // 重标 MS2 缺失时的说明
}
pub struct XicTrace { pub theoretical_mz: f64, pub data_points: Vec<XicDataPoint>, pub is_heavy: bool, /* ion_label/type/number/charge */ }
pub struct XicDataPoint { pub retention_time_min: f64, pub scan_number: u32, pub intensity: f64, pub observed_mz: Option<f64> }
```

SILAC 重标类型在 `core/label.rs`，按 K/R 残基计移量：

```rust
pub enum LabelType {
    Silac { heavy_k_delta: f64, heavy_r_delta: f64 }, // 标准 K+8.014199, R+10.008269
    Custom { residue_deltas: Vec<(char, f64)> },
}
// total_heavy_delta(seq, label) 累加序列内 K/R 移量；
// compute_heavy_precursor_mz = light_mz + delta / charge
```

隔离窗在 `core/spectrum.rs`——中心加上下偏移，DDA 窄（1-3 Da）DIA 宽（如 25 Da）：

```rust
pub struct IsolationWindow { pub target_mz: f64, pub lower_offset: f64, pub upper_offset: f64 }
// 覆盖区间 = [target_mz - lower_offset, target_mz + upper_offset]
```

## 4. XIC 提取伪代码（extract_xic_unified）

核心思想：**不全文件流式扫描**，而是先用内存索引 `list_scan_meta`（亚毫秒）规划出 ~30 个 scan，再用 `read_spectrum` 做 O(1) 定点 seek——把大 DIA 文件单次注释从 ~240s 降到 <1s。

```text
fn extract_xic_unified(reader, file, target_scan, peptide, charge, precursor_mz, mods, params, ms1_window_da):
  # 1 读目标 scan（必须 MS2），取 RT 与隔离窗
  target = reader.read_spectrum(file, target_scan)
  target_window = target.precursors[0].isolation_window
  is_dia = (lower_offset + upper_offset) > 1.0

  # 2 生成 b/y 目标离子；按 SILAC 算重标镜像（仅当 total_heavy_delta > 1e-6）
  light_ions   = build_target_ions(peptide, mods, charge)   # 前体>=3 时含 2+ 碎片
  heavy_ions   = compute_heavy_target_ions(light_ions, peptide, label)  # b 取前缀/y 取后缀 K,R
  heavy_prec_mz = compute_heavy_precursor_mz(precursor_mz, charge, peptide, label)

  # 3 一次性取全部 scan 元数据（内存 ScanIndex，零额外 I/O）
  all_meta = reader.list_scan_meta(file)

  # 4 规划 light MS2：同一隔离窗 + 目标前后 +/-n_cycles
  light = [m in all_meta if m.ms_level==2 and same_isolation_window(target_window, m.window)]
  light = sort_by_scan(light)[target +/- n_cycles]            # 窗不存在(DDA无窗)则全收

  # 5 规划 heavy MS2（仅 DIA+SILAC）：宽窗包含 heavy_prec_mz，取 RT 最近为中心 +/-n_cycles
  heavy = plan_heavy_window(all_meta, heavy_prec_mz) if is_dia and heavy_ions else []

  # 6 规划 MS1：覆盖上面 MS2 的 RT 区间 [lo, hi]
  ms1 = [m in all_meta if m.ms_level==1 and lo <= m.rt <= hi]

  # 7 仅对这 ~30 个 scan 定点读 + 抽强度
  for s in light: ms2_light[s] = [extract_intensity(ion.mz, peaks, tol, rule) for ion in light_ions]
  if is_dia:   for s in heavy: ms2_heavy[s] = extract over heavy_ions   # DIA 走宽窗包含
  elif heavy_ions:                                                      # DDA+SILAC
      center = find_dda_heavy_scan(candidates, target_rt, heavy_prec_mz, ppm)  # 选重前体所在 MS2
      ms2_heavy = extract over center +/- n_cycles（守卫：该 scan 确实选中重前体）
  for s in ms1: light_int = extract_intensity(precursor_mz); heavy_int = extract_intensity(heavy_prec_mz)
                raw_ms1 = trim_peaks_to_window(peaks, precursor_mz, ms1_window_da)  # 削峰减体积

  # 8/9 组装：按总强度降序 -> 截 top_n -> 去全零 trace；重标只取 top 标签
  return XicUnifiedResult { xic_data, raw_scans, ion_metadata }
```

要点：`extract_intensity` 用 `partition_point` 二分定位再按 `IntensityRule`（`MaxInWindow`/`SumInWindow`/`NearestPeak`）取值，返回 `(intensity, observed_mz)`；`same_isolation_window` 以「中心差 < 1.0 Da 且窗宽相对差 < 20%」判定同窗。轻重肽共洗脱但前体 m/z 不同，故重碎片必须从重前体所在的 MS2 取，**绝不复用轻 scan**——DIA 用宽窗包含定位，DDA 用 `find_dda_heavy_scan` 按 ppm + 最近 RT 定位。重标离子找不到对应 MS2 时不报错，置 `heavy_warning` 说明。

几个工程细节：`n_cycles` 默认 5，故轻 MS2 约 11 张，叠加重 MS2 与窗内 MS1，总规划约 30 次定点读，正是把「全文件流式」换成「索引规划 + O(1) seek」的收益来源。MS1 原始峰用 `trim_peaks_to_window` 只留前体附近一段；窗宽取 `ms1_window_da` 与「最大重标移量 + 5」的较大者，保证浏览器端重算 SILAC 时轻、重两簇都落在窗内。碎片色谱先按总强度降序、截 `top_n`、再去掉全零 trace；重标 trace 只为入选 top 的标签生成，二者一一对应便于叠图比对。

## 5. DIA 前体提取伪代码（extract_dia_precursors）

```text
fn extract_dia_precursors(spectra, extractor, config):
  (ms1, ms2) = separate_by_ms_level(spectra)
  if ms2.empty: Err(NoMs2Spectra)

  # 检测采集模式：MS2 隔离窗宽 (lower+upper) 的中位数 > 阈值(默认 5 Da) -> DIA
  mode = config.acquisition_mode ?? detect_acquisition_mode(spectra, dia_threshold_da)
  if mode in {DDA, Unknown}: return 原样 MS2（0 提取）       # 旁路，不动谱图
  if ms1.empty: Err(NoMs1Spectra)

  # MS2 关联最近 MS1：三级回退 source_scan -> scan_order -> rt_nearest
  idxs = correlate_ms1_ms2(ms1, ms2)
  for (spec, idx) in zip(ms2, idxs):
      iw = spec.isolation_window
      spec.precursors = extractor.extract(ms1[idx], iw)    # 隔离窗内同位素分析 -> 候选前体
  return DiaExtractionResult { detected_mode, enhanced_spectra, stats }
```

`IsotopePatternExtractor::extract`（`isotope.rs`）在隔离窗内做同位素簇搜索，定候选前体**电荷**：

```text
peaks = 窗内峰按强度降序
for seed in peaks (未用):
  for z in min_charge..=max_charge:            # 默认 2..=5
      delta = 1.003355 / z                     # 13C-12C 间距 / 电荷
      向后找 M+0(seed-delta)，向前串 M+1,M+2...（容差 0.01 Da，强度不得 > 前峰 1.5x）
  选峰数最多的簇；并列时取更低电荷
  单同位素 m/z = 簇内最小 m/z；charge = z
去重(同电荷 0.01 Da 内留强者) -> 按强度降序
```

`correlate_ms1_ms2` 为每张 MS2 选一张 MS1，采用三级回退（`correlation.rs`）：(1) `source_scan`——mzML 里前体显式引用的 MS1 scan，最可靠；(2) `scan_order`——退而取「scan 号小于本 MS2 的最大 MS1」，即采集序上紧邻的前一张 MS1；(3) `rt_nearest`——再退而按保留时间绝对差最近。`extract_single_spectrum_precursors` 会把命中的级别名（`source_scan`/`scan_order`/`rt_nearest`）一并返回，便于诊断关联是否可信。

产物经 `into_pseudo_spectra`（多前体谱按前体拆成多张伪谱，单前体谱直接迁移）或 `into_enhanced_spectra`（多前体保留）落入 MCP 的 `dia_cache`（键为 `run_id`），**缓存供搜索**：后续 `run_search(dia_run_id=...)` 直接取增强谱，绕过「DIA 不能直接搜」的守卫。`ExtractionStats` 同时记录 MS1/MS2 数、提取前体总数、每张 MS2 平均前体数与电荷分布，写入运行元数据供审计。DDA 或 `Unknown` 模式下本 crate 不动任何谱图、提取计数为 0，保证误判采集模式时也不会破坏数据。

## 6. 约定与不变量

- **单位**：m/z 与质量移量用 `f64`（Da）；电荷 `XicData.charge` 用 `i32`、碎片电荷用 `u32`；scan 编号 `u32` 且 1-based。保留时间字段 `retention_time_min` 原样透传自 `Spectrum::retention_time_min`（`core` 标注为分钟），**不做单位换算**（注：`xic/lib.rs` 的文档注释把它写作 seconds 是过时标签，值即源谱所载）。
- **索引**：色谱规划只读 `list_scan_meta` 的内存索引，按 `scan_number` seek；同位素间距固定 `C13_C12_MASS_DIFF = 1.003355`；SILAC 标准移量 K=8.014199、R=10.008269。
- **确定性**：所有排序显式 `sort_by_key(scan_number)` 或 `partial_cmp(...).unwrap_or(Equal)`，不依赖 `HashMap` 迭代序（`charge_distribution` 仅用于统计计数，不参与排序）；相同输入产出逐位一致。
- **错误处理**：`XicError`（`UnsupportedFormat`/`ScanNotFound`/`NoIsolationWindow`/`NoMatchingCycles`/`InvalidPeptide`/`SpectrumIo`）与 `DiaExtractionError`（`NoMs1Spectra`/`NoMs2Spectra`/`ScanNotFound`/`NoIsolationWindow`/…）皆 `thiserror` 定义；库代码数据路径不 `unwrap/expect`，用 `Result`/`Option`；重标 MS2 缺失走 `heavy_warning` 软告警而非报错。

---

继续阅读：[README](README.md) 总目录 -- 上溯 [L2 架构](L2-architecture.md)。
