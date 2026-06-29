# L4 — dia-extraction（DIA 母离子提取）

承接 [L3-xic-dia](L3-xic-dia.md) 与 [L2](L2-architecture.md)。本篇只讲 `crates/dia-extraction` 一个 crate：对外 API 与代码骨架。所有签名、常量、默认值均以源码 `src/{lib,config,error,detection,correlation,extractor,isotope}.rs` 为准，不臆造。

## 1. 用途 + 位置 + 依赖

DIA 谱图的宽隔离窗里混着多个共碎裂前体，搜索引擎无法用「一谱一前体」直接搜。本 crate 从 MS1 同位素模式反推候选前体（m/z + 电荷），写回每张 MS2 的 `precursors`。

- 输入：一批已在内存的 `Spectrum`（MS1 + MS2），本 crate 不读盘、不联网、无 LLM。
- 输出：`DiaExtractionResult`（增强后的 MS2 + 统计），由 mcp-server 缓存为 `run_id`，供 `run_search(dia_run_id=...)` 取用。
- 依赖：仅 `protein-copilot-core`（`Spectrum`/`PrecursorInfo`/`IsolationWindow`/`AcquisitionMode`/`MsLevel`）+ `thiserror` + `serde` + `schemars` + `tracing`。无 tokio——纯同步内存计算。

提取算法以 trait object（`&dyn PrecursorExtractor`）注入，编排逻辑与具体打分解耦：换一种同位素判定只需另写一个实现，`extract_dia_precursors` 的流程不动。读盘、缓存、生成 `run_id` 都在 mcp-server 层完成，本 crate 收到的永远是内存里的 `&[Spectrum]`，因而单元测试无需任何 mzML 文件。

## 2. 对外 API（lib.rs 重导出）

两个自由函数 + 一个 trait + 一个默认实现：

```rust
pub fn extract_dia_precursors(
    spectra: &[Spectrum],
    extractor: &dyn PrecursorExtractor,
    config: &DiaExtractionConfig,
) -> Result<DiaExtractionResult, DiaExtractionError>;

pub fn extract_single_spectrum_precursors(
    spectra: &[Spectrum],
    scan_number: u32,
    extractor: &dyn PrecursorExtractor,
) -> Result<SingleSpectrumExtractionResult, DiaExtractionError>;

pub trait PrecursorExtractor: Send + Sync {
    fn extract(&self, ms1: &Spectrum, isolation_window: &IsolationWindow) -> Vec<PrecursorInfo>;
}
```

detection / correlation 子模块导出的辅助函数：

```rust
pub fn detect_acquisition_mode(spectra: &[Spectrum], threshold_da: f64) -> AcquisitionMode;
pub fn separate_by_ms_level(spectra: &[Spectrum]) -> (Vec<&Spectrum>, Vec<&Spectrum>);
pub fn correlate_ms1_ms2(ms1: &[&Spectrum], ms2: &[&Spectrum]) -> Vec<Option<usize>>;
pub fn correlate_single_with_method(ms1: &[&Spectrum], ms2: &Spectrum) -> (Option<usize>, &'static str);
```

配置与结果字段（config.rs，逐字段）：

```rust
pub struct DiaExtractionConfig {
    pub acquisition_mode: Option<AcquisitionMode>,  // None = 自动检测
    pub dia_threshold_da: f64,                      // 默认 5.0
}
pub struct DiaExtractionResult {
    pub detected_mode: AcquisitionMode,
    pub enhanced_spectra: Vec<Spectrum>,
    pub stats: ExtractionStats,
}
pub struct ExtractionStats {
    pub ms1_count: u32,
    pub ms2_count: u32,
    pub total_precursors_extracted: u32,
    pub avg_precursors_per_ms2: f64,
    pub charge_distribution: HashMap<i32, u32>,
}
pub struct SingleSpectrumExtractionResult {
    pub ms2_scan: u32,
    pub ms1_scan_used: u32,
    pub correlation_method: String,                 // source_scan | scan_order | rt_nearest
    pub isolation_window: Option<IsolationWindow>,
    pub precursors: Vec<PrecursorInfo>,
}
```

`DiaExtractionResult` 有三个取数方法：`into_enhanced_spectra(self)`（按窗多前体一谱）、`expand_to_pseudo_spectra(&self)` 与 `into_pseudo_spectra(self)`（拆成一谱一前体的伪谱，兼容「一谱一前体」引擎）。错误枚举 `DiaExtractionError`：`NoMs1Spectra`、`NoMs2Spectra`、`InvalidIsolationWindow{detail}`、`ExtractionFailed{detail}`、`ScanNotFound{scan}`、`NoIsolationWindow{scan}`。

## 3. 关键常量与默认

- `IsotopePatternExtractor`（isotope.rs，default 实现）：`min_charge=2`、`max_charge=5`（电荷范围 `2..=5`）、`isotope_tolerance_da=0.01`、`min_isotope_peaks=2`、`min_intensity_ratio=0.1`。
- 同位素间距常量：`const C13_C12_MASS_DIFF: f64 = 1.003_355;`，第 z 电荷态相邻峰间距 `= C13_C12_MASS_DIFF / z`。
- DIA 判定阈值：`dia_threshold_da` 默认 `5.0`，比较的是 MS2 窗宽 `lower_offset + upper_offset` 的中位数。
- 关联三级回退顺序：`source_scan -> scan_order -> rt_nearest`（都不命中则 `none`）。

几个默认值的含义：向前串峰时若某峰强度 > 前一峰的 `1.5` 倍即停（同位素包络应随质量增大而递减，强度反弹说明撞上了另一簇）；任一峰强度低于 `min_intensity_ratio * seed_int` 也停（滤掉噪声尾巴）。一个种子峰可能本身是 M+1（大肽里 M+1 常高于 M+0），故算法先向后探一格找单同位素峰 M+0，最终上报的 `mz` 取簇内最小 m/z 而非种子峰。多电荷态都成簇时选峰数最多者，峰数并列再取更低电荷——低电荷态在 DIA 中更常见，作为先验更稳。

## 4. 源码骨架（简化）

采集模式：取所有 MS2 窗宽 `lower+upper` 的中位数与阈值比较（detection.rs）：

```rust
pub fn detect_acquisition_mode(spectra: &[Spectrum], threshold_da: f64) -> AcquisitionMode {
    let widths: Vec<f64> = /* MS2 每个前体 isolation_window 的 lower_offset + upper_offset */;
    if widths.is_empty() { return AcquisitionMode::Unknown; }
    if median_f64(&widths) > threshold_da { AcquisitionMode::DIA } else { AcquisitionMode::DDA }
}
```

MS2 关联 MS1：三级回退（correlation.rs）：

```rust
pub fn correlate_single_with_method(ms1: &[&Spectrum], ms2: &Spectrum) -> (Option<usize>, &'static str) {
    // L1 source_scan：前体元数据里显式引用的 MS1 scan
    if let Some(src) = ms2.precursors.first().and_then(|p| p.source_scan) {
        if let Some(i) = ms1.iter().position(|m| m.scan_number == src) { return (Some(i), "source_scan"); }
    }
    // L2 scan_order：scan_number 最大且 < MS2 的 MS1
    if let Some((i, _)) = ms1.iter().enumerate()
        .filter(|(_, m)| m.scan_number < ms2.scan_number)
        .max_by_key(|(_, m)| m.scan_number) { return (Some(i), "scan_order"); }
    // L3 rt_nearest：保留时间最近（跳过 NaN）
    let by_rt = /* min_by |RT(MS1) - RT(MS2)| */;
    match by_rt { Some(i) => (Some(i), "rt_nearest"), None => (None, "none") }
}
```

同位素模式提取骨架（isotope.rs，`PrecursorExtractor::extract`）：

```rust
fn extract(&self, ms1: &Spectrum, iw: &IsolationWindow) -> Vec<PrecursorInfo> {
    if self.min_charge < 1 || self.max_charge < self.min_charge { return Vec::new(); }
    let (low, high) = (iw.target_mz - iw.lower_offset, iw.target_mz + iw.upper_offset);
    let mut peaks = /* 窗内峰 (idx, mz, intensity) */;
    peaks.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(Equal));   // 强度降序
    let mut used = vec![false; ms1.mz_array.len()];
    for &(seed, seed_mz, seed_int) in &peaks {
        if used[seed] { continue; }
        for z in self.min_charge..=self.max_charge {   // 2..=5
            let delta = C13_C12_MASS_DIFF / z as f64;   // 1.003355 / z
            // 向后找 M+0 (seed_mz - delta)；向前串 M+1, M+2 ...
            //   容差 isotope_tolerance_da(0.01)；峰强不得 > 前峰 1.5 倍；
            //   且 >= min_intensity_ratio(0.1) * seed_int
        }
        // 选峰数最多的簇；峰数并列取更低电荷；mono m/z = 簇内最小 mz；charge = z
    }
    // 去重(同电荷 0.01 Da 内留强者) -> 按强度降序
}
```

`extract_dia_precursors` 的编排：`separate_by_ms_level` 分组；MS2 为空报 `NoMs2Spectra`；模式取 `config.acquisition_mode` 或自动检测；DDA/Unknown 时原样回传（0 提取，旁路）；DIA 且 MS1 为空报 `NoMs1Spectra`；否则 `correlate_ms1_ms2` 后逐 MS2 在其隔离窗内 `extractor.extract`。无隔离窗的 MS2 会被计数并 `tracing::warn!`，不报错。

## 5. 调用链（mcp-server）

```text
extract_dia_precursors (MCP tool, tools.rs)
  | validate_file_path -> get_or_create_reader -> read_all(path)   // 需全部谱图
  | IsotopePatternExtractor{min_charge,max_charge} + DiaExtractionConfig{acquisition_mode}
  | run_dia_extraction(spectra, extractor, config) -> DiaExtractionResult
  | output_mode=="multi" -> into_enhanced_spectra() | 否则 -> into_pseudo_spectra()
  | dia_cache.insert(run_id = UUID, spectra) -> 返回 DiaExtractionOutput{ run_id, stats }
  +-> run_search(dia_run_id = run_id) 从 dia_cache 取增强谱再搜（见 L2 第 3 节旁路）

extract_spectrum_precursors (MCP tool, tools.rs)
  | get_or_create_reader -> read_spectrum(path, scan) O(1) seek   // 须为 MS2
  | list_scan_meta -> 选 RT +/- 1.0 min 内的 MS1 -> 逐张 read_spectrum
  | extract_single_spectrum_precursors(spectra, scan, extractor)
  +-> SingleSpectrumExtractionResult{ ms1_scan_used, correlation_method, precursors }
```

要点：批量提取走 `read_all`（要全谱建 MS1/MS2 关联）；单谱探查走 `read_spectrum` 的 O(1) 定点 seek + `list_scan_meta` 内存索引，只读目标 MS2 与邻近 MS1，符合 spectrum-io 读取规范。两个工具都用 `IsotopePatternExtractor::default()` 再按入参覆盖电荷范围（`min_charge>=1`、`max_charge>=min_charge`，否则 `INVALID_PARAMS`）。

## 6. 测试入口

```bash
cargo test -p protein-copilot-dia-extraction --offline
```

单元测试随模块：detection(6) + correlation(5) + isotope(8) + lib(10)，共 29 个，覆盖 DDA/DIA 判定、三级回退、z=2/z=3 同位素簇、DDA 旁路、伪谱展开与四类错误路径。

返回 [README](README.md)。
