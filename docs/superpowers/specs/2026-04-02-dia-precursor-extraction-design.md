# DIA Precursor Extraction — Design Spec

## Problem

当前 ProteinCopilot 对 DIA 数据的支持仅限于解析 mzML 中的 `IsolationWindow`，但搜索引擎只使用 `precursors.first()` 的单一 precursor_mz，忽略了 DIA 宽隔离窗口内包含多个共碎裂母离子的事实。这导致 DIA 数据的搜索结果不准确。

## Approach

创建独立的 `dia-extraction` crate，从 MS1 谱图中通过同位素模式检测提取候选母离子，填充到 DIA MS2 谱图的 `precursors` vec 中，作为搜索引擎的前置预处理步骤。同时暴露为 MCP tool `extract_dia_precursors` 供用户调用。

## Architecture

### Data Flow

```
mzML/mgf 文件
    ↓
spectrum-io（读取 MS1+MS2，解析 spectrumRef）
    ↓ Vec<Spectrum>
dia-extraction（NEW）
    ├─ 自动检测 DDA/DIA（隔离窗口宽度中位数 > 5 Da → DIA）
    ├─ 分离 MS1 / MS2 谱图
    ├─ 关联 MS1 ↔ MS2（spectrumRef → 扫描顺序 → RT）
    ├─ MS1 同位素模式检测 → 候选 precursor
    ├─ 按隔离窗口范围过滤
    └─ 输出增强后的谱图
    ↓
    ├─ 模式 A: 多 precursor 格式（内部引擎优化）
    └─ 模式 B: pseudo-spectra 拆分（外部引擎通用）
    ↓
search-engine（遍历所有 precursors 匹配 → PSM 结果，FDR 控制）
```

### New Crate: `crates/dia-extraction`

#### Core Trait

```rust
pub trait PrecursorExtractor: Send + Sync {
    fn extract(
        &self,
        ms1: &Spectrum,
        isolation_window: &IsolationWindow,
    ) -> Vec<PrecursorInfo>;
}
```

`PrecursorExtractor` trait 抽象算法，内置 `IsotopePatternExtractor` 实现，预留外部算法扩展。

#### Built-in Implementation: IsotopePatternExtractor

配置参数：
- `min_charge: i32` — 最小电荷态（默认 2）
- `max_charge: i32` — 最大电荷态（默认 5）
- `isotope_tolerance_da: f64` — 同位素间距容差（默认 0.01 Da）
- `min_isotope_peaks: usize` — 最少同位素峰数（默认 2）
- `min_intensity_ratio: f64` — 最小强度比（默认 0.1）

算法步骤：

1. **筛选 MS1 峰**：只保留 m/z 在隔离窗口 `[target - lower, target + upper]` 范围内的峰，按强度降序排序。
2. **构建同位素簇**：对每个峰，尝试 `z = min_charge..=max_charge`，计算期望间距 `Δ = 1.00335 / z`（¹³C-¹²C 质量差 / 电荷），向右寻找匹配峰（容差 `isotope_tolerance_da`），检查强度递减趋势。
3. **构建 PrecursorInfo**：mz = 单同位素峰，charge = 推断电荷态，intensity = 单同位素峰强度。
4. **去重与排序**：合并相近候选（Δm/z < 0.01 Da），按强度降序排序，标记已使用峰避免重复。

#### Main Entry Function

```rust
pub fn extract_dia_precursors(
    spectra: &[Spectrum],
    extractor: &dyn PrecursorExtractor,
    config: &DiaExtractionConfig,
) -> Result<DiaExtractionResult, DiaExtractionError>
```

#### DiaExtractionConfig

```rust
pub struct DiaExtractionConfig {
    pub acquisition_mode: Option<AcquisitionMode>,  // None = 自动检测
    pub dia_threshold_da: f64,                       // DIA 判断阈值（默认 5.0 Da）
}
```

#### DiaExtractionResult

```rust
pub struct DiaExtractionResult {
    pub detected_mode: AcquisitionMode,
    pub enhanced_spectra: Vec<Spectrum>,  // MS2 谱图（precursors 已填充）
    pub stats: ExtractionStats,
}

pub struct ExtractionStats {
    pub ms1_count: u32,
    pub ms2_count: u32,
    pub total_precursors_extracted: u32,
    pub avg_precursors_per_ms2: f64,
    pub charge_distribution: HashMap<i32, u32>,
}

impl DiaExtractionResult {
    /// 拆分为 pseudo-spectra（每个 spectrum 只有 1 个 precursor）
    pub fn expand_to_pseudo_spectra(&self) -> Vec<Spectrum>;
}
```

### MS1 ↔ MS2 关联策略（三级回退）

1. **spectrumRef**（最可靠）：mzML 的 `<precursor spectrumRef="...">` 属性直接指向对应 MS1 scan。需要在 spectrum-io mzML reader 中解析此属性。
2. **扫描顺序**（回退）：对每个 MS2，向前找 scan_number 最近的 MS1 谱图。假设标准 DIA 采集模式（MS1 → MS2×N → MS1 → ...）。
3. **保留时间**（最后回退）：找 `retention_time_sec` 最接近的 MS1 谱图。适用于不规则采集或 scan_number 不连续的情况。

### DDA/DIA 自动检测

分析所有 MS2 谱图的隔离窗口宽度（`lower_offset + upper_offset`），计算中位数。中位数 > `dia_threshold_da`（默认 5.0 Da）→ DIA，否则 → DDA。用户可通过 `DiaExtractionConfig.acquisition_mode` 覆盖。

## Changes to Existing Crates

### crates/core/src/spectrum.rs

- `PrecursorInfo` 新增字段：`pub source_scan: Option<u32>` — spectrumRef 对应的 MS1 scan number
- 新增枚举：`pub enum AcquisitionMode { DDA, DIA, Unknown }`

### crates/core/src/search_params.rs

- `SearchParams` 新增字段：`pub acquisition_mode: Option<AcquisitionMode>` — 用户可指定采集模式

### crates/spectrum-io/src/mzml.rs

- 解析 `<precursor spectrumRef="...">` 属性，提取 scan number 存入 `PrecursorInfo.source_scan`
- spectrumRef 格式通常为 `"controllerType=0 controllerNumber=1 scan=1234"`，需要解析出 scan number

### crates/search-engine/src/matching.rs

- `match_spectrum()` 改为遍历所有 `spectrum.precursors`（而非只取 `first()`）
- 返回值改为 `Vec<PeptideMatch>`，收集所有 precursor 的匹配结果，由下游 FDR 控制筛选

### crates/search-engine/src/simple_engine.rs

- 主搜索循环增加 MS level 过滤：只搜索 `MsLevel::MS2` 谱图

### crates/mcp-server/src/main.rs

- 新增 MCP Tool：`extract_dia_precursors`
  - 输入：`file_path: String`, `output_mode: "multi" | "pseudo"`, `min_charge: Option<i32>`, `max_charge: Option<i32>`
  - 输出：`acquisition_mode`, `ms1_count`, `ms2_count`, `total_precursors_extracted`, `avg_precursors_per_ms2`
  - 增强后的谱图缓存在内存中，供后续 `run_search` 使用

## MCP Tool: extract_dia_precursors

### Description

从 DIA 质谱数据中提取候选母离子。读取 mzML 文件中的 MS1 和 MS2 谱图，通过同位素模式检测在 MS1 中识别可能的母离子，将候选母离子信息填充到对应 DIA MS2 谱图中。输出增强后的谱图可直接用于 `run_search`。

### Input Schema

```json
{
  "file_path": "string (required) — mzML 文件路径",
  "output_mode": "string (optional, default: 'pseudo') — 'multi'(多precursor) | 'pseudo'(拆分)",
  "min_charge": "integer (optional, default: 2) — 最小电荷态",
  "max_charge": "integer (optional, default: 5) — 最大电荷态",
  "acquisition_mode": "string (optional) — 'DDA' | 'DIA'，覆盖自动检测"
}
```

### Output Schema

```json
{
  "detected_mode": "DDA | DIA",
  "ms1_count": "u32",
  "ms2_count": "u32",
  "total_precursors_extracted": "u32",
  "avg_precursors_per_ms2": "f64",
  "charge_distribution": "HashMap<i32, u32>",
  "run_id": "string — 缓存 ID，供 run_search 引用"
}
```

## Testing

### Unit Tests (dia-extraction crate)

- 同位素簇检测：已知 2+、3+、4+ 模式 → 正确电荷态和 m/z
- 噪声峰不误检：随机峰不构成有效同位素簇
- 隔离窗口过滤：窗口外的同位素簇不被提取
- MS1-MS2 关联：三种策略各自正确性
- DDA/DIA 自动检测：窄窗口 → DDA，宽窗口 → DIA
- 边界情况：空 MS1、空隔离窗口、无匹配峰

### Integration Tests

- 合成 DIA mzML → 完整提取流程 → 验证候选 precursor 正确
- pseudo-spectra 展开 → 验证每个 spectrum 只有 1 个 precursor
- DIA 提取 → 搜索引擎 → 验证端到端流程
- DDA 数据直通 → 不执行提取，保持原有行为
- MCP Tool 调用 → 验证输入输出格式

## Error Handling

- `DiaExtractionError::NoMs1Spectra` — DIA 数据但无 MS1 谱图
- `DiaExtractionError::NoMs2Spectra` — 无 MS2 谱图可处理
- `DiaExtractionError::InvalidIsolationWindow` — 隔离窗口参数异常
- `DiaExtractionError::SpectrumIoError` — 底层谱图读取错误

## Out of Scope

- Averagine 模型拟合（后续扩展，通过 PrecursorExtractor trait）
- 外部 DIA 工具输出导入（如 DIA-NN 结果）
- 色谱维度的 feature detection（RT 方向的峰检测）
- DIA-specific FDR 计算
- 谱库搜索（spectral library search）
