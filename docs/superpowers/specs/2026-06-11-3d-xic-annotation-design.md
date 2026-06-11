# 3D MS2 标注可视化设计文档

> **日期**: 2026-06-11
> **状态**: 设计确认（经可视化伴侣 4 轮迭代确认呈现形态，终端确认实现范围）
> **方案**: 在现有 `extract_xic` 工具上新增 `view=3d` 参数，复用 XIC 提取与单谱标注逻辑，生成新的 3D + 逐谱图标注 HTML

---

## 1. 需求概述

### 1.1 核心目标

现有 XIC 可视化是 **2D 折线图**：X=保留时间、Y=强度，每条曲线是一个选定的碎片离子。它只回答「匹配了哪几根碎片离子、它们怎么洗脱」，**看不到整张 MS2 谱图的全貌**。

本特性提供一个 **3D 视角**，让用户不仅看到匹配的碎片离子数量，还能看到 **整张 MS2 谱图随保留时间的演变**，并对窗口内每一张谱图做 b/y 离子标注：

1. **3D 总览**：把 RT 窗口内每一张 MS2 谱图画成竖直柱子（stick），沿保留时间排成 3D 瀑布。三轴 = **保留时间 RT × m/z × 强度**。灰柱 = 普通峰，彩色粗柱 = 匹配到的目标肽段 b/y 碎片离子。
2. **逐谱图 b/y 标注明细**：在 3D 下方，依次列出窗口内**每一张** MS2 谱图，每张给出：
   - scan 号、保留时间（鉴定到肽段的谱图高亮）
   - **总谱峰数**（该谱图全部峰数目）
   - 匹配情况（匹配到几根 / 共几根目标离子）+ 匹配离子彩色标签
   - 经典的 b/y 离子标注棒状图（灰柱 = 普通峰，彩色柱 + 文字 = 匹配的 b/y 离子）

### 1.2 用户场景

- **鉴定可信度判断**：在整张 MS2 谱图背景下观察匹配碎片是否显著、是否被干扰峰淹没。
- **跨谱图共洗脱核查**（DIA）：窗口内多张 MS2 谱图逐张标注，直观看到目标碎片离子随 RT 的出现/消失。
- **谱图质量评估**：每张谱图的总峰数 + 匹配根数一目了然。

### 1.3 视觉形态确认

经可视化伴侣 4 轮迭代（concept → v2 柱状 + 峰数图示 → v3 逐谱图明细表 → v4 逐谱图 b/y 标注），最终确认形态为 **v4**：上方 3D 柱状总览 + 下方逐谱图 b/y 标注棒状图列表（每张含总谱峰数与匹配情况）。

---

## 2. 设计决策汇总

| 决策项 | 选择 | 理由 |
|--------|------|------|
| 交付形态 | 扩展现有 `extract_xic`，新增 `view: "standard" \| "3d"` 参数 | 用户选择；复用现有定位/提取参数与逻辑，零新增工具表面 |
| 默认行为 | `view` 默认 `standard`，保持现有 2D XIC 不变 | 向后兼容 |
| 3D HTML 输出 | `view="3d"` 生成独立自包含 HTML（默认 `output/xic3d_scan{N}.html`） | 与现有 `xic_scan{N}.html` 并列，互不影响 |
| 数据源 | 复用 `extract_xic_unified` 的 `raw_scans.ms2_scans`（窗口内每张完整峰列表） | 数据已现成，无需新增提取通路 |
| 逐谱图标注 | 复用 `search-engine::annotate::annotate_spectrum` 逐张匹配 | 避免重复实现 b/y 匹配；确定性逻辑留在 Rust |
| 谱图范围 | RT 窗口内（±n_cycles，默认 ±5）的**全部** MS2 谱图 | 用户选择「完整呈现」 |
| 3D 总览降密 | 可选 `max_peaks_per_scan_3d`（默认 top-200 按强度），仅作用于 3D 总览 | 控制 3D 点密度/文件大小；逐谱图标注图不受影响 |
| 逐谱图标注图 | 始终包含该谱图全部峰，且必含匹配峰 | 忠实呈现总谱峰，匹配峰不被裁剪 |
| 计算位置 | b/y 匹配、ppm、计数全部在 Rust；JS 仅渲染 | 项目铁律 2.1：确定性逻辑全在 Rust |
| SILAC heavy | v1 不做（仅 light/目标肽跨谱图标注） | YAGNI；可后续扩展 |
| 文件格式 | 仅 mzML | 与现有 XIC 一致（需 MS1+MS2 + 隔离窗口） |
| JS 可视化 | Plotly.js（CDN / Embedded 可配置） | 沿用现有 XIC/unified 技术栈 |

---

## 3. 架构设计

### 3.1 系统分层与数据流

```text
┌──────────────────────────────────────────────────────────────┐
│  MCP Tool Layer (mcp-server)                                   │
│  extract_xic(view = "standard" | "3d")                         │
│    view="standard" → 现有路径：render_xic_html (2D)            │
│    view="3d"       → 新路径（见下）                            │
└──────────────────────────────────────────────────────────────┘
                              │ view="3d"
                              ▼
┌──────────────────────────────────────────────────────────────┐
│  xic crate (已有)                                              │
│  extract_xic_unified(...) → XicUnifiedResult {                 │
│      xic_data, raw_scans: { ms2_scans: Vec<RawScan> }, ... }   │
│  RawScan { scan_number, retention_time_min, mz_array,          │
│            intensity_array }   ← 窗口内每张 MS2 完整峰列表      │
└──────────────────────────────────────────────────────────────┘
                              │ raw_scans.ms2_scans
                              ▼
┌──────────────────────────────────────────────────────────────┐
│  report crate (新增 xic3d 模块)                                │
│  build_xic3d_data(ms2_scans, peptide, charge, precursor_mz,    │
│                   mods, tol, target_scan)                      │
│    for each RawScan:                                           │
│      Spectrum { peaks + 合成目标前体 }                         │
│        → search_engine::annotate::annotate_spectrum(...)       │
│        → SpectrumAnnotation (峰带注释, b/y, matched/total)     │
│    ⇒ Xic3dData { scans: Vec<Ms2ScanAnnotation> }               │
│                              │                                 │
│  render_xic_3d(&Xic3dData, &Path, PlotlyMode)                  │
│    → 注入 JSON 到 templates/xic3d.html (Plotly)                │
└──────────────────────────────────────────────────────────────┘
                              │
                              ▼
                  output/xic3d_scan{N}.html  (自包含)
```

### 3.2 模块归属

- 新增数据结构与 builder、渲染逻辑均放在 **report crate**（已依赖 `search-engine` 与 `xic`，见 `report/src/unified_visualize.rs`、`lib.rs`）。
- mcp-server 仅做参数解析、调用 `extract_xic_unified` + `report` 的 builder/render，保持工具薄。

---

## 4. 数据结构（report crate，全部 `Serialize + Deserialize`）

```rust
// report/src/xic3d_types.rs
use protein_copilot_core::search_params::MassTolerance;
use protein_copilot_search_engine::annotate::SpectrumAnnotation;
use serde::{Deserialize, Serialize};

/// 完整的 3D MS2 标注视图数据。
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Xic3dData {
    pub peptide_sequence: String,
    pub charge: i32,
    pub precursor_mz: f64,
    pub target_scan: u32,
    pub source_file: String,
    pub mz_tolerance: MassTolerance,
    /// 窗口内每张 MS2 谱图的标注（按 scan/RT 升序）。
    pub scans: Vec<Ms2ScanAnnotation>,
}

/// 单张 MS2 谱图在 3D 视图中的标注条目。
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ms2ScanAnnotation {
    pub scan_number: u32,
    pub retention_time_min: f64,
    /// 该谱图的总峰数（= 完整 mz_array 长度）。
    pub total_peaks: usize,
    /// 是否为鉴定到目标肽段的那张谱图。
    pub is_target: bool,
    /// 复用单谱标注结果：含全部峰（可选 b/y 注释）、b_ions/y_ions、matched_ions/total_ions。
    pub annotation: SpectrumAnnotation,
}
```

> 复用 `SpectrumAnnotation`（`search-engine::annotate`）的现有字段：`peaks: Vec<AnnotatedPeak>`、`b_ions`/`y_ions: Vec<TheoreticalIon>`、`matched_ions`、`total_ions`、`scan_number`、`retention_time_min` 等，避免重复定义。

---

## 5. 确定性计算（Rust）

### 5.1 builder：`build_xic3d_data`

```rust
// report/src/xic3d_types.rs (或 xic3d_build.rs)
pub fn build_xic3d_data(
    ms2_scans: &[protein_copilot_xic::RawScan],
    peptide_sequence: &str,
    charge: i32,
    precursor_mz: f64,
    modifications: &[Modification],
    fragment_tolerance: &MassTolerance,
    target_scan: u32,
    source_file: &str,
) -> Result<Xic3dData, ReportError>;
```

实现要点：
1. 对每张 `RawScan` 构造 `Spectrum`：
   - `scan_number / retention_time_min / mz_array / intensity_array` 来自 RawScan；
   - `ms_level = MsLevel::Ms2`；
   - `precursors = vec![PrecursorInfo { mz: precursor_mz, charge: Some(charge), intensity: None, isolation_window: None, source_scan: None }]`（合成目标前体，供 annotate 计算前体 ppm）。
2. 调 `annotate_spectrum(&spectrum, peptide_sequence, charge, fragment_tolerance, modifications, vec![], false, false)` 得 `SpectrumAnnotation`。
3. `total_peaks = spectrum.mz_array.len()`；`is_target = (scan_number == target_scan)`。
4. 收集为 `Vec<Ms2ScanAnnotation>`，按 scan 升序。

> builder 产出的 `scans[].annotation.peaks` 始终是**完整峰列表**（供逐谱图标注图全量显示，也作为 3D 总览的数据源）。3D 总览的 top-N 降密属于**纯显示优化**，在渲染层处理（见 §5.2），不改动任何统计数字。
>
> `ms2_scans` 为空时，返回明确的 `ReportError`（窗口内无 MS2 谱图）；单张谱图无峰或 annotate 出错（如非标准残基），返回带上下文的 `ReportError`，由 mcp-server 转成清晰的 MCP 错误。

### 5.2 渲染：`render_xic_3d`

```rust
// report/src/xic3d_visualize.rs
const TEMPLATE: &str = include_str!("../templates/xic3d.html");

pub fn render_xic_3d(
    data: &Xic3dData,
    output_path: &Path,
    plotly_mode: PlotlyMode,
    max_peaks_per_scan_3d: Option<usize>, // 仅 3D 总览的 display-only 降密；缺省 = 200
) -> Result<(), ReportError>;
```

- 序列化 `Xic3dData` → `escape_json_for_html`（复用 `report::escape_json_for_html`，防 `</script>` 注入）→ 注入模板占位符 `window.__XIC3D_DATA__ = {...};`。
- `max_peaks_per_scan_3d` 注入为模板 JS 常量（display-only 降密，匹配峰恒保留），不改动数据本身。
- Plotly 源按 `PlotlyMode`（Cdn / Embedded），与 `unified_visualize.rs` 一致。
- 自动创建父目录；写出自包含 HTML。
- `ReportGenerator::render_xic_3d(...)` 在 `report/src/lib.rs` 暴露。

### 5.3 模板：`templates/xic3d.html`

- 复用现有 xic/unified 模板的深色主题与布局。
- JS（消费 `window.__XIC3D_DATA__`）：
  - **3D 总览**：遍历 `scans`，每张谱图的每个峰画竖直 stick（`scatter3d` mode lines，`(rt,mz,0)→(rt,mz,intensity)`，null 分隔）；匹配峰（peak.annotation 非空）按 b/y 上色加粗；普通峰灰色。可选 top-N 降密。
  - **逐谱图列表**：遍历 `scans`，每张渲染标题（scan/RT/总谱峰数/匹配 n/N/匹配离子标签，target 高亮）+ 一个 2D 标注棒状图（灰柱全部峰，彩色柱 + 文字标注匹配 b/y 离子）。

---

## 6. MCP 工具变更

`extract_xic` 输入新增字段：

```rust
/// 视图模式：standard（默认，现有 2D XIC）或 3d（新 3D + 逐谱图标注）。
#[serde(default)]
view: Option<XicView>,            // enum { Standard, ThreeD }，缺省 = Standard

/// 仅 view=3d：3D 总览每张谱图最多绘制的峰数（按强度 top-N）。默认 200。
#[serde(default)]
max_peaks_per_scan_3d: Option<usize>,
```

流程分支（`extract_xic` 现有逻辑之后）：
- `Standard` → 现有 `render_xic_html(&xic_data, ...)`，输出 `xic_scan{N}.html`（不变）。
- `ThreeD` → `build_xic3d_data(unified_result.raw_scans.ms2_scans, ...)` + `render_xic_3d(...)`，输出默认 `xic3d_scan{N}.html`。

返回结果 `ExtractXicResult` 在 3D 模式复用 `output_path`，并将 `ms2_scan_count` 设为窗口内 MS2 谱图张数；新增 `annotated_scan_count`（成功标注的谱图数）与 `target_matched_ions`（目标 scan 匹配到的碎片离子数）。`summary` 文案相应描述 3D 视图。

---

## 7. 边界与错误处理

- **窗口内无 MS2 谱图**：返回明确错误（提示调大 n_cycles 或检查 scan/RT）。
- **非标准残基 / 无法生成理论离子**：annotate 报错，转成面向用户的 MCP 错误（说明哪个序列、如何修复）。
- **DDA 数据**：窗口内 MS2 谱图可能较少，逐张列出同样适用（复用现有隔离窗口匹配逻辑）。
- **大谱图**：3D 总览用 `max_peaks_per_scan_3d` 降密；逐谱图标注图保留全部峰（数量大时浏览器仍可承受，~十几张 × 数百峰）。
- **文件大小**：预计几十 KB 量级，可接受。

---

## 8. 测试计划

### 8.1 单元测试（report crate）

- `build_xic3d_data`：
  - `total_peaks` 等于输入峰数；
  - `matched_ions`/`total_ions` 计数正确（构造已知匹配的合成谱）；
  - `is_target` 正确标记目标 scan；
  - 空窗口（无 ms2_scans）→ 返回 `ReportError`（明确「窗口内无 MS2」）；
  - 非标准残基 → 返回错误。
- `render_xic_3d`：
  - 产物 HTML 含注入数据 `__XIC3D_DATA__`；
  - 含 Plotly 源（CDN 字符串 / Embedded）；
  - JSON 经 `escape_json_for_html` 转义（无裸 `</script>`）；
  - `max_peaks_per_scan_3d` 注入为模板常量（断言其值出现在 HTML 中）；
  - 自动创建父目录。

### 8.2 集成测试（tests/ 或 mcp-server）

- `extract_xic` with `view="3d"` 跑 `tests/fixtures` 中的 mzML：
  - 成功生成 HTML 文件；
  - scan 张数与窗口（±n_cycles）一致；
  - 命中并标记 target scan；
  - `view` 缺省时行为与现状一致（回归保护）。

---

## 9. 范围边界（v1）

**做**：
- `extract_xic view=3d` 生成 3D 总览 + 逐谱图 b/y 标注 HTML（light / 目标肽）。
- 复用现有 XIC 提取与单谱标注。

**不做（YAGNI，可后续）**：
- SILAC heavy 跨谱图 3D 标注。
- 表格行 ↔ 3D 联动高亮、镜像（理论 vs 实测）图等交互增强。
- mzML 之外的格式。

---

## 10. 涉及文件

| 文件 | 变更 |
|------|------|
| `crates/report/src/xic3d_types.rs` | 新增：`Xic3dData`、`Ms2ScanAnnotation`、`build_xic3d_data` |
| `crates/report/src/xic3d_visualize.rs` | 新增：`render_xic_3d` |
| `crates/report/templates/xic3d.html` | 新增：3D + 逐谱图标注 Plotly 模板 |
| `crates/report/src/lib.rs` | 新增模块声明 + `ReportGenerator::render_xic_3d` |
| `crates/mcp-server/src/tools.rs` | `ExtractXicInput` 增 `view` / `max_peaks_per_scan_3d`；`extract_xic` 分支 3D 路径 |
| `tests/` | 集成测试：`view=3d` |
