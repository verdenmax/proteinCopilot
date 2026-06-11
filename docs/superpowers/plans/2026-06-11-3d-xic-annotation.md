# 3D MS2 标注视图 实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在现有 `extract_xic` MCP 工具上新增 `view=3d` 模式，生成「3D 柱状 MS2 总览 + 逐谱图 b/y 离子标注（含总谱峰数）」的自包含 HTML。

**Architecture:** 复用 `extract_xic_unified` 已提取的 `raw_scans.ms2_scans`（目标隔离窗口内 ±n_cycles 的完整 MS2 峰列表），对每张谱图调用 `annotate_spectrum` 做确定性 b/y 匹配，组装为 `Xic3dData`，再由 report crate 的新模板 `xic3d.html` 用 Plotly 渲染。所有匹配/计数在 Rust 完成，JS 仅渲染。

**Tech Stack:** Rust（report / search-engine / xic / mcp-server crates）、serde、Plotly.js（CDN，自包含 HTML）、thiserror、tracing；测试用 `cargo test` + tempfile。

> 设计文档：`docs/superpowers/specs/2026-06-11-3d-xic-annotation-design.md`

---

## 文件结构（创建/修改）

| 文件 | 职责 | 变更 |
|------|------|------|
| `crates/report/src/xic3d_types.rs` | `Xic3dData` / `Ms2ScanAnnotation` 数据结构 | 创建 |
| `crates/report/src/error.rs` | 新增 `EmptyMs2Window` / `AnnotationError` 错误变体 | 修改 |
| `crates/report/src/xic3d_build.rs` | `build_xic3d_data`：逐谱图 b/y 标注（确定性） | 创建 |
| `crates/report/templates/xic3d.html` | 3D + 逐谱图标注 Plotly 模板 | 创建 |
| `crates/report/src/xic3d_visualize.rs` | `render_xic_3d`：注入 JSON 渲染 HTML | 创建 |
| `crates/report/src/lib.rs` | 模块声明 + `ReportGenerator::render_xic_3d` | 修改 |
| `crates/mcp-server/src/tools.rs` | `ExtractXicInput` 增 `view`/`max_peaks_per_scan_3d`；`extract_xic` 分支 3D 路径；`ExtractXicResult` 增字段 | 修改 |

## 任务总览

- **Task 1**：report crate — `Xic3dData` / `Ms2ScanAnnotation` 类型 + `ReportError` 新变体（serde round-trip 测试）
- **Task 2**：report crate — `build_xic3d_data` 构建器（核心确定性逻辑 + 角落测试）
- **Task 3**：report crate — `templates/xic3d.html` 模板（创建）
- **Task 4**：report crate — `render_xic_3d` + `ReportGenerator::render_xic_3d` + lib 接线（渲染测试）
- **Task 5**：mcp-server — `extract_xic` 接入 `view=3d`（参数/分支/返回字段 + 构建与 clippy 验证 + 冒烟）

---

### Task 1: report crate — `Xic3dData` / `Ms2ScanAnnotation` 类型

**Files:**
- Create: `crates/report/src/xic3d_types.rs`
- Modify: `crates/report/src/lib.rs`（新增模块声明）
- Test: `crates/report/src/xic3d_types.rs`（同文件 `#[cfg(test)] mod tests`）

- [ ] **Step 1: 声明模块 + 写失败测试**

修改 `crates/report/src/lib.rs`，在 `pub mod unified_visualize;` 之后新增一行：

```rust
pub mod xic3d_types;
```

创建 `crates/report/src/xic3d_types.rs`，先只放测试（类型尚未定义，故意编译失败）：

```rust
//! Data structures for the 3D MS2 annotation view (`extract_xic view=3d`).
//!
//! [`Xic3dData`] bundles, for one identified peptide, the per-scan b/y
//! annotation of every MS2 spectrum in the target isolation window's
//! ±n_cycles RT range. Built by [`crate::xic3d_build::build_xic3d_data`]
//! and rendered by [`crate::xic3d_visualize::render_xic_3d`].

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};

    #[test]
    fn xic3d_data_serde_round_trip() {
        let data = Xic3dData {
            peptide_sequence: "PEPTIDEK".to_string(),
            charge: 2,
            precursor_mz: 450.25,
            target_scan: 1852,
            source_file: "sample.mzML".to_string(),
            mz_tolerance: MassTolerance { value: 20.0, unit: ToleranceUnit::Ppm },
            scans: vec![],
        };
        let json = serde_json::to_string(&data).unwrap();
        let back: Xic3dData = serde_json::from_str(&json).unwrap();
        assert_eq!(back.peptide_sequence, "PEPTIDEK");
        assert_eq!(back.target_scan, 1852);
        assert_eq!(back.charge, 2);
        assert!(back.scans.is_empty());
    }
}
```

- [ ] **Step 2: 跑测试确认编译失败**

Run: `cargo test -p protein-copilot-report xic3d_data_serde_round_trip`
Expected: 编译错误 `cannot find type Xic3dData in this scope`

- [ ] **Step 3: 写最小实现（在文件顶部、测试模块之前插入）**

```rust
use protein_copilot_core::search_params::MassTolerance;
use protein_copilot_search_engine::annotate::SpectrumAnnotation;
use serde::{Deserialize, Serialize};

/// Complete data for the 3D MS2 annotation view.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Xic3dData {
    /// Identified peptide sequence.
    pub peptide_sequence: String,
    /// Precursor charge state.
    pub charge: i32,
    /// Real precursor m/z (PSM-derived, not DIA window center).
    pub precursor_mz: f64,
    /// Scan number that identified the peptide.
    pub target_scan: u32,
    /// Source spectrum file name.
    pub source_file: String,
    /// Fragment m/z tolerance used for b/y matching.
    pub mz_tolerance: MassTolerance,
    /// Per-scan annotations (ascending by scan number).
    pub scans: Vec<Ms2ScanAnnotation>,
}

/// One MS2 spectrum's annotation entry within the 3D view.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ms2ScanAnnotation {
    /// Scan number (1-based).
    pub scan_number: u32,
    /// Retention time in minutes.
    pub retention_time_min: f64,
    /// Total number of peaks in this spectrum (= full `mz_array` length).
    pub total_peaks: usize,
    /// Whether this is the scan that identified the peptide.
    pub is_target: bool,
    /// Reused single-spectrum annotation: peaks (optionally b/y annotated),
    /// `b_ions`/`y_ions`, `matched_ions`/`total_ions`, etc.
    pub annotation: SpectrumAnnotation,
}
```

- [ ] **Step 4: 跑测试确认通过**

Run: `cargo test -p protein-copilot-report xic3d_data_serde_round_trip`
Expected: PASS（1 passed）

- [ ] **Step 5: 提交**

```bash
git add crates/report/src/xic3d_types.rs crates/report/src/lib.rs
git commit -m "feat(report): add Xic3dData/Ms2ScanAnnotation types for 3D MS2 view

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: report crate — `build_xic3d_data` 构建器

**Files:**
- Modify: `crates/report/src/error.rs`（新增 2 个错误变体）
- Create: `crates/report/src/xic3d_build.rs`
- Modify: `crates/report/src/lib.rs`（新增模块声明）
- Test: `crates/report/src/xic3d_build.rs`（同文件 `#[cfg(test)] mod tests`）

- [ ] **Step 1: 写失败测试 + 声明模块**

修改 `crates/report/src/lib.rs`，在 `pub mod xic3d_types;` 之后新增：

```rust
pub mod xic3d_build;
```

创建 `crates/report/src/xic3d_build.rs`，先只放测试（构建器未定义，编译失败）：

```rust
//! Builder for [`Xic3dData`]: per-scan b/y annotation of every MS2 spectrum
//! in the target isolation window's ±n_cycles RT range.

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::search_params::ToleranceUnit;
    use protein_copilot_search_engine::matching::{
        generate_b_ions_with_charge, generate_y_ions_with_charge,
    };
    use protein_copilot_xic::RawScan;

    fn tol() -> MassTolerance {
        MassTolerance { value: 20.0, unit: ToleranceUnit::Ppm }
    }

    /// A RawScan whose peaks include the first 3 b- and 3 y-ions (1+) of `pep`,
    /// plus two non-matching peaks. Peaks sorted ascending by m/z.
    fn raw_scan_with_matches(scan: u32, rt: f64, pep: &str) -> RawScan {
        let b = generate_b_ions_with_charge(pep, &[], 1);
        let y = generate_y_ions_with_charge(pep, &[], 1);
        let mut mz: Vec<f64> = Vec::new();
        mz.extend(b.iter().take(3).copied());
        mz.extend(y.iter().take(3).copied());
        mz.push(150.0);
        mz.push(1234.5);
        mz.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let intensity = vec![1000.0; mz.len()];
        RawScan { scan_number: scan, retention_time_min: rt, mz_array: mz, intensity_array: intensity }
    }

    #[test]
    fn build_annotates_each_scan_and_marks_target() {
        let pep = "PEPTIDEK";
        let scans = vec![
            raw_scan_with_matches(100, 19.9, pep),
            raw_scan_with_matches(101, 20.0, pep),
            raw_scan_with_matches(102, 20.1, pep),
        ];
        let data = build_xic3d_data(&scans, pep, 2, 460.0, &[], &tol(), 101, "s.mzML").unwrap();
        assert_eq!(data.scans.len(), 3);
        assert_eq!(data.scans[0].total_peaks, scans[0].mz_array.len());
        assert_eq!(data.scans.iter().filter(|s| s.is_target).count(), 1);
        assert!(data.scans[1].is_target);
        assert!(
            data.scans[1].annotation.matched_ions >= 6,
            "expected >=6 matched, got {}",
            data.scans[1].annotation.matched_ions
        );
    }

    #[test]
    fn build_marks_matched_peaks() {
        let pep = "PEPTIDEK";
        let scans = vec![raw_scan_with_matches(100, 20.0, pep)];
        let data = build_xic3d_data(&scans, pep, 2, 460.0, &[], &tol(), 100, "s.mzML").unwrap();
        let annotated = data.scans[0]
            .annotation
            .peaks
            .iter()
            .filter(|p| p.annotation.is_some())
            .count();
        assert!(annotated >= 6, "expected >=6 annotated peaks, got {annotated}");
    }

    #[test]
    fn build_empty_window_errors() {
        let err = build_xic3d_data(&[], "PEPTIDEK", 2, 460.0, &[], &tol(), 100, "s.mzML")
            .unwrap_err();
        assert!(matches!(err, ReportError::EmptyMs2Window { scan: 100 }));
    }

    #[test]
    fn build_scan_with_no_peaks_errors() {
        let scans = vec![RawScan {
            scan_number: 100,
            retention_time_min: 20.0,
            mz_array: vec![],
            intensity_array: vec![],
        }];
        let err = build_xic3d_data(&scans, "PEPTIDEK", 2, 460.0, &[], &tol(), 100, "s.mzML")
            .unwrap_err();
        assert!(matches!(err, ReportError::AnnotationError { scan: 100, .. }));
    }
}
```

- [ ] **Step 2: 跑测试确认编译失败**

Run: `cargo test -p protein-copilot-report xic3d_build`
Expected: 编译错误 `cannot find function build_xic3d_data` / `ReportError::EmptyMs2Window`

- [ ] **Step 3a: 新增错误变体**

修改 `crates/report/src/error.rs`，在 `ReportError` 枚举中 `EmptyResult` 变体之后插入：

```rust
    /// No MS2 scans were found in the RT window around the target scan.
    #[error("no MS2 scans found in the RT window around scan {scan}")]
    EmptyMs2Window {
        /// Target scan number.
        scan: u32,
    },

    /// Per-scan annotation failed.
    #[error("annotation failed for scan {scan}: {detail}")]
    AnnotationError {
        /// Scan number that failed.
        scan: u32,
        /// Underlying error detail.
        detail: String,
    },
```

- [ ] **Step 3b: 写构建器实现（在 `xic3d_build.rs` 顶部、测试模块之前插入）**

```rust
use protein_copilot_core::search_params::{MassTolerance, Modification};
use protein_copilot_core::spectrum::{MsLevel, PrecursorInfo, Spectrum};
use protein_copilot_search_engine::annotate::annotate_spectrum;
use protein_copilot_xic::RawScan;

use crate::error::ReportError;
use crate::xic3d_types::{Ms2ScanAnnotation, Xic3dData};

/// Builds [`Xic3dData`] from raw MS2 scans by annotating each spectrum against
/// the identified peptide's theoretical b/y ions.
///
/// `ms2_scans` are the full peak lists captured by `extract_xic_unified`
/// (target isolation window, ±n_cycles), already sorted by scan number.
/// Returns [`ReportError::EmptyMs2Window`] if `ms2_scans` is empty, or
/// [`ReportError::AnnotationError`] if a spectrum fails to annotate.
#[allow(clippy::too_many_arguments)]
pub fn build_xic3d_data(
    ms2_scans: &[RawScan],
    peptide_sequence: &str,
    charge: i32,
    precursor_mz: f64,
    modifications: &[Modification],
    fragment_tolerance: &MassTolerance,
    target_scan: u32,
    source_file: &str,
) -> Result<Xic3dData, ReportError> {
    if ms2_scans.is_empty() {
        return Err(ReportError::EmptyMs2Window { scan: target_scan });
    }

    let mut scans = Vec::with_capacity(ms2_scans.len());
    for raw in ms2_scans {
        let spectrum = Spectrum::new(
            raw.scan_number,
            MsLevel::MS2,
            raw.retention_time_min,
            vec![PrecursorInfo {
                mz: precursor_mz,
                charge: Some(charge),
                intensity: None,
                isolation_window: None,
                source_scan: None,
            }],
            raw.mz_array.clone(),
            raw.intensity_array.clone(),
        )
        .map_err(|e| ReportError::AnnotationError {
            scan: raw.scan_number,
            detail: e.to_string(),
        })?;

        let annotation = annotate_spectrum(
            &spectrum,
            peptide_sequence,
            charge,
            fragment_tolerance,
            modifications,
            Vec::new(),
            false,
            false,
        )
        .map_err(|e| ReportError::AnnotationError {
            scan: raw.scan_number,
            detail: e.to_string(),
        })?;

        scans.push(Ms2ScanAnnotation {
            scan_number: raw.scan_number,
            retention_time_min: raw.retention_time_min,
            total_peaks: raw.mz_array.len(),
            is_target: raw.scan_number == target_scan,
            annotation,
        });
    }

    Ok(Xic3dData {
        peptide_sequence: peptide_sequence.to_string(),
        charge,
        precursor_mz,
        target_scan,
        source_file: source_file.to_string(),
        mz_tolerance: fragment_tolerance.clone(),
        scans,
    })
}
```

- [ ] **Step 4: 跑测试确认通过**

Run: `cargo test -p protein-copilot-report xic3d_build`
Expected: PASS（4 passed）

- [ ] **Step 5: 提交**

```bash
git add crates/report/src/xic3d_build.rs crates/report/src/error.rs crates/report/src/lib.rs
git commit -m "feat(report): build_xic3d_data — per-scan b/y annotation of window MS2

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: report crate — `templates/xic3d.html` 模板

**Files:**
- Create: `crates/report/templates/xic3d.html`

模板内容正确性由 Task 4 的渲染测试验证（注入数据后断言）。本任务只创建文件并做结构检查。

- [ ] **Step 1: 创建模板文件**

创建 `crates/report/templates/xic3d.html`，完整内容如下（占位符：`__PEPTIDE_PLACEHOLDER__`、`__PLOTLY_SRC__`、`/*__XIC3D_JSON__*/`、`/*__MAX_PEAKS_3D__*/`）：

```html
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>3D MS2 — __PEPTIDE_PLACEHOLDER__</title>
<script src="__PLOTLY_SRC__"></script>
<style>
  * { margin:0; padding:0; box-sizing:border-box; }
  body { font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif; background:#0a0e1a; color:#e0e6f0; padding:20px; }
  .container { max-width:1400px; margin:0 auto; }
  .header { background:linear-gradient(135deg,#1a2744,#0f1829); border:1px solid #2a3a5c; border-radius:12px; padding:20px; margin-bottom:16px; }
  h1 { font-size:20px; margin-bottom:6px; }
  .subtitle { color:#7a8aaa; font-size:13px; }
  .meta-grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:10px; margin-top:14px; }
  .meta-card { background:#0f1829; border:1px solid #2a3a5c; border-radius:8px; padding:10px 14px; }
  .meta-label { font-size:11px; color:#7a8aaa; text-transform:uppercase; }
  .meta-value { font-size:17px; font-weight:600; margin-top:3px; }
  .panel { background:#0f1829; border:1px solid #2a3a5c; border-radius:12px; padding:16px; margin-bottom:16px; }
  .panel-title { font-size:15px; font-weight:600; margin-bottom:12px; }
  #plot3d { width:100%; height:520px; }
  .scan-card { border:1px solid #2a3a5c; border-radius:10px; margin-bottom:12px; overflow:hidden; }
  .scan-card.target { box-shadow:0 0 0 2px #3b82f6 inset; }
  .scan-head { display:flex; justify-content:space-between; align-items:center; gap:8px; flex-wrap:wrap; padding:9px 14px; background:#13203a; border-bottom:1px solid #2a3a5c; }
  .scan-head .left { font-size:14px; }
  .scan-head .right { font-size:12px; color:#9fb0d0; }
  .badge { font-size:10px; color:#3b82f6; border:1px solid #3b82f6; border-radius:6px; padding:0 5px; margin-left:6px; }
  .chip { display:inline-block; padding:1px 7px; border-radius:10px; font-size:11px; color:#fff; margin-left:3px; }
  .spec-plot { width:100%; height:190px; padding:4px 8px; }
</style>
</head>
<body>
<div class="container">
  <div class="header">
    <h1 id="pep"></h1>
    <div class="subtitle" id="sub"></div>
    <div class="meta-grid">
      <div class="meta-card"><div class="meta-label">Charge</div><div class="meta-value" id="m-charge"></div></div>
      <div class="meta-card"><div class="meta-label">Precursor m/z</div><div class="meta-value" id="m-mz"></div></div>
      <div class="meta-card"><div class="meta-label">谱图数</div><div class="meta-value" id="m-scans"></div></div>
      <div class="meta-card"><div class="meta-label">鉴定谱图</div><div class="meta-value" id="m-target"></div></div>
      <div class="meta-card"><div class="meta-label">m/z 容差</div><div class="meta-value" id="m-tol"></div></div>
    </div>
  </div>

  <div class="panel">
    <div class="panel-title">① 3D 柱状 MS2 总览（拖动旋转；灰柱=普通峰，红=b 离子，蓝=y 离子）</div>
    <div id="plot3d"></div>
  </div>

  <div class="panel">
    <div class="panel-title">② 逐谱图 b/y 标注（每张：总谱峰数 + 匹配情况 + 标注棒状图）</div>
    <div id="scan-list"></div>
  </div>
</div>

<script>
var data = /*__XIC3D_JSON__*/;
var MAX_PEAKS_3D = /*__MAX_PEAKS_3D__*/;
(function(){
  var B_COLOR = '#e74c3c', Y_COLOR = '#3498db', GRAY = '#7f8a97';
  function ionLabel(a){ var t = (a.ion_type === 'B' ? 'b' : 'y'); var s = t + a.ion_number; if (a.charge > 1) s += '^' + a.charge + '+'; return s; }
  function ionColor(a){ return a.ion_type === 'B' ? B_COLOR : Y_COLOR; }

  function header(){
    document.getElementById('pep').textContent = '3D MS2 — ' + data.peptide_sequence;
    document.getElementById('sub').textContent = data.source_file + ' · 目标 scan ' + data.target_scan;
    document.getElementById('m-charge').textContent = data.charge + '+';
    document.getElementById('m-mz').textContent = data.precursor_mz.toFixed(4);
    document.getElementById('m-scans').textContent = data.scans.length;
    document.getElementById('m-target').textContent = 'scan ' + data.target_scan;
    document.getElementById('m-tol').textContent = data.mz_tolerance.value + ' ' + data.mz_tolerance.unit;
  }

  function build3d(){
    var gX=[],gY=[],gZ=[], bX=[],bY=[],bZ=[], yX=[],yY=[],yZ=[];
    data.scans.forEach(function(s){
      var rt = s.retention_time_min;
      var peaks = s.annotation.peaks;
      var matched = peaks.filter(function(p){ return p.annotation; });
      var others = peaks.filter(function(p){ return !p.annotation; });
      if (MAX_PEAKS_3D && others.length > MAX_PEAKS_3D) {
        others = others.slice().sort(function(a,b){ return b.intensity - a.intensity; }).slice(0, MAX_PEAKS_3D);
      }
      others.forEach(function(p){ gX.push(rt,rt,null); gY.push(p.mz,p.mz,null); gZ.push(0,p.intensity,null); });
      matched.forEach(function(p){
        if (p.annotation.ion_type === 'B') { bX.push(rt,rt,null); bY.push(p.mz,p.mz,null); bZ.push(0,p.intensity,null); }
        else { yX.push(rt,rt,null); yY.push(p.mz,p.mz,null); yZ.push(0,p.intensity,null); }
      });
    });
    var traces = [
      { type:'scatter3d', mode:'lines', x:gX, y:gY, z:gZ, name:'其他 MS2 峰', line:{ color:GRAY, width:2 }, opacity:0.55, hoverinfo:'skip' },
      { type:'scatter3d', mode:'lines', x:bX, y:bY, z:bZ, name:'b 离子', line:{ color:B_COLOR, width:6 } },
      { type:'scatter3d', mode:'lines', x:yX, y:yY, z:yZ, name:'y 离子', line:{ color:Y_COLOR, width:6 } }
    ];
    var layout = {
      margin:{ l:0, r:0, t:0, b:0 }, paper_bgcolor:'rgba(0,0,0,0)', font:{ color:'#9fb0d0', size:10 },
      showlegend:true, legend:{ x:0, y:1, font:{ size:10 }, bgcolor:'rgba(0,0,0,0)' },
      scene:{
        xaxis:{ title:'RT (min)', gridcolor:'rgba(120,140,180,0.2)' },
        yaxis:{ title:'m/z', gridcolor:'rgba(120,140,180,0.2)' },
        zaxis:{ title:'强度', gridcolor:'rgba(120,140,180,0.2)' },
        camera:{ eye:{ x:1.75, y:1.45, z:0.85 } }
      }
    };
    Plotly.newPlot('plot3d', traces, layout, { displayModeBar:false, responsive:true });
  }

  function plotScan(divId, ann){
    var gx=[],gy=[];
    ann.peaks.forEach(function(p){ if (!p.annotation) { gx.push(p.mz,p.mz,null); gy.push(0,p.intensity,null); } });
    var traces = [{ type:'scatter', mode:'lines', x:gx, y:gy, line:{ color:'#5a6b88', width:1 }, hoverinfo:'skip', showlegend:false }];
    var lx=[],ly=[],lt=[],lc=[];
    ann.peaks.forEach(function(p){
      if (p.annotation) {
        traces.push({ type:'scatter', mode:'lines', x:[p.mz,p.mz], y:[0,p.intensity], line:{ color:ionColor(p.annotation), width:2.5 }, hoverinfo:'skip', showlegend:false });
        lx.push(p.mz); ly.push(p.intensity); lt.push(ionLabel(p.annotation)); lc.push(ionColor(p.annotation));
      }
    });
    traces.push({ type:'scatter', mode:'markers+text', x:lx, y:ly, text:lt, textposition:'top center', textfont:{ size:10, color:'#cdd3da' }, marker:{ size:4, color:lc }, showlegend:false, hovertemplate:'%{text}<br>m/z %{x:.3f}<br>强度 %{y:.0f}<extra></extra>' });
    var layout = {
      margin:{ l:50, r:8, t:8, b:30 }, paper_bgcolor:'rgba(0,0,0,0)', plot_bgcolor:'rgba(0,0,0,0)', font:{ color:'#9fb0d0', size:9 },
      xaxis:{ title:'m/z', gridcolor:'rgba(120,140,180,0.12)', zeroline:false },
      yaxis:{ title:'强度', gridcolor:'rgba(120,140,180,0.12)', zeroline:false, rangemode:'tozero' }
    };
    Plotly.newPlot(divId, traces, layout, { displayModeBar:false, responsive:true });
  }

  function buildList(){
    var html = data.scans.map(function(s, i){
      var ann = s.annotation;
      var matched = ann.peaks.filter(function(p){ return p.annotation; });
      var chips = matched.map(function(p){ return '<span class="chip" style="background:' + ionColor(p.annotation) + '">' + ionLabel(p.annotation) + '</span>'; }).join('');
      var tgt = s.is_target ? '<span class="badge">鉴定谱图</span>' : '';
      return '<div class="scan-card' + (s.is_target ? ' target' : '') + '">'
        + '<div class="scan-head"><div class="left"><b>scan ' + s.scan_number + '</b>' + tgt + ' <span style="color:#7a8aaa">RT ' + s.retention_time_min.toFixed(2) + ' min</span></div>'
        + '<div class="right">总谱峰数 <b>' + s.total_peaks + '</b> · 匹配 <b>' + ann.matched_ions + '/' + ann.total_ions + '</b> ' + chips + '</div></div>'
        + '<div class="spec-plot" id="spec' + i + '"></div></div>';
    }).join('');
    document.getElementById('scan-list').innerHTML = html;
    data.scans.forEach(function(s, i){ plotScan('spec' + i, s.annotation); });
  }

  header();
  build3d();
  buildList();
})();
</script>
</body>
</html>
```

- [ ] **Step 2: 结构检查（确认 4 个占位符与关键元素存在）**

Run:
```bash
grep -c -e '__PEPTIDE_PLACEHOLDER__' -e '__PLOTLY_SRC__' -e '/\*__XIC3D_JSON__\*/' -e '/\*__MAX_PEAKS_3D__\*/' crates/report/templates/xic3d.html
grep -c -e 'id="plot3d"' -e 'id="scan-list"' crates/report/templates/xic3d.html
```
Expected: 第一条输出 `4`，第二条输出 `2`

- [ ] **Step 3: 提交**

```bash
git add crates/report/templates/xic3d.html
git commit -m "feat(report): add xic3d.html template (3D overview + per-scan b/y)

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 4: report crate — `render_xic_3d` + `ReportGenerator::render_xic_3d`

**Files:**
- Create: `crates/report/src/xic3d_visualize.rs`
- Modify: `crates/report/src/lib.rs`（模块声明 + `ReportGenerator` 方法）
- Test: `crates/report/src/xic3d_visualize.rs`（同文件 `#[cfg(test)] mod tests`）

- [ ] **Step 1: 声明模块 + 写失败测试**

修改 `crates/report/src/lib.rs`，在 `pub mod xic3d_build;` 之后新增：

```rust
pub mod xic3d_visualize;
```

创建 `crates/report/src/xic3d_visualize.rs`，先只放测试（渲染函数未定义，编译失败）：

```rust
//! 3D MS2 annotation visualization — renders [`Xic3dData`] to a self-contained
//! HTML file with Plotly.js (3D overview + per-scan b/y annotated spectra).

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};
    use protein_copilot_search_engine::matching::{
        generate_b_ions_with_charge, generate_y_ions_with_charge,
    };
    use protein_copilot_xic::{PlotlyMode, RawScan};

    fn sample_data() -> crate::xic3d_types::Xic3dData {
        let pep = "PEPTIDEK";
        let b = generate_b_ions_with_charge(pep, &[], 1);
        let y = generate_y_ions_with_charge(pep, &[], 1);
        let mut mz: Vec<f64> = Vec::new();
        mz.extend(b.iter().take(3).copied());
        mz.extend(y.iter().take(3).copied());
        mz.push(180.0);
        mz.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let intensity = vec![1000.0; mz.len()];
        let scans = vec![
            RawScan { scan_number: 100, retention_time_min: 19.9, mz_array: mz.clone(), intensity_array: intensity.clone() },
            RawScan { scan_number: 101, retention_time_min: 20.0, mz_array: mz, intensity_array: intensity },
        ];
        let tol = MassTolerance { value: 20.0, unit: ToleranceUnit::Ppm };
        crate::xic3d_build::build_xic3d_data(&scans, pep, 2, 460.0, &[], &tol, 101, "sample.mzML").unwrap()
    }

    #[test]
    fn render_creates_html_with_injected_data() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("xic3d.html");
        render_xic_3d(&sample_data(), &out, PlotlyMode::Cdn, Some(150)).unwrap();
        let html = std::fs::read_to_string(&out).unwrap();
        assert!(html.contains("PEPTIDEK"), "peptide missing");
        assert!(html.contains("plotly-2.35.2"), "plotly cdn missing");
        assert!(html.contains("\"scan_number\":101"), "scan data missing");
        assert!(html.contains("var MAX_PEAKS_3D = 150;"), "max peaks const missing");
        assert!(!html.contains("__XIC3D_JSON__"), "JSON placeholder not replaced");
        assert!(!html.contains("__PLOTLY_SRC__"), "plotly placeholder not replaced");
        assert!(!html.contains("__MAX_PEAKS_3D__"), "max peaks placeholder not replaced");
        assert!(!html.contains("__PEPTIDE_PLACEHOLDER__"), "peptide placeholder not replaced");
    }

    #[test]
    fn render_defaults_max_peaks_when_none() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("xic3d2.html");
        render_xic_3d(&sample_data(), &out, PlotlyMode::Cdn, None).unwrap();
        let html = std::fs::read_to_string(&out).unwrap();
        assert!(html.contains("var MAX_PEAKS_3D = 200;"), "default max peaks (200) missing");
    }

    #[test]
    fn render_creates_parent_dirs() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("a").join("b").join("xic3d.html");
        render_xic_3d(&sample_data(), &out, PlotlyMode::Cdn, None).unwrap();
        assert!(out.exists());
    }

    #[test]
    fn render_escapes_angle_brackets_in_data() {
        // A source_file containing angle brackets must be escaped in the
        // injected JSON (defends against </script> breakout / HTML injection).
        let pep = "PEPTIDEK";
        let b = generate_b_ions_with_charge(pep, &[], 1);
        let mut mz: Vec<f64> = b.iter().take(2).copied().collect();
        mz.push(180.0);
        mz.sort_by(|a, c| a.partial_cmp(c).unwrap());
        let intensity = vec![1000.0; mz.len()];
        let scans = vec![RawScan { scan_number: 1, retention_time_min: 1.0, mz_array: mz, intensity_array: intensity }];
        let tol = MassTolerance { value: 20.0, unit: ToleranceUnit::Ppm };
        let data = crate::xic3d_build::build_xic3d_data(&scans, pep, 2, 460.0, &[], &tol, 1, "a<script>b.mzML").unwrap();
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("x.html");
        render_xic_3d(&data, &out, PlotlyMode::Cdn, None).unwrap();
        let html = std::fs::read_to_string(&out).unwrap();
        assert!(html.contains("a\\u003cscript\\u003eb.mzML"), "angle brackets in data not escaped");
    }
}
```

- [ ] **Step 2: 跑测试确认编译失败**

Run: `cargo test -p protein-copilot-report xic3d_visualize`
Expected: 编译错误 `cannot find function render_xic_3d in this scope`

- [ ] **Step 3a: 写渲染实现（在 `xic3d_visualize.rs` 顶部、测试模块之前插入）**

```rust
use std::fs;
use std::path::Path;

use protein_copilot_xic::PlotlyMode;

use crate::error::ReportError;
use crate::xic3d_types::Xic3dData;

const PLOTLY_CDN: &str = "https://cdn.plot.ly/plotly-2.35.2.min.js";
const TEMPLATE: &str = include_str!("../templates/xic3d.html");

/// Default number of non-matched peaks per scan drawn in the 3D overview.
const DEFAULT_MAX_PEAKS_PER_SCAN_3D: usize = 200;

/// Renders [`Xic3dData`] into a standalone HTML file with Plotly.js charts:
/// a 3D MS2 overview plus per-scan b/y annotated spectra.
///
/// `max_peaks_per_scan_3d` caps the number of non-matched peaks drawn per scan
/// in the 3D overview (display-only declutter; matched b/y peaks are always
/// kept, and the per-scan annotated spectra always show all peaks). Defaults
/// to 200 when `None`.
pub fn render_xic_3d(
    data: &Xic3dData,
    output_path: &Path,
    plotly_mode: PlotlyMode,
    max_peaks_per_scan_3d: Option<usize>,
) -> Result<(), ReportError> {
    let json = serde_json::to_string(data)
        .map_err(|e| ReportError::SerializationError(e.to_string()))?;
    let json = crate::escape_json_for_html(&json);

    let plotly_src = match plotly_mode {
        PlotlyMode::Cdn | PlotlyMode::Embedded => PLOTLY_CDN.to_string(),
    };

    let max_peaks = max_peaks_per_scan_3d.unwrap_or(DEFAULT_MAX_PEAKS_PER_SCAN_3D);

    let escaped_peptide = data
        .peptide_sequence
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;");

    let html = TEMPLATE
        .replace("/*__XIC3D_JSON__*/", &json)
        .replace("/*__MAX_PEAKS_3D__*/", &max_peaks.to_string())
        .replace("__PLOTLY_SRC__", &plotly_src)
        .replace("__PEPTIDE_PLACEHOLDER__", &escaped_peptide);

    if let Some(parent) = output_path.parent() {
        fs::create_dir_all(parent).map_err(|e| ReportError::IoError {
            path: parent.to_path_buf(),
            detail: e.to_string(),
        })?;
    }

    fs::write(output_path, &html).map_err(|e| ReportError::IoError {
        path: output_path.to_path_buf(),
        detail: e.to_string(),
    })?;

    Ok(())
}
```

- [ ] **Step 3b: 在 `ReportGenerator` 暴露方法**

修改 `crates/report/src/lib.rs`，在 `impl ReportGenerator { ... }` 内 `render_unified` 方法之后插入：

```rust
    /// Renders 3D MS2 annotation data as a self-contained HTML file.
    pub fn render_xic_3d(
        data: &crate::xic3d_types::Xic3dData,
        output_path: &Path,
        plotly_mode: protein_copilot_xic::PlotlyMode,
        max_peaks_per_scan_3d: Option<usize>,
    ) -> Result<(), ReportError> {
        xic3d_visualize::render_xic_3d(data, output_path, plotly_mode, max_peaks_per_scan_3d)
    }
```

- [ ] **Step 4: 跑测试确认通过**

Run: `cargo test -p protein-copilot-report xic3d`
Expected: PASS（全部 `xic3d_*` 测试通过：types 1 项 + build 4 项 + visualize 4 项）

- [ ] **Step 5: 提交**

```bash
git add crates/report/src/xic3d_visualize.rs crates/report/src/lib.rs
git commit -m "feat(report): render_xic_3d — Plotly HTML for 3D MS2 + per-scan b/y

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 5: mcp-server — `extract_xic` 接入 `view=3d`

**Files:**
- Modify: `crates/mcp-server/src/tools.rs`
  - 新增本地 `XicView` 枚举（约 `ExtractXicInput` 定义之前，line ~446）
  - `ExtractXicInput` 新增 `view` / `max_peaks_per_scan_3d` 字段（`run_id` 字段后，line ~513）
  - `ExtractXicResult` 新增 `annotated_scan_count` / `target_matched_ions` 字段（`summary` 字段后，line ~530）
  - `extract_xic` 方法：插入 3D 分支（`let unified_result = …?;` 之后、`let xic_data = …;` 之前，line ~3138）
  - `extract_xic` 标准路径构造 `ExtractXicResult` 处补两个 `None`（line ~3170）
  - 更新 `extract_xic` 工具 description（line ~2996）

> 本任务为机械接线：核心逻辑已在 Task 2/4 单元测试覆盖。验证门槛 = 构建 + clippy（warnings 视为 error）+ 全量测试。

- [ ] **Step 1: 新增 `XicView` 枚举**

在 `ExtractXicInput` 结构体定义之前插入：

```rust
/// View mode for the `extract_xic` tool.
#[derive(Debug, Clone, Copy, Deserialize, schemars::JsonSchema)]
enum XicView {
    /// Standard 2D XIC line chart (default).
    #[serde(rename = "standard")]
    Standard,
    /// 3D MS2 overview + per-scan b/y annotated spectra.
    #[serde(rename = "3d")]
    ThreeD,
}
```

- [ ] **Step 2: `ExtractXicInput` 新增字段**

在 `ExtractXicInput` 的 `run_id: Option<String>,` 字段之后插入：

```rust
    /// View mode: standard 2D XIC (default) or 3D MS2 annotation view.
    #[serde(default)]
    #[schemars(
        description = "View mode: 'standard' (default, 2D XIC line chart) or '3d' (3D MS2 overview RT x m/z x intensity + per-scan b/y annotated spectra with total peak counts)."
    )]
    view: Option<XicView>,

    /// (3D only) Max non-matched peaks per scan drawn in the 3D overview.
    #[serde(default)]
    #[schemars(
        description = "Only for view=3d: max non-matched peaks per scan drawn in the 3D overview (display declutter; matched b/y always kept). Default 200."
    )]
    max_peaks_per_scan_3d: Option<usize>,
```

- [ ] **Step 3: `ExtractXicResult` 新增字段**

在 `ExtractXicResult` 的 `summary: String,` 字段之后插入：

```rust
    /// (3D mode only) Number of MS2 scans annotated in the window.
    #[serde(skip_serializing_if = "Option::is_none")]
    annotated_scan_count: Option<usize>,
    /// (3D mode only) Matched fragment ions in the target scan.
    #[serde(skip_serializing_if = "Option::is_none")]
    target_matched_ions: Option<u32>,
```

- [ ] **Step 4: 插入 3D 分支**

在 `extract_xic` 方法中，紧接 `let unified_result = protein_copilot_xic::extract::extract_xic_unified(…)?;`（即 `.map_err(…)?;` 这一行）之后、`let xic_data = unified_result.xic_data;` 之前，插入：

```rust
        // 3D view: annotate every MS2 scan in the window and render 3D HTML.
        if matches!(input.view, Some(XicView::ThreeD)) {
            let plotly_mode = input
                .plotly_mode
                .unwrap_or(protein_copilot_xic::PlotlyMode::Cdn);
            let out_path = input
                .output_path
                .clone()
                .map(PathBuf::from)
                .unwrap_or_else(|| PathBuf::from(format!("output/xic3d_scan{}.html", resolved_scan)));

            let data = protein_copilot_report::xic3d_build::build_xic3d_data(
                &unified_result.raw_scans.ms2_scans,
                &peptide,
                charge,
                precursor_mz,
                &modifications,
                &params.mz_tolerance,
                resolved_scan,
                &file_path.to_string_lossy(),
            )
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?;

            protein_copilot_report::xic3d_visualize::render_xic_3d(
                &data,
                &out_path,
                plotly_mode,
                input.max_peaks_per_scan_3d,
            )
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?;

            let annotated = data.scans.len();
            let target_matched = data
                .scans
                .iter()
                .find(|s| s.is_target)
                .map(|s| s.annotation.matched_ions);

            tracing::info!(scans = annotated, "completed (3d)");
            return Ok(Json(ExtractXicResult {
                output_path: out_path.to_string_lossy().to_string(),
                ms2_scan_count: annotated,
                light_trace_count: 0,
                heavy_trace_count: 0,
                has_ms1_xic: false,
                annotated_scan_count: Some(annotated),
                target_matched_ions: target_matched,
                summary: format!(
                    "3D MS2 view for {} ({}+): {} MS2 scans annotated in window",
                    peptide, charge, annotated
                ),
            }));
        }
```

- [ ] **Step 5: 标准路径补 `None` 字段**

在标准路径返回 `Ok(Json(ExtractXicResult { … }))`（line ~3170）中，`summary,` 之后补两行：

```rust
            annotated_scan_count: None,
            target_matched_ions: None,
```

- [ ] **Step 6: 更新工具 description**

把 `extract_xic` 的 `description = "…"`（line ~2996）末尾追加一句（合并进同一字符串）：

```
" Set view='3d' for a 3D MS2 overview (RT x m/z x intensity sticks) plus per-scan b/y annotated spectra with total peak counts (output: xic3d_scan{N}.html)."
```

- [ ] **Step 7: 构建 + clippy + 全量测试**

Run:
```bash
cargo build -p protein-copilot-mcp-server
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
```
Expected: 全部成功，无警告，新旧测试全部通过

- [ ] **Step 8（可选冒烟，需真实 mzML）**

如手头有 mzML 及一个 MS2 scan，可手动验证生成文件：通过 MCP 客户端调用 `extract_xic`，参数 `{"file_path":"<.mzML>","scan_number":<N>,"peptide_sequence":"<SEQ>","charge":<z>,"precursor_mz":<mz>,"view":"3d"}`，确认在 `output/xic3d_scan<N>.html` 生成且浏览器可打开（上方 3D + 下方逐谱图标注）。

- [ ] **Step 9: 提交**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "feat(mcp-server): add view=3d to extract_xic (3D MS2 + per-scan b/y)

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## 完成标准

- [ ] `cargo test --workspace` 全绿；`cargo clippy --workspace --all-targets -- -D warnings` 无警告。
- [ ] `extract_xic view=3d` 生成 `xic3d_scan{N}.html`：上方 3D 柱状 MS2 总览（灰=普通峰、红=b、蓝=y），下方逐谱图列出总谱峰数 + b/y 匹配标注棒状图，鉴定谱图高亮。
- [ ] `view` 缺省时行为与现状完全一致（回归保护）。
