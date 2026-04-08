# XIC 可视化设计文档

> **日期**: 2026-04-08
> **状态**: 设计确认（经 7 轮独立 review，v2 修订完成）
> **方案**: 方案 C — 分层架构 + Plotly.js 交互式可视化

---

## 1. 需求概述

### 1.1 核心目标

实现碎片离子级别的 XIC（Extracted Ion Chromatogram）可视化，重点支持 DIA 数据场景。用户可以：

1. 查看某个 PSM 对应的碎片离子 XIC 色谱曲线
2. 同时查看 MS1 前体离子 XIC 作为参考
3. 在轻标/重标模式下对比碎片离子 XIC（MVP: SILAC）
4. 交互式地浏览、缩放、切换离子显隐

### 1.2 用户场景

- **DIA 数据分析**：DIA 模式下每个隔离窗口在每个 cycle 都有 MS2 谱图，碎片离子 XIC 曲线平滑且信息丰富
- **轻重标对比**：SILAC 标记实验中，对比轻/重标肽段的碎片离子色谱行为
- **定性验证**：通过 XIC 共洗脱（co-elution）确认 PSM 鉴定的可靠性

---

## 2. 设计决策汇总

| 决策项 | 选择 | 理由 |
|--------|------|------|
| 数据源 | MS1 + MS2 | MS2 碎片离子 XIC 是核心，MS1 前体 XIC 作为参考 |
| 重标类型 | MVP: SILAC；v2: TMT | SILAC 影响碎片离子 m/z，TMT 仅影响 reporter 区段，复杂度不同 |
| MCP Tool | 独立 tool `extract_xic`；v2: annotate 可选嵌入 | 模块化，灵活组合 |
| 碎片离子选择 | 优先目标 scan 已匹配离子，再按强度 top-N，UI 可切换 | 已匹配离子比纯强度排序更可靠（避免 DIA 干扰） |
| 轻重标布局 | 叠加 / 上下对比，可切换 | 两种模式适用不同分析场景 |
| 交互级别 | tooltip + 缩放平移 + 离子开关切换 | Plotly 原生支持，开发量可控 |
| RT 范围 | 前后 N 个同隔离窗口 cycle（默认 ±5） | DIA 场景下更直观 |
| 质量偏移 | 自动计算 + 用户可覆盖 | 根据肽段序列 + SILAC 标记自动推算 |
| m/z 提取窗口 | 复用 fragment_tolerance，可覆盖 | 保持一致性 |
| 强度提取规则 | 窗口内最高峰（max），缺失 scan 零填充 | 简单可靠，对齐所有 trace 到相同 scan 集合 |
| HTML 输出 | 独立 HTML（MVP）；合并模式 v2 | 先做好独立功能 |
| JS 可视化 | Plotly.js (CDN / 嵌入可配置) | 科学可视化标准库，交互功能完备 |
| 文件格式 | MVP: 仅 mzML | MGF 无 MS1/隔离窗口信息，不适合 XIC |
| 修饰处理 | 传入完整修饰列表（fixed + variable + applied） | 碎片离子 m/z 依赖实际修饰状态 |
| 碎片离子生成 | 复用 `matching.rs` 公开 API | 避免第三次重复实现 |
| 文件解析 | 1.5-pass 流式提取（Pass 0 定位 + Pass 1 流式遍历） | 避免全量加载，同时保证前后 scan 完整 |

---

## 3. 架构设计

### 3.1 系统分层

```text
┌─────────────────────────────────────────────────────┐
│  MCP Tool Layer (mcp-server)                        │
│  ┌──────────────┐  ┌────────────────────────────┐   │
│  │ extract_xic  │  │ annotate_spectrum (扩展)    │   │
│  │ (独立 tool)  │  │ + optional XIC 嵌入        │   │
│  └──────┬───────┘  └─────────────┬──────────────┘   │
│         │                        │                   │
├─────────┼────────────────────────┼───────────────────┤
│  Library Crates                  │                   │
│  ┌──────▼───────┐  ┌────────────▼──────────────┐    │
│  │  crates/xic  │  │  crates/report (扩展)     │    │
│  │              │  │  + xic.html template       │    │
│  │ - MS1 XIC    │  │  + render_xic_html()      │    │
│  │ - MS2 XIC    │  │  + render_combined_html() │    │
│  │ - Heavy calc │  └───────────────────────────┘    │
│  └──────────────┘                                    │
│         ▲                                            │
│  ┌──────┴───────┐  ┌──────────────────┐             │
│  │ spectrum-io  │  │      core        │             │
│  │ (读取谱图)   │  │ (数据结构)       │             │
│  └──────────────┘  └──────────────────┘             │
└─────────────────────────────────────────────────────┘
```

### 3.2 新增 crate: `crates/xic`

纯 library crate，负责所有 XIC 提取的确定性计算。

```rust
// crates/xic/src/lib.rs

pub mod extract;    // XIC 数据提取
pub mod heavy;      // 重标质量偏移计算
pub mod error;      // 错误类型

// 核心数据结构
/// 单条 XIC 色谱曲线
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XicTrace {
    /// 离子标识（如 "y5¹⁺", "b3²⁺"）
    pub ion_label: String,
    /// 离子类型
    pub ion_type: IonType,
    /// 离子编号
    pub ion_number: u32,
    /// 电荷态
    pub charge: u32,
    /// 理论 m/z
    pub theoretical_mz: f64,
    /// 数据点: (retention_time_sec, intensity)
    pub data_points: Vec<XicDataPoint>,
    /// 是否为重标离子
    pub is_heavy: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XicDataPoint {
    pub retention_time_sec: f64,
    pub scan_number: u32,
    pub intensity: f64,
}

/// 完整 XIC 数据包（传给前端渲染）
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XicData {
    /// 目标肽段序列
    pub peptide_sequence: String,
    /// 目标 scan 的 RT
    pub target_rt_sec: f64,
    /// 目标 scan number
    pub target_scan: u32,
    /// 电荷态
    pub charge: i32,
    /// 真实前体 m/z（来自 PSM 或用户输入，非 DIA 隔离窗口中心）
    pub precursor_mz: f64,
    /// MS1 前体 XIC（轻标）
    pub ms1_precursor_xic: Option<XicTrace>,
    /// MS1 前体 XIC（重标，如有）
    pub ms1_heavy_precursor_xic: Option<XicTrace>,
    /// MS2 碎片离子 XIC（轻标）
    pub fragment_xic_traces: Vec<XicTrace>,
    /// MS2 碎片离子 XIC（重标）
    pub heavy_fragment_xic_traces: Vec<XicTrace>,
    /// 提取参数
    pub extraction_params: ExtractionParams,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExtractionParams {
    pub mz_tolerance: MassTolerance,
    pub n_cycles: u32,           // 前后 cycle 数
    pub top_n_ions: usize,       // 默认展示 top-N 离子
    pub label_type: Option<LabelType>,
    pub intensity_rule: IntensityRule,
}

/// 强度提取规则
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default)]
pub enum IntensityRule {
    /// 窗口内最高峰（默认）
    #[default]
    MaxInWindow,
    /// 窗口内峰面积之和
    SumInWindow,
    /// 最近峰
    NearestPeak,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum LabelType {
    /// SILAC: K+8.014 / R+10.008
    Silac {
        heavy_k_delta: f64,  // 默认 8.014199
        heavy_r_delta: f64,  // 默认 10.008269
    },
    /// 自定义质量偏移（如非标准 SILAC 标签）
    Custom {
        residue_deltas: Vec<(char, f64)>,
    },
}
```

**Review 修正说明**：
- `precursor_mz` 字段明确为"真实前体 m/z"（来自 PSM 或用户输入），而非 DIA 隔离窗口中心
- `ms1_precursor_xic` 改为 `Option`——DIA 模式下 MS1 中可能找不到前体
- 移除 TMT（v2 再支持），移除 `ScanAnnotationRef`（点击跳转 MS2 移出 MVP）
- 新增 `IntensityRule` 枚举，明确强度提取策略

### 3.3 XIC 提取算法

#### 1.5-pass 提取策略

**关键设计**：完全单次遍历不可行（需要 target scan 的 RT 和 isolation window 才能筛选其他 scan）。
采用 "1.5-pass" 策略：
- **Pass 0**（廉价）：`read_spectrum(target_scan)` 获取 RT、isolation window（仅读一个 scan）
- **Pass 1**（流式）：带着已知的 target RT + isolation window 遍历全文件，同时提取所有目标离子强度

```
输入: spectrum_file (mzML), target_scan, peptide_sequence, charge,
      modifications (完整修饰列表: fixed + applied variable),
      precursor_mz (真实前体 m/z), extraction_params

输出: XicData

步骤:
1. Pass 0 — 获取 target scan 信息（单次 read_spectrum 调用）:
   a. read_spectrum(target_scan) → target_rt, target_isolation_window
   b. 计算 RT 范围估计（用于粗筛 MS1 scan）

2. 前置计算（无需文件 I/O）:
   a. 生成理论碎片离子 m/z 列表（复用 matching::generate_b_ions_with_charge / generate_y_ions_with_charge）
   b. 如有 SILAC: 计算重标碎片离子 m/z（heavy::compute_heavy_shifts）
   c. 合并所有目标 m/z 列表（轻标 + 重标 + 前体 + 重标前体）

3. Pass 1 — 流式遍历文件（轻标 + 重标同时提取）:
   对每个 spectrum:
   a. 如果是 MS1 且 RT 在范围内:
      - 提取前体 m/z ± tolerance 内的最高峰 → ms1_precursor_xic
      - 如有重标: 提取重标前体 m/z → ms1_heavy_precursor_xic
   b. 如果是 MS2 且 isolation_window 匹配:
      - 对每个目标碎片离子 m/z: 提取 tolerance 内最高峰
      - 同时提取轻标和重标碎片离子强度
      - 记录 (scan_number, RT, intensity) 数据点

4. 后处理:
   a. 按 scan_number 排序，取 target_scan 前后 N 个同窗口 cycle
   b. 零填充：确保每条 trace 在相同的 scan 集合上对齐
   c. 选择 top-N 离子：优先选目标 scan 中已匹配的碎片离子，再按总强度排序
```

#### 隔离窗口匹配

```rust
/// 判断两个 isolation window 是否覆盖相同区域
fn same_isolation_window(a: &IsolationWindow, b: &IsolationWindow) -> bool {
    // 比较完整窗口区间，而非仅 target_mz
    let a_lo = a.target_mz - a.lower_offset;
    let a_hi = a.target_mz + a.upper_offset;
    let b_lo = b.target_mz - b.lower_offset;
    let b_hi = b.target_mz + b.upper_offset;
    // 区间中心差 < 1 Da 且宽度差 < 20%
    let center_close = ((a_lo + a_hi) / 2.0 - (b_lo + b_hi) / 2.0).abs() < 1.0;
    let width_a = a_hi - a_lo;
    let width_b = b_hi - b_lo;
    let width_close = width_a > 0.0 && ((width_a - width_b).abs() / width_a) < 0.2;
    center_close && width_close
}
```

#### 强度提取规则

```rust
/// 在实验谱图的 mz_array / intensity_array 中提取目标 m/z 的强度
fn extract_intensity(
    target_mz: f64,
    exp_mz: &[f64],
    exp_int: &[f64],
    tolerance: &MassTolerance,
    rule: IntensityRule,
) -> f64 {
    match rule {
        IntensityRule::MaxInWindow => {
            // 找 tolerance 窗口内的最高峰
            // binary search + neighbor check
        }
        IntensityRule::SumInWindow => { /* 窗口内求和 */ }
        IntensityRule::NearestPeak => { /* 最近峰 */ }
    }
}
// 如果窗口内无峰，返回 0.0（零填充）
```

### 3.4 重标质量偏移计算（MVP: SILAC only）

SILAC 影响碎片离子 m/z（含标记氨基酸的碎片离子质量增加），这是 XIC 核心需求。

#### SILAC 碎片离子偏移

```
对于 b-ion (N-terminal fragment):
  heavy_mz = light_mz + (count_K_in_prefix × 8.014199 + count_R_in_prefix × 10.008269) / charge

对于 y-ion (C-terminal fragment):
  heavy_mz = light_mz + (count_K_in_suffix × 8.014199 + count_R_in_suffix × 10.008269) / charge
```

```rust
/// 计算 SILAC 重标碎片离子的 m/z
/// 
/// 与搜索引擎中的修饰应用逻辑类似：
/// 遍历肽段序列，对包含 K 或 R 的碎片离子累加 delta
pub fn compute_heavy_fragment_mz(
    light_fragments: &[FragmentIon],
    peptide_sequence: &str,
    label: &SilacLabel,
) -> Vec<FragmentIon> {
    // 对每个 b-ion: 检查 prefix[1..=n] 中有多少 K/R
    // 对每个 y-ion: 检查 suffix[len-n..] 中有多少 K/R
    // heavy_mz = light_mz + (count_K * delta_K + count_R * delta_R) / charge
}

/// SILAC 重标前体 m/z
pub fn compute_heavy_precursor_mz(
    light_mz: f64,
    charge: i32,
    peptide_sequence: &str,
    label: &SilacLabel,
) -> f64 {
    let count_k = peptide_sequence.chars().filter(|&c| c == 'K').count() as f64;
    let count_r = peptide_sequence.chars().filter(|&c| c == 'R').count() as f64;
    light_mz + (count_k * label.heavy_k_delta + count_r * label.heavy_r_delta) / charge as f64
}
```

**v2 预留：TMT**
TMT 的 reporter ions (126–134 Da) 在低 m/z 区，提取逻辑完全不同（不影响碎片离子 m/z），留待 v2 专门设计。

### 3.5 HTML 可视化模板

#### 页面布局

```
┌──────────────────────────────────────────────────┐
│  Header: 肽段信息 (序列, 电荷, m/z, score)       │
├──────────────────────────────────────────────────┤
│  [模式切换] 叠加 | 上下对比    [标记] 轻 | 重    │
├──────────────────────────────────────────────────┤
│  MS1 前体 XIC (参考)                              │
│  ┌──────────────────────────────────────────┐    │
│  │  ──── light precursor                    │    │
│  │  ---- heavy precursor                    │    │
│  │  ▼ target RT                             │    │
│  └──────────────────────────────────────────┘    │
├──────────────────────────────────────────────────┤
│  MS2 碎片离子 XIC (主体)                          │
│  ┌──────────────────────────────────────────┐    │
│  │  ── y5¹⁺  ── y4¹⁺  ── b3¹⁺  ...       │    │
│  │  (每条 trace 不同颜色，可切换显示/隐藏)    │    │
│  │  ▼ target RT                             │    │
│  └──────────────────────────────────────────┘    │
├──────────────────────────────────────────────────┤
│  [离子开关面板]                                    │
│  ☑ y5¹⁺  ☑ y4¹⁺  ☑ b3¹⁺  ☐ y3¹⁺  ☐ b2¹⁺ ... │
├──────────────────────────────────────────────────┤
│  提取参数摘要 + 数据质量指标                       │
└──────────────────────────────────────────────────┘
```

**MVP 移除**：点击 RT → 展示 MS2 annotation 功能（静态 HTML 无法动态加载数据，v2 考虑 WebSocket 方案）。

#### Plotly.js 集成

- **加载方式可配置**：CDN（默认）或嵌入 plotly-basic.min.js
- **子图布局**：`Plotly.newPlot()` with `subplots` 实现 MS1 + MS2 上下排列
- **Legend 交互**：点击 legend 切换离子显示/隐藏（Plotly 原生支持）
- **Hover**：显示 RT、m/z、intensity、scan number

#### 轻重标切换

- **叠加模式**：同一子图中，轻标实线 + 重标虚线，颜色对应
- **上下对比模式**：3 个子图（MS1 / Light MS2 / Heavy MS2），共享 x 轴
- **切换按钮**：在 HTML 中用 JS 控制 Plotly layout 更新

---

## 4. MCP Tool 设计

### 4.1 `extract_xic` Tool（独立）

```rust
#[tool(name = "extract_xic")]
struct ExtractXicInput {
    /// 谱图文件路径（仅支持 mzML）
    file_path: String,
    /// 目标 scan number
    scan_number: u32,
    /// 肽段序列
    peptide_sequence: String,
    /// 电荷态
    charge: i32,
    /// 真实前体 m/z（对 DIA 数据，不能用隔离窗口中心，需来自 PSM 或用户指定）
    precursor_mz: f64,
    /// 完整修饰列表（fixed + 该 PSM 实际应用的 variable modifications）
    modifications: Option<Vec<Modification>>,
    /// 前后 cycle 数（默认 5）
    n_cycles: Option<u32>,
    /// 展示 top-N 碎片离子（默认 6）
    top_n_ions: Option<usize>,
    /// 重标类型（可选）
    label_type: Option<LabelType>,
    /// m/z 提取容差（默认复用 fragment_tolerance）
    extraction_tolerance: Option<MassTolerance>,
    /// 强度提取规则（默认 MaxInWindow）
    intensity_rule: Option<IntensityRule>,
    /// Plotly 加载模式
    plotly_mode: Option<PlotlyMode>,  // Cdn | Embedded
    /// 输出路径
    output_path: Option<String>,
    /// 也可从已有搜索结果获取上下文（自动填入 peptide, charge, mods, precursor_mz）
    /// MVP 限制：仅支持单文件搜索的 run_id（PSM 无 source_file 字段）
    run_id: Option<String>,
}
```

### 4.2 `annotate_spectrum` 扩展（v2）

在现有 `AnnotateSpectrumInput` 中新增可选字段（**v2 功能**，MVP 不实现）：

```rust
/// 是否同时生成 XIC（嵌入到 annotation HTML 中）
include_xic: Option<bool>,       // 默认 false
/// XIC 相关参数（仅当 include_xic=true 时使用）
xic_n_cycles: Option<u32>,
xic_label_type: Option<LabelType>,
```

---

## 5. 文件结构

```
crates/
├── xic/                          ← 新增 crate
│   ├── Cargo.toml
│   └── src/
│       ├── lib.rs                ← 公开 API + 数据结构
│       ├── extract.rs            ← XIC 提取核心逻辑
│       ├── heavy.rs              ← 重标质量偏移计算
│       └── error.rs              ← XicError
├── report/
│   ├── src/
│   │   ├── visualize.rs          ← 扩展: render_xic_html()
│   │   └── ...
│   └── templates/
│       ├── annotation.html       ← 现有
│       └── xic.html              ← 新增: XIC 可视化模板
└── mcp-server/
    └── src/
        └── tools.rs              ← 新增 extract_xic tool
```

---

## 6. 数据流（1.5-pass 提取）

```
用户请求 extract_xic(file, scan=1234, peptide="PEPTIDEK", charge=3, precursor_mz=500.25)
  │
  ▼
MCP Tool: extract_xic
  │
  ├─ Pass 0: spectrum-io::read_spectrum(scan=1234)
  │     → target_rt, target_isolation_window
  │
  ├─ 前置计算（无 I/O）
  │     ├─ matching::generate_b_ions_with_charge(peptide, mods, charge)
  │     ├─ matching::generate_y_ions_with_charge(peptide, mods, charge)
  │     ├─ [SILAC] heavy::compute_heavy_fragment_mz(light_ions, peptide, label)
  │     └─ 合并目标 m/z 列表（轻 + 重 + 前体）
  │
  ├─ Pass 1: spectrum-io::for_each_spectrum() 流式遍历
  │     ├─ MS1 (RT 范围内) → 提取前体 XIC + 重标前体 XIC
  │     └─ MS2 (同窗口) → 提取所有目标离子强度 (轻+重同时)
  │
  ├─ 后处理
  │     ├─ 取 target_scan 前后 N 个同窗口 cycle
  │     ├─ 零填充缺失 scan
  │     └─ 选择 top-N 离子（优先已匹配离子）
  │
  ├─ 组装 XicData JSON
  │
  ├─ report crate: render_xic_html(xic_data) → HTML 文件
  │
  └─ 返回 XicResult { output_path, summary }
```

---

## 7. 边界条件与限制

### 7.1 已识别的边界条件

- **DDA 数据**：同一前体可能只有 1-2 个 MS2 scan，XIC 曲线稀疏 → 显示警告
- **大文件**：单次流式遍历，峰值内存 = 目标离子数 × scan 数（每点仅几十字节），可接受
- **隔离窗口不匹配**：overlapping DIA windows → 用完整区间比较（中心差 < 1 Da 且宽度差 < 20%）
- **重标肽段未检测到**：重标 XIC 全为 0 → 在 UI 中显示提示
- **DIA 前体 m/z**：隔离窗口中心 ≠ 真实前体 m/z → 要求用户输入或从 PSM 获取
- **MGF 文件**：无 MS1/隔离窗口信息 → 返回明确错误提示（需 mzML）
- **修饰缺失**：碎片离子 m/z 依赖修饰 → 未传入修饰时使用未修饰序列并发出警告
- **前体不在 MS1 中**：DIA 模式下 MS1 稀疏或前体弱 → `ms1_precursor_xic = None`
- **空 XIC**：目标离子在所有 scan 中均无匹配峰 → 该 trace 全零，top-N 不选入

### 7.2 MVP 范围

**包含**:
- MS2 碎片离子 XIC 提取与可视化（DIA + DDA）
- MS1 前体 XIC 参考（可选，可能为空）
- SILAC 重标支持（K+8.014199, R+10.008269）
- Plotly.js 交互（CDN + 嵌入，tooltip + 缩放 + legend 切换）
- 独立 MCP tool `extract_xic`
- 自动 top-N 离子选择（优先已匹配离子）
- 完整修饰支持（fixed + variable）
- 碎片离子复用 `matching.rs` API
- 仅 mzML 文件格式
- 1.5-pass 流式文件遍历（Pass 0 定位 + Pass 1 提取）

**推迟 (v2)**:
- TMT reporter ion XIC（提取逻辑完全不同）
- annotate_spectrum XIC 嵌入模式
- 点击 RT → 内嵌 MS2 annotation（静态 HTML 限制）
- MGF 格式支持（缺少必要数据）

---

## 8. 碎片离子生成复用策略

**核心原则**：碎片离子 m/z 计算不允许第三次重复实现。

已有实现：
1. `matching.rs` — `generate_b_ions_with_charge()` / `generate_y_ions_with_charge()`（搜索打分用）
2. `annotate.rs` — `generate_b_entries()` / `generate_y_entries()`（谱图标注用）

XIC crate 策略：
- **直接依赖 `search-engine` crate**，调用 `matching.rs` 中的公开 API
- 需要确保 `generate_b_ions_with_charge` / `generate_y_ions_with_charge` 是 `pub` 可见的
- XIC 只做"给定 m/z 列表 + 谱图数据 → 提取强度"，不做离子计算

```rust
// xic crate 中的调用方式
use protein_copilot_search_engine::matching::{
    generate_b_ions_with_charge,
    generate_y_ions_with_charge,
};

fn build_target_mz_list(
    peptide: &str,
    mods: &[Modification],
    charge: i32,
) -> Vec<TargetIon> {
    let max_frag_charge = if charge >= 3 { 2 } else { 1 };
    let b_ions = generate_b_ions_with_charge(peptide, mods, max_frag_charge);
    let y_ions = generate_y_ions_with_charge(peptide, mods, max_frag_charge);
    // 转换为 TargetIon { label, mz, ion_type, ion_number, charge }
}
```

---

## 9. 实现前置条件

### 9.1 spectrum-io 流式遍历 API（必需）

当前 `SpectrumReader` 仅有 `read_all()`（全量加载）和 `read_spectrum()`（单 scan）。
XIC 需要遍历所有谱图但不需要全部驻留内存。

**新增 API**：
```rust
/// 流式遍历谱图，每个 spectrum 调用 callback
fn for_each_spectrum<F>(&self, path: &Path, f: F) -> Result<(), SpectrumIoError>
where
    F: FnMut(Spectrum) -> ControlFlow<()>;
```

这也是 DIA 大文件场景的通用优化，非 XIC 特有需求。

### 9.2 HTML 数据注入改进

当前 `render_annotation_html()` 直接注入原始 JSON 到 `<script>` 中。XIC 数据量更大，改用安全模式：

```html
<!-- 旧方式（annotation.html 沿用，不改） -->
<script>const data = {{DATA_PLACEHOLDER}};</script>

<!-- 新方式（xic.html） -->
<script type="application/json" id="xic-data">{{DATA_PLACEHOLDER}}</script>
<script>
  const data = JSON.parse(document.getElementById('xic-data').textContent);
</script>
```

### 9.3 `run_id` 模式限制

当前 `Psm` 无 `source_file` 字段，`run_id` 模式仅取 `input_files[0]`。
MVP 限制：`run_id` 仅支持单文件搜索。多文件搜索需 v2 为 PSM 添加 source file 追踪。

### 9.4 位点特异性可变修饰（已知限制）

当前 `Modification` 类型仅描述修饰类别（name、residues、position class），不记录肽段中的具体修饰位点。
`mod_delta_fragment()` 对碎片中所有匹配残基应用修饰，对固定修饰和 SILAC 标记正确，
但对位点特异性可变修饰（如"第 2 个 M 氧化而非第 1 个"）会产生错误的碎片离子 m/z。

**MVP 影响**：低。MVP 核心场景为 SILAC（K/R 全标记）和固定修饰（C 全标记），不受影响。
**v2 修复**：引入 `AppliedModification { position_index: usize, mass_delta: f64 }` 表示 PSM 级别的位点特异性修饰，
并在碎片离子生成 API 中支持。
