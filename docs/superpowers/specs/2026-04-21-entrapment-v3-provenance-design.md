# Entrapment v3 设计文档 — 碎片离子溯源 + 修饰感知比对 + Mirror Plot 可视化

> **日期**：2026-04-21
> **状态**：Design Approved
> **前置**：entrapment v1+v2 已合并到 master（L0-L4 分级 + Levenshtein + k-mer + SubstitutionType）

---

## 1. 概述

### 1.1 目标

v3 的核心目标是**碎片离子溯源（Fragment Ion Provenance）**：对每条 trap PSM 的实际 MS2 谱图中的碎片离子，判断它到底来源于 trap 肽段还是共洗脱的 target 同源肽段。

具体实现三个子系统：

1. **碎片离子溯源引擎**（A）— 逐峰分类为 trap-only / target-only / shared / unassigned，计算 shared_ratio，判定嵌合谱
2. **修饰感知比对**（B）— 将修饰质量差纳入 delta_mass 计算，提升分级精度
3. **Mirror Plot 可视化**（C）— trap 碎片在上、target 碎片在下的交互式 Plotly.js 图表

### 1.2 非目标

- ML 模型训练与特征提取（独立模块）
- 定量层面分析（L/H ratio 提取等，独立模块）
- c/z 离子、a 离子、中性丢失等高级离子类型（v3 只做 b/y）
- 修饰感知的编辑距离算法改造（v3 只在 delta_mass 计算中叠加修饰质量差）

---

## 2. 架构

### 2.1 数据流

```text
┌──────────┐   ┌──────────┐   ┌──────────────┐   ┌───────────────────┐
│ search   │   │ target   │   │   mzML 谱图   │   │ entrapment YAML  │
│ results  │   │ FASTA    │   │  (v3 新增!)    │   │   config         │
│ parquet  │   │          │   │               │   │                  │
└────┬─────┘   └────┬─────┘   └──────┬────────┘   └────────┬─────────┘
     │              │                │                      │
     ▼              ▼                │                      ▼
┌─────────────────────────┐          │          ┌──────────────────┐
│ v1+v2: classify         │          │          │ config + mods    │
│  L0/L1/L2/L3/L4        │          │          │ (v3: mod masses) │
│  + SubstitutionType     │          │          └────────┬─────────┘
│  + edit distance        │          │                   │
│  + mod-aware Δm (v3)   │◄─────────────────────────────┘
└───────────┬─────────────┘          │
            │                        │
            │ classified PSMs        │
            ▼                        ▼
┌─────────────────────────────────────────────────────────────────┐
│ v3 NEW: Fragment Provenance Engine                              │
│                                                                  │
│  对每条 L2/L3/L4 trap PSM:                                       │
│                                                                  │
│  1. 读取实际 MS2 谱图 (mzML, scan_number)                        │
│  2. 生成 trap 肽段理论碎片 (b/y ions, charge 1..max)             │
│  3. 生成 target 同源肽段理论碎片 (best_target_peptide)            │
│  4. 逐峰匹配:                                                    │
│     observed peak → 匹配 trap 理论峰?  → TrapOnly                │
│                   → 匹配 target 理论峰? → TargetOnly              │
│                   → 两者都匹配?         → Shared                  │
│                   → 都不匹配?           → Unassigned              │
│  5. 计算 shared_ratio = shared / (trap_matched + target_matched) │
│  6. 判断 chimera: shared_ratio > config.chimera_threshold        │
│                                                                  │
└───────────────────────┬──────────────────────────────────────────┘
                        │
            ┌───────────┼───────────┐
            ▼           ▼           ▼
┌──────────────┐ ┌───────────┐ ┌──────────────────┐
│ TSV output   │ │ 单条 HTML │ │ entrapment 报告  │
│ + provenance │ │ mirror    │ │ 嵌入 mirror plot │
│ columns      │ │ plot      │ │ + chimera 统计   │
└──────────────┘ └───────────┘ └──────────────────┘
```

### 2.2 模块结构

```text
crates/entrapment-analysis/src/
├── provenance.rs        ← v3 新增：碎片离子溯源引擎
│   ├── PeakProvenance enum
│   ├── AnnotatedPeak struct
│   ├── FragmentProvenance struct
│   ├── trace_provenance() — 单条 PSM 溯源
│   └── trace_batch() — 批量溯源（L2/L3/L4）
│
├── mirror_plot.rs       ← v3 新增：Mirror Plot HTML 生成
│   ├── render_mirror_plot() — 生成单条 PSM 的 Plotly.js mirror HTML
│   └── embed_mirror_plots() — 将 mirror plots 嵌入 entrapment 报告
│
├── mod_parser.rs        ← v3 新增：修饰序列解析
│   ├── parse_modified_sequence() — 解析 UniMod 格式修饰
│   └── ModificationInfo struct
│
├── similarity.rs        ← v3 修改：修饰感知 delta_mass
│   ├── hamming_diff() — 修饰 Δm 叠加
│   └── classify_single() — 传入修饰信息
│
├── types.rs             ← v3 修改：UnifiedPsm 扩展
│   └── modifications: Vec<(usize, f64)>
│
├── loader/              ← v3 修改：DIA-NN loader 解析 Modified.Sequence
│   └── diann_parquet.rs
│
├── output.rs            ← v3 修改：TSV 新增 5 列
├── report.rs            ← v3 修改：HTML 报告含 mirror plot + chimera 统计
└── config.rs            ← v3 修改：新增 provenance 配置段
```

### 2.3 依赖关系

```text
entrapment-analysis ──→ core (数据结构)
                    ──→ spectrum-io (v3 新增：读取 mzML 谱图)
                    ──→ search-engine (复用 b/y 离子生成 + within_tolerance)

entrapment-cli ──→ entrapment-analysis
mcp-server     ──→ entrapment-analysis
```

v3 新增对 `spectrum-io` 和 `search-engine` 的依赖（碎片离子生成和匹配逻辑复用）。

---

## 3. 核心数据结构

### 3.1 碎片离子分类标签

```rust
/// 每个实测峰的来源分类
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PeakProvenance {
    TrapOnly,      // 只匹配 trap 肽段的理论碎片
    TargetOnly,    // 只匹配 target 同源肽段的理论碎片
    Shared,        // 两者都匹配（共享碎片）
    Unassigned,    // 两者都不匹配（噪声/其他来源）
}
```

### 3.2 单峰溯源结果

```rust
/// 单个峰的溯源结果
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnnotatedPeak {
    pub observed_mz: f64,
    pub observed_intensity: f64,
    pub provenance: PeakProvenance,
    pub trap_ion: Option<String>,      // e.g. "b4+1", "y7+2"
    pub target_ion: Option<String>,    // e.g. "b4+1", "y6+1"
    pub delta_mz_trap: Option<f64>,    // 与 trap 理论峰的偏差 (ppm)
    pub delta_mz_target: Option<f64>,  // 与 target 理论峰的偏差 (ppm)
}
```

### 3.3 PSM 级别溯源摘要

```rust
/// 一条 PSM 的碎片离子溯源结果
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FragmentProvenance {
    // 输入标识
    pub scan_number: u32,
    pub spectrum_file: String,
    pub trap_peptide: String,
    pub target_peptide: Option<String>,

    // 逐峰结果
    pub peaks: Vec<AnnotatedPeak>,

    // 统计摘要
    pub total_peaks: u32,
    pub trap_only_count: u32,
    pub target_only_count: u32,
    pub shared_count: u32,
    pub unassigned_count: u32,

    // 关键指标
    pub shared_ratio: f64,            // shared / (trap_matched + target_matched)
    pub trap_explained_ratio: f64,    // (trap_only + shared) / total_assigned
    pub is_chimera: bool,             // shared_ratio > config.chimera_threshold
}
```

---

## 4. 修饰感知比对

### 4.1 修饰解析

从 DIA-NN 的 `Modified.Sequence` 列解析 UniMod 格式修饰：

```text
输入: "AAAC(UniMod:4)DFK"
输出: peptide = "AAACDFK", modifications = [(3, 57.021464)]  // Cys carbamidomethyl
```

新增 `mod_parser.rs` 模块：

```rust
pub struct ModificationInfo {
    pub position: usize,  // 0-based residue index
    pub delta_mass: f64,  // 修饰质量差 (Da)
    pub name: Option<String>,  // UniMod name if known
}

/// 解析 Modified.Sequence 为 stripped sequence + modifications
pub fn parse_modified_sequence(modified_seq: &str) -> (String, Vec<ModificationInfo>);
```

### 4.2 Delta Mass 叠加

在 `similarity.rs` 的 `hamming_diff()` 和 Levenshtein 比对中：

```text
现有（v2）:
  delta_mass = Σ residue_mass(trap[i]) - Σ residue_mass(target[i])

v3 升级:
  trap_mass   = Σ residue_mass(trap[i])   + Σ mod_delta(trap modifications)
  target_mass = Σ residue_mass(target[i]) + Σ mod_delta(target modifications)
  delta_mass  = trap_mass - target_mass
```

target 肽段通常无修饰（来自 in-silico digest），因此 target_mod_delta = 0。效果是 trap 上的修饰质量会加入到 delta_mass 中，使分级更准确。

### 4.3 向后兼容

- `UnifiedPsm.modifications` 为空时，行为与 v2 完全一致
- Generic TSV loader 不解析修饰（字段为空），仅 DIA-NN loader 填充
- 配置文件无需修改

---

## 5. Mirror Plot 可视化

### 5.1 布局

```text
          trap 碎片 (向上)
    ▲  b3   b4    y8    b5    y7    y5    y4    b8    y3
    │  ██   ██    ██    ██    ██    ██    ██    ██    ██
────┼──────────────────────────────────────────────────── m/z →
    │  ██   ██    ██    ██    ██    ██         ██    ██
    ▼  b4   b5    y8    y6    y7    y4         b9    y3
          target 碎片 (向下)

颜色：■ shared (黄)  ■ trap-only (橙)  ■ target-only (青)
```

### 5.2 实现

使用 Plotly.js bar chart：
- trap 碎片强度为正值（向上），target 碎片强度取负值（向下）
- 每个 bar 按 PeakProvenance 着色
- Hover 显示：离子标签、m/z、强度、Δppm
- 标题栏显示：trap 肽段 vs target 肽段、Δm、shared_ratio、chimera 判定

### 5.3 两种输出模式

1. **单条 HTML**（`annotate_provenance` MCP tool）— 独立 HTML 文件，类似现有 `annotate_spectrum`
2. **嵌入报告**（批量模式）— 在 entrapment HTML 报告的 PSM 表格中，点击 L2/L3/L4 行展开内嵌 mirror plot

---

## 6. 输出扩展

### 6.1 TSV 新增列

| 列名 | 类型 | 说明 |
|------|------|------|
| `shared_ratio` | f64 | shared / (trap_matched + target_matched)，无溯源时为空 |
| `trap_only_ions` | u32 | trap 独有的匹配峰数 |
| `target_only_ions` | u32 | target 独有的匹配峰数 |
| `shared_ions` | u32 | trap 和 target 共享的匹配峰数 |
| `is_chimera` | bool | shared_ratio > chimera_threshold |

### 6.2 HTML 报告增强

- 统计面板新增：chimera 计数、shared_ratio 分布直方图
- PSM 表格新增 5 列
- L2/L3/L4 行可展开显示内嵌 mirror plot（Plotly.js 懒加载）

---

## 7. MCP Tools

### 7.1 新增 Tools

| Tool | 输入 | 输出 | 说明 |
|------|------|------|------|
| `trace_fragment_provenance` | classified_file, mzml_dir, config_file, target_fasta | provenance TSV + 统计摘要 JSON | 批量模式：对 L2/L3/L4 做碎片溯源 |
| `annotate_provenance` | classified_file, mzml_path, scan_number, config_file, target_fasta | mirror plot HTML 路径 | 单条交互：生成单条 PSM 的 mirror plot |

### 7.2 现有 Tool 修改

`classify_entrapment_hits` 新增可选参数：
- `mzml_dir: Option<String>` — 提供时自动对 L2/L3/L4 做碎片溯源
- 输出中新增 provenance 列和 chimera 统计

---

## 8. 配置扩展

```yaml
# entrapment.yaml — v3 新增 provenance 段
similarity:
  # ... v1+v2 字段不变 ...

provenance:
  fragment_tolerance_ppm: 20.0      # 碎片匹配容差 (ppm)
  max_fragment_charge: 2            # 理论碎片最大电荷态
  chimera_threshold: 0.3            # shared_ratio > 此值判定为 chimera
  min_peaks_for_analysis: 6         # 峰数少于此值跳过溯源
  levels_to_trace: ["L2", "L3", "L4"]  # 需要做溯源的 Level
```

所有字段均有 serde 默认值，旧配置文件无需修改。

---

## 9. 算法细节

### 9.1 碎片离子溯源算法

```text
fn trace_provenance(psm, spectrum, config) -> FragmentProvenance:
  1. 生成 trap 理论碎片:
     trap_ions = generate_b_ions(psm.peptide, charge=1..max_charge)
              ∪ generate_y_ions(psm.peptide, charge=1..max_charge)

  2. 生成 target 理论碎片 (如有 best_target_peptide):
     target_ions = generate_b_ions(best_target, charge=1..max_charge)
               ∪ generate_y_ions(best_target, charge=1..max_charge)

  3. 对谱图中每个实测峰:
     match_trap   = find_closest_within_tolerance(peak.mz, trap_ions, tolerance)
     match_target = find_closest_within_tolerance(peak.mz, target_ions, tolerance)

     provenance = match (match_trap, match_target):
       (Some, Some) → Shared       (记录两个离子标签)
       (Some, None) → TrapOnly     (记录 trap 离子标签)
       (None, Some) → TargetOnly   (记录 target 离子标签)
       (None, None) → Unassigned

  4. 统计:
     assigned = trap_only + target_only + shared
     shared_ratio = if assigned > 0 { shared / assigned } else { 0.0 }
     is_chimera = shared_ratio > config.chimera_threshold
```

### 9.2 L4 特殊处理

L4 PSM 没有 `best_target_peptide`（v2 未找到匹配），此时：
- `target_ions` 为空
- 所有匹配 trap 理论碎片的峰 → TrapOnly
- 所有不匹配的峰 → Unassigned
- `shared_ratio = 0.0`，`is_chimera = false`
- 仍然有价值：可以看到 trap 肽段对谱图的解释程度

### 9.3 修饰在碎片离子生成中的应用

trap 肽段如果有修饰信息，生成理论碎片时应用修饰质量：
- 修饰位置 i 之前的 b 离子不受影响
- 修饰位置 i 及之后的 b 离子 += mod_delta
- y 离子对称处理

使用现有 `search-engine` 的 `generate_b_ions` / `generate_y_ions` 函数，传入修饰作为 fixed_mods。

---

## 10. 测试策略

### 10.1 单元测试

- `provenance.rs`：已知谱图 + 已知肽段对 → 验证逐峰分类
- `mod_parser.rs`：各种 UniMod 格式解析（单修饰、多修饰、N-term 修饰）
- `mirror_plot.rs`：HTML 生成正确性（包含 Plotly 数据）
- `similarity.rs`：修饰感知 delta_mass 计算正确性

### 10.2 集成测试

- 端到端：parquet + FASTA + mzML → classified TSV (含 provenance 列) + HTML 报告 (含 mirror plot)
- 向后兼容：不提供 mzML 时，行为与 v2 完全一致

### 10.3 测试数据

- 使用现有 HeLa 混合物种数据（131K PSMs）的 mzML 子集
- 选取若干 L2/L3 PSM 手动验证碎片分类

---

## 11. 实施顺序建议

| 阶段 | 内容 | 依赖 |
|------|------|------|
| T1 | 修饰解析器 (`mod_parser.rs`) | 无 |
| T2 | UnifiedPsm 扩展 + DIA-NN loader 修饰解析 | T1 |
| T3 | 修饰感知 delta_mass（改造 similarity.rs） | T2 |
| T4 | 碎片离子溯源引擎 (`provenance.rs`) | 无（可与 T1-T3 并行） |
| T5 | 配置扩展 (`config.rs` provenance 段) | 无（可并行） |
| T6 | Mirror Plot 渲染 (`mirror_plot.rs`) | T4 |
| T7 | TSV 输出扩展 | T4 |
| T8 | HTML 报告增强（chimera 统计 + 嵌入 mirror plot） | T6, T7 |
| T9 | MCP Tools（trace_fragment_provenance + annotate_provenance） | T4, T6 |
| T10 | CLI 扩展（entrapment-cli 新增 mzml_dir 参数） | T4, T7 |
| T11 | 集成测试 + 回归测试 | 全部 |
