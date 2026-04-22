# v4 Multi-Target Fragment Ion Provenance — Design Spec

> **日期**：2026-04-22
> **状态**：Approved
> **分支**：待创建 `feat/entrapment-v4-multi-target-provenance`
> **前置**：v3 provenance（已合并 master）

## 1. 问题陈述

v3 provenance 仅对比 **trap vs 1 个 best_target**，输出 shared_ratio + is_chimeric 布尔值。但在 DIA 宽窗口采集模式下，一个 MS2 谱图中可能包含来自多个共洗脱 target 肽段的碎片离子。

**v4 目标**：对每个 L2/L3 trap PSM，找出所有共洗脱 target（轻标 + 重标），将每个观测碎片离子归属到具体的 target 来源，生成 per-PSM 可视化报告。

**关键洞察**：轻标和重标的共洗脱差异本身就是判别证据——如果轻标 target 在同一 DIA 窗口但重标不在，则共享碎片只可能来自轻标形式。

## 2. 设计决策

| 决策 | 选择 | 理由 |
|------|------|------|
| 共洗脱 target 来源 | DIA-NN 搜索结果 | 已有 q-value 过滤，数量可控 |
| 共洗脱判定 | RT 洗脱窗口有交集 | 更全面，不遗漏边缘共洗脱 |
| DIA 窗口匹配 | 同一个 isolation window | 从 mzML 读取窗口边界，精确匹配 |
| SILAC 参数 | YAML config 配置 | K+8.014, R+10.008 可配置 |
| 报告生成范围 | 所有 L2/L3 trap PSMs | 约 100 个，全自动生成 |
| 架构方案 | 方案 A：索引 + 溯源 | CoElutionIndex 查询高效，结构清晰 |

## 3. 架构设计

### 3.1 处理管道

```
输入: L2/L3 trap PSM (e.g. STTTGHLIYK, RT=35.2 min)
  │
  ├─ Step 1: 查找共洗脱 targets
  │   ├─ RT 交叉检查: target.[RT.Start, RT.Stop] ∩ trap.[RT.Start, RT.Stop] ≠ ∅
  │   ├─ DIA 窗口过滤: target.precursor_mz 落在同一个 isolation window
  │   └─ 重标配对: heavy_mz = light_mz + (K×8.014 + R×10.008) / charge
  │
  ├─ Step 2: 碎片离子多路匹配
  │   ├─ 为 trap + 每个 candidate 生成理论 b/y 离子
  │   ├─ 对每个观测峰: 匹配所有候选的理论离子
  │   └─ 分类: TrapOnly / TargetOnly(指定来源) / Shared(列出所有匹配) / Unassigned
  │
  └─ Step 3: per-PSM HTML 报告
      ├─ Section 1: 共洗脱 Target 候选表
      ├─ Section 2: 镜像谱图 (Trap ↑ vs Targets ↓, 颜色区分来源)
      └─ Section 3: 碎片离子归属表 (每峰每个匹配候选)
```

### 3.2 核心数据结构

```rust
/// 共洗脱候选 target
struct CoElutingCandidate {
    peptide: String,
    protein_ids: Vec<String>,
    precursor_mz: f64,           // 轻标 precursor m/z
    charge: i32,
    rt_start: f64,               // 洗脱窗口起点 (min)
    rt_stop: f64,                // 洗脱窗口终点 (min)
    label_form: LabelForm,       // Light / Heavy
    modifications: Vec<(usize, f64)>,
}

/// 轻标/重标标识
enum LabelForm {
    Light,
    Heavy {
        precursor_mz_heavy: f64,
        residue_deltas: Vec<(usize, f64)>,  // 哪些残基加了多少 Da
    },
}

/// 多 target 溯源结果（per-PSM）
struct MultiTargetProvenance {
    trap_psm: UnifiedPsm,
    candidates: Vec<CoElutingCandidate>,
    annotated_peaks: Vec<MultiAnnotatedPeak>,
    scan_number: u32,
}

/// 观测峰的多路归属
struct MultiAnnotatedPeak {
    mz_observed: f64,
    intensity: f64,
    trap_ion: Option<String>,             // e.g. "b3+1"
    target_matches: Vec<TargetIonMatch>,   // 所有匹配到的 target 离子
}

/// 单个 target 离子匹配
struct TargetIonMatch {
    candidate_index: usize,       // 指向 candidates[i]
    ion_label: String,            // e.g. "y5+1"
    delta_ppm: f64,               // 匹配误差
}
```

### 3.3 共洗脱索引 (CoElutionIndex)

```rust
struct CoElutionIndex {
    /// 按 Run 分组 → 按 RT.Start 排序的 target PSMs
    by_run: HashMap<String, Vec<TargetEntry>>,
    /// DIA isolation windows (从 mzML 提取)
    dia_windows: Vec<(f64, f64, f64)>,   // (center, low, high)
    /// SILAC 配置
    silac: Option<SilacConfig>,
}

impl CoElutionIndex {
    /// 从 DIA-NN parquet 的 target PSMs 构建
    fn build(
        all_psms: &[UnifiedPsm],
        groups: &[PsmGroup],    // target/trap 标记
        dia_windows: &[(f64, f64, f64)],
        silac: Option<&SilacConfig>,
    ) -> Self;

    /// 查找与指定 trap PSM 共洗脱的所有 targets（轻 + 重）
    fn find_co_eluting(
        &self,
        trap: &UnifiedPsm,
        run: &str,
    ) -> Vec<CoElutingCandidate>;
}
```

**查询算法**：
1. 取 `by_run[run]`（已按 `RT.Start` 排序）
2. 二分查找 `RT.Start ≤ trap.RT.Stop` 的起始位置
3. 向后扫描直到 `target.RT.Start > trap.RT.Stop`
4. 过滤：`target.RT.Stop ≥ trap.RT.Start`（窗口交集）
5. 过滤：`target.precursor_mz` 落在同一个 DIA isolation window
6. 对每个轻标 target → 计算 `heavy_mz` → 检查 heavy 是否在同一/不同窗口 → 生成 Heavy 候选

**复杂度**：O(log N + k)，其中 k 为 RT 窗口内的 target 数

### 3.4 DIA Isolation Window 提取

从 mzML 的 `SpectrumIndex` 中提取 DIA isolation windows：
- 遍历所有 MS2 scan 的 `isolation_window: Option<(target_mz, lower, upper)>`
- 去重得到 DIA scheme（所有唯一窗口）
- 判断两个 precursor_mz 是否在同一窗口：`low ≤ mz ≤ high`

### 3.5 重标碎片离子生成

对每个 Heavy 候选：
- 获取轻标序列的理论 b/y 离子（复用现有 `generate_theoretical_ions()`）
- 对每个 b_n 离子：`heavy_mz = light_mz + cumulative_delta(prefix[0..n]) / charge`
- 对每个 y_n 离子：`heavy_mz = light_mz + cumulative_delta(suffix[len-n..]) / charge`
- `cumulative_delta` 根据 prefix/suffix 中的 K/R 残基数量计算

## 4. 可视化报告设计

### 4.1 Per-PSM HTML 报告

每个 L2/L3 trap PSM 生成一个独立 HTML 文件，包含 3 个 section：

**Section 1: 共洗脱 Target 候选表**
| 列 | 说明 |
|----|------|
| # | 颜色标识（与镜像图颜色对应） |
| Peptide | target 肽段序列 |
| Protein | 蛋白 accession |
| Label | 🔵 Light / 🔴 Heavy |
| m/z | precursor m/z |
| z | 电荷态 |
| RT range | 洗脱窗口 |
| DIA window | isolation window 范围 |
| Matched ions | 匹配的碎片离子数 / 总理论离子数 |

**Section 2: 镜像谱图**（Plotly.js 交互式）
- X 轴：m/z
- 上半（正向）：Trap 碎片离子
- 下半（负向）：所有 Target 碎片离子叠加
- 颜色编码：每个 candidate 一种颜色，Shared 用紫色
- Hover：显示 ion label、delta ppm、来源 target

**Section 3: 碎片离子归属表**
| 列 | 说明 |
|----|------|
| Obs. m/z | 观测值 |
| Intensity | 强度 |
| Trap Ion | trap 匹配的离子标签 |
| Origin | TrapOnly / TargetOnly / Shared / Shared×N / Unassigned |
| Target Matches | 颜色标识 + peptide + ion_label + (light/heavy) |
| Δppm | 匹配误差 |

**Footer**：Summary 统计（总峰数、各分类计数、结论）

### 4.2 汇总报告 (provenance_summary.html)

所有溯源 PSMs 的概览表，包含：
- Peptide、Level、Run、#Candidates、#Shared peaks、shared_ratio、is_chimeric
- 可筛选/排序
- 点击跳转到 per-PSM HTML

## 5. 配置扩展

### 5.1 YAML config 新增

```yaml
provenance:
  fragment_tolerance_ppm: 20.0
  max_fragment_charge: 2
  levels_to_trace: ["L2", "L3"]
  rt_tolerance_min: 0.5
  # v4 新增
  silac:
    heavy_k_delta: 8.014199      # ¹³C₆¹⁵N₂-Lys
    heavy_r_delta: 10.008269     # ¹³C₆¹⁵N₄-Arg
    enable_heavy_search: true
  generate_per_psm_reports: true
  max_co_eluting_candidates: 20  # 防止密集区域候选爆炸
```

### 5.2 Rust 类型

```rust
struct SilacConfig {
    heavy_k_delta: f64,          // default 8.014199
    heavy_r_delta: f64,          // default 10.008269
    enable_heavy_search: bool,   // default true
}

// ProvenanceConfig 扩展:
struct ProvenanceConfig {
    // ... 已有字段 ...
    silac: Option<SilacConfig>,
    generate_per_psm_reports: bool,      // default true
    max_co_eluting_candidates: usize,    // default 20
}
```

## 6. 输出结构

```
output/entrapment/
├── classified.tsv                 # v1/v2（不变）
├── razor_errors.tsv               # v1（不变）
├── run_metadata.json              # v1（不变）
├── entrapment_report.html         # v1/v2（不变）
├── provenance_summary.html        # v4 NEW: 所有溯源 PSMs 汇总
└── provenance/                    # v4 NEW: per-PSM HTML
    ├── STTTGHLIYK_z2_Rep1.html
    ├── AELTALAPSTMK_z2_Rep1.html
    └── ...
```

## 7. CLI & MCP 集成

### 7.1 CLI

已有的 `--mzml-dir` 参数即可触发 v4 溯源。当 config 中包含 `silac` 配置时自动启用多 target 溯源。

### 7.2 MCP Tool

新增单 PSM MCP tool：
```
trace_psm_provenance(
    peptide, run, results_file, mzml_dir, config_file
) → MultiTargetProvenance + HTML path
```

## 8. 文件变更范围

### 8.1 新增文件
| 文件 | 内容 |
|------|------|
| `crates/entrapment-analysis/src/coelution.rs` | CoElutionIndex 构建 + 查询 |
| `crates/entrapment-analysis/src/multi_provenance.rs` | 多 target 碎片匹配逻辑 |
| `crates/entrapment-analysis/src/multi_report.rs` | per-PSM HTML 报告渲染 |

### 8.2 修改文件
| 文件 | 变更 |
|------|------|
| `config.rs` | 新增 SilacConfig, ProvenanceConfig 扩展 |
| `types.rs` | 新增 CoElutingCandidate, MultiTargetProvenance 等类型 |
| `lib.rs` | 新增 trace_multi_target_provenance() 入口 |
| `output.rs` | provenance_summary.html 输出 |
| `entrapment-cli/src/main.rs` | 调用多 target 溯源 |
| `mcp-server/src/tools.rs` | 新增 trace_psm_provenance tool |

### 8.3 不变文件
- `provenance.rs`（v3 单 target 溯源保留，供 v4 复用 `generate_theoretical_ions()`）
- `mirror_plot.rs`（v3 简单镜像图保留）
- `similarity.rs`、`digest.rs`、`tagger.rs`、`levenshtein.rs`（v1/v2 核心逻辑不变）

## 9. 测试策略

### 9.1 单元测试
- `coelution.rs`：RT 交集判定、DIA 窗口匹配、重标 m/z 计算
- `multi_provenance.rs`：多路碎片匹配（Shared 跨多 target、TrapOnly、Unassigned）
- `multi_report.rs`：HTML 渲染正确性（候选表、归属表结构）

### 9.2 集成测试
- 使用 mini parquet + mini mzML 构建共洗脱场景
- 验证 per-PSM 报告生成完整性

### 9.3 真实数据验证
- HeLa mixed-species (hela-mix-2da) Rep1 550-600
- 对 L2/L3 trap PSMs 生成报告，人工检查碎片归属合理性

## 10. 复用 v3 基础设施

| v3 组件 | v4 复用方式 |
|---------|-----------|
| `generate_theoretical_ions()` | 为每个 candidate 生成理论离子 |
| `amino_acid_mass()` | 残基质量查表 |
| `parse_modified_sequence()` | UniMod 修饰解析 |
| `find_by_rt()` | RT-based scan lookup |
| `IndexedMzMLReader` | 读取 MS2 谱图 |
| `trace_provenance_batch()` 的文件分组逻辑 | 按 Run 分组处理 |
