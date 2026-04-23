# Entrapment Analysis（陷阱库分析）

对搜索引擎报告中的 trap（陷阱库）命中进行系统分类，判断每个命中到底是
真正的误鉴定，还是因序列同源性导致的"假阳性"。

## 核心概念

**Entrapment Database**（陷阱库）是将已知不应存在于样本中的蛋白序列加入搜索数据库，
用来独立于 target-decoy 策略评估 FDR 的方法。然而，某些 trap 命中并非真正的误鉴定：

- 人源蛋白中存在与酵母等物种高度同源的保守肽段
- L/I（亮氨酸/异亮氨酸）同分异构体在质谱中完全不可区分
- 近等重氨基酸替换（如 D→N）的质量差异小于仪器分辨能力

### L0–L4 分级体系

| 级别 | 名称 | 含义 | MS 可分辨 |
|------|------|------|-----------|
| **L0** | Razor Error | 肽段序列完全相同，存在于 target 和 trap 蛋白中（razor 分配错误） | ❌ |
| **L1** | L/I Isomer | 仅 L↔I 替换，单同位素质量相同（113.084064 Da） | ❌ |
| **L2** | Near-isobaric | 氨基酸替换的质量差 < fragment tolerance（如 D→N，Δm=0.9840 Da） | ❌ (MS1) / ⚠️ (MS2)¹ |
| **L3** | Distinguishable Homolog | 有同源性但质量可区分（Hamming distance ≤ max_mismatches） | ✅ |
| **L4** | True Trap | 无显著同源性的真正陷阱命中 | ✅ |

**分类优先级**：L0 → L1 → L2/L3 → L4（依次检查，首次命中即返回）

## 安装与构建

```bash
# 构建 CLI
cargo build --release -p protein-copilot-entrapment-cli

# 二进制位于
target/release/protein-copilot-entrapment-cli
```

## 快速开始

### 1. 准备 YAML 配置文件

```yaml
version: 1

target:
  rules:
    - type: Fasta
      path: human_swissprot.fasta

trap:
  rules:
    - type: AccessionContains
      any_of: ["_YEAST", "_ECOLI", "_DICDI"]

conflict_resolution: prefer_target
unmatched: trap

similarity:
  max_mismatches: 2
  delta_mz_threshold_da: 1.0
  max_missed_cleavages: 2
```

### 2. 运行分析

```bash
entrapment analyze \
  --results search_report.parquet \
  --config entrapment.yaml \
  --target-fasta human_swissprot.fasta \
  --mzml-dir /path/to/mzml/files \
  --out output/entrapment
```

### 3. 查看输出

输出目录包含 4 个文件：

| 文件 | 内容 |
|------|------|
| `classified.tsv` | 所有 PSM 的分类结果（含 level, best_target, provenance 等） |
| `entrapment_report.html` | 交互式 HTML 报告（饼图 + 柱状图 + Δm 直方图 + 溯源统计 + 可筛选表格） |
| `razor_errors.tsv` | L0 级别的 razor 分配错误（仅 trap 组） |
| `run_metadata.json` | 运行元数据（输入文件 SHA256、配置快照、计数统计） |

### 4. 检查单个肽段

```bash
entrapment inspect \
  --peptide ELTALAPSTMK \
  --target-fasta human_swissprot.fasta
```

## 配置参考

### 分类规则（`target.rules` / `trap.rules`）

每个规则使用 `type` 字段指定类型：

| 类型 | 参数 | 说明 |
|------|------|------|
| `AccessionContains` | `any_of: [...]` | 蛋白 accession 包含任一子串 |
| `AccessionRegex` | `pattern: "..."` | 蛋白 accession 匹配正则表达式 |
| `Fasta` | `path: "..."` | 蛋白 accession 出现在指定 FASTA 文件中 |
| `AccessionList` | `path: "..."` | 蛋白 accession 出现在纯文本列表中（每行一个） |

### 冲突解决（`conflict_resolution`）

当蛋白同时匹配 target 和 trap 规则时：

- `prefer_target`（默认）：分类为 target
- `prefer_trap`：分类为 trap
- `mark_ambiguous`：标记为 ambiguous

### 未匹配策略（`unmatched`）

当蛋白既不匹配 target 也不匹配 trap 规则时：

- `ignore`（默认）：视为 target
- `trap`：视为 trap
- `target`：视为 target
- `error`：报错退出

### 相似性参数（`similarity`）

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `max_mismatches` | 2 | Hamming 距离上限（超过则直接归 L4） |
| `delta_mz_threshold_da` | 1.0 | L2 判定阈值：Δm < 此值则为 near-isobaric |
| `max_missed_cleavages` | 2 | FASTA 酶切时允许的漏切位点数 |

## 输入格式

### DIA-NN Parquet

自动识别以下列（大小写敏感）：

| 列名 | 必需 | 说明 |
|------|------|------|
| `Stripped.Sequence` | ✅ | 去修饰肽段序列 |
| `Protein.Ids` | ✅ | 蛋白 accession（`;` 分隔多蛋白） |
| `Modified.Sequence` | | 带 UniMod 修饰的序列（v3 修饰解析） |
| `Precursor.Charge` | | 电荷态 |
| `Precursor.Mz` | | 前体 m/z |
| `RT` | | 保留时间 apex（分钟） |
| `RT.Start` | | 洗脱窗口起点（分钟，v3 RT-based scan lookup） |
| `RT.Stop` | | 洗脱窗口终点（分钟，v3 RT-based scan lookup） |
| `Q.Value` | | q-value |
| `Run` | | 运行名称 |

### 通用 TSV

支持自定义列名，必须包含 `peptide` 和 `protein` 列。

## MCP Tools（3 个）

当通过 MCP Server 使用时，提供以下工具：

| Tool | 功能 |
|------|------|
| `classify_entrapment_hits` | 运行完整分类流程（含溯源），返回统计摘要 + HTML 报告 |
| `analyze_entrapment_stats` | 从已分类 TSV 生成统计分析 |
| `find_similar_targets` | 查找单个肽段在 target 库中的最相似序列 |
| `annotate_provenance` | 对单个 trap PSM 进行碎片离子溯源分析（v3） |

## 示例

参见 `examples/hela-mix-2da-entrapment.yaml` — HeLa + 混合物种陷阱库实验的完整配置。

## 典型输出

```
=== Entrapment Analysis Summary ===
Total PSMs:     131159
  Target:       130625
  Trap:         534
  Ambiguous:    0

Trap PSM breakdown by discriminability level:
  L0 (razor error):         0
  L1 (L/I isomer):          39
  L2 (near-isobaric):       8
  L3 (distinguishable):     92
  L4 (true trap):           395

Provenance traced for 62 PSMs
```

## 版本历程

### v2（✅ 已完成）
- **编辑距离**：Levenshtein 替代 Hamming，支持不等长肽段比较
- **Q/K 检测**：SubstitutionType::QKSubstitution 注释（Δm=36.4 mDa）
- **等重二肽替换**：IsobaricDipeptide 检测（GG↔N, AG↔Q, AD↔EG）
- **k-mer 倒排索引**：pigeonhole 预筛加速跨长搜索

### v3（✅ 已完成）
- **碎片离子溯源**：b/y 离子匹配分类（TrapOnly / TargetOnly / Shared / Unassigned）
- **嵌合谱检测**：shared_ratio > threshold 标记 chimeric
- **UniMod 修饰解析**：DIA-NN Modified.Sequence → 位置 + delta mass
- **RT-based scan lookup**：DIA-NN 无 scan_number 时通过 RT + precursor m/z 查找 MS2
- **镜像图可视化**：trap vs target 碎片离子对比 HTML
- **容错设计**：缺失 mzML 文件跳过而非中断

### v4: 多目标碎片溯源 (Multi-Target Fragment Provenance)

v4 对每个 L2/L3 trap PSM 自动查找所有共洗脱的 target 肽段（轻标 + 重标 SILAC），将每个观测碎片离子归属到具体的 target 来源。

**功能特性：**
- **共洗脱索引（CoElutionIndex）**：基于 RT 窗口交叉 + DIA 隔离窗口匹配，O(log N + k) 查询
- **轻重标搜索**：同时查找 light 和 SILAC heavy 形式的共洗脱 target
- **多目标碎片匹配**：每个观测峰可归属到多个 target 的理论碎片离子
- **DIA 双扫描镜像（Dual-Scan Mirror）**：轻标和重标前体落在不同隔离窗口，分别从对应的 MS2 扫描读取观测谱图
  - 轻标镜像：使用轻标前体 m/z 通过 `find_by_rt` 定位 MS2 扫描
  - 重标镜像：使用 SILAC 偏移后的重标前体 m/z 定位不同的 MS2 扫描
  - 重标镜像的 trap 理论离子使用 SILAC delta 偏移（`shift_ions_heavy`）
  - 验证 `heavy_scan != light_scan`，避免重复读取同一扫描
- **Per-PSM HTML 报告**：
  - 信息头：轻/重标前体 m/z、电荷、轻标扫描号、重标扫描号、谱图文件名
  - 候选表：共洗脱 target 肽段（含前体 m/z、蛋白 ID、标记形式）
  - 轻标镜像谱图（Light Mirror）：归一化强度、离子标注、噪声过滤
  - 重标镜像谱图（Heavy Mirror）：独立扫描号、SILAC 偏移 trap 离子
  - Trap 碎片归属表：仅展示 trap 碎片离子的 target 匹配

**数据结构：**
- `MirrorData`：每个镜像的独立数据（scan_number + annotated_peaks + 统计计数）
- `MultiTargetProvenance`：`light: MirrorData` + `heavy: Option<MirrorData>` + 候选列表

**前置条件：**
- `--mzml-dir` 参数指向 mzML 文件目录
- config 中配置 `provenance.silac` 块（可选，启用重标搜索）

**输出文件：**
- `provenance_summary.html` — 所有溯源 PSMs 的汇总表
- `provenance/` 目录 — 每个 PSM 一份独立的 HTML 报告

**新增配置参数（`provenance` 块下）：**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `silac.heavy_k_delta` | 8.014199 | 重标 Lysine delta mass (Da) |
| `silac.heavy_r_delta` | 10.008269 | 重标 Arginine delta mass (Da) |
| `silac.enable_heavy_search` | true | 是否搜索重标候选 |
| `generate_per_psm_reports` | true | 是否生成 per-PSM HTML 报告 |
| `max_co_eluting_candidates` | 20 | 每个 trap PSM 的最大候选数 |

**镜像谱图展示优化：**
- 噪声峰过滤：<5% 相对强度的峰以极淡色（`#f0f0f0`）背景展示
- 匹配峰标注：每个峰顶标注 b/y 离子名称（字号 10、颜色跟随离子）
- 柱宽统一 1.2 Da（匹配峰）/ 0.5 Da（未匹配）/ 0.3 Da（噪声）
- 无边框纯色柱形，视觉清晰
