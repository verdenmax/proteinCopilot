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
| **L2** | Near-isobaric | 氨基酸替换的质量差 < fragment tolerance（如 D→N，Δm=0.9840 Da） | ❌ (MS1) / ⚠️ (MS2) |
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
  --out output/entrapment
```

### 3. 查看输出

输出目录包含 4 个文件：

| 文件 | 内容 |
|------|------|
| `classified.tsv` | 所有 PSM 的分类结果（含 level, best_target 等） |
| `entrapment_report.html` | 交互式 HTML 报告（饼图 + 柱状图 + Δm 直方图 + 可筛选表格） |
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
| `Precursor.Charge` | | 电荷态 |
| `Precursor.Mz` | | 前体 m/z |
| `RT` | | 保留时间（分钟） |
| `Q.Value` | | q-value |
| `Run` | | 运行名称 |
| `Precursor.Id` | | scan 编号 |

### 通用 TSV

支持自定义列名，必须包含 `peptide` 和 `protein` 列。

## MCP Tools（3 个）

当通过 MCP Server 使用时，提供以下工具：

| Tool | 功能 |
|------|------|
| `classify_entrapment_hits` | 运行完整分类流程，返回统计摘要 + HTML 报告 |
| `analyze_entrapment_stats` | 从已分类 TSV 生成统计分析 |
| `find_similar_targets` | 查找单个肽段在 target 库中的最相似序列 |

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
  L2 (near-isobaric):       6
  L3 (distinguishable):     71
  L4 (true trap):           418
```

## v2 路线图

- **L1.5 级别**：Q↔K 近等重替换（Δm=36.4 mDa）+ 等重双肽替换（GG↔N, AG↔Q, AD↔EG）
- **编辑距离**：替代 Hamming 距离，支持不等长肽段比较
- **共洗脱碎片追踪**：利用 mzML 原始数据验证轻重标 XIC 一致性
- **修饰感知比对**：在比较时考虑 PTM 位置差异
