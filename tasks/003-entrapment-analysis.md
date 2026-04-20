# PRD: Entrapment Analysis — 陷阱库命中分类与同源性分析工具

> **文件名**：`003-entrapment-analysis.md`
> **版本**：1.0
> **创建日期**：2026-04-20
> **状态**：Design Approved — 待实施
> **前置**：无（独立工具，可选集成到 proteinCopilot）

---

## 1. 概述

### 1.1 背景与动机

在蛋白质组学中，entrapment（陷阱库）策略是评估搜索引擎 FDR 控制质量的常用方法：将非目标物种的蛋白数据库混入搜索库，命中陷阱库的 PSM 理论上应全为假阳性。然而实际分析发现，由于以下原因，大量陷阱库命中其实是"真实的人源肽段被错误归属"：

1. **Razor 归属错误**：跨物种保守蛋白（actin、EF1A、enolase 等）的肽段序列在人类和陷阱物种间完全一致，搜索引擎 razor 规则将其归给陷阱物种
2. **L/I 同分异构**：Leu(L) 和 Ile(I) 分子量完全相同（113.08406 Da），质谱无法区分
3. **近同源肽段**：仅 1-2 个氨基酸差异的保守家族肽段，在宽隔离窗 DIA 中可能被共洗脱信号"污染"

这些"不可分"或"难以分辨"的命中会导致 entrapment FDR 评估**系统性偏高**（过于悲观），干扰后续 ML 模型训练。

### 1.2 目标

构建一个**通用的陷阱库命中分类分析工具**，能够：

1. 对搜索结果中的 trap PSM 进行同源性分级（L0-L4）
2. 识别并报告 razor 归属错误
3. 产出干净的"真陷阱"负样本集（L3+L4）
4. 提供交互式 HTML 报告供人工审查
5. 独立运行（CLI），也可集成到 proteinCopilot（MCP tool）

### 1.3 非目标（v2 或更后）

- **不等长序列比对**：v1 的 L2/L3 仅比较等长肽段（Hamming distance），遗漏了不等长的同源物。v2 应升级为 edit distance（Levenshtein）或序列比对算法，使 L2/L3 能检测 indel 导致的不等长同源 → v2
- **Q/K 近等质量检测**：Q↔K 替换 Δm=36.4 mDa，MS1 在高电荷态（z≥3）时难以区分；需在 L1 中增加 near-isobaric 子类 → v2
- **等质量二肽替换**：GG↔N (114.043 Da)、AG↔Q (128.059 Da)、AD↔EG (186.064 Da) 等不等长但等质量的替换，MS2 理论可分但实际困难 → v2
- 共洗脱碎片溯源可视化（需 mzML，复杂度高）→ v2
- 修饰感知的序列比对 → v2
- ML 模型训练与特征提取 → 独立模块
- 定量层面分析（L/H ratio 提取等）→ 独立模块

---

## 2. 分级体系（Discriminability Level）

每条 trap PSM 被分为 5 个 Level：

| Level | 名称 | 定义 | 质谱可分性 | 处理建议 |
|-------|------|------|-----------|---------|
| **L0** | razor-error | stripped sequence 精确存在于 target 库 | 完全不可分 | 标记剔除 |
| **L1** | LI-isomer | 仅 L↔I 互换，其余位置完全相同 | 质谱不可分（Δm=0） | 标记剔除 |
| **L1.5** *(v2)* | near-isobaric-sub | Q↔K 单替换（Δm=36.4 mDa, z≥3 时 MS1 不可分）；或等质量二肽替换 GG↔N, AG↔Q, AD↔EG（不等长但等质量） | MS1 高电荷不可分；MS2 弱可分 | 重点观察 |
| **L2** | near-identical | 1 处非 L/I 差异 且 Δm < 阈值（默认 1.0 Da） | 弱可分 | 重点观察 |
| **L3** | homolog | 1-2 处差异 且 Δm ≥ 阈值，或 2-mismatch | 理论可分 | 保留进 ML |
| **L4** | true-trap | target 库中无 ≤2 mismatch 的对应肽段 | 理想负样本 | 保留进 ML |

### 2.1 判定算法

```text
fn classify(trap_pep, target_peptides, config):
    1. L0: exact string match in target digest
    2. L1: L/I-normalized match (replace L↔I then compare)
    3. [v1: 等长 only] Scan all target peptides of same length:
       - compute hamming distance + delta_mass
       - best = min by hamming
       [v2: 升级为 edit distance，扫描 len±2 范围的 target 肽段]
    4. if best.hamming == 1 && |delta_mass| < threshold → L2
    5. if best.hamming in {1,2} → L3
    6. else → L4
```

### 2.2 可配置参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `max_mismatches` | 2 | 最大允许的氨基酸差异数 |
| `delta_mz_threshold_da` | 1.0 | L2 vs L3 的 Δ mass 分界线 |
| `require_tryptic_ends` | true | 要求 target 对应肽两端也是胰酶切位点 |
| `max_missed_cleavages` | 2 | in-silico digest 允许的 missed cleavage |

---

## 3. 输入格式

### 3.1 搜索结果适配

通过 `ResultLoader` trait 支持多种搜索引擎输出：

| 引擎 | 格式 | v1 支持 | 备注 |
|------|------|---------|------|
| DIA-NN | Parquet | ✅ | 首要支持，可复用 result-import crate |
| 通用 | TSV | ✅ | 用户映射列名 |
| pFind | .spectra | Phase 7 | 按需 |
| MSFragger | pepXML/TSV | Phase 7 | 按需 |

### 3.2 target/trap 规则配置（YAML）

支持三种规则类型混合使用：

```yaml
# entrapment.yaml
version: 1

target:
  rules:
    - type: accession_contains
      any_of: ["_HUMAN"]
  fasta:
    - path: ./dbs/human_sp.fasta    # 用于同源性扫描
  # accession_list: ./target_ids.txt  # 可选白名单

trap:
  rules:
    - type: accession_contains
      any_of: ["_YEAST", "_ECOLI", "_YARLI", "_DICDI", "_YERPE"]

conflict_resolution: prefer_target   # target & trap 均匹配时
unmatched: ignore                    # 不匹配任一规则时: ignore | trap | target | error

similarity:
  max_mismatches: 2
  delta_mz_threshold_da: 1.0
  require_tryptic_ends: true
  max_missed_cleavages: 2
```

### 3.3 规则类型详述

| 类型 | 语法 | 说明 |
|------|------|------|
| `accession_contains` | `any_of: ["_HUMAN"]` | protein accession 包含任一子串 |
| `accession_regex` | `pattern: "^sp\\|[OPQ]"` | 正则匹配 |
| `fasta` | `path: ./db.fasta` | 从 FASTA 文件提取 accession 集合 |
| `accession_list` | `./ids.txt` | 显式白名单文件（一行一个） |

### 3.4 冲突解决

当一条 PSM 的 Protein.Ids 列同时匹配 target 和 trap 规则时：
- `prefer_target`：归为 target（推荐，保守策略）
- `prefer_trap`：归为 trap
- `mark_ambiguous`：单独标记为 ambiguous，不参与后续分级

---

## 4. 输出

### 4.1 Classified TSV / Parquet

每 PSM 一行，包含以下列：

| 列 | 类型 | 说明 |
|----|------|------|
| `peptide` | String | Stripped sequence |
| `charge` | i32 | 前体电荷 |
| `precursor_mz` | f64 | 前体 m/z |
| `retention_time` | f64 | 保留时间（min） |
| `scan_number` | u32 | 扫描号 |
| `spectrum_file` | String | 谱图文件名 |
| `protein_ids` | String | 蛋白 accession（分号分隔） |
| `group` | String | "target" / "trap" / "ambiguous" |
| `level` | String | "L0" / "L1" / "L2" / "L3" / "L4" |
| `best_target_peptide` | Option<String> | 最相似的 target 肽段序列 |
| `best_target_protein` | Option<String> | 对应的 target 蛋白 accession |
| `mismatches` | u8 | 差异数（L4 为 255 或 null） |
| `delta_mass_da` | f64 | 单同位素质量差（Da） |
| `diff_positions` | Option<String> | 差异位置描述 `[pos:X→Y,...]` |
| `q_value` | f64 | 原始搜索 q-value |

### 4.2 HTML 交互报告

使用 Plotly.js 生成交互式图表：

- **统计面板**：5 个 Level 的计数和占比
- **饼图**：Level 分布
- **蛋白家族 Bar Chart**：按保守蛋白家族聚类（actin、EF1A、enolase、HSP 等）
- **Δm 分布直方图**：L2/L3 的 delta mass 分布
- **可筛选 PSM 明细表**：按 Level / Protein / Δm 过滤排序，点击可展开详情

### 4.3 Razor Error List

专门的 TSV，列出所有 L0 PSM 的重新分配建议：

```
peptide       current_razor     suggested_razor     reason
STTTGHLIYK    EF1A_YARLI   →   EF1A1_HUMAN         exact_match
GYSFTTTAER    ACT1_DICDI   →   ACTB_HUMAN          exact_match
```

### 4.4 Run Metadata JSON

```json
{
  "tool_version": "0.1.0",
  "run_timestamp": "2026-04-20T12:00:00Z",
  "input_file": "hela_mix.parquet",
  "input_file_sha256": "abc123...",
  "config_snapshot": { ... },
  "target_fasta": "human_sp.fasta",
  "target_fasta_sha256": "def456...",
  "total_psms": 2421,
  "trap_psms": 2421,
  "level_counts": { "L0": 128, "L1": 45, "L2": 89, "L3": 312, "L4": 1847 }
}
```

---

## 5. 架构

### 5.1 Workspace Crate 结构

```text
crates/
├── entrapment-analysis/          ← 新增 lib crate (核心)
│   ├── src/
│   │   ├── lib.rs                EntrapmentAnalyzer, 公共 API
│   │   ├── config.rs             YAML 解析 (serde_yaml)
│   │   ├── loader/
│   │   │   ├── mod.rs            trait ResultLoader + UnifiedPsm
│   │   │   ├── diann_parquet.rs  DIA-NN parquet 加载
│   │   │   └── generic_tsv.rs    通用 TSV（用户映射列名）
│   │   ├── tagger.rs             target/trap 标记
│   │   ├── digest.rs             tryptic in-silico digest
│   │   ├── similarity.rs         brute-force L0-L4 分级
│   │   ├── report/
│   │   │   ├── mod.rs            HTML 生成 (Plotly.js)
│   │   │   └── charts.rs         图表数据生成
│   │   └── error.rs
│   ├── tests/
│   │   ├── fixtures/             mini parquet + fasta
│   │   └── integration.rs
│   └── Cargo.toml
│
├── entrapment-cli/               ← 新增 bin crate (薄壳)
│   ├── src/main.rs               clap: analyze / report / inspect
│   └── Cargo.toml
│
└── mcp-server/                   ← 扩展
    └── src/tools/entrapment.rs   3 个 MCP tool
```

### 5.2 核心 API

```rust
pub struct EntrapmentAnalyzer {
    config: EntrapmentConfig,
    target_digest: Vec<DigestedPeptide>,  // in-silico digest 结果
}

impl EntrapmentAnalyzer {
    pub fn new(config: EntrapmentConfig, fasta_path: &Path) -> Result<Self>;
    pub fn classify(&self, psm: &UnifiedPsm) -> ClassifiedPsm;
    pub fn classify_all(&self, psms: &[UnifiedPsm]) -> Vec<ClassifiedPsm>;
    pub fn summary(&self, classified: &[ClassifiedPsm]) -> EntrapmentSummary;
}
```

### 5.3 依赖关系

```text
entrapment-cli ──→ entrapment-analysis ──→ core (数据结构)
mcp-server ──→ entrapment-analysis
```

不依赖 spectrum-io、search-engine 等（v1 不需要 mzML）。

---

## 6. CLI 接口

```bash
# 完整分析流程
entrapment analyze \
    --results hela_mix.parquet \
    --format diann-parquet \
    --config entrapment.yaml \
    --target-fasta human_sp.fasta \
    --out ./output

# 从已分级结果重新生成报告
entrapment report \
    --classified ./output/psms_classified.parquet \
    --out ./output/report.html

# 查看单条肽段详情
entrapment inspect \
    --classified ./output/psms_classified.parquet \
    --peptide DGFLLDGFPR
```

---

## 7. MCP Tools

```text
classify_entrapment_hits(results_file, format?, config_file, target_fasta, output_dir?)
    → EntrapmentSummary { counts_by_level, output_files, top_protein_families }

analyze_entrapment_stats(classified_file)
    → DetailedStats { level_distribution, protein_family_clusters, delta_mass_histogram }

find_similar_targets(peptide, target_fasta, max_mismatches?)
    → Vec<SimilarityHit> { target_peptide, protein, mismatches, delta_mass_da, level }
```

---

## 8. 实施阶段

| Phase | 内容 | 交付物 |
|-------|------|--------|
| **P1 · 核心骨架** | 数据结构 + config 解析 + digest + similarity (L0/L1) + 单元测试 | lib crate 可编译 |
| **P2 · Loader** | DIA-NN parquet loader + tagger + generic TSV loader | 可加载真实数据 |
| **P3 · 分级完善** | L2/L3/L4 分级 + diff_positions + razor-error 导出 + 集成测试 | classify_all 跑通 |
| **P4 · CLI** | `analyze` / `inspect` 子命令 + Parquet/TSV 输出 | 独立可执行 |
| **P5 · HTML 报告** | 交互式 Plotly 报告 + 可筛选表格 + `report` 子命令 | 完整分析流程 |
| **P6 · MCP Tools** | 3 个 MCP tool 集成到 mcp-server | LLM 可调用 |
| **P7 · 扩展 Loader** | pFind / MSFragger 适配器（按需） | 多引擎支持 |
| **[v2]** | 不等长序列比对（edit distance）· Q/K 近等质量检测 · 等质量二肽替换 (GG↔N, AG↔Q, AD↔EG) · L1.5 子级实现 · 共洗脱溯源可视化 · 修饰感知比对 | — |

---

## 10. v2 扩展

### 10.1 不等长序列比对（L2/L3 升级）

**问题**：v1 的 L2/L3 仅使用 Hamming distance 比较等长肽段，遗漏了因 indel（插入/删除）导致长度不同但高度同源的 target 肽段。

**现状（v1）**：
- 数据结构：`by_length: HashMap<usize, Vec<TargetPeptide>>` —— 按长度分桶
- 查询：只搜索 `peptides_of_length(trap_peptide.len())`
- 不等长同源物直接落入 L4（误判为"真陷阱"）

**升级方案（v2）**：
- 将 Hamming distance 替换为 **Levenshtein edit distance**，支持插入/删除/替换操作
- 扩展查询范围：搜索 `len-2..=len+2` 的 target 肽段
- edit distance ≤ max_mismatches 的匹配进入 L2/L3 判定
- 需要计算对齐后的 Δm（考虑插入/删除残基的质量）
- 性能考虑：edit distance O(mn) 比 hamming O(n) 慢，可能需要预筛（如 k-mer 过滤）

### 10.2 Q/K 近等质量替换

**问题**：Q (Gln, 128.0586 Da) 和 K (Lys, 128.0950 Da) 差异仅 36.4 mDa。

| 电荷 z | Δ(m/z) | Orbitrap @m/z 1000 (R=60k) 可分性 |
|--------|--------|-----------------------------------|
| +2 | 18.2 mDa | ✅ 可分 |
| +3 | 12.1 mDa | ⚠️ 边界 |
| +4 | 9.1 mDa | ❌ 困难 |
| +5 | 7.3 mDa | ❌ 不可分 |

**实现方案**：
- 在 `classify_single()` 中增加 Q↔K 检查步骤（L0/L1 之后、hamming 之前）
- 条件：hamming==1 且差异位是 Q↔K 替换 → 归入 L1.5
- 可选：结合 PSM 的 charge 信息判断实际可分性

### 10.3 等质量二肽替换

**问题**：某些单氨基酸与二肽组合具有完全相同的残基质量：

| 单残基 | 二肽 | 共享质量 (Da) |
|--------|------|--------------|
| N (Asn) | GG | 114.04293 |
| Q (Gln) | AG / GA | 128.05858 |
| — | AD / DA = EG / GE | 186.06405 |

**影响**：
- 肽段长度不同（trap 和 target 长度差 1），v1 等长比较完全遗漏
- MS2 碎片系列不同（多/少一个切割位点），理论可分但实际困难
- 搜索引擎评分差异可能很小，razor 错误风险高

**实现方案**：
- 在 `TargetDigestIndex` 中增加不等长查询：对每个 trap 肽段，也搜索 len±1 的 target 肽段
- 新增 `isobaric_dipeptide_check()`：检查差异位是否属于已知等质量替换对
- 匹配则归入 L1.5（near-isobaric-sub）

---

## 9. 测试策略

### 9.1 单元测试
- `config.rs`：YAML 解析各种规则组合、冲突解决、边界情况
- `digest.rs`：tryptic digest 正确性（missed cleavage、N/C 端）
- `similarity.rs`：每种 Level 的判定逻辑（精确匹配、L/I 互换、hamming distance、Δm 阈值）
- `tagger.rs`：target/trap 标记规则匹配

### 9.2 集成测试
- 使用 mini 数据集（~10 条 PSM + ~100 蛋白 FASTA）端到端测试
- 验证 6 条已知肽段（STTTGHLIYK 等）的分级结果与手动分析一致

### 9.3 Fixtures
- `tests/fixtures/mini_report.parquet`：从真实 DIA-NN 结果截取的小子集
- `tests/fixtures/mini_human.fasta`：仅包含 EF1A1、ACTB、ACTG、KAD2、GUAA、ENOA 的 mini FASTA
