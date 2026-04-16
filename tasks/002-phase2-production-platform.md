# PRD: ProteinCopilot Phase 2 — 从原型到生产级蛋白质组学平台

> **文件名**：`002-phase2-production-platform.md`
> **版本**：1.0
> **创建日期**：2026-04-16
> **状态**：Planning
> **前置**：001-mvp-proteomics-search-platform.md (✅ 已完成)

---

## 1. Phase 2 概述

### 1.1 从 MVP 到生产平台

MVP 阶段验证了 **Rust MCP Server + LLM Agent** 的全流程可行性：19 个 MCP 工具、11 个 crate、510+ 测试、DDA/DIA 双模式支持。但 MVP 使用的 SimpleSearchEngine（b/y 离子匹配）仅适合功能验证，无法用于真实科研。

Phase 2 的目标是将 ProteinCopilot 从**功能原型**升级为**可用于真实科研的生产级平台**，具备与 MaxQuant、FragPipe、Spectronaut 等成熟工具可比的核心能力。

### 1.2 竞品对标

| 能力 | MaxQuant | FragPipe | Spectronaut | ProteinCopilot MVP | Phase 2 目标 |
|------|----------|----------|-------------|-------------------|-------------|
| 搜索引擎 | Andromeda | MSFragger | Pulsar | SimpleSearch ❌ | pFind ✅ |
| 蛋白推断 | ✅ | ✅ | ✅ | ❌ | ✅ Parsimony |
| 多文件批处理 | ✅ | ✅ | ✅ | 单文件 ❌ | ✅ 并行队列 |
| 定量分析 | ✅ LFQ/TMT/SILAC | ✅ | ✅ DIA-LFQ | ❌ | ✅ LFQ + TMT |
| 蛋白 FDR | ✅ | ✅ | ✅ | PSM 级 ⚠️ | ✅ 多级 FDR |
| MBR | ✅ | ✅ IonQuant | ✅ | ❌ | ⚠️ Phase 3 |
| PTM 定位 | ✅ | ✅ PTMProphet | ✅ | ❌ | ✅ 概率评分 |
| 标准格式导出 | ✅ | ✅ | ✅ | TSV/JSON ⚠️ | ✅ mzTab |
| QC | ❌ | ❌ | ✅ | ❌ | ✅ |
| AI 驱动 | ❌ | ❌ | ❌ | ✅ 🎯 | ✅ 强化 🎯 |
| DIA 支持 | 有限 | ✅ | ✅ (核心) | ✅ 基础 | ✅ 增强 |

### 1.3 Phase 2 设计原则

1. **搜索引擎为核心**：pFind 对接是最高优先级，没有生产级搜索引擎一切都是空谈
2. **批处理优先**：真实实验 = 几十到几百个文件，单文件模式无法满足需求
3. **分析深度递进**：PSM → 肽段 → 蛋白推断 → 定量，每一层依赖上一层
4. **多进程架构**：搜索引擎（pFind/MSFragger）天然是独立进程，Rust 做编排和结果聚合
5. **AI 差异化**：这是我们相对竞品的核心优势，每个新模块都要有 AI 解释能力

---

## 2. MVP 延迟项（必须在 Phase 2 完成）

这些是 MVP 文档中明确标记为延迟的功能：

| ID | 延迟项 | 原始优先级 | Phase 2 归属 |
|----|--------|-----------|-------------|
| FR-3.2~3.4 | pFind adapter 完整实现 | MVP P0 (未完成) | **M2.1** |
| FR-4.6 | 搜索结果比较工具 | MVP P2 | **M2.8** |
| NG-1 | 独立 FDR 模块 (已有 fdr crate) | Phase 2 | **M2.3** 增强 |
| NG-2 | 质量控制模块 | Phase 2 | **M2.6** |
| NG-3 | 蛋白推断 (Parsimony) | Phase 2 | **M2.3** |
| NG-4 | MSFragger / Comet adapter | Phase 2 | **M2.9** |
| NG-5 | 搜索失败诊断 Agent | Phase 2 | **M2.7** |
| pFind .spectra | pFind 结果解析器 | 阻塞于样本文件 | **M2.1** |

---

## 3. Phase 2 功能需求

### 优先级定义

| 层级 | 含义 | 标准 |
|------|------|------|
| **Tier 1** | 生产可用的最低要求 | 没有这些功能，平台无法用于真实科研 |
| **Tier 2** | 竞品对标的核心能力 | 缺少会被用户认为是"不完整"的平台 |
| **Tier 3** | 差异化与生态完善 | 提升用户体验和竞争力 |

---

### Tier 1：生产可用（必须完成）

#### M2.1：pFind 搜索引擎对接 ⭐ 最高优先级

**目标**：对接真实 pFind 搜索引擎，替代 SimpleSearchEngine 用于生产。

**背景**：pFind 是 ICR 自研的蛋白质搜索引擎，支持 open search、多引擎打分。当前仅有 stub adapter。

| ID | 需求 | 优先级 | 说明 |
|----|------|--------|------|
| M2.1.1 | pFind 二进制发现与健康检查 | P0 | 自动发现 pFind 安装路径，检查版本 |
| M2.1.2 | SearchParams → pFind cfg 文件转换 | P0 | 将标准化参数翻译为 pFind 配置文件格式 |
| M2.1.3 | pFind 进程管理（启动/监控/终止） | P0 | 使用 `tokio::process::Command` 管理 pFind 子进程 |
| M2.1.4 | pFind .spectra 结果解析器 | P0 | 解析 pFind 输出文件，转换为标准 SearchResult |
| M2.1.5 | pFind 日志实时流解析 | P1 | 解析 pFind 控制台输出，提取进度信息 |
| M2.1.6 | pFind 远程调用支持 (SSH) | P2 | 通过 SSH 在远程服务器/集群上运行 pFind |

**架构决策**：
- pFind 是独立进程（C++ 二进制），通过 `std::process::Command` 启动
- Rust 不直接链接 pFind，而是生成配置文件 → 启动进程 → 解析结果文件
- 进度监控通过解析 pFind 的 stdout/stderr 日志

**涉及 crate**：
- `search-engine` — `adapters/pfind.rs` 从 stub 升级为完整实现
- `result-import` — `pfind.rs` 实现 .spectra 解析

---

#### M2.2：多文件批处理与并行调度

**目标**：支持一次分析多个谱图文件，并行执行搜索，统一汇总结果。

**背景**：真实蛋白质组学实验通常产生 10-200 个 raw/mzML 文件。MVP 的单文件模式无法满足需求。

| ID | 需求 | 优先级 | 说明 |
|----|------|--------|------|
| M2.2.1 | 多文件输入支持 | P0 | `run_search` 接受文件列表，分配到多个搜索任务 |
| M2.2.2 | 并行搜索调度器 | P0 | 控制并发数（默认 CPU 核心数），队列管理 |
| M2.2.3 | 进度聚合 | P0 | 每个文件独立进度 + 总体进度百分比 |
| M2.2.4 | 结果合并 | P0 | 多文件 PSM 合并，统一 FDR 计算 |
| M2.2.5 | 单文件失败不阻塞 | P1 | 某个文件搜索失败时，其他文件继续，最终报告部分结果 |
| M2.2.6 | 批量任务持久化 | P1 | 批处理状态写入磁盘，支持中断恢复 |

**架构决策**：
- **进程级并行**：每个文件启动独立 pFind 进程（pFind 本身不支持多文件合并搜索）
- **Rust 做编排**：`tokio::task::spawn` 管理并发，`Semaphore` 控制最大进程数
- **不使用分布式**：Phase 2 仅支持单机多进程，分布式（多节点）留到 Phase 3

**新 crate**：`batch` (或集成到 `search-engine`)
- `BatchScheduler` — 任务队列、并发控制、进度聚合
- `BatchResult` — 多文件合并结果

**新 MCP 工具**：
- `run_batch_search` — 提交多文件批处理任务
- `get_batch_status` — 查询批处理总体进度

---

#### M2.3：蛋白推断与多级 FDR

**目标**：从 PSM 级别结果推断蛋白质组，实现蛋白级别 FDR 控制。

**背景**：MVP 的 FDR 仅在 PSM 级别。蛋白推断（Parsimony）是蛋白质组学分析的标准步骤，所有竞品都有此功能。

| ID | 需求 | 优先级 | 说明 |
|----|------|--------|------|
| M2.3.1 | Parsimony 蛋白推断算法 | P0 | 最小蛋白集合覆盖所有肽段 |
| M2.3.2 | 蛋白分组 (Protein Groups) | P0 | 同源蛋白归组，主蛋白 + 亚组成员 |
| M2.3.3 | Razor 肽段分配 | P0 | 共享肽段分配给最"丰富"的蛋白组 |
| M2.3.4 | 蛋白级 FDR | P0 | 基于 picked-protein approach 的蛋白 FDR |
| M2.3.5 | 肽段级 FDR | P0 | 独立于 PSM 的肽段级 q-value |
| M2.3.6 | 序列覆盖率计算 | P1 | 每个蛋白的氨基酸序列覆盖百分比 |
| M2.3.7 | 推断过程可视化 | P2 | 蛋白-肽段映射关系图 (HTML) |

**新 crate**：`protein-inference`
- `ParsimonyEngine` — 贪心集合覆盖
- `ProteinGroup` — 蛋白组数据结构
- `RazorAssigner` — 共享肽段分配
- `ProteinFdr` — picked-protein FDR

**新 MCP 工具**：
- `infer_proteins` — 执行蛋白推断
- `get_protein_groups` — 查看蛋白分组详情

---

#### M2.4：Rescoring / Percolator 集成

**目标**：使用半监督机器学习方法（Percolator 或 mokapot）重新评分 PSM，显著提高鉴定灵敏度。

**背景**：Percolator 是蛋白质组学的标准重评分工具，通常能提高 10-30% 的 PSM 鉴定数。MaxQuant、FragPipe 都集成了类似功能。

| ID | 需求 | 优先级 | 说明 |
|----|------|--------|------|
| M2.4.1 | PSM 特征提取 | P0 | 从搜索结果提取 Percolator 需要的特征向量 |
| M2.4.2 | Percolator/mokapot 进程调用 | P0 | 生成 PIN 文件 → 调用 Percolator → 解析输出 |
| M2.4.3 | 重评分结果回写 | P0 | 将 Percolator q-value 回写到 SearchResult |
| M2.4.4 | 自动判断是否需要重评分 | P1 | 基于初始搜索质量，AI 建议是否启用 Percolator |

**架构决策**：
- Percolator 是 Python/C++ 工具，通过进程调用
- 也可选 mokapot (纯 Python)，更易安装
- 特征提取在 Rust 中完成（确定性），仅调用外部工具做 ML 部分

**涉及 crate**：`fdr` (扩展) 或新 crate `rescoring`

---

### Tier 2：竞品对标（重要功能）

#### M2.5：定量分析

**目标**：从鉴定结果到蛋白丰度定量，支持 Label-free、TMT、SILAC 三种主流方式。

| ID | 需求 | 优先级 | 说明 |
|----|------|--------|------|
| M2.5.1 | Spectral Counting 定量 | P0 | 基于 PSM 计数的最简定量 |
| M2.5.2 | MS1 强度提取 (LFQ 基础) | P0 | 从 MS1 提取前体离子色谱峰面积 |
| M2.5.3 | TMT Reporter Ion 提取 | P1 | 从 MS2/MS3 提取 TMT 标记离子强度 |
| M2.5.4 | SILAC H/L 比值计算 | P1 | 基于 XIC 峰面积比计算 H/L ratio |
| M2.5.5 | 蛋白丰度汇总 | P1 | 从肽段定量到蛋白定量 (Top3/iBAQ/MaxLFQ) |
| M2.5.6 | 归一化 | P1 | Median normalization, Quantile normalization |
| M2.5.7 | 缺失值填充 | P2 | 基于正态分布的 imputation |

**新 crate**：`quantitation`
- `SpectralCounter` — PSM 计数定量
- `Ms1Quantifier` — MS1 XIC 峰面积
- `TmtQuantifier` — TMT reporter 离子
- `SilacQuantifier` — SILAC H/L ratio
- `ProteinQuantifier` — Top3, iBAQ, MaxLFQ 蛋白汇总

**新 MCP 工具**：
- `quantify` — 执行定量分析
- `normalize` — 数据归一化

**依赖**：M2.3 (蛋白推断) 必须先完成，因为定量需要蛋白分组

---

#### M2.6：质量控制模块

**目标**：自动化评估谱图和搜索质量，生成诊断报告。

| ID | 需求 | 优先级 | 说明 |
|----|------|--------|------|
| M2.6.1 | 谱图级 QC 指标 | P0 | TIC 分布、MS1/MS2 比例、峰数分布、空谱率 |
| M2.6.2 | 搜索级 QC 指标 | P0 | ID 率 vs 期望值、质量偏差分布、电荷分布异常 |
| M2.6.3 | 保留时间稳定性 | P1 | RT 偏差分析，批次效应检测 |
| M2.6.4 | 自动评级 | P1 | 将 QC 结果划分为 Good/Acceptable/Poor，AI 解释 |
| M2.6.5 | QC 报告 HTML | P2 | 交互式 HTML 报告，Plotly.js 图表 |

**新 crate**：`qc`
**新 MCP 工具**：`assess_quality`, `get_qc_report`

---

#### M2.7：PTM 定位评分

**目标**：为带修饰的 PSM 提供修饰位点的概率评分。

**背景**：磷酸化等修饰搜索中，修饰可能在多个 Ser/Thr/Tyr 位点之间。用户需要知道"修饰在哪个位点最可能"。

| ID | 需求 | 优先级 | 说明 |
|----|------|--------|------|
| M2.7.1 | 修饰位点概率评分 (类似 Ascore) | P0 | 基于碎片离子的位点区分度打分 |
| M2.7.2 | 位点定位分类 | P0 | Class I (>0.75), Class II (0.5-0.75), Class III (<0.5) |
| M2.7.3 | 多位点修饰处理 | P1 | 同一肽段多个修饰的联合定位 |
| M2.7.4 | 定位结果集成到导出 | P0 | PSM.tsv 增加 localization_probability 列 |

**涉及 crate**：`search-engine` (扩展 matching.rs) 或新 crate `ptm-scoring`

---

### Tier 3：差异化与生态完善

#### M2.8：搜索结果比较工具

**目标**：比较两次搜索的 PSM/肽段/蛋白 overlap，支持参数优化迭代。

| ID | 需求 | 优先级 | 说明 |
|----|------|--------|------|
| M2.8.1 | 两次搜索的 PSM overlap | P0 | 共有/独有 PSM 统计 |
| M2.8.2 | 肽段/蛋白级 overlap | P0 | Venn 图数据 |
| M2.8.3 | 参数差异对比 | P1 | 高亮两次搜索的参数差异 |
| M2.8.4 | HTML 可视化 | P2 | 交互式 Venn 图 + 差异表 |

**新 MCP 工具**：`compare_results`

---

#### M2.9：多搜索引擎支持

**目标**：扩展搜索引擎 adapter，支持 MSFragger 和 Comet。

| ID | 需求 | 优先级 | 说明 |
|----|------|--------|------|
| M2.9.1 | MSFragger adapter | P1 | SearchEngineAdapter 实现 |
| M2.9.2 | Comet adapter | P2 | SearchEngineAdapter 实现 |
| M2.9.3 | 多引擎结果合并 | P2 | iProphet 风格的多引擎评分整合 |

---

#### M2.10：标准格式导出

**目标**：支持蛋白质组学社区标准格式，方便与其他工具互操作和论文发表。

| ID | 需求 | 优先级 | 说明 |
|----|------|--------|------|
| M2.10.1 | mzTab 导出 | P1 | PSI mzTab 1.0 格式 (蛋白/肽段/PSM) |
| M2.10.2 | mzIdentML 导出 | P2 | PSI mzIdentML 1.2 格式 |
| M2.10.3 | MaxQuant 格式兼容 | P2 | proteinGroups.txt / evidence.txt 风格 |

**涉及 crate**：`report` (扩展)

---

#### M2.11：搜索失败诊断 Agent

**目标**：当搜索结果质量差时，AI 自动分析原因并建议修复方案。

| ID | 需求 | 优先级 | 说明 |
|----|------|--------|------|
| M2.11.1 | 失败模式分类 | P0 | 低 ID 率、质量偏差大、修饰异常、数据库不匹配 |
| M2.11.2 | 诊断决策树 | P0 | 结构化规则判断失败原因 |
| M2.11.3 | 自动建议参数调整 | P1 | 基于诊断结果推荐修改的参数 |
| M2.11.4 | 一键重搜 | P2 | 应用建议参数直接重新搜索 |

**涉及**：`.github/agents/failure-diagnosis.agent.md`, `.github/prompts/diagnosis.prompt.md`

---

#### M2.12：AI 增强功能

**目标**：强化 AI 作为核心差异化能力，覆盖分析全流程。

| ID | 需求 | 优先级 | 说明 |
|----|------|--------|------|
| M2.12.1 | 多轮分析对话 | P1 | 在一次会话中迭代搜索：调参 → 重搜 → 比较 |
| M2.12.2 | 实验设计建议 | P2 | 基于用户描述的实验目标，建议搜索策略 |
| M2.12.3 | 文献关联 | P2 | 搜索到的蛋白关联 UniProt/PubMed 信息 |
| M2.12.4 | 结果摘要 Markdown 报告 | P1 | 生成可直接用于论文方法部分的段落 |

---

## 4. 架构演进

### 4.1 多进程架构

```text
┌─────────────────────────────────────────────────┐
│              MCP Server (Rust)                    │
│  ┌───────────────────────────────────────────┐  │
│  │         Batch Scheduler (tokio)           │  │
│  │  ┌─────┐ ┌─────┐ ┌─────┐ ┌─────┐       │  │
│  │  │Job 1│ │Job 2│ │Job 3│ │Job N│ ...    │  │
│  │  └──┬──┘ └──┬──┘ └──┬──┘ └──┬──┘       │  │
│  │     │       │       │       │            │  │
│  │     ▼       ▼       ▼       ▼            │  │
│  │  ┌──────────────────────────────────┐    │  │
│  │  │     Semaphore (max_workers=N)    │    │  │
│  │  └──────────────────────────────────┘    │  │
│  └───────────────────────────────────────────┘  │
│     │       │       │       │                    │
│     ▼       ▼       ▼       ▼                    │
│  ┌─────┐ ┌─────┐ ┌─────┐ ┌─────┐               │
│  │pFind│ │pFind│ │pFind│ │pFind│  (OS 子进程)    │
│  └─────┘ └─────┘ └─────┘ └─────┘               │
└─────────────────────────────────────────────────┘
```

**决策**：
- 搜索引擎作为 **OS 子进程** 运行（pFind/MSFragger 是独立二进制）
- Rust MCP Server 使用 **tokio Semaphore** 控制并发数
- 每个文件一个搜索任务，互不干扰
- 失败任务不阻塞其他任务

### 4.2 新 Crate 规划

| Crate | 类型 | 依赖 | 说明 |
|-------|------|------|------|
| `protein-inference` | lib | core, fdr | Parsimony + 蛋白分组 + Razor 肽段 |
| `quantitation` | lib | core, xic, protein-inference | 定量分析 (LFQ/TMT/SILAC) |
| `qc` | lib | core, spectrum-io | 质量控制指标计算 |
| `ptm-scoring` | lib | core, search-engine | PTM 定位概率评分 |
| `batch` | lib | core, search-engine | 批处理调度器 |
| `rescoring` | lib | core, fdr | Percolator/mokapot 集成 |

### 4.3 MCP 工具扩展（19 → ~30 个）

| 新工具 | 所属模块 | 说明 |
|--------|----------|------|
| `run_batch_search` | batch | 多文件批处理搜索 |
| `get_batch_status` | batch | 批处理进度查询 |
| `infer_proteins` | protein-inference | 蛋白推断 |
| `get_protein_groups` | protein-inference | 蛋白分组详情 |
| `quantify` | quantitation | 定量分析 |
| `normalize` | quantitation | 数据归一化 |
| `assess_quality` | qc | 质量评估 |
| `get_qc_report` | qc | QC 报告 |
| `rescore` | rescoring | Percolator 重评分 |
| `compare_results` | report | 结果比较 |
| `localize_ptms` | ptm-scoring | PTM 定位评分 |

### 4.4 数据流（Phase 2 完整管线）

```text
谱图文件 (mzML/mgf × N)
    │
    ▼
┌──────────────┐
│ read_spectra │ × N (并行)
└──────┬───────┘
       │
       ▼
┌──────────────────┐     ┌─────────────────┐
│ recommend_params │ ──→ │ 👤 用户确认参数 │
└──────────────────┘     └────────┬────────┘
                                  │
                                  ▼
                       ┌──────────────────┐
                       │ run_batch_search │ (并行 pFind 进程)
                       └────────┬─────────┘
                                │
                                ▼
                       ┌──────────────────┐
                       │     rescore      │ (Percolator, 可选)
                       └────────┬─────────┘
                                │
                                ▼
                       ┌──────────────────┐
                       │ infer_proteins   │ (Parsimony + 蛋白 FDR)
                       └────────┬─────────┘
                                │
                        ┌───────┴───────┐
                        ▼               ▼
               ┌──────────────┐ ┌──────────────┐
               │   quantify   │ │ localize_ptms│
               └──────┬───────┘ └──────┬───────┘
                      │                │
                      ▼                ▼
               ┌──────────────────────────────┐
               │     generate_summary         │
               │     assess_quality           │
               │     export_results (mzTab)   │
               └──────────────────────────────┘
                                │
                                ▼
                       🤖 AI 结果解读 → 👤 用户
```

---

## 5. 非功能需求

### 5.1 性能

| 指标 | 目标 | 说明 |
|------|------|------|
| 多文件并行 | 支持 N 个 pFind 进程并行 (N = CPU 核数) | Semaphore 控制 |
| 批处理 100 文件 | < 2 小时 (取决于 pFind 速度) | 单文件约 1 分钟 |
| 蛋白推断 (10万 PSM) | < 30 秒 | Parsimony 算法优化 |
| 定量 (1000 蛋白 × 100 文件) | < 5 分钟 | XIC 提取为瓶颈 |
| 内存 (100 文件合并) | < 8 GB | 流式聚合，不全部加载 |

### 5.2 可靠性

| 指标 | 目标 |
|------|------|
| 单文件失败不影响批处理 | 错误隔离，继续其他文件 |
| 批处理可中断恢复 | 状态持久化到磁盘 |
| 搜索引擎崩溃处理 | 超时 + 重试 (最多 2 次) |
| 数据完整性 | 所有中间结果可审计 |

### 5.3 可扩展性

| 指标 | 目标 |
|------|------|
| 新搜索引擎 | 实现 `SearchEngineAdapter` trait 即可 |
| 新定量方法 | 实现 `Quantifier` trait |
| 新 QC 指标 | 注册到 QC 指标集合 |
| 自定义数据库 | 支持用户 FASTA + 污染物合并 |

---

## 6. 实施计划

### 6.1 依赖关系图

```text
M2.1 (pFind)  ←── M2.2 (批处理) ←── M2.5 (定量)
     │                  │                  │
     │                  ▼                  │
     │         M2.3 (蛋白推断) ──────────→ │
     │                  │                  │
     ▼                  ▼                  ▼
M2.4 (Rescoring)    M2.7 (PTM)     M2.6 (QC)
                                         │
                                         ▼
                                   M2.8-M2.12
                                 (Tier 3 功能)
```

### 6.2 建议执行顺序

| 阶段 | Milestone | 依赖 | 估计复杂度 |
|------|-----------|------|-----------|
| **Phase 2a** | M2.1 pFind 对接 | 需要 pFind 二进制 | 高 |
| **Phase 2a** | M2.3 蛋白推断 | 无 (可并行) | 中 |
| **Phase 2b** | M2.2 批处理 | M2.1 | 高 |
| **Phase 2b** | M2.4 Rescoring | M2.1 | 中 |
| **Phase 2c** | M2.5 定量 | M2.3 | 高 |
| **Phase 2c** | M2.6 QC | 无 | 中 |
| **Phase 2c** | M2.7 PTM 定位 | M2.1 | 中 |
| **Phase 2d** | M2.8~M2.12 | 前置完成 | 各项低-中 |

### 6.3 Phase 2a 先行任务（可立即开始）

即使没有 pFind 二进制，以下任务可以立即开始：

1. ✅ **M2.3 蛋白推断** — 纯算法，仅依赖 PSM 结果（SimpleSearchEngine 可产出）
2. ✅ **M2.6 QC 模块** — 谱图质量评估不依赖搜索引擎
3. ✅ **M2.10 mzTab 导出** — 格式转换，不依赖搜索引擎
4. ✅ **M2.8 结果比较** — 已有 SearchResult 结构
5. ✅ **M2.11 诊断 Agent** — AI prompt 定义

---

## 7. 测试策略

### 7.1 单元测试

每个新 crate 必须有完整单元测试：
- 蛋白推断：已知蛋白-肽段映射 → 验证 Parsimony 结果
- 定量：已知 XIC 峰面积 → 验证 LFQ/TMT 计算
- QC：已知异常数据 → 验证异常检测
- PTM 定位：已知碎片离子 → 验证 Ascore 计算

### 7.2 集成测试

- pFind 端到端：mzML → pFind 搜索 → 结果解析 → FDR → 蛋白推断 → 定量
- 批处理：3 个文件 → 并行搜索 → 合并结果
- 全管线：read_spectra → run_batch_search → rescore → infer_proteins → quantify → export_results

### 7.3 基准测试

- 与 MaxQuant/FragPipe 对比：相同数据、相同 FASTA、相同 FDR 阈值下的鉴定数
- 蛋白推断一致性：与 ProteinProphet 结果对比
- 定量精度：与 MaxQuant LFQ 强度相关性

---

## 8. Phase 3 预览（Phase 2 之后）

| 功能 | 说明 |
|------|------|
| Match-Between-Runs (MBR) | 跨运行 ID 转移，显著提高低丰度蛋白鉴定 |
| 差异表达分析 | 统计学显著性检验，火山图，GO 富集 |
| 分布式计算 | 多节点并行搜索 (Slurm/PBS 集群) |
| Web Dashboard | 实时监控 UI，交互式结果浏览 |
| Spectral Library | 从搜索结果构建谱图库，用于 library search |
| DIA-NN 风格 library-free search | 不依赖谱图库的 DIA 鉴定 |
| Glycoproteomics | 糖基化修饰搜索 |
| Cross-linking MS | 交联质谱分析 |

---

## 变更记录

| 日期 | 变更内容 | 原因 |
|------|----------|------|
| 2026-04-16 | Phase 2 初始规划 | MVP 完成，进入下一阶段 |

