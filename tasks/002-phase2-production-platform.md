# PRD: ProteinCopilot Phase 2 — 从原型到生产级蛋白质组学平台

> **文件名**：`002-phase2-production-platform.md`
> **版本**：1.0
> **创建日期**：2026-04-16
> **状态**：In Progress — M2.4 ✅ 蛋白推断完成
> **前置**：001-mvp-proteomics-search-platform.md (✅ 已完成)

---

## 1. Phase 2 概述

### 1.1 从 MVP 到生产平台

MVP 阶段验证了 **Rust MCP Server + LLM Agent** 的全流程可行性：19 个 MCP 工具、11 个 crate、510+ 测试、DDA/DIA 双模式支持。但 MVP 使用的 SimpleSearchEngine（b/y 离子匹配）仅适合功能验证，无法用于真实科研。

Phase 2 的目标是将 ProteinCopilot 从**功能原型**升级为**可用于真实科研的生产级平台**，具备与 MaxQuant、FragPipe、Spectronaut 等成熟工具可比的核心能力。

### 1.2 竞品对标

#### 1.2.1 传统 DDA 为主平台

| 能力 | MaxQuant | FragPipe | ProteinCopilot MVP | Phase 2 目标 |
|------|----------|----------|-------------------|-------------|
| 搜索引擎 | Andromeda | MSFragger | SimpleSearch ❌ | pFind ✅ |
| 蛋白推断 | ✅ Razor | ✅ ProteinProphet | ❌ | ✅ Parsimony |
| 多文件批处理 | ✅ | ✅ | 单文件 ❌ | ✅ 并行队列 |
| 定量分析 | ✅ LFQ/TMT/SILAC | ✅ IonQuant | ❌ | ✅ LFQ + TMT |
| 蛋白 FDR | ✅ | ✅ | PSM 级 ⚠️ | ✅ 多级 FDR |
| MBR | ✅ | ✅ IonQuant | ❌ | ⚠️ Phase 3 |
| PTM 定位 | ✅ | ✅ PTMProphet | ❌ | ✅ 概率评分 |
| 标准格式导出 | ✅ | ✅ | TSV/JSON ⚠️ | ✅ mzTab |
| QC | ❌ | ❌ | ❌ | ✅ |
| AI 驱动 | ❌ | ❌ | ✅ 🎯 | ✅ 强化 🎯 |

#### 1.2.2 DIA 专用平台

| 能力 | DIA-NN | AlphaDIA | Spectronaut | ProteinCopilot MVP | Phase 2 目标 |
|------|--------|----------|-------------|-------------------|-------------|
| Library-free DIA | ✅ (核心优势) | ✅ AlphaPeptDeep | ✅ directDIA | ❌ | ⚠️ Phase 3 |
| 深度学习预测 | ✅ RT + fragment | ✅ AlphaPeptDeep | ✅ | ❌ | ⚠️ Phase 3 |
| DIA 定量 | ✅ (极高精度) | ✅ | ✅ | ❌ | ✅ 基础 XIC |
| 处理速度 | ⚡ 极快 (C++) | ⚡ 快 (Python+Numba) | 快 | 快 (Rust) | ✅ Rust 优势 |
| 多文件批处理 | ✅ | ✅ | ✅ | 单文件 ❌ | ✅ |
| MBR | ✅ (跨 run 转移) | ✅ | ✅ | ❌ | ⚠️ Phase 3 |
| 蛋白推断 | ✅ | ✅ | ✅ | ❌ | ✅ |
| 结果导入 | — | — | — | ✅ Parquet 导入 | ✅ 增强 |
| AI 驱动 | ❌ | ❌ | ❌ | ✅ 🎯 | ✅ 🎯 |

> **关键洞察**：
> - **DIA-NN** (Demichev, 2020) 已成为 DIA 领域事实标准，library-free 模式通过深度学习预测 RT 和碎片离子强度，实现无需谱图库的 DIA 搜索。C++ 实现，单文件性能极强。
> - **AlphaDIA** (Mann lab, 2024) 是 AlphaPept 生态的 DIA 模块，基于 AlphaPeptDeep 深度学习框架，代表最前沿架构方向。Python + Numba 实现。
> - **Spectronaut** (Biognosys) 是商业 DIA 平台标杆，directDIA 模式 + 深度学习，但闭源且昂贵。
> - **ProteinCopilot 的 DIA 策略**：
>   1. MVP 已支持 DIA 谱图读取、前体提取、DIA-NN 结果导入（Parquet 格式）
>   2. Phase 2 通过 pFind 对接实现 pseudo-DDA 搜索（已有 `extract_dia_precursors` 工具）
>   3. Phase 3 考虑 library-free DIA（需深度学习 RT/fragment 预测模型）
>   4. **AI 编排是核心差异化**——DIA-NN/AlphaDIA 都没有自然语言交互和智能参数推荐

#### 1.2.3 新兴工具

| 工具 | 语言 | 特点 | 与 ProteinCopilot 的关系 |
|------|------|------|------------------------|
| **Sage** | Rust | 纯 Rust 搜索引擎，极快，开源 | 架构最近，可作为第三方引擎集成 |
| **MSBooster** | Java | 深度学习辅助重评分 | M2.4 重评分可参考 |
| **AlphaPept** | Python | Mann lab 全流程平台 | 生态竞品，AI 方面我们有优势 |
| **ionbot** | Python | 开放搜索 + 深度学习 | M2.9 多引擎候选 |

### 1.3 Phase 2 设计原则

1. **搜索引擎编排为核心**：不重新发明搜索算法，集成最好的开源引擎（Sage/DIA-NN/pFind）
2. **Sage 优先**：Rust 原生，一个集成 = DDA+DIA+LFQ+TMT+FDR 五个能力（见竞品分析）
3. **批处理优先**：真实实验 = 几十到几百个文件，单文件模式无法满足需求
4. **分析深度递进**：PSM → 肽段 → 蛋白推断 → 定量，每一层依赖上一层
5. **混合架构**：Sage 用 library 集成，pFind/DIA-NN 用子进程编排，统一 `SearchEngineAdapter` trait
6. **AI 差异化**：这是我们相对竞品的核心优势，每个新模块都要有 AI 解释能力

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

---

#### M2.0：Sage 搜索引擎集成 ⭐⭐ 最高优先级

**目标**：集成 Sage (Rust 原生搜索引擎)，一个集成获得 DDA+DIA+LFQ+TMT+FDR 五大能力。

**背景**：Sage (lazear/sage) 是纯 Rust 开源搜索引擎，支持 fragment indexing、wide-window DIA、LFQ/TMT 定量、LDA+KDE FDR、chimeric 谱匹配。可作为 Rust library crate 直接集成，无需子进程。

**新 crate**：无（集成到 `search-engine` crate）

**详细任务**：

| ID | 任务 | 优先级 | 涉及文件 | 验收标准 | 依赖 |
|----|------|--------|----------|----------|------|
| M2.0.1 | **添加 sage-core crate 依赖** | P0 | `crates/search-engine/Cargo.toml` | `cargo build` 通过，sage-core 可导入 | 无 |
| M2.0.2 | **SearchParams → SageConfig 转换器** | P0 | `search-engine/src/adapters/sage.rs` (新) | 所有字段正确映射：enzyme、mods、tolerance、FASTA path、DIA mode | M2.0.1 |
| M2.0.2a | → enzyme 枚举映射 | P0 | 同上 | Trypsin/LysC/GluC/AspN/Chymotrypsin/TrypsinP/NonSpecific 全映射 | — |
| M2.0.2b | → modification 映射 | P0 | 同上 | fixed_mods/variable_mods → Sage static_mods/variable_mods JSON 格式 | — |
| M2.0.2c | → tolerance 单位转换 | P0 | 同上 | ppm/Da → Sage [lower, upper] 格式 | — |
| M2.0.2d | → DIA 模式配置 | P0 | 同上 | `wide_window: true`, `chimera: true` 当 acquisition_mode=DIA | — |
| M2.0.3 | **SageAdapter 实现 SearchEngineAdapter trait** | P0 | `search-engine/src/adapters/sage.rs` | `search()` 调用 sage-core API，返回标准 `SearchResult` | M2.0.2 |
| M2.0.3a | → Sage 结果 → PSM 转换 | P0 | 同上 | scan_number, peptide, charge, mz, score, q_value 全映射 | — |
| M2.0.3b | → Sage 结果 → Protein 转换 | P0 | 同上 | protein accession, coverage, peptide_count 映射 | — |
| M2.0.3c | → Sage decoy 标记映射 | P0 | 同上 | Sage `is_decoy` → 我们的 `Psm.is_decoy` | — |
| M2.0.3d | → 错误处理包装 | P0 | 同上 | Sage 错误 → `CoreError::SearchEngineFailed` | — |
| M2.0.4 | **注册 SageAdapter 到 EngineRegistry** | P0 | `search-engine/src/registry.rs` | `create_adapter("sage")` 返回 SageAdapter 实例 | M2.0.3 |
| M2.0.5 | **health_check 实现** | P0 | `search-engine/src/adapters/sage.rs` | 验证 sage-core 版本，返回 `HealthStatus::Healthy` | M2.0.1 |
| M2.0.6 | **engine_info 返回 Sage 元数据** | P0 | 同上 | name="sage", version 从 sage-core crate 版本获取 | M2.0.1 |
| M2.0.7 | **MCP 工具适配：run_search 支持 Sage** | P0 | `mcp-server/src/tools.rs` | `run_search` 自动检测/允许指定 engine="sage" | M2.0.4 |
| M2.0.8 | **param-recommend 适配** | P1 | `param-recommend/src/*.rs` | 当检测到 Sage 可用时，推荐 Sage 作为引擎（DDA/DIA 均可） | M2.0.4 |
| M2.0.9 | **LFQ 定量包装** | P1 | `search-engine/src/adapters/sage.rs` | 提取 Sage LFQ 结果，映射到标准定量数据结构 | M2.0.3 |
| M2.0.10 | **TMT 定量包装** | P1 | 同上 | 提取 Sage TMT reporter ion 结果 | M2.0.3 |
| M2.0.11 | **DIA chimeric 模式测试** | P1 | `integration-tests/` | DIA mzML → Sage wide_window 搜索 → 验证多肽匹配 | M2.0.3 |
| M2.0.12 | **单元测试：config 转换** | P0 | `search-engine/src/adapters/sage.rs` | ≥10 个测试用例覆盖所有参数组合 | M2.0.2 |
| M2.0.13 | **集成测试：端到端搜索** | P0 | `integration-tests/` | mzML + FASTA → Sage 搜索 → PSM/Protein 结果验证 | M2.0.7 |

---

#### M2.1：pFind 搜索引擎对接

**目标**：对接真实 pFind 搜索引擎，替代 SimpleSearchEngine 用于生产。

**背景**：pFind 是 ICR 自研的蛋白质搜索引擎，支持 open search、多引擎打分。当前仅有 stub adapter。

**涉及 crate**：`search-engine` (adapters/pfind.rs), `result-import` (pfind.rs)

**详细任务**：

| ID | 任务 | 优先级 | 涉及文件 | 验收标准 | 依赖 |
|----|------|--------|----------|----------|------|
| M2.1.1 | **pFind 二进制路径发现** | P0 | `search-engine/src/adapters/pfind.rs` | 按优先级搜索：config → 环境变量 PFIND_HOME → PATH → 默认路径 | 无 |
| M2.1.1a | → 配置文件路径读取 | P0 | 同上 | 从 `.proteincopilot/config.toml` 读取 `pfind.binary_path` | — |
| M2.1.1b | → 环境变量发现 | P0 | 同上 | `$PFIND_HOME/bin/pFind3` 或 `$PFIND_PATH` | — |
| M2.1.1c | → PATH 搜索 | P0 | 同上 | `which pFind3` 等价逻辑 | — |
| M2.1.2 | **pFind 版本检测与健康检查** | P0 | 同上 | 执行 `pFind3 --version`，解析版本号，返回 HealthStatus | M2.1.1 |
| M2.1.3 | **SearchParams → pFind .cfg 文件生成** | P0 | `search-engine/src/adapters/pfind_config.rs` (新) | 生成完整 pFind 配置文件 | 无 |
| M2.1.3a | → [basic] section 生成 | P0 | 同上 | enzyme, missed_cleavages, min/max peptide length | — |
| M2.1.3b | → [modification] section 生成 | P0 | 同上 | fixed_mods/variable_mods → pFind mod 格式 (name@residue) | — |
| M2.1.3c | → [tolerance] section 生成 | P0 | 同上 | precursor_tolerance, fragment_tolerance (ppm/Da) | — |
| M2.1.3d | → [database] section 生成 | P0 | 同上 | FASTA path, decoy strategy (reverse/shuffle) | — |
| M2.1.3e | → [spectrum] section 生成 | P0 | 同上 | input file paths, spectrum format | — |
| M2.1.3f | → [output] section 生成 | P0 | 同上 | output directory, output format | — |
| M2.1.4 | **pFind 子进程启动与管理** | P0 | `search-engine/src/adapters/pfind.rs` | tokio::process::Command 启动 pFind，监控 exit code | M2.1.1 |
| M2.1.4a | → 进程启动 | P0 | 同上 | spawn pFind with cfg file path arg | — |
| M2.1.4b | → stdout/stderr 日志捕获 | P0 | 同上 | 实时转发到 tracing::info/error | — |
| M2.1.4c | → 超时处理 | P0 | 同上 | 可配置超时（默认 30min），超时 kill 进程 | — |
| M2.1.4d | → 退出码处理 | P0 | 同上 | exit 0 = 成功, 非 0 = 返回 SearchEngineFailed error | — |
| M2.1.5 | **pFind 日志进度解析** | P1 | `search-engine/src/adapters/pfind_progress.rs` (新) | 从 pFind stdout 解析 "Searching... X%" 进度信息 | M2.1.4 |
| M2.1.5a | → 日志格式正则匹配 | P1 | 同上 | 解析 pFind 3.x 的进度行格式 | — |
| M2.1.5b | → 进度回调通道 | P1 | 同上 | `tokio::sync::watch` 通道推送进度到 `get_search_status` | — |
| M2.1.6 | **pFind .spectra 结果解析器** | P0 | `result-import/src/pfind.rs` | 解析 pFind .spectra 输出 → 标准 `Vec<Psm>` | 无 |
| M2.1.6a | → 文件格式解析 | P0 | 同上 | 解析 .spectra 列格式：scan, peptide, charge, score, etc. | — |
| M2.1.6b | → 修饰名称映射 | P0 | 同上 | pFind 修饰名 → 标准 Modification 结构 | — |
| M2.1.6c | → decoy 标记识别 | P0 | 同上 | 识别 REV_ / DECOY_ 前缀蛋白 | — |
| M2.1.7 | **PFindAdapter 完整实现** | P0 | `search-engine/src/adapters/pfind.rs` | search() 完整流程：生成 cfg → 启动进程 → 等待 → 解析结果 | M2.1.3-M2.1.6 |
| M2.1.8 | **单元测试：cfg 生成** | P0 | `search-engine/src/adapters/pfind_config.rs` | ≥8 个用例，覆盖所有参数组合 | M2.1.3 |
| M2.1.9 | **集成测试：pFind 端到端** | P0 | `integration-tests/` | 需要 pFind 二进制可用，标注 `#[ignore]` 供有环境时执行 | M2.1.7 |
| M2.1.10 | **pFind 远程调用 (SSH)** | P2 | `search-engine/src/adapters/pfind_remote.rs` (新) | SSH 连接 → 上传 cfg → 远程执行 → 下载结果 | M2.1.7 |

---

#### M2.2：DIA-NN 搜索引擎对接

**目标**：编排 DIA-NN CLI，实现一键 DIA 分析。在已有 Parquet 结果导入基础上，增加进程管理。

**背景**：DIA-NN 是 DIA 领域事实标准。MVP 已实现 `import_search_results(diann_parquet)` 导入 DIA-NN 结果。Phase 2 将其升级为全自动编排：参数生成 → 进程管理 → 结果导入。

**涉及 crate**：`search-engine` (新 adapter), `result-import` (已有 diann.rs)

**详细任务**：

| ID | 任务 | 优先级 | 涉及文件 | 验收标准 | 依赖 |
|----|------|--------|----------|----------|------|
| M2.2.1 | **DIA-NN 二进制路径发现** | P0 | `search-engine/src/adapters/diann.rs` (新) | config → $DIANN_PATH → PATH → 默认路径 | 无 |
| M2.2.2 | **DIA-NN 版本检测** | P0 | 同上 | 执行 `diann --version` 或解析 banner 输出 | M2.2.1 |
| M2.2.3 | **SearchParams → DIA-NN CLI 参数转换** | P0 | `search-engine/src/adapters/diann_config.rs` (新) | 完整 CLI 参数映射 | 无 |
| M2.2.3a | → 基础参数映射 | P0 | 同上 | `--fasta`, `--threads`, `--qvalue`, `--out` | — |
| M2.2.3b | → 消化酶映射 | P0 | 同上 | `--cut K*,R*` (Trypsin) 等 → DIA-NN `--cut` 参数格式 | — |
| M2.2.3c | → 修饰映射 | P0 | 同上 | `--var-mod UniMod:35,15.994915,M` 格式 | — |
| M2.2.3d | → tolerance 映射 | P0 | 同上 | `--mass-acc`, `--mass-acc-ms1` (ppm) | — |
| M2.2.3e | → DIA 特有参数 | P0 | 同上 | `--gen-spec-lib`, `--predictor`, `--min-pep-len`, `--max-pep-len` | — |
| M2.2.3f | → Library-based 模式 | P1 | 同上 | `--lib library.tsv` 当提供谱图库时 | — |
| M2.2.3g | → plexDIA 参数 | P2 | 同上 | `--channel` 参数用于 SILAC/mTRAQ DIA | — |
| M2.2.4 | **DIA-NN 子进程管理** | P0 | `search-engine/src/adapters/diann.rs` | spawn → monitor → collect | M2.2.1 |
| M2.2.4a | → 进程启动 | P0 | 同上 | 构建完整命令行，spawn 子进程 | — |
| M2.2.4b | → 日志捕获 | P0 | 同上 | DIA-NN stdout 实时转发到 tracing | — |
| M2.2.4c | → 进度解析 | P1 | 同上 | 解析 `[INFO] Running file X of N` 等进度信息 | — |
| M2.2.4d | → 超时与 kill | P0 | 同上 | 可配置超时（默认 60min DIA 搜索较慢） | — |
| M2.2.5 | **DIA-NN 结果收集** | P0 | 同上 | 搜索完成后调用已有 `import_diann_parquet()` 导入 report.parquet | M2.2.4 |
| M2.2.5a | → 结果文件定位 | P0 | 同上 | 根据 `--out` 参数定位 report.tsv / report.parquet | — |
| M2.2.5b | → 定量矩阵收集 | P1 | 同上 | 收集 pg_matrix.tsv / pr_matrix.tsv 到输出目录 | — |
| M2.2.5c | → 谱图库保存 | P1 | 同上 | 保存 report-lib.tsv 供后续 library-based 搜索复用 | — |
| M2.2.6 | **DiannAdapter 实现 SearchEngineAdapter** | P0 | 同上 | 完整 search() 流程 | M2.2.3-M2.2.5 |
| M2.2.7 | **注册到 EngineRegistry** | P0 | `search-engine/src/registry.rs` | `create_adapter("diann")` 返回 DiannAdapter | M2.2.6 |
| M2.2.8 | **AI 自动引擎选择增强** | P1 | `param-recommend/src/*.rs` | 检测到 DIA 数据时推荐 DIA-NN；DDA 推荐 Sage/pFind | M2.2.7 |
| M2.2.9 | **单元测试：CLI 参数生成** | P0 | `search-engine/src/adapters/diann_config.rs` | ≥8 个测试用例 | M2.2.3 |
| M2.2.10 | **集成测试** | P0 | `integration-tests/` | 需要 DIA-NN 二进制，标注 `#[ignore]` | M2.2.6 |

---

#### M2.3：多文件批处理与并行调度

**目标**：支持一次分析多个谱图文件，并行执行搜索，统一汇总结果。

**背景**：真实蛋白质组学实验通常产生 10-200 个 raw/mzML 文件。MVP 的单文件模式无法满足需求。

**新 crate**：`batch` 或集成到 `search-engine`

**详细任务**：

| ID | 任务 | 优先级 | 涉及文件 | 验收标准 | 依赖 |
|----|------|--------|----------|----------|------|
| M2.3.1 | **BatchConfig 数据结构** | P0 | `core/src/batch.rs` (新) | max_workers, input_files, engine, params, output_dir | 无 |
| M2.3.2 | **BatchScheduler 核心调度器** | P0 | `search-engine/src/batch.rs` (新) | tokio::Semaphore 控制并发，per-file 任务分发 | M2.3.1 |
| M2.3.2a | → 任务队列 | P0 | 同上 | Vec<FileTask> → 按顺序/优先级分发到 worker | — |
| M2.3.2b | → Semaphore 并发控制 | P0 | 同上 | max_workers 默认 = CPU 核心数，可配置 | — |
| M2.3.2c | → 单文件搜索任务封装 | P0 | 同上 | `async fn run_single_file(adapter, params, file) → Result<SearchResult>` | — |
| M2.3.3 | **进度聚合** | P0 | 同上 | 每文件独立进度 + 总体 completed/total 百分比 | M2.3.2 |
| M2.3.3a | → BatchProgress 数据结构 | P0 | `core/src/batch.rs` | `{ total, completed, failed, running, per_file: Vec<FileProgress> }` | — |
| M2.3.3b | → watch 通道推送进度 | P0 | `search-engine/src/batch.rs` | `tokio::sync::watch::Sender<BatchProgress>` | — |
| M2.3.4 | **结果合并** | P0 | `search-engine/src/batch.rs` | 多文件 PSM 合并 → 统一 FDR 重新计算 | M2.3.2 |
| M2.3.4a | → PSM 列表合并 | P0 | 同上 | 保留 file_name 元数据，合并到单一 Vec<Psm> | — |
| M2.3.4b | → 全局 FDR 重新计算 | P0 | 同上 | 对合并后 PSM 重新计算 q-value（使用 fdr crate） | — |
| M2.3.4c | → Peptide/Protein 重新聚合 | P0 | 同上 | 跨文件 unique peptide 统计，protein coverage 更新 | — |
| M2.3.5 | **单文件失败不阻塞** | P1 | 同上 | 某文件失败 → 记录 error → 继续其他文件 → 最终报告 partial result | M2.3.2 |
| M2.3.5a | → 错误收集器 | P1 | 同上 | `BatchResult { results: Vec<FileResult>, errors: Vec<(PathBuf, Error)> }` | — |
| M2.3.5b | → partial success 策略 | P1 | 同上 | 当 ≥1 文件成功时返回 Ok(BatchResult)，全部失败时 Err | — |
| M2.3.6 | **批处理状态持久化** | P1 | `search-engine/src/batch.rs` | 批处理状态写入 `.proteincopilot/batch/{run_id}.json` | M2.3.2 |
| M2.3.6a | → 状态文件序列化 | P1 | 同上 | 每个文件完成时写入/更新状态文件 | — |
| M2.3.6b | → 中断恢复 | P2 | 同上 | 从状态文件恢复，跳过已完成文件 | M2.3.6a |
| M2.3.7 | **MCP 工具：run_batch_search** | P0 | `mcp-server/src/tools.rs` | 新 MCP tool，接受 file list + params | M2.3.2 |
| M2.3.8 | **MCP 工具：get_batch_status** | P0 | `mcp-server/src/tools.rs` | 查询批处理进度 | M2.3.3 |
| M2.3.9 | **单元测试** | P0 | `search-engine/src/batch.rs` | Semaphore 并发测试，结果合并测试 | M2.3.2 |
| M2.3.10 | **集成测试：3 文件并行** | P0 | `integration-tests/` | 3 个小 mzML → 并行搜索 → 合并结果 → FDR | M2.3.7 |

---

#### M2.4：蛋白推断与多级 FDR ✅ 已完成

**目标**：从 PSM 级别结果推断蛋白质组，实现蛋白级别 FDR 控制。

**状态**：✅ 全部完成（2026-04-16）

**新 crate**：`protein-inference`（mapper + parsimony + razor + coverage）

**已完成任务**：

| ID | 任务 | 状态 | 测试数 |
|----|------|------|--------|
| M2.4.1 | ProteinGroup + InferenceResult 数据结构 | ✅ | 验证测试 |
| M2.4.2 | 肽段-蛋白映射（I/L 等价 + decoy 分离 + q-value 过滤） | ✅ | 9 |
| M2.4.3 | Parsimony 算法（贪心集合覆盖 + 子集合并 + 不可区分分组） | ✅ | 11 |
| M2.4.4 | Razor 肽段分配（按 unique 数 → 分数 → 字母序） | ✅ | 10 |
| M2.4.5 | 肽段级 FDR（best PSM score → TDC q-value） | ✅ | 8 |
| M2.4.6 | 蛋白级 FDR（picked-protein 配对竞争） | ✅ | 10 |
| M2.4.7 | 序列覆盖率（FASTA 定位 + I/L 归一化 + 多出现点） | ✅ | 12 |
| M2.4.8 | MCP 工具 `infer_proteins` | ✅ | — |
| M2.4.10-12 | 单元 + 集成测试（14 个端到端场景） | ✅ | 14 |
| | **合计** | | **74+ 新测试** |

**未实现（降级到后续）**：
- M2.4.9 `get_protein_groups` 查看详情工具（可用 `infer_proteins` 结果替代）
- M2.4.13 推断过程可视化 HTML（P2 优先级）

---

#### M2.5：Rescoring / Percolator 集成

**目标**：使用半监督机器学习方法重新评分 PSM，提高 10-30% 鉴定灵敏度。

**背景**：Percolator 是标准重评分工具。Sage 也内置了 LDA+KDE rescoring。

**涉及 crate**：`fdr` (扩展) 或新 crate `rescoring`

**详细任务**：

| ID | 任务 | 优先级 | 涉及文件 | 验收标准 | 依赖 |
|----|------|--------|----------|----------|------|
| M2.5.1 | **PSM 特征提取器** | P0 | `rescoring/src/features.rs` (新) | 从 PSM 提取 Percolator 特征向量 | 无 |
| M2.5.1a | → 基础特征 | P0 | 同上 | score, delta_score, mass_error_ppm, charge, peptide_length | — |
| M2.5.1b | → 碎片匹配特征 | P0 | 同上 | matched_ions, matched_intensity_ratio, longest_b_series, longest_y_series | — |
| M2.5.1c | → 修饰特征 | P0 | 同上 | n_variable_mods, is_oxidized, etc. | — |
| M2.5.1d | → 消化特征 | P0 | 同上 | missed_cleavages, n_enzymatic_termini | — |
| M2.5.2 | **PIN 文件生成器** | P0 | `rescoring/src/pin.rs` (新) | PSM list → Percolator PIN format (TSV with header) | M2.5.1 |
| M2.5.3 | **Percolator 进程调用** | P0 | `rescoring/src/percolator.rs` (新) | spawn percolator → wait → parse output POUT file | M2.5.2 |
| M2.5.3a | → 二进制发现 | P0 | 同上 | config → $PERCOLATOR_PATH → PATH | — |
| M2.5.3b | → 进程管理 | P0 | 同上 | `percolator -X pout.xml pin.tsv` 调用 | — |
| M2.5.3c | → POUT 结果解析 | P0 | 同上 | 解析 Percolator XML/TSV output → 更新 q-value | — |
| M2.5.4 | **mokapot 替代方案** | P1 | `rescoring/src/mokapot.rs` (新) | Python mokapot 调用作为 Percolator 替代 | M2.5.2 |
| M2.5.5 | **重评分结果回写** | P0 | `rescoring/src/lib.rs` | Percolator q-value → 更新 SearchResult.psms[i].q_value | M2.5.3 |
| M2.5.6 | **AI 自动判断是否需要重评分** | P1 | `param-recommend/src/*.rs` | 初始 ID 率 < 30% 时建议启用 Percolator | M2.5.5 |
| M2.5.7 | **MCP 工具：rescore** | P0 | `mcp-server/src/tools.rs` | 新 MCP tool，输入 run_id → 重评分后的 SearchResult | M2.5.5 |
| M2.5.8 | **单元测试** | P0 | `rescoring/src/features.rs` | 特征提取正确性 | M2.5.1 |
| M2.5.9 | **集成测试** | P0 | `integration-tests/` | 需要 Percolator 二进制，标注 `#[ignore]` | M2.5.3 |

---

### Tier 2：竞品对标（重要功能）

---

#### M2.6：定量分析

**目标**：从鉴定结果到蛋白丰度定量，支持 Label-free、TMT、SILAC。

**新 crate**：`quantitation`

**详细任务**：

| ID | 任务 | 优先级 | 涉及文件 | 验收标准 | 依赖 |
|----|------|--------|----------|----------|------|
| M2.6.1 | **Quantifier trait 定义** | P0 | `core/src/quantitation.rs` (新) | `trait Quantifier { fn quantify(&self, result, spectra) → QuantResult }` | 无 |
| M2.6.1a | → QuantResult struct | P0 | 同上 | `{ peptide_quant: Vec<PeptideQuant>, protein_quant: Vec<ProteinQuant> }` | — |
| M2.6.1b | → PeptideQuant struct | P0 | 同上 | sequence, intensity, ratio (optional), file_name | — |
| M2.6.1c | → ProteinQuant struct | P0 | 同上 | accession, abundance, method, samples | — |
| M2.6.2 | **Spectral Counting** | P0 | `quantitation/src/spectral_count.rs` (新) | 每蛋白 PSM 计数 → 归一化频率 | 无 |
| M2.6.2a | → raw count | P0 | 同上 | protein_accession → psm_count | — |
| M2.6.2b | → NSAF 归一化 | P1 | 同上 | Normalized Spectral Abundance Factor | — |
| M2.6.3 | **MS1 强度提取 (LFQ)** | P0 | `quantitation/src/ms1_lfq.rs` (新) | 从 MS1 谱图提取 precursor XIC 面积 | 无 |
| M2.6.3a | → precursor XIC 提取 | P0 | 同上 | 利用已有 xic crate，提取 MS1 级色谱峰 | — |
| M2.6.3b | → 色谱峰面积积分 | P0 | 同上 | 梯形法 / 高斯拟合积分 | — |
| M2.6.3c | → 同位素包络求和 | P1 | 同上 | monoisotopic + M+1 + M+2 求和 | — |
| M2.6.4 | **TMT Reporter Ion 提取** | P1 | `quantitation/src/tmt.rs` (新) | 从 MS2/MS3 提取 TMT 标记离子强度 | 无 |
| M2.6.4a | → reporter ion m/z 定义 | P1 | 同上 | TMT-6plex, TMT-10plex, TMT-16plex, TMTpro 通道定义 | — |
| M2.6.4b | → MS2 reporter 提取 | P1 | 同上 | 在 reporter m/z ± tolerance 内找峰取强度 | — |
| M2.6.4c | → MS3 reporter 提取 | P2 | 同上 | 支持 SPS-MS3 模式 | — |
| M2.6.4d | → isotope purity correction | P2 | 同上 | 基于 TMT 批次校准矩阵 | — |
| M2.6.5 | **SILAC H/L 比值** | P1 | `quantitation/src/silac.rs` (新) | 基于已有 XIC heavy/light 对比（xic crate） | 无 |
| M2.6.5a | → heavy/light XIC 配对 | P1 | 同上 | 利用 xic crate 的 heavy label 功能 | — |
| M2.6.5b | → H/L ratio 计算 | P1 | 同上 | area_heavy / area_light | — |
| M2.6.5c | → ratio 归一化 | P1 | 同上 | median ratio normalization | — |
| M2.6.6 | **蛋白丰度汇总** | P1 | `quantitation/src/protein_quant.rs` (新) | 从肽段定量到蛋白定量 | M2.6.3 或 M2.6.4 |
| M2.6.6a | → Top3 法 | P1 | 同上 | 取前 3 个最强肽段的平均强度 | — |
| M2.6.6b | → iBAQ 法 | P1 | 同上 | 总强度 / 理论可观测肽段数 | — |
| M2.6.6c | → MaxLFQ 法 | P2 | 同上 | 最小化 ratio 差异的方法 (复杂) | — |
| M2.6.7 | **归一化** | P1 | `quantitation/src/normalize.rs` (新) | 数据归一化 | M2.6.6 |
| M2.6.7a | → Median normalization | P1 | 同上 | 每个样本除以中位数 | — |
| M2.6.7b | → Quantile normalization | P2 | 同上 | 分位数归一化 | — |
| M2.6.8 | **缺失值填充** | P2 | `quantitation/src/impute.rs` (新) | 基于正态分布低端采样 | M2.6.6 |
| M2.6.9 | **MCP 工具：quantify** | P0 | `mcp-server/src/tools.rs` | 新 MCP tool | M2.6.2-M2.6.6 |
| M2.6.10 | **MCP 工具：normalize** | P1 | `mcp-server/src/tools.rs` | 新 MCP tool | M2.6.7 |
| M2.6.11 | **单元测试** | P0 | 各子模块 | 每个定量方法 ≥5 个用例 | M2.6.2-M2.6.7 |

---

#### M2.7：质量控制模块

**目标**：自动化评估谱图和搜索质量，生成诊断报告。

**新 crate**：`qc`

**详细任务**：

| ID | 任务 | 优先级 | 涉及文件 | 验收标准 | 依赖 |
|----|------|--------|----------|----------|------|
| M2.7.1 | **QcMetric trait** | P0 | `qc/src/lib.rs` (新) | `trait QcMetric { fn name(); fn compute(data) → QcValue; fn assess() → QcGrade }` | 无 |
| M2.7.1a | → QcValue enum | P0 | 同上 | Numeric(f64), Distribution(Vec<f64>), Categorical(String) | — |
| M2.7.1b | → QcGrade enum | P0 | 同上 | Good, Acceptable, Poor, Critical | — |
| M2.7.2 | **谱图级 QC 指标** | P0 | `qc/src/spectrum_qc.rs` (新) | 基于 SpectrumSummary 计算 | 无 |
| M2.7.2a | → TIC 分布分析 | P0 | 同上 | TIC 随 RT 的变化曲线，检测异常下降 | — |
| M2.7.2b | → MS1/MS2 比例 | P0 | 同上 | 期望 MS2 >> MS1，异常则 DDA trigger 可能有问题 | — |
| M2.7.2c | → 空谱率 | P0 | 同上 | 无峰谱图占比，> 10% 告警 | — |
| M2.7.2d | → 峰数分布 | P0 | 同上 | median peaks per spectrum, 过低说明信号差 | — |
| M2.7.2e | → 前体质量分布 | P0 | 同上 | m/z 分布是否在预期范围 | — |
| M2.7.3 | **搜索级 QC 指标** | P0 | `qc/src/search_qc.rs` (新) | 基于 SearchResult 计算 | 无 |
| M2.7.3a | → ID 率 vs 期望值 | P0 | 同上 | 实际 ID 率 vs 样本类型期望值 (HeLa ~40-50%, plasma ~15-25%) | — |
| M2.7.3b | → 质量偏差分布 | P0 | 同上 | delta_mass_ppm 的 mean/std/偏斜度 | — |
| M2.7.3c | → 电荷分布异常检测 | P0 | 同上 | 期望 2+ 主导，5+/6+ 占比高说明可能有问题 | — |
| M2.7.3d | → 修饰频率异常 | P0 | 同上 | Oxidation(M) > 20% 提示样品制备问题 | — |
| M2.7.3e | → decoy 比例检查 | P0 | 同上 | target:decoy 比应接近 1:0 在高分段 | — |
| M2.7.4 | **RT 稳定性分析** | P1 | `qc/src/rt_qc.rs` (新) | 保留时间偏差，批次效应 | 无 |
| M2.7.4a | → RT 偏差计算 | P1 | 同上 | 同一肽段跨文件的 RT 标准差 | — |
| M2.7.4b | → RT 漂移检测 | P1 | 同上 | RT 随文件序号的系统性偏移 | — |
| M2.7.5 | **自动评级 (QcReport)** | P1 | `qc/src/report.rs` (新) | 汇总所有指标 → 总体 Good/Acceptable/Poor | M2.7.2-M2.7.4 |
| M2.7.5a | → 评级规则引擎 | P1 | 同上 | 各指标权重 + 阈值 → 综合评分 | — |
| M2.7.5b | → AI 解释文本生成 | P1 | 同上 | 为每个 Poor 指标生成人类可读解释 | — |
| M2.7.6 | **QC HTML 报告** | P2 | `report/src/qc_visualize.rs` (新) | Plotly.js 交互图表 | M2.7.5 |
| M2.7.6a | → TIC 图 | P2 | 同上 | TIC vs RT 曲线 | — |
| M2.7.6b | → 质量偏差直方图 | P2 | 同上 | delta_ppm distribution | — |
| M2.7.6c | → 电荷分布饼图 | P2 | 同上 | charge state distribution | — |
| M2.7.7 | **MCP 工具：assess_quality** | P0 | `mcp-server/src/tools.rs` | 新 MCP tool | M2.7.5 |
| M2.7.8 | **MCP 工具：get_qc_report** | P1 | `mcp-server/src/tools.rs` | 返回 HTML 报告路径 | M2.7.6 |
| M2.7.9 | **单元测试** | P0 | 各子模块 | 已知好/差数据的正确评级 | M2.7.2-M2.7.5 |

---

#### M2.8：PTM 定位评分

**目标**：为带修饰的 PSM 提供修饰位点的概率评分（类似 Ascore）。

**涉及 crate**：新 crate `ptm-scoring` 或 `search-engine` 扩展

**详细任务**：

| ID | 任务 | 优先级 | 涉及文件 | 验收标准 | 依赖 |
|----|------|--------|----------|----------|------|
| M2.8.1 | **位点候选枚举** | P0 | `ptm-scoring/src/candidates.rs` (新) | 给定肽段 + 修饰类型 → 枚举所有可能位点 | 无 |
| M2.8.1a | → 单修饰位点枚举 | P0 | 同上 | PEPTMS + Phospho → P1, P4, P5 (S/T/Y sites) | — |
| M2.8.1b | → 多修饰联合枚举 | P1 | 同上 | 2 个 Phospho → C(n,2) 种组合 | — |
| M2.8.2 | **Ascore 算法实现** | P0 | `ptm-scoring/src/ascore.rs` (新) | 基于碎片离子的位点区分度打分 | M2.8.1 |
| M2.8.2a | → site-determining ions 识别 | P0 | 同上 | 找到能区分不同位点分配的 b/y 离子 | — |
| M2.8.2b | → 离子匹配对比 | P0 | 同上 | 对每个候选位点，计算 matched / total site-determining ions | — |
| M2.8.2c | → Ascore 概率计算 | P0 | 同上 | -10 * log10(P) 其中 P 是随机匹配概率 | — |
| M2.8.3 | **概率转换与分类** | P0 | `ptm-scoring/src/lib.rs` | Ascore → 概率 → Class I/II/III | M2.8.2 |
| M2.8.3a | → Class I | P0 | 同上 | probability > 0.75 (localized) | — |
| M2.8.3b | → Class II | P0 | 同上 | 0.5 < probability ≤ 0.75 (ambiguous) | — |
| M2.8.3c | → Class III | P0 | 同上 | probability ≤ 0.5 (not localized) | — |
| M2.8.4 | **结果集成到 PSM 导出** | P0 | `report/src/lib.rs` | PSM.tsv 增加 `localization_probability`, `localization_class` 列 | M2.8.3 |
| M2.8.5 | **MCP 工具：localize_ptms** | P0 | `mcp-server/src/tools.rs` | 新 MCP tool | M2.8.3 |
| M2.8.6 | **单元测试** | P0 | `ptm-scoring/src/ascore.rs` | 已知谱图 → 已知位点概率 | M2.8.2 |
| M2.8.7 | **多位点修饰测试** | P1 | 同上 | 2-phospho peptide 正确评分 | M2.8.1b |

---

### Tier 3：差异化与生态完善

---

#### M2.9：搜索结果比较工具

**目标**：比较两次搜索的 PSM/肽段/蛋白 overlap，支持参数优化迭代。

**详细任务**：

| ID | 任务 | 优先级 | 涉及文件 | 验收标准 | 依赖 |
|----|------|--------|----------|----------|------|
| M2.9.1 | **PSM 级 overlap 计算** | P0 | `report/src/compare.rs` (新) | 两个 SearchResult → shared / A-only / B-only PSM counts | 无 |
| M2.9.1a | → PSM 匹配键定义 | P0 | 同上 | (scan_number, peptide_sequence, charge) 作为唯一键 | — |
| M2.9.1b | → overlap 统计 | P0 | 同上 | shared, a_only, b_only, jaccard_index | — |
| M2.9.2 | **肽段级 overlap** | P0 | 同上 | unique peptide sequences → Venn data | M2.9.1 |
| M2.9.3 | **蛋白级 overlap** | P0 | 同上 | protein accessions → Venn data | M2.9.1 |
| M2.9.4 | **参数差异对比** | P1 | 同上 | 两个 SearchParams → diff (哪些参数不同) | 无 |
| M2.9.5 | **Score 分布对比** | P1 | 同上 | 两次搜索的 score histogram 数据 | M2.9.1 |
| M2.9.6 | **HTML Venn 图 + 差异表** | P2 | `report/src/compare_visualize.rs` (新) | Plotly.js Venn + 表格 | M2.9.1-M2.9.5 |
| M2.9.7 | **MCP 工具：compare_results** | P0 | `mcp-server/src/tools.rs` | 新 MCP tool | M2.9.1-M2.9.3 |
| M2.9.8 | **单元测试** | P0 | `report/src/compare.rs` | 已知 overlap 数据验证 | M2.9.1 |

---

#### M2.10：多搜索引擎支持（扩展）

**目标**：在 Sage/pFind/DIA-NN 基础上，扩展更多引擎和多引擎结果合并。

**详细任务**：

| ID | 任务 | 优先级 | 涉及文件 | 验收标准 | 依赖 |
|----|------|--------|----------|----------|------|
| M2.10.1 | **MSFragger adapter** | P1 | `search-engine/src/adapters/msfragger.rs` (新) | SearchParams → .params 文件, 子进程管理, pepXML 解析 | M2.0 或 M2.1 (参考模式) |
| M2.10.1a | → MSFragger .params 生成 | P1 | 同上 | SearchParams 字段映射 | — |
| M2.10.1b | → 进程管理 | P1 | 同上 | java -jar MSFragger.jar 调用 | — |
| M2.10.1c | → pepXML 结果解析 | P1 | 同上 | pepXML → Vec<Psm> | — |
| M2.10.2 | **Comet adapter** | P2 | `search-engine/src/adapters/comet.rs` (新) | 类似 MSFragger 模式 | M2.10.1 |
| M2.10.3 | **多引擎结果合并** | P2 | `search-engine/src/multi_engine.rs` (新) | iProphet 风格多引擎 score 整合 | M2.0+M2.1 |
| M2.10.3a | → 结果对齐 | P2 | 同上 | 按 (scan, peptide, charge) 对齐不同引擎的 PSM | — |
| M2.10.3b | → 联合评分 | P2 | 同上 | 多引擎 score → 联合 probability | — |
| M2.10.4 | **单元测试** | P1 | 各 adapter | 配置生成正确性 | 各 adapter |

---

#### M2.11：标准格式导出

**目标**：支持蛋白质组学社区标准格式。

**详细任务**：

| ID | 任务 | 优先级 | 涉及文件 | 验收标准 | 依赖 |
|----|------|--------|----------|----------|------|
| M2.11.1 | **mzTab 1.0 导出** | P1 | `report/src/mztab.rs` (新) | PSI mzTab 格式：MTD + PRT + PEP + PSM sections | M2.4 (蛋白推断) |
| M2.11.1a | → MTD (metadata) section | P1 | 同上 | search engine, database, modifications, tolerance | — |
| M2.11.1b | → PRT (protein) section | P1 | 同上 | accession, description, coverage, abundance | — |
| M2.11.1c | → PEP (peptide) section | P1 | 同上 | sequence, modifications, charge, retention time | — |
| M2.11.1d | → PSM section | P1 | 同上 | spectrum ref, peptide, score, q-value | — |
| M2.11.2 | **mzIdentML 1.2 导出** | P2 | `report/src/mzidentml.rs` (新) | XML 格式，符合 PSI 规范 | M2.4 |
| M2.11.3 | **MaxQuant 兼容格式** | P2 | `report/src/maxquant_compat.rs` (新) | proteinGroups.txt + evidence.txt 风格 | M2.4 + M2.6 |
| M2.11.4 | **MCP 工具更新** | P1 | `mcp-server/src/tools.rs` | export_results 增加 format 参数: tsv/mztab/mzidentml | M2.11.1 |
| M2.11.5 | **mzTab 格式验证测试** | P1 | `report/tests/` | 输出通过 jmzTab validator | M2.11.1 |

---

#### M2.12：搜索失败诊断 Agent

**目标**：当搜索结果质量差时，AI 自动分析原因并建议修复方案。

**详细任务**：

| ID | 任务 | 优先级 | 涉及文件 | 验收标准 | 依赖 |
|----|------|--------|----------|----------|------|
| M2.12.1 | **失败模式分类规则** | P0 | `param-recommend/src/diagnosis.rs` (新) | 结构化规则：低 ID 率 / 高 mass error / 修饰异常 / DB 不匹配 | 无 |
| M2.12.1a | → 低 ID 率诊断 | P0 | 同上 | ID rate < 15% → 检查 DB, enzyme, tolerance | — |
| M2.12.1b | → 高质量偏差诊断 | P0 | 同上 | median_ppm > 10 → 建议校准或放宽 tolerance | — |
| M2.12.1c | → 修饰异常诊断 | P0 | 同上 | Oxidation > 30% → 样品制备问题 | — |
| M2.12.1d | → 数据库不匹配 | P0 | 同上 | 0 hit proteins → 检查物种是否正确 | — |
| M2.12.2 | **诊断决策树** | P0 | 同上 | 结构化 if/else 规则树 → `DiagnosisReport` struct | M2.12.1 |
| M2.12.3 | **自动参数调整建议** | P1 | 同上 | 基于诊断结果推荐具体参数修改 | M2.12.2 |
| M2.12.4 | **一键重搜** | P2 | `mcp-server/src/tools.rs` | 应用建议参数 → 自动重新搜索 | M2.12.3 |
| M2.12.5 | **Agent 定义** | P0 | `.github/agents/failure-diagnosis.agent.md` (新) | AI agent 描述，使用 diagnose 工具结合 LLM 解释 | M2.12.2 |
| M2.12.6 | **Prompt 模板** | P0 | `.github/prompts/diagnosis.prompt.md` (新) | 诊断流程 prompt | M2.12.5 |

---

#### M2.13：AI 增强功能

**目标**：强化 AI 作为核心差异化能力。

**详细任务**：

| ID | 任务 | 优先级 | 涉及文件 | 验收标准 | 依赖 |
|----|------|--------|----------|----------|------|
| M2.13.1 | **多轮分析对话 Prompt** | P1 | `.github/prompts/iterative-search.prompt.md` (新) | 调参 → 重搜 → 比较 流程模板 | M2.9 |
| M2.13.2 | **实验设计建议 Prompt** | P2 | `.github/prompts/experiment-design.prompt.md` (新) | 根据实验目标建议搜索策略 | 无 |
| M2.13.3 | **文献关联** | P2 | `report/src/literature.rs` (新) | 蛋白 accession → UniProt API → function description | 无 |
| M2.13.3a | → UniProt REST API 调用 | P2 | 同上 | 批量查询蛋白 function/GO annotation | — |
| M2.13.3b | → 缓存机制 | P2 | 同上 | 本地缓存 UniProt 查询结果 | — |
| M2.13.4 | **论文方法段落生成** | P1 | `.github/prompts/methods-section.prompt.md` (新) | 从 run_metadata → 生成 Methods section 文本 | 无 |
| M2.13.5 | **AI 引擎推荐增强** | P1 | `param-recommend/src/*.rs` | DDA → Sage/pFind; DIA → DIA-NN/Sage; 磷酸化 → pFind open | M2.0-M2.2 |

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
| `rescoring` | lib | core, fdr | Percolator/mokapot 集成 |

> 注：batch 功能集成到 `search-engine` crate，Sage/DIA-NN adapter 也在 `search-engine/src/adapters/`。

### 4.3 MCP 工具扩展（19 → ~32 个）

| 新工具 | 所属模块 | 说明 |
|--------|----------|------|
| `run_batch_search` | search-engine (batch) | 多文件批处理搜索 |
| `get_batch_status` | search-engine (batch) | 批处理进度查询 |
| `infer_proteins` | protein-inference | 蛋白推断 |
| `get_protein_groups` | protein-inference | 蛋白分组详情 |
| `quantify` | quantitation | 定量分析 |
| `normalize` | quantitation | 数据归一化 |
| `assess_quality` | qc | 质量评估 |
| `get_qc_report` | qc | QC 报告 |
| `rescore` | rescoring | Percolator 重评分 |
| `compare_results` | report | 结果比较 |
| `localize_ptms` | ptm-scoring | PTM 定位评分 |
| `diagnose_search` | param-recommend | 搜索失败诊断 |
| `export_mztab` | report | mzTab 格式导出 |

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
│ (AI: 自动选引擎) │     │ (引擎 + 参数)   │
└──────────────────┘     └────────┬────────┘
                                  │
                    ┌─────────────┼─────────────┐
                    ▼             ▼              ▼
            ┌───────────┐ ┌───────────┐ ┌───────────┐
            │   Sage    │ │  DIA-NN   │ │  pFind    │
            │ (library) │ │ (子进程)  │ │ (子进程)  │
            └─────┬─────┘ └─────┬─────┘ └─────┬─────┘
                  └─────────────┼─────────────┘
                                │ 标准 SearchResult
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
               │     compare_results          │
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
M2.0 (Sage) ──┬──→ M2.3 (批处理) ──→ M2.6 (定量, 部分由 Sage 覆盖)
               │         │
M2.1 (pFind) ─┤    M2.4 (蛋白推断) ──→ M2.6 (定量, 需要蛋白分组)
               │         │
M2.2 (DIA-NN)─┘    M2.5 (Rescoring)    M2.8 (PTM)
                         │
                    M2.7 (QC, 独立)
                         │
                    M2.9-M2.13 (Tier 3)
```

### 6.2 建议执行顺序

| 阶段 | Milestone | 依赖 | 复杂度 | 任务数 |
|------|-----------|------|--------|--------|
| **Phase 2a** | **M2.0** Sage 集成 | 无 | 高 | 13 |
| **Phase 2a** | **M2.4** 蛋白推断 | 无 (可并行) | 高 | 13 |
| **Phase 2a** | **M2.7** QC 模块 | 无 (可并行) | 中 | 9 |
| **Phase 2b** | **M2.1** pFind 对接 | 需要 pFind 二进制 | 高 | 10 |
| **Phase 2b** | **M2.2** DIA-NN 对接 | 需要 DIA-NN 二进制 | 高 | 10 |
| **Phase 2b** | **M2.3** 批处理 | M2.0 或 M2.1 | 高 | 10 |
| **Phase 2c** | **M2.5** Rescoring | M2.0 或 M2.1 | 中 | 9 |
| **Phase 2c** | **M2.6** 定量 | M2.4 | 高 | 11 |
| **Phase 2c** | **M2.8** PTM 定位 | M2.0 或 M2.1 | 中 | 7 |
| **Phase 2d** | **M2.9** 结果比较 | 任一引擎 | 低 | 8 |
| **Phase 2d** | **M2.10** 多引擎扩展 | M2.0+M2.1 | 中 | 4 |
| **Phase 2d** | **M2.11** 标准导出 | M2.4 | 中 | 5 |
| **Phase 2d** | **M2.12** 失败诊断 | 任一引擎 | 低 | 6 |
| **Phase 2d** | **M2.13** AI 增强 | 前置完成 | 低 | 5 |

**总计：~120 个细分任务**

### 6.3 Phase 2a 先行任务（可立即开始）

无需外部引擎二进制即可启动：

1. ✅ **M2.0 Sage 集成** — Rust crate 依赖，无需额外安装
2. ✅ **M2.4 蛋白推断** — 纯算法，仅依赖 PSM 结果（SimpleSearchEngine 可产出）
3. ✅ **M2.7 QC 模块** — 谱图质量评估不依赖搜索引擎
4. ✅ **M2.9 结果比较** — 已有 SearchResult 结构
5. ✅ **M2.11 mzTab 导出** — 格式转换
6. ✅ **M2.12 诊断 Agent** — AI prompt 定义

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
| Library-free DIA | DIA-NN/AlphaDIA 风格，深度学习预测 RT + fragment ion + 离子迁移率 |
| 深度学习谱图预测 | 集成 Prosit/AlphaPeptDeep 模型，用于 library 生成和重评分 |
| AlphaDIA 集成 | Mann lab DIA 引擎 adapter，AlphaPeptDeep 生态对接 |
| 差异表达分析 | 统计学显著性检验，火山图，GO 富集 |
| 分布式计算 | 多节点并行搜索 (Slurm/PBS 集群) |
| Web Dashboard | 实时监控 UI，交互式结果浏览 |
| Spectral Library | 从搜索结果构建谱图库，用于 library search |
| Glycoproteomics | 糖基化修饰搜索 |
| Cross-linking MS | 交联质谱分析 |
| 离子迁移率 (IM) 支持 | timsTOF 数据，4D 蛋白组学 |

---

## 变更记录

| 日期 | 变更内容 | 原因 |
|------|----------|------|
| 2026-04-16 | Phase 2 初始规划 | MVP 完成，进入下一阶段 |
| 2026-04-16 | 补充 DIA-NN/AlphaDIA/Sage 竞品分析 | 竞品对标遗漏 DIA 领域工具 |
| 2026-04-16 | 新增 M2.0 Sage 集成 (P0), M2.2 DIA-NN 对接 (P1) | 基于竞品分析调整引擎策略 |
| 2026-04-16 | 所有 Milestone 细化为 ~120 个具体任务 | 任务粒度细化，可直接用于开发 |

