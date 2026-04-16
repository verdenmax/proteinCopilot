# ProteinCopilot 竞品深度分析与战略定位

> 2026-04-16 · 基于 DIA-NN / AlphaDIA / Sage / FragPipe / Spectronaut 调研

---

## 1. 三大开源竞品功能拆解

### 1.1 DIA-NN（Demichev, Nature Methods 2020; v2.0+ 2024–2025）

**定位**：DIA 领域事实标准，学术界使用最广泛的 DIA 分析工具。

| 模块 | 功能 | 技术细节 |
|------|------|----------|
| **Library-free 搜索** | 从 FASTA 直接生成预测谱图库 | 深度学习预测 RT + fragment ion 强度，无需 DDA 实验建库 |
| **Cranium 神经网络** | 核心 DL 引擎 | C++ 实现，peptide → (RT, fragment intensities) 预测 |
| **QuantUMS 定量** | 不确定性最小化定量 | 比 MaxLFQ 更精确，CV 通常 3-5% |
| **MBR** | Match-Between-Runs | 从 study 数据构建经验库，跨 run ID 转移 |
| **plexDIA** | 多标记 DIA | SILAC / mTRAQ / dimethyl 标记的 DIA 分析 |
| **仪器支持** | 全平台 | Orbitrap, timsTOF (PASEF, Slice-PASEF), Astral, SCIEX |
| **输出格式** | Parquet + TSV matrix | 适合下游统计分析，我们已支持 Parquet 导入 |
| **CLI** | 清晰参数化 | `--f` `--fasta` `--lib` `--gen-spec-lib` `--predictor` `--threads` |

**DIA-NN 的核心护城河**：
1. **深度学习谱图预测** — 从序列直接预测 MS2 谱图，消除了对 DDA 实验建库的依赖
2. **极致性能** — C++ 实现，号称 1000 runs/hour
3. **社区生态** — 引用量极高，几乎所有 DIA 论文都用 DIA-NN

**DIA-NN CLI 集成接口**（我们最关心的）：
```bash
# Library-free DIA 搜索
diann --f sample.mzML --fasta uniprot.fasta --gen-spec-lib --predictor \
      --out report.tsv --threads 16 --qvalue 0.01

# Library-based 搜索
diann --f sample.mzML --lib library.tsv --out report.tsv --threads 8

# 输出文件
# report.tsv          — PSM 级结果
# report.pr_matrix.tsv — precursor 定量矩阵
# report.pg_matrix.tsv — protein group 定量矩阵
# report.gg_matrix.tsv — gene group 矩阵
# report-lib.tsv       — 生成的谱图库（可复用）
```

---

### 1.2 AlphaDIA（Mann lab, Nature Biotechnology 2025）

**定位**：下一代 feature-free DIA 分析，深度学习驱动，AlphaPept 生态核心组件。

| 模块 | 功能 | 技术细节 |
|------|------|----------|
| **Feature-free 分析** | 跳过传统 peak picking | 直接在原始数据上用 ML 检测模式 |
| **Transfer Learning** | 端到端迁移学习 | DL 模型持续适配特定仪器/实验条件 |
| **AlphaPeptDeep** | 谱图库预测 | 从序列预测 RT + fragment + 离子迁移率 |
| **生态整合** | 模块化设计 | AlphaBase + AlphaRaw + AlphaTims + directLFQ |
| **离子迁移率** | 原生 IM 支持 | timsTOF 4D 数据原生处理 |
| **平台** | 跨平台 | Python + Numba, Docker, GUI/CLI |

**AlphaDIA 的核心创新**：
1. **Feature-free** — 不依赖传统色谱峰提取，直接在原始信号上做 ML 匹配
2. **Transfer Learning** — 模型在分析过程中持续优化，适配当前数据
3. **结果**：dimethyl HeLa 数据上，unique peptide +48%, protein group +25%（vs 传统方法）

**与 ProteinCopilot 的关系**：
- AlphaDIA 是 Python 生态，通过子进程 + TOML 配置驱动
- 可作为 Phase 3 的高级 DIA 引擎候选
- 优先级低于 DIA-NN（DIA-NN 更成熟、社区更大）

---

### 1.3 Sage（Lazear, J. Proteome Res. 2023）

**定位**：纯 Rust 搜索引擎，极快，开源，架构与 ProteinCopilot 最匹配。

| 模块 | 功能 | 技术细节 |
|------|------|----------|
| **Fragment Indexing** | MSFragger 级别速度 | 宽/窄 mass tolerance 都支持（>500 Da precursor window） |
| **DIA 支持** | Wide-window 搜索 | WWA/PRM/DIA 模式，2025 beta 支持 diaPASEF |
| **定量** | LFQ + TMT/iTRAQ | 内置 feature detection + RT alignment + isotope pattern |
| **FDR** | 多级 FDR | 内置 LDA + KDE，picked peptide/protein FDR |
| **Chimeric** | 嵌合谱匹配 | 单张谱图多肽段鉴定（DIA 关键能力） |
| **RT 预测** | 内置模型 | 训练回归模型，跨实验 RT 对齐 |
| **云原生** | S3/GCS/Azure | 直接流式读写云存储 |
| **配置** | JSON 文件 | 所有参数 JSON 化，与 MCP JSON-RPC 天然兼容 |
| **输出** | CSV / Parquet | Percolator/Mokapot 兼容 |

**Sage 的核心优势（对 ProteinCopilot 而言）**：
1. **同为 Rust** — 可直接作为 library crate 集成，无需子进程
2. **JSON 配置** — 与我们的 MCP JSON-RPC 协议天然兼容
3. **DDA + DIA 双模** — 一个引擎覆盖两种模式
4. **内置定量** — LFQ + TMT，省去单独开发定量模块
5. **内置 FDR** — LDA + KDE，比简单 TDC 更先进
6. **活跃维护** — MIT 开源，社区活跃

**Sage JSON 配置示例**：
```json
{
  "database": {
    "fasta": "human.fasta",
    "enzyme": { "missed_cleavages": 2, "min_len": 7, "max_len": 50 },
    "static_mods": { "C": 57.0215 },
    "variable_mods": { "M": [15.9949] }
  },
  "precursor_tol": { "ppm": [-20, 20] },
  "fragment_tol": { "ppm": [-20, 20] },
  "report_psms": 1,
  "wide_window": true,  // DIA 模式
  "chimera": true        // 嵌合谱匹配
}
```

---

## 2. ProteinCopilot 现有 DIA 能力盘点

| 能力 | 状态 | 所在 crate |
|------|------|-----------|
| DIA mzML 读取 + 隔离窗检测 | ✅ 已实现 | `spectrum-io` / `dia-extraction` |
| DIA 前体提取 (pseudo-DDA) | ✅ 已实现 | `dia-extraction` (同位素模式 + MS1 关联) |
| DIA-NN Parquet 结果导入 | ✅ 已实现 | `result-import/src/diann.rs` |
| pFind .spectra 结果导入 | ✅ 已实现 | `result-import/src/pfind.rs` |
| DIA 感知的 XIC 提取 | ✅ 已实现 | `xic` (cycle-based extraction) |
| DIA 感知的谱图注释 | ✅ 已实现 | `search-engine/src/annotate.rs` |
| SimpleSearchEngine (DDA) | ✅ 已实现 | `search-engine/src/simple_engine.rs` |
| pFind adapter (stub) | ⚠️ 框架已有 | `search-engine/src/adapters/pfind.rs` |
| DIA-NN 进程编排 | ❌ 未实现 | — |
| Sage 集成 | ❌ 未实现 | — |
| Library-free DIA | ❌ 未实现 | — |
| 蛋白推断 | ❌ 未实现 | — |
| 定量分析 | ❌ 未实现 | — |

---

## 3. 战略定位分析

### 3.1 我们不应该做什么

| 不该做 | 原因 |
|--------|------|
| ❌ 重新实现 DIA 深度学习预测 | DIA-NN / AlphaPeptDeep 已经做到极致，我们没有 ML 团队 |
| ❌ 重新实现 fragment indexing | Sage 已经是 Rust 实现，直接用 |
| ❌ 重新实现 library-free search | 这是 5+ 年研究的结果，不是工程问题 |
| ❌ 试图在 DIA 搜索算法上超越 DIA-NN | 他们的护城河太深 |

### 3.2 我们应该做什么

#### 🎯 核心差异化：AI 编排层（没有竞品有这个）

```
用户: "帮我分析这批 timsTOF DIA 数据"

ProteinCopilot AI 编排:
  1. read_spectra → 检测到 DIA 模式, timsTOF, 120 min gradient
  2. AI 推荐: "检测到 DIA-PASEF 数据, 建议使用 DIA-NN library-free 模式"
  3. 自动配置 DIA-NN CLI 参数
  4. run_diann → 启动 DIA-NN 子进程
  5. 结果导入 → import_search_results(diann_parquet)
  6. generate_summary → AI 解释结果
  7. "鉴定了 6,832 蛋白 (1% FDR)，比典型 HeLa 高 15%，数据质量优秀"
```

**DIA-NN / AlphaDIA / Sage 都没有这种自然语言交互和智能编排能力。**

#### 🔧 引擎集成优先级

| 优先级 | 引擎 | 集成方式 | 理由 |
|--------|------|----------|------|
| **P0** | **Sage** | Rust library crate 依赖 | 同语言，最深集成，DDA+DIA，内置定量+FDR |
| **P1** | **pFind** | 子进程 + .spectra/.cfg 文件 | ICR 自研，已有 adapter 框架 |
| **P1** | **DIA-NN** | 子进程 + CLI 参数 + Parquet 导入 | DIA 标准，已有结果导入 |
| **P2** | **MSFragger** | 子进程 + .params 文件 | FragPipe 生态，DDA 标准 |
| **P3** | **AlphaDIA** | 子进程 + TOML 配置 | 前沿但还不够成熟 |

#### 📊 功能开发优先级重排

基于竞品分析，Phase 2 优先级应调整为：

| # | 任务 | 理由 | 依赖 |
|---|------|------|------|
| 1 | **Sage 集成** | 一个引擎解决 DDA+DIA+定量+FDR，Rust 原生 | 无 |
| 2 | **DIA-NN 进程编排** | DIA 标准，已有 Parquet 导入基础 | 无 |
| 3 | **蛋白推断** | 所有竞品都有，缺失则不可用 | Sage 或 pFind |
| 4 | **多引擎结果比较** | AI 解释 Sage vs DIA-NN 差异 — 独特功能 | 1 + 2 |
| 5 | **pFind 对接** | ICR 自研引擎 | 无 |
| 6 | **QC 模块** | AI 驱动质量评估 — 差异化 | 1 或 2 |
| 7 | **批处理** | 多文件并行 | 1 或 2 |
| 8 | **定量分析** | Sage 内置 LFQ+TMT，只需包装 | 1 |

---

## 4. Sage 集成方案（推荐 P0）

### 4.1 为什么 Sage 应该是最高优先级

```
ProteinCopilot (Rust) ─── 直接 crate 依赖 ──→ Sage (Rust)
                                                   │
                            ┌─────────────────────┤
                            │                     │
                     DDA 搜索 + FDR        DIA wide-window
                            │                     │
                     LFQ 定量              TMT 定量
                            │                     │
                     chimeric 匹配        RT 预测
```

- **零序列化开销**：直接调用 Rust API，不需要子进程 + JSON/Parquet 序列化
- **一个集成 = 5 个功能**：搜索 + DIA + LFQ + TMT + FDR
- **JSON 配置天然兼容**：Sage 的 JSON config 可直接映射到我们的 `SearchParams`
- **错误处理一致**：Rust Result 类型，不需要解析外部进程 stderr

### 4.2 集成架构

```rust
// crates/search-engine/src/adapters/sage.rs

use sage_core::{scoring::Scorer, database::IndexedDatabase, ...};

pub struct SageAdapter {
    config: SageConfig,
}

impl SearchEngineAdapter for SageAdapter {
    async fn search(&self, params: &SearchParams, input_files: &[PathBuf])
        -> Result<SearchResult, CoreError>
    {
        // 1. 将 SearchParams → Sage JSON config
        let sage_config = self.to_sage_config(params)?;
        // 2. 调用 Sage 库 API（非子进程）
        let results = sage_core::run_search(&sage_config, input_files)?;
        // 3. 将 Sage 结果 → 标准 SearchResult
        self.convert_results(results)
    }
}
```

### 4.3 Sage 带来的"免费"功能

| 我们原来需要自己开发的 | Sage 直接提供 |
|----------------------|--------------|
| M2.3 蛋白推断 (部分) | Picked protein FDR |
| M2.4 重评分 | 内置 LDA + KDE rescoring |
| M2.5 LFQ 定量 | Feature detection + RT alignment + LFQ |
| M2.5 TMT 定量 | MS2/MS3 reporter ion 提取 |
| M2.7 chimeric 谱匹配 | `chimera: true` 配置 |
| DIA 搜索 | `wide_window: true` |

---

## 5. DIA-NN 集成方案（推荐 P1）

### 5.1 集成路径

```
ProteinCopilot
  │
  ├── 已有: import_search_results(diann_parquet) ← 结果导入
  │
  └── 需要新增:
      ├── DiaNN CLI 参数生成 (SearchParams → diann args)
      ├── 子进程编排 (spawn + progress monitoring)
      ├── 输出结果收集 (report.tsv + matrices)
      └── AI 参数推荐增强 (DIA 模式自动选择 DIA-NN)
```

### 5.2 实现细节

```rust
// crates/search-engine/src/adapters/diann.rs

pub struct DiannAdapter {
    binary_path: PathBuf,  // diann 可执行文件路径
}

impl SearchEngineAdapter for DiannAdapter {
    async fn search(&self, params: &SearchParams, input_files: &[PathBuf])
        -> Result<SearchResult, CoreError>
    {
        // 1. 构建 CLI 参数
        let args = self.build_cli_args(params, input_files)?;
        // e.g.: ["--f", "sample.mzML", "--fasta", "human.fasta",
        //        "--gen-spec-lib", "--predictor", "--threads", "16"]

        // 2. 启动子进程
        let child = tokio::process::Command::new(&self.binary_path)
            .args(&args)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()?;

        // 3. 监控进度 (解析 stderr 输出)
        // 4. 等待完成
        // 5. 导入 Parquet 结果 → 标准 SearchResult
        let result = import_diann_parquet(&output_path)?;
        Ok(result)
    }
}
```

### 5.3 DIA-NN 关键 CLI 参数映射

| SearchParams 字段 | DIA-NN CLI 参数 | 说明 |
|-------------------|----------------|------|
| `database_path` | `--fasta` | FASTA 文件 |
| `precursor_tolerance` | `--mass-acc` | 前体质量偏差 (ppm) |
| `fragment_tolerance` | `--mass-acc-ms1` | MS1 质量偏差 |
| `enzyme` | `--cut` + `--missed-cleavages` | 消化酶 |
| `variable_modifications` | `--var-mod` | 可变修饰 |
| `fixed_modifications` | `--fixed-mod` | 固定修饰 |
| `min_peptide_length` | `--min-pep-len` | 最短肽段 |
| `max_peptide_length` | `--max-pep-len` | 最长肽段 |
| `acquisition_mode: DIA` | `--gen-spec-lib --predictor` | Library-free 模式 |

---

## 6. ProteinCopilot 独特价值主张

### 6.1 没有竞品拥有的能力

| 能力 | DIA-NN | AlphaDIA | Sage | FragPipe | **ProteinCopilot** |
|------|--------|----------|------|----------|-------------------|
| 自然语言交互 | ❌ | ❌ | ❌ | ❌ | ✅ |
| AI 参数推荐 + 解释 | ❌ | ❌ | ❌ | ❌ | ✅ |
| AI 结果解释 | ❌ | ❌ | ❌ | ❌ | ✅ |
| 多引擎智能选择 | ❌ | ❌ | ❌ | 部分 | ✅ |
| 交互式谱图注释 | ❌ | ❌ | ❌ | ❌ | ✅ |
| AI 驱动 QC | ❌ | ❌ | ❌ | ❌ | ✅ |
| MCP 协议生态 | ❌ | ❌ | ❌ | ❌ | ✅ |

### 6.2 用户体验差异

**DIA-NN 用户**：
```
1. 手动选择参数（enzyme, mods, tolerance）
2. 手动构建 CLI 命令
3. 等待完成
4. 手动用 R/Python 分析 report.tsv
5. 手动判断结果质量
```

**ProteinCopilot 用户**：
```
用户: "分析这批 HeLa DIA 数据，用 human FASTA"
AI:   "检测到 DIA 数据，推荐 DIA-NN library-free。开始分析..."
      [自动: 参数推荐 → DIA-NN 搜索 → 结果导入 → QC → 解释]
AI:   "完成！6,832 蛋白 (1% FDR)，鉴定率 42%。
       ⚠️ 注意：Oxidation(M) 占比 18%（偏高），建议检查样品制备流程。
       需要我生成谱图注释或 XIC 验证吗？"
```

### 6.3 战略总结

```
ProteinCopilot 的定位不是「又一个搜索引擎」

而是：「搜索引擎的智能编排层」

                    ┌─── Sage (Rust native, DDA+DIA)
用户 ←→ AI 编排 ←→ ├─── DIA-NN (DIA 标准)
                    ├─── pFind (ICR 自研)
                    └─── MSFragger (DDA 标准)
```

**核心竞争力 = 搜索引擎集成 × AI 编排 × 交互体验**

不是取代 DIA-NN，而是让 DIA-NN 更好用。
不是重写 Sage，而是让 Sage 会说话。
