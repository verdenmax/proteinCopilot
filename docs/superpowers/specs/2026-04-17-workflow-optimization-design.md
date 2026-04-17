# 设计文档：ProteinCopilot 工具链优化与 Agent/Prompt 补全

> **日期**：2026-04-17
> **状态**：Draft
> **范围**：A（Agent/Prompt 更新）+ B（代码层优化）

---

## 1. 问题陈述

审计发现 ProteinCopilot 的 21 个 MCP 工具在**单独使用时功能正确**，但在**组合编排时存在 3 个断裂点和多处摩擦**，同时 Agent/Prompt 层未覆盖 Phase 2 新增的 5 个工具。

### 断裂点

| # | 位置 | 问题 |
|---|------|------|
| P0-1 | recommend_params → run_search | 类型断裂：返回 `AiDecision<SearchParams>`，run_search 只接受裸 `SearchParams`；database_path 是占位符 |
| P0-2 | DIA 缓存 | 内存上限 10，静默淘汰，无检查接口 |
| P0-3 | run_search 输入校验 | 文件/数据库路径在 tokio::spawn 内才校验，错误延迟 |

### Agent/Prompt 缺失

| 缺失 | 影响 |
|------|------|
| Agent 缺 5 个工具 | LLM 不知道 infer_proteins / list_databases / download_database / get_database_info 存在 |
| 缺 4 个 Prompt | 蛋白推断、数据库管理、Sage 搜索、批量搜索无引导 |

---

## 2. 设计方案

### 2.1 代码层：新增 `prepare_search` 复合工具

**目标**：将 7 步手动流程缩减为 4 步，不改动现有接口。

**接口定义**：

```rust
/// MCP Tool: prepare_search
///
/// 组合 recommend_params + 数据库自动解析，返回可直接用于 run_search 的参数。
/// LLM 展示给用户确认后，将 params 传给 run_search 即可。
struct PrepareSearchInput {
    /// 谱图文件路径列表
    input_files: Vec<String>,
    /// 用户提示（实验类型等），可选
    user_hints: Option<UserHints>,
    /// 目标物种（用于自动查找/下载数据库），可选
    organism: Option<String>,
    /// 直接指定数据库路径（优先于 organism 自动解析），可选
    database_path: Option<String>,
    /// 搜索引擎，可选（默认 SimpleSearch）
    engine: Option<String>,
}

struct PrepareSearchOutput {
    /// 推荐的搜索参数（database_path 已填充真实路径）
    params: SearchParams,
    /// 推荐理由
    reasoning: String,
    /// 推荐置信度
    confidence: f64,
    /// 备选方案
    alternatives: Vec<String>,
    /// 谱图摘要（供 LLM 展示）
    spectra_summary: SpectraSummary,
    /// 数据库信息（名称、路径、序列数）
    database_info: Option<DatabaseInfo>,
}
```

**内部逻辑**：
1. 调用 spectrum-io 读取谱图摘要
2. 调用 param-recommend 推荐参数
3. 如果提供了 organism：
   - 查询 fasta-db 本地缓存
   - 如果已缓存，直接使用路径
   - 如果未缓存，自动下载并使用新路径
4. 数据库路径优先级：用户在 input 中显式传 `database_path` > organism 自动解析 > 默认空（报错提示用户指定）
5. 返回完整的 PrepareSearchOutput

**工作流变化**：

```
之前（7 步）：
  read_spectra → recommend_params → list_databases → download_database
  → 手动提取 .decision → 手动替换 database_path → run_search

之后（4 步）：
  prepare_search(files, organism="human")
  → LLM 展示参数给用户确认
  → run_search(params, files)
  → get_search_status / generate_summary
```

**不变的接口**：recommend_params、run_search、list_databases、download_database 全部保留，高级用户仍可单独调用。

---

### 2.2 代码层：DIA 缓存持久化

**当前状态**：`DiaCache` 是 `HashMap<String, Vec<ProcessedSpectrum>>` 内存缓存，上限 10。

**改进设计**：

1. **磁盘溢出**：缓存超限时，将最旧的条目序列化到临时目录 `{output_dir}/dia_cache/{dia_run_id}/`
   - 使用 bincode 序列化 `Vec<ProcessedSpectrum>`
   - 搜索时先查内存，miss 则从磁盘加载
   - 服务退出时清理临时文件

2. **新增 `get_dia_cache_status` 工具**：

```rust
/// MCP Tool: get_dia_cache_status
struct DiaCacheStatusInput {
    dia_run_id: String,
}

struct DiaCacheStatusOutput {
    /// 缓存是否存在
    exists: bool,
    /// 在内存还是磁盘
    location: String, // "memory" | "disk" | "not_found"
    /// 谱图数量
    spectrum_count: Option<usize>,
    /// 提取时间
    extracted_at: Option<String>,
}
```

3. **搜索前自动检查**：`run_search` 在使用 `dia_run_id` 时，自动尝试从磁盘恢复缓存，失败才报错。

---

### 2.3 代码层：输入预校验

**改动位置**：`tools.rs` 中 `run_search` 的两个代码路径（DIA 和 file-based）。

**在 `tokio::spawn` 之前增加同步校验**：

```rust
// 在 spawn 前同步校验
for f in &input_files {
    if !Path::new(f).exists() {
        return Err(mcp_err(ErrorCode::INVALID_PARAMS,
            format!("Input file not found: {}", f)));
    }
}
if !Path::new(&params.database_path).exists() {
    return Err(mcp_err(ErrorCode::INVALID_PARAMS,
        format!("Database file not found: {}", params.database_path)));
}
```

**效果**：文件不存在时立即返回 MCP 错误，不创建 run_id，LLM 立刻得到明确反馈。

---

### 2.4 Agent 更新：proteomics-search.agent.md

**补齐 tools 声明**（新增 5 个）：
```yaml
tools:
  # ... existing 16 tools ...
  - prepare_search      # NEW: 复合工具
  - get_dia_cache_status # NEW: DIA 缓存检查
  - infer_proteins
  - list_databases
  - download_database
  - get_database_info
```

**新增工作流段落**：

1. **快速搜索工作流**（使用 prepare_search）：
   - 一步完成参数推荐 + 数据库解析
   - 用户确认后直接 run_search

2. **蛋白推断工作流**：
   - 搜索完成后调用 infer_proteins
   - 解释 Parsimony 最小集 + Razor 肽段
   - 三级 FDR 结果展示

3. **数据库管理工作流**：
   - list_databases 查看缓存
   - download_database 按物种下载
   - get_database_info 查看详情

4. **check_engine 使用说明**：
   - 搜索前调用确认引擎可用
   - 展示 all_engines 列表

5. **DIA 缓存检查**（更新 DIA 工作流）：
   - 提取后用 get_dia_cache_status 确认
   - 搜索前再次确认缓存存在

**修复 DIA 阈值**：
- 统一说明：5 Da 是自动检测启发式阈值，<1 Th 是 DDA/DIA 定义级区分

---

### 2.5 新增 Prompt：protein-inference.prompt.md

**内容要点**：
- Parsimony 算法说明（最小蛋白集覆盖所有肽段）
- Razor 肽段分配（共享肽段归属到证据最多的蛋白）
- 多级 FDR：PSM 1% → 肽段 1% → 蛋白质 1%
- 调用示例：`infer_proteins(run_id, fasta_path, fdr_threshold)`
- 结果解读：蛋白组数、序列覆盖率、unique vs shared 肽段

### 2.6 新增 Prompt：database-management.prompt.md

**内容要点**：
- 物种 → 数据库名映射（human_swissprot, mouse_swissprot 等）
- 三步流程：list_databases → 判断 → download_database
- cRAP 污染物数据库说明
- Decoy 数据库生成（target-decoy 策略）
- 错误处理（网络超时、磁盘空间）

### 2.7 新增 Prompt：sage-search.prompt.md

**内容要点**：
- Sage 引擎特点（rayon 并行、LDA rescoring、三级 FDR）
- 使用方式：`run_search(params, files)` 中 `engine: "Sage"`
- 与 SimpleSearch 的差异对比
- Sage 特有参数（peptide mass range 500-5000、min_ion_index 2）
- 结果 extra 字段说明（hyperscore、discriminant_score、spectrum_q/peptide_q/protein_q）

### 2.8 新增 Prompt：batch-search.prompt.md

**内容要点**：
- 多文件输入：run_search 接受文件列表
- 建议按实验分组搜索（相同参数的文件一起）
- 跨文件结果聚合策略
- 与 infer_proteins 联合使用（多文件 PSM 合并后推断蛋白）
- 性能预期（文件数 × 谱图数 → 估算时间）

---

## 3. 不改动的部分

| 项目 | 理由 |
|------|------|
| recommend_params 接口 | 高级用户仍需单独调用 |
| run_search 接口 | 保持向后兼容 |
| annotate_spectrum / extract_xic 重复 I/O | 低优先级，未来可加共享缓存 |
| 开发类 agent 领域化 | 范围外，独立任务 |
| pFind .spectra 导入 | Phase 3 功能 |

---

## 4. 实施顺序

1. **prepare_search 工具** — 最高价值，解决核心断裂点
2. **输入预校验** — 改动小，立即生效
3. **DIA 缓存持久化 + get_dia_cache_status** — 中等改动
4. **Agent 更新** — 文档工作，依赖新工具名确定
5. **4 个新 Prompt** — 文档工作，可并行

---

## 5. 测试策略

- prepare_search：集成测试（mock 数据库下载 + 真实谱图推荐）
- 输入预校验：单元测试（不存在的文件路径 → 同步错误）
- DIA 缓存：单元测试（内存→磁盘溢出→恢复）+ 集成测试
- Agent/Prompt：人工验证（LLM 调用流程是否通顺）
