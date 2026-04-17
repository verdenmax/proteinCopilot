# 失败诊断功能设计文档

> 日期: 2026-04-17
> 状态: 已通过
> 方案: 方案 C（混合）— 轻量增强现有工具 + 新增专用诊断工具

## 1. 概述

### 问题
当前 ProteinCopilot 的错误处理存在以下不足：
1. **搜索失败时**：`get_search_status` 只返回 `"Failed: {error_string}"`，无结构化分类、失败阶段、修复建议
2. **结果质量差时**：搜索"成功"但鉴定率极低，LLM 缺乏诊断数据来分析原因
3. **没有重试方案**：用户不知道该调整哪些参数

### 目标
实现全面诊断系统：搜索失败诊断 + 结果质量评估 + 参数调优建议 + 重试方案生成。

### 设计原则
- **Rust 做确定性计算**：异常检测规则、阶段指标收集、建议生成
- **LLM 做推理解释**：通过 Prompt 引导 LLM 解读结构化诊断数据
- **不在 Rust 中调用 LLM API**

---

## 2. 架构

### 分层设计

| 层 | 职责 | 实现 |
|---|---|---|
| **Layer 1: 结构化错误** | 错误分类、失败阶段、每个错误附带建议 | Rust `ErrorCategory` 枚举 |
| **Layer 2: 搜索诊断报告** | 每阶段指标、异常检测、参数匹配度 | Rust `SearchDiagnostics` |
| **Layer 3: 智能建议** | LLM 驱动的根因分析、参数调整、重试方案 | Prompt + Agent |

### 数据流

```
run_search → 搜索引擎(收集 SearchDiagnostics) → RunState(存储)
                                                     ↓
get_search_status ← error_category + has_diagnostics (轻量预览)
                                                     ↓
diagnose_search ← 完整 SearchDiagnostics (按需获取)
                                                     ↓
LLM(通过 failure-diagnosis.prompt.md 解读) → 用户
```

---

## 3. 核心数据结构

### 3.1 SearchDiagnostics

放在 `crates/core/src/diagnostics.rs`。

```rust
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SearchDiagnostics {
    /// 错误分类（搜索失败时）
    pub error_category: Option<ErrorCategory>,
    /// 失败发生的阶段名
    pub failure_stage: Option<String>,
    /// 错误详情
    pub error_detail: Option<String>,
    /// 各阶段指标
    pub stages: Vec<DiagnosticStage>,
    /// 检测到的异常
    pub anomalies: Vec<SearchAnomaly>,
    /// 修复/优化建议
    pub suggestions: Vec<DiagnosticSuggestion>,
    /// 搜索总耗时（秒）
    pub total_elapsed_sec: f64,
}
```

### 3.2 ErrorCategory

```rust
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum ErrorCategory {
    /// 输入数据问题：文件损坏、格式错误、无 MS2 谱图
    InputData,
    /// 搜索参数问题：容差不合理、酶不匹配
    Parameters,
    /// 数据库问题：物种不匹配、FASTA 格式错误、数据库过小
    Database,
    /// 搜索引擎问题：内部错误、资源不足(OOM)
    Engine,
}
```

### 3.3 DiagnosticStage

```rust
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct DiagnosticStage {
    /// 阶段名称: "file_reading", "fasta_parsing", "digestion", "matching", "fdr_calculation"
    pub name: String,
    /// 状态: "completed", "failed", "skipped"
    pub status: String,
    /// 该阶段耗时（秒）
    pub elapsed_sec: f64,
    /// 处理的项目数量（谱图数、蛋白数、候选肽段数等）
    pub items_processed: Option<u64>,
    /// 总项目数量
    pub items_total: Option<u64>,
    /// 附加信息
    pub detail: Option<String>,
}
```

### 3.4 SearchAnomaly

```rust
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SearchAnomaly {
    /// 严重程度: "warning", "error"
    pub severity: String,
    /// 异常分类
    pub category: AnomalyCategory,
    /// 人类可读的异常描述
    pub message: String,
    /// 相关指标名称
    pub metric_name: Option<String>,
    /// 指标实际值
    pub metric_value: Option<f64>,
    /// 预期范围（人类可读）
    pub expected_range: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum AnomalyCategory {
    /// PSM 鉴定率 < 10%
    LowIdentificationRate,
    /// FDR 分布异常（decoy 过多或异常分布）
    HighFdr,
    /// 无 decoy 命中（FDR 不可计算）
    NoDecoyHits,
    /// 前体容差过窄导致匹配候选不足
    NarrowTolerance,
    /// 前体容差过宽导致 FDR 不可靠
    WideTolerance,
    /// 谱图碎片离子数过少
    LowSpectraQuality,
    /// 物种可能不匹配
    DatabaseMismatch,
    /// 搜索耗时异常（匹配阶段 > 90% 总时间）
    SlowSearch,
}
```

### 3.5 DiagnosticSuggestion

```rust
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct DiagnosticSuggestion {
    /// 优先级 1-5（1 最高）
    pub priority: u8,
    /// 建议的操作
    pub action: String,
    /// 建议的理由
    pub reason: String,
    /// 具体参数调整建议（字段名 → 新值）
    pub param_changes: Option<HashMap<String, serde_json::Value>>,
}
```

### 3.6 SearchDiagnostics 构建方法

```rust
impl SearchDiagnostics {
    pub fn new() -> Self { /* 空诊断 */ }
    pub fn begin_stage(&mut self, name: &str) { /* 记录开始时间 */ }
    pub fn end_stage(&mut self, items_processed: Option<u64>) { /* 计算耗时 */ }
    pub fn fail_stage(&mut self, detail: &str) { /* 标记阶段失败 */ }
    pub fn set_error(&mut self, category: ErrorCategory, detail: &str) { /* 设置错误信息 */ }
    pub fn finalize(&mut self) { /* 运行异常检测规则 */ }
}
```

---

## 4. 异常检测规则

在 `SearchDiagnostics::finalize()` 中运行，纯确定性逻辑：

| 异常类型 | 触发条件 | 自动建议 |
|---------|---------|---------|
| LowIdentificationRate | PSM FDR 1% 鉴定率 < 10% | 检查物种、酶、修饰、容差 |
| NoDecoyHits | decoy PSM 数量 = 0 | 数据库可能缺少 decoy 序列 |
| HighFdr | 1% FDR 下 PSM < 50 | 放宽 FDR 或调整搜索参数 |
| NarrowTolerance | 前体容差 < 5 ppm 且鉴定率低 | 尝试 10-20 ppm |
| WideTolerance | 前体容差 > 50 ppm 或 > 0.5 Da | 缩小容差 |
| LowSpectraQuality | 中位碎片离子数 < 10 | 数据质量不足 |
| DatabaseMismatch | 鉴定率 < 5% 且数据库蛋白数正常 | 确认物种 |
| SlowSearch | 匹配阶段耗时 > 总时间 90% | 缩小搜索空间或使用 Sage |

`finalize()` 需要接收 `SearchResult`（用于计算鉴定率等指标）和 `SearchParams`（用于判断容差范围）。

---

## 5. MCP 工具变更

### 5.1 get_search_status 增强

在 `SearchProgress` 中新增 2 个可选字段：

```rust
pub struct SearchProgress {
    // ... 现有字段保持不变 ...
    /// 错误分类（仅失败时有值）
    pub error_category: Option<ErrorCategory>,
    /// 是否有详细诊断数据可用
    pub has_diagnostics: bool,
}
```

兼容性：这两个字段带默认值，不影响现有消费者。

### 5.2 新增 diagnose_search 工具

```
名称: diagnose_search
描述: 获取搜索运行的诊断报告，包括阶段指标、异常检测和修复建议。
      支持失败搜索（分析错误原因）和成功搜索（评估结果质量）。
      搜索完成后调用（无论成功/失败），需先通过 get_search_status 确认搜索已结束。

输入: { run_id: String }
输出: SearchDiagnostics (JSON)
```

---

## 6. 搜索引擎集成

### 6.1 SimpleSearch 集成

在 `SimpleSearchEngine::search()` 的各阶段插入诊断回调：
- file_reading: 谱图读取阶段
- fasta_parsing: FASTA 解析阶段
- digestion: 蛋白酶切阶段
- matching: 谱图-肽段匹配阶段
- fdr_calculation: FDR 计算阶段

### 6.2 Sage 集成

在 `SageAdapter::search()` 的各阶段插入诊断回调：
- file_reading: 谱图读取/转换
- fasta_parsing: FASTA → IndexedDatabase
- matching: rayon 并行匹配
- rescoring: LDA 重打分
- fdr_calculation: 三级 FDR

### 6.3 SearchEngineAdapter trait 变更

在 `search()` 方法签名中添加 `&mut SearchDiagnostics` 参数：

```rust
pub trait SearchEngineAdapter: Send + Sync {
    async fn search(
        &self,
        params: &SearchParams,
        input_files: &[PathBuf],
        progress_cb: Option<&ProgressCallback>,
        diagnostics: &mut SearchDiagnostics,  // 新增
    ) -> Result<SearchResult, CoreError>;
}
```

---

## 7. RunState 扩展

```rust
struct RunState {
    progress: SearchProgress,
    result: Option<SearchResult>,
    handle: Option<tokio::task::JoinHandle<()>>,
    diagnostics: Option<SearchDiagnostics>,  // 新增
    params: Option<SearchParams>,            // 新增（用于异常检测规则）
}
```

---

## 8. Prompt 设计

### 8.1 failure-diagnosis.prompt.md

```markdown
---
mode: agent
description: "搜索失败诊断 — 分析搜索失败原因、评估结果质量、提供参数调优建议"
---

# 搜索诊断

## 使用时机
- 搜索失败时（status = "Failed"）
- 搜索成功但 generate_summary 显示异常指标时
- 用户主动要求分析搜索质量

## 诊断流程

### 1. 获取诊断数据
调用 diagnose_search(run_id) 获取 SearchDiagnostics

### 2. 解读诊断（按场景）

#### 搜索失败
- error_category 确定大方向（InputData/Parameters/Database/Engine）
- failure_stage 定位失败阶段
- stages[] 展示"搜索走了多远"
- suggestions[] 按 priority 排序展示修复方案
- 如有 param_changes → 提供修改后的参数供重试

#### 搜索成功但质量异常
- anomalies[] 列出检测到的异常
- 对每个异常：说明含义 + 影响
- suggestions[] 提供优化建议
- 询问是否用调整后参数重新搜索

#### 搜索成功且正常
- 展示各阶段耗时（性能概况）
- 确认无异常
- 建议下一步（蛋白推断/报告导出）

### 3. 领域参考值
- HeLa DDA 鉴定率: 15-40%
- FDR 1% PSM 数量: > 1000（通常）
- Sage 搜索: 1-5 min / SimpleSearch: 5-30 min
```

---

## 9. Agent 集成

### 9.1 proteomics-search.agent.md 更新

新增诊断工具声明：
```yaml
tools:
  - diagnose_search  # 新增
```

新增诊断工作流章节：
```markdown
## 搜索诊断工作流

### 自动诊断
1. 搜索完成后（无论成功/失败），检查 get_search_status 的 has_diagnostics 字段
2. 如果 has_diagnostics = true 且 (status = "Failed" 或 error_category 有值)：
   - 自动调用 diagnose_search(run_id)
3. 基于诊断数据向用户解释原因和建议

### 手动诊断
1. 搜索成功后，用户觉得结果不理想时可手动请求诊断
2. 调用 diagnose_search(run_id) 评估结果质量
3. 展示异常和优化建议

### 重试搜索
- 如果 suggestions 中包含 param_changes，向用户展示调整后的参数
- **必须用户确认后**才能使用新参数调用 run_search
- 保留原始 run_id 供对比
```

### 9.2 决策边界更新

| 操作 | 自动/手动 |
|------|----------|
| 调用 diagnose_search | ✅ 可自动执行 |
| 解读诊断结果 | ✅ LLM 自动解读并展示 |
| 参数调整后重试 | ⚠️ 必须用户确认 |

---

## 10. 修改范围汇总

| 文件 | 操作 | 说明 |
|------|------|------|
| `crates/core/src/diagnostics.rs` | 创建 | SearchDiagnostics、DiagnosticStage、SearchAnomaly、DiagnosticSuggestion、ErrorCategory、AnomalyCategory |
| `crates/core/src/lib.rs` | 修改 | 新增 `pub mod diagnostics` |
| `crates/core/src/progress.rs` | 修改 | SearchProgress 增加 error_category + has_diagnostics |
| `crates/core/src/engine.rs` | 修改 | SearchEngineAdapter::search() 签名增加 diagnostics 参数 |
| `crates/search-engine/src/simple_engine.rs` | 修改 | 各阶段插入诊断回调 |
| `crates/search-engine/src/adapters/sage/mod.rs` | 修改 | 各阶段插入诊断回调 |
| `crates/mcp-server/src/tools.rs` | 修改 | RunState 扩展、get_search_status 增强、diagnose_search 工具、finalize 调用 |
| `.github/prompts/failure-diagnosis.prompt.md` | 创建 | 诊断解读 Prompt |
| `.github/agents/proteomics-search.agent.md` | 修改 | 新增 diagnose_search 工具 + 诊断工作流 |
