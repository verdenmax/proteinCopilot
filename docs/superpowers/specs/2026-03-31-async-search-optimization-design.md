# Post-MVP: 异步搜索优化 — 设计文档

> **日期**：2026-03-31
> **状态**：待实施
> **范围**：mcp-server, search-engine (core trait), Agent 指令

---

## 1. 问题陈述

当前异步搜索已实现基础架构（tokio::spawn + RunCache + get_search_status），但存在以下不足：

1. **无真实进度** — progress_pct 只有 0%（Running）和 100%（Completed），中间无更新
2. **无取消能力** — 长搜索无法中止，只能等待完成或失败
3. **无历史查询** — 服务重启后所有搜索记录丢失，无法回溯之前的运行
4. **Agent 轮询策略不明确** — 缺少关于轮询间隔和超时的指令

## 2. 设计目标

| # | 目标 | 验收标准 |
|---|------|---------|
| G1 | 阶段级进度上报 | get_search_status 返回当前阶段名称和阶段进度 |
| G2 | 即时取消搜索 | cancel_search(run_id) 在 1 秒内终止搜索任务 |
| G3 | 历史持久化与查询 | 服务重启后仍可查询已完成搜索的元数据 |
| G4 | Agent 轮询指导 | Agent 指令明确轮询策略和取消场景 |
| G5 | pFind 兼容 | trait 设计预留 pFind SSH 远程取消和进度上报能力 |

## 3. 架构变更

### 3.1 SearchEngineAdapter trait 扩展

**文件**：`crates/core/src/engine.rs`

注意：`ProgressCallback` 引用 `SearchProgress`。当前 `SearchProgress` 定义在 `search-engine` crate 中。为避免循环依赖，将 `SearchProgress` 移到 `core` crate（因为它是共享数据结构），然后 `search-engine` crate re-export 它。

```rust
/// 进度回调类型：搜索引擎在内部按阶段调用此回调更新进度。
pub type ProgressCallback = Box<dyn Fn(SearchProgress) + Send + Sync>;

pub trait SearchEngineAdapter: Send + Sync {
    /// 执行搜索。on_progress 回调用于报告阶段进度。
    async fn search(
        &self,
        params: &SearchParams,
        input_files: &[PathBuf],
        on_progress: ProgressCallback,
    ) -> Result<SearchResult, CoreError>;

    fn engine_info(&self) -> EngineInfo;

    async fn health_check(&self) -> Result<HealthStatus, CoreError>;

    /// 取消一个正在运行的搜索。
    /// 默认实现为空操作（依赖 JoinHandle::abort() 终止本地 task）。
    /// pFind adapter 可覆写为 SSH kill 远程进程。
    async fn cancel(&self, _run_id: Uuid) -> Result<(), CoreError> {
        Ok(())
    }
}
```

**变更影响**：
- `SimpleSearchEngine::search()` — 增加 `on_progress` 参数，在 4 个阶段回调
- `PFindAdapter::search()` — 桩实现增加参数（暂不回调）
- `PFindAdapter::cancel()` — 桩实现返回 "not implemented"
- mcp-server `run_search` — 构造回调闭包，捕获 RunCache 写入进度

### 3.2 SearchProgress 扩展

**文件**：`crates/search-engine/src/progress.rs`

```rust
pub struct SearchProgress {
    pub run_id: Uuid,
    pub status: String,           // "Running" | "Completed" | "Failed: ..." | "Cancelled"
    pub stage: Option<String>,    // 新增：当前阶段名 e.g. "Matching spectra (450/1000)"
    pub progress_pct: Option<f64>,
    pub elapsed_sec: f64,
    pub estimated_remaining_sec: Option<f64>,
}
```

**搜索阶段定义**（SimpleSearchEngine）：

| 阶段 | stage 值 | progress_pct 范围 |
|------|---------|------------------|
| 读取 FASTA | `"Reading FASTA database"` | 0% → 5% |
| 酶切消化 | `"Digesting proteins"` | 5% → 15% |
| 谱图匹配 | `"Matching spectra (N/M)"` | 15% → 90% |
| 统计聚合 | `"Aggregating results"` | 90% → 100% |

谱图匹配阶段的进度按已处理谱图数 / 总谱图数线性插值。

### 3.3 取消搜索

**新增 MCP Tool**：`cancel_search`

```
输入：{ run_id: String }
输出：{ run_id, status: "Cancelled", message }
```

**实现流程**：
1. 验证 run_id 存在且状态为 "Running"
2. 调用 `engine.cancel(run_id)` — SimpleSearchEngine 为空操作，pFind 未来走 SSH kill
3. 调用 `handle.abort()` — 终止 tokio task
4. 更新 RunState：status = "Cancelled", progress_pct = None
5. PanicGuard 的 Drop 逻辑增加：若状态已是 "Cancelled" 则不覆盖为 "Failed: task panicked"

**RunState 扩展**：

```rust
struct RunState {
    progress: SearchProgress,
    result: Option<SearchResult>,
    handle: Option<JoinHandle<()>>,  // 新增：用于 abort
}
```

### 3.4 历史持久化

**存储位置**：`~/.protein-copilot/history/`

**文件格式**：每次搜索完成/失败/取消后写入 `{run_id}.json`：

```json
{
    "run_id": "uuid",
    "status": "Completed",
    "stage": null,
    "created_at": "2026-03-31T...",
    "elapsed_sec": 62.3,
    "params_used": { ... },
    "engine_info": { ... },
    "input_files": ["..."],
    "summary": {
        "total_psms": 1234,
        "psms_at_1pct_fdr": 1100,
        "identification_rate": 0.32,
        "protein_groups_at_1pct_fdr": 456
    }
}
```

注意：**不持久化完整 SearchResult**（PSM 列表可能很大），只保存元数据和摘要统计。完整结果仍在内存缓存中，用于 `generate_summary` 和 `export_results`。

**启动加载**：MCP Server 启动时扫描 `~/.protein-copilot/history/`，将历史记录加载到只读索引中。

**FIFO 策略**：历史文件超过 `MAX_HISTORY`（默认 500）时，按 created_at 删除最旧的。

**新增 MCP Tool**：`list_searches`

```
输入：{ status_filter?: String, limit?: u32 }
输出：{ searches: Vec<SearchHistoryEntry> }
```

返回所有搜索记录（内存中的活跃搜索 + 磁盘历史），按时间倒序排列。

### 3.5 MCP Tool 变更汇总

| Tool | 变更类型 | 说明 |
|------|---------|------|
| `run_search` | 修改 | 构造 progress callback，存储 JoinHandle |
| `get_search_status` | 修改 | 返回新增的 stage 字段 |
| `generate_summary` | 不变 | 通过 run_id 从内存缓存取结果 |
| `export_results` | 不变 | 通过 run_id 从内存缓存取结果 |
| `cancel_search` | **新增** | 取消正在运行的搜索 |
| `list_searches` | **新增** | 查询搜索历史（活跃 + 持久化） |

### 3.6 Agent 指令更新

**文件**：`.github/agents/proteomics-search.agent.md`

更新搜索步骤指令：

```markdown
Step 4: Execute Search
  - Call run_search(input_files, database_path) → returns run_id
  - Report to user: "Search started (run_id: xxx)"

Step 5: Monitor Progress
  - Poll get_search_status(run_id) every 5-10 seconds
  - Report stage changes to user: "Reading database...", "Matching spectra (300/1000)..."
  - If user requests cancellation, call cancel_search(run_id)
  - If status is "Completed", proceed to Step 6
  - If status starts with "Failed", report error and suggest next steps

Step 6: Generate Results
  - Call generate_summary(run_id) → SearchResultSummary
  ...
```

新增取消场景：
```markdown
Cancellation:
  - If user says "stop", "cancel", or "abort" during search, call cancel_search(run_id)
  - Confirm cancellation: "Search cancelled. Would you like to start a new search?"
```

新增历史查询场景：
```markdown
History:
  - If user asks "what searches have I run?", call list_searches()
  - Display recent searches with status, duration, and key metrics
```

## 4. 数据流

```text
用户: "搜索这批数据"
  │
  ├─ Agent: run_search(files, db) → {run_id, status: "Running"}
  │
  ├─ Agent: get_search_status(run_id)
  │  → {stage: "Reading FASTA database", progress_pct: 0.03, elapsed: 1.2s}
  │
  ├─ Agent: get_search_status(run_id)
  │  → {stage: "Matching spectra (300/1000)", progress_pct: 0.42, elapsed: 15.3s}
  │
  ├─ Agent: get_search_status(run_id)
  │  → {status: "Completed", progress_pct: 1.0, elapsed: 35.1s}
  │
  ├─ Agent: generate_summary(run_id=xxx)
  │  → SearchResultSummary
  │
  └─ 历史自动写入 ~/.protein-copilot/history/{run_id}.json
```

取消流程：
```text
用户: "停止搜索"
  │
  ├─ Agent: cancel_search(run_id)
  │  → engine.cancel(run_id)  // SimpleEngine: no-op; pFind: SSH kill
  │  → handle.abort()         // 终止 tokio task
  │  → {status: "Cancelled"}
  │
  └─ 历史写入（status: "Cancelled"）
```

## 5. 错误处理

| 场景 | 处理 |
|------|------|
| 取消不存在的 run_id | 返回 INVALID_PARAMS，"run_id not found" |
| 取消已完成的搜索 | 返回 INVALID_PARAMS，"search already completed" |
| 取消已取消的搜索 | 返回 INVALID_PARAMS，"search already cancelled" |
| 历史目录不可写 | tracing::warn 记录，不阻塞搜索流程 |
| 历史文件损坏 | 跳过该文件，tracing::warn |
| JoinHandle abort 后 PanicGuard 触发 | 检查状态是否已是 "Cancelled"，不覆盖 |

## 6. 测试策略

| 测试类型 | 覆盖内容 |
|---------|---------|
| 单元测试 | SearchProgress 新字段序列化；历史文件读写；FIFO 淘汰 |
| 集成测试 | run_search → 轮询 get_search_status（验证阶段变化）→ generate_summary |
| 集成测试 | run_search → cancel_search → 验证状态为 "Cancelled" |
| 集成测试 | 重启后 list_searches 返回之前的历史记录 |
| e2e 测试 | 更新现有 e2e_integration.rs 验证进度回调和取消 |

## 7. 不包含（YAGNI）

- 进度推送（MCP 协议不支持 server → client push，只能轮询）
- 搜索重试/恢复
- 完整 SearchResult 持久化（PSM 列表太大，只持久化摘要）
- 搜索结果比较（Phase 2 功能）
