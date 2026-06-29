# L4 - mcp-server（bin crate 代码组织与实现骨架）

承接 [L3-mcp-server](L3-mcp-server.md) 与 [L2](L2-architecture.md)。L3 已列全部 27 个工具清单、调用链与异步搜索模型，本篇不复述，只聚焦 `crates/mcp-server` 的**代码组织**与**实现骨架** —— main / tools / history 三文件各干什么、状态如何注入、后台搜索如何编排。所有签名、常量、路径逐一核对源码（`src/{main.rs, tools.rs, history.rs}`），禁止臆造。

## 1. 用途 + 位置 + 依赖

- 位置：`crates/mcp-server`，package 名 `protein-copilot-mcp-server`；`[[bin]] name = "protein-copilot-mcp"`，`path = "src/main.rs"`。
- 用途：workspace 的 MCP 服务 bin crate（另有 `entrapment-cli` 一个独立 CLI bin），把全部 12 个确定性库 crate 组装成 MCP 工具，经 rmcp 的 stdio transport 暴露给 LLM 客户端。`edition` / `version` / `rust-version` 均继承 workspace。
- 依赖（`Cargo.toml`）：12 个库 —— `core`、`spectrum-io`、`param-recommend`、`search-engine`、`report`、`dia-extraction`、`xic`、`result-import`、`fasta-db`、`protein-inference`、`fdr`、`entrapment-analysis`；外加 `rmcp = "1.3"`（features `server` + `transport-io`）、`tokio`、`schemars`、`uuid`、`chrono`、`tracing(-subscriber)`、`dirs = "6"`、`lru = "0.12"`、`bincode = "1"`、`csv`。dev：`tempfile`。

## 2. 代码组织

| 文件 | 行数 | 职责 |
|------|------|------|
| `main.rs` | 62 | `#[tokio::main]` 入口：初始化 tracing（`RUST_LOG` 默认 `info`，`PROTEIN_LOG_JSON=1` 切 JSON，均写 stderr）-> `ProteinCopilotServer::new()` -> `serve(stdio())` -> `service.waiting()`，错误即 `exit(1)`。 |
| `tools.rs` | 4551 | 27 个 `#[rmcp::tool(...)]` 处理器、全部输入/输出结构体、三类缓存、`ProteinCopilotServer` 状态结构、错误助手。 |
| `history.rs` | 174 | 搜索历史持久化：`SearchHistoryEntry` + `history_dir / save_entry / load_all / evict_oldest`。 |

工具注册靠两个宏：`#[rmcp::tool_router] impl ProteinCopilotServer` 收集所有 `#[rmcp::tool(name = .., description = ..)]` 方法到 `tool_router`；`#[rmcp::tool_handler] impl ServerHandler` 提供 `get_info()`。每个处理器签名同构 —— `&self` + `Parameters(input): Parameters<XxxInput>` -> `Result<Json<T>, ErrorData>`，输入结构体派生 `Deserialize + schemars::JsonSchema`（自动生成 JSON Schema）。两个真实签名：

```rust
#[rmcp::tool(name = "read_spectra", description = "Read a mass spectrometry file ...")]
fn read_spectra(&self, Parameters(input): Parameters<ReadSpectraInput>)
    -> Result<Json<SpectrumSummary>, ErrorData> { ... }

#[rmcp::tool(name = "run_search", description = "Run a proteomics database search ...")]
async fn run_search(&self, Parameters(input): Parameters<RunSearchInput>)
    -> Result<Json<SearchStarted>, ErrorData> { ... }

#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ReadSpectraInput { file_path: String }     // scan_number: u32（GetSpectrumInput）
```

`#[rmcp::tool_handler] impl ServerHandler` 仅实现 `get_info()`，返回带 `instructions` 的 `ServerInfo`，提示客户端"先 read_spectra，再 recommend_params / run_search / generate_summary"。文件类工具共享两个校验器：`validate_file_path`（非空且存在）与 `validate_scan_number`（>= 1，1-based 索引），任一不满足即返回 `INVALID_PARAMS`。少数无失败路径的工具（如 `check_engine`）直接返回裸 `Json<T>` 而非 `Result`。27 个工具按读谱 / 推参 / 搜索 / 报告 / DIA / XIC / 结果导入 / 蛋白推断 / 数据库 / 诊断 / entrapment 等类别组织，完整名字与 `description` 见 L3。

## 3. 运行时状态

状态全部挂在 `ProteinCopilotServer` 字段上（无全局变量），三类缓存 + 历史各有容量与落盘策略：

| 状态 | 类型 | 容量常量 | 溢出策略 / 路径 |
|------|------|---------|----------------|
| `run_cache` | `Arc<Mutex<OrderedRunCache>>` | `MAX_CACHE_SIZE = 100` | FIFO `evict_if_full()`，跳过 `status == "Running"` 条目 |
| `dia_cache` | `Arc<Mutex<OrderedDiaCache>>` | `MAX_DIA_CACHE_SIZE = 10` | 超量 `bincode` 落盘 `.proteincopilot/dia_cache/{uuid}.bin`，仍可取回 |
| `reader_cache` | `Arc<Mutex<lru::LruCache<PathBuf, Arc<dyn SpectrumReader>>>>` | LRU `NonZeroUsize(8)` | LRU 淘汰；`get_or_create_reader()` 先 canonicalize 再查缓存 |
| 历史 | JSON 文件 | `MAX_HISTORY = 500` | FIFO，目录 `~/.protein-copilot/history/{run_id}.json` |

其它常量：`MAX_FASTA_SIZE = 512 MiB`、`DIA_ISOLATION_WINDOW_THRESHOLD_DA = 5.0`、`RT_AUTO_LOOKUP_TOLERANCE_MIN = 0.5`、`FDR_1PCT_THRESHOLD = 0.01`。所有锁均用 `lock().unwrap_or_else(|e| e.into_inner())` 实现**锁中毒自愈** —— 某线程 panic 毒化 Mutex 后，下次取回内部数据继续服务而非级联 panic。`list_searches` 合并 `history::load_all()`（磁盘）与 `run_cache` 活跃条目。

`OrderedDiaCache` 内部用 `entries: HashMap<Uuid, Vec<Spectrum>>` + `order: Vec<Uuid>` + `extracted_at` 三件套维护 FIFO 顺序与抽取时间戳；`status()` 返回三态 `Memory{spectrum_count, extracted_at}` / `Disk{extracted_at}` / `NotFound`，`get_dia_cache_status` 据此告知谱图在内存还是已落盘。落盘若失败（目录不可写），该条目**回退保留在内存**以避免数据丢失并发 warn 日志。`get_or_create_reader` 先 `canonicalize`（失败回退原路径）再查 `reader_cache`，命中即 `Arc::clone` 复用同一 `IndexedMzMLReader`。

## 4. 异步搜索模型

`run_search` 立即返回、后台执行，三步走 + 两个查询工具：

```
run_search:
  + 先 resolve/validate params + resolve_engine    (fail-fast；失败不动任何缓存)
  + run_id = Uuid::new_v4()
  + run_cache.evict_if_full(); insert RunState{ status="Running", progress_pct=0.0 }
  + tokio::spawn:
      + PanicGuard            (Drop -> 仍为 "Running" 则置 "Failed: task panicked")
      + engine.search_with_spectra(.., on_progress, &mut diagnostics).await
      + 单锁写回 result + "Completed"（若已 "Cancelled" 不覆盖）+ history::save_entry
  -> Json<SearchStarted{ run_id, status, message }>

get_search_status(run_id) -> Json<SearchProgress>   (读 run_cache，克隆 progress)
cancel_search(run_id)     -> handle.take().abort(); status="Cancelled"  (须先为 "Running")
```

DIA 分支额外接受 `dia_run_id`：在校验通过后才从 `dia_cache.remove()` 取走谱图（失败不消耗缓存）。每个处理器进入即开 `tracing::info_span!("mcp_tool", name = ...)` 并 `info!("started"/"completed")`，结构化记录文件、run_id、引擎、容差等字段。

进度回调 `on_progress: ProgressCallback`（即 `Box<dyn Fn(SearchProgress) + Send + Sync>`）仅当条目仍 `"Running"` 时写回 `stage / progress_pct / elapsed_sec`，避免覆盖已取消的状态。`PanicGuard{ run_id, cache, start }` 以 RAII 在任务异常 unwind 时兜底置错。任务收尾的**单次加锁**里一次性写入 `SearchResult`、跑异常检测 `diagnostics.finalize(..)`、并落历史条目（含 `total_psms / psms_at_1pct_fdr / identification_rate / protein_groups`）；若期间已被 `cancel_search` 置为 `"Cancelled"`，则保留取消状态不覆盖。

## 5. 简化源码片段

下列三段均截取自源码、仅省略无关行：

main 的 stdio 启动（faithful 简化）：

```rust
#[tokio::main]
async fn main() {
    let env_filter = EnvFilter::try_from_default_env()
        .unwrap_or_else(|_| EnvFilter::new("info"));            // RUST_LOG，默认 info
    tracing_subscriber::registry().with(env_filter).with(
        fmt::layer().with_writer(std::io::stderr)               // 日志走 stderr，stdout 留给协议
            .with_span_events(fmt::format::FmtSpan::CLOSE)
            .with_timer(fmt::time::uptime())).init();
    let server = ProteinCopilotServer::new();
    let service = match server.serve(stdio()).await {           // rmcp stdio transport
        Ok(s) => s,
        Err(e) => { tracing::error!("Failed to start MCP server: {e}"); std::process::exit(1); }
    };
    if let Err(e) = service.waiting().await { std::process::exit(1); }   // 阻塞至连接关闭
}
```

一个工具处理器骨架（解析 -> 校验 -> 委托库 -> 结构化结果/错误）：

```rust
fn read_spectra(&self, Parameters(input): Parameters<ReadSpectraInput>)
    -> Result<Json<SpectrumSummary>, ErrorData> {
    let _span = tracing::info_span!("mcp_tool", name = "read_spectra").entered();
    validate_file_path(&input.file_path)?;
    let path = Path::new(&input.file_path);
    let reader = self.get_or_create_reader(path)?;              // 复用 LRU(8) 缓存的索引读取器
    let summary = reader.read_summary(path)
        .map_err(|e| mcp_core_err(CoreError::from(e)))?;        // 库错误 -> 结构化 ErrorData
    Ok(Json(summary))
}
```

锁中毒自愈 + 后台兜底守卫：

```rust
let mut cache = self.run_cache.lock().unwrap_or_else(|e| e.into_inner());  // 中毒自愈

impl Drop for PanicGuard {                       // 后台任务异常退出兜底
    fn drop(&mut self) {
        if let Ok(mut cache) = self.cache.lock() {
            if let Some(s) = cache.get_mut(&self.run_id) {
                if s.progress.status == "Running" {              // 不覆盖 Cancelled / Failed
                    s.progress.status = "Failed: task panicked".to_string();
                }
            }
        }
    }
}
```

## 6. 设计约束

- **不调 LLM**：本 crate 是纯确定性工具层，工具体只做"解析 -> 校验 -> 委托库 crate -> 返回 JSON"，FDR / 打分 / 质量偏差等数值计算全在库内完成，server 不含任何模型调用。
- **无全局可变状态**：没有 `static mut` / `lazy_static` / `OnceLock`；全部状态实例化于 `ProteinCopilotServer`，经 `new()` 构造并注入（依赖注入），run_id 用 `Uuid::new_v4()` 生成而非全局计数器；并发安全靠 `Arc<Mutex<..>>` 与 LRU。
- **结构化错误**：统一返回 `rmcp::ErrorData`。`mcp_err(code, msg)` 构造 `INVALID_PARAMS` 类错误；`mcp_core_err(CoreError)` 走 `INTERNAL_ERROR` 并附 `err.suggestion()` 修复建议。工具输入输出均为强类型结构体（派生 `serde` + `schemars::JsonSchema`），无自由文本；返回值用 `Json<T>` 包装以保证 schema 化输出。

这条"bin 薄、库厚"的边界让所有数值与算法都落在可单元测试、可复现的库 crate 内：LLM 客户端只负责编排工具调用顺序，server 把每一步固化为带 JSON Schema 的结构化输入输出，二者职责清晰隔离 —— 这也是 L2 中"确定性内核 + 智能编排层"约束在最上层的落地。

## 7. 测试入口

```bash
cargo test -p protein-copilot-mcp-server --offline
```

| 测试二进制 | 通过数 |
|-----------|--------|
| `unittests src/main.rs`（含 `tools::tests` 5 + `history::tests` 2） | 7 |
| **合计** | **7** |

覆盖：失败引擎不消耗 DIA 缓存、不留孤儿 run、锁中毒自愈、pFind 导入路径穿越拦截 / 输入排序、历史 serde 往返与目录创建 —— 正是三类缓存与后台搜索容错路径的回归护栏。

返回目录 [README](README.md)。
