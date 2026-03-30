# MCP Tools 参考

ProteinCopilot MCP Server 提供 8 个工具，通过 JSON-RPC over stdio 暴露给 LLM。

> **当前搜索引擎**：MVP 使用内置的 **SimpleSearchEngine**（基于 b/y 离子匹配的简化搜索），
> 后续将接入 pFind 作为生产级搜索引擎。SimpleSearch 足以验证完整流程，但搜索质量和性能不如专业引擎。

## 启动方式

```bash
cargo run --release -p protein-copilot-mcp-server
```

或通过 `.mcp.json` 配置自动发现。

---

## read_spectra

读取质谱文件，返回统计摘要。

**输入**：
```json
{
  "file_path": "/data/sample.mgf"
}
```

**输出**：`SpectrumSummary`

**使用场景**：分析数据特征，为参数推荐提供输入。

---

## get_spectrum

按 scan 号读取单张谱图。

**输入**：
```json
{
  "file_path": "/data/sample.mgf",
  "scan_number": 42
}
```

**输出**：`Spectrum`（含 mz_array、intensity_array、precursors）

---

## recommend_params

基于谱图特征推荐搜索参数。可直接传文件路径 + 数据库路径。

**输入**：
```json
{
  "file_path": "/data/sample.mgf",
  "database_path": "/data/human.fasta",
  "hints": {
    "experiment_type": "phosphorylation"
  }
}
```

**输出**：`AiDecision<SearchParams>`（含 decision、confidence、explanation）

> **注意**：`database_path` 会自动注入到推荐结果的 `decision` 中，LLM 无需手动修改。

---

## list_presets

列出所有内置搜索参数预设。

**输入**：无

**输出**：`{ "presets": [...] }`（5 个预设）

---

## run_search

执行蛋白质数据库搜索。支持两种调用方式：

**简单模式（LLM 推荐）**— 自动推荐参数：
```json
{
  "input_files": ["/data/sample.mgf"],
  "database_path": "/data/human.fasta"
}
```

**高级模式** — 传入 recommend_params 返回的参数：
```json
{
  "params": { ... },
  "input_files": ["/data/sample.mgf"]
}
```

结果自动缓存，后续用 `run_id` 引用。

**输出**：`{run_id, status: "Running", message: "..."}`（立即返回，搜索在后台执行）

> **异步模式**：run_search 不阻塞。搜索完成后用 `get_search_status(run_id)` 确认，再用 `generate_summary(run_id)` 获取结果。

---

## get_search_status

查询搜索进度。搜索完成后 status 变为 "Completed"。

**输入**：
```json
{
  "run_id": "9f71e493-..."
}
```

**输出**：`SearchProgress`
```json
{
  "run_id": "...",
  "status": "Completed",
  "progress_pct": 1.0,
  "elapsed_sec": 5.2
}
```

status 值：`Running` / `Completed` / `Failed: <reason>`

---

## check_engine

检查搜索引擎状态。

**输入**：无

**输出**：`{ "engine": {...}, "status": "Healthy" }`

---

## generate_summary

从搜索结果生成 FDR 过滤后的统计摘要。支持通过 `run_id` 引用缓存结果。

**输入**：
```json
{
  "run_id": "7ab6d7d4-df4d-4aa0-..."
}
```

或直接传 `{"result": {...}}`

**输出**：`SearchResultSummary`

---

## export_results

将搜索结果导出为文件。支持通过 `run_id` 引用缓存结果。

**输入**：
```json
{
  "run_id": "7ab6d7d4-df4d-4aa0-...",
  "output_dir": "./output"
}
```

**输出**：导出文件列表

---

## LLM 完整工作流

**最简模式（一步搜索）：**
```
① run_search(input_files, database_path)       →  {run_id, "Running"}
② get_search_status(run_id)                    →  轮询直到 "Completed"
③ generate_summary(run_id)                     →  统计摘要
④ export_results(run_id, output_dir)           →  导出文件
```

**标准模式（分步控制）：**
```
① read_spectra(file_path)                      →  数据摘要
② recommend_params(file_path, database_path)   →  推荐参数
③ run_search(decision, input_files)            →  {run_id, "Running"}
④ get_search_status(run_id)                    →  轮询直到 "Completed"
⑤ generate_summary(run_id)                     →  统计摘要
⑥ export_results(run_id, output_dir)           →  导出文件
```

LLM 全程只需传简单参数（路径、run_id），无需构造复杂 JSON。
搜索不会阻塞 — LLM 可以在等待期间与用户交互。
