# L2 — 系统架构

承接 [L1](L1-overview.md)。本篇讲清 4 层职责、15 个 crate 如何分层依赖、一次搜索的数据流、以及贯穿全局的设计原则。

## 1. 四层职责

```text
+------------------------------------------------------------+
| 用户交互层   MCP Client（Copilot CLI / Claude Desktop）      |  自然语言进、结果出
+------------------------------------------------------------+
| AI 编排层    .github/agents + prompts + LLM                 |  意图、推参理由、解释、诊断
+------------------------------------------------------------+
| MCP Tool 层  crates/mcp-server（bin）                        |  27 工具，组装确定性能力
+------------------------------------------------------------+
| 计算/引擎层  12 个 library crate + 搜索引擎 adapter           |  确定性：解析/打分/FDR/推断
+------------------------------------------------------------+
```

确定性逻辑全在下两层（Rust），推理在上两层（LLM）。两者不混：Rust 不调 LLM，LLM 不算 FDR。

## 2. Workspace 与 crate 依赖

15 个 crate（12 库 + 2 bin + 1 集成测试）。`mcp-server`（bin `protein-copilot-mcp`）组装全部 12 个库；`entrapment-cli`（bin `entrapment`）是独立 CLI；`core` 是共享底座，无内部依赖。

```text
              +-------- core --------+   fasta-db   (两者无内部依赖)
              | spectrum  search_    |
              | params/result/ai_    |
              | decision/label/engine|
              +----------------------+
        +-----------+-----------+
   spectrum-io     fdr      dia-extraction
        |           |
  param-recommend  protein-inference
        |
   search-engine --(fdr, param-recommend, spectrum-io)
        +--------------+---------------+--------------+
       xic        result-import       report     entrapment-analysis
                                                       |
                                                  entrapment-cli
        +------------- mcp-server 依赖全部库 -------------+
```

依赖边（非 dev）：
- `core`、`fasta-db`：无内部依赖
- `spectrum-io`/`fdr`/`dia-extraction` -> core
- `param-recommend` -> core, spectrum-io
- `protein-inference` -> core, fdr
- `search-engine` -> core, fdr, param-recommend, spectrum-io
- `xic` -> core, search-engine, spectrum-io
- `result-import` -> core, search-engine, spectrum-io
- `report` -> core, search-engine, xic
- `entrapment-analysis` -> core, search-engine, spectrum-io
- `entrapment-cli` -> entrapment-analysis
- `mcp-server` -> 全部 12 个库

`search-engine` 与 `report` 没有环：search-engine 仅在 **dev-dependencies** 引 report。

## 3. 一次 DDA 搜索的数据流

```text
read_spectra(file) --> SpectrumSummary           [spectrum-io]
recommend_params(summary) --> AiDecision<Params>  [param-recommend]
run_search(params) --> run_id（后台执行）          [mcp-server 异步]
   +- digest > match > score > target-decoy FDR   [search-engine + fdr]
get_search_status(run_id) --> Completed
generate_summary(run_id) --> 1% FDR 统计           [report]
export_results(run_id) --> psm/peptide/protein.tsv + result.json + run_metadata.json
infer_proteins(run_id) --> parsimony + razor + 蛋白FDR + 覆盖率  [protein-inference]
```

DIA 旁路：`extract_dia_precursors(file)` -> 缓存 -> `run_search(dia_run_id=...)`。
导入旁路：`import_search_results(DIA-NN/pFind/json)` -> run_id -> 同样可注释/XIC/汇总。

## 4. 搜索引擎通过 Adapter 抽象

`core::engine::SearchEngineAdapter`（async trait）统一 `search / engine_info / health_check`；`search-engine` 内 `EngineRegistry` 注册具体实现（SimpleSearch、Sage；pFind 预留）。不同引擎原始输出各自 adapter 内解析，对外统一为 `SearchResult`（PSM/肽/蛋白三级 + score/q-value/谱图引用/修饰）。

## 5. 贯穿原则

- **可复现/可审计**：每次运行有 `run_id`，落 `run_metadata.json`（参数、输入摘要、引擎版本、时间）。
- **结构化优先**：工具 I/O 皆 JSON Schema 可描述；AI 决策统一 `AiDecision`。
- **确定性**：排序/聚合不依赖 HashMap 迭代序；FDR target-decoy 默认 1%。
- **错误处理**：每 crate 自有 `thiserror` 错误；MCP 工具返回结构化错误（码 + 描述 + 建议）；库代码不 unwrap/expect。
- **可观测**：tracing 记录每次工具调用的输入摘要/耗时/输出摘要。

## 6. 异步搜索与缓存

`run_search` 立即返回 `run_id`，搜索在后台跑；`get_search_status` 轮询，`cancel_search` 取消。run / dia 结果进有锁缓存（中毒自愈），供注释、XIC、导入、推断复用。

## 7. 往下读

各子系统内部见 L3：[谱图IO](L3-spectrum-io.md)、[搜索](L3-search-engine.md)、[FDR+推断](L3-fdr-protein.md)、[entrapment](L3-entrapment.md)、[XIC+DIA](L3-xic-dia.md)、[导入](L3-result-import.md)、[MCP](L3-mcp-server.md)。
