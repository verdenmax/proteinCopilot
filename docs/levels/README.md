# ProteinCopilot 分层文档（L1-L4）

按"深度递进"组织的文档，从整体认知一路读到 crate 核心源码。建议按 L1 -> L2 -> L3 -> L4 顺序阅读，也可按需直达。

## 阅读顺序

| 级别 | 关注 | 适合 |
|------|------|------|
| **L1** 总览 | 这是什么、解决什么问题、有哪些能力 | 第一次接触项目 |
| **L2** 架构 | 分层、15 个 crate 如何协作、数据流 | 想理解整体设计 |
| **L3** 子系统 | 单个子系统内部模块、数据流、伪代码 | 要改/读某一类功能 |
| **L4** crate | 单 crate 的结构体、核心函数、源码片段 | 要动具体代码 |

## L1 — 项目总览
- [L1-overview.md](L1-overview.md)

## L2 — 系统架构
- [L2-architecture.md](L2-architecture.md)

## L3 — 子系统
- [L3-spectrum-io.md](L3-spectrum-io.md) — 谱图读取与解析
- [L3-search-engine.md](L3-search-engine.md) — 搜索引擎调度与匹配
- [L3-fdr-protein.md](L3-fdr-protein.md) — FDR + 蛋白推断
- [L3-entrapment.md](L3-entrapment.md) — entrapment 同源性分级
- [L3-xic-dia.md](L3-xic-dia.md) — XIC 提取 + DIA
- [L3-result-import.md](L3-result-import.md) — 外部结果导入
- [L3-mcp-server.md](L3-mcp-server.md) — MCP Server / 工具层

## L4 — 逐 crate 核心路径
core, spectrum-io, param-recommend, search-engine, fdr, protein-inference,
xic, dia-extraction, result-import, report, fasta-db, entrapment-analysis,
entrapment-cli, mcp-server, integration-tests

> 约定：纯 ASCII 图，无 mermaid/图片；术语用蛋白组学领域词；代码块用简化片段，完整逻辑以源码为准。
