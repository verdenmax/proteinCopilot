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

确定性计算层（库）：
- [L4-core.md](L4-core.md) — 共享数据结构与领域 trait
- [L4-spectrum-io.md](L4-spectrum-io.md) — 谱图读取解析
- [L4-param-recommend.md](L4-param-recommend.md) — 参数推荐
- [L4-search-engine.md](L4-search-engine.md) — 酶切/匹配/打分/引擎调度
- [L4-fdr.md](L4-fdr.md) — target-decoy FDR / q-value
- [L4-protein-inference.md](L4-protein-inference.md) — parsimony / razor / 覆盖率
- [L4-xic.md](L4-xic.md) — 碎片离子 XIC 提取
- [L4-dia-extraction.md](L4-dia-extraction.md) — DIA 母离子提取
- [L4-result-import.md](L4-result-import.md) — 外部结果导入
- [L4-report.md](L4-report.md) — 摘要 + 导出 + 可视化
- [L4-fasta-db.md](L4-fasta-db.md) — FASTA 数据库管理
- [L4-entrapment-analysis.md](L4-entrapment-analysis.md) — entrapment 分级 + 报告

可执行与测试：
- [L4-mcp-server.md](L4-mcp-server.md) — MCP Server（bin，27 工具）
- [L4-entrapment-cli.md](L4-entrapment-cli.md) — entrapment 命令行（bin）
- [L4-integration-tests.md](L4-integration-tests.md) — 端到端测试 + fixtures

工程层：
- [L4-workspace.md](L4-workspace.md) — workspace 约定、构建/测试、新增 crate

> 约定：纯 ASCII 图，无 mermaid/图片；术语用蛋白组学领域词；代码块用简化片段，完整逻辑以源码为准。
