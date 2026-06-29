# L1–L4 分层文档体系设计

> **日期**：2026-06-29
> **目标**：为 proteinCopilot 补充按"深度递进"组织的四级文档（L1 总览 → L2 架构 → L3 子系统 → L4 crate 核心路径），方便从认识项目到读懂源码逐层深入。

## 1. 背景与问题

现有文档（`docs/architecture.md`、`development.md`、`mcp-tools.md`、32 篇 spec）按主题/特性散布，缺少一条"从不懂到读源码"的连续路径。新成员或外部读者无法按深度递进理解：先有整体认知，再看架构，再看子系统数据流，最后落到具体 crate 的核心代码。

本设计建立统一的 **L1–L4 深度递进** 文档体系填补此缺口。

## 2. 目标与非目标

**目标**
- 四级深度：L1 项目总览 → L2 架构 → L3 子系统 → L4 单 crate 核心路径
- 覆盖全部 16 个 crate（L4 逐 crate 精讲核心路径）
- 中文叙述，**少配图**（仅必要的纯 ASCII 框图），重伪代码与从源码简化的片段
- 每篇独立成稿，可增量编写与审查
- 与现有 `docs/architecture.md` 等不重复，提供索引导航

**非目标**
- 不替换现有 spec / 用户指南 / HTML 可视化
- 不引入 mermaid 或图片渲染（纯 ASCII）
- 不做 API 自动生成（rustdoc 另行覆盖）

## 3. 目录结构

```
docs/levels/
├── README.md                 # 四级索引导航 + 阅读顺序
├── L1-overview.md            # 项目定位、解决的问题、能力、导览
├── L2-architecture.md        # 分层、16 crate 依赖、数据流、设计原则
├── L3-spectrum-io.md         # 谱图读取与解析子系统
├── L3-search-engine.md       # 搜索引擎调度与匹配子系统
├── L3-fdr-protein.md         # FDR + 蛋白推断子系统
├── L3-entrapment.md          # entrapment 分析子系统
├── L3-xic-dia.md             # XIC 提取 + DIA 子系统
├── L3-result-import.md       # 外部结果导入子系统
├── L3-mcp-server.md          # MCP Server / 工具层子系统
└── L4-<crate>.md ×16         # 逐 crate 核心路径精讲
```

L4 覆盖（16）：core, spectrum-io, param-recommend, search-engine, report, fdr,
dia-extraction, xic, result-import, protein-inference, fasta-db, entrapment-analysis,
entrapment-cli, mcp-server, integration-tests（测试约定）, 以及 workspace 级 1 篇收尾。

## 4. 各级内容规范

**L1（1 篇，~600–900 字）**：一句话定位、要解决的领域问题、核心能力清单、Rust 确定性 + LLM 编排分层理念、文档导览。0–1 张 ASCII 总览图。

**L2（1 篇，~1500–2500 字）**：四层职责表、16 crate 依赖（ASCII 树）、典型请求数据流（读谱→推参→搜索→FDR→报告）、核心原则（确定性/LLM 分层、adapter、可复现）。

**L3（7 篇，每篇 ~1000–1800 字）**：子系统职责、内部模块边界、关键数据结构、主流程伪代码、跨 crate 交互、错误处理。

**L4（16 篇，每篇 ~800–1500 字）**：crate 用途与依赖、关键结构体/trait、核心函数签名、1–3 段简化源码片段、调用链、测试入口。

**通用**：中文；标题层级一致；交叉引用相对链接；ASCII 优先；代码块标 `rust`/`text`。

## 5. 数据流（贯穿示例，用于 L1/L2）

```
mzML/mgf ─▶ spectrum-io ─▶ param-recommend ─▶ search-engine(+adapters)
                                                   │
                                fdr ◀── PSM ◀──────┘
                                 │
        protein-inference ◀──────┴──▶ report ─▶ MCP 结构化结果 ─▶ LLM 解释
```

## 6. 实施方式

- 单一 PR，分步增量写：先 README+L1+L2，再 7 篇 L3，再 16 篇 L4。
- 每篇基于真实源码（核对签名/路径），引用 `crates/...` 相对位置。
- 文档变更不需构建/测试；仅校对内链有效。
- 每完成一批提交一次，便于审查。

## 7. 验收

- 25 篇文档 + README 全部存在且内链可达
- 内容与当前源码一致（结构体/函数签名核对）
- 纯 ASCII，无 mermaid/图片
- README 提供清晰阅读顺序
