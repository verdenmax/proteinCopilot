# L1 — 项目总览

> 一句话：**ProteinCopilot 是一个 Rust workspace + MCP Server + Agent/Skill 驱动的蛋白质组学质谱智能搜索平台**。Rust 做所有确定性计算，LLM（经 MCP Client）做意图理解、参数推荐、结果解释。

## 1. 它解决什么问题

蛋白质组学质谱分析链条长且专业：原始谱图读取 -> 选搜索参数 -> 跑搜索引擎 -> 控 FDR -> 蛋白推断 -> 解释结果。传统做法需要专家手工配参、对着一堆 TSV 看结果。ProteinCopilot 把"确定性计算"封装成 MCP 工具，把"判断与解释"交给 LLM，让用户用自然语言完成全流程，并得到可解释、可复现、可审计的结果。

## 2. 核心理念：确定性与推理严格分层

```text
用户(自然语言)
   |
   v
MCP Client + LLM        <- 推理：意图理解、参数推荐理由、结果解释、失败诊断
   |  调用 MCP 工具
   v
Rust MCP Server         <- 确定性：谱图解析、酶切打分、FDR、蛋白推断、报告
   |  通过 Adapter
   v
搜索引擎 (SimpleSearch / Sage / pFind 预留)
```

- **Rust 只算不猜**：FDR、打分、质量偏差等数值计算全部在 Rust 内，禁止交给 LLM。
- **LLM 只猜不算**：推荐参数前必须先调工具拿数据特征，解释结果前必须先拿统计摘要，不"凭空"给结论。
- 所有 AI 决策输出结构化（`AiDecision`：decision/confidence/explanation/alternatives/evidence）。

## 3. 能力清单

从质谱文件到蛋白推断的完整流程：

```text
质谱文件 (mgf/mzML) + FASTA
  01 spectrum-io        读取解析 -> SpectrumSummary（索引随机访问）
  02 param-recommend    推荐参数 -> AiDecision<SearchParams>
  03 search-engine      酶切->匹配->打分 -> SearchResult（SimpleSearch + Sage）
  04 fdr                target-decoy q-value（PSM/肽/蛋白 1%）
  05 protein-inference  parsimony + razor + 蛋白级 FDR + 覆盖率
  06 report             统计摘要 + TSV/JSON 导出
  07 xic                碎片离子 XIC 提取 + Plotly 可视化（含 SILAC）
  08 dia-extraction     DIA 母离子提取（宽窗同位素分析）
  09 result-import      外部结果导入（DIA-NN / pFind / custom JSON）
  10 fasta-db           FASTA 数据库管理（UniProt 注册表 + 缓存）
  11 entrapment         陷阱库命中分类（L0-L4 同源性分级 + HTML 报告）
```

- **格式**：mgf、mzML（DDA + DIA 自动检测）。
- **引擎**：SimpleSearch（MVP）、Sage（生产级）、pFind（预留 adapter）。
- **接口**：27 个 MCP Tool，覆盖读谱、推参、搜索（异步 run_id）、状态、注释、XIC、导入、推断、数据库、entrapment。

## 4. 技术栈

Rust 2021（1.85+）；tokio 异步；serde 序列化；thiserror/anyhow；tracing 可观测；quick-xml/base64/flate2 解析 mzML；rmcp/JSON-RPC 走 MCP；sage-core 集成 Sage。

## 5. 怎么读下去

- 想懂整体设计 -> **L2 架构**：15 crate 依赖与数据流。
- 想改某类功能 -> **L3 子系统**：谱图IO / 搜索 / FDR+推断 / entrapment / XIC+DIA / 导入 / MCP。
- 想动具体代码 -> **L4 crate**：结构体、核心函数、源码片段。

文档总入口见 [README.md](README.md)。
