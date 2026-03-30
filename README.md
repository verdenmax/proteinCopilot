# ProteinCopilot

AI 驱动的蛋白质组学质谱搜索与结果解释平台。

## 功能

从质谱文件到搜索报告的完整流程：

```text
质谱文件 (mgf/mzML) + FASTA 蛋白数据库
        │
  ① spectrum-io    读取解析 → SpectrumSummary
  ② param-recommend 推荐参数 → AiDecision<SearchParams>
  ③ search-engine   酶切→匹配→打分 → SearchResult
  ④ report          统计摘要 + TSV/JSON 导出
```

**支持格式**：mgf、mzML（DDA + DIA）
**搜索引擎**：内置 SimpleSearch（MVP）、pFind adapter 预留
**输出文件**：psm.tsv、peptide.tsv、protein.tsv、result.json、run_metadata.json

## 快速测试

```bash
# 读取谱图文件
cargo run -p protein-copilot-spectrum-io --example read_spectra -- <file.mgf|mzML>

# 完整搜索流程（谱图 → 参数推荐 → 搜索 → 报告导出）
cargo run --release -p protein-copilot-search-engine --example full_search -- \
  <spectrum.mgf|mzML> <database.fasta> [output_dir]
```

## 项目结构

```text
crates/
├── core/              共享领域模型（Spectrum, SearchParams, SearchResult 等）
├── spectrum-io/       谱图文件解析（mgf/mzML streaming 读取）
├── param-recommend/   参数推荐规则引擎（确定性，不调 LLM）
├── search-engine/     搜索引擎（SimpleSearch + pFind adapter 预留）
├── report/            报告生成（摘要 + TSV/JSON 导出）
└── mcp-server/        MCP Server 二进制（8 tools，stdio transport）

.github/
├── agents/proteomics-search.agent.md   蛋白搜索助手 Agent
├── prompts/basic-search.prompt.md      基础搜索 Skill
└── prompts/result-interpretation.prompt.md  结果解读 Skill
```

## MCP Tools（8 个）

| Tool | 功能 |
|------|------|
| `read_spectra` | 读取谱图文件 → 统计摘要 |
| `get_spectrum` | 按 scan 读取单张谱图 |
| `recommend_params` | 推荐搜索参数 + 解释 |
| `list_presets` | 列出内置预设 |
| `run_search` | 执行数据库搜索 |
| `check_engine` | 检查引擎状态 |
| `generate_summary` | FDR 过滤统计摘要 |
| `export_results` | 导出 TSV/JSON 文件 |

## 架构原则

- **确定性/LLM 分层**：Rust 做所有计算，LLM 做意图理解和结果解释
- **MCP 协议**：所有能力通过 MCP tools 暴露给 LLM（M1.6 实现中）
- **可测试**：299 个单元/集成测试，0 clippy warnings
- **可审计**：每次搜索生成 run_id + 完整参数 + 引擎版本记录

## 当前进度

| 里程碑 | 状态 |
|--------|------|
| M1.1 core | ✅ 共享类型 + 验证 + trait |
| M1.2 spectrum-io | ✅ mgf/mzML 解析 |
| M1.3 param-recommend | ✅ 规则引擎 + 5 个预设 |
| M1.4 search-engine | ✅ SimpleSearch + pFind 预留 |
| M1.5 report | ✅ 摘要 + TSV/JSON 导出 |
| M1.6 mcp-server | ✅ 8 MCP tools + Agent + Skill |
| M1.7 integration | ✅ 端到端测试 + 文档 |

详细计划：`tasks/001-mvp-proteomics-search-platform.md`
架构设计：`docs/architecture.md`

## License

MIT
