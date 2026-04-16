# ProteinCopilot

AI 驱动的蛋白质组学质谱搜索与结果解释平台。

## 功能

从质谱文件到蛋白推断的完整流程：

```text
质谱文件 (mgf/mzML) + FASTA 蛋白数据库
        │
  ① spectrum-io        读取解析 → SpectrumSummary（支持索引随机访问）
  ② param-recommend    推荐参数 → AiDecision<SearchParams>
  ③ search-engine      酶切→匹配→打分 → SearchResult
  ④ protein-inference  蛋白推断（parsimony + razor + 蛋白级 FDR + 序列覆盖率）
  ⑤ report             统计摘要 + TSV/JSON 导出
  ⑥ xic                碎片离子 XIC 提取 + Plotly.js 可视化
  ⑦ result-import      外部搜索结果导入（DIA-NN / custom JSON）
```

**支持格式**：mgf、mzML（DDA + DIA，自动检测采集模式）
**搜索引擎**：内置 SimpleSearch（MVP）、pFind adapter 预留
**输出文件**：psm.tsv、peptide.tsv、protein.tsv、result.json、run_metadata.json

**蛋白推断**：`infer_proteins(run_id)` → parsimony 最小蛋白集 + razor 肽段分配 + 蛋白级 FDR + 序列覆盖率
**DIA 工作流**：`extract_dia_precursors(file)` → 缓存提取结果 → `run_search(dia_run_id=...)` → 端到端搜索
**单谱图检查**：`extract_spectrum_precursors(file, scan)` → 查看单张 MS2 的母离子提取详情
**外部结果导入**：`import_search_results(parquet/json, mzML)` → RT 匹配扫描号 → 可直接注释/XIC
**XIC 可视化**：`extract_xic(run_id, scan)` → 碎片离子色谱图 HTML（支持 SILAC 轻重标记）
**FASTA 管理**：`list_databases` / `download_database` → 内置 UniProt 数据库注册表 + 自动缓存

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
├── core/                共享领域模型（Spectrum, SearchParams, SearchResult, ProteinGroup 等）
├── spectrum-io/         谱图文件解析（mgf/mzML streaming + indexed 随机访问）
├── param-recommend/     参数推荐规则引擎（确定性，不调 LLM）
├── search-engine/       搜索引擎（SimpleSearch + pFind adapter 预留）
├── dia-extraction/      DIA 前体离子提取（同位素模式检测 + MS1↔MS2 关联）
├── fdr/                 FDR 计算（PSM/肽段/蛋白 三级 + decoy 生成 + picked-protein）
├── protein-inference/   蛋白推断（mapper + parsimony + razor + 序列覆盖率）
├── xic/                 XIC 碎片离子色谱图提取与 Plotly.js HTML 可视化
├── result-import/       外部搜索结果导入（DIA-NN parquet / custom JSON / UnimodDb）
├── fasta-db/            FASTA 数据库管理（内置注册表 + HTTPS 下载 + 缓存）
├── report/              报告生成（摘要 + TSV/JSON 导出）
├── integration-tests/   集成测试（端到端流水线验证）
└── mcp-server/          MCP Server 二进制（20 tools，stdio transport）

.github/
├── agents/proteomics-search.agent.md   蛋白搜索助手 Agent
├── prompts/basic-search.prompt.md      基础搜索 Skill
└── prompts/result-interpretation.prompt.md  结果解读 Skill
```

## MCP Tools（20 个）

| Tool | 功能 |
|------|------|
| `read_spectra` | 读取谱图文件 → 统计摘要 |
| `get_spectrum` | 按 scan 读取单张谱图 |
| `recommend_params` | 推荐搜索参数 + 解释 |
| `list_presets` | 列出内置预设 |
| `run_search` | 异步执行数据库搜索（立即返回 run_id） |
| `get_search_status` | 查询搜索进度（阶段 + 百分比） |
| `cancel_search` | 取消正在运行的搜索 |
| `check_engine` | 检查引擎状态 |
| `generate_summary` | FDR 过滤统计摘要 |
| `export_results` | 导出 TSV/JSON 文件 |
| `list_searches` | 列出搜索历史（活跃 + 持久化） |
| `annotate_spectrum` | 谱图注释（DIA: 标注+XIC+SILAC 统一视图；DDA: 标注 only） |
| `extract_dia_precursors` | DIA MS1 前体离子提取（同位素模式检测） |
| `extract_spectrum_precursors` | 单张 MS2 谱图母离子提取（调试用） |
| `extract_xic` | 碎片离子 XIC 色谱图（支持 SILAC 轻重标记） |
| `import_search_results` | 导入外部搜索结果（DIA-NN / custom JSON） |
| `infer_proteins` | 蛋白推断（parsimony + razor + 蛋白级 FDR + 序列覆盖率） |
| `list_databases` | 列出内置 FASTA 数据库（UniProt 物种库） |
| `download_database` | 下载 FASTA 数据库到本地缓存 |
| `get_database_info` | 查询已缓存数据库的详细信息 |

## 架构原则

- **确定性/LLM 分层**：Rust 做所有计算，LLM 做意图理解和结果解释
- **MCP 协议**：所有能力通过 MCP tools 暴露给 LLM
- **三级 FDR**：PSM → 肽段 → 蛋白质，各级独立 FDR 控制
- **DDA + DIA 支持**：自动检测采集模式，DIA 数据通过 MS1 同位素模式提取前体离子后搜索
- **外部结果导入**：DIA-NN parquet / 自定义 JSON → RT 匹配 mzML 扫描号 → 标准 SearchResult
- **可测试**：600+ 个单元/集成测试，0 clippy warnings
- **可审计**：每次搜索生成 run_id + 完整参数 + 引擎版本记录

## 当前进度

| 里程碑 | 状态 |
|--------|------|
| M1.1 core | ✅ 共享类型 + 验证 + trait |
| M1.2 spectrum-io | ✅ mgf/mzML 解析 + indexed 随机访问 |
| M1.3 param-recommend | ✅ 规则引擎 + 5 个预设 |
| M1.4 search-engine | ✅ SimpleSearch + pFind 预留 |
| M1.5 report | ✅ 摘要 + TSV/JSON 导出 |
| M1.6 mcp-server | ✅ 20 MCP tools + Agent + Skill |
| M1.7 integration | ✅ 端到端测试 + 文档 |
| Post-MVP | ✅ 异步搜索 + 历史持久化 + 谱图注释 + FW-1/2/3/4/6 |
| DIA 支持 | ✅ DIA 前体提取 + 搜索集成 + 端到端工作流 |
| XIC 可视化 | ✅ 碎片离子 XIC + SILAC 轻重标记 + Plotly.js HTML |
| 统一标注+XIC | ✅ 标注+XIC 合并视图 + 客户端 SILAC + 逐离子 L/H 开关 |
| 外部结果导入 | ✅ DIA-NN parquet + custom JSON + RT 扫描匹配 + UnimodDb |
| Biology Audit | ✅ 全部审计，单位统一，score 方向规范化 |
| FASTA 管理 | ✅ 内置 UniProt 注册表 + HTTPS 下载 + 本地缓存 |
| **蛋白推断** | ✅ **parsimony + razor + 三级 FDR + 序列覆盖率 + MCP tool** |

详细计划：`tasks/001-mvp-proteomics-search-platform.md`
架构设计：`docs/architecture.md`
架构演示：`docs/architecture.html`

## License

MIT
