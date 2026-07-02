# ProteinCopilot

AI 驱动的蛋白质组学质谱搜索与结果解释平台。

## 功能

从质谱文件到蛋白推断的完整流程：

```text
质谱文件 (mgf/mzML) + FASTA 蛋白数据库
        │
  ① spectrum-io        读取解析 → SpectrumSummary（支持索引随机访问）
  ② param-recommend    推荐参数 → AiDecision<SearchParams>
  ③ search-engine      酶切→匹配→打分 → SearchResult（SimpleSearch + Sage）
  ④ protein-inference  蛋白推断（parsimony + razor + 蛋白级 FDR + 序列覆盖率）
  ⑤ report             统计摘要 + TSV/JSON 导出
  ⑥ xic                碎片离子 XIC 提取 + Plotly.js 可视化
  ⑦ result-import      外部搜索结果导入（DIA-NN / custom JSON）
  ⑧ fasta-db           FASTA 数据库管理（UniProt 注册表 + 缓存）
  ⑨ diagnostics        搜索失败诊断 + 质量异常检测 + 修复建议
  ⑩ entrapment         陷阱库命中分类（L0-L4 同源性分级 + HTML 报告）
```

**支持格式**：mgf、mzML（DDA + DIA，自动检测采集模式）
**搜索引擎**：内置 SimpleSearch（MVP）、Sage（生产级，sage-core 库集成）、pFind adapter 预留
**输出文件**：psm.tsv、peptide.tsv、protein.tsv、result.json、run_metadata.json

**蛋白推断**：`infer_proteins(run_id)` → parsimony 最小蛋白集 + razor 肽段分配 + 蛋白级 FDR + 序列覆盖率
**DIA 工作流**：`extract_dia_precursors(file)` → 缓存提取结果 → `run_search(dia_run_id=...)` → 端到端搜索
**单谱图检查**：`extract_spectrum_precursors(file, scan)` → 查看单张 MS2 的母离子提取详情
**外部结果导入**：`import_search_results(parquet/json, mzML)` → RT 匹配扫描号 → 可直接注释/XIC
**XIC 可视化**：`extract_xic(run_id, scan)` → 碎片离子色谱图 HTML（支持 SILAC 轻重标记）
**Sage 搜索引擎**：在 `run_search` 中指定 `engine: "Sage"` 即可使用 sage-core 进行生产级蛋白组学搜索（rayon 并行打分 + LDA rescoring）
**FASTA 管理**：`list_databases` / `download_database` → 内置 UniProt 数据库注册表 + 自动缓存
**搜索诊断**：`diagnose_search(run_id)` → 阶段耗时 + 7 条异常检测规则 + 分级修复建议
**陷阱库分析**：`classify_entrapment_hits` → L0-L4 同源性分级 + Levenshtein edit distance 跨长匹配 + SubstitutionType 注释（Q/K、等质量二肽等）+ HTML 交互报告 + CLI 独立工具

## 安装与使用

ProteinCopilot 以**单个 MCP 二进制** `protein-copilot-mcp` 发布。用户在自己的 MCP 客户端中登记该二进制即可使用，无需阅读源码。

### 安装二进制

**方式一：npm（推荐，无需安装 Rust/手动下载）**

```bash
# 直接运行（无需安装）
npx protein-copilot-mcp --list-tools

# 或全局安装
npm install -g protein-copilot-mcp
protein-copilot-mcp --list-tools
```

`npx` 会自动下载对应平台的预编译二进制并运行，支持 Linux (x64 glibc/musl)、macOS (Intel/Apple Silicon)、Windows (x64)。

**方式二：下载预编译二进制**

从 [GitHub Releases](https://github.com/verdenmax/proteinCopilot/releases) 下载对应平台的压缩包并解压：

| 平台 | 资产 |
|------|------|
| Linux x86_64（静态链接，推荐） | `protein-copilot-mcp-x86_64-unknown-linux-musl.tar.gz` |
| Linux x86_64（glibc） | `protein-copilot-mcp-x86_64-unknown-linux-gnu.tar.gz` |
| macOS Intel | `protein-copilot-mcp-x86_64-apple-darwin.tar.gz` |
| macOS Apple Silicon | `protein-copilot-mcp-aarch64-apple-darwin.tar.gz` |
| Windows x86_64 | `protein-copilot-mcp-x86_64-pc-windows-msvc.zip` |

```bash
# Linux/macOS 安装示例
tar xzf protein-copilot-mcp-*.tar.gz
sudo cp protein-copilot-mcp-*/protein-copilot-mcp /usr/local/bin/
chmod +x /usr/local/bin/protein-copilot-mcp

# Windows 安装示例（PowerShell）
Expand-Archive protein-copilot-mcp-*.zip -DestinationPath .
# 将 protein-copilot-mcp.exe 放到 PATH 中的目录，如 C:\Users\<你>\bin\
```

**方式三：用 Cargo 从源码安装**（需 Rust 1.85+）

```bash
cargo install --git https://github.com/verdenmax/proteinCopilot.git \
  -p protein-copilot-mcp-server
# 安装后二进制名为 protein-copilot-mcp
```

**方式四：Docker**

```bash
docker build -t protein-copilot-mcp .
docker run -i --rm protein-copilot-mcp        # 通过 stdio 提供 MCP 服务
```

### 配置 MCP 客户端

以下配置均使用 `npx protein-copilot-mcp`（无需手动安装二进制）。如果手动安装了二进制，将 `"command"` 改为 `"<path>/protein-copilot-mcp"` 并去掉 `"args"`。如果用 Docker，将 `"command"` 替换为 `"docker"`，`"args"` 设置为 `["run", "-i", "--rm", "protein-copilot-mcp"]`。

#### Claude Code

**项目级**（在项目根目录创建 `.mcp.json`，可提交到 git）：

```json
{
  "mcpServers": {
    "protein-copilot": {
      "command": "npx",
      "args": ["-y", "protein-copilot-mcp"],
      "timeout": 300,
      "env": {
        "RUST_LOG": "info"
      }
    }
  }
}
```

> 如果你手动安装了二进制（npm install -g / Cargo / 下载），将 `"command"` 改为 `"<path>/protein-copilot-mcp"` 并去掉 `"args"`。

**用户级**（全局生效，写入 `~/.claude.json` 的 `"mcpServers"` 字段）：

```json
{
  "mcpServers": {
    "protein-copilot": {
      "command": "<path>/protein-copilot-mcp",
      "timeout": 300,
      "env": {
        "RUST_LOG": "info"
      }
    }
  }
}
```

#### Claude Desktop

在 Claude Desktop 的配置文件（`~/Library/Application Support/Claude/claude_desktop_config.json`（macOS）或 `%APPDATA%\Claude\claude_desktop_config.json`（Windows））中添加：

```json
{
  "mcpServers": {
    "protein-copilot": {
      "command": "npx",
      "args": ["-y", "protein-copilot-mcp"],
      "env": {
        "RUST_LOG": "info"
      }
    }
  }
}
```

#### GitHub Copilot (VS Code)

在项目根目录的 `.vscode/mcp.json` 中配置（支持工作区级共享）：

```json
{
  "servers": {
    "protein-copilot": {
      "command": "npx",
      "args": ["-y", "protein-copilot-mcp"],
      "env": {
        "RUST_LOG": "info"
      }
    }
  }
}
```

> **注意**：需要 VS Code 1.99+、Copilot Agent mode 开启，且 `chat.mcp.discovery.enabled` 设为 `true`。配置完成后在 Agent mode 下点击 tools 图标即可看到 `protein-copilot` 工具。

#### OpenAI Codex CLI

在 `~/.codex/config.toml` 中添加（需 Codex CLI ≥ 0.43）：

```toml
[mcp_servers.protein-copilot]
command = "npx"
args = ["-y", "protein-copilot-mcp"]
env = { RUST_LOG = "info" }
startup_timeout_ms = 300_000
```

### 验证安装

```bash
# 在终端直接查看工具列表，确认二进制可用
protein-copilot-mcp --list-tools          # 文本目录：每个工具的摘要
protein-copilot-mcp --list-tools --json   # 完整 JSON Schema
protein-copilot-mcp --help                # 用法

# 测试 MCP 通信（启动 stdio server 后发送 initialize 请求）
echo '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{},"clientInfo":{"name":"test","version":"1.0"}}}' | protein-copilot-mcp
```

完整接口契约见 [`docs/mcp-tools.md`](docs/mcp-tools.md)；分层文档见 [`docs/levels/`](docs/levels/README.md)。

### 典型用法（自然语言驱动）

接入客户端后，用自然语言即可走完整流程——LLM 会按需调用下列工具：

```text
你：帮我分析 /data/sample.mzML，数据库用人类的。
  -> read_spectra(/data/sample.mzML)              # 谱图特征摘要
  -> prepare_search(organism="human")             # 推荐参数 + 解析 FASTA
  -> run_search(...)                              # 立即返回 run_id（后台搜索）

你：好了吗？
  -> get_search_status(run_id)                    # 阶段 + 进度%（必要时 diagnose_search）

你：给我结果。
  -> generate_summary(run_id)                      # 1% FDR 统计摘要
  -> infer_proteins(run_id)                        # parsimony + razor + 蛋白级 FDR
  -> export_results(run_id)                        # psm/peptide/protein.tsv + json

你：把扫描 12345 的肽段画出来。
  -> annotate_spectrum(run_id, scan_number=12345)  # b/y 离子标注 HTML（绝对路径返回）
  -> extract_xic(run_id, scan_number=12345)        # 碎片离子 XIC（支持 SILAC）
```

生成的 HTML/TSV 默认写入 `./output/`，可用 `PROTEIN_OUTPUT_DIR` 改到任意目录；返回的路径均为绝对路径。

### 环境变量

| 变量 | 作用 | 默认 |
|------|------|------|
| `RUST_LOG` | 日志级别（写 stderr） | `info` |
| `PROTEIN_LOG_JSON` | 设 `1` 时日志输出 JSON | 关 |
| `PROTEIN_OUTPUT_DIR` | 生成文件的基目录 | `./output` |

## 快速测试

```bash
# 读取谱图文件
cargo run -p protein-copilot-spectrum-io --example read_spectra -- <file.mgf|mzML>

# 完整搜索流程（谱图 → 参数推荐 → 搜索 → 报告导出）
cargo run --release -p protein-copilot-search-engine --example full_search -- \
  <spectrum.mgf|mzML> <database.fasta> [output_dir]

# 陷阱库分析（DIA-NN parquet → L0-L4 分级 → HTML 报告）
cargo run --release -p protein-copilot-entrapment-cli -- analyze \
  --results <report.parquet> --config <config.yaml> --target-fasta <human.fasta>
```

## 项目结构

```text
crates/
├── core/                共享领域模型（Spectrum, SearchParams, SearchResult, ProteinGroup, SearchDiagnostics 等）
├── spectrum-io/         谱图文件解析（mgf/mzML streaming + indexed 随机访问）
├── param-recommend/     参数推荐规则引擎（确定性，不调 LLM）
├── search-engine/       搜索引擎（SimpleSearch + Sage adapter + pFind 预留）
├── dia-extraction/      DIA 前体离子提取（同位素模式检测 + MS1↔MS2 关联）
├── fdr/                 FDR 计算（PSM/肽段/蛋白 三级 + decoy 生成 + picked-protein）
├── protein-inference/   蛋白推断（mapper + parsimony + razor + 序列覆盖率）
├── xic/                 XIC 碎片离子色谱图提取与 Plotly.js HTML 可视化
├── result-import/       外部搜索结果导入（DIA-NN parquet / custom JSON / UnimodDb）
├── fasta-db/            FASTA 数据库管理（内置注册表 + HTTPS 下载 + 缓存）
├── report/              报告生成（摘要 + TSV/JSON 导出）
├── integration-tests/   集成测试（端到端流水线验证）
├── entrapment-analysis/ 陷阱库分析（L0-L4 分级 + 配置 + 加载器 + 报告）
├── entrapment-cli/      陷阱库分析 CLI 工具
└── mcp-server/          MCP Server 二进制（28 tools，stdio transport）

.github/
├── agents/proteomics-search.agent.md     蛋白搜索助手 Agent（25 tools 完整工作流）
├── prompts/basic-search.prompt.md        基础搜索 Skill
├── prompts/failure-diagnosis.prompt.md   搜索失败诊断 Skill
├── prompts/sage-search.prompt.md         Sage 引擎搜索 Skill
├── prompts/protein-inference.prompt.md   蛋白推断 Skill
├── prompts/database-management.prompt.md FASTA 管理 Skill
├── prompts/batch-search.prompt.md        批处理搜索 Skill
└── prompts/result-interpretation.prompt.md  结果解读 Skill
    (+5 more prompts: dia, spectrum-annotation, prd-creation, task-*)
```

## MCP Tools（27 个）

> 下表为速查摘要。**完整接口契约**（每个工具的输入参数、必填/默认、类型定义、输出结构）见自动生成的 [`docs/mcp-tools.md`](docs/mcp-tools.md)（由 `scripts/gen_mcp_tools_doc.py` 从二进制 `tools/list` 生成）。也可在终端直接查看：`protein-copilot-mcp --list-tools`（文本目录）或 `--list-tools --json`（完整 JSON Schema）。

| Tool | 功能 |
|------|------|
| `read_spectra` | 读取谱图文件 → 统计摘要 |
| `get_spectrum` | 按 scan 读取单张谱图 |
| `recommend_params` | 推荐搜索参数 + 解释 |
| `list_presets` | 列出内置预设 |
| `prepare_search` | 组合操作：参数推荐 + 验证 + 准备（recommend→search 桥接） |
| `run_search` | 异步执行数据库搜索（立即返回 run_id） |
| `get_search_status` | 查询搜索进度（阶段 + 百分比 + 诊断标记） |
| `cancel_search` | 取消正在运行的搜索 |
| `diagnose_search` | 搜索诊断报告（阶段耗时 + 异常检测 + 修复建议） |
| `check_engine` | 检查引擎状态 |
| `generate_summary` | FDR 过滤统计摘要 |
| `export_results` | 导出 TSV/JSON 文件 |
| `list_searches` | 列出搜索历史（活跃 + 持久化） |
| `annotate_spectrum` | 谱图注释（DIA: 标注+XIC+SILAC 统一视图；DDA: 标注 only） |
| `extract_dia_precursors` | DIA MS1 前体离子提取（同位素模式检测） |
| `extract_spectrum_precursors` | 单张 MS2 谱图母离子提取（调试用） |
| `get_dia_cache_status` | DIA 提取缓存状态（内存/磁盘溢出统计） |
| `extract_xic` | 碎片离子 XIC 色谱图（支持 SILAC 轻重标记） |
| `import_search_results` | 导入外部搜索结果（DIA-NN / custom JSON） |
| `infer_proteins` | 蛋白推断（parsimony + razor + 蛋白级 FDR + 序列覆盖率） |
| `list_databases` | 列出内置 FASTA 数据库（UniProt 物种库） |
| `download_database` | 下载 FASTA 数据库到本地缓存 |
| `get_database_info` | 查询已缓存数据库的详细信息 |
| `classify_entrapment_hits` | 运行陷阱库分类流程（L0-L4 分级 + HTML 报告） |
| `analyze_entrapment_stats` | 从已分类 TSV 生成统计分析 |
| `find_similar_targets` | 查找肽段在 target 库中的最相似序列 |
| `annotate_provenance` | 标注命中溯源（谱图证据 + 多靶 SILAC 链路） |

## 架构原则

- **确定性/LLM 分层**：Rust 做所有计算，LLM 做意图理解和结果解释
- **MCP 协议**：所有能力通过 MCP tools 暴露给 LLM
- **三级 FDR**：PSM → 肽段 → 蛋白质，各级独立 FDR 控制
- **DDA + DIA 支持**：自动检测采集模式，DIA 数据通过 MS1 同位素模式提取前体离子后搜索
- **外部结果导入**：DIA-NN parquet / 自定义 JSON → RT 匹配 mzML 扫描号 → 标准 SearchResult
- **搜索诊断**：结构化错误分类 + 7 条异常检测规则 + 分级修复建议（确定性，不依赖 LLM）
- **可测试**：795 个单元/集成测试，0 clippy warnings
- **可审计**：每次搜索生成 run_id + 完整参数 + 引擎版本 + 诊断报告

## 当前进度

| 里程碑 | 状态 |
|--------|------|
| M1.1 core | ✅ 共享类型 + 验证 + trait |
| M1.2 spectrum-io | ✅ mgf/mzML 解析 + indexed 随机访问 |
| M1.3 param-recommend | ✅ 规则引擎 + 5 个预设 |
| M1.4 search-engine | ✅ SimpleSearch + pFind 预留 |
| M1.5 report | ✅ 摘要 + TSV/JSON 导出 |
| M1.6 mcp-server | ✅ 25 MCP tools + Agent + 12 Skill Prompts |
| M1.7 integration | ✅ 端到端测试 + 文档 |
| Post-MVP | ✅ 异步搜索 + 历史持久化 + 谱图注释 + FW-1/2/3/4/6 |
| DIA 支持 | ✅ DIA 前体提取 + 搜索集成 + 端到端工作流 |
| XIC 可视化 | ✅ 碎片离子 XIC + SILAC 轻重标记 + Plotly.js HTML |
| 统一标注+XIC | ✅ 标注+XIC 合并视图 + 客户端 SILAC + 逐离子 L/H 开关 |
| 外部结果导入 | ✅ DIA-NN parquet + custom JSON + RT 扫描匹配 + UnimodDb |
| Biology Audit | ✅ 全部审计，单位统一，score 方向规范化 |
| FASTA 管理 | ✅ 内置 UniProt 注册表 + HTTPS 下载 + 本地缓存 |
| **蛋白推断** | ✅ **parsimony + razor + 三级 FDR + 序列覆盖率 + MCP tool** |
| **Sage 集成** | ✅ **sage-core v0.15.0 库集成 + rayon 并行 + LDA rescoring** |
| **工作流优化** | ✅ **prepare_search 桥接 + DIA 缓存溢出 + Agent 工作流更新** |
| **搜索诊断** | ✅ **错误分类 + 阶段指标 + 7 条异常检测 + 修复建议 + diagnose_search tool** |
| **RT 二分查找** | ✅ **ScanIndex + PCIX v2 缓存 + O(log N) find_by_rt + collect_ms2_info 零 I/O** |
| **陷阱库分析** | ✅ **L0-L4 同源性分级 + DIA-NN parquet 加载 + FASTA 酶切索引 + HTML 报告 + CLI** |
| **陷阱库 v2** | ✅ **Levenshtein edit distance + k-mer 预筛 + SubstitutionType 注释 + 3 新输出列 + mDa 显示** |

详细计划：`tasks/001-mvp-proteomics-search-platform.md`
Phase 2 计划：`tasks/002-phase2-production-platform.md`
陷阱库分析：`docs/entrapment-analysis.md`
架构设计：`docs/architecture.md`
架构演示：`docs/architecture.html`
MCP 工具接口契约：[`docs/mcp-tools.md`](docs/mcp-tools.md)（自动生成）
分层文档 L1-L4：[`docs/levels/`](docs/levels/README.md)

## 大文件性能优化

处理大型 mzML 文件（>1GB）时，ProteinCopilot 使用三层索引加速：

1. **PCIX v2 磁盘缓存**（`.mzml.idx`）— 首次打开后自动生成，后续毫秒级加载（46B/entry，含 RT + ms_level + 隔离窗口）
2. **SIMD 字节扫描**（首次构建）— 使用 `memchr` 加速全文件扫描，提取完整元数据
3. **O(log N) RT 二分查找**（`find_by_rt()`）— 按保留时间 + 前体 m/z 定位 MS2 scan，用于谱图标注、XIC 提取、SILAC 重标谱匹配

### 性能数据（7.5GB mzML，SSD）

| 操作 | 耗时 |
|------|------|
| PCIX v2 缓存加载 | <1 ms |
| RT 二分查找（单次） | ~5 ms |
| 字节扫描首次构建 | ~5 s |
| collect_ms2_info（从索引） | <1 ms |

### MCP 超时配置

对于 8GB+ 的大文件，建议在 `.mcp.json` 中增加超时时间：

```json
{
  "mcpServers": {
    "protein-copilot": {
      "timeout": 300
    }
  }
}
```

默认超时 60 秒可能不足以完成首次索引构建。设置 `timeout: 300`（5 分钟）可避免超时错误。

## License

MIT
