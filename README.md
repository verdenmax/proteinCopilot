# ProteinCopilot

AI 驱动的蛋白质组学质谱搜索与结果解释平台。

## 架构

- **Rust Workspace + 多 MCP Server**：确定性计算（谱图解析、搜索引擎调度、FDR 计算等）
- **Agent.md + Skill Prompt**：LLM 编排层（意图理解、参数推荐、结果解释）
- **MCP Client（Copilot CLI / Claude Desktop）**：用户交互层

```text
用户 <─> MCP Client (Copilot CLI)
              │
              ├── agents / prompts    ← AI 编排
              └── MCP Servers (Rust)  ← 确定性计算
                    ├── mcp-spectrum-io
                    ├── mcp-param-recommend
                    ├── mcp-search-engine (pFind)
                    ├── mcp-report
                    └── core
```

## 状态

🚧 项目初始化阶段 — 参见 `tasks/001-mvp-proteomics-search-platform.md`

## License

TBD
