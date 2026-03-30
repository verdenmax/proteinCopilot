# 配置指南

## 当前状态（MVP）

MVP 阶段使用内置的 **SimpleSearchEngine**，无需外部配置即可运行。

### MCP Server 启动

```bash
# 直接启动（使用 .mcp.json 自动发现）
cargo run --release -p protein-copilot-mcp-server

# 指定日志级别
RUST_LOG=info cargo run --release -p protein-copilot-mcp-server
```

### 环境变量

| 变量 | 说明 | 默认值 |
|------|------|--------|
| `RUST_LOG` | 日志级别（error/warn/info/debug/trace） | 无（不输出日志） |

### .mcp.json

项目根目录的 `.mcp.json` 用于 Copilot CLI 自动发现 MCP Server：

```json
{
  "mcpServers": {
    "protein-copilot": {
      "command": "cargo",
      "args": ["run", "--release", "-p", "protein-copilot-mcp-server"],
      "env": { "RUST_LOG": "info" }
    }
  }
}
```

---

## 未来配置（pFind 接入后）

接入 pFind 搜索引擎后，需要额外的 SSH 和引擎配置。以下是预留的配置结构：

### SSH 配置（SshConfig）

```json
{
  "host": "compute-01.lab.edu",
  "port": 22,
  "user": "researcher",
  "key_path": "/home/user/.ssh/id_rsa",
  "pfind_executable": "/opt/pfind/bin/pfind",
  "work_dir": "/tmp/pfind_work"
}
```

### 计划的 AppConfig 结构

```rust
struct AppConfig {
    ssh_config: SshConfig,        // SSH 连接配置
    pfind_config: PFindConfig,    // pFind 引擎配置
    data_dirs: Vec<PathBuf>,      // 数据目录列表
}
```

配置文件将通过 `--config` 命令行参数或 `PROTEIN_COPILOT_CONFIG` 环境变量加载。

> **注意**：AppConfig 和 SSH 配置将在 pFind adapter 实现时一并完成。参见 tasks/001 M1.4 "推迟" 部分。
