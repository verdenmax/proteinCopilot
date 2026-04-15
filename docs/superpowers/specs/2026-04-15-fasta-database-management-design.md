# FASTA Database Management — Design Spec

> **Date**: 2026-04-15
> **Status**: Approved
> **Implements**: FR-6.2 (内置常用数据库列表) + FR-6.3 (自动下载并缓存)

## Problem

当前 FASTA 数据库管理仅支持手动指定本地路径（FR-6.1）。用户必须自己从 UniProt 下载 FASTA 文件，知道正确的 URL 和格式。需要：

1. 内置常用物种的 Swiss-Prot 数据库注册表
2. 一键下载并缓存到本地
3. Agent 可根据用户意图自动建议合适的数据库

## Architecture

### New Crate: `crates/fasta-db`

独立 lib crate，职责：
- 管理内置数据库注册表（物种 → UniProt URL 映射）
- 下载 FASTA 文件到本地缓存
- 跟踪已下载数据库的元数据（蛋白数量、文件大小、下载时间）
- 提供 cRAP 污染物库

**不包含**：decoy 生成（由搜索引擎负责）、FASTA 解析（已在 search-engine/fasta.rs）。

### Dependencies

```
fasta-db
├── reqwest (HTTPS 下载, with rustls)
├── serde / serde_json (registry 序列化)
├── tokio (async I/O)
├── tracing (日志)
└── protein-copilot-core (错误类型复用)
```

## Built-in Database Registry

硬编码在 Rust 常量中，不从远程加载。

| ID | 物种 | Taxonomy ID | 类型 | 预估大小 |
|----|------|-------------|------|----------|
| `human_swissprot` | Homo sapiens | 9606 | Swiss-Prot | ~12 MB |
| `mouse_swissprot` | Mus musculus | 10090 | Swiss-Prot | ~9 MB |
| `ecoli_swissprot` | Escherichia coli (K12) | 83333 | Swiss-Prot | ~2 MB |
| `yeast_swissprot` | Saccharomyces cerevisiae | 559292 | Swiss-Prot | ~3 MB |
| `arabidopsis_swissprot` | Arabidopsis thaliana | 3702 | Swiss-Prot | ~6 MB |
| `crap` | Contaminants (cRAP) | — | GPM cRAP | ~50 KB |

### Download URL Pattern

UniProt REST API (HTTPS):
```
https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(reviewed:true)+AND+(model_organism:{taxonomy_id})
```

cRAP: 从 GPM 官方 URL 或内嵌静态文件。

## Cache Directory

位置：项目根目录下 `.proteincopilot/databases/`

```
.proteincopilot/databases/
├── registry.json              ← 已下载库的元数据
├── human_swissprot.fasta      ← 下载的 FASTA 文件
├── mouse_swissprot.fasta
├── crap.fasta
└── ...
```

### registry.json Schema

```json
{
  "version": 1,
  "databases": {
    "human_swissprot": {
      "file_name": "human_swissprot.fasta",
      "downloaded_at": "2026-04-15T12:00:00Z",
      "file_size_bytes": 12345678,
      "protein_count": 20422,
      "sha256": "abc123..."
    }
  }
}
```

## Data Structures

### DatabaseEntry (内置注册表项)

```rust
pub struct DatabaseEntry {
    pub id: &'static str,           // "human_swissprot"
    pub species: &'static str,      // "Homo sapiens"
    pub taxonomy_id: u32,           // 9606
    pub db_type: &'static str,      // "Swiss-Prot"
    pub description: &'static str,  // "Human reviewed proteome"
    pub url: &'static str,          // UniProt REST URL
}
```

### CachedDatabase (已下载的元数据)

```rust
pub struct CachedDatabase {
    pub id: String,
    pub file_name: String,
    pub downloaded_at: DateTime<Utc>,
    pub file_size_bytes: u64,
    pub protein_count: u64,
    pub sha256: String,
}
```

### DatabaseStatus (列表返回)

```rust
pub struct DatabaseStatus {
    pub id: String,
    pub species: String,
    pub db_type: String,
    pub description: String,
    pub status: DownloadStatus,  // Available | Downloaded { info } | Outdated
}

pub enum DownloadStatus {
    Available,                     // 未下载，可下载
    Downloaded { cached: CachedDatabase },  // 已下载
}
```

## MCP Tools

### 1. `list_databases`

**输入**: 无必选参数，可选 `cache_dir` 覆盖默认路径。

**输出**: `Vec<DatabaseStatus>` — 所有内置库 + 缓存状态。

**逻辑**:
1. 加载内置注册表（6 个条目）
2. 读取 `registry.json`（如存在）
3. 合并：对每个内置库，检查是否已缓存
4. 返回完整列表

### 2. `download_database`

**输入**:
- `database_id: String` — 必选，如 "human_swissprot"
- `cache_dir: Option<String>` — 可选覆盖
- `force: Option<bool>` — 是否强制重新下载

**输出**: `DownloadResult { path, protein_count, file_size_bytes, duration_sec }`

**逻辑**:
1. 查找内置注册表中的 URL
2. 如已缓存且非 force，返回已有路径
3. 创建缓存目录
4. 通过 reqwest 流式下载到临时文件
5. 计算 SHA256、统计蛋白数量
6. 原子重命名到目标路径
7. 更新 registry.json
8. 返回结果

**错误处理**:
- 网络失败 → 保留已有缓存，返回错误
- 无效 database_id → 返回可用 ID 列表
- 磁盘空间不足 → 清理临时文件，返回错误

### 3. `get_database_info`

**输入**: `database_id: String`

**输出**: `DatabaseInfo { id, species, path, protein_count, file_size, downloaded_at, sha256, first_5_accessions }`

**逻辑**:
1. 读取 registry.json
2. 如未下载，返回错误（建议先 download）
3. 读取 FASTA 文件前 5 个 entry 的 accession（供用户确认正确性）
4. 返回详细信息

## Agent Integration

在 `proteomics-search.agent.md` 和 `basic-search.prompt.md` 中补充指导：

- 当用户提到物种名（如"人"、"小鼠"、"大肠杆菌"）时，先调用 `list_databases` 检查缓存
- 如未下载，建议用户下载并等待确认
- 下载完成后，将路径设为 `SearchParams.database_path`
- 如用户未指定物种，提示选择

## cRAP Contaminant Database

cRAP 库特殊处理：
- 可内嵌为 Rust 静态字符串（~50 KB，可接受）
- 或从 GPM 下载：`https://www.thegpm.org/crap/crap.fasta`
- `list_databases` 中标记为独立条目
- 用户可选择与物种库合并或单独使用

## Error Types

```rust
pub enum FastaDbError {
    NetworkError { url: String, detail: String },
    IoError { path: PathBuf, source: std::io::Error },
    RegistryError { detail: String },
    UnknownDatabase { id: String, available: Vec<String> },
    DownloadFailed { id: String, detail: String },
}
```

## Testing Strategy

1. **单元测试**: 注册表查找、registry.json 序列化/反序列化、蛋白计数
2. **集成测试**: 真实 HTTP 下载小库（ecoli ~2MB），验证 FASTA 有效性
3. **Mock 测试**: 使用 mockito 模拟 HTTP 响应，测试错误处理
4. **缓存测试**: 验证重复下载跳过、force 覆盖、registry.json 原子更新

## Non-Goals

- 不做 decoy 拼接（搜索引擎负责）
- 不做定期自动更新（用户手动 force 刷新）
- 不做 TrEMBL 支持（文件太大，Swiss-Prot 足够 MVP）
- 不做自定义 URL 注册（用户可直接指定本地路径）
