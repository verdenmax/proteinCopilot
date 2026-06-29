# L4 — fasta-db crate

回溯 [L2](L2-architecture.md)。本篇只聚焦 `crates/fasta-db` 一个 crate 的函数级 API、常量与骨架，签名/常量/内置 DB ID 均按源码核验，完整逻辑以源码为准。

## 1. 用途 + 位置 + 依赖

`protein-copilot-fasta-db` 负责 FASTA 蛋白库的"注册表 + 下载 + 缓存"三件事：内置一份常见蛋白组学库的静态清单（UniProt Swiss-Prot 主流模式生物 + cRAP 污染库），通过 HTTPS 流式下载到本地缓存目录，并用 `registry.json` 追踪元数据。

无任何内部（workspace）依赖；外部仅 `reqwest`（HTTPS）、`sha2`（SHA256）、`tokio`、`chrono`、`serde`/`serde_json`、`schemars`、`thiserror`、`tracing`。唯一上游是 mcp-server：

```
mcp-server tools: list_databases | download_database | get_database_info | prepare_search
   -> protein-copilot-fasta-db (lib API)
      -> registry (内置清单, 只读静态)
      +  cache    (CacheManager + registry.json 原子写)
      +  download (reqwest 流式 + SHA256 + protein 计数)
```

模块导出见 `lib.rs:7-14`：`cache` / `download` / `error` / `registry`，并 re-export `CacheManager` / `CacheRegistry` / `CachedDatabase` / `FastaDbError` / `DatabaseEntry`。整个 crate 不含 LLM、不做打分，纯粹是"取库 + 落盘 + 记账"的确定性基础设施，全部失败路径都走 `Result<_, FastaDbError>`，库代码无 `unwrap/expect`。

## 2. 对外 API + 内置库清单

| 函数 / 类型 | 签名（源码原样） | 位置 |
|---|---|---|
| `list_databases` | `fn list_databases(cache_dir: &Path) -> Result<Vec<DatabaseStatus>, FastaDbError>` | lib.rs:66 |
| `download_database` | `async fn download_database(database_id: &str, cache_dir: &Path, force: bool) -> Result<DownloadDatabaseResult, FastaDbError>` | lib.rs:102 |
| `get_database_info` | `fn get_database_info(database_id: &str, cache_dir: &Path) -> Result<DatabaseInfo, FastaDbError>` | lib.rs:155 |
| `registry::all_databases` | `fn all_databases() -> &'static [DatabaseEntry]` | registry.rs:66 |
| `registry::all_database_ids` | `fn all_database_ids() -> Vec<&'static str>` | registry.rs:71 |
| `registry::get_database` | `fn get_database(id: &str) -> Option<&'static DatabaseEntry>` | registry.rs:76 |
| `download::download_fasta` | `async fn download_fasta(url: &str, dest_path: &Path) -> Result<DownloadResult, FastaDbError>` | download.rs:84 |
| `CacheManager::save_entry` | `fn save_entry(&self, entry: &CachedDatabase) -> Result<(), FastaDbError>` | cache.rs:104 |
| `CacheManager::load_registry` | `fn load_registry(&self) -> Result<CacheRegistry, FastaDbError>` | cache.rs:75 |
| `DatabaseStatus` | `struct { id, species, db_type, description, status: DownloadStatus }` | lib.rs:21 |
| `DownloadStatus` | `enum #[serde(tag="state")] { Available, Downloaded { file_name, file_size_bytes, protein_count, downloaded_at } }` | lib.rs:32 |
| `DatabaseInfo` | `struct { id, species, db_type, path, protein_count, file_size_bytes, downloaded_at, sha256, first_accessions: Vec<String> }` | lib.rs:44 |
| `DownloadDatabaseResult` | `struct { id, path, protein_count, file_size_bytes }` | lib.rs:58 |
| `DatabaseEntry` | `struct { id, species, taxonomy_id: u32, db_type, description, url: &'static str }` | registry.rs:5 |
| `FastaDbError` | `enum { NetworkError, IoError, RegistryError, UnknownDatabase, DownloadFailed }` | error.rs:6 |

内置清单 `BUILTIN_DATABASES`（registry.rs:14-63，共 6 条，URL 取自 UniProt REST / theGPM）：

| id | species | taxonomy_id | db_type | URL 来源 |
|---|---|---|---|---|
| `human_swissprot` | Homo sapiens | 9606 | Swiss-Prot | rest.uniprot.org (reviewed+9606) |
| `mouse_swissprot` | Mus musculus | 10090 | Swiss-Prot | rest.uniprot.org (reviewed+10090) |
| `ecoli_swissprot` | Escherichia coli (K12) | 83333 | Swiss-Prot | rest.uniprot.org (reviewed+83333) |
| `yeast_swissprot` | Saccharomyces cerevisiae | 559292 | Swiss-Prot | rest.uniprot.org (reviewed+559292) |
| `arabidopsis_swissprot` | Arabidopsis thaliana | 3702 | Swiss-Prot | rest.uniprot.org (reviewed+3702) |
| `crap` | Contaminants | 0 | cRAP | ftp.thegpm.org/fasta/cRAP/crap.fasta |

未知 id 返回 `UnknownDatabase { id, available }`，`available` 即 `all_database_ids()`，错误信息直接列出全部合法 id。缓存目录内的文件布局固定：每个库一份 `{cache_dir}/{id}.fasta`（`CacheManager::fasta_path`，cache.rs:59），共享一份 `{cache_dir}/registry.json`（`registry_path`，cache.rs:54）；`get_database_info` 的 `first_accessions` 由 `read_first_accessions` 读 FASTA 前 5 行 `>` 头、按空白切出首个 token 得到（lib.rs:193-218）。

## 3. 安全 / 健壮性

- 下载超时（download.rs:18-23）：`CONNECT_TIMEOUT = 30s`（TCP/TLS 建连），`READ_TIMEOUT = 60s`（两次收字节之间的空闲上限，per-read 而非总时长，大库仍能下完）。`build_download_client` 把二者注入 `reqwest::Client`，慢服务器不会挂死。
- 流式 + 大小上限：响应体逐 chunk 消费，绝不整体进内存；累计超 `MAX_FASTA_DOWNLOAD_BYTES = 10 * 1024^3`（10 GiB，download.rs:15）即清理 `.part` 并报 `DownloadFailed`，防恶意/错配服务器撑爆磁盘。protein 计数与 SHA256 都在流上增量完成（`StreamingProteinCounter` 跨 chunk 记 `prev_byte`，保证行首 `>` 即便落在 chunk 边界也不漏计）。
- `.part` 原子落盘：先写 `dest.fasta.part`，`flush` 后 `tokio::fs::rename` 原子改名到目标；任一步出错（网络/IO/超限/空响应）都 best-effort `remove_file(.part)`，绝不留半截文件冒充成品。空响应（`total_bytes == 0`）也判 `DownloadFailed`。
- 注册表原子写（cache.rs:104-158）：`save_entry` 读-改-写，写入带 `pid + 单调序号` 的唯一 tmp 文件（`REGISTRY_TMP_SEQ` 防并发撞名），`write_all -> flush -> sync_all -> rename`，同盘 rename 原子，读者绝不会看到半截 JSON；失败清理 tmp。
- 损坏自愈（cache.rs:75-96）：`load_registry` 遇 JSON 解析失败只 `warn` 并返回空注册表（崩溃截断可自恢复，磁盘上的 FASTA 仍在、重存即重建索引）；只有读 IO 错误才是硬错误。
- SHA256：下载即算 64 位十六进制摘要存入 `CachedDatabase.sha256`，`get_database_info` 回显，供完整性核对。

## 4. 简化源码片段

带超时的下载客户端（download.rs:27-36）：

```rust
fn build_download_client() -> Result<reqwest::Client, FastaDbError> {
    reqwest::Client::builder()
        .connect_timeout(CONNECT_TIMEOUT)   // 30s 建连
        .read_timeout(READ_TIMEOUT)         // 60s 空闲, per-read
        .build()
        .map_err(|e| FastaDbError::NetworkError { url: String::new(), detail: format!("...: {e}") })
}
```

流式下载 + 大小上限 + 原子改名（download.rs:129-191，节选）：

```rust
loop {
    let chunk = match response.chunk().await {
        Ok(Some(c)) => c,
        Ok(None)    => break,
        Err(e)      => { let _ = remove_file(&part_path).await; return Err(NetworkError{..}); }
    };
    total_bytes += chunk.len() as u64;
    if total_bytes > MAX_FASTA_DOWNLOAD_BYTES {           // 10 GiB 硬上限
        let _ = remove_file(&part_path).await;
        return Err(DownloadFailed { id: url.into(), detail: "exceeds maximum ...".into() });
    }
    file.write_all(&chunk).await?;                        // 出错同样清 .part
    hasher.update(&chunk);                                // SHA256 增量
    counter.update(&chunk);                               // protein 计数增量
}
// ... flush + 空响应检查 ...
tokio::fs::rename(&part_path, dest_path).await?;          // 原子改名为成品
```

注册表原子写（cache.rs:118-130 + write_atomic 135-158，节选）：

```rust
let tmp_path = self.cache_dir.join(format!(
    "registry.json.{}.{}.tmp",
    std::process::id(),
    REGISTRY_TMP_SEQ.fetch_add(1, Ordering::Relaxed),   // 并发不撞名
));
// write_atomic: create -> write_all -> flush -> sync_all -> rename(tmp, dest)
if let Err(e) = Self::write_atomic(&tmp_path, &self.registry_path(), json.as_bytes()) {
    let _ = std::fs::remove_file(&tmp_path);             // 失败清 tmp
    return Err(e);
}
```

## 5. 调用链

- `list_databases` 工具（tools.rs:3765-3778）-> `fasta_db::list_databases(&cache_dir)`：合并内置清单与 `registry.json`，每条标 `Available` 或 `Downloaded{..}`（须 registry 有项且 `.fasta` 实际存在）。
- `download_database` 工具（tools.rs:3785-3801）-> `fasta_db::download_database(id, &cache_dir, force)`，`force = input.force.unwrap_or(false)`；`force=false` 且已缓存则直接回缓存路径，否则 `download_fasta` 后 `save_entry`。
- `get_database_info` 工具（tools.rs:3808-3819）-> `fasta_db::get_database_info(id, &cache_dir)`：未下载报 `DownloadFailed` 提示先 download，已下载则读前 5 条 accession（`read_first_accessions`）+ SHA256 回显。
- `prepare_search`（tools.rs:3620 起）：`organism_to_database_id` 把 `human/mouse/E.coli/yeast/arabidopsis`（含中文 人/小鼠/大肠杆菌）模糊映射成 id，再 `list_databases` 判已缓存、`download_database(.., false)` 取路径写入 `params.database_path`。
- 缓存目录由 mcp-server 的 `default_cache_dir` 决定：入参覆盖，否则默认 `.proteincopilot/databases`（tools.rs:863-869）。三个工具均为只读或幂等：重复 `list/info` 无副作用，`download` 在 `force=false` 时命中缓存即 O(1) 返回路径，整体可安全重试。

## 6. 测试入口

```
cargo test -p protein-copilot-fasta-db --offline
```

16 个单元测试全过（cache 5 + download 6 + registry 3 + lib 2），覆盖：流式计数/SHA256 跨 chunk 边界与 one-shot 一致、客户端超时可构建、注册表往返与"不留 tmp"、损坏自愈、`is_cached` 状态、`list_databases` 全 `Available`、未缓存 `get_database_info` 报错。`tests/integration.rs` 的 2 个端到端用例（下载 E.coli、未知库报错）标 `#[ignore = "requires network access"]`，`--offline` 下计为 ignored，需真实网络时另跑 `-- --ignored`。

---

回到 [README](README.md) 选择其它层级或子系统。
