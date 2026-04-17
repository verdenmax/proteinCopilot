# Workflow Optimization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix 3 workflow break-points (recommend→search disconnect, DIA cache eviction, database_path pre-validation), add 2 new MCP tools (`prepare_search`, `get_dia_cache_status`), update the proteomics-search agent with 5 missing tools + 5 workflow sections, and create 4 new prompt files.

**Architecture:** New `prepare_search` composite tool in `tools.rs` combines spectrum reading + param recommendation + database auto-resolution. DIA cache gets disk spillover via bincode serialization. Agent and prompts updated to cover all 23 tools (21 existing + 2 new).

**Tech Stack:** Rust, rmcp (MCP framework), serde/bincode, tokio, protein_copilot_* crates

---

## File Structure

| Action | File | Responsibility |
|--------|------|----------------|
| Modify | `crates/mcp-server/src/tools.rs` | Add `prepare_search` tool, `get_dia_cache_status` tool, database_path pre-validation, DIA disk spillover |
| Modify | `crates/mcp-server/Cargo.toml` | Add `bincode` dependency for DIA cache serialization |
| Modify | `.github/agents/proteomics-search.agent.md` | Add 7 tools to declaration + 5 workflow sections |
| Create | `.github/prompts/protein-inference.prompt.md` | Protein inference workflow |
| Create | `.github/prompts/database-management.prompt.md` | FASTA database management workflow |
| Create | `.github/prompts/sage-search.prompt.md` | Sage engine specific guidance |
| Create | `.github/prompts/batch-search.prompt.md` | Multi-file batch search workflow |

---

### Task 1: Database path pre-validation in `run_search`

**Files:**
- Modify: `crates/mcp-server/src/tools.rs:1268-1274` (file-based path)
- Modify: `crates/mcp-server/src/tools.rs:1019-1025` (DIA path)

- [ ] **Step 1: Add database_path existence check to file-based path**

In `crates/mcp-server/src/tools.rs`, after line 1274 (`params.validate()`), insert database file existence check:

```rust
// After: params.validate().map_err(...)?;
// Insert before: let run_id = Uuid::new_v4();

// Validate database file exists before spawning background task
let db_path = Path::new(&params.database_path);
if !db_path.exists() {
    return Err(mcp_err(
        ErrorCode::INVALID_PARAMS,
        format!(
            "Database file not found: {}. Use list_databases to see available databases, \
             or download_database to fetch one.",
            params.database_path
        ),
    ));
}
```

- [ ] **Step 2: Add database_path existence check to DIA path**

In `crates/mcp-server/src/tools.rs`, after line 1025 (`params.validate()`), insert the same check:

```rust
// After: params.validate().map_err(...)?;
// Insert before: let dia_spectra = { ...

let db_path = Path::new(&params.database_path);
if !db_path.exists() {
    return Err(mcp_err(
        ErrorCode::INVALID_PARAMS,
        format!(
            "Database file not found: {}. Use list_databases to see available databases, \
             or download_database to fetch one.",
            params.database_path
        ),
    ));
}
```

- [ ] **Step 3: Verify build compiles**

Run: `cargo build --workspace 2>&1 | tail -3`
Expected: `Finished` with no errors

- [ ] **Step 4: Run existing tests**

Run: `cargo test --workspace 2>&1 | grep -E "^test result|FAILED"`
Expected: All pass, 0 failed

- [ ] **Step 5: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "fix: validate database_path exists before spawning search task

Previously, an invalid database_path would only be detected inside the
tokio::spawn block, requiring the LLM to poll get_search_status to
discover the error. Now it returns an immediate MCP error with actionable
guidance to use list_databases/download_database.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: DIA cache disk spillover + `get_dia_cache_status` tool

**Files:**
- Modify: `crates/mcp-server/Cargo.toml` (add bincode)
- Modify: `crates/mcp-server/src/tools.rs:634-668` (OrderedDiaCache)
- Modify: `crates/mcp-server/src/tools.rs:751-758` (ProteinCopilotServer struct)
- Modify: `crates/mcp-server/src/tools.rs:2157-2162` (extract_dia_precursors cache insert)
- Modify: `crates/mcp-server/src/tools.rs:1029-1043` (run_search DIA cache lookup)

- [ ] **Step 1: Add bincode dependency**

In `crates/mcp-server/Cargo.toml`, add to `[dependencies]`:

```toml
bincode = "1"
```

- [ ] **Step 2: Add disk spillover to OrderedDiaCache**

Replace the `OrderedDiaCache` implementation (lines 634-668) with:

```rust
struct OrderedDiaCache {
    entries: HashMap<Uuid, Vec<Spectrum>>,
    order: Vec<Uuid>,
    spill_dir: PathBuf,
    extracted_at: HashMap<Uuid, chrono::DateTime<chrono::Utc>>,
}

const MAX_DIA_CACHE_SIZE: usize = 10;

impl OrderedDiaCache {
    fn new() -> Self {
        let spill_dir = PathBuf::from(".proteincopilot/dia_cache");
        Self {
            entries: HashMap::new(),
            order: Vec::new(),
            spill_dir,
            extracted_at: HashMap::new(),
        }
    }

    fn remove(&mut self, id: &Uuid) -> Option<Vec<Spectrum>> {
        // Try memory first
        if let Some(spectra) = self.entries.remove(id) {
            self.order.retain(|x| x != id);
            self.extracted_at.remove(id);
            return Some(spectra);
        }
        // Try disk
        let path = self.spill_dir.join(format!("{}.bin", id));
        if path.exists() {
            match std::fs::read(&path) {
                Ok(data) => {
                    let _ = std::fs::remove_file(&path);
                    self.extracted_at.remove(id);
                    match bincode::deserialize(&data) {
                        Ok(spectra) => return Some(spectra),
                        Err(e) => {
                            tracing::warn!("Failed to deserialize DIA cache {}: {}", id, e);
                            return None;
                        }
                    }
                }
                Err(e) => {
                    tracing::warn!("Failed to read DIA cache file {}: {}", id, e);
                    return None;
                }
            }
        }
        None
    }

    fn insert(&mut self, id: Uuid, spectra: Vec<Spectrum>) {
        // Spill oldest to disk when at capacity
        while self.order.len() >= MAX_DIA_CACHE_SIZE {
            if let Some(oldest) = self.order.first().copied() {
                self.order.remove(0);
                if let Some(old_spectra) = self.entries.remove(&oldest) {
                    self.spill_to_disk(oldest, &old_spectra);
                }
            }
        }
        self.extracted_at.insert(id, chrono::Utc::now());
        self.entries.insert(id, spectra);
        self.order.push(id);
    }

    fn spill_to_disk(&self, id: Uuid, spectra: &[Spectrum]) {
        if let Err(e) = std::fs::create_dir_all(&self.spill_dir) {
            tracing::warn!("Failed to create DIA spill dir: {}", e);
            return;
        }
        let path = self.spill_dir.join(format!("{}.bin", id));
        match bincode::serialize(spectra) {
            Ok(data) => {
                if let Err(e) = std::fs::write(&path, &data) {
                    tracing::warn!("Failed to write DIA cache to disk {}: {}", id, e);
                }
            }
            Err(e) => {
                tracing::warn!("Failed to serialize DIA cache {}: {}", id, e);
            }
        }
    }

    fn status(&self, id: &Uuid) -> DiaCacheLocation {
        if self.entries.contains_key(id) {
            let count = self.entries[id].len();
            let ts = self.extracted_at.get(id).copied();
            DiaCacheLocation::Memory { spectrum_count: count, extracted_at: ts }
        } else {
            let path = self.spill_dir.join(format!("{}.bin", id));
            if path.exists() {
                let ts = self.extracted_at.get(id).copied();
                DiaCacheLocation::Disk { extracted_at: ts }
            } else {
                DiaCacheLocation::NotFound
            }
        }
    }
}

enum DiaCacheLocation {
    Memory { spectrum_count: usize, extracted_at: Option<chrono::DateTime<chrono::Utc>> },
    Disk { extracted_at: Option<chrono::DateTime<chrono::Utc>> },
    NotFound,
}
```

- [ ] **Step 3: Add `get_dia_cache_status` input/output structs**

Add after the existing tool input/output structs (around line 600):

```rust
#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct GetDiaCacheStatusInput {
    /// The dia_run_id returned by extract_dia_precursors
    dia_run_id: String,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct DiaCacheStatusOutput {
    /// Whether the cached extraction exists
    exists: bool,
    /// Where the cache is stored: "memory", "disk", or "not_found"
    location: String,
    /// Number of spectra (only available if in memory)
    spectrum_count: Option<usize>,
    /// When the extraction was performed
    extracted_at: Option<String>,
}
```

- [ ] **Step 4: Add `get_dia_cache_status` tool handler**

Add the tool method in the `#[rmcp::tool_router]` impl block, after `extract_spectrum_precursors`:

```rust
#[rmcp::tool(
    name = "get_dia_cache_status",
    description = "Check if a DIA extraction result is still cached and available for use with run_search. Returns cache location (memory/disk/not_found) and spectrum count. Call this before run_search(dia_run_id=...) to verify the extraction hasn't been evicted."
)]
fn get_dia_cache_status(
    &self,
    Parameters(input): Parameters<GetDiaCacheStatusInput>,
) -> Result<Json<DiaCacheStatusOutput>, ErrorData> {
    let dia_uuid = Uuid::parse_str(&input.dia_run_id)
        .map_err(|_| mcp_err(ErrorCode::INVALID_PARAMS, "invalid dia_run_id format"))?;

    let cache = self
        .dia_cache
        .lock()
        .map_err(|_| mcp_err(ErrorCode::INTERNAL_ERROR, "DIA cache lock is poisoned"))?;

    let output = match cache.status(&dia_uuid) {
        DiaCacheLocation::Memory { spectrum_count, extracted_at } => DiaCacheStatusOutput {
            exists: true,
            location: "memory".to_string(),
            spectrum_count: Some(spectrum_count),
            extracted_at: extracted_at.map(|t| t.to_rfc3339()),
        },
        DiaCacheLocation::Disk { extracted_at } => DiaCacheStatusOutput {
            exists: true,
            location: "disk".to_string(),
            spectrum_count: None,
            extracted_at: extracted_at.map(|t| t.to_rfc3339()),
        },
        DiaCacheLocation::NotFound => DiaCacheStatusOutput {
            exists: false,
            location: "not_found".to_string(),
            spectrum_count: None,
            extracted_at: None,
        },
    };

    Ok(Json(output))
}
```

- [ ] **Step 5: Add chrono and bincode imports**

At the top of `tools.rs`, ensure these imports are present (add if missing):

```rust
use std::path::PathBuf;  // already imported
// bincode is used via bincode::serialize/deserialize — no use statement needed
```

Also add `chrono` to `crates/mcp-server/Cargo.toml` dependencies if not already present (check first — it may already be there).

- [ ] **Step 6: Verify Spectrum implements Serialize + Deserialize**

Check that `protein_copilot_core::spectrum::Spectrum` derives `Serialize` and `Deserialize` (required for bincode). If not, add `#[derive(Serialize, Deserialize)]` to the struct in `crates/core/src/spectrum.rs`.

- [ ] **Step 7: Build and test**

Run: `cargo build --workspace 2>&1 | tail -3`
Expected: `Finished` with no errors

Run: `cargo test --workspace 2>&1 | grep -E "^test result|FAILED"`
Expected: All pass, 0 failed

- [ ] **Step 8: Commit**

```bash
git add crates/mcp-server/Cargo.toml crates/mcp-server/src/tools.rs
git commit -m "feat: DIA cache disk spillover and get_dia_cache_status tool

- OrderedDiaCache now spills evicted entries to disk via bincode
- Entries are restored from disk when accessed via run_search(dia_run_id)
- New get_dia_cache_status tool lets Agent verify cache before searching
- Prevents silent data loss when cache exceeds 10 entries

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: `prepare_search` composite tool

**Files:**
- Modify: `crates/mcp-server/src/tools.rs` (add input/output structs + tool handler)

- [ ] **Step 1: Add organism-to-database-id mapping**

Add helper function in `tools.rs` (near `default_cache_dir` around line 576):

```rust
/// Maps common organism names/keywords to database IDs.
fn organism_to_database_id(organism: &str) -> Option<&'static str> {
    let lower = organism.to_lowercase();
    // Check exact IDs first
    match lower.as_str() {
        "human_swissprot" | "mouse_swissprot" | "ecoli_swissprot"
        | "yeast_swissprot" | "arabidopsis_swissprot" | "crap" => {
            return Some(match lower.as_str() {
                "human_swissprot" => "human_swissprot",
                "mouse_swissprot" => "mouse_swissprot",
                "ecoli_swissprot" => "ecoli_swissprot",
                "yeast_swissprot" => "yeast_swissprot",
                "arabidopsis_swissprot" => "arabidopsis_swissprot",
                "crap" => "crap",
                _ => unreachable!(),
            });
        }
        _ => {}
    }
    // Fuzzy keyword matching
    if lower.contains("human") || lower.contains("人") || lower.contains("homo sapiens") || lower.contains("9606") {
        Some("human_swissprot")
    } else if lower.contains("mouse") || lower.contains("小鼠") || lower.contains("mus musculus") || lower.contains("10090") {
        Some("mouse_swissprot")
    } else if lower.contains("ecoli") || lower.contains("e.coli") || lower.contains("大肠杆菌") || lower.contains("escherichia") {
        Some("ecoli_swissprot")
    } else if lower.contains("yeast") || lower.contains("酵母") || lower.contains("saccharomyces") {
        Some("yeast_swissprot")
    } else if lower.contains("arabidopsis") || lower.contains("拟南芥") {
        Some("arabidopsis_swissprot")
    } else if lower.contains("contaminant") || lower.contains("污染") || lower.contains("crap") {
        Some("crap")
    } else {
        None
    }
}
```

- [ ] **Step 2: Add PrepareSearchInput and PrepareSearchOutput structs**

Add after the existing input/output struct definitions:

```rust
#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct PrepareSearchInput {
    /// Paths to spectrum files (.mgf or .mzML)
    input_files: Vec<String>,
    /// Optional user hints (experiment_type, instrument_type, enzyme)
    #[serde(default, deserialize_with = "deserialize_hints")]
    #[schemars(with = "Option<UserHints>")]
    hints: Option<UserHints>,
    /// Target organism for auto database resolution (e.g. "human", "mouse", "E.coli", "小鼠").
    /// Maps to built-in database IDs. If not provided and database_path is empty,
    /// returns an error asking user to specify.
    organism: Option<String>,
    /// Direct FASTA database path. Takes priority over organism auto-resolution.
    database_path: Option<String>,
    /// Search engine to use: "Sage" or "SimpleSearch". Default: "SimpleSearch".
    engine: Option<String>,
    /// Override cache directory for database downloads. Default: .proteincopilot/databases
    cache_dir: Option<String>,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct PrepareSearchOutput {
    /// Recommended search parameters with real database_path filled in.
    /// Pass this to run_search after user confirms.
    params: SearchParams,
    /// Explanation of why these parameters were recommended.
    reasoning: String,
    /// Confidence score (0.0 to 1.0).
    confidence: f64,
    /// Alternative approaches the user might consider.
    alternatives: Vec<String>,
    /// Evidence supporting the recommendation.
    evidence: Vec<String>,
    /// Summary of the input spectra.
    spectra_summary: SpectrumSummary,
    /// Database info if auto-resolved. None if user provided database_path directly.
    database_info: Option<PreparedDatabaseInfo>,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
struct PreparedDatabaseInfo {
    /// Database ID (e.g. "human_swissprot")
    id: String,
    /// Local file path
    path: String,
    /// Number of protein sequences
    protein_count: u64,
    /// Whether this was freshly downloaded
    freshly_downloaded: bool,
}
```

- [ ] **Step 3: Implement `prepare_search` tool handler**

Add in the `#[rmcp::tool_router]` impl block:

```rust
#[rmcp::tool(
    name = "prepare_search",
    description = "Prepare a proteomics search by combining spectrum analysis, parameter recommendation, and database resolution in one step. Returns recommended SearchParams with a real database_path filled in. Present the result to the user for confirmation, then pass params to run_search. Supports organism-based auto database lookup (e.g. organism='human' auto-downloads UniProt Human Swiss-Prot). If database_path is provided directly, it takes priority over organism."
)]
async fn prepare_search(
    &self,
    Parameters(input): Parameters<PrepareSearchInput>,
) -> Result<Json<PrepareSearchOutput>, ErrorData> {
    // 1. Validate input_files not empty
    if input.input_files.is_empty() {
        return Err(mcp_err(
            ErrorCode::INVALID_PARAMS,
            "input_files is empty — provide at least one spectrum file path",
        ));
    }

    // 2. Validate all input files exist
    for file_str in &input.input_files {
        let p = Path::new(file_str);
        if !p.exists() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!("Input file does not exist: {file_str}"),
            ));
        }
    }

    // 3. Read spectra summary from first file
    let first_file = input.input_files.first().unwrap();
    let path = Path::new(first_file);
    let info = protein_copilot_spectrum_io::detect_format(path)
        .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
    let reader = protein_copilot_spectrum_io::create_reader(&info);
    let summary = reader
        .read_summary(path)
        .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

    // 4. Recommend parameters
    let decision = ParamRecommender
        .recommend(&summary, input.hints.as_ref())
        .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

    let mut params = decision.decision;

    // 5. Set engine if specified
    if let Some(ref engine) = input.engine {
        params.engine = Some(engine.clone());
    }

    // 6. Resolve database path (priority: database_path > organism > error)
    let mut database_info = None;

    if let Some(ref db_path) = input.database_path {
        // User provided direct path
        params.database_path = db_path.clone();
        if !Path::new(db_path).exists() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!(
                    "Database file not found: {}. Check the path or use organism parameter \
                     for auto-download.",
                    db_path
                ),
            ));
        }
    } else if let Some(ref organism) = input.organism {
        // Auto-resolve from organism
        let db_id = organism_to_database_id(organism).ok_or_else(|| {
            mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!(
                    "Cannot map organism '{}' to a database. Supported: human, mouse, \
                     E.coli, yeast, arabidopsis, crap. Or provide database_path directly.",
                    organism
                ),
            )
        })?;

        let cache_dir = default_cache_dir(&input.cache_dir);
        let dbs = protein_copilot_fasta_db::list_databases(&cache_dir)
            .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?;

        let status = dbs.iter().find(|d| d.id == db_id);
        let (db_path, protein_count, freshly_downloaded) = match status {
            Some(protein_copilot_fasta_db::DatabaseStatus {
                status: protein_copilot_fasta_db::DownloadStatus::Downloaded { file_name, protein_count, .. },
                ..
            }) => {
                let full_path = cache_dir.join(file_name).to_string_lossy().to_string();
                (full_path, *protein_count, false)
            }
            _ => {
                // Not cached — download
                let result = protein_copilot_fasta_db::download_database(db_id, &cache_dir, false)
                    .await
                    .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, format!("Database download failed: {e}")))?;
                (result.path.clone(), result.protein_count, true)
            }
        };

        params.database_path = db_path.clone();
        database_info = Some(PreparedDatabaseInfo {
            id: db_id.to_string(),
            path: db_path,
            protein_count,
            freshly_downloaded,
        });
    }
    // If neither database_path nor organism provided, params.database_path stays as default
    // (empty string from recommend_params), and validate() below will catch it

    // 7. Validate final params
    params
        .validate()
        .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

    // 8. Validate database file exists
    let db_path = Path::new(&params.database_path);
    if !db_path.exists() {
        return Err(mcp_err(
            ErrorCode::INVALID_PARAMS,
            format!(
                "Database file not found: {}. Provide database_path or organism parameter.",
                params.database_path
            ),
        ));
    }

    Ok(Json(PrepareSearchOutput {
        params,
        reasoning: decision.explanation,
        confidence: decision.confidence,
        alternatives: decision.alternatives,
        evidence: decision.evidence,
        spectra_summary: summary,
        database_info,
    }))
}
```

- [ ] **Step 4: Build and test**

Run: `cargo build --workspace 2>&1 | tail -3`
Expected: `Finished` with no errors

Run: `cargo test --workspace 2>&1 | grep -E "^test result|FAILED"`
Expected: All pass, 0 failed

Run: `cargo clippy --workspace 2>&1 | tail -3`
Expected: 0 warnings

- [ ] **Step 5: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "feat: add prepare_search composite MCP tool

Combines spectrum reading + parameter recommendation + database
auto-resolution into a single tool call. Supports organism-based
database lookup (e.g. organism='human' auto-downloads UniProt
Human Swiss-Prot). Returns ready-to-use SearchParams for user
confirmation before calling run_search.

Reduces standard workflow from 7 manual steps to 4.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 4: Update proteomics-search.agent.md

**Files:**
- Modify: `.github/agents/proteomics-search.agent.md`

- [ ] **Step 1: Update tools declaration**

Replace the tools list (lines 3-19) with:

```yaml
tools:
  - read_spectra
  - get_spectrum
  - recommend_params
  - list_presets
  - run_search
  - get_search_status
  - cancel_search
  - check_engine
  - generate_summary
  - export_results
  - list_searches
  - annotate_spectrum
  - extract_xic
  - import_search_results
  - extract_dia_precursors
  - extract_spectrum_precursors
  - prepare_search
  - get_dia_cache_status
  - infer_proteins
  - list_databases
  - download_database
  - get_database_info
```

- [ ] **Step 2: Add quick search workflow section**

After the "标准工作流程" section (after line 80), add:

```markdown
## 快速搜索工作流（推荐）

使用 `prepare_search` 复合工具，一步完成参数推荐 + 数据库解析：

### Step 1：准备搜索
- 调用 `prepare_search(input_files=[...], organism="human")`
- 工具自动完成：读取谱图摘要、推荐参数、查找/下载 FASTA 数据库
- 向用户展示推荐参数、置信度、推荐理由

### Step 2：确认参数（必须）
- **等待用户确认或修改**
- 如果用户要求修改，调整返回的 params 对象

### Step 3：执行搜索
- 将确认后的 params 传给 `run_search(params=..., input_files=[...])`
- 后续流程同标准工作流 Step 4.5 - Step 6

### 何时使用快速搜索 vs 标准搜索
- **快速搜索**：适合有明确物种信息的常规搜索
- **标准搜索**：需要精细控制参数、使用自定义数据库、或分步调试时
```

- [ ] **Step 3: Add protein inference workflow section**

After the quick search section, add:

```markdown
## 蛋白推断工作流

搜索完成后，将 PSM 结果聚合到蛋白质水平：

### Step 1：确认搜索完成
- `get_search_status(run_id)` 确认 status = "Completed"

### Step 2：执行蛋白推断
- 调用 `infer_proteins(run_id=xxx, fasta_path=xxx)`
- 可选参数：`fdr_threshold`（默认 0.01）、`min_peptides`（默认 1）

### Step 3：解读结果
- 报告蛋白组数量（protein groups）
- 解释 Parsimony 原则：最小蛋白集覆盖所有肽段
- 区分 unique peptides（仅属于一个蛋白）和 shared peptides（多个蛋白共享）
- Razor 肽段：共享肽段归属到证据最多的蛋白质
- 序列覆盖率：matched peptides 覆盖蛋白序列的百分比

### 领域知识
- 典型 HeLa 样品：3000-6000 蛋白组（取决于分析深度）
- 蛋白 FDR 1% 是标准阈值，发表级别可用 0.1%
- unique peptides ≥ 2 的蛋白鉴定更可靠
```

- [ ] **Step 4: Add database management workflow section**

```markdown
## 数据库管理

### 查看可用数据库
- 调用 `list_databases()` 查看所有内置数据库及缓存状态
- 内置数据库：human_swissprot, mouse_swissprot, ecoli_swissprot, yeast_swissprot, arabidopsis_swissprot, crap

### 下载数据库
- 调用 `download_database(database_id="human_swissprot")` 下载并缓存
- 支持 force=true 强制重新下载
- 下载后可用 `get_database_info(database_id=xxx)` 查看详情（蛋白数量、文件大小、SHA256）

### 自动解析
- 使用 `prepare_search(organism="human")` 时自动处理数据库查找和下载
- 支持中英文物种名：human/人/Homo sapiens, mouse/小鼠, E.coli/大肠杆菌 等

### cRAP 污染物数据库
- `database_id="crap"` 是 Common Repository of Adventitious Proteins
- 包含角蛋白、胰蛋白酶自切等常见污染物
- 建议在正式搜索数据库中包含 cRAP 序列
```

- [ ] **Step 5: Add check_engine usage section**

```markdown
## 搜索引擎管理

### 检查引擎状态
- 调用 `check_engine(engine="Sage")` 确认引擎可用
- 返回引擎名称、版本、健康状态，以及所有已注册引擎列表
- 支持的引擎：Sage（生产级，推荐）、SimpleSearch（内置 MVP）

### 引擎选择指南
- **Sage**：rayon 并行打分、LDA rescoring、三级 FDR（spectrum/peptide/protein），适合生产使用
- **SimpleSearch**：内置简化引擎，适合快速测试和小规模搜索
- 引擎通过 `run_search(params={...engine: "Sage"...})` 或 `prepare_search(engine="Sage")` 指定
```

- [ ] **Step 6: Update DIA workflow with cache check**

In the existing "DIA Data Workflow" section (lines 159-163), update to:

```markdown
### DIA Data Workflow
1. Use `read_spectra` to check if data is DIA (wide isolation windows, median > 5 Da)
2. Call `extract_dia_precursors` to extract candidate precursors from MS1
3. **Call `get_dia_cache_status(dia_run_id=xxx)` to verify cache is available**
4. Use the returned run_id with `run_search(dia_run_id=xxx)` to search the extracted spectra

**注意**：DIA 缓存内存上限为 10 条，超出后自动写入磁盘。使用 `get_dia_cache_status`
确认缓存存在后再调用 `run_search`，避免"not found"错误。

### DIA 检测标准
- **自动检测阈值**：中位隔离窗口宽度 > 5 Da → 判定为 DIA 数据
- 这是启发式阈值，用于自动模式选择
- 仪器级定义中，DDA 使用窄窗口（通常 < 2 Th），DIA 使用宽窗口（通常 10-25 Da）
```

- [ ] **Step 7: Update decision boundary table**

Add new entries to the decision boundary table:

```markdown
| 下载数据库 | ✅ 可自动执行（prepare_search 自动处理） |
| 蛋白推断 | ✅ 搜索完成后可自动执行，但结果需展示给用户 |
| 检查引擎状态 | ✅ 可自动执行 |
| 检查 DIA 缓存 | ✅ 可自动执行 |
```

- [ ] **Step 8: Fix DIA threshold in spectrum-annotation.prompt.md**

In `.github/prompts/spectrum-annotation.prompt.md`, line 37, change:

```markdown
- DIA 检测标准：隔离窗口宽度 > 1 Th 判定为 DIA，≤ 1 Th 为 DDA。
```

to:

```markdown
- DIA 自动检测：中位隔离窗口宽度 > 5 Da 判定为 DIA。仪器级别 DDA 窗口通常 < 2 Th，DIA 窗口通常 10-25 Da。
```

- [ ] **Step 9: Commit**

```bash
git add .github/agents/proteomics-search.agent.md .github/prompts/spectrum-annotation.prompt.md
git commit -m "docs: update proteomics-search agent with all 22 tools and workflows

- Add 6 new tools: prepare_search, get_dia_cache_status, infer_proteins,
  list_databases, download_database, get_database_info
- Add quick search workflow (prepare_search composite tool)
- Add protein inference workflow (Parsimony + Razor)
- Add database management section
- Add engine management section (check_engine usage)
- Update DIA workflow with cache status check
- Fix DIA threshold documentation conflict

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 5: Create protein-inference.prompt.md

**Files:**
- Create: `.github/prompts/protein-inference.prompt.md`

- [ ] **Step 1: Create file**

```markdown
---
mode: agent
description: "蛋白推断工作流 — 从 PSM 搜索结果到蛋白质水平鉴定，支持 Parsimony 最小集和 Razor 肽段"
---

# 蛋白推断

将肽段-谱图匹配（PSM）结果聚合到蛋白质水平，确定样品中存在哪些蛋白质。

## 输入要求
- 已完成搜索的 run_id（从 `run_search` 获得，status = "Completed"）
- FASTA 数据库路径（与搜索使用的相同）

## 流程

1. 确认搜索已完成：`get_search_status(run_id)` → status = "Completed"
2. 调用 `infer_proteins(run_id=xxx, fasta_path=xxx)`
   - 可选参数：
     - `fdr_threshold`：FDR 阈值，默认 0.01（1%）
     - `min_peptides`：最少肽段数，默认 1
3. 解读推断结果：
   - **蛋白组数**：通过 FDR 阈值的蛋白组数量
   - **Unique peptides**：仅属于一个蛋白的肽段，是蛋白鉴定的最强证据
   - **Shared peptides**：被多个蛋白共享的肽段
   - **Razor peptides**：共享肽段中，归属到证据最多蛋白的那一份
   - **序列覆盖率**：匹配肽段覆盖蛋白序列的百分比

## 算法说明

### Parsimony（最小蛋白集）
- 找到能解释所有鉴定肽段的**最少蛋白质数量**
- 一个肽段可能匹配多个蛋白（同源蛋白、亚型）
- Parsimony 消除冗余：如果蛋白 A 的所有肽段都被蛋白 B 包含，则 A 是冗余的

### Razor 肽段分配
- 共享肽段归属到拥有最多 unique peptides 的蛋白质
- 每个共享肽段只计算一次（不重复计数）
- 这是 MaxQuant / Proteome Discoverer 使用的标准方法

### 多级 FDR
- **PSM FDR 1%** → **肽段 FDR 1%** → **蛋白 FDR 1%**
- 每一级独立过滤，逐级收紧
- 发表级别分析通常使用更严格的蛋白 FDR（如 0.1%）

## 结果质量评估

| 指标 | 参考范围（HeLa 标准样品） |
|------|--------------------------|
| 蛋白组数 | 3,000 - 6,000+ |
| Unique peptides/蛋白 中位数 | ≥ 2 |
| 序列覆盖率中位数 | 15-30% |
| 1-peptide 蛋白占比 | < 30%（过高提示数据深度不足） |

## 适用场景
- 标准蛋白质组学实验的蛋白水平报告
- 比较不同样品的蛋白鉴定差异
- 验证目标蛋白是否被鉴定到
```

- [ ] **Step 2: Commit**

```bash
git add .github/prompts/protein-inference.prompt.md
git commit -m "docs: add protein inference prompt

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 6: Create database-management.prompt.md

**Files:**
- Create: `.github/prompts/database-management.prompt.md`

- [ ] **Step 1: Create file**

```markdown
---
mode: agent
description: "FASTA 数据库管理 — 查看、下载和管理蛋白质序列数据库"
---

# 数据库管理

管理蛋白质序列数据库（FASTA 格式），支持内置数据库自动下载和本地缓存。

## 内置数据库

| ID | 物种 | 来源 | 关键词（自动匹配） |
|----|------|------|---------------------|
| human_swissprot | 人 (Homo sapiens) | UniProt Swiss-Prot | human, 人, 人类, homo sapiens, 9606 |
| mouse_swissprot | 小鼠 (Mus musculus) | UniProt Swiss-Prot | mouse, 小鼠, mus musculus, 10090 |
| ecoli_swissprot | 大肠杆菌 (E. coli) | UniProt Swiss-Prot | ecoli, e.coli, 大肠杆菌, escherichia |
| yeast_swissprot | 酵母 (S. cerevisiae) | UniProt Swiss-Prot | yeast, 酵母, saccharomyces |
| arabidopsis_swissprot | 拟南芥 (A. thaliana) | UniProt Swiss-Prot | arabidopsis, 拟南芥 |
| crap | 污染物 | cRAP | contaminant, 污染, crap |

## 流程

### 查看数据库状态
1. 调用 `list_databases()` 查看所有数据库
2. 每个数据库显示：Available（可下载）或 Downloaded（已缓存，含文件大小和蛋白数量）

### 下载数据库
1. 调用 `download_database(database_id="human_swissprot")`
2. 从 UniProt 通过 HTTPS 下载，自动解析 FASTA 统计蛋白数量
3. 缓存到本地目录，下次使用无需重新下载
4. 返回本地文件路径，可直接用作搜索的 `database_path`
5. 如需更新：`download_database(database_id="human_swissprot", force=true)`

### 查看数据库详情
- 调用 `get_database_info(database_id="human_swissprot")` 查看：
  - 蛋白序列数量
  - 文件大小
  - SHA256 校验和
  - 下载时间
  - 前 5 个蛋白 accession（验证正确性）

### 自动数据库解析（推荐）
- 使用 `prepare_search(organism="human")` 时自动处理：
  1. 检查本地缓存
  2. 未缓存则自动下载
  3. 填充到搜索参数的 database_path
- 支持中英文物种名和 NCBI Taxonomy ID

## cRAP 污染物数据库

- **Common Repository of Adventitious Proteins**
- 包含实验室常见污染物：人角蛋白、胰蛋白酶自切产物、BSA 等
- 建议在搜索时将 cRAP 序列合并到物种数据库中
- 高占比的 cRAP 鉴定可能提示样品制备问题

## 适用场景
- 首次使用时下载所需数据库
- 检查已缓存数据库是否需要更新
- 使用自定义 FASTA 数据库时，直接提供 database_path 即可
```

- [ ] **Step 2: Commit**

```bash
git add .github/prompts/database-management.prompt.md
git commit -m "docs: add database management prompt

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 7: Create sage-search.prompt.md

**Files:**
- Create: `.github/prompts/sage-search.prompt.md`

- [ ] **Step 1: Create file**

```markdown
---
mode: agent
description: "Sage 搜索引擎 — 使用 sage-core 库进行高性能蛋白质数据库搜索"
---

# Sage 搜索引擎

使用 sage-core（v0.15.0）进行生产级蛋白质数据库搜索。Sage 是一个高性能开源搜索引擎，
内置于 ProteinCopilot 作为库调用（非子进程）。

## 使用方式

在搜索参数中指定引擎：

```
# 方式 1：通过 prepare_search
prepare_search(input_files=[...], organism="human", engine="Sage")

# 方式 2：通过 run_search
run_search(params={...engine: "Sage"...}, input_files=[...])
```

## Sage vs SimpleSearch 对比

| 特性 | Sage | SimpleSearch |
|------|------|--------------|
| 打分算法 | Hyperscore + LDA rescoring | 简化打分 |
| 并行化 | rayon 多线程 | 单线程 |
| FDR 计算 | 三级（spectrum/peptide/protein） | 基础 target-decoy |
| 适用场景 | 生产分析 | 快速测试、小规模数据 |
| 数据规模 | 大规模（万级谱图） | 小规模（千级谱图） |

## Sage 搜索流程

1. **参数转换**：ProteinCopilot SearchParams → Sage Parameters
   - 消化酶映射（Trypsin → KR|P 正则规则）
   - 修饰映射（名称 → 质量偏移）
   - 质量容差映射（Da/ppm → SageTolerance）
2. **数据库构建**：FASTA → IndexedDatabase（内存中构建 + 索引）
3. **打分**：rayon 并行谱图匹配，计算 hyperscore
4. **LDA Rescoring**：线性判别分析重打分，综合多维特征
5. **FDR 计算**：
   - Spectrum-level q-value（基于 discriminant_score）
   - Peptide-level q-value（最佳 PSM 代表）
   - Protein-level q-value（picked-protein 方法）

## 结果字段说明

搜索结果的 `extra` 字段包含 Sage 特有信息：

| 字段 | 说明 |
|------|------|
| `hyperscore` | Sage 原始打分（越高越好） |
| `discriminant_score` | LDA 重打分后的综合分数（用于排序） |
| `spectrum_q` | 谱图水平 q-value |
| `peptide_q` | 肽段水平 q-value |
| `protein_q` | 蛋白水平 q-value |
| `delta_hyperscore` | 最佳与次佳 PSM 的 hyperscore 差值 |
| `matched_intensity_pct` | 匹配碎片离子强度占总强度的百分比 |
| `poisson` | 随机匹配概率（Poisson 模型） |

## Sage 特有参数（内部默认值）

| 参数 | 默认值 | 说明 |
|------|--------|------|
| 肽段质量范围 | 500-5000 Da | 过滤过短/过长肽段 |
| min_ion_index | 2 | 跳过前 2 个碎片离子（通常噪声高） |
| max_variable_mods | 2 | 每条肽段最多 2 个可变修饰 |
| chimera | false | 不启用嵌合谱图处理 |

这些参数在 SageAdapter 中硬编码为合理默认值，无需用户调整。

## 引擎健康检查

```
check_engine(engine="Sage")
```

返回引擎名称、版本、健康状态，以及所有已注册引擎列表。

## 适用场景
- 生产级蛋白质组学搜索
- 大规模数据集（>10,000 谱图）
- 需要多级 FDR 控制的发表级分析
- 需要高性能并行处理时
```

- [ ] **Step 2: Commit**

```bash
git add .github/prompts/sage-search.prompt.md
git commit -m "docs: add Sage search engine prompt

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 8: Create batch-search.prompt.md

**Files:**
- Create: `.github/prompts/batch-search.prompt.md`

- [ ] **Step 1: Create file**

```markdown
---
mode: agent
description: "批量搜索 — 多文件质谱数据的批量搜索策略"
---

# 批量搜索

处理包含多个质谱文件的实验，如多样品比较、技术重复、分级分离等场景。

## 使用方式

`run_search` 的 `input_files` 参数接受文件列表：

```
run_search(
  params={...},
  input_files=["sample1.mzML", "sample2.mzML", "sample3.mzML"]
)
```

或通过 `prepare_search`：

```
prepare_search(
  input_files=["sample1.mzML", "sample2.mzML", "sample3.mzML"],
  organism="human"
)
```

## 流程

### 1. 数据概览
- 调用 `read_spectra` 检查**第一个文件**的数据特征
- 确认所有文件来自相同实验条件（相同仪器、相同采集模式）

### 2. 参数推荐
- `prepare_search` 或 `recommend_params` 基于第一个文件推荐参数
- 同一批次的文件通常使用**相同的搜索参数**

### 3. 执行搜索
- 将所有文件一次性传给 `run_search`
- 搜索引擎内部合并处理所有谱图
- 返回单个 `run_id`，统一管理

### 4. 蛋白推断（推荐）
- 多文件搜索后，调用 `infer_proteins(run_id=xxx, fasta_path=xxx)` 聚合蛋白结果
- 跨文件的 PSM 被合并后再做蛋白推断，提高覆盖率

### 5. 结果导出
- `generate_summary(run_id)` 显示合并统计
- `export_results(run_id)` 导出包含所有文件结果的 TSV/JSON

## 分组策略

| 场景 | 策略 |
|------|------|
| 同一样品的技术重复 | 合并到一个 run_search 调用 |
| 分级分离（fractionation） | 合并到一个 run_search 调用 |
| 不同实验条件 | 分别搜索，各自独立 run_id |
| 不同物种 | 必须分别搜索（不同数据库） |
| DDA + DIA 混合 | 必须分别处理（不同工作流） |

## 性能预期

- **SimpleSearch**：~100 谱图/秒（单线程）
- **Sage**：~1000-5000 谱图/秒（多线程，取决于 CPU 核数）
- 10 个文件 × 10,000 谱图/文件 = 100,000 谱图 → Sage 约 20-100 秒

## 注意事项
- 所有文件必须存在且可读（搜索前同步校验）
- 参数推荐基于第一个文件的谱图特征
- 建议同一批次文件使用相同的仪器参数
- 跨文件的蛋白推断需要 `infer_proteins`，不会自动执行

## 适用场景
- 多样品比较蛋白质组学实验
- 分级分离样品（SCX、高 pH RP 等）
- 技术或生物学重复实验
- 时间序列实验（不同时间点的样品）
```

- [ ] **Step 2: Commit**

```bash
git add .github/prompts/batch-search.prompt.md
git commit -m "docs: add batch search prompt

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 9: Final verification and format

**Files:**
- All modified files

- [ ] **Step 1: Run full build**

Run: `cargo build --workspace 2>&1 | tail -3`
Expected: `Finished` with no errors

- [ ] **Step 2: Run full test suite**

Run: `cargo test --workspace 2>&1 | grep -E "^test result|FAILED"`
Expected: All pass, 0 failed

- [ ] **Step 3: Run clippy**

Run: `cargo clippy --workspace 2>&1 | tail -3`
Expected: 0 warnings

- [ ] **Step 4: Run cargo fmt**

Run: `cargo fmt --all`

- [ ] **Step 5: Final commit if formatting changed**

```bash
git add -A
git diff --cached --stat  # only commit if changes exist
git commit -m "style: apply cargo fmt

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

- [ ] **Step 6: Verify git log**

Run: `git log --oneline -10`
Expected: All task commits present, clean history
