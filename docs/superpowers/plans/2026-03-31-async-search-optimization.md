# Async Search Optimization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add real-time stage-based progress reporting, search cancellation, and disk-persisted search history to the MCP server.

**Architecture:** Extend `SearchEngineAdapter` trait with `on_progress` callback and `cancel()` method. Store `JoinHandle` in `RunState` for abort. Persist completed search metadata as JSON files in `~/.protein-copilot/history/`. Add `cancel_search` and `list_searches` MCP tools. Update Agent instructions for polling and cancellation workflows.

**Tech Stack:** Rust, tokio (JoinHandle::abort), serde_json (history persistence), async-trait, rmcp

**Spec:** `docs/superpowers/specs/2026-03-31-async-search-optimization-design.md`

---

## File Map

| Action | File | Responsibility |
|--------|------|----------------|
| Modify | `crates/core/src/engine.rs` | Add `ProgressCallback` type, `on_progress` param to `search()`, default `cancel()` method |
| Modify | `crates/search-engine/src/progress.rs` | Add `stage` field to `SearchProgress` |
| Modify | `crates/search-engine/src/simple_engine.rs` | Call `on_progress` at each search phase |
| Modify | `crates/search-engine/src/adapters/pfind.rs` | Update stub signature, add `cancel()` stub |
| Modify | `crates/search-engine/src/lib.rs` | Re-export new types |
| Create | `crates/mcp-server/src/history.rs` | JSON file persistence for search history |
| Modify | `crates/mcp-server/src/tools.rs` | Store JoinHandle, add `cancel_search`/`list_searches` tools, construct progress callback |
| Modify | `crates/mcp-server/src/main.rs` | Register history module |
| Modify | `.github/agents/proteomics-search.agent.md` | Polling strategy, cancel/history instructions |
| Modify | `crates/search-engine/tests/e2e_integration.rs` | Update e2e tests for new search() signature |

---

### Task 1: Extend SearchProgress with `stage` field

**Files:**
- Modify: `crates/search-engine/src/progress.rs:8-20`

- [ ] **Step 1: Add `stage` field to SearchProgress**

```rust
// crates/search-engine/src/progress.rs — replace entire file
//! Search progress tracking.

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use uuid::Uuid;

/// Progress information for a running search.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SearchProgress {
    /// Unique identifier for the search run.
    pub run_id: Uuid,
    /// Current status: "Running", "Completed", "Failed: ...", or "Cancelled".
    pub status: String,
    /// Current search stage, e.g. "Matching spectra (300/1000)".
    pub stage: Option<String>,
    /// Progress percentage (0.0 to 1.0), `None` if indeterminate.
    pub progress_pct: Option<f64>,
    /// Elapsed time in seconds.
    pub elapsed_sec: f64,
    /// Estimated remaining time in seconds, `None` if unknown.
    pub estimated_remaining_sec: Option<f64>,
}
```

- [ ] **Step 2: Fix all compile errors from the new field**

Every place that constructs `SearchProgress` must now include `stage: None` (or `stage: Some(...)` as appropriate). The main place is `crates/mcp-server/src/tools.rs` around line 497. Search for `SearchProgress {` across the workspace and add the field.

Run: `cargo build --workspace 2>&1 | head -40`
Expected: Compilation errors pointing to missing `stage` field. Fix each one by adding `stage: None`.

- [ ] **Step 3: Run tests**

Run: `cargo test --workspace 2>&1 | tail -5`
Expected: All 310 tests pass.

- [ ] **Step 4: Commit**

```bash
git add -A && git commit -m "feat(search-engine): add stage field to SearchProgress

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: Extend SearchEngineAdapter trait with `on_progress` and `cancel()`

**Files:**
- Modify: `crates/core/src/engine.rs:14-95`

- [ ] **Step 1: Write test for the new ProgressCallback type**

Add to the existing test module in `crates/core/src/engine.rs`:

```rust
#[test]
fn progress_callback_is_send_sync() {
    // Ensure ProgressCallback can be sent across threads
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<Box<dyn Fn(protein_copilot_search_engine::SearchProgress) + Send + Sync>>();
}
```

Run: `cargo test -p protein-copilot-core -- progress_callback 2>&1`
Expected: FAIL — `ProgressCallback` type not defined yet.

- [ ] **Step 2: Add imports, ProgressCallback type alias, and update trait**

In `crates/core/src/engine.rs`, add the `uuid` import and `SearchProgress` re-definition concern. Since `SearchProgress` lives in `search-engine` (which depends on `core`), we cannot import it in `core`. Instead, define the callback as accepting a generic progress struct, or move `SearchProgress` to `core`.

**Best approach: move SearchProgress to core.** It's a shared data structure.

Create `crates/core/src/progress.rs`:

```rust
//! Search progress tracking (shared data structure).

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use uuid::Uuid;

/// Progress information for a running search.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SearchProgress {
    /// Unique identifier for the search run.
    pub run_id: Uuid,
    /// Current status: "Running", "Completed", "Failed: ...", or "Cancelled".
    pub status: String,
    /// Current search stage, e.g. "Matching spectra (300/1000)".
    pub stage: Option<String>,
    /// Progress percentage (0.0 to 1.0), `None` if indeterminate.
    pub progress_pct: Option<f64>,
    /// Elapsed time in seconds.
    pub elapsed_sec: f64,
    /// Estimated remaining time in seconds, `None` if unknown.
    pub estimated_remaining_sec: Option<f64>,
}

/// Callback type for progress reporting from search engines.
pub type ProgressCallback = Box<dyn Fn(SearchProgress) + Send + Sync>;

/// A no-op progress callback for cases where progress is not needed.
pub fn noop_progress() -> ProgressCallback {
    Box::new(|_| {})
}
```

Add to `crates/core/src/lib.rs`:

```rust
pub mod progress;
```

- [ ] **Step 3: Update search-engine to re-export from core instead of defining its own**

In `crates/search-engine/src/progress.rs`, replace the content with:

```rust
//! Re-export progress types from core.
pub use protein_copilot_core::progress::*;
```

In `crates/search-engine/src/lib.rs`, the existing `pub use progress::SearchProgress;` will still work via re-export.

- [ ] **Step 4: Update SearchEngineAdapter trait**

In `crates/core/src/engine.rs`, update the imports and trait:

```rust
use crate::progress::{ProgressCallback, SearchProgress};
use uuid::Uuid;
```

Update the trait:

```rust
#[async_trait::async_trait]
pub trait SearchEngineAdapter: Send + Sync {
    /// Execute a search with the given parameters against input spectrum files.
    /// The `on_progress` callback is invoked at each search stage.
    async fn search(
        &self,
        params: &SearchParams,
        input_files: &[PathBuf],
        on_progress: ProgressCallback,
    ) -> Result<SearchResult, CoreError>;

    fn engine_info(&self) -> EngineInfo;

    async fn health_check(&self) -> Result<HealthStatus, CoreError>;

    /// Cancel a running search. Default is no-op (relies on JoinHandle::abort).
    /// pFind adapter can override to SSH kill the remote process.
    async fn cancel(&self, _run_id: Uuid) -> Result<(), CoreError> {
        Ok(())
    }
}
```

- [ ] **Step 5: Fix all callers — SimpleSearchEngine, PFindAdapter, tests**

Update `crates/search-engine/src/simple_engine.rs` (line 155-162):

```rust
#[async_trait::async_trait]
impl SearchEngineAdapter for SimpleSearchEngine {
    async fn search(
        &self,
        params: &SearchParams,
        input_files: &[PathBuf],
        _on_progress: ProgressCallback,
    ) -> Result<SearchResult, CoreError> {
        self.run_search(params, input_files)
            .map_err(CoreError::from)
    }
    // ... engine_info and health_check unchanged
}
```

Update `crates/search-engine/src/adapters/pfind.rs` — add `_on_progress: ProgressCallback` param and `cancel()` stub:

```rust
async fn search(
    &self,
    _params: &SearchParams,
    _input_files: &[PathBuf],
    _on_progress: ProgressCallback,
) -> Result<SearchResult, CoreError> {
    Err(CoreError::SearchEngineError { ... }) // existing error
}

async fn cancel(&self, _run_id: Uuid) -> Result<(), CoreError> {
    Err(CoreError::SearchEngineError {
        engine: "pFind".to_string(),
        detail: "cancel not yet implemented".to_string(),
        suggestion: "pFind remote cancellation requires SSH integration".to_string(),
    })
}
```

Update all test files that call `engine.search(&params, &files)` to add `noop_progress()`:

- `crates/search-engine/src/simple_engine.rs` (tests)
- `crates/search-engine/tests/e2e_integration.rs`
- `crates/search-engine/tests/integration.rs`
- `crates/mcp-server/src/tools.rs` (run_search)

Use: `engine.search(&params, &files, noop_progress()).await`

- [ ] **Step 6: Build and test**

Run: `cargo build --workspace && cargo test --workspace 2>&1 | tail -5`
Expected: All tests pass.

- [ ] **Step 7: Commit**

```bash
git add -A && git commit -m "feat(core): add ProgressCallback and cancel() to SearchEngineAdapter trait

Move SearchProgress from search-engine to core crate.
Add on_progress param to search() and default cancel() method.
All existing callers updated with noop_progress().

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: Implement stage-based progress in SimpleSearchEngine

**Files:**
- Modify: `crates/search-engine/src/simple_engine.rs:40-162`

- [ ] **Step 1: Write test for progress callbacks**

Add to `crates/search-engine/src/simple_engine.rs` test module:

```rust
#[tokio::test]
async fn search_reports_progress_stages() {
    use std::sync::{Arc, Mutex};
    use protein_copilot_core::progress::SearchProgress;

    let stages: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
    let stages_clone = Arc::clone(&stages);

    let on_progress: ProgressCallback = Box::new(move |p: SearchProgress| {
        if let Some(stage) = p.stage {
            stages_clone.lock().unwrap().push(stage);
        }
    });

    let engine = SimpleSearchEngine::new();
    let params = test_params(); // existing test helper
    let result = engine
        .search(&params, &[test_mgf_path()], on_progress)
        .await;
    assert!(result.is_ok());

    let recorded = stages.lock().unwrap();
    assert!(recorded.iter().any(|s| s.contains("Reading FASTA")));
    assert!(recorded.iter().any(|s| s.contains("Digesting")));
    assert!(recorded.iter().any(|s| s.contains("Matching")));
    assert!(recorded.iter().any(|s| s.contains("Aggregating")));
}
```

Run: `cargo test -p protein-copilot-search-engine -- search_reports_progress 2>&1`
Expected: FAIL — stages list is empty because `run_search` doesn't call progress yet.

- [ ] **Step 2: Refactor `run_search` to accept and call `on_progress`**

Update `run_search` signature and add progress callbacks at each stage:

```rust
fn run_search(
    &self,
    params: &SearchParams,
    input_files: &[PathBuf],
    on_progress: &dyn Fn(SearchProgress),
) -> Result<SearchResult, SearchEngineError> {
    let start = Instant::now();
    let run_id = Uuid::new_v4();

    // Helper to build progress
    let report = |stage: &str, pct: f64| {
        on_progress(SearchProgress {
            run_id,
            status: "Running".to_string(),
            stage: Some(stage.to_string()),
            progress_pct: Some(pct),
            elapsed_sec: start.elapsed().as_secs_f64(),
            estimated_remaining_sec: None,
        });
    };

    // Step 1: Validate
    params.validate().map_err(|e| SearchEngineError::InvalidParams { detail: e.to_string() })?;
    if input_files.is_empty() {
        return Err(SearchEngineError::NoInputSpectra);
    }

    // Step 2: Read FASTA
    report("Reading FASTA database", 0.02);
    let fasta_path = Path::new(&params.database_path);
    let proteins = parse_fasta(fasta_path)?;

    // Step 3: Digest
    report("Digesting proteins", 0.08);
    let mut all_peptides: Vec<DigestedPeptide> = Vec::new();
    for protein in &proteins {
        let peptides = digest(&protein.sequence, &protein.accession, &params.enzyme, params.missed_cleavages);
        all_peptides.extend(peptides);
    }
    if all_peptides.is_empty() {
        return Err(SearchEngineError::ExecutionError { detail: format!("no candidate peptides generated from {} proteins", proteins.len()) });
    }

    // Step 4: Read spectra
    report("Reading spectra", 0.15);
    let mut all_spectra: Vec<Spectrum> = Vec::new();
    for file_path in input_files {
        let info = protein_copilot_spectrum_io::detect_format(file_path).map_err(|e| SearchEngineError::IoError { detail: e.to_string() })?;
        let reader = protein_copilot_spectrum_io::create_reader(&info);
        let spectra = reader.read_all(file_path).map_err(|e| SearchEngineError::IoError { detail: e.to_string() })?;
        all_spectra.extend(spectra);
    }
    if all_spectra.is_empty() {
        return Err(SearchEngineError::NoInputSpectra);
    }

    // Step 5: Match spectra
    let total = all_spectra.len();
    let mut psms: Vec<Psm> = Vec::new();
    for (i, spectrum) in all_spectra.iter().enumerate() {
        if i % 50 == 0 || i == total - 1 {
            let pct = 0.15 + 0.75 * (i as f64 / total as f64);
            report(&format!("Matching spectra ({}/{})", i + 1, total), pct);
        }
        if let Some(m) = match_spectrum(spectrum, &all_peptides, &params.precursor_tolerance, &params.fragment_tolerance, &params.fixed_modifications) {
            psms.push(build_psm(spectrum, &m, &params.fixed_modifications));
        }
    }

    // Step 6: Aggregate
    report("Aggregating results", 0.92);
    let peptides = aggregate_peptides(&psms);
    let protein_results = aggregate_proteins(&psms, &proteins);

    let duration = start.elapsed().as_secs_f64();
    let summary = build_summary(&psms, all_spectra.len() as u64, duration);

    let engine_info = self.engine_info();
    let mut metadata = RunMetadata::new(params.clone(), engine_info.clone(), input_files.to_vec());
    metadata.run_id = run_id;
    metadata.status = RunStatus::Completed;
    metadata.duration_sec = Some(duration);

    Ok(SearchResult { run_id, engine_info, params_used: params.clone(), psms, peptides, proteins: protein_results, summary, metadata })
}
```

Update the trait impl to forward the callback:

```rust
async fn search(&self, params: &SearchParams, input_files: &[PathBuf], on_progress: ProgressCallback) -> Result<SearchResult, CoreError> {
    self.run_search(params, input_files, &*on_progress).map_err(CoreError::from)
}
```

- [ ] **Step 3: Run tests**

Run: `cargo test --workspace 2>&1 | tail -10`
Expected: All tests pass including the new `search_reports_progress_stages`.

- [ ] **Step 4: Commit**

```bash
git add -A && git commit -m "feat(search-engine): report stage-based progress during search

SimpleSearchEngine now calls on_progress at 4 stages:
Reading FASTA → Digesting → Matching spectra (N/M) → Aggregating.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 4: Implement cancel_search with JoinHandle::abort

**Files:**
- Modify: `crates/mcp-server/src/tools.rs`

- [ ] **Step 1: Extend RunState to store JoinHandle**

Update the `RunState` struct (around line 187):

```rust
struct RunState {
    progress: SearchProgress,
    result: Option<SearchResult>,
    handle: Option<tokio::task::JoinHandle<()>>,
}
```

Fix all `RunState` constructions to include `handle: None`, then in `run_search` store `handle: Some(handle)` after `tokio::spawn`.

- [ ] **Step 2: Add CancelSearchInput and implement cancel_search tool**

Add the input struct near the other input structs:

```rust
#[derive(Debug, Deserialize, JsonSchema)]
struct CancelSearchInput {
    /// Run ID of the search to cancel.
    run_id: String,
}
```

Add the tool inside the `#[rmcp::tool_router]` impl:

```rust
#[rmcp::tool(
    name = "cancel_search",
    description = "Cancel a running search. The search task is immediately terminated and status is set to Cancelled."
)]
fn cancel_search(
    &self,
    Parameters(input): Parameters<CancelSearchInput>,
) -> Result<Json<SearchProgress>, ErrorData> {
    let id = Uuid::parse_str(&input.run_id)
        .map_err(|_| mcp_err(ErrorCode::INVALID_PARAMS, "invalid run_id format"))?;
    let mut cache = self.run_cache.lock()
        .map_err(|_| mcp_err(ErrorCode::INTERNAL_ERROR, "cache lock failed"))?;
    let state = cache.get_mut(&id)
        .ok_or_else(|| mcp_err(ErrorCode::INVALID_PARAMS, format!("run_id {id} not found")))?;

    if state.progress.status != "Running" {
        return Err(mcp_err(
            ErrorCode::INVALID_PARAMS,
            format!("search is not running (status: {})", state.progress.status),
        ));
    }

    // Abort the tokio task
    if let Some(handle) = state.handle.take() {
        handle.abort();
    }

    state.progress.status = "Cancelled".to_string();
    state.progress.stage = Some("Cancelled by user".to_string());
    state.progress.progress_pct = None;

    Ok(Json(state.progress.clone()))
}
```

- [ ] **Step 3: Update PanicGuard to respect Cancelled status**

In the `PanicGuard::Drop` impl (around line 257-268), change the condition:

```rust
impl Drop for PanicGuard {
    fn drop(&mut self) {
        if let Ok(mut cache) = self.cache.lock() {
            if let Some(state) = cache.get_mut(&self.run_id) {
                if state.progress.status == "Running" {
                    // Only overwrite if still Running (not Cancelled or already Failed)
                    state.progress.status = "Failed: task panicked".to_string();
                    state.progress.elapsed_sec = self.start.elapsed().as_secs_f64();
                    state.progress.progress_pct = None;
                }
            }
        }
    }
}
```

- [ ] **Step 4: Update run_search to construct progress callback and store handle**

In `run_search`, construct the `on_progress` callback that writes to the cache, and store the JoinHandle:

```rust
// Construct progress callback that updates cache
let progress_cache = Arc::clone(&self.run_cache);
let progress_run_id = run_id;
let on_progress: ProgressCallback = Box::new(move |p: SearchProgress| {
    if let Ok(mut cache) = progress_cache.lock() {
        if let Some(state) = cache.get_mut(&progress_run_id) {
            if state.progress.status == "Running" {
                state.progress.stage = p.stage;
                state.progress.progress_pct = p.progress_pct;
                state.progress.elapsed_sec = p.elapsed_sec;
                state.progress.estimated_remaining_sec = p.estimated_remaining_sec;
            }
        }
    }
});

// ... spawn and store handle
let handle = tokio::spawn(async move {
    // ... existing search logic with on_progress
    let search_result = engine.search(&params, &files, on_progress).await;
    // ... existing cache update logic
});

// Store handle in cache
if let Ok(mut cache) = self.run_cache.lock() {
    if let Some(state) = cache.get_mut(&run_id) {
        state.handle = Some(handle);
    }
}
```

- [ ] **Step 5: Build and test**

Run: `cargo build --workspace && cargo test --workspace 2>&1 | tail -5`
Expected: All tests pass.

- [ ] **Step 6: Commit**

```bash
git add -A && git commit -m "feat(mcp-server): add cancel_search tool with JoinHandle::abort

Store JoinHandle in RunState. cancel_search aborts the task and
sets status to Cancelled. PanicGuard respects Cancelled status.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 5: Implement search history persistence

**Files:**
- Create: `crates/mcp-server/src/history.rs`
- Modify: `crates/mcp-server/src/tools.rs`

- [ ] **Step 1: Create history module**

Create `crates/mcp-server/src/history.rs`:

```rust
//! Search history persistence — JSON files in ~/.protein-copilot/history/

use std::fs;
use std::path::{Path, PathBuf};

use chrono::{DateTime, Utc};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use uuid::Uuid;

use protein_copilot_core::engine::EngineInfo;
use protein_copilot_core::search_params::SearchParams;

const MAX_HISTORY: usize = 500;

/// Summary metadata persisted for each completed/failed/cancelled search.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SearchHistoryEntry {
    pub run_id: Uuid,
    pub status: String,
    pub created_at: DateTime<Utc>,
    pub elapsed_sec: f64,
    pub engine_info: EngineInfo,
    pub input_files: Vec<PathBuf>,
    pub params_used: SearchParams,
    pub total_psms: Option<u64>,
    pub psms_at_1pct_fdr: Option<u64>,
    pub identification_rate: Option<f64>,
    pub protein_groups: Option<u64>,
}

/// Returns the history directory path, creating it if needed.
pub fn history_dir() -> Option<PathBuf> {
    let home = dirs::home_dir()?;
    let dir = home.join(".protein-copilot").join("history");
    fs::create_dir_all(&dir).ok()?;
    Some(dir)
}

/// Save a history entry to disk.
pub fn save_entry(entry: &SearchHistoryEntry) {
    let Some(dir) = history_dir() else {
        tracing::warn!("cannot determine history directory; skipping persistence");
        return;
    };
    let path = dir.join(format!("{}.json", entry.run_id));
    match serde_json::to_string_pretty(entry) {
        Ok(json) => {
            if let Err(e) = fs::write(&path, json) {
                tracing::warn!("failed to write history {}: {e}", path.display());
            }
        }
        Err(e) => tracing::warn!("failed to serialize history: {e}"),
    }
    evict_oldest(&dir);
}

/// Load all history entries from disk.
pub fn load_all() -> Vec<SearchHistoryEntry> {
    let Some(dir) = history_dir() else { return Vec::new() };
    let mut entries = Vec::new();
    let Ok(read_dir) = fs::read_dir(&dir) else { return entries };

    for entry in read_dir.flatten() {
        let path = entry.path();
        if path.extension().and_then(|e| e.to_str()) != Some("json") {
            continue;
        }
        match fs::read_to_string(&path) {
            Ok(content) => match serde_json::from_str::<SearchHistoryEntry>(&content) {
                Ok(e) => entries.push(e),
                Err(err) => tracing::warn!("corrupt history file {}: {err}", path.display()),
            },
            Err(e) => tracing::warn!("cannot read {}: {e}", path.display()),
        }
    }
    entries.sort_by(|a, b| b.created_at.cmp(&a.created_at));
    entries
}

/// FIFO eviction: delete oldest entries when over MAX_HISTORY.
fn evict_oldest(dir: &Path) {
    let mut entries = load_all();
    while entries.len() > MAX_HISTORY {
        if let Some(oldest) = entries.pop() {
            let path = dir.join(format!("{}.json", oldest.run_id));
            let _ = fs::remove_file(path);
        }
    }
}
```

Add `dirs` dependency to `crates/mcp-server/Cargo.toml`:

```toml
dirs = "5"
```

Add module declaration in `crates/mcp-server/src/tools.rs` or create a separate `mod history;` in a `lib.rs`. Since mcp-server is a bin crate, add at the top of `main.rs`:

```rust
mod history;
mod tools;
```

And move `tools.rs` to be a module of main. Or simpler: add `pub mod history;` as a submodule. The exact approach depends on the current structure — use `mod history;` at the top of `tools.rs`.

- [ ] **Step 2: Add `list_searches` MCP tool**

Add input struct and tool in `tools.rs`:

```rust
#[derive(Debug, Deserialize, JsonSchema)]
struct ListSearchesInput {
    /// Filter by status (e.g. "Completed", "Failed"). Optional.
    #[serde(default)]
    status_filter: Option<String>,
    /// Maximum results to return. Default 20.
    #[serde(default)]
    limit: Option<u32>,
}

#[rmcp::tool(
    name = "list_searches",
    description = "List recent search runs with their status, duration, and key metrics. Includes both active searches and completed history."
)]
fn list_searches(
    &self,
    Parameters(input): Parameters<ListSearchesInput>,
) -> Json<Vec<history::SearchHistoryEntry>> {
    let limit = input.limit.unwrap_or(20) as usize;
    let mut entries = history::load_all();

    // Merge active runs from cache
    if let Ok(cache) = self.run_cache.lock() {
        for (id, state) in cache.iter() {
            if !entries.iter().any(|e| e.run_id == *id) {
                entries.push(history::SearchHistoryEntry {
                    run_id: *id,
                    status: state.progress.status.clone(),
                    created_at: Utc::now(), // approximate
                    elapsed_sec: state.progress.elapsed_sec,
                    engine_info: EngineInfo { name: "SimpleSearch".into(), version: "0.1.0".into(), supported_features: vec![] },
                    input_files: vec![],
                    params_used: SearchParams::default_placeholder(),
                    total_psms: None,
                    psms_at_1pct_fdr: None,
                    identification_rate: None,
                    protein_groups: None,
                });
            }
        }
    }

    if let Some(ref filter) = input.status_filter {
        entries.retain(|e| e.status.starts_with(filter.as_str()));
    }
    entries.sort_by(|a, b| b.created_at.cmp(&a.created_at));
    entries.truncate(limit);
    Json(entries)
}
```

Note: `OrderedRunCache` needs an `iter()` method. Add it:

```rust
fn iter(&self) -> impl Iterator<Item = (&Uuid, &RunState)> {
    self.map.iter()
}
```

- [ ] **Step 3: Persist history on search completion/failure/cancellation**

In the `tokio::spawn` block inside `run_search`, after updating the cache on completion/failure, call `history::save_entry()`. Similarly in `cancel_search`.

- [ ] **Step 4: Build and test**

Run: `cargo build --workspace && cargo test --workspace 2>&1 | tail -5`
Expected: All tests pass.

- [ ] **Step 5: Commit**

```bash
git add -A && git commit -m "feat(mcp-server): add search history persistence and list_searches tool

Search metadata persisted as JSON in ~/.protein-copilot/history/.
FIFO eviction at 500 entries. list_searches merges active + history.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 6: Update Agent and Skill instructions

**Files:**
- Modify: `.github/agents/proteomics-search.agent.md`

- [ ] **Step 1: Add cancel_search and list_searches to tools list**

Update the frontmatter tools list to include the 2 new tools:

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
```

- [ ] **Step 2: Update search workflow with polling instructions**

Find the search execution step and update to include:

```markdown
Step 4: Execute Search
  - Call `run_search(input_files, database_path)` → returns `{run_id, status: "Running"}`
  - Report to user: "搜索已提交 (run_id: xxx)"

Step 5: Monitor Progress
  - Poll `get_search_status(run_id)` every 5-10 seconds
  - Report stage changes to user when stage field changes:
    - "正在读取蛋白数据库..."
    - "正在消化蛋白序列..."
    - "正在匹配谱图 (300/1000)..."
    - "正在聚合结果..."
  - If user says "停止", "取消", or "cancel", call `cancel_search(run_id)`
  - If status is "Completed", proceed to Step 6
  - If status starts with "Failed", report error details and suggest next steps
  - If status is "Cancelled", confirm to user and ask if they want to start a new search
```

- [ ] **Step 3: Add history query instructions**

Add a new section:

```markdown
## 历史查询

当用户询问"之前搜索过什么"或"搜索历史"时：
  - Call `list_searches(limit=10)` to get recent searches
  - Display as a table: run_id (shortened), status, duration, PSMs found, identification rate
  - User can then ask to re-export or review specific results
```

- [ ] **Step 4: Commit**

```bash
git add -A && git commit -m "docs(agent): update proteomics-search agent with polling, cancel, and history instructions

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 7: Integration tests and final verification

**Files:**
- Modify: `crates/search-engine/tests/e2e_integration.rs`

- [ ] **Step 1: Add progress tracking e2e test**

```rust
#[tokio::test]
async fn scenario_progress_tracking() {
    use std::sync::{Arc, Mutex};
    use protein_copilot_core::progress::{SearchProgress, ProgressCallback};

    let stages: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
    let stages_clone = Arc::clone(&stages);
    let on_progress: ProgressCallback = Box::new(move |p: SearchProgress| {
        if let Some(ref stage) = p.stage {
            let mut s = stages_clone.lock().unwrap();
            if s.last().map(|l| l != stage).unwrap_or(true) {
                s.push(stage.clone());
            }
        }
    });

    let file_info = detect_format(&mgf_path()).unwrap();
    let summary = create_reader(&file_info).read_summary(&mgf_path()).unwrap();
    let mut params = ParamRecommender.recommend(&summary, None).unwrap().decision;
    params.database_path = fasta_path().to_string_lossy().to_string();

    let engine = SimpleSearchEngine::new();
    let result = engine.search(&params, &[mgf_path()], on_progress).await.unwrap();
    assert!(!result.psms.is_empty());

    let recorded = stages.lock().unwrap();
    assert!(recorded.len() >= 4, "Expected at least 4 stages, got: {:?}", *recorded);
    assert!(recorded[0].contains("FASTA"), "First stage should be FASTA reading");
    assert!(recorded.iter().any(|s| s.contains("Matching")));
    assert!(recorded.last().unwrap().contains("Aggregating"));
}
```

- [ ] **Step 2: Update existing e2e tests to use noop_progress**

Every existing `engine.search(&params, &[...]).await` call needs the third argument. Add `protein_copilot_core::progress::noop_progress` import and use `noop_progress()` as the third arg.

- [ ] **Step 3: Run full verification**

```bash
cargo build --workspace
cargo test --workspace
cargo clippy --workspace
cargo fmt --check
```

Expected: All pass, 0 warnings, 0 format issues.

- [ ] **Step 4: Final commit**

```bash
git add -A && git commit -m "test: add progress tracking e2e test, update existing tests for new signature

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```
