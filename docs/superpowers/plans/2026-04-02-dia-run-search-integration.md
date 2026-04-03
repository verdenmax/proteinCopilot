# DIA run_search Integration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wire `dia_cache` into `run_search` so DIA extraction → search is end-to-end.

**Architecture:** Add `OrderedDiaCache::remove()`, add `search_with_spectra()` to the search engine trait + SimpleSearchEngine, and add `dia_run_id` parameter to `RunSearchInput` to branch between file-based and cache-based search paths.

**Tech Stack:** Rust, rmcp, uuid, tokio

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `crates/mcp-server/src/tools.rs` | Modify | `remove()` on cache, `dia_run_id` field, wire into run_search |
| `crates/core/src/engine.rs` | Modify | Add `search_with_spectra()` to trait |
| `crates/search-engine/src/simple_engine.rs` | Modify | Implement `search_with_spectra()` |

---

### Task 1: Add `remove()` to OrderedDiaCache + `dia_run_id` to RunSearchInput

**Files:**
- Modify: `crates/mcp-server/src/tools.rs`

- [ ] **Step 1: Add `remove()` method to OrderedDiaCache**

In `crates/mcp-server/src/tools.rs`, add a `remove` method to the `impl OrderedDiaCache` block (near the existing `insert` method):

```rust
fn remove(&mut self, id: &Uuid) -> Option<Vec<Spectrum>> {
    if let Some(spectra) = self.entries.remove(id) {
        self.order.retain(|x| x != id);
        Some(spectra)
    } else {
        None
    }
}
```

- [ ] **Step 2: Add `dia_run_id` field to RunSearchInput**

In the `RunSearchInput` struct, add after `hints`:

```rust
    /// Optional run_id from extract_dia_precursors. When provided, uses cached
    /// DIA-extracted spectra instead of reading from input_files.
    dia_run_id: Option<String>,
```

- [ ] **Step 3: Update extract_dia_precursors output message**

Find the message in the `extract_dia_precursors` tool method. Change from:

```rust
message: format!(
    "DIA extraction complete. {} precursors extracted from {} MS2 spectra. \
     Results cached as run_id '{}'.",
    result.stats.total_precursors_extracted, result.stats.ms2_count, run_id
),
```

To:

```rust
message: format!(
    "DIA extraction complete. {} precursors extracted from {} MS2 spectra. \
     Pass dia_run_id=\"{}\" to run_search to search these spectra.",
    result.stats.total_precursors_extracted, result.stats.ms2_count, run_id
),
```

- [ ] **Step 4: Build and test**

```bash
cd /home/verden/pfind/2026-spring/code/proteinCopilot
cargo build && cargo test --quiet && cargo clippy -- -D warnings
```

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(mcp-server): add dia_run_id to RunSearchInput and cache remove

- OrderedDiaCache::remove() for moving spectra out of cache
- dia_run_id field on RunSearchInput for DIA pipeline integration
- Updated extract_dia_precursors message to guide user

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: Add `search_with_spectra()` to search engine

**Files:**
- Modify: `crates/core/src/engine.rs`
- Modify: `crates/search-engine/src/simple_engine.rs`

- [ ] **Step 1: Add method to SearchEngineAdapter trait**

In `crates/core/src/engine.rs`, add to the `SearchEngineAdapter` trait after the existing `search` method:

```rust
    /// Execute a search with pre-loaded spectra (e.g., from DIA extraction cache).
    /// Default implementation returns an error; engines must opt-in.
    async fn search_with_spectra(
        &self,
        params: &SearchParams,
        spectra: Vec<Spectrum>,
        on_progress: ProgressCallback,
    ) -> Result<SearchResult, CoreError> {
        let _ = (params, spectra, on_progress);
        Err(CoreError::SearchError(
            "this engine does not support search_with_spectra".to_string(),
        ))
    }
```

Add the import if not already present:

```rust
use crate::spectrum::Spectrum;
```

- [ ] **Step 2: Refactor SimpleSearchEngine to extract core matching logic**

In `crates/search-engine/src/simple_engine.rs`, the current `run_search` method reads files then does matching. Refactor by extracting the matching-and-scoring logic (everything AFTER file reading) into a private method:

```rust
/// Core search logic operating on pre-loaded spectra.
fn run_search_on_spectra(
    &self,
    params: &SearchParams,
    all_spectra: Vec<Spectrum>,
    on_progress: &dyn Fn(SearchProgress),
) -> Result<SearchResult, SearchEngineError> {
    // Move everything from run_search AFTER the file-reading loop here.
    // This starts at the MS2 filter: "let ms2_spectra: Vec<&Spectrum> = ..."
    // Through to the final "Ok(SearchResult { ... })" return.
    // The parameter validation + file reading stays in run_search.
}
```

Then simplify `run_search` to: validate params → read files → call `run_search_on_spectra`.

**IMPORTANT:** Read the actual current `run_search` code carefully. The split point is after the file-reading loop ends (after `all_spectra.extend(spectra);`). Everything before that stays in `run_search`. Everything from the MS2 filter onward moves to `run_search_on_spectra`.

- [ ] **Step 3: Implement `search_with_spectra` on SimpleSearchEngine**

Add the trait method implementation. Since `SimpleSearchEngine` implements `SearchEngineAdapter`, add:

```rust
async fn search_with_spectra(
    &self,
    params: &SearchParams,
    spectra: Vec<Spectrum>,
    on_progress: ProgressCallback,
) -> Result<SearchResult, CoreError> {
    let progress_fn = move |p: SearchProgress| {
        on_progress(p);
    };

    self.run_search_on_spectra(params, spectra, &progress_fn)
        .map_err(|e| CoreError::SearchError(e.to_string()))
}
```

Note: `run_search_on_spectra` needs to include parameter validation since it may be called directly without file-reading. Add `params.validate()` at the start of `run_search_on_spectra` if not already there.

- [ ] **Step 4: Run tests**

```bash
cd /home/verden/pfind/2026-spring/code/proteinCopilot
cargo test --quiet
cargo clippy -- -D warnings
```

All existing tests must pass — the refactoring is purely structural.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(search-engine): add search_with_spectra for pre-loaded spectra

- Extract run_search_on_spectra() from run_search() for code reuse
- Add search_with_spectra() to SearchEngineAdapter trait (default: error)
- SimpleSearchEngine implements it by delegating to run_search_on_spectra

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: Wire dia_run_id path in MCP tool's run_search

**Files:**
- Modify: `crates/mcp-server/src/tools.rs`

- [ ] **Step 1: Add DIA branch in run_search tool method**

In the `run_search` tool method, find the section where `input_files` are validated and the search is spawned. Add a branch: if `dia_run_id` is provided, extract spectra from cache and use `search_with_spectra`; otherwise use existing file-based path.

Before the input_files empty check (around line 636), add:

```rust
// DIA mode: use cached spectra from extract_dia_precursors
if let Some(ref run_id_str) = input.dia_run_id {
    let run_id = Uuid::parse_str(run_id_str).map_err(|_| {
        mcp_err(ErrorCode::INVALID_PARAMS, "invalid dia_run_id format")
    })?;

    let dia_spectra = {
        let mut cache = self.dia_cache.lock().map_err(|_| {
            mcp_err(ErrorCode::INTERNAL_ERROR, "DIA cache lock is poisoned")
        })?;
        cache.remove(&run_id).ok_or_else(|| {
            mcp_err(
                ErrorCode::INVALID_PARAMS,
                &format!(
                    "dia_run_id '{}' not found in cache (may have been evicted or already used)",
                    run_id_str
                ),
            )
        })?
    };

    // Use the DIA spectra for parameter recommendation if no params provided
    // (simplified: require params or database_path for DIA mode)

    let final_params = /* resolve params same as existing code */;

    let search_run_id = Uuid::new_v4();
    let engine = self.registry.get_engine(&final_params.engine_name())
        .ok_or_else(|| mcp_err(ErrorCode::INTERNAL_ERROR, "search engine not found"))?;
    let run_cache = self.run_cache.clone();

    tokio::spawn(async move {
        // ... similar to existing spawn block but calls
        // engine.search_with_spectra(&final_params, dia_spectra, callback).await
    });

    // Return run_id for status polling (same as existing)
}
```

**IMPORTANT:** Read the actual existing `run_search` method carefully to understand:
- How `final_params` is resolved (auto-recommendation flow)
- How the tokio task is spawned
- How the run_cache is populated with results
- How the PanicGuard works

The DIA path should mirror this structure but:
1. Skip `input_files` validation
2. Skip `read_summary` (we already have spectra)
3. Call `engine.search_with_spectra()` instead of `engine.search()`
4. Still use the same run_cache, PanicGuard, and result handling

For auto-recommendation when `params` is None in DIA mode: require that the user provides params explicitly (or at minimum `database_path`). If params is None and no database_path, return an error telling the user to provide params.

- [ ] **Step 2: Handle input_files validation for DIA mode**

The existing `input_files.is_empty()` check should be skipped when `dia_run_id` is provided. Adjust the validation:

```rust
if input.dia_run_id.is_none() && input.input_files.is_empty() {
    return Err(mcp_err(
        ErrorCode::INVALID_PARAMS,
        "input_files is empty — provide at least one spectrum file path, or use dia_run_id",
    ));
}
```

- [ ] **Step 3: Build and test**

```bash
cd /home/verden/pfind/2026-spring/code/proteinCopilot
cargo build
cargo test --quiet
cargo clippy -- -D warnings
```

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "feat(mcp-server): wire dia_run_id into run_search for end-to-end DIA

- When dia_run_id is provided, reads spectra from dia_cache instead of files
- Calls search_with_spectra() for pre-loaded DIA spectra
- Spectra moved out of cache on use (freed after search)
- Clear error messages for invalid/missing/evicted run_id

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```
