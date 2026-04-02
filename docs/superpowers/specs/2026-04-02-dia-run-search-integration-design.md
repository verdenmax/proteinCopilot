# DIA run_search Integration Design

## Problem

`extract_dia_precursors` MCP tool extracts candidate precursors from DIA data and caches
enhanced spectra in `OrderedDiaCache`. But `run_search` only reads spectra from file paths —
it has no way to consume the cached DIA spectra. The DIA pipeline is broken at this handoff.

## Solution

Add an optional `dia_run_id` parameter to `RunSearchInput`. When provided, `run_search`
reads spectra from `dia_cache` instead of from files. The cached spectra are **moved** (not
cloned) out of the cache to avoid double memory usage.

## Changes

### 1. `RunSearchInput` (mcp-server/src/tools.rs)

Add field:

```rust
/// Optional run_id from extract_dia_precursors. When provided, uses cached
/// DIA-extracted spectra instead of reading from input_files.
dia_run_id: Option<String>,
```

When `dia_run_id` is `Some`, `input_files` may be empty (not required).

### 2. `run_search` tool method (mcp-server/src/tools.rs)

Before the file-reading step, check for `dia_run_id`:

```rust
let spectra = if let Some(ref run_id_str) = input.dia_run_id {
    let run_id = Uuid::parse_str(run_id_str)
        .map_err(|_| mcp_err(ErrorCode::INVALID_PARAMS, "invalid dia_run_id format"))?;
    let mut cache = self.dia_cache.lock()
        .map_err(|_| mcp_err(ErrorCode::INTERNAL_ERROR, "DIA cache lock is poisoned"))?;
    cache.remove(&run_id)
        .ok_or_else(|| mcp_err(ErrorCode::INVALID_PARAMS,
            &format!("dia_run_id '{}' not found in cache (may have been evicted)", run_id_str)))?
} else {
    // existing file-reading logic
};
```

Key decisions:
- **`remove`** instead of clone — moves spectra out of cache, frees memory
- Parse UUID with clear error on invalid format
- Clear error message if run_id not found (may have been evicted by FIFO)

### 3. `OrderedDiaCache` — add `remove` method

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

### 4. `extract_dia_precursors` output message

Update to mention the `dia_run_id` parameter name:

```
"DIA extraction complete. N precursors extracted from M MS2 spectra.
 Pass dia_run_id='<id>' to run_search to search these spectra."
```

## Not Changed

- `SearchEngineAdapter` trait — unchanged
- `simple_engine.rs` — unchanged (already handles multi-precursor via `match_spectrum_all`)
- `SearchParams` — unchanged
- `dia-extraction` crate — unchanged

## User Workflow (LLM Agent)

```
1. read_spectra(file_path) → see DIA characteristics
2. extract_dia_precursors(file_path) → get dia_run_id
3. run_search(dia_run_id=..., database_path=...) → search results
4. generate_summary(run_id) → interpret results
```
