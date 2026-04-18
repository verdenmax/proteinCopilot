# mzML Index Optimization Design

## Problem

When `annotate_spectrum` is called on a large (8GB) mzML file for the first time, `IndexedMzMLReader::open()` may take 30-60+ seconds to build a scan index, exceeding the MCP Client's request timeout (typically 60s). The server completes the work but the client reports a timeout error.

### Root Cause Chain

1. `annotate_spectrum` calls `get_or_create_reader(path)` (tools.rs:1021)
2. First call: `IndexedMzMLReader::open(path)` → `build_index_from_native_mzml(path)?`
3. If no native `<indexList>`: falls back to `build_index_by_scanning(path)`
4. Fallback reads entire 8GB file line-by-line (`BufReader::read_line()`)
5. ~30-60s elapsed → MCP Client timeout fires
6. Server continues and completes; client reports error `-32001: Request timed out`

### Why Subsequent Calls Are Fast

`get_or_create_reader` caches the `IndexedMzMLReader` in an LRU cache (capacity 8). After the first open, `read_spectrum(scan)` is O(1): seek to byte offset → parse single `<spectrum>` node. Even the SILAC heavy scan search (±50 scans, ~100 reads) takes <1s with cached reader.

## Solution: Three-Layer Index Resolution + Faster Fallback

### 1. Disk-Persisted Index Cache (`.mzml.idx`)

**New index load priority in `IndexedMzMLReader::open()`:**

```
1. Disk cache (.mzml.idx sidecar file) → milliseconds
2. Native mzML <indexList>             → milliseconds (reads EOF 4KB)
3. Byte-scan fallback (accelerated)    → 3-8s for 8GB file
```

After steps 2 or 3 succeed, the result is written to disk as `.mzml.idx` so future opens (including across MCP Server restarts) are instant.

**Sidecar file location:** Same directory as the mzML file, with `.idx` appended:
- `/data/experiment.mzML` → `/data/experiment.mzML.idx`

**Binary format:**

```
Offset  Size   Field                Description
0       4      magic                b"PCIX" (ProteinCopilot IndeX)
4       1      version              Format version (1)
5       8      source_file_size     u64 LE, mzML file size in bytes
13      8      source_file_mtime    u64 LE, mzML file modification time (Unix epoch secs)
21      4      entry_count          u32 LE, number of scan entries
25      N×12   entries              Array of (scan_number: u32 LE, byte_offset: u64 LE)
```

**Size estimate:** 100,000 spectra × 12 bytes = 1.2 MB + 25 bytes header ≈ 1.2 MB.

**Staleness check:** When loading `.idx`, compare `source_file_size` and `source_file_mtime` against the actual mzML file. If either differs, discard and rebuild.

**Write-failure tolerance:** If writing `.idx` fails (e.g., read-only filesystem), log a warning and continue — disk caching is an optimization, not a requirement.

### 2. Accelerated Fallback Scan

Current `build_index_by_scanning` uses `BufReader::read_line()` which has overhead from:
- Per-line String allocation
- UTF-8 validation of entire lines (including base64-encoded peak data)
- Cannot handle `<spectrum` tags split across lines (though rare)

**Replacement: block-based byte scanning**

```rust
fn build_index_by_byte_scan(path: &Path) -> Result<ScanIndex, SpectrumIoError> {
    let file = File::open(path)?;
    let mut reader = BufReader::with_capacity(64 * 1024, file);
    let needle = b"<spectrum ";
    let mut offsets = HashMap::new();
    let mut global_pos: u64 = 0;
    let mut fallback_scan: u32 = 0;

    loop {
        let buf = reader.fill_buf()?;
        if buf.is_empty() { break; }
        
        // Search for needle in current buffer
        let mut search_start = 0;
        while let Some(pos) = memchr::memmem::find(&buf[search_start..], needle) {
            let abs_pos = global_pos + (search_start + pos) as u64;
            fallback_scan += 1;
            
            // Extract id attribute from the tag (scan forward to '>')
            let tag_slice = &buf[search_start + pos..];
            if let Some(scan) = extract_scan_from_tag_bytes(tag_slice, fallback_scan) {
                offsets.insert(scan, abs_pos);
            }
            search_start += pos + needle.len();
        }
        
        // Consume buffer, keeping overlap for cross-boundary matches
        let consumed = buf.len().saturating_sub(needle.len());
        global_pos += consumed as u64;
        reader.consume(consumed);
    }
    
    Ok(ScanIndex::new(offsets, IndexSource::BuiltFromScan))
}
```

Key improvements:
- **No UTF-8 conversion**: operates on raw bytes
- **No per-line allocation**: searches in-place within the BufReader buffer
- **`memchr::memmem`**: SIMD-accelerated byte pattern search
- **Cross-boundary safety**: keeps `needle.len()` overlap between buffer fills

**Expected performance:** 8GB file from ~30-60s → ~3-8s (limited by disk I/O, not CPU).

### 3. MCP Client Timeout Configuration

The `.mcp.json` file supports a `timeout` field (in seconds) for each MCP server. Document the recommended configuration:

```json
{
  "mcpServers": {
    "protein-copilot": {
      "command": "cargo",
      "args": ["run", "--release", "-p", "protein-copilot-mcp-server"],
      "timeout": 300,
      "env": {
        "RUST_LOG": "info"
      }
    }
  }
}
```

This should be documented in:
- `README.md` installation/configuration section
- `.mcp.json` in the repository (add timeout field)

## Implementation Scope

### Changes to `crates/spectrum-io`

1. **New module: `src/disk_cache.rs`**
   - `save_index(path: &Path, index: &ScanIndex, file_size: u64, mtime: u64) -> Result<()>`
   - `load_index(path: &Path, expected_size: u64, expected_mtime: u64) -> Result<Option<ScanIndex>>`
   - Binary format read/write with the PCIX format described above

2. **Modified: `src/index.rs`**
   - Add `build_index_by_byte_scan()` using `memchr::memmem` for SIMD-accelerated scanning
   - Keep old `build_index_by_scanning()` as reference (can remove after validation)

3. **Modified: `src/indexed_mzml.rs`**
   - Update `IndexedMzMLReader::open()` to use three-layer resolution:
     1. Try disk cache → 2. Try native index → 3. Byte-scan fallback
   - After steps 2 or 3, persist to disk cache

4. **Modified: `Cargo.toml`**
   - Add `memchr` dependency for SIMD byte scanning

### Changes to `crates/mcp-server`

5. **No code changes needed** — the server already uses `get_or_create_reader()` which delegates to `IndexedMzMLReader::open()`. All optimizations are transparent.

### Configuration and Documentation

6. **Modified: `.mcp.json`** — add `timeout: 300`
7. **Modified: `README.md`** — document timeout configuration for large files

## Testing Strategy

- Unit tests for disk cache: write, read, staleness detection, corrupt file handling
- Unit tests for byte-scan: compare results with existing line-scan on small.mzml fixture
- Unit test for cross-buffer-boundary `<spectrum` tag detection
- Integration test: IndexedMzMLReader::open() with disk cache enabled
- Performance: not testable in CI, but add a benchmark target for local testing

## Non-Goals

- Async/streaming index building (too complex, changes API semantics)
- Memory-mapped I/O (portability concerns, less predictable memory usage)
- Index server or shared index across processes
