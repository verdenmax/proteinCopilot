# Full-Stack Observability Tracing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add comprehensive structured tracing to all 15 crates — 28 MCP tool spans, ~35 library function spans, 23 hot loop progress points, and subscriber enhancement with JSON output support.

**Architecture:** Three layers of instrumentation — L1 (MCP tool entry/exit spans), L2 (library function `#[instrument]`/`info_span!`), L3 (hot loop progress with rate/ETA). All output to stderr via `tracing-subscriber` with `FmtSpan::CLOSE` for auto-duration. `PROTEIN_LOG_JSON=1` toggles JSON format.

**Tech Stack:** `tracing 0.1`, `tracing-subscriber 0.3` (features: `env-filter`, `json`), `#[instrument]`, `info_span!`, `tracing::info!`

**Spec:** `docs/superpowers/specs/2026-04-28-observability-tracing-design.md`

---

## Task 1: Infrastructure — Dependencies & Subscriber Enhancement

**Files:**
- Modify: `Cargo.toml` (workspace root, add `json` feature to tracing-subscriber)
- Modify: `crates/fdr/Cargo.toml` (add tracing dep)
- Modify: `crates/report/Cargo.toml` (add tracing dep)
- Modify: `crates/param-recommend/Cargo.toml` (add tracing dep)
- Modify: `crates/core/Cargo.toml` (add tracing dep)
- Modify: `crates/mcp-server/Cargo.toml` (add json feature)
- Modify: `crates/mcp-server/src/main.rs` (enhance subscriber)

- [ ] **Step 1: Add tracing dependency to 4 crates missing it**

In `crates/fdr/Cargo.toml`, add under `[dependencies]`:
```toml
tracing = { workspace = true }
```

In `crates/report/Cargo.toml`, add under `[dependencies]`:
```toml
tracing = { workspace = true }
```

In `crates/param-recommend/Cargo.toml`, add under `[dependencies]`:
```toml
tracing = { workspace = true }
```

In `crates/core/Cargo.toml`, add under `[dependencies]`:
```toml
tracing = { workspace = true }
```

- [ ] **Step 2: Add `json` feature to workspace tracing-subscriber**

In root `Cargo.toml`, change:
```toml
tracing-subscriber = { version = "0.3", features = ["env-filter"] }
```
to:
```toml
tracing-subscriber = { version = "0.3", features = ["env-filter", "json"] }
```

- [ ] **Step 3: Enhance subscriber in main.rs**

Replace the subscriber block in `crates/mcp-server/src/main.rs` (lines 17-21):

```rust
// Initialize tracing (respects RUST_LOG env var)
// PROTEIN_LOG_JSON=1 switches to structured JSON output
let env_filter = tracing_subscriber::EnvFilter::try_from_default_env()
    .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info"));

let use_json = std::env::var("PROTEIN_LOG_JSON").map_or(false, |v| v == "1");

if use_json {
    use tracing_subscriber::layer::SubscriberExt;
    use tracing_subscriber::util::SubscriberInitExt;
    tracing_subscriber::registry()
        .with(env_filter)
        .with(
            tracing_subscriber::fmt::layer()
                .json()
                .with_writer(std::io::stderr)
                .with_span_events(tracing_subscriber::fmt::format::FmtSpan::CLOSE)
                .with_target(true)
                .with_timer(tracing_subscriber::fmt::time::uptime()),
        )
        .init();
} else {
    use tracing_subscriber::layer::SubscriberExt;
    use tracing_subscriber::util::SubscriberInitExt;
    tracing_subscriber::registry()
        .with(env_filter)
        .with(
            tracing_subscriber::fmt::layer()
                .with_writer(std::io::stderr)
                .with_span_events(tracing_subscriber::fmt::format::FmtSpan::CLOSE)
                .with_target(true)
                .with_timer(tracing_subscriber::fmt::time::uptime()),
        )
        .init();
}
```

- [ ] **Step 4: Verify build**

Run: `cargo build --workspace 2>&1 | tail -5`
Expected: build succeeds, no errors

- [ ] **Step 5: Verify subscriber works**

Run: `cargo test -p protein-copilot-mcp-server --lib -- --nocapture 2>&1 | head -20`
Expected: tracing output appears on stderr with uptime timestamps

- [ ] **Step 6: Commit**

```bash
git add -A && git commit -m "feat(tracing): infrastructure — deps + enhanced subscriber

Add tracing dependency to fdr, report, param-recommend, core crates.
Add json feature to tracing-subscriber for PROTEIN_LOG_JSON=1 support.
Enhance subscriber with FmtSpan::CLOSE (auto-duration), uptime timer,
and JSON/text format toggle via PROTEIN_LOG_JSON env var.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 2: L1 — MCP Tool Entry/Exit Spans (28 tools)

**Files:**
- Modify: `crates/mcp-server/src/tools.rs`

For each of the 28 `#[rmcp::tool]` methods, add an `info_span!` at the method entry and `tracing::info!` for start/completion. The pattern is the same for all tools — adapt the fields to each tool's context.

- [ ] **Step 1: Add spans to quick read tools (read_spectra, get_spectrum, recommend_params, list_presets, check_engine, list_databases, get_database_info, list_searches)**

For each tool, wrap the body in a span. Example for `read_spectra`:

```rust
fn read_spectra(
    &self,
    Parameters(input): Parameters<ReadSpectraInput>,
) -> Result<Json<SpectrumSummary>, ErrorData> {
    let span = tracing::info_span!("mcp_tool", name = "read_spectra");
    let _enter = span.enter();
    tracing::info!(file = %input.file_path, "started");

    validate_file_path(&input.file_path)?;
    let path = Path::new(&input.file_path);
    let reader = self.get_or_create_reader(path)?;
    let summary = reader
        .read_summary(path)
        .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

    tracing::info!(
        ms1 = summary.ms1_count,
        ms2 = summary.ms2_count,
        total = summary.total_spectra,
        "completed"
    );
    Ok(Json(summary))
}
```

Apply the same pattern to all 8 quick tools. For each tool, log the most relevant input (file path, database_id, run_id) on start and the most relevant output metric on completion.

- [ ] **Step 2: Add spans to search tools (run_search, get_search_status, cancel_search, prepare_search)**

For `run_search`, the span should be created inside the spawned async task so it covers the entire search duration. Log `engine`, `files`, `database` on start, and `psms_1pct`, `proteins`, `total_sec` on completion.

For `prepare_search`, log `files`, `organism`/`database_path` on start, and recommended params on completion.

- [ ] **Step 3: Add spans to spectrum annotation tools (annotate_spectrum, extract_xic, extract_spectrum_precursors, extract_dia_precursors, get_dia_cache_status)**

For `annotate_spectrum`, log `scan`, `peptide`, `charge` on start, and `matched_ions` on completion.

For `extract_xic`, log `peptide`, `scan`, `precursor_mz` on start, and output path on completion.

For `extract_dia_precursors`, log `file`, `output_mode` on start, and `candidates`, `ms2_count` on completion.

- [ ] **Step 4: Add spans to result tools (generate_summary, export_results, import_search_results, infer_proteins)**

For `generate_summary`, log `run_id` on start, `psms_1pct`, `id_rate` on completion.

For `export_results`, log `run_id`, `output_dir` on start, file counts on completion.

For `import_search_results`, log `result_file`, `format` on start, `psm_count`, `matched` on completion.

For `infer_proteins`, log `run_id` on start, `protein_groups`, `coverage` on completion.

- [ ] **Step 5: Add spans to entrapment tools (classify_entrapment_hits, analyze_entrapment_stats, annotate_provenance, find_similar_targets)**

For `classify_entrapment_hits`, log `results_file`, `target_fasta` on start, `trap_count`, `L0-L4` distribution on completion.

- [ ] **Step 6: Add spans to database tools (download_database)**

Log `database_id` on start, `path`, `protein_count` on completion.

- [ ] **Step 7: Add spans to diagnostic tools (diagnose_search)**

Log `run_id` on start, anomaly count on completion.

- [ ] **Step 8: Verify all 28 tools have spans**

Run: `grep -c 'info_span!.*mcp_tool' crates/mcp-server/src/tools.rs`
Expected: 28 (or close — some tools may share entry points)

Run: `cargo build -p protein-copilot-mcp-server 2>&1 | tail -3`
Expected: build succeeds

- [ ] **Step 9: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "feat(tracing): L1 — add entry/exit spans to all 28 MCP tools

Each tool now logs start (with key inputs) and completion (with key
metrics) at info level. Spans use mcp_tool{name=...} for filtering.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 3: L2 — spectrum-io Library Spans (8 spans)

**Files:**
- Modify: `crates/spectrum-io/src/indexed_mzml.rs`
- Modify: `crates/spectrum-io/src/index.rs`
- Modify: `crates/spectrum-io/src/disk_cache.rs`
- Modify: `crates/spectrum-io/src/mzml.rs`
- Modify: `crates/spectrum-io/src/lib.rs`

- [ ] **Step 1: Add span to IndexedMzMLReader::open()**

In `indexed_mzml.rs`, at the start of the `open()` method, add:

```rust
pub fn open(path: &Path) -> Result<Self, SpectrumIoError> {
    let file_size = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
    let span = tracing::info_span!("open_reader",
        file = %path.display(),
        file_size_mb = file_size / (1024 * 1024),
        index_source = tracing::field::Empty,
        scan_count = tracing::field::Empty,
    );
    let _enter = span.enter();
    // ... existing code ...
    // After index is built, record:
    span.record("index_source", match index.source() {
        IndexSource::NativeMzMLIndex => "native_index",
        IndexSource::DiskCache => "disk_cache",
        IndexSource::BuiltFromScan => "byte_scan",
    });
    span.record("scan_count", index.len() as u64);
    tracing::info!("indexed reader opened");
    // ...
}
```

- [ ] **Step 2: Add span to build_index_by_byte_scan()**

In `index.rs`, at the start of `build_index_by_byte_scan()`:

```rust
pub fn build_index_by_byte_scan(path: &Path) -> Result<ScanIndex, SpectrumIoError> {
    let file_size = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
    let span = tracing::info_span!("byte_scan_index",
        file = %path.display(),
        file_size_mb = file_size / (1024 * 1024),
        scan_count = tracing::field::Empty,
    );
    let _enter = span.enter();
    // ... existing code ...
    // Before return:
    span.record("scan_count", entries.len() as u64);
    tracing::info!("byte scan complete");
    Ok(ScanIndex::from_meta(entries, IndexSource::BuiltFromScan))
}
```

- [ ] **Step 3: Add spans to disk_cache load/save**

In `disk_cache.rs`, add `tracing::info!` at key points:

```rust
// In load function:
tracing::info!(path = %idx_path.display(), "loading disk cache");
// After successful load:
tracing::info!(scans = index.len(), "disk cache loaded");

// In save function:
tracing::info!(path = %idx_path.display(), scans = index.len(), "saving disk cache");
```

- [ ] **Step 4: Add span to for_each_spectrum and read_summary**

In `mzml.rs`, add `info_span!` to `read_summary()` and `for_each_spectrum()`:

```rust
fn read_summary(&self, path: &Path) -> Result<SpectrumSummary, SpectrumIoError> {
    let span = tracing::info_span!("read_summary", file = %path.display());
    let _enter = span.enter();
    // ... existing code ...
    tracing::info!(ms1 = summary.ms1_count, ms2 = summary.ms2_count, "summary read");
    Ok(summary)
}
```

- [ ] **Step 5: Verify build + tests**

Run: `cargo test -p protein-copilot-spectrum-io --lib 2>&1 | tail -5`
Expected: all 102 tests pass

- [ ] **Step 6: Commit**

```bash
git add crates/spectrum-io/
git commit -m "feat(tracing): L2 — spectrum-io library spans

Add info spans to IndexedMzMLReader::open (index_source, scan_count),
build_index_by_byte_scan, disk_cache load/save, read_summary, and
for_each_spectrum. Logs index source decision context at info level.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 4: L2 — search-engine Library Spans (15 spans)

**Files:**
- Modify: `crates/search-engine/src/simple_engine.rs`
- Modify: `crates/search-engine/src/adapters/sage/mod.rs`
- Modify: `crates/search-engine/src/fasta.rs`
- Modify: `crates/search-engine/src/digest.rs`
- Modify: `crates/search-engine/src/matching.rs`
- Modify: `crates/search-engine/src/annotate.rs`

- [ ] **Step 1: Add stage spans to SimpleSearchEngine::search()**

Wrap each stage (FASTA parse, digest, match, score) in `info_span!`. Keep existing `Instant::now()` for `SearchProgress` callback but also add span-based logging:

```rust
// Before FASTA parse:
let fasta_span = tracing::info_span!("parse_fasta",
    path = %params.database_path,
    protein_count = tracing::field::Empty,
);
{
    let _enter = fasta_span.enter();
    // ... existing parse_fasta code ...
    fasta_span.record("protein_count", proteins.len() as u64);
    tracing::info!("FASTA parsed");
}

// Before digest:
let digest_span = tracing::info_span!("digest",
    protein_count = proteins.len(),
    peptide_count = tracing::field::Empty,
    enzyme = ?params.enzyme,
);
{
    let _enter = digest_span.enter();
    // ... existing digest code ...
    digest_span.record("peptide_count", peptides.len() as u64);
    tracing::info!("digestion complete");
}

// Before matching:
let match_span = tracing::info_span!("match_spectra",
    spectrum_count = spectra.len(),
    matched = tracing::field::Empty,
);
{
    let _enter = match_span.enter();
    // ... existing matching code ...
    match_span.record("matched", matched_count as u64);
    tracing::info!("matching complete");
}
```

- [ ] **Step 2: Add spans to SageSearchEngine::search()**

Same pattern as SimpleSearch — wrap the main search call and result processing in spans. Log `engine="Sage"`, spectrum_count, and result metrics.

- [ ] **Step 3: Add spans to fasta.rs parse_fasta()**

```rust
pub fn parse_fasta(path: &Path) -> Result<Vec<FastaEntry>, SearchEngineError> {
    let span = tracing::info_span!("parse_fasta", path = %path.display());
    let _enter = span.enter();
    // ... existing code ...
    tracing::info!(proteins = entries.len(), "FASTA parsed");
    Ok(entries)
}
```

- [ ] **Step 4: Add span to digest.rs digest()**

```rust
pub fn digest(proteins: &[FastaEntry], params: &SearchParams) -> Vec<DigestedPeptide> {
    let span = tracing::info_span!("digest",
        proteins = proteins.len(),
        enzyme = ?params.enzyme,
        missed_cleavages = params.missed_cleavages,
    );
    let _enter = span.enter();
    // ... existing code ...
    tracing::info!(peptides = result.len(), "digestion complete");
    result
}
```

- [ ] **Step 5: Add spans to annotate.rs functions**

Add `#[tracing::instrument(skip(spectrum))]` to `annotate_spectrum_impl` and related public functions.

- [ ] **Step 6: Verify build + tests**

Run: `cargo test -p protein-copilot-search-engine --lib 2>&1 | tail -5`
Expected: all tests pass

- [ ] **Step 7: Commit**

```bash
git add crates/search-engine/
git commit -m "feat(tracing): L2 — search-engine library spans

Add info spans to SimpleSearch/Sage search stages (FASTA parse, digest,
match, score), parse_fasta, digest, and annotate functions. Each stage
logs input size and output metrics at info level.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 5: L2 — xic, dia-extraction, result-import, fdr, report, param-recommend, protein-inference, entrapment Spans (~25 spans)

**Files:**
- Modify: `crates/xic/src/extract.rs`
- Modify: `crates/dia-extraction/src/lib.rs`, `correlation.rs`, `detection.rs`
- Modify: `crates/result-import/src/lib.rs`, `diann.rs`, `pfind.rs`, `scan_matcher.rs`
- Modify: `crates/fdr/src/calculation.rs`
- Modify: `crates/report/src/lib.rs`
- Modify: `crates/param-recommend/src/lib.rs`
- Modify: `crates/protein-inference/src/mapper.rs`, `parsimony.rs`, `razor.rs`, `coverage.rs`
- Modify: `crates/entrapment-analysis/src/lib.rs`, `digest.rs`, `levenshtein.rs`

- [ ] **Step 1: xic crate — add spans to extract_xic_unified and helpers**

In `extract.rs`, add `info_span!` to `extract_xic_unified()`:

```rust
pub fn extract_xic_unified(...) -> Result<XicUnifiedResult, XicError> {
    let span = tracing::info_span!("extract_xic_unified",
        peptide = %peptide_sequence,
        precursor_mz = precursor_mz,
        scans_planned = tracing::field::Empty,
    );
    let _enter = span.enter();
    // ... existing code ...
    span.record("scans_planned", planned_scans.len() as u64);
    tracing::info!("XIC extraction complete");
    // ...
}
```

- [ ] **Step 2: dia-extraction crate — add spans to main entry and correlation**

In `lib.rs`, add span to `extract_dia_precursors()`:
```rust
let span = tracing::info_span!("extract_dia_precursors",
    file = %file_path.display(),
    mode = tracing::field::Empty,
    ms2_count = tracing::field::Empty,
    candidates = tracing::field::Empty,
);
```

In `detection.rs`, log detected mode:
```rust
tracing::info!(mode = %detected_mode, "acquisition mode detected");
```

In `correlation.rs`, log correlation stats:
```rust
tracing::info!(ms1 = ms1_count, ms2 = ms2_count, "MS1-MS2 correlation complete");
```

- [ ] **Step 3: result-import crate — add spans to format detection, parse, and scan matching**

In `lib.rs`:
```rust
tracing::info!(file = %path.display(), format = %detected, "format detected");
```

In `diann.rs` / `pfind.rs`:
```rust
let span = tracing::info_span!("parse_results", format = "diann_parquet", file = %path.display());
// After parse:
tracing::info!(rows = row_count, valid_psms = psm_count, "results parsed");
```

In `scan_matcher.rs`:
```rust
let span = tracing::info_span!("match_scans", psm_count = psms.len());
// After matching:
tracing::info!(matched = matched, unmatched = unmatched, "scan matching complete");
```

- [ ] **Step 4: fdr crate — add spans to calculate_fdr**

In `calculation.rs`:
```rust
pub fn calculate_fdr(psms: &mut [Psm]) -> FdrResult {
    let span = tracing::info_span!("calculate_fdr", psm_count = psms.len());
    let _enter = span.enter();
    // ... existing code ...
    tracing::info!(
        target = target_count,
        decoy = decoy_count,
        psms_at_1pct = result.psms_at_1pct_fdr,
        "FDR calculated"
    );
    result
}
```

- [ ] **Step 5: report crate — add spans to generate_summary and export functions**

In `lib.rs`:
```rust
pub fn generate_summary(...) -> SearchResultSummary {
    let span = tracing::info_span!("generate_summary");
    let _enter = span.enter();
    // ... existing code ...
    tracing::info!(
        id_rate = format!("{:.1}%", summary.identification_rate * 100.0),
        median_ppm = format!("{:.1}", summary.median_delta_mass_ppm),
        psms_1pct = summary.psms_at_1pct_fdr,
        "summary generated"
    );
    summary
}
```

- [ ] **Step 6: param-recommend crate — add span**

In `lib.rs`:
```rust
tracing::info!(
    enzyme = ?recommended.enzyme,
    precursor_tol = ?recommended.precursor_tolerance,
    confidence = confidence,
    "parameters recommended"
);
```

- [ ] **Step 7: protein-inference crate — add spans to all 4 public functions**

```rust
// mapper.rs
let span = tracing::info_span!("build_peptide_protein_map", psm_count = psms.len());

// parsimony.rs
let span = tracing::info_span!("run_parsimony", peptide_count, protein_count);

// razor.rs
let span = tracing::info_span!("assign_razor_peptides", peptide_count);

// coverage.rs
let span = tracing::info_span!("calculate_coverage", protein_count);
```

- [ ] **Step 8: entrapment-analysis crate — add spans to classify_all, digest, find_similar**

```rust
// lib.rs - classify entry
let span = tracing::info_span!("classify_entrapment",
    psm_count, trap_count = tracing::field::Empty);

// digest.rs - target digest
let span = tracing::info_span!("digest_target_fasta",
    fasta = %path.display(), proteins = tracing::field::Empty);

// levenshtein.rs - find_similar
let span = tracing::info_span!("find_similar",
    query = %peptide, max_mismatches);
```

- [ ] **Step 9: Verify full workspace build + tests**

Run: `cargo test --workspace 2>&1 | grep 'test result' | awk '{s+=$4; f+=$6} END {print s " passed, " f " failed"}'`
Expected: 894 passed, 0 failed

- [ ] **Step 10: Commit**

```bash
git add crates/xic/ crates/dia-extraction/ crates/result-import/ crates/fdr/ crates/report/ crates/param-recommend/ crates/protein-inference/ crates/entrapment-analysis/
git commit -m "feat(tracing): L2 — library spans for 8 crates

Add info spans to xic (extract_xic_unified), dia-extraction (main +
correlation + detection), result-import (format detect + parse + scan
match), fdr (calculate_fdr), report (generate_summary + export),
param-recommend (recommend_params), protein-inference (mapper +
parsimony + razor + coverage), entrapment-analysis (classify +
digest + find_similar). All log input sizes and output metrics.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 6: L3 — Hot Loop Progress (23 locations)

**Files:** Same files as Tasks 3-5, plus additional loop bodies.

The pattern for all hot loops is identical:

```rust
let total = items.len();
let progress_interval = N; // varies per loop, see table below
let start = std::time::Instant::now();

for (i, item) in items.iter().enumerate() {
    // ... existing processing ...

    if (i + 1) % progress_interval == 0 || i + 1 == total {
        let elapsed = start.elapsed().as_secs_f64();
        let rate = (i + 1) as f64 / elapsed;
        let eta = if rate > 0.0 { (total - i - 1) as f64 / rate } else { 0.0 };
        tracing::info!(
            progress = i + 1,
            total = total,
            pct = format!("{:.1}%", (i + 1) as f64 / total as f64 * 100.0),
            rate = format!("{:.0}/s", rate),
            eta_sec = format!("{:.0}", eta),
            "processing"  // message varies per loop
        );
    }
}
```

- [ ] **Step 1: search-engine hot loops**

Add progress logging to:

| Location | Loop content | Interval | Message |
|----------|-------------|----------|---------|
| `simple_engine.rs` spectrum matching loop | `for spectrum in spectra` | 500 | `"matching spectra"` |
| `simple_engine.rs` scoring loop (if separate) | PSM scoring | 5000 | `"scoring PSMs"` |
| `fasta.rs` parse loop | protein parsing | 5000 | `"parsing FASTA"` |
| `digest.rs` digest loop | protein digestion | 1000 | `"digesting proteins"` |

- [ ] **Step 2: spectrum-io hot loops**

| Location | Loop content | Interval | Message |
|----------|-------------|----------|---------|
| `mzml.rs` for_each_spectrum | streaming spectra | 1000 | `"streaming spectra"` |
| `index.rs` build_index_by_byte_scan | byte scanning | 5000 | `"scanning for spectra"` |

- [ ] **Step 3: xic, dia-extraction hot loops**

| Location | Loop content | Interval | Message |
|----------|-------------|----------|---------|
| `xic/extract.rs` spectrum read loop | O(1) seek reads | per scan | `"extracting XIC"` |
| `dia-extraction/correlation.rs` | MS1-MS2 correlation | 1000 | `"correlating MS1-MS2"` |
| `dia-extraction/lib.rs` isotope extraction | isotope patterns | 500 | `"extracting isotope patterns"` |

- [ ] **Step 4: result-import, fdr, report hot loops**

| Location | Loop content | Interval | Message |
|----------|-------------|----------|---------|
| `result-import/scan_matcher.rs` | scan matching loop | 1000 | `"matching scans"` |
| `fdr/calculation.rs` | FDR scan loop | 5000 | `"calculating FDR"` |
| `report/lib.rs` export_tsv | TSV writing | 5000 | `"exporting TSV"` |

- [ ] **Step 5: protein-inference, entrapment hot loops**

| Location | Loop content | Interval | Message |
|----------|-------------|----------|---------|
| `protein-inference/mapper.rs` | peptide-protein mapping | 5000 | `"building peptide map"` |
| `protein-inference/parsimony.rs` | parsimony iterations | per iter | `"parsimony iteration"` |
| `protein-inference/razor.rs` | razor assignment | 1000 | `"assigning razor peptides"` |
| `protein-inference/coverage.rs` | coverage calculation | 500 | `"calculating coverage"` |
| `entrapment-analysis/lib.rs` classify_all | PSM classification | 500 | `"classifying PSMs"` |
| `entrapment-analysis/digest.rs` | target digestion | 1000 | `"digesting target FASTA"` |
| `entrapment-analysis` provenance | provenance tracing | 50 | `"tracing provenance"` |

- [ ] **Step 6: Verify full workspace tests**

Run: `cargo test --workspace 2>&1 | grep 'test result' | awk '{s+=$4; f+=$6} END {print s " passed, " f " failed"}'`
Expected: 894 passed, 0 failed

- [ ] **Step 7: Commit**

```bash
git add -A
git commit -m "feat(tracing): L3 — hot loop progress for 23 locations

Add progress logging with rate/ETA to all hot loops: spectrum matching
(every 500), FASTA parsing (5000), digestion (1000), XIC extraction
(per scan), DIA correlation (1000), scan matching (1000), FDR
calculation (5000), protein inference (varies), entrapment
classification (500). All at info level with progress/total/pct/rate.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 7: Migrate Instant::now() + Add Anomaly Warnings

**Files:**
- Modify: `crates/search-engine/src/simple_engine.rs`
- Modify: `crates/search-engine/src/adapters/sage/mod.rs`
- Modify: `crates/mcp-server/src/tools.rs`

- [ ] **Step 1: Migrate manual Instant::now() in simple_engine.rs**

Keep the `start` Instant for `SearchProgress` callback (it feeds `elapsed_sec` to MCP `get_search_status`), but remove duplicate manual timing logs — the `FmtSpan::CLOSE` on spans now auto-reports duration.

Review each `Instant::now()` usage:
- `simple_engine.rs:185` — keep for `SearchProgress.elapsed_sec`, add span wrapping
- `simple_engine.rs:477` — if this is a sub-stage timer, convert to `info_span!`
- `sage/mod.rs:123` — same: keep for progress callback, add span

- [ ] **Step 2: Add anomaly warnings to search-engine**

In the spectrum matching loop (after progress logging):

```rust
// Track rate for anomaly detection
if i > progress_interval * 2 {
    let recent_rate = progress_interval as f64 / batch_elapsed;
    if recent_rate < overall_rate * 0.5 {
        tracing::warn!(
            batch_rate = format!("{:.0}/s", recent_rate),
            normal_rate = format!("{:.0}/s", overall_rate),
            scan_range = format!("{}-{}", i - progress_interval, i),
            "slow batch detected"
        );
    }
}
```

- [ ] **Step 3: Add anomaly warning for low identification rate**

In `report/lib.rs` generate_summary, after computing `identification_rate`:

```rust
if summary.identification_rate < 0.10 {
    tracing::warn!(
        id_rate = format!("{:.1}%", summary.identification_rate * 100.0),
        "low identification rate — check search parameters or database"
    );
}
```

- [ ] **Step 4: Verify**

Run: `cargo test --workspace 2>&1 | grep 'test result' | awk '{s+=$4; f+=$6} END {print s " passed, " f " failed"}'`
Expected: 894 passed, 0 failed

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(tracing): migrate Instant::now + add anomaly warnings

Convert manual timing to span-based where possible, keeping Instant for
SearchProgress callbacks. Add slow-batch detection warning in spectrum
matching and low-identification-rate warning in report generation.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 8: Verification — End-to-End Tracing Output

**Files:** None (verification only)

- [ ] **Step 1: Run full workspace tests**

```bash
cd /home/verden/pfind/2026-spring/code/proteinCopilot
cargo test --workspace 2>&1 | grep 'test result' | awk '{s+=$4; f+=$6} END {print s " passed, " f " failed"}'
```
Expected: 894 passed, 0 failed

- [ ] **Step 2: Verify tracing output with a real search**

```bash
RUST_LOG=info cargo run -p protein-copilot-mcp-server 2>tracing_output.log &
# Send a read_spectra call via MCP client, then check:
grep 'mcp_tool' tracing_output.log | head -5
grep 'index_source' tracing_output.log | head -3
```

Verify the output shows tool entry/exit, index source decisions, and stage durations.

- [ ] **Step 3: Verify JSON mode**

```bash
PROTEIN_LOG_JSON=1 RUST_LOG=info cargo run -p protein-copilot-mcp-server 2>json_output.log &
# Verify valid JSON lines:
head -3 json_output.log | python3 -m json.tool
```

- [ ] **Step 4: Verify no stdout pollution**

```bash
RUST_LOG=debug cargo run -p protein-copilot-mcp-server 1>stdout.log 2>stderr.log &
# stdout should only have MCP JSON-RPC, no tracing:
grep -c 'INFO\|DEBUG\|WARN' stdout.log  # should be 0
```

- [ ] **Step 5: Clean up temp files**

```bash
rm -f tracing_output.log json_output.log stdout.log stderr.log
```

