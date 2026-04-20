# Entrapment Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a reusable entrapment analysis library + CLI + MCP tools that classify trap PSMs by homology level (L0-L4) against a target proteome, producing structured outputs and interactive HTML reports.

**Architecture:** Two new crates: `entrapment-analysis` (lib) for core logic (config, loading, tagging, digest, similarity, reporting) and `entrapment-cli` (bin) for the CLI wrapper. Three new MCP tools added to the existing `mcp-server` crate. The library depends on `search-engine` for FASTA parsing, in-silico digest, and chemistry functions; on `core` for shared types; and adds `serde_yaml`, `arrow`, `parquet` for its own I/O needs.

**Tech Stack:** Rust, serde + serde_yaml, arrow 54 + parquet 54, Plotly.js (HTML template), clap 4 (CLI), rmcp (MCP tools)

**Spec:** `tasks/003-entrapment-analysis.md`

---

## File Structure

### New: `crates/entrapment-analysis/` (lib crate)

| File | Responsibility |
|------|---------------|
| `Cargo.toml` | Crate metadata and dependencies |
| `src/lib.rs` | Public API: `EntrapmentAnalyzer`, re-exports |
| `src/error.rs` | `EntrapmentError` enum (thiserror) |
| `src/config.rs` | `EntrapmentConfig` - YAML parsing + validation |
| `src/types.rs` | `UnifiedPsm`, `ClassifiedPsm`, `DiscriminabilityLevel`, `EntrapmentSummary` |
| `src/tagger.rs` | `Tagger` - applies target/trap rules to PSMs |
| `src/digest.rs` | Thin wrapper: digest target FASTA into HashMap<len, Vec<(seq, accession)>> |
| `src/similarity.rs` | `classify_single()` - L0/L1/L2/L3/L4 classification |
| `src/loader/mod.rs` | `ResultLoader` trait + format detection |
| `src/loader/diann_parquet.rs` | DIA-NN parquet loader |
| `src/loader/generic_tsv.rs` | Generic TSV loader |
| `src/output.rs` | Write classified PSMs to TSV + razor_errors.tsv + run_metadata.json |
| `src/report.rs` | Generate interactive HTML report (Plotly.js) |
| `templates/entrapment_report.html` | HTML template with Plotly.js |

### New: `crates/entrapment-cli/` (bin crate)

| File | Responsibility |
|------|---------------|
| `Cargo.toml` | Crate metadata, clap dependency |
| `src/main.rs` | CLI: `analyze`, `report`, `inspect` subcommands |

### Modified files

| File | Change |
|------|--------|
| `Cargo.toml` (root) | Add workspace deps for new crates + serde_yaml + clap + csv |
| `crates/mcp-server/Cargo.toml` | Add `protein-copilot-entrapment-analysis` dep |
| `crates/mcp-server/src/tools.rs` | Add 3 entrapment MCP tools |

---

## Task 1: Crate Scaffolding & Error Types

**Files:**
- Create: `crates/entrapment-analysis/Cargo.toml`
- Create: `crates/entrapment-analysis/src/lib.rs`
- Create: `crates/entrapment-analysis/src/error.rs`
- Modify: `Cargo.toml` (root workspace)

- [ ] **Step 1: Create `crates/entrapment-analysis/Cargo.toml`**

```toml
[package]
name = "protein-copilot-entrapment-analysis"
version.workspace = true
edition.workspace = true
rust-version.workspace = true
license.workspace = true
description = "Entrapment (trap-database) hit classification and homology analysis"

[dependencies]
protein-copilot-core = { workspace = true }
protein-copilot-search-engine = { workspace = true }

serde = { workspace = true }
serde_json = { workspace = true }
serde_yaml = "0.9"
thiserror = { workspace = true }
tracing = { workspace = true }
chrono = { workspace = true }
sha2 = { workspace = true }
regex = "1"
csv = "1"
arrow = { version = "54", default-features = false, features = ["prettyprint"] }
parquet = { version = "54", default-features = false, features = ["arrow", "zstd", "snap", "lz4"] }

[dev-dependencies]
tokio = { workspace = true }
tempfile = "3"
```

- [ ] **Step 2: Create `crates/entrapment-analysis/src/error.rs`**

```rust
//! Error types for the entrapment analysis crate.

use std::path::PathBuf;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum EntrapmentError {
    #[error("config error: {detail}")]
    ConfigError { detail: String },

    #[error("config file I/O error at {path}: {detail}")]
    ConfigIoError { path: PathBuf, detail: String },

    #[error("FASTA error at {path}: {detail}")]
    FastaError { path: PathBuf, detail: String },

    #[error("loader error for {format}: {detail}")]
    LoaderError { format: String, detail: String },

    #[error("I/O error at {path}: {detail}")]
    IoError { path: PathBuf, detail: String },

    #[error("output error: {detail}")]
    OutputError { detail: String },

    #[error("report error: {detail}")]
    ReportError { detail: String },
}
```

- [ ] **Step 3: Create `crates/entrapment-analysis/src/lib.rs`**

```rust
//! Entrapment analysis - classify trap-database PSM hits by homology to target proteome.
//!
//! Provides L0-L4 discriminability levels for each trap PSM, identifying
//! razor attribution errors, L/I isomers, near-identical homologs, and true trap hits.

pub mod config;
pub mod error;

pub use error::EntrapmentError;
```

- [ ] **Step 4: Add workspace deps to root `Cargo.toml`**

Add these lines in the `[workspace.dependencies]` section (after the existing internal crates around line 49):

```toml
protein-copilot-entrapment-analysis = { path = "crates/entrapment-analysis" }

# YAML config
serde_yaml = "0.9"

# CSV I/O
csv = "1"

# CLI
clap = { version = "4", features = ["derive"] }
```

- [ ] **Step 5: Verify it compiles**

Run: `cargo check -p protein-copilot-entrapment-analysis 2>&1 | tail -5`
Expected: `Finished` with no errors.

- [ ] **Step 6: Commit**

```bash
git add crates/entrapment-analysis/ Cargo.toml Cargo.lock
git commit -m "feat(entrapment): scaffold crate with error types"
```

---

## Task 2: Config & YAML Parsing

**Files:**
- Create: `crates/entrapment-analysis/src/config.rs`
- Modify: `crates/entrapment-analysis/src/lib.rs` (already exports config)

- [ ] **Step 1: Create config module with types and tests**

Create `crates/entrapment-analysis/src/config.rs` with:

- `EntrapmentConfig` struct: version, target (GroupConfig), trap (GroupConfig), conflict_resolution, unmatched, similarity (SimilarityConfig)
- `GroupConfig` struct: rules (Vec<Rule>), fasta (Vec<FastaRef>), accession_list (Option<PathBuf>)
- `Rule` enum (serde tagged): AccessionContains { any_of }, AccessionRegex { pattern }, Fasta { path }, AccessionList { path }
- `ConflictResolution` enum: PreferTarget, PreferTrap, MarkAmbiguous
- `UnmatchedPolicy` enum: Ignore, Trap, Target, Error
- `SimilarityConfig` struct: max_mismatches (default 2), delta_mz_threshold_da (default 1.0), require_tryptic_ends (default true), max_missed_cleavages (default 2)
- `EntrapmentConfig::from_yaml(path)` and `from_yaml_str(yaml)` methods
- `validate()` that checks version==1 and at least one rule per group

Tests:
- `test_parse_minimal_config` - defaults for conflict_resolution, unmatched, similarity
- `test_parse_full_config` - all fields explicit
- `test_reject_bad_version` - version != 1
- `test_reject_empty_rules` - no rules in target

See Task 2 in the detailed code blocks below.

- [ ] **Step 2: Run tests**

Run: `cargo test -p protein-copilot-entrapment-analysis 2>&1 | tail -15`
Expected: 4 tests passed.

- [ ] **Step 3: Commit**

```bash
git add crates/entrapment-analysis/src/config.rs
git commit -m "feat(entrapment): add YAML config parsing with validation"
```

---

## Task 3: Core Types - UnifiedPsm, ClassifiedPsm, DiscriminabilityLevel

**Files:**
- Create: `crates/entrapment-analysis/src/types.rs`
- Modify: `crates/entrapment-analysis/src/lib.rs`

- [ ] **Step 1: Create types module**

Create `crates/entrapment-analysis/src/types.rs` with:

- `DiscriminabilityLevel` enum: L0, L1, L2, L3, L4 with `as_str()` and `Display`
- `UnifiedPsm` struct: peptide, charge, precursor_mz, retention_time, scan_number, spectrum_file, protein_ids, q_value
- `PsmGroup` enum: Target, Trap, Ambiguous with `Display`
- `ClassifiedPsm` struct: psm (UnifiedPsm), group, level, best_target_peptide, best_target_protein, mismatches, delta_mass_da, diff_positions
- `LevelCounts` struct: l0-l4 counts with `total()` and `increment(level)`
- `EntrapmentSummary` struct: total_psms, target_psms, trap_psms, ambiguous_psms, level_counts, top_razor_families
- `RazorFamily` struct: family, count, example_peptide, example_trap_protein, example_target_protein

Tests:
- `test_level_display`
- `test_level_counts`
- `test_psm_group_display`

- [ ] **Step 2: Update lib.rs exports**

Add `pub mod types;` and re-export key types.

- [ ] **Step 3: Run tests and commit**

Run: `cargo test -p protein-copilot-entrapment-analysis 2>&1 | tail -15`
Expected: 7 tests.

```bash
git add crates/entrapment-analysis/src/types.rs crates/entrapment-analysis/src/lib.rs
git commit -m "feat(entrapment): add core types - UnifiedPsm, ClassifiedPsm, levels"
```

---

## Task 4: Tagger - Target/Trap Classification

**Files:**
- Create: `crates/entrapment-analysis/src/tagger.rs`
- Modify: `crates/entrapment-analysis/src/lib.rs`

- [ ] **Step 1: Create tagger module**

Create `crates/entrapment-analysis/src/tagger.rs` with:

- `Tagger` struct with target/trap matchers and accession sets
- `Matcher` enum: Contains(Vec<String>), Regex(Regex)
- `Tagger::new(config)` - builds matchers from rules, loads accession lists
- `Tagger::tag(protein_ids) -> Result<PsmGroup>` - semicolon-separated protein IDs, first-match wins
- Conflict resolution and unmatched policy handling
- FASTA accession loading via `parse_fasta` from search-engine crate
- Accession list file loading

Tests:
- `test_simple_contains_rules` - target/trap/unmatched
- `test_conflict_prefer_target`
- `test_conflict_mark_ambiguous`
- `test_unmatched_error`
- `test_semicolon_separated_proteins`
- `test_regex_rule`

- [ ] **Step 2: Run tests and commit**

Run: `cargo test -p protein-copilot-entrapment-analysis 2>&1 | tail -15`
Expected: 13 tests.

```bash
git add crates/entrapment-analysis/src/tagger.rs crates/entrapment-analysis/src/lib.rs
git commit -m "feat(entrapment): add tagger for target/trap classification"
```

---

## Task 5: Digest Index - Target FASTA to Peptide Lookup

**Files:**
- Create: `crates/entrapment-analysis/src/digest.rs`
- Modify: `crates/entrapment-analysis/src/lib.rs`

- [ ] **Step 1: Create digest module**

Create `crates/entrapment-analysis/src/digest.rs` with:

- `TargetPeptide` struct: sequence, protein_accession, neutral_mass
- `TargetDigestIndex` struct containing:
  - `by_length: HashMap<usize, Vec<TargetPeptide>>` for hamming scan
  - `exact_set: HashSet<String>` for L0 lookup
  - `normalized_set: HashSet<String>` for L1 lookup (L/I normalized)
  - `exact_to_protein: HashMap<String, String>` for protein accession
  - `normalized_to_original: HashMap<String, (String, String)>` for original seq + protein
- `normalize_li(seq) -> String` - replace all I with L
- `TargetDigestIndex::from_fasta(path, sim_config)` - digest with Trypsin, configurable missed cleavages
- Query methods: `has_exact()`, `exact_protein()`, `has_normalized()`, `normalized_match()`, `peptides_of_length()`
- `len()`, `is_empty()`
- `#[cfg(test)] empty_for_test()` constructor

Tests:
- `test_normalize_li` - verify I->L replacement
- `test_normalize_li_empty`

- [ ] **Step 2: Run tests and commit**

Run: `cargo test -p protein-copilot-entrapment-analysis 2>&1 | tail -15`

```bash
git add crates/entrapment-analysis/src/digest.rs crates/entrapment-analysis/src/lib.rs
git commit -m "feat(entrapment): add target FASTA digest index"
```

---

## Task 6: Similarity Classification - L0/L1/L2/L3/L4

**Files:**
- Create: `crates/entrapment-analysis/src/similarity.rs`
- Modify: `crates/entrapment-analysis/src/lib.rs`

- [ ] **Step 1: Create similarity module**

Create `crates/entrapment-analysis/src/similarity.rs` with:

- `hamming_diff(a, b) -> Option<(u8, f64, String)>` - hamming distance + delta mass + diff positions string
  - Returns None if lengths differ
  - Uses `residue_mass()` from search-engine for mass delta calculation
  - Diff format: `[pos:X->Y,pos2:A->B]`
- `classify_single(psm, group, index, config) -> ClassifiedPsm`:
  1. If not Trap group, return L4 immediately
  2. L0: `index.has_exact(seq)` -> return L0
  3. L1: `index.has_normalized(seq)` and NOT exact -> return L1
  4. L2/L3/L4: brute-force scan `index.peptides_of_length(len)`:
     - Skip if mismatches > max_mismatches
     - Skip if all differences are only L/I (already handled by L1)
     - Track best match (lowest mismatches, then lowest abs delta_mass)
     - L2: 1 mismatch AND abs(delta_mass) < threshold
     - L3: 1-2 mismatches AND not L2
  5. No match within threshold -> L4

Tests:
- `test_hamming_identical`
- `test_hamming_one_mismatch` - verify diff string format
- `test_hamming_two_mismatches`
- `test_hamming_different_length` - returns None
- `test_hamming_mass_difference` - D vs N substitution, verify delta ~0.984 Da
- `test_classify_target_psm_gets_l4`

- [ ] **Step 2: Run tests and commit**

Run: `cargo test -p protein-copilot-entrapment-analysis 2>&1 | tail -15`

```bash
git add crates/entrapment-analysis/src/similarity.rs crates/entrapment-analysis/src/digest.rs crates/entrapment-analysis/src/lib.rs
git commit -m "feat(entrapment): add L0-L4 similarity classification"
```

---

## Task 7: Result Loaders - DIA-NN Parquet & Generic TSV

**Files:**
- Create: `crates/entrapment-analysis/src/loader/mod.rs`
- Create: `crates/entrapment-analysis/src/loader/diann_parquet.rs`
- Create: `crates/entrapment-analysis/src/loader/generic_tsv.rs`
- Modify: `crates/entrapment-analysis/src/lib.rs`

- [ ] **Step 1: Create loader module**

`src/loader/mod.rs`:
- `ResultFormat` enum: DiannParquet, GenericTsv
- `ResultFormat::from_path(path)` - detect from extension (.parquet, .tsv/.txt)
- `load_psms(path, format, tsv_column_map)` - dispatch to format-specific loader

`src/loader/diann_parquet.rs`:
- Load columns: Stripped.Sequence, Precursor.Charge, Precursor.Mz, RT, Q.Value, Run, Protein.Ids
- Uses arrow Float64Array/Float32Array/Int32Array/StringArray
- Pattern follows `crates/result-import/src/diann.rs`

`src/loader/generic_tsv.rs`:
- `TsvColumnMap` struct with serde defaults for column names
- Uses `csv` crate with configurable delimiter
- Optional columns (scan_number, spectrum_file) gracefully handled

Tests:
- `test_load_tsv` - create temp TSV, load, verify fields

- [ ] **Step 2: Run tests and commit**

Run: `cargo test -p protein-copilot-entrapment-analysis 2>&1 | tail -15`

```bash
git add crates/entrapment-analysis/src/loader/
git commit -m "feat(entrapment): add DIA-NN parquet and generic TSV loaders"
```

---

## Task 8: EntrapmentAnalyzer - Public API & Output

**Files:**
- Create: `crates/entrapment-analysis/src/output.rs`
- Modify: `crates/entrapment-analysis/src/lib.rs`

- [ ] **Step 1: Create output module**

`src/output.rs`:
- `RunMetadata` struct (Serialize): tool_version, run_timestamp, input/fasta file + SHA256, config snapshot, total/trap PSMs, level_counts
- `write_classified_tsv(psms, path)` - TSV with all ClassifiedPsm fields
- `write_razor_errors_tsv(psms, path)` - only L0 trap PSMs: peptide, current_razor, suggested_razor, reason
- `write_run_metadata(metadata, path)` - pretty-printed JSON
- `file_sha256(path)` - compute SHA-256 hex digest

- [ ] **Step 2: Implement EntrapmentAnalyzer in lib.rs**

Update `src/lib.rs` to add:
- `EntrapmentAnalyzer` struct: config, tagger, index
- `EntrapmentAnalyzer::new(config, fasta_path)` - builds tagger + digest index
- `classify(psm) -> Result<ClassifiedPsm>` - tag + classify_single
- `classify_all(psms) -> Result<Vec<ClassifiedPsm>>`
- `summary(classified) -> EntrapmentSummary` - counts by level/group, top razor families
- Helper `extract_family_name(accession)` - parse UniProt format

Test: `test_extract_family_name`

- [ ] **Step 3: Run tests and commit**

Run: `cargo test -p protein-copilot-entrapment-analysis 2>&1 | tail -15`

```bash
git add crates/entrapment-analysis/src/output.rs crates/entrapment-analysis/src/lib.rs
git commit -m "feat(entrapment): add EntrapmentAnalyzer API and output writers"
```

---

## Task 9: CLI - `entrapment-cli` Crate

**Files:**
- Create: `crates/entrapment-cli/Cargo.toml`
- Create: `crates/entrapment-cli/src/main.rs`
- Modify: `Cargo.toml` (root)

- [ ] **Step 1: Create CLI crate**

`Cargo.toml`:
- bin name: `entrapment`
- deps: entrapment-analysis, clap (derive), serde_json, tracing, tracing-subscriber, chrono

`src/main.rs` with clap derive:
- `Cli` with `Commands` subcommand enum
- `Commands::Analyze` - results, format, config, target-fasta, out
- `Commands::Report` - classified (TSV/Parquet path), out (HTML path)
- `Commands::Inspect` - peptide, target-fasta, config (optional)
- `run_analyze()` - load config, detect format, load PSMs, build analyzer, classify, write outputs, print summary
- `run_report()` - read classified TSV, regenerate HTML report
- `run_inspect()` - build index, classify dummy PSM, print result

- [ ] **Step 2: Add workspace dep and verify build**

Add to root `Cargo.toml`:
```toml
protein-copilot-entrapment-cli = { path = "crates/entrapment-cli" }
```

Run: `cargo build -p protein-copilot-entrapment-cli 2>&1 | tail -5`
Expected: `Finished`

- [ ] **Step 3: Commit**

```bash
git add crates/entrapment-cli/ Cargo.toml Cargo.lock
git commit -m "feat(entrapment): add CLI with analyze, report, and inspect subcommands"
```

---

## Task 10: HTML Report - Interactive Plotly Dashboard

**Files:**
- Create: `crates/entrapment-analysis/templates/entrapment_report.html`
- Create: `crates/entrapment-analysis/src/report.rs`
- Modify: `crates/entrapment-analysis/src/lib.rs`
- Modify: `crates/entrapment-cli/src/main.rs`

- [ ] **Step 1: Create HTML template**

Create `crates/entrapment-analysis/templates/entrapment_report.html`:
- Self-contained HTML with Plotly.js CDN
- `/*__REPORT_DATA__*/` placeholder replaced with JSON
- Sections: summary stats panel, pie chart (L0-L4), bar chart (razor families), delta-mass histogram (L2/L3), searchable PSM table
- ~300 lines of HTML + JS

- [ ] **Step 2: Create report module**

`src/report.rs`:
- `ReportData` struct (Serialize): summary + psms (flattened PsmRow)
- `PsmRow` struct - flat version of ClassifiedPsm for JS consumption
- `render_report(summary, classified, output_path)` - serialize to JSON, escape for HTML, replace placeholder, write file
- Pattern follows `crates/report/src/visualize.rs`

- [ ] **Step 3: Add report generation to CLI analyze**

After `write_run_metadata` in `run_analyze()`, add:
```rust
protein_copilot_entrapment_analysis::report::render_report(&summary, &classified, &out_dir.join("entrapment_report.html"))?;
```

- [ ] **Step 4: Verify and commit**

Run: `cargo check -p protein-copilot-entrapment-analysis -p protein-copilot-entrapment-cli 2>&1 | tail -5`

```bash
git add crates/entrapment-analysis/src/report.rs crates/entrapment-analysis/templates/ crates/entrapment-cli/src/main.rs
git commit -m "feat(entrapment): add interactive HTML report with Plotly.js"
```

---

## Task 11: MCP Tools - Integrate into mcp-server

**Files:**
- Modify: `crates/mcp-server/Cargo.toml`
- Modify: `crates/mcp-server/src/tools.rs`

- [ ] **Step 1: Add entrapment dependency**

In `crates/mcp-server/Cargo.toml`, add:
```toml
protein-copilot-entrapment-analysis = { workspace = true }
```

- [ ] **Step 2: Add MCP tool input types and implementations**

In `crates/mcp-server/src/tools.rs`:

Input types:
- `ClassifyEntrapmentHitsInput`: results_file, format?, config_file, target_fasta, output_dir?
- `AnalyzeEntrapmentStatsInput`: classified_file (path to classified TSV)
- `FindSimilarTargetsInput`: peptide, target_fasta, max_mismatches?

Tool implementations in `#[rmcp::tool_router]` block:
- `classify_entrapment_hits` - full pipeline: load config + PSMs + analyzer, classify all, write outputs, return summary
- `analyze_entrapment_stats` - read classified TSV, compute summary stats (level distribution, protein family clusters, delta-mass histogram data), return DetailedStats
- `find_similar_targets` - build index, classify single peptide, return result

- [ ] **Step 3: Verify and commit**

Run: `cargo check -p protein-copilot-mcp-server 2>&1 | tail -5`

```bash
git add crates/mcp-server/Cargo.toml crates/mcp-server/src/tools.rs
git commit -m "feat(entrapment): add MCP tools - classify, analyze_stats, find_similar"
```

---

## Task 12: Integration Test with Real Data

**Files:**
- Create: `tests/entrapment_integration.rs`

- [ ] **Step 1: Write integration test**

Test the 6 known peptides from the analysis:

| Peptide | Expected Level | Reason |
|---------|---------------|--------|
| STTTGHLIYK | L0 | Exact match in EF1A1_HUMAN |
| GYSFTTTAER | L0 | Exact match in ACTB_HUMAN |
| ELTALAPSTMK | L1 | L/I isomer of EITALAPSTMK |
| HPFPGPGIAIR | L1 | L/I isomer of HPFPGPGLAIR |
| DGFLLDGFPR | L2 | 1mm from NGFLLDGFPR, delta ~0.98 Da |
| IGSEVYHNLK | L3 | 1mm from IGAEVYHNLK, delta ~16 Da |

Test gracefully skips if `human_swissprot.fasta` is not present.

- [ ] **Step 2: Run and commit**

Run: `cargo test --test entrapment_integration 2>&1 | tail -15`

```bash
git add tests/entrapment_integration.rs
git commit -m "test(entrapment): add integration test with 6 known peptides"
```

---

## Task 13: Full Build & Test Verification

- [ ] **Step 1: Full workspace build**

Run: `cargo build --workspace 2>&1 | tail -10`
Expected: `Finished`

- [ ] **Step 2: Run all tests**

Run: `cargo test --workspace 2>&1 | tail -20`
Expected: All pass.

- [ ] **Step 3: Run clippy**

Run: `cargo clippy --workspace -- -D warnings 2>&1 | tail -10`
Expected: No warnings.

- [ ] **Step 4: Run rustfmt**

Run: `cargo fmt --all -- --check 2>&1 | tail -10`
Expected: No issues.

- [ ] **Step 5: Final fixup commit if needed**

```bash
git add -A
git commit -m "chore(entrapment): fix clippy warnings and formatting"
```
