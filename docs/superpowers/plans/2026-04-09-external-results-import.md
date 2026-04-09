# External Search Results Import — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create a `result-import` lib crate that parses external search results (custom JSON, DIA-NN parquet, pFind skeleton), matches PSMs to mzML scans via RT + isolation window, converts to standard `SearchResult`, and integrates with the existing MCP tool ecosystem via a new `import_search_results` tool.

**Architecture:** Three-layer pipeline — format parsers → scan matcher → converter. The `result-import` crate is a pure library; the MCP server wires it into the tool router and run cache. All RT values are normalized to seconds at parse time.

**Tech Stack:** Rust, `arrow`/`parquet` (DIA-NN), `quick-xml` (unimod.xml), `regex` (Modified.Sequence parsing), `serde_json` (custom JSON), existing `protein-copilot-core` types.

**Design Spec:** `docs/superpowers/specs/2026-04-09-external-results-import-design.md`

---

## File Structure

```
crates/result-import/
  Cargo.toml                     ← New lib crate
  src/
    lib.rs                       ← Public types (ImportedPsm, ImportResult, MatchReport, etc.), ResultParser trait, import() orchestrator
    error.rs                     ← ResultImportError enum (thiserror)
    unimod.rs                    ← UnimodDb: builtin table + XML parser + to_modification()
    custom_json.rs               ← CustomJsonParser: hela.json format
    diann.rs                     ← DiannParser: report.parquet format
    pfind.rs                     ← PFindParser: skeleton with ResultParser trait impl
    scan_matcher.rs              ← ScanMatcher: RT + isolation window matching
    converter.rs                 ← ImportedPsm → core::Psm → SearchResult builder

crates/mcp-server/
  Cargo.toml                     ← Add result-import dependency
  src/tools.rs                   ← Add import_search_results tool (~120 lines)

Cargo.toml                       ← Add result-import to workspace members + dependencies
```

---

### Task 1: Crate Scaffold + Core Types

**Files:**
- Create: `crates/result-import/Cargo.toml`
- Create: `crates/result-import/src/lib.rs`
- Create: `crates/result-import/src/error.rs`
- Modify: `Cargo.toml` (workspace root)

- [ ] **Step 1: Create `Cargo.toml` for result-import crate**

Create `crates/result-import/Cargo.toml`:

```toml
[package]
name = "protein-copilot-result-import"
version.workspace = true
edition.workspace = true
rust-version.workspace = true
license.workspace = true
description = "Import external proteomics search results (DIA-NN, custom JSON, pFind)"

[dependencies]
protein-copilot-core = { workspace = true }
protein-copilot-spectrum-io = { workspace = true }

serde = { workspace = true }
serde_json = { workspace = true }
thiserror = { workspace = true }
tracing = { workspace = true }
uuid = { workspace = true }
chrono = { workspace = true }
quick-xml = { workspace = true }
regex = "1"
arrow = { version = "54", default-features = false, features = ["prettyprint"] }
parquet = { version = "54", default-features = false, features = ["arrow"] }

[dev-dependencies]
tokio = { workspace = true }
```

- [ ] **Step 2: Create `error.rs`**

Create `crates/result-import/src/error.rs`:

```rust
//! Error types for the result-import crate.

use std::path::PathBuf;

/// Errors that can occur during external result import.
#[derive(Debug, thiserror::Error)]
pub enum ResultImportError {
    #[error("file not found: {path}")]
    FileNotFound { path: PathBuf },

    #[error("format detection failed for {path}: expected one of .json, .parquet, .spectra")]
    FormatDetectionFailed { path: PathBuf },

    #[error("unknown Unimod record_id: {0} — check unimod.xml or use builtin database")]
    UnknownUnimodId(u32),

    #[error("mzML file not found for raw name '{raw_name}' in directory {dir} — available files: {available}")]
    MzmlNotFound {
        raw_name: String,
        dir: PathBuf,
        available: String,
    },

    #[error("missing required column '{column}' in parquet file — expected columns: {expected}")]
    MissingColumn { column: String, expected: String },

    #[error("JSON parse error: {0}")]
    JsonError(#[from] serde_json::Error),

    #[error("XML parse error: {0}")]
    XmlError(#[from] quick_xml::Error),

    #[error("Parquet read error: {0}")]
    ParquetError(#[from] parquet::errors::ParquetError),

    #[error("Arrow error: {0}")]
    ArrowError(#[from] arrow::error::ArrowError),

    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),

    #[error("invalid modification position {position} for sequence of length {seq_len}")]
    InvalidModPosition { position: usize, seq_len: usize },

    #[error("spectrum-io error: {0}")]
    SpectrumIo(String),

    #[error("{0}")]
    Other(String),
}
```

- [ ] **Step 3: Create `lib.rs` with core types**

Create `crates/result-import/src/lib.rs`:

```rust
//! Import external proteomics search results and match to mzML scans.
//!
//! Supports:
//! - Custom JSON (hela.json format)
//! - DIA-NN report.parquet
//! - pFind .spectra (skeleton)

pub mod error;
pub mod unimod;
pub mod custom_json;
pub mod diann;
pub mod pfind;
pub mod scan_matcher;
pub mod converter;

use std::collections::HashMap;
use std::path::Path;

use protein_copilot_core::search_params::Modification;
use serde::{Deserialize, Serialize};

pub use error::ResultImportError;

/// A PSM imported from an external search result file.
///
/// RT is always in seconds (converted from source format at parse time).
/// `matched_scan` is `None` until scan matching is performed.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ImportedPsm {
    pub sequence: String,
    pub charge: i32,
    pub precursor_mz: f64,
    /// Retention time in seconds (converted from source format at parse time).
    pub rt_sec: f64,
    pub modifications: Vec<Modification>,
    /// Search engine score (e.g. DIA-NN Q.Value). `None` if not available.
    pub score: Option<f64>,
    /// q-value for FDR control. `None` if not available.
    pub q_value: Option<f64>,
    pub protein_accessions: Vec<String>,
    /// Raw file name without extension (used for mzML association).
    pub raw_name: String,
    /// Filled by scan matcher.
    pub matched_scan: Option<u32>,
    /// RT delta in seconds between PSM and matched MS2 scan.
    pub rt_delta_sec: Option<f64>,
}

/// Result of the import operation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ImportResult {
    pub run_id: String,
    pub match_report: MatchReport,
    pub imported_psm_count: usize,
    pub unique_peptides: usize,
    pub protein_count: usize,
    pub raw_files: Vec<String>,
}

/// Scan matching quality report.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MatchReport {
    pub total_psms: usize,
    pub matched: usize,
    pub unmatched: usize,
    pub median_rt_delta_sec: f64,
    pub max_rt_delta_sec: f64,
    pub per_file: HashMap<String, FileMatchStats>,
}

/// Per-file matching statistics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileMatchStats {
    pub total: usize,
    pub matched: usize,
    pub ms2_count: usize,
}

/// Trait for format-specific parsers.
pub trait ResultParser: Send + Sync {
    fn parse(
        &self,
        path: &Path,
        unimod: &unimod::UnimodDb,
    ) -> Result<Vec<ImportedPsm>, ResultImportError>;
}

/// Supported import formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImportFormat {
    CustomJson,
    DiannParquet,
    PFindSpectra,
}

/// Detect format from file extension.
pub fn detect_format(path: &Path) -> Result<ImportFormat, ResultImportError> {
    match path.extension().and_then(|e| e.to_str()) {
        Some("json") => Ok(ImportFormat::CustomJson),
        Some("parquet") => Ok(ImportFormat::DiannParquet),
        Some("spectra") => Ok(ImportFormat::PFindSpectra),
        _ => Err(ResultImportError::FormatDetectionFailed {
            path: path.to_path_buf(),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn detect_format_json() {
        assert_eq!(
            detect_format(Path::new("test.json")).unwrap(),
            ImportFormat::CustomJson
        );
    }

    #[test]
    fn detect_format_parquet() {
        assert_eq!(
            detect_format(Path::new("report.parquet")).unwrap(),
            ImportFormat::DiannParquet
        );
    }

    #[test]
    fn detect_format_spectra() {
        assert_eq!(
            detect_format(Path::new("result.spectra")).unwrap(),
            ImportFormat::PFindSpectra
        );
    }

    #[test]
    fn detect_format_unknown() {
        assert!(detect_format(Path::new("data.csv")).is_err());
    }
}
```

- [ ] **Step 4: Add to workspace**

Add to root `Cargo.toml` workspace dependencies section:

```toml
protein-copilot-result-import = { path = "crates/result-import" }
```

- [ ] **Step 5: Run `cargo check -p protein-copilot-result-import` and `cargo test -p protein-copilot-result-import`**

Expected: 4 tests pass (detect_format tests).

- [ ] **Step 6: Commit**

```bash
git add -A && git commit -m "feat(result-import): scaffold crate with core types and error handling

Add ImportedPsm, ImportResult, MatchReport, ResultParser trait,
detect_format(), and ResultImportError enum."
```

---

### Task 2: UnimodDb — Builtin Table + XML Parser

**Files:**
- Create: `crates/result-import/src/unimod.rs`

**Reference:** Unimod XML at `/home/verden/pfind/2025-fall/code/ms2-met/unimod.xml` — key structure: `<umod:mod title="..." record_id="..."><umod:delta mono_mass="..."/><umod:specificity site="..." .../>`.

- [ ] **Step 1: Write tests for UnimodDb**

Add to `crates/result-import/src/unimod.rs`:

```rust
//! Unimod modification database: maps record_id → name + mono_mass.
//!
//! Provides a builtin table of ~20 common modifications and an XML parser
//! for the full unimod.xml database.

use std::collections::HashMap;
use std::path::Path;

use protein_copilot_core::search_params::{ModPosition, Modification};

use crate::ResultImportError;

/// A single Unimod modification entry.
#[derive(Debug, Clone)]
pub struct UnimodEntry {
    pub record_id: u32,
    pub title: String,
    pub mono_mass: f64,
    /// Residues this modification can occur on (empty = any).
    pub residues: Vec<char>,
}

/// Unimod modification database.
pub struct UnimodDb {
    entries: HashMap<u32, UnimodEntry>,
}

impl UnimodDb {
    /// Create a database with ~20 common builtin modifications.
    pub fn builtin() -> Self {
        let mut entries = HashMap::new();
        let common = vec![
            (1, "Acetyl", 42.010565, vec![]),
            (4, "Carbamidomethyl", 57.021464, vec!['C']),
            (5, "Carbamyl", 43.005814, vec![]),
            (7, "Deamidated", 0.984016, vec!['N', 'Q']),
            (21, "Phospho", 79.966331, vec!['S', 'T', 'Y']),
            (28, "Gln->pyro-Glu", -17.026549, vec!['Q']),
            (27, "Glu->pyro-Glu", -18.010565, vec!['E']),
            (34, "Methyl", 14.015650, vec![]),
            (35, "Oxidation", 15.994915, vec!['M', 'W', 'H']),
            (36, "Dimethyl", 28.031300, vec![]),
            (37, "Trimethyl", 42.046950, vec![]),
            (39, "Dehydrated", -18.010565, vec![]),
            (40, "Formyl", 27.994915, vec![]),
            (121, "GlyGly", 114.042927, vec!['K']),
            (188, "Label:13C(6)", 6.020129, vec!['K', 'R']),
            (199, "Label:13C(6)15N(2)", 8.014199, vec!['K']),
            (259, "Label:13C(6)15N(4)", 10.008269, vec!['R']),
            (214, "Label:2H(4)", 4.025107, vec![]),
            (267, "Label:13C(6)15N(1)", 7.017165, vec![]),
            (268, "iTRAQ4plex", 144.102063, vec![]),
            (737, "TMT6plex", 229.162932, vec![]),
            (738, "TMTpro", 304.207146, vec![]),
        ];
        for (id, title, mass, residues) in common {
            entries.insert(
                id,
                UnimodEntry {
                    record_id: id,
                    title: title.to_string(),
                    mono_mass: mass,
                    residues,
                },
            );
        }
        Self { entries }
    }

    /// Parse the full Unimod XML database.
    pub fn from_xml(path: &Path) -> Result<Self, ResultImportError> {
        use quick_xml::events::Event;
        use quick_xml::Reader;

        if !path.exists() {
            return Err(ResultImportError::FileNotFound {
                path: path.to_path_buf(),
            });
        }

        let xml_bytes = std::fs::read(path)?;
        let mut reader = Reader::from_reader(xml_bytes.as_slice());
        reader.config_mut().trim_text(true);

        let mut entries = HashMap::new();
        let mut buf = Vec::new();

        // State for current <umod:mod> element
        let mut current_id: Option<u32> = None;
        let mut current_title: Option<String> = None;
        let mut current_mass: Option<f64> = None;
        let mut current_residues: Vec<char> = Vec::new();

        loop {
            match reader.read_event_into(&mut buf) {
                Ok(Event::Start(ref e) | Event::Empty(ref e)) => {
                    let local = e.local_name();
                    let local_str = std::str::from_utf8(local.as_ref()).unwrap_or("");

                    if local_str == "mod" {
                        // Reset state
                        current_id = None;
                        current_title = None;
                        current_mass = None;
                        current_residues.clear();

                        for attr in e.attributes().flatten() {
                            let key = std::str::from_utf8(attr.key.as_ref()).unwrap_or("");
                            let val = std::str::from_utf8(&attr.value).unwrap_or("");
                            match key {
                                "record_id" => current_id = val.parse().ok(),
                                "title" => current_title = Some(val.to_string()),
                                _ => {}
                            }
                        }
                    } else if local_str == "delta" {
                        for attr in e.attributes().flatten() {
                            let key = std::str::from_utf8(attr.key.as_ref()).unwrap_or("");
                            let val = std::str::from_utf8(&attr.value).unwrap_or("");
                            if key == "mono_mass" {
                                current_mass = val.parse().ok();
                            }
                        }
                    } else if local_str == "specificity" {
                        for attr in e.attributes().flatten() {
                            let key = std::str::from_utf8(attr.key.as_ref()).unwrap_or("");
                            let val = std::str::from_utf8(&attr.value).unwrap_or("");
                            if key == "site" && val.len() == 1 {
                                let ch = val.chars().next().unwrap();
                                if ch.is_ascii_uppercase()
                                    && !current_residues.contains(&ch)
                                {
                                    current_residues.push(ch);
                                }
                            }
                        }
                    }
                }
                Ok(Event::End(ref e)) => {
                    let local = std::str::from_utf8(e.local_name().as_ref()).unwrap_or("");
                    if local == "mod" {
                        if let (Some(id), Some(title), Some(mass)) =
                            (current_id, current_title.take(), current_mass)
                        {
                            entries.insert(
                                id,
                                UnimodEntry {
                                    record_id: id,
                                    title,
                                    mono_mass: mass,
                                    residues: current_residues.clone(),
                                },
                            );
                        }
                    }
                }
                Ok(Event::Eof) => break,
                Err(e) => return Err(ResultImportError::XmlError(e)),
                _ => {}
            }
            buf.clear();
        }

        tracing::info!("loaded {} Unimod modifications from {}", entries.len(), path.display());
        Ok(Self { entries })
    }

    /// Look up a modification by Unimod record_id.
    pub fn get(&self, record_id: u32) -> Option<&UnimodEntry> {
        self.entries.get(&record_id)
    }

    /// Convert a Unimod record_id at a given position to a core::Modification.
    ///
    /// `position` is 1-based. The residue at that position in `sequence` is
    /// used to set the `residues` field of the resulting Modification.
    pub fn to_modification(
        &self,
        record_id: u32,
        position: usize,
        sequence: &str,
    ) -> Result<Modification, ResultImportError> {
        let entry = self
            .entries
            .get(&record_id)
            .ok_or(ResultImportError::UnknownUnimodId(record_id))?;

        let residue = if position >= 1 && position <= sequence.len() {
            sequence.chars().nth(position - 1).unwrap_or('X')
        } else {
            return Err(ResultImportError::InvalidModPosition {
                position,
                seq_len: sequence.len(),
            });
        };

        Ok(Modification {
            name: entry.title.clone(),
            mass_delta: entry.mono_mass,
            residues: vec![residue],
            position: ModPosition::Anywhere,
        })
    }

    /// Number of entries in the database.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the database is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builtin_has_common_mods() {
        let db = UnimodDb::builtin();
        assert!(db.len() >= 20);

        // Oxidation
        let ox = db.get(35).expect("Oxidation should be builtin");
        assert_eq!(ox.title, "Oxidation");
        assert!((ox.mono_mass - 15.994915).abs() < 0.001);

        // Carbamidomethyl
        let cam = db.get(4).expect("Carbamidomethyl should be builtin");
        assert_eq!(cam.title, "Carbamidomethyl");
        assert!((cam.mono_mass - 57.021464).abs() < 0.001);

        // Phospho
        let phos = db.get(21).expect("Phospho should be builtin");
        assert_eq!(phos.title, "Phospho");
        assert!((phos.mono_mass - 79.966331).abs() < 0.001);
    }

    #[test]
    fn to_modification_oxidation_on_m() {
        let db = UnimodDb::builtin();
        let m = db.to_modification(35, 5, "PEPTMIDE").unwrap();
        assert_eq!(m.name, "Oxidation");
        assert!((m.mass_delta - 15.994915).abs() < 0.001);
        assert_eq!(m.residues, vec!['M']);
        assert_eq!(m.position, ModPosition::Anywhere);
    }

    #[test]
    fn to_modification_unknown_id_errors() {
        let db = UnimodDb::builtin();
        assert!(db.to_modification(99999, 1, "PEPTIDE").is_err());
    }

    #[test]
    fn to_modification_invalid_position_errors() {
        let db = UnimodDb::builtin();
        assert!(db.to_modification(35, 0, "PEPTIDE").is_err());
        assert!(db.to_modification(35, 100, "PEPTIDE").is_err());
    }

    #[test]
    fn from_xml_loads_real_unimod() {
        let xml_path = Path::new("/home/verden/pfind/2025-fall/code/ms2-met/unimod.xml");
        if !xml_path.exists() {
            eprintln!("skipping XML test: unimod.xml not found");
            return;
        }
        let db = UnimodDb::from_xml(xml_path).unwrap();
        assert!(db.len() > 1000, "expected >1000 mods, got {}", db.len());

        let ox = db.get(35).expect("Oxidation");
        assert_eq!(ox.title, "Oxidation");
        assert!((ox.mono_mass - 15.994915).abs() < 0.001);
    }
}
```

- [ ] **Step 2: Run tests**

Run: `cargo test -p protein-copilot-result-import -- unimod`

Expected: 5 tests pass (4 unit + 1 XML if file exists).

- [ ] **Step 3: Commit**

```bash
git add -A && git commit -m "feat(result-import): add UnimodDb with builtin table and XML parser

22 common modifications hardcoded, full unimod.xml parser via quick-xml.
to_modification() maps record_id + position → core::Modification."
```

---

### Task 3: CustomJsonParser (hela.json)

**Files:**
- Create: `crates/result-import/src/custom_json.rs`

**Reference:** hela.json format — array of `{ sequence, charge, modify, rt, precursor_mz, raw_title, protein_names }`. `modify` is `[[position, unimod_id]]` (1-based position). `rt` is in minutes.

- [ ] **Step 1: Create test fixture**

Create `crates/result-import/src/custom_json.rs`:

```rust
//! Parser for custom JSON format (hela.json style).
//!
//! Input format: JSON array of objects with fields:
//! - `sequence`: peptide sequence
//! - `charge`: charge state
//! - `modify`: `[[position, unimod_id], ...]` (1-based position)
//! - `rt`: retention time in minutes
//! - `precursor_mz`: precursor m/z
//! - `raw_title`: raw file name (without extension)
//! - `protein_names`: array of protein accessions

use std::path::Path;

use serde::Deserialize;

use crate::unimod::UnimodDb;
use crate::{ImportedPsm, ResultImportError, ResultParser};

/// Raw JSON record matching the hela.json schema.
#[derive(Debug, Deserialize)]
struct RawJsonPsm {
    sequence: String,
    charge: i32,
    /// Modifications as [[position, unimod_id], ...]. May be empty or absent.
    #[serde(default)]
    modify: Vec<Vec<u32>>,
    /// Retention time in minutes.
    rt: f64,
    precursor_mz: f64,
    raw_title: String,
    #[serde(default)]
    protein_names: Vec<String>,
}

/// Parser for the custom JSON format.
pub struct CustomJsonParser;

impl ResultParser for CustomJsonParser {
    fn parse(
        &self,
        path: &Path,
        unimod: &UnimodDb,
    ) -> Result<Vec<ImportedPsm>, ResultImportError> {
        if !path.exists() {
            return Err(ResultImportError::FileNotFound {
                path: path.to_path_buf(),
            });
        }

        let data = std::fs::read_to_string(path)?;
        let raw_psms: Vec<RawJsonPsm> = serde_json::from_str(&data)?;

        let mut psms = Vec::with_capacity(raw_psms.len());
        let mut mod_errors = 0u64;

        for (i, raw) in raw_psms.into_iter().enumerate() {
            let mut modifications = Vec::new();
            for pair in &raw.modify {
                if pair.len() != 2 {
                    tracing::warn!("PSM #{i}: invalid modify entry (expected [pos, id]): {pair:?}");
                    continue;
                }
                let position = pair[0] as usize;
                let unimod_id = pair[1];
                match unimod.to_modification(unimod_id, position, &raw.sequence) {
                    Ok(m) => modifications.push(m),
                    Err(e) => {
                        mod_errors += 1;
                        if mod_errors <= 5 {
                            tracing::warn!(
                                "PSM #{i} ({}): modification error: {e}",
                                raw.sequence
                            );
                        }
                    }
                }
            }

            psms.push(ImportedPsm {
                sequence: raw.sequence,
                charge: raw.charge,
                precursor_mz: raw.precursor_mz,
                rt_sec: raw.rt * 60.0, // minutes → seconds
                modifications,
                score: None,
                q_value: None,
                protein_accessions: raw.protein_names,
                raw_name: raw.raw_title,
                matched_scan: None,
                rt_delta_sec: None,
            });
        }

        if mod_errors > 5 {
            tracing::warn!("... and {} more modification errors suppressed", mod_errors - 5);
        }
        tracing::info!(
            "parsed {} PSMs from custom JSON: {}",
            psms.len(),
            path.display()
        );
        Ok(psms)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn sample_json() -> &'static str {
        r#"[
            {
                "sequence": "AADLLDDVSQK",
                "charge": 2,
                "modify": [],
                "rt": 12.648,
                "precursor_mz": 588.2993,
                "raw_title": "hela_Rep1",
                "protein_names": ["sp|P12345|TEST_HUMAN"]
            },
            {
                "sequence": "PEPTMIDER",
                "charge": 3,
                "modify": [[5, 35]],
                "rt": 25.5,
                "precursor_mz": 372.185,
                "raw_title": "hela_Rep1",
                "protein_names": ["sp|P67890|TEST2_HUMAN"]
            },
            {
                "sequence": "ACDEFGHIK",
                "charge": 2,
                "modify": [[2, 4]],
                "rt": 40.0,
                "precursor_mz": 510.234,
                "raw_title": "hela_Rep2",
                "protein_names": ["sp|P11111|TEST3_HUMAN"]
            }
        ]"#
    }

    #[test]
    fn parse_custom_json_basic() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.json");
        std::fs::write(&path, sample_json()).unwrap();

        let db = UnimodDb::builtin();
        let parser = CustomJsonParser;
        let psms = parser.parse(&path, &db).unwrap();

        assert_eq!(psms.len(), 3);

        // First PSM: no modifications, RT converted to seconds
        assert_eq!(psms[0].sequence, "AADLLDDVSQK");
        assert_eq!(psms[0].charge, 2);
        assert!((psms[0].rt_sec - 12.648 * 60.0).abs() < 0.01);
        assert!(psms[0].modifications.is_empty());
        assert_eq!(psms[0].raw_name, "hela_Rep1");

        // Second PSM: Oxidation at position 5 (M)
        assert_eq!(psms[1].sequence, "PEPTMIDER");
        assert_eq!(psms[1].modifications.len(), 1);
        assert_eq!(psms[1].modifications[0].name, "Oxidation");
        assert!((psms[1].modifications[0].mass_delta - 15.994915).abs() < 0.001);
        assert_eq!(psms[1].modifications[0].residues, vec!['M']);

        // Third PSM: Carbamidomethyl at position 2 (C)
        assert_eq!(psms[2].sequence, "ACDEFGHIK");
        assert_eq!(psms[2].modifications.len(), 1);
        assert_eq!(psms[2].modifications[0].name, "Carbamidomethyl");
        assert_eq!(psms[2].modifications[0].residues, vec!['C']);
        assert_eq!(psms[2].raw_name, "hela_Rep2");
    }

    #[test]
    fn parse_custom_json_rt_conversion() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.json");
        std::fs::write(&path, r#"[{"sequence":"PEPTIDE","charge":2,"modify":[],"rt":1.0,"precursor_mz":400.0,"raw_title":"test","protein_names":[]}]"#).unwrap();

        let db = UnimodDb::builtin();
        let psms = CustomJsonParser.parse(&path, &db).unwrap();
        // 1 minute = 60 seconds
        assert!((psms[0].rt_sec - 60.0).abs() < 0.01);
    }

    #[test]
    fn parse_custom_json_file_not_found() {
        let db = UnimodDb::builtin();
        let result = CustomJsonParser.parse(Path::new("/nonexistent.json"), &db);
        assert!(result.is_err());
    }
}
```

- [ ] **Step 2: Add `tempfile` dev-dependency to `Cargo.toml`**

Add to `crates/result-import/Cargo.toml` under `[dev-dependencies]`:

```toml
tempfile = "3"
```

- [ ] **Step 3: Run tests**

Run: `cargo test -p protein-copilot-result-import -- custom_json`

Expected: 3 tests pass.

- [ ] **Step 4: Commit**

```bash
git add -A && git commit -m "feat(result-import): add CustomJsonParser for hela.json format

Parses JSON array of PSMs, converts modify [[pos, unimod_id]] via UnimodDb,
RT minutes → seconds. Graceful handling of modification errors."
```

---

### Task 4: ScanMatcher — RT + Isolation Window Matching

**Files:**
- Create: `crates/result-import/src/scan_matcher.rs`

**Key algorithm:** Pre-scan mzML for all MS2 `(scan, rt_sec, isolation_window)`, sort by RT, binary search for candidates within `rt_tolerance_sec`, filter by `precursor_mz ∈ isolation_window`, pick closest RT.

- [ ] **Step 1: Implement ScanMatcher**

Create `crates/result-import/src/scan_matcher.rs`:

```rust
//! Scan matching: associates imported PSMs with mzML MS2 scans.
//!
//! Algorithm:
//! 1. Pre-scan mzML to collect (scan_number, rt_sec, isolation_window) for all MS2 spectra
//! 2. Sort by RT
//! 3. For each PSM: binary search for RT-proximate candidates, filter by isolation window,
//!    pick the closest RT match

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use protein_copilot_core::spectrum::{MsLevel, Spectrum};
use protein_copilot_spectrum_io::reader::SpectrumReader;

use crate::{FileMatchStats, ImportedPsm, MatchReport, ResultImportError};

/// MS2 spectrum info extracted from mzML for scan matching.
#[derive(Debug, Clone)]
struct Ms2Info {
    scan_number: u32,
    rt_sec: f64,
    /// (target_mz, lower_offset, upper_offset)
    isolation_window: Option<(f64, f64, f64)>,
}

/// Scan matcher configuration.
pub struct ScanMatcherConfig {
    pub rt_tolerance_sec: f64,
    pub mzml_dir: PathBuf,
}

/// Match imported PSMs to mzML MS2 scans.
///
/// PSMs are matched by:
/// 1. `raw_name` → mzML file (raw_name + ".mzML" in `mzml_dir`)
/// 2. RT proximity (within `rt_tolerance_sec`)
/// 3. precursor_mz falls within the MS2's isolation window
///
/// Returns the mutated PSMs (with `matched_scan` and `rt_delta_sec` filled)
/// and a `MatchReport` with quality statistics.
pub fn match_scans(
    psms: &mut [ImportedPsm],
    config: &ScanMatcherConfig,
    reader_factory: &dyn Fn(&Path) -> Result<Arc<dyn SpectrumReader>, ResultImportError>,
) -> Result<MatchReport, ResultImportError> {
    // Group PSMs by raw_name
    let mut groups: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, psm) in psms.iter().enumerate() {
        groups.entry(psm.raw_name.clone()).or_default().push(i);
    }

    let mut per_file = HashMap::new();
    let mut all_rt_deltas: Vec<f64> = Vec::new();
    let mut total_matched = 0usize;
    let mut total_unmatched = 0usize;

    for (raw_name, indices) in &groups {
        let mzml_path = config.mzml_dir.join(format!("{raw_name}.mzML"));
        if !mzml_path.exists() {
            // Try lowercase extension
            let mzml_path_lower = config.mzml_dir.join(format!("{raw_name}.mzml"));
            if !mzml_path_lower.exists() {
                let available = list_mzml_files(&config.mzml_dir);
                return Err(ResultImportError::MzmlNotFound {
                    raw_name: raw_name.clone(),
                    dir: config.mzml_dir.clone(),
                    available,
                });
            }
        }
        let actual_path = if mzml_path.exists() {
            mzml_path
        } else {
            config.mzml_dir.join(format!("{raw_name}.mzml"))
        };

        let reader = reader_factory(&actual_path)?;

        // Pre-scan all MS2 spectra
        let ms2_infos = collect_ms2_info(&*reader)?;
        let ms2_count = ms2_infos.len();

        // Sort by RT for binary search
        let mut sorted_ms2 = ms2_infos;
        sorted_ms2.sort_by(|a, b| a.rt_sec.partial_cmp(&b.rt_sec).unwrap_or(std::cmp::Ordering::Equal));

        let mut file_matched = 0usize;
        let mut file_unmatched = 0usize;

        for &idx in indices {
            let psm = &psms[idx];
            if let Some((scan, delta)) =
                find_best_match(&sorted_ms2, psm.rt_sec, psm.precursor_mz, config.rt_tolerance_sec)
            {
                let psm_mut = &mut psms[idx];
                psm_mut.matched_scan = Some(scan);
                psm_mut.rt_delta_sec = Some(delta);
                all_rt_deltas.push(delta.abs());
                file_matched += 1;
            } else {
                file_unmatched += 1;
            }
        }

        total_matched += file_matched;
        total_unmatched += file_unmatched;

        per_file.insert(
            raw_name.clone(),
            FileMatchStats {
                total: indices.len(),
                matched: file_matched,
                ms2_count,
            },
        );
    }

    all_rt_deltas.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median_rt_delta = if all_rt_deltas.is_empty() {
        0.0
    } else {
        all_rt_deltas[all_rt_deltas.len() / 2]
    };
    let max_rt_delta = all_rt_deltas.last().copied().unwrap_or(0.0);

    Ok(MatchReport {
        total_psms: psms.len(),
        matched: total_matched,
        unmatched: total_unmatched,
        median_rt_delta_sec: median_rt_delta,
        max_rt_delta_sec: max_rt_delta,
        per_file,
    })
}

/// Collect (scan, rt, isolation_window) for all MS2 spectra from a reader.
fn collect_ms2_info(
    reader: &dyn SpectrumReader,
) -> Result<Vec<Ms2Info>, ResultImportError> {
    let summary = reader.read_all_spectra().map_err(|e| ResultImportError::SpectrumIo(e.to_string()))?;
    let mut infos = Vec::new();
    for spec in &summary {
        if spec.ms_level == MsLevel::MS2 {
            let isolation = spec.precursors.first().and_then(|p| {
                p.isolation_window.as_ref().map(|w| {
                    (w.target_mz, w.lower_offset, w.upper_offset)
                })
            });
            infos.push(Ms2Info {
                scan_number: spec.scan_number,
                rt_sec: spec.retention_time_sec,
                isolation_window: isolation,
            });
        }
    }
    Ok(infos)
}

/// Find the best matching MS2 for a given PSM.
///
/// Returns `(scan_number, rt_delta_sec)` or `None` if no match found.
fn find_best_match(
    sorted_ms2: &[Ms2Info],
    psm_rt_sec: f64,
    psm_precursor_mz: f64,
    rt_tolerance_sec: f64,
) -> Option<(u32, f64)> {
    if sorted_ms2.is_empty() {
        return None;
    }

    // Binary search for the closest RT
    let insert_pos = sorted_ms2
        .partition_point(|m| m.rt_sec < psm_rt_sec - rt_tolerance_sec);

    let mut best: Option<(u32, f64)> = None;

    for ms2 in sorted_ms2[insert_pos..].iter() {
        let delta = ms2.rt_sec - psm_rt_sec;
        if delta > rt_tolerance_sec {
            break; // past the RT window
        }
        if delta.abs() > rt_tolerance_sec {
            continue;
        }

        // Check isolation window
        if let Some((target, lower, upper)) = ms2.isolation_window {
            let low = target - lower;
            let high = target + upper;
            if psm_precursor_mz < low || psm_precursor_mz > high {
                continue; // precursor_mz outside isolation window
            }
        }
        // else: no isolation window info → accept based on RT only (DDA fallback)

        match &best {
            None => best = Some((ms2.scan_number, delta)),
            Some((_, best_delta)) => {
                if delta.abs() < best_delta.abs() {
                    best = Some((ms2.scan_number, delta));
                }
            }
        }
    }

    best
}

/// List .mzML files in a directory for error messages.
fn list_mzml_files(dir: &Path) -> String {
    match std::fs::read_dir(dir) {
        Ok(entries) => {
            let files: Vec<String> = entries
                .filter_map(|e| e.ok())
                .filter(|e| {
                    e.path()
                        .extension()
                        .is_some_and(|ext| ext.eq_ignore_ascii_case("mzml"))
                })
                .map(|e| e.file_name().to_string_lossy().to_string())
                .collect();
            if files.is_empty() {
                "none".to_string()
            } else {
                files.join(", ")
            }
        }
        Err(_) => "directory not readable".to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_ms2_infos() -> Vec<Ms2Info> {
        vec![
            Ms2Info { scan_number: 10, rt_sec: 100.0, isolation_window: Some((500.0, 12.5, 12.5)) },
            Ms2Info { scan_number: 20, rt_sec: 200.0, isolation_window: Some((600.0, 12.5, 12.5)) },
            Ms2Info { scan_number: 30, rt_sec: 300.0, isolation_window: Some((500.0, 12.5, 12.5)) },
            Ms2Info { scan_number: 40, rt_sec: 400.0, isolation_window: Some((700.0, 12.5, 12.5)) },
            Ms2Info { scan_number: 50, rt_sec: 500.0, isolation_window: Some((500.0, 25.0, 25.0)) }, // wide DIA window
        ]
    }

    #[test]
    fn find_best_match_exact_rt_and_mz() {
        let ms2s = make_ms2_infos();
        let result = find_best_match(&ms2s, 100.0, 500.0, 30.0);
        assert_eq!(result, Some((10, 0.0)));
    }

    #[test]
    fn find_best_match_within_tolerance() {
        let ms2s = make_ms2_infos();
        // RT=198, mz=605 → should match scan 20 (RT=200, window 587.5–612.5)
        let result = find_best_match(&ms2s, 198.0, 605.0, 30.0);
        assert_eq!(result.unwrap().0, 20);
        assert!((result.unwrap().1 - 2.0).abs() < 0.01);
    }

    #[test]
    fn find_best_match_mz_outside_window() {
        let ms2s = make_ms2_infos();
        // RT=100, mz=550 → scan 10 has window 487.5–512.5, mz=550 is outside
        let result = find_best_match(&ms2s, 100.0, 550.0, 30.0);
        assert!(result.is_none());
    }

    #[test]
    fn find_best_match_rt_outside_tolerance() {
        let ms2s = make_ms2_infos();
        // RT=150, tolerance=30 → nearest scan 10 (RT=100) is 50s away
        let result = find_best_match(&ms2s, 150.0, 500.0, 30.0);
        assert!(result.is_none());
    }

    #[test]
    fn find_best_match_wide_dia_window() {
        let ms2s = make_ms2_infos();
        // scan 50: RT=500, window 475–525 (wide DIA)
        let result = find_best_match(&ms2s, 502.0, 520.0, 30.0);
        assert_eq!(result.unwrap().0, 50);
    }

    #[test]
    fn find_best_match_picks_closest_rt() {
        // Two scans at similar RT, both covering the mz
        let ms2s = vec![
            Ms2Info { scan_number: 1, rt_sec: 100.0, isolation_window: Some((500.0, 25.0, 25.0)) },
            Ms2Info { scan_number: 2, rt_sec: 105.0, isolation_window: Some((500.0, 25.0, 25.0)) },
        ];
        // PSM at RT=103 → closer to scan 2 (105)
        let result = find_best_match(&ms2s, 103.0, 500.0, 30.0);
        assert_eq!(result.unwrap().0, 2);
    }

    #[test]
    fn find_best_match_no_isolation_window_fallback() {
        // DDA: no isolation window → match on RT only
        let ms2s = vec![
            Ms2Info { scan_number: 1, rt_sec: 100.0, isolation_window: None },
        ];
        let result = find_best_match(&ms2s, 105.0, 999.0, 30.0);
        assert_eq!(result.unwrap().0, 1);
    }

    #[test]
    fn find_best_match_empty_ms2_list() {
        let result = find_best_match(&[], 100.0, 500.0, 30.0);
        assert!(result.is_none());
    }
}
```

- [ ] **Step 2: Run tests**

Run: `cargo test -p protein-copilot-result-import -- scan_matcher`

Expected: 8 tests pass.

- [ ] **Step 3: Commit**

```bash
git add -A && git commit -m "feat(result-import): add ScanMatcher with RT + isolation window matching

Binary search for RT-proximate MS2 candidates, filter by isolation window,
pick closest RT. DIA-aware: wide windows allow multiple PSM matches.
8 unit tests covering exact match, tolerance, DIA, edge cases."
```

---

### Task 5: Converter — ImportedPsm → SearchResult

**Files:**
- Create: `crates/result-import/src/converter.rs`

**Key:** Convert matched `ImportedPsm`s to `core::Psm`, aggregate to `PeptideResult`/`ProteinResult`, build `SearchResult` with proper metadata. Uses `search-engine::chemistry::peptide_mass` and `peptide_mz` for `calculated_mz`.

- [ ] **Step 1: Add `protein-copilot-search-engine` dependency**

Add to `crates/result-import/Cargo.toml` `[dependencies]`:

```toml
protein-copilot-search-engine = { workspace = true }
```

- [ ] **Step 2: Implement converter**

Create `crates/result-import/src/converter.rs`:

```rust
//! Converts matched ImportedPsms to core::SearchResult.

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use protein_copilot_core::engine::EngineInfo;
use protein_copilot_core::run_metadata::RunMetadata;
use protein_copilot_core::search_params::{
    DecoyStrategy, Enzyme, MassTolerance, SearchParams, ToleranceUnit,
};
use protein_copilot_core::search_result::{
    PeptideResult, ProteinResult, SearchResult, SearchResultSummary,
};
use protein_copilot_core::search_result::Psm;
use protein_copilot_search_engine::chemistry::{peptide_mass, peptide_mz};
use uuid::Uuid;

use crate::{ImportResult, ImportedPsm, MatchReport};

/// Convert matched ImportedPsms into a standard SearchResult.
///
/// Only PSMs with `matched_scan.is_some()` are included.
/// Returns `(SearchResult, ImportResult)`.
pub fn build_search_result(
    psms: &[ImportedPsm],
    match_report: MatchReport,
    format_name: &str,
    input_files: Vec<PathBuf>,
) -> (SearchResult, ImportResult) {
    let run_id = Uuid::new_v4();

    // Convert matched PSMs to core::Psm
    let core_psms: Vec<Psm> = psms
        .iter()
        .filter(|p| p.matched_scan.is_some())
        .map(|p| to_core_psm(p))
        .collect();

    // Aggregate peptides
    let peptides = aggregate_peptides(&core_psms);

    // Aggregate proteins
    let proteins = aggregate_proteins(&core_psms);

    // Build summary
    let unique_peptides: HashSet<&str> = core_psms.iter().map(|p| p.peptide_sequence.as_str()).collect();
    let unique_proteins: HashSet<&str> = core_psms
        .iter()
        .flat_map(|p| p.protein_accessions.iter().map(|a| a.as_str()))
        .collect();

    let scores: Vec<f64> = core_psms.iter().map(|p| p.score).collect();
    let delta_ppms: Vec<f64> = core_psms.iter().map(|p| p.delta_mass_ppm).collect();
    let median_score = median(&scores);
    let median_delta = median(&delta_ppms);

    let mut charge_dist: HashMap<i32, u64> = HashMap::new();
    let mut mod_dist: HashMap<String, u64> = HashMap::new();
    for psm in &core_psms {
        *charge_dist.entry(psm.charge).or_default() += 1;
        for m in &psm.modifications {
            *mod_dist.entry(m.name.clone()).or_default() += 1;
        }
    }

    // PSMs at 1% FDR: if we have q-values, filter; otherwise count all
    let has_qvalues = core_psms.iter().any(|p| p.q_value.is_some());
    let psms_at_fdr = if has_qvalues {
        core_psms.iter().filter(|p| p.q_value.is_some_and(|q| q <= 0.01)).count() as u64
    } else {
        core_psms.len() as u64
    };
    let peptides_at_fdr = if has_qvalues {
        core_psms
            .iter()
            .filter(|p| p.q_value.is_some_and(|q| q <= 0.01))
            .map(|p| p.peptide_sequence.as_str())
            .collect::<HashSet<_>>()
            .len() as u64
    } else {
        unique_peptides.len() as u64
    };

    let summary = SearchResultSummary {
        total_spectra_searched: core_psms.len() as u64,
        total_psms: core_psms.len() as u64,
        psms_at_1pct_fdr: psms_at_fdr,
        unique_peptides_at_1pct_fdr: peptides_at_fdr,
        protein_groups_at_1pct_fdr: unique_proteins.len() as u64,
        median_score,
        median_delta_mass_ppm: median_delta,
        identification_rate: 1.0, // all imported PSMs are "identified"
        modification_distribution: mod_dist,
        charge_distribution: charge_dist,
        search_duration_sec: 0.0,
    };

    let engine_info = EngineInfo {
        name: "imported".to_string(),
        version: format_name.to_string(),
        supported_features: vec![],
    };

    // Build a placeholder SearchParams (external results don't have search params)
    let params = SearchParams {
        enzyme: Enzyme::Trypsin,
        missed_cleavages: 2,
        fixed_modifications: vec![],
        variable_modifications: vec![],
        precursor_tolerance: MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        },
        fragment_tolerance: MassTolerance {
            value: 0.02,
            unit: ToleranceUnit::Da,
        },
        database_path: "imported".to_string(),
        decoy_strategy: DecoyStrategy::None,
        acquisition_mode: None,
        max_variable_modifications: 3,
        min_peptide_length: 7,
        max_peptide_length: 50,
    };

    let mut metadata = RunMetadata::new(params.clone(), engine_info.clone(), input_files);
    // Override the auto-generated run_id to match our chosen one
    metadata.run_id = run_id;

    let raw_files: Vec<String> = psms
        .iter()
        .map(|p| p.raw_name.clone())
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();

    let search_result = SearchResult {
        run_id,
        engine_info,
        params_used: params,
        psms: core_psms,
        peptides,
        proteins,
        summary,
        metadata,
    };

    let import_result = ImportResult {
        run_id: run_id.to_string(),
        match_report,
        imported_psm_count: search_result.psms.len(),
        unique_peptides: unique_peptides.len(),
        protein_count: unique_proteins.len(),
        raw_files,
    };

    (search_result, import_result)
}

fn to_core_psm(imported: &ImportedPsm) -> Psm {
    let scan = imported.matched_scan.unwrap_or(0);

    // Calculate theoretical m/z
    let mod_mass: f64 = imported.modifications.iter().map(|m| m.mass_delta).sum();
    let calculated_mz = peptide_mass(&imported.sequence)
        .map(|neutral| peptide_mz(neutral + mod_mass, imported.charge))
        .unwrap_or(imported.precursor_mz); // fallback if non-standard AA

    let delta_ppm = if calculated_mz > 0.0 {
        (imported.precursor_mz - calculated_mz) / calculated_mz * 1e6
    } else {
        0.0
    };

    Psm {
        spectrum_scan: scan,
        peptide_sequence: imported.sequence.clone(),
        modifications: imported.modifications.clone(),
        charge: imported.charge,
        precursor_mz: imported.precursor_mz,
        calculated_mz,
        delta_mass_ppm: delta_ppm,
        score: imported.score.unwrap_or(0.0),
        q_value: imported.q_value,
        protein_accessions: imported.protein_accessions.clone(),
        is_decoy: false,
    }
}

fn aggregate_peptides(psms: &[Psm]) -> Vec<PeptideResult> {
    let mut map: HashMap<&str, (Vec<&str>, f64, u64, Option<f64>)> = HashMap::new();
    for psm in psms {
        let entry = map.entry(psm.peptide_sequence.as_str()).or_insert_with(|| {
            (Vec::new(), f64::MIN, 0, None)
        });
        for acc in &psm.protein_accessions {
            if !entry.0.contains(&acc.as_str()) {
                entry.0.push(acc.as_str());
            }
        }
        if psm.score > entry.1 {
            entry.1 = psm.score;
        }
        entry.2 += 1;
        if let Some(q) = psm.q_value {
            entry.3 = Some(entry.3.map_or(q, |prev: f64| prev.min(q)));
        }
    }

    map.into_iter()
        .map(|(seq, (proteins, best_score, count, q))| PeptideResult {
            sequence: seq.to_string(),
            protein_accessions: proteins.into_iter().map(|s| s.to_string()).collect(),
            best_score: if best_score == f64::MIN { 0.0 } else { best_score },
            q_value: q,
            psm_count: count,
        })
        .collect()
}

fn aggregate_proteins(psms: &[Psm]) -> Vec<ProteinResult> {
    let mut map: HashMap<&str, HashSet<&str>> = HashMap::new();
    for psm in psms {
        for acc in &psm.protein_accessions {
            map.entry(acc.as_str())
                .or_default()
                .insert(psm.peptide_sequence.as_str());
        }
    }

    map.into_iter()
        .map(|(acc, peptides)| ProteinResult {
            accession: acc.to_string(),
            description: String::new(),
            coverage: 0.0, // cannot compute without protein sequence
            peptide_count: peptides.len() as u64,
            unique_peptide_count: peptides.len() as u64,
        })
        .collect()
}

fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    sorted[sorted.len() / 2]
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::search_params::ModPosition;
    use protein_copilot_core::search_params::Modification;

    fn sample_psms() -> Vec<ImportedPsm> {
        vec![
            ImportedPsm {
                sequence: "PEPTIDE".to_string(),
                charge: 2,
                precursor_mz: 400.19,
                rt_sec: 600.0,
                modifications: vec![],
                score: Some(0.001),
                q_value: Some(0.001),
                protein_accessions: vec!["P12345".to_string()],
                raw_name: "test".to_string(),
                matched_scan: Some(100),
                rt_delta_sec: Some(1.5),
            },
            ImportedPsm {
                sequence: "PEPTMIDE".to_string(),
                charge: 3,
                precursor_mz: 310.5,
                rt_sec: 1200.0,
                modifications: vec![Modification {
                    name: "Oxidation".to_string(),
                    mass_delta: 15.994915,
                    residues: vec!['M'],
                    position: ModPosition::Anywhere,
                }],
                score: Some(0.005),
                q_value: Some(0.005),
                protein_accessions: vec!["P12345".to_string()],
                raw_name: "test".to_string(),
                matched_scan: Some(200),
                rt_delta_sec: Some(0.5),
            },
            ImportedPsm {
                sequence: "PEPTIDE".to_string(),
                charge: 2,
                precursor_mz: 400.19,
                rt_sec: 1800.0,
                modifications: vec![],
                score: None,
                q_value: None,
                protein_accessions: vec!["P67890".to_string()],
                raw_name: "test".to_string(),
                matched_scan: None, // unmatched — should be excluded
                rt_delta_sec: None,
            },
        ]
    }

    #[test]
    fn build_search_result_excludes_unmatched() {
        let psms = sample_psms();
        let report = MatchReport {
            total_psms: 3,
            matched: 2,
            unmatched: 1,
            median_rt_delta_sec: 1.0,
            max_rt_delta_sec: 1.5,
            per_file: HashMap::new(),
        };
        let (result, import) = build_search_result(&psms, report, "custom_json", vec![]);
        assert_eq!(result.psms.len(), 2); // only matched
        assert_eq!(import.imported_psm_count, 2);
    }

    #[test]
    fn build_search_result_calculates_mz() {
        let psms = sample_psms();
        let report = MatchReport {
            total_psms: 3, matched: 2, unmatched: 1,
            median_rt_delta_sec: 0.0, max_rt_delta_sec: 0.0,
            per_file: HashMap::new(),
        };
        let (result, _) = build_search_result(&psms, report, "custom_json", vec![]);
        for psm in &result.psms {
            assert!(psm.calculated_mz > 0.0, "calculated_mz should be positive");
            assert!(psm.delta_mass_ppm.is_finite(), "delta_mass_ppm should be finite");
        }
    }

    #[test]
    fn build_search_result_aggregates_peptides() {
        let psms = sample_psms();
        let report = MatchReport {
            total_psms: 3, matched: 2, unmatched: 1,
            median_rt_delta_sec: 0.0, max_rt_delta_sec: 0.0,
            per_file: HashMap::new(),
        };
        let (result, _) = build_search_result(&psms, report, "custom_json", vec![]);
        assert_eq!(result.peptides.len(), 2); // PEPTIDE and PEPTMIDE
    }

    #[test]
    fn build_search_result_engine_info() {
        let psms = sample_psms();
        let report = MatchReport {
            total_psms: 3, matched: 2, unmatched: 1,
            median_rt_delta_sec: 0.0, max_rt_delta_sec: 0.0,
            per_file: HashMap::new(),
        };
        let (result, _) = build_search_result(&psms, report, "diann_parquet", vec![]);
        assert_eq!(result.engine_info.name, "imported");
        assert_eq!(result.engine_info.version, "diann_parquet");
    }
}
```

- [ ] **Step 3: Run tests**

Run: `cargo test -p protein-copilot-result-import -- converter`

Expected: 4 tests pass.

- [ ] **Step 4: Commit**

```bash
git add -A && git commit -m "feat(result-import): add converter ImportedPsm → SearchResult

Converts matched PSMs to core::Psm with calculated_mz from peptide_mass().
Aggregates peptides and proteins. Builds full SearchResult for run_cache."
```

---

### Task 6: DIA-NN Parquet Parser

**Files:**
- Create: `crates/result-import/src/diann.rs`

**Reference:** DIA-NN report.parquet columns — `Modified.Sequence` (e.g. `_AAAC(UniMod:4)DM(UniMod:35)K_`), `Precursor.Charge`, `Precursor.Mz`, `RT` (minutes), `Q.Value`, `Run`, `Protein.Names`.

- [ ] **Step 1: Implement DiannParser**

Create `crates/result-import/src/diann.rs`:

```rust
//! Parser for DIA-NN report.parquet files.
//!
//! Key columns:
//! - `Modified.Sequence`: `_AAAC(UniMod:4)DM(UniMod:35)K_`
//! - `Precursor.Charge`: charge state
//! - `Precursor.Mz`: precursor m/z
//! - `RT`: retention time in minutes
//! - `Q.Value`: FDR q-value
//! - `Run`: raw file name
//! - `Protein.Names`: protein accessions (semicolon-separated)

use std::path::Path;
use std::sync::Arc;

use arrow::array::{AsArray, Float32Array, Float64Array, Int32Array, StringArray};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use regex::Regex;

use protein_copilot_core::search_params::{ModPosition, Modification};

use crate::unimod::UnimodDb;
use crate::{ImportedPsm, ResultImportError, ResultParser};

pub struct DiannParser {
    /// Maximum Q.Value to include. PSMs above this threshold are filtered out.
    pub filter_qvalue: Option<f64>,
    /// Optional: only import PSMs from this specific Run.
    pub run_filter: Option<String>,
}

impl DiannParser {
    pub fn new() -> Self {
        Self {
            filter_qvalue: Some(0.01),
            run_filter: None,
        }
    }
}

/// Parse DIA-NN `Modified.Sequence` like `_AAAC(UniMod:4)DM(UniMod:35)K_`
///
/// Returns `(clean_sequence, Vec<(position_1based, unimod_id)>)`.
fn parse_modified_sequence(modified_seq: &str) -> (String, Vec<(usize, u32)>) {
    let re = Regex::new(r"\(UniMod:(\d+)\)").unwrap();
    let mut clean = String::new();
    let mut mods = Vec::new();

    // Strip leading/trailing underscores
    let trimmed = modified_seq.trim_matches('_');

    let mut pos = 0usize; // 1-based position in clean sequence
    let mut chars = trimmed.char_indices().peekable();
    let bytes = trimmed.as_bytes();

    let mut i = 0;
    while i < trimmed.len() {
        if bytes[i] == b'(' {
            // Find closing paren
            if let Some(j) = trimmed[i..].find(')') {
                let inside = &trimmed[i..i + j + 1];
                if let Some(cap) = re.captures(inside) {
                    if let Ok(id) = cap[1].parse::<u32>() {
                        // This mod applies to the previous residue
                        mods.push((pos, id));
                    }
                }
                i += j + 1;
                continue;
            }
        }

        let ch = trimmed.as_bytes()[i] as char;
        if ch.is_ascii_uppercase() {
            pos += 1;
            clean.push(ch);
        }
        i += 1;
    }

    (clean, mods)
}

impl ResultParser for DiannParser {
    fn parse(
        &self,
        path: &Path,
        unimod: &UnimodDb,
    ) -> Result<Vec<ImportedPsm>, ResultImportError> {
        if !path.exists() {
            return Err(ResultImportError::FileNotFound {
                path: path.to_path_buf(),
            });
        }

        let file = std::fs::File::open(path)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
        let reader = builder.build()?;

        let mut psms = Vec::new();
        let mut filtered_count = 0u64;

        for batch_result in reader {
            let batch = batch_result?;
            let schema = batch.schema();

            // Get column arrays
            let mod_seq_col = get_string_column(&batch, &schema, "Modified.Sequence")?;
            let charge_col = get_int_column(&batch, &schema, "Precursor.Charge")?;
            let mz_col = get_float_column(&batch, &schema, "Precursor.Mz")?;
            let rt_col = get_float_column(&batch, &schema, "RT")?;
            let qvalue_col = get_float_column(&batch, &schema, "Q.Value")?;
            let run_col = get_string_column(&batch, &schema, "Run")?;
            let protein_col = get_string_column_optional(&batch, &schema, "Protein.Names");

            for row in 0..batch.num_rows() {
                // Q.Value filter
                let qvalue = get_f64(&qvalue_col, row);
                if let Some(max_q) = self.filter_qvalue {
                    if qvalue > max_q {
                        filtered_count += 1;
                        continue;
                    }
                }

                let run = get_str(&run_col, row).to_string();

                // Run filter
                if let Some(ref filter) = self.run_filter {
                    if &run != filter {
                        continue;
                    }
                }

                let mod_seq = get_str(&mod_seq_col, row);
                let (sequence, mod_positions) = parse_modified_sequence(mod_seq);

                let mut modifications = Vec::new();
                for (pos, id) in &mod_positions {
                    match unimod.to_modification(*id, *pos, &sequence) {
                        Ok(m) => modifications.push(m),
                        Err(e) => {
                            tracing::debug!("DIA-NN mod conversion: {e}");
                        }
                    }
                }

                let proteins = protein_col
                    .as_ref()
                    .map(|col| {
                        get_str(col, row)
                            .split(';')
                            .map(|s| s.trim().to_string())
                            .filter(|s| !s.is_empty())
                            .collect()
                    })
                    .unwrap_or_default();

                psms.push(ImportedPsm {
                    sequence,
                    charge: get_i32(&charge_col, row),
                    precursor_mz: get_f64(&mz_col, row),
                    rt_sec: get_f64(&rt_col, row) * 60.0, // minutes → seconds
                    modifications,
                    score: Some(qvalue),
                    q_value: Some(qvalue),
                    protein_accessions: proteins,
                    raw_name: run,
                    matched_scan: None,
                    rt_delta_sec: None,
                });
            }
        }

        if filtered_count > 0 {
            tracing::info!("filtered {filtered_count} PSMs above Q.Value threshold");
        }
        tracing::info!("parsed {} PSMs from DIA-NN parquet: {}", psms.len(), path.display());
        Ok(psms)
    }
}

// Column helper functions

fn get_string_column(
    batch: &arrow::record_batch::RecordBatch,
    schema: &arrow::datatypes::SchemaRef,
    name: &str,
) -> Result<Arc<StringArray>, ResultImportError> {
    let idx = schema.index_of(name).map_err(|_| ResultImportError::MissingColumn {
        column: name.to_string(),
        expected: "Modified.Sequence, Precursor.Charge, Precursor.Mz, RT, Q.Value, Run".to_string(),
    })?;
    Ok(Arc::new(
        batch
            .column(idx)
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| ResultImportError::Other(format!("column '{name}' is not String type")))?
            .clone(),
    ))
}

fn get_string_column_optional(
    batch: &arrow::record_batch::RecordBatch,
    schema: &arrow::datatypes::SchemaRef,
    name: &str,
) -> Option<Arc<StringArray>> {
    let idx = schema.index_of(name).ok()?;
    batch
        .column(idx)
        .as_any()
        .downcast_ref::<StringArray>()
        .map(|a| Arc::new(a.clone()))
}

fn get_int_column(
    batch: &arrow::record_batch::RecordBatch,
    schema: &arrow::datatypes::SchemaRef,
    name: &str,
) -> Result<Arc<dyn arrow::array::Array>, ResultImportError> {
    let idx = schema.index_of(name).map_err(|_| ResultImportError::MissingColumn {
        column: name.to_string(),
        expected: "Modified.Sequence, Precursor.Charge, Precursor.Mz, RT, Q.Value, Run".to_string(),
    })?;
    Ok(batch.column(idx).clone())
}

fn get_float_column(
    batch: &arrow::record_batch::RecordBatch,
    schema: &arrow::datatypes::SchemaRef,
    name: &str,
) -> Result<Arc<dyn arrow::array::Array>, ResultImportError> {
    let idx = schema.index_of(name).map_err(|_| ResultImportError::MissingColumn {
        column: name.to_string(),
        expected: "Modified.Sequence, Precursor.Charge, Precursor.Mz, RT, Q.Value, Run".to_string(),
    })?;
    Ok(batch.column(idx).clone())
}

fn get_str(col: &Arc<StringArray>, row: usize) -> &str {
    col.value(row)
}

fn get_i32(col: &Arc<dyn arrow::array::Array>, row: usize) -> i32 {
    if let Some(a) = col.as_any().downcast_ref::<Int32Array>() {
        return a.value(row);
    }
    // DIA-NN might use other int types
    if let Some(a) = col.as_any().downcast_ref::<arrow::array::Int64Array>() {
        return a.value(row) as i32;
    }
    if let Some(a) = col.as_any().downcast_ref::<arrow::array::Int16Array>() {
        return a.value(row) as i32;
    }
    0
}

fn get_f64(col: &Arc<dyn arrow::array::Array>, row: usize) -> f64 {
    if let Some(a) = col.as_any().downcast_ref::<Float64Array>() {
        return a.value(row);
    }
    if let Some(a) = col.as_any().downcast_ref::<Float32Array>() {
        return a.value(row) as f64;
    }
    0.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_modified_sequence_no_mods() {
        let (seq, mods) = parse_modified_sequence("_PEPTIDE_");
        assert_eq!(seq, "PEPTIDE");
        assert!(mods.is_empty());
    }

    #[test]
    fn parse_modified_sequence_single_mod() {
        let (seq, mods) = parse_modified_sequence("_PEPTM(UniMod:35)IDE_");
        assert_eq!(seq, "PEPTMIDE");
        assert_eq!(mods, vec![(5, 35)]); // M is at position 5
    }

    #[test]
    fn parse_modified_sequence_multiple_mods() {
        let (seq, mods) = parse_modified_sequence("_AAAC(UniMod:4)DM(UniMod:35)K_");
        assert_eq!(seq, "AAACDMK");
        assert_eq!(mods, vec![(4, 4), (6, 35)]); // C at pos 4, M at pos 6
    }

    #[test]
    fn parse_modified_sequence_nterm_mod() {
        let (seq, mods) = parse_modified_sequence("_(UniMod:1)PEPTIDE_");
        assert_eq!(seq, "PEPTIDE");
        // N-term mod at position 0 (before first AA)
        assert_eq!(mods, vec![(0, 1)]);
    }

    #[test]
    fn parse_modified_sequence_bare() {
        let (seq, mods) = parse_modified_sequence("PEPTIDE");
        assert_eq!(seq, "PEPTIDE");
        assert!(mods.is_empty());
    }

    #[test]
    fn parse_real_diann_parquet() {
        let path = Path::new("/home/verden/pfind/2025-fall/code/report.parquet");
        if !path.exists() {
            eprintln!("skipping parquet test: report.parquet not found");
            return;
        }
        let db = UnimodDb::builtin();
        let mut parser = DiannParser::new();
        parser.filter_qvalue = Some(0.001); // strict filter for speed
        let psms = parser.parse(path, &db).unwrap();
        assert!(!psms.is_empty(), "should parse some PSMs");
        // Verify RT was converted to seconds
        for psm in &psms[..5.min(psms.len())] {
            assert!(psm.rt_sec > 60.0, "RT should be in seconds (>60), got {}", psm.rt_sec);
        }
        tracing::info!("parsed {} PSMs from real DIA-NN parquet", psms.len());
    }
}
```

- [ ] **Step 2: Run tests**

Run: `cargo test -p protein-copilot-result-import -- diann`

Expected: 5 unit tests pass + 1 integration test if parquet file exists.

- [ ] **Step 3: Commit**

```bash
git add -A && git commit -m "feat(result-import): add DiannParser for DIA-NN report.parquet

Parses Modified.Sequence with UniMod regex, RT minutes → seconds,
Q.Value filtering, flexible column types (f32/f64, i32/i64)."
```

---

### Task 7: pFind Parser Skeleton

**Files:**
- Create: `crates/result-import/src/pfind.rs`

- [ ] **Step 1: Create skeleton**

Create `crates/result-import/src/pfind.rs`:

```rust
//! Parser for pFind .spectra result files (skeleton).
//!
//! pFind results include scan numbers directly, so scan matching
//! is not required. Implementation pending sample file availability.

use std::path::Path;

use crate::unimod::UnimodDb;
use crate::{ImportedPsm, ResultImportError, ResultParser};

/// Parser for pFind .spectra result files.
pub struct PFindParser;

impl ResultParser for PFindParser {
    fn parse(
        &self,
        path: &Path,
        _unimod: &UnimodDb,
    ) -> Result<Vec<ImportedPsm>, ResultImportError> {
        if !path.exists() {
            return Err(ResultImportError::FileNotFound {
                path: path.to_path_buf(),
            });
        }

        // TODO: Implement when sample .spectra file is available.
        // pFind results include scan numbers, so scan matching is not needed.
        Err(ResultImportError::Other(
            "pFind .spectra parser not yet implemented — awaiting sample file".to_string(),
        ))
    }
}

/// Detect whether a file is a pFind .spectra file.
pub fn detect(path: &Path) -> bool {
    path.extension().is_some_and(|e| e == "spectra")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detect_spectra_extension() {
        assert!(detect(Path::new("result.spectra")));
        assert!(!detect(Path::new("result.json")));
    }

    #[test]
    fn parse_returns_not_implemented() {
        let db = UnimodDb::builtin();
        let result = PFindParser.parse(Path::new("/nonexistent.spectra"), &db);
        assert!(result.is_err());
    }
}
```

- [ ] **Step 2: Run tests**

Run: `cargo test -p protein-copilot-result-import -- pfind`

Expected: 2 tests pass.

- [ ] **Step 3: Commit**

```bash
git add -A && git commit -m "feat(result-import): add pFind parser skeleton

Skeleton implementation with ResultParser trait. Returns error until
sample .spectra file is available for development."
```

---

### Task 8: MCP Tool Integration — import_search_results

**Files:**
- Modify: `crates/mcp-server/Cargo.toml`
- Modify: `crates/mcp-server/src/tools.rs`
- Modify: `Cargo.toml` (workspace root, if not done in Task 1)

- [ ] **Step 1: Add result-import dependency to mcp-server**

Add to `crates/mcp-server/Cargo.toml` `[dependencies]`:

```toml
protein-copilot-result-import = { workspace = true }
```

- [ ] **Step 2: Add import types and tool to tools.rs**

Add these imports to the top of `crates/mcp-server/src/tools.rs`:

```rust
use protein_copilot_result_import::{
    self, ImportFormat, ImportResult, ImportedPsm, ResultParser,
    custom_json::CustomJsonParser, diann::DiannParser, pfind::PFindParser,
    unimod::UnimodDb, scan_matcher::{ScanMatcherConfig, match_scans},
    converter::build_search_result,
};
```

Add the input struct in the tool input types section:

```rust
/// Input for the import_search_results tool.
#[derive(Debug, Deserialize, schemars::JsonSchema)]
struct ImportSearchResultsInput {
    /// Path to external search result file (.json, .parquet, .spectra).
    result_file: String,
    /// Result file format. 'auto' detects from extension.
    #[serde(default = "default_import_format")]
    format: String,
    /// Directory containing mzML files. File association: raw_name + '.mzML'.
    mzml_dir: String,
    /// Path to unimod.xml. If not provided, uses builtin modification database.
    #[serde(default)]
    unimod_path: Option<String>,
    /// RT tolerance in seconds for scan matching. Default: 30.
    #[serde(default = "default_rt_tolerance")]
    rt_tolerance_sec: f64,
    /// Q-value threshold for filtering (DIA-NN). Default: 0.01.
    #[serde(default = "default_filter_qvalue")]
    filter_qvalue: f64,
    /// Optional: only import PSMs from this specific run/raw_title.
    #[serde(default)]
    run_filter: Option<String>,
}

fn default_import_format() -> String { "auto".to_string() }
fn default_rt_tolerance() -> f64 { 30.0 }
fn default_filter_qvalue() -> f64 { 0.01 }
```

Add the tool method in the `#[rmcp::tool_handler] impl ServerHandler for ProteinCopilotServer` block (at the end, before the closing brace of the tool handler impl):

```rust
    /// Import external search results and match to mzML scans.
    #[rmcp::tool(
        name = "import_search_results",
        description = "Import external search results (DIA-NN, custom JSON, pFind) and match to mzML scans. Returns a run_id for use with annotate_spectrum, extract_xic, and generate_summary."
    )]
    async fn import_search_results(
        &self,
        Parameters(input): Parameters<ImportSearchResultsInput>,
    ) -> Result<Json<ImportResult>, ErrorData> {
        let result_path = PathBuf::from(&input.result_file);
        let mzml_dir = PathBuf::from(&input.mzml_dir);

        if !result_path.exists() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!("result file not found: {}", result_path.display()),
            ));
        }
        if !mzml_dir.is_dir() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                format!("mzml_dir is not a directory: {}", mzml_dir.display()),
            ));
        }

        // Load Unimod database
        let unimod = if let Some(ref xml_path) = input.unimod_path {
            UnimodDb::from_xml(Path::new(xml_path))
                .map_err(|e| mcp_err(ErrorCode::INVALID_PARAMS, format!("unimod.xml error: {e}")))?
        } else {
            UnimodDb::builtin()
        };

        // Detect format
        let format = match input.format.as_str() {
            "auto" => protein_copilot_result_import::detect_format(&result_path)
                .map_err(|e| mcp_err(ErrorCode::INVALID_PARAMS, e.to_string()))?,
            "custom_json" => ImportFormat::CustomJson,
            "diann_parquet" => ImportFormat::DiannParquet,
            "pfind_spectra" => ImportFormat::PFindSpectra,
            other => {
                return Err(mcp_err(
                    ErrorCode::INVALID_PARAMS,
                    format!("unknown format: '{other}'. Supported: auto, custom_json, diann_parquet, pfind_spectra"),
                ));
            }
        };

        // Parse
        let mut psms: Vec<ImportedPsm> = match format {
            ImportFormat::CustomJson => {
                CustomJsonParser.parse(&result_path, &unimod)
                    .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?
            }
            ImportFormat::DiannParquet => {
                let mut parser = DiannParser::new();
                parser.filter_qvalue = Some(input.filter_qvalue);
                parser.run_filter = input.run_filter.clone();
                parser.parse(&result_path, &unimod)
                    .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?
            }
            ImportFormat::PFindSpectra => {
                PFindParser.parse(&result_path, &unimod)
                    .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?
            }
        };

        if psms.is_empty() {
            return Err(mcp_err(
                ErrorCode::INVALID_PARAMS,
                "no PSMs parsed from the result file (check format and filters)",
            ));
        }

        // Scan matching
        let reader_cache = Arc::clone(&self.reader_cache);
        let config = ScanMatcherConfig {
            rt_tolerance_sec: input.rt_tolerance_sec,
            mzml_dir: mzml_dir.clone(),
        };
        let match_report = match_scans(&mut psms, &config, &|path| {
            self.get_or_create_reader(path)
                .map_err(|_| protein_copilot_result_import::ResultImportError::SpectrumIo(
                    format!("failed to open {}", path.display()),
                ))
        }).map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?;

        // Convert to SearchResult
        let format_name = match format {
            ImportFormat::CustomJson => "custom_json",
            ImportFormat::DiannParquet => "diann_parquet",
            ImportFormat::PFindSpectra => "pfind_spectra",
        };
        let input_files = vec![result_path];
        let (search_result, import_result) =
            build_search_result(&psms, match_report, format_name, input_files);

        // Store in run_cache
        {
            let mut cache = self
                .run_cache
                .lock()
                .map_err(|_| mcp_err(ErrorCode::INTERNAL_ERROR, "run cache lock poisoned"))?;
            cache.evict_if_full();
            cache.insert(
                search_result.run_id,
                RunState {
                    progress: SearchProgress {
                        run_id: search_result.run_id,
                        status: "Completed".to_string(),
                        stage: Some("Imported".to_string()),
                        progress_pct: Some(100.0),
                        elapsed_sec: 0.0,
                        estimated_remaining_sec: None,
                    },
                    result: Some(search_result),
                    handle: None,
                },
            );
        }

        Ok(Json(import_result))
    }
```

- [ ] **Step 3: Fix `get_or_create_reader` visibility**

The `get_or_create_reader` method on `ProteinCopilotServer` takes `&self` and returns `Result<Arc<dyn SpectrumReader>, ErrorData>`. In the import tool, we need to pass it as a closure. Verify the existing method signature works or adapt the closure.

Check if `get_or_create_reader` is already a `&self` method (it should be based on prior work). The closure adapts the error type from `ErrorData` to `ResultImportError`.

- [ ] **Step 4: Run `cargo check -p protein-copilot-mcp-server`**

Expected: compiles without errors.

- [ ] **Step 5: Run full workspace tests**

Run: `cargo test --workspace`

Expected: all existing tests pass + new result-import tests pass.

- [ ] **Step 6: Commit**

```bash
git add -A && git commit -m "feat(mcp-server): add import_search_results tool

Wires result-import crate into MCP server. Supports auto format detection,
Unimod builtin/XML, scan matching, and stores result in run_cache for
use with annotate_spectrum, extract_xic, and generate_summary."
```

---

### Task 9: Full Verification

- [ ] **Step 1: Run clippy**

Run: `cargo clippy --workspace -- -D warnings`

Expected: no warnings.

- [ ] **Step 2: Run full test suite**

Run: `cargo test --workspace`

Expected: all tests pass.

- [ ] **Step 3: Verify tool appears in MCP server tool list**

The `import_search_results` tool should appear when the server reports its capabilities via `#[rmcp::tool(...)]`.

- [ ] **Step 4: Final commit if any fixes needed**

```bash
git add -A && git commit -m "fix: address clippy warnings and test issues"
```
