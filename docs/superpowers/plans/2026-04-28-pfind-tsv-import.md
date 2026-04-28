# pFind TSV Import — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Support pFind result TSV files in both `import_search_results` (spectrum annotation) and `classify_entrapment_hits` (entrapment analysis) MCP tools.

**Architecture:** Add a dedicated pFind TSV parser to both `result-import` and `entrapment-analysis` crates. Auto-detect format by header column names (`PeptideSequence`, `ScanNo`, `FileName`). Since pFind TSV already contains scan numbers, the `import_search_results` path skips RT-based scan matching entirely.

**Tech Stack:** Rust, csv crate, existing `ImportedPsm`/`UnifiedPsm` types, MCP Server tools layer.

---

## File Structure

| Action | File | Responsibility |
|--------|------|---------------|
| Create | `crates/result-import/src/pfind_tsv.rs` | Parse pFind TSV → `Vec<ImportedPsm>` with `ResultParser` trait |
| Create | `crates/entrapment-analysis/src/loader/pfind_tsv.rs` | Parse pFind TSV → `Vec<UnifiedPsm>` |
| Create | `tests/fixtures/pfind_sample.tsv` | 10-row fixture from real data |
| Modify | `crates/result-import/Cargo.toml` | Add `csv = "1"` dependency |
| Modify | `crates/result-import/src/lib.rs` | Add `PFindTsv` variant to `ImportFormat`, update `detect_format`, add `pub mod pfind_tsv` |
| Modify | `crates/entrapment-analysis/src/loader/mod.rs` | Add `PFindTsv` variant to `ResultFormat`, update `from_path` and `load_psms` |
| Modify | `crates/mcp-server/src/tools.rs` | Wire `pfind_tsv` format in `import_search_results` and `classify_entrapment_hits` |

## Dependency Graph

```
T1 (fixture) ──────────────────────────┐
T2 (result-import parser) ─────────────┤
T3 (entrapment-analysis parser) ───────┼── T5 (MCP Server) ── T6 (verify)
T4 (auto-detect in both crates) ───────┘
```

- T1 is independent (fixture file)
- T2, T3 are independent of each other (different crates)
- T4 depends on T2 and T3 (format detection references the parsers)
- T5 depends on T2, T3, T4 (MCP server wires everything)
- T6 depends on T5 (end-to-end verification)

---

### Task 1: Create Test Fixture

**Files:**
- Create: `tests/fixtures/pfind_sample.tsv`

- [ ] **Step 1: Create the fixture file**

Extract the header + first 10 data rows from `output/intersection_non_human.tsv`. Include variety: rows with modifications, empty modifications, different charge states, different raw file names.

```bash
# Header + selected rows covering: no mods (rows 2,5,6), single mod (rows 4,14), double mod (find one), different charges (1,2,3,4)
head -1 output/intersection_non_human.tsv > tests/fixtures/pfind_sample.tsv
# Row with no mod, charge 2
sed -n '2p' output/intersection_non_human.tsv >> tests/fixtures/pfind_sample.tsv
# Row with no mod, charge 3
sed -n '3p' output/intersection_non_human.tsv >> tests/fixtures/pfind_sample.tsv
# Row with single Carbamidomethyl, charge 3
sed -n '4p' output/intersection_non_human.tsv >> tests/fixtures/pfind_sample.tsv
# Row with no mod, charge 4
sed -n '7p' output/intersection_non_human.tsv >> tests/fixtures/pfind_sample.tsv
# Row with no mod, charge 2
sed -n '9p' output/intersection_non_human.tsv >> tests/fixtures/pfind_sample.tsv
# Row with single Carbamidomethyl
sed -n '14p' output/intersection_non_human.tsv >> tests/fixtures/pfind_sample.tsv
# Row with Oxidation
grep 'Oxidation' output/intersection_non_human.tsv | head -1 >> tests/fixtures/pfind_sample.tsv
# Row with Acetyl
grep 'Acetyl' output/intersection_non_human.tsv | head -1 >> tests/fixtures/pfind_sample.tsv
# Row with double Carbamidomethyl
grep -P 'Carbamidomethyl.*Carbamidomethyl' output/intersection_non_human.tsv | head -1 >> tests/fixtures/pfind_sample.tsv
# Row with charge 1
awk -F'\t' '$10==1' output/intersection_non_human.tsv | head -1 >> tests/fixtures/pfind_sample.tsv
```

- [ ] **Step 2: Verify the fixture**

```bash
wc -l tests/fixtures/pfind_sample.tsv   # expect 11-12 lines (header + 10-11 data rows)
head -2 tests/fixtures/pfind_sample.tsv  # verify header and first row
```

- [ ] **Step 3: Commit**

```bash
git add tests/fixtures/pfind_sample.tsv
git commit -m "test: add pFind TSV fixture for parser tests"
```

---

### Task 2: Implement `result-import` pFind TSV Parser

**Files:**
- Modify: `crates/result-import/Cargo.toml` (add `csv = "1"`)
- Create: `crates/result-import/src/pfind_tsv.rs`
- Modify: `crates/result-import/src/lib.rs` (add module, update `ImportFormat`, update `detect_format`)

- [ ] **Step 1: Add csv dependency**

In `crates/result-import/Cargo.toml`, add under `[dependencies]`:

```toml
csv = "1"
```

- [ ] **Step 2: Write the parser module `pfind_tsv.rs`**

Create `crates/result-import/src/pfind_tsv.rs` with:

```rust
//! Parser for pFind result TSV files.
//!
//! pFind TSV files contain tab-separated PSM results with columns:
//! FileName, PeptideSequence, Modifications, PepMass, PredRT, CleavageType,
//! ProNCTerm, Proteins, MH+, Charge, ScanNo, RawScore, DeltaMassPPM,
//! DeltaRT(Min), FinalScore, QValue.
//!
//! Key features:
//! - ScanNo is already present → no RT-based scan matching needed
//! - Modifications use pFind format: "pos,Name[Residue];" (semicolon-separated)
//! - Proteins have trailing "/" separator

use std::path::Path;

use protein_copilot_core::search_params::{ModPosition, Modification};

use crate::unimod::UnimodDb;
use crate::{ImportedPsm, ResultImportError, ResultParser};

/// Proton mass in Da (used for MH+ → precursor_mz conversion).
const PROTON_MASS: f64 = 1.007276;

/// Required header columns for pFind TSV format detection.
const PFIND_REQUIRED_HEADERS: [&str; 3] = ["PeptideSequence", "ScanNo", "FileName"];

/// Parser for pFind result TSV files.
pub struct PFindTsvParser;

impl ResultParser for PFindTsvParser {
    fn parse(
        &self,
        path: &Path,
        _unimod: &UnimodDb,
    ) -> Result<Vec<ImportedPsm>, ResultImportError> {
        let _span = tracing::info_span!("parse_pfind_tsv",
            file = %path.display(),
        )
        .entered();

        if !path.exists() {
            return Err(ResultImportError::FileNotFound {
                path: path.to_path_buf(),
            });
        }

        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(path)
            .map_err(|e| ResultImportError::IoError(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("failed to open {}: {e}", path.display()),
            )))?;

        let headers = rdr.headers().map_err(|e| ResultImportError::Other(
            format!("failed to read TSV headers from {}: {e}", path.display()),
        ))?.clone();

        // Find column indices
        let col = |name: &str| -> Result<usize, ResultImportError> {
            headers.iter().position(|h| h == name).ok_or_else(|| {
                ResultImportError::Other(format!(
                    "missing required column '{}' in pFind TSV '{}'",
                    name,
                    path.display()
                ))
            })
        };

        let file_name_idx = col("FileName")?;
        let peptide_idx = col("PeptideSequence")?;
        let mods_idx = col("Modifications")?;
        let proteins_idx = col("Proteins")?;
        let mh_plus_idx = col("MH+")?;
        let charge_idx = col("Charge")?;
        let scan_no_idx = col("ScanNo")?;
        let final_score_idx = col("FinalScore")?;
        let q_value_idx = col("QValue")?;
        let pred_rt_idx = col("PredRT")?;

        let mut psms = Vec::new();
        let mut row_num: usize = 1;

        for result in rdr.records() {
            row_num += 1;
            let record = result.map_err(|e| ResultImportError::Other(
                format!("failed to parse row {} in {}: {e}", row_num, path.display()),
            ))?;

            let sequence = record.get(peptide_idx).unwrap_or("").trim().to_string();
            if sequence.is_empty() {
                continue;
            }

            let raw_name = record.get(file_name_idx).unwrap_or("").trim().to_string();
            let charge: i32 = record.get(charge_idx)
                .unwrap_or("0")
                .trim()
                .parse()
                .unwrap_or(0);
            let mh_plus: f64 = record.get(mh_plus_idx)
                .unwrap_or("0")
                .trim()
                .parse()
                .unwrap_or(0.0);
            let scan_no: u32 = record.get(scan_no_idx)
                .unwrap_or("0")
                .trim()
                .parse()
                .unwrap_or(0);
            let final_score: f64 = record.get(final_score_idx)
                .unwrap_or("0")
                .trim()
                .parse()
                .unwrap_or(0.0);
            let q_value: f64 = record.get(q_value_idx)
                .unwrap_or("1")
                .trim()
                .parse()
                .unwrap_or(1.0);
            let pred_rt: f64 = record.get(pred_rt_idx)
                .unwrap_or("0")
                .trim()
                .parse()
                .unwrap_or(0.0);

            // MH+ → precursor_mz
            let precursor_mz = if charge > 0 {
                (mh_plus + (charge as f64 - 1.0) * PROTON_MASS) / charge as f64
            } else {
                mh_plus
            };

            // Parse modifications
            let mod_str = record.get(mods_idx).unwrap_or("").trim();
            let modifications = parse_pfind_modifications(mod_str);

            // Parse proteins (split by "/" separator, remove trailing empty)
            let protein_str = record.get(proteins_idx).unwrap_or("").trim();
            let protein_accessions = parse_pfind_proteins(protein_str);

            psms.push(ImportedPsm {
                sequence,
                charge,
                precursor_mz,
                rt_min: pred_rt,
                modifications,
                score: Some(final_score),
                q_value: Some(q_value),
                protein_accessions,
                raw_name,
                matched_scan: Some(scan_no), // pFind TSV has scan numbers
                rt_delta_min: None,
            });
        }

        tracing::info!(psm_count = psms.len(), "parsed pFind TSV");
        Ok(psms)
    }
}

/// Parse pFind modification string: "pos,Name[Residue];pos2,Name2[Residue2];"
///
/// Returns a Vec of core `Modification` structs.
/// Position is 1-based in pFind format (converted to match internal use).
///
/// Known modification mass deltas (builtin lookup):
/// - Carbamidomethyl: 57.021464 Da
/// - Oxidation: 15.994915 Da
/// - Acetyl: 42.010565 Da
/// - Phospho: 79.966331 Da
/// - Deamidated: 0.984016 Da
pub fn parse_pfind_modifications(mod_str: &str) -> Vec<Modification> {
    if mod_str.is_empty() {
        return Vec::new();
    }

    let mut mods = Vec::new();
    for part in mod_str.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        // Format: "pos,Name[Residue]"
        let Some((pos_str, name_residue)) = part.split_once(',') else {
            tracing::debug!(entry = part, "skipping malformed modification entry");
            continue;
        };

        let _position: usize = match pos_str.trim().parse() {
            Ok(p) => p,
            Err(_) => {
                tracing::debug!(entry = part, "skipping modification with unparseable position");
                continue;
            }
        };

        // Parse "Name[Residue]" or "Name[ProteinN-term]"
        let name_residue = name_residue.trim();
        let (name, residues, position) = if let Some(bracket_start) = name_residue.find('[') {
            let name = &name_residue[..bracket_start];
            let residue_str = name_residue[bracket_start + 1..]
                .trim_end_matches(']')
                .trim();

            if residue_str == "ProteinN-term" {
                (name.to_string(), vec![], ModPosition::ProteinNTerm)
            } else if residue_str == "AnyN-term" {
                (name.to_string(), vec![], ModPosition::AnyNTerm)
            } else if residue_str == "ProteinC-term" {
                (name.to_string(), vec![], ModPosition::ProteinCTerm)
            } else if residue_str == "AnyC-term" {
                (name.to_string(), vec![], ModPosition::AnyCTerm)
            } else {
                let residues: Vec<char> = residue_str.chars().collect();
                (name.to_string(), residues, ModPosition::Anywhere)
            }
        } else {
            (name_residue.to_string(), vec![], ModPosition::Anywhere)
        };

        let mass_delta = mod_name_to_mass(&name);

        mods.push(Modification {
            name,
            mass_delta,
            residues,
            position,
        });
    }

    mods
}

/// Parse pFind protein string with "/" separator.
///
/// pFind format: "sp|P50475|SYAC_RAT/" (single) or
/// "sp|P50475|SYAC_RAT/sp|Q12345|TEST_HUMAN/" (multiple).
/// Splits by "/" and removes empty entries.
pub fn parse_pfind_proteins(protein_str: &str) -> Vec<String> {
    protein_str
        .split('/')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect()
}

/// Look up mass delta for common modification names.
fn mod_name_to_mass(name: &str) -> f64 {
    match name {
        "Carbamidomethyl" => 57.021464,
        "Oxidation" => 15.994915,
        "Acetyl" => 42.010565,
        "Phospho" => 79.966331,
        "Deamidated" => 0.984016,
        "Methyl" => 14.015650,
        "Dimethyl" => 28.031300,
        "Trimethyl" => 42.046950,
        "GlyGly" => 114.042927,
        "Formyl" => 27.994915,
        "Dehydrated" | "Glu->pyro-Glu" => -18.010565,
        "Gln->pyro-Glu" => -17.026549,
        "Carbamyl" => 43.005814,
        other => {
            tracing::warn!(modification = other, "unknown pFind modification name, using mass_delta=0.0");
            0.0
        }
    }
}

/// Detect whether a TSV file is in pFind format by reading its header line.
///
/// Returns `true` if the header contains all three required columns:
/// `PeptideSequence`, `ScanNo`, and `FileName`.
pub fn detect(path: &Path) -> bool {
    let Ok(mut rdr) = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
    else {
        return false;
    };

    let Ok(headers) = rdr.headers() else {
        return false;
    };

    PFIND_REQUIRED_HEADERS
        .iter()
        .all(|required| headers.iter().any(|h| h == *required))
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::search_params::ModPosition;

    #[test]
    fn parse_empty_modifications() {
        let mods = parse_pfind_modifications("");
        assert!(mods.is_empty());
    }

    #[test]
    fn parse_single_modification() {
        let mods = parse_pfind_modifications("10,Carbamidomethyl[C];");
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].name, "Carbamidomethyl");
        assert!((mods[0].mass_delta - 57.021464).abs() < 1e-6);
        assert_eq!(mods[0].residues, vec!['C']);
        assert_eq!(mods[0].position, ModPosition::Anywhere);
    }

    #[test]
    fn parse_double_modification() {
        let mods = parse_pfind_modifications("1,Carbamidomethyl[C];5,Carbamidomethyl[C];");
        assert_eq!(mods.len(), 2);
        assert_eq!(mods[0].name, "Carbamidomethyl");
        assert_eq!(mods[1].name, "Carbamidomethyl");
    }

    #[test]
    fn parse_oxidation() {
        let mods = parse_pfind_modifications("5,Oxidation[M];");
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].name, "Oxidation");
        assert!((mods[0].mass_delta - 15.994915).abs() < 1e-6);
        assert_eq!(mods[0].residues, vec!['M']);
    }

    #[test]
    fn parse_nterm_modification() {
        let mods = parse_pfind_modifications("0,Acetyl[ProteinN-term];");
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].name, "Acetyl");
        assert!((mods[0].mass_delta - 42.010565).abs() < 1e-6);
        assert!(mods[0].residues.is_empty());
        assert_eq!(mods[0].position, ModPosition::ProteinNTerm);
    }

    #[test]
    fn parse_proteins_single() {
        let proteins = parse_pfind_proteins("sp|P50475|SYAC_RAT/");
        assert_eq!(proteins, vec!["sp|P50475|SYAC_RAT"]);
    }

    #[test]
    fn parse_proteins_multiple() {
        let proteins = parse_pfind_proteins("sp|P50475|SYAC_RAT/sp|Q12345|TEST_HUMAN/");
        assert_eq!(
            proteins,
            vec!["sp|P50475|SYAC_RAT", "sp|Q12345|TEST_HUMAN"]
        );
    }

    #[test]
    fn parse_proteins_empty() {
        let proteins = parse_pfind_proteins("");
        assert!(proteins.is_empty());
    }

    #[test]
    fn precursor_mz_from_mh_plus() {
        // MH+ = 1012.470123, Charge = 2
        let mh_plus = 1012.470123;
        let charge = 2;
        let precursor_mz = (mh_plus + (charge as f64 - 1.0) * PROTON_MASS) / charge as f64;
        assert!((precursor_mz - 506.738700).abs() < 0.001);
    }

    #[test]
    fn parse_fixture_file() {
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/fixtures/pfind_sample.tsv");
        if !fixture.exists() {
            return; // skip if fixture not yet created
        }
        let parser = PFindTsvParser;
        let unimod = crate::unimod::UnimodDb::builtin();
        let psms = parser.parse(&fixture, &unimod).unwrap();
        assert!(!psms.is_empty(), "should parse at least one PSM");
        // All PSMs should have matched_scan set
        for psm in &psms {
            assert!(psm.matched_scan.is_some(), "pFind PSMs should have scan numbers");
        }
    }

    #[test]
    fn detect_pfind_format() {
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/fixtures/pfind_sample.tsv");
        if !fixture.exists() {
            return;
        }
        assert!(detect(&fixture));
    }

    #[test]
    fn detect_non_pfind_tsv() {
        // A TSV with different headers should not be detected as pFind
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("other.tsv");
        std::fs::write(&path, "peptide\tcharge\tprotein_ids\nPEPTIDEK\t2\tP12345\n").unwrap();
        assert!(!detect(&path));
    }
}
```

- [ ] **Step 3: Register the module in `lib.rs`**

In `crates/result-import/src/lib.rs`:

1. Add `pub mod pfind_tsv;` after the existing `pub mod pfind;` line.

2. Add `PFindTsv` variant to the `ImportFormat` enum:

```rust
pub enum ImportFormat {
    CustomJson,
    DiannParquet,
    PFindSpectra,
    PFindTsv,
}
```

3. Update `detect_format()` to handle `.tsv` extension by probing the header:

```rust
pub fn detect_format(path: &Path) -> Result<ImportFormat, ResultImportError> {
    match path.extension().and_then(|e| e.to_str()) {
        Some("json") => Ok(ImportFormat::CustomJson),
        Some("parquet") => Ok(ImportFormat::DiannParquet),
        Some("spectra") => Ok(ImportFormat::PFindSpectra),
        Some("tsv") | Some("txt") => {
            if pfind_tsv::detect(path) {
                Ok(ImportFormat::PFindTsv)
            } else {
                Err(ResultImportError::FormatDetectionFailed {
                    path: path.to_path_buf(),
                })
            }
        }
        _ => Err(ResultImportError::FormatDetectionFailed {
            path: path.to_path_buf(),
        }),
    }
}
```

- [ ] **Step 4: Run tests**

```bash
cargo test -p protein-copilot-result-import 2>&1
```

Expected: all existing tests pass + new `pfind_tsv` tests pass.

- [ ] **Step 5: Commit**

```bash
git add crates/result-import/
git commit -m "feat(result-import): add pFind TSV parser with auto-detection

- PFindTsvParser implements ResultParser trait
- Parse modifications: pos,Name[Residue]; format
- Parse proteins: /-separated with trailing / removal
- MH+ → precursor_mz conversion
- ScanNo directly used as matched_scan (skip RT matching)
- Header-based auto-detection (PeptideSequence+ScanNo+FileName)
- Unit tests for mods, proteins, precursor_mz, fixture parsing"
```

---

### Task 3: Implement `entrapment-analysis` pFind TSV Loader

**Files:**
- Create: `crates/entrapment-analysis/src/loader/pfind_tsv.rs`
- Modify: `crates/entrapment-analysis/src/loader/mod.rs`

- [ ] **Step 1: Create `loader/pfind_tsv.rs`**

```rust
//! Loader for pFind result TSV files into [`UnifiedPsm`] records.
//!
//! pFind TSV files contain: FileName, PeptideSequence, Modifications, PepMass,
//! PredRT, CleavageType, ProNCTerm, Proteins, MH+, Charge, ScanNo, RawScore,
//! DeltaMassPPM, DeltaRT(Min), FinalScore, QValue.

use std::path::Path;

use crate::error::EntrapmentError;
use crate::types::UnifiedPsm;

/// Proton mass in Da.
const PROTON_MASS: f64 = 1.007276;

/// Required header columns for pFind TSV detection.
const PFIND_REQUIRED_HEADERS: [&str; 3] = ["PeptideSequence", "ScanNo", "FileName"];

/// Detect whether a TSV file is in pFind format by checking header columns.
pub fn detect(path: &Path) -> bool {
    let Ok(mut rdr) = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
    else {
        return false;
    };

    let Ok(headers) = rdr.headers() else {
        return false;
    };

    PFIND_REQUIRED_HEADERS
        .iter()
        .all(|required| headers.iter().any(|h| h == *required))
}

/// Load a pFind result TSV file into a vector of [`UnifiedPsm`].
pub fn load_pfind_tsv(path: &Path) -> Result<Vec<UnifiedPsm>, EntrapmentError> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .map_err(|e| EntrapmentError::IoError {
            path: path.to_path_buf(),
            detail: e.to_string(),
        })?;

    let headers = rdr.headers().map_err(|e| EntrapmentError::LoaderError {
        format: "pFind TSV".to_string(),
        detail: format!("failed to read headers from '{}': {e}", path.display()),
    })?.clone();

    let col = |name: &str| -> Option<usize> {
        headers.iter().position(|h| h == name)
    };

    let peptide_idx = col("PeptideSequence").ok_or_else(|| EntrapmentError::LoaderError {
        format: "pFind TSV".to_string(),
        detail: format!("missing 'PeptideSequence' column in '{}'", path.display()),
    })?;

    let charge_idx = col("Charge");
    let mh_plus_idx = col("MH+");
    let pred_rt_idx = col("PredRT");
    let scan_no_idx = col("ScanNo");
    let file_name_idx = col("FileName");
    let proteins_idx = col("Proteins");
    let q_value_idx = col("QValue");
    let mods_idx = col("Modifications");

    let mut psms = Vec::new();

    for result in rdr.records() {
        let record = result.map_err(|e| EntrapmentError::LoaderError {
            format: "pFind TSV".to_string(),
            detail: format!("failed to parse row in '{}': {e}", path.display()),
        })?;

        let peptide = record.get(peptide_idx).unwrap_or("").trim().to_string();
        if peptide.is_empty() {
            continue;
        }

        let charge = charge_idx
            .and_then(|i| record.get(i))
            .and_then(|s| s.trim().parse::<i32>().ok());

        let mh_plus = mh_plus_idx
            .and_then(|i| record.get(i))
            .and_then(|s| s.trim().parse::<f64>().ok());

        let precursor_mz = match (mh_plus, charge) {
            (Some(mh), Some(z)) if z > 0 => {
                Some((mh + (z as f64 - 1.0) * PROTON_MASS) / z as f64)
            }
            (Some(mh), _) => Some(mh),
            _ => None,
        };

        let retention_time = pred_rt_idx
            .and_then(|i| record.get(i))
            .and_then(|s| s.trim().parse::<f64>().ok());

        let scan_number = scan_no_idx
            .and_then(|i| record.get(i))
            .and_then(|s| s.trim().parse::<u32>().ok());

        let spectrum_file = file_name_idx
            .and_then(|i| record.get(i))
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty());

        // Parse protein IDs: split by "/" and rejoin with ";"
        let protein_ids = proteins_idx
            .and_then(|i| record.get(i))
            .map(|s| {
                s.trim()
                    .split('/')
                    .map(|p| p.trim())
                    .filter(|p| !p.is_empty())
                    .collect::<Vec<_>>()
                    .join(";")
            })
            .unwrap_or_default();

        let q_value = q_value_idx
            .and_then(|i| record.get(i))
            .and_then(|s| s.trim().parse::<f64>().ok());

        // Parse modifications: "pos,Name[Residue];" → (0-based position, delta_mass)
        let modifications = mods_idx
            .and_then(|i| record.get(i))
            .map(|s| parse_pfind_mods_to_tuples(s.trim()))
            .unwrap_or_default();

        psms.push(UnifiedPsm {
            peptide,
            charge,
            precursor_mz,
            retention_time,
            rt_start: None,
            rt_stop: None,
            scan_number,
            spectrum_file,
            protein_ids,
            q_value,
            modifications,
        });
    }

    tracing::info!(
        "loaded {} PSMs from pFind TSV: {}",
        psms.len(),
        path.display()
    );
    Ok(psms)
}

/// Parse pFind modification string into (0-based position, delta_mass) tuples.
///
/// Input: "10,Carbamidomethyl[C];5,Oxidation[M];"
/// Output: [(9, 57.021464), (4, 15.994915)]
///
/// Position is converted from 1-based (pFind) to 0-based.
fn parse_pfind_mods_to_tuples(mod_str: &str) -> Vec<(usize, f64)> {
    if mod_str.is_empty() {
        return Vec::new();
    }

    let mut mods = Vec::new();
    for part in mod_str.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }

        let Some((pos_str, name_residue)) = part.split_once(',') else {
            continue;
        };

        let position: usize = match pos_str.trim().parse::<usize>() {
            Ok(p) if p > 0 => p - 1, // 1-based → 0-based
            Ok(p) => p,               // 0 stays 0 (N-term)
            Err(_) => continue,
        };

        // Extract modification name from "Name[Residue]"
        let name = name_residue
            .trim()
            .split('[')
            .next()
            .unwrap_or(name_residue.trim());

        let mass_delta = pfind_mod_mass(name);
        mods.push((position, mass_delta));
    }

    mods
}

/// Look up mass delta for common pFind modification names.
fn pfind_mod_mass(name: &str) -> f64 {
    match name {
        "Carbamidomethyl" => 57.021464,
        "Oxidation" => 15.994915,
        "Acetyl" => 42.010565,
        "Phospho" => 79.966331,
        "Deamidated" => 0.984016,
        "Methyl" => 14.015650,
        "Dimethyl" => 28.031300,
        "Trimethyl" => 42.046950,
        "GlyGly" => 114.042927,
        "Formyl" => 27.994915,
        "Dehydrated" | "Glu->pyro-Glu" => -18.010565,
        "Gln->pyro-Glu" => -17.026549,
        "Carbamyl" => 43.005814,
        _ => {
            tracing::warn!(modification = name, "unknown pFind modification, using delta_mass=0.0");
            0.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_mods_empty() {
        assert!(parse_pfind_mods_to_tuples("").is_empty());
    }

    #[test]
    fn test_parse_mods_single() {
        let mods = parse_pfind_mods_to_tuples("10,Carbamidomethyl[C];");
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].0, 9); // 10 → 9 (0-based)
        assert!((mods[0].1 - 57.021464).abs() < 1e-6);
    }

    #[test]
    fn test_parse_mods_double() {
        let mods = parse_pfind_mods_to_tuples("3,Carbamidomethyl[C];8,Carbamidomethyl[C];");
        assert_eq!(mods.len(), 2);
        assert_eq!(mods[0].0, 2);
        assert_eq!(mods[1].0, 7);
    }

    #[test]
    fn test_parse_mods_nterm() {
        let mods = parse_pfind_mods_to_tuples("0,Acetyl[ProteinN-term];");
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].0, 0);
        assert!((mods[0].1 - 42.010565).abs() < 1e-6);
    }

    #[test]
    fn test_detect_pfind() {
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/fixtures/pfind_sample.tsv");
        if !fixture.exists() {
            return;
        }
        assert!(detect(&fixture));
    }

    #[test]
    fn test_load_fixture() {
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/fixtures/pfind_sample.tsv");
        if !fixture.exists() {
            return;
        }
        let psms = load_pfind_tsv(&fixture).unwrap();
        assert!(!psms.is_empty());
        for psm in &psms {
            assert!(!psm.peptide.is_empty());
            assert!(psm.scan_number.is_some());
            assert!(psm.spectrum_file.is_some());
        }
    }
}
```

- [ ] **Step 2: Update `loader/mod.rs`**

Add `pub mod pfind_tsv;` to the module declarations.

Add `PFindTsv` variant to `ResultFormat`:

```rust
pub enum ResultFormat {
    DiannParquet,
    GenericTsv,
    PFindTsv,
}
```

Update `ResultFormat::from_path()` — for `.tsv`/`.txt`, probe header before defaulting to `GenericTsv`:

```rust
"tsv" | "txt" => {
    if pfind_tsv::detect(path) {
        Ok(Self::PFindTsv)
    } else {
        Ok(Self::GenericTsv)
    }
}
```

Update `load_psms()` to handle the new variant:

```rust
ResultFormat::PFindTsv => pfind_tsv::load_pfind_tsv(path),
```

- [ ] **Step 3: Run tests**

```bash
cargo test -p protein-copilot-entrapment-analysis 2>&1
```

Expected: all existing tests pass + new `pfind_tsv` tests pass.

- [ ] **Step 4: Commit**

```bash
git add crates/entrapment-analysis/
git commit -m "feat(entrapment-analysis): add pFind TSV loader

- load_pfind_tsv() parses pFind TSV into Vec<UnifiedPsm>
- Modification parsing: pos,Name[Residue]; → (0-based, delta_mass)
- Protein IDs: /-separated → ;-separated
- MH+ → precursor_mz conversion
- Auto-detection via header probe in ResultFormat::from_path
- Unit tests for mods, detection, fixture loading"
```

---

### Task 4: Wire pFind TSV into MCP Server

**Files:**
- Modify: `crates/mcp-server/src/tools.rs`

- [ ] **Step 1: Update `import_search_results` tool**

In `crates/mcp-server/src/tools.rs`, make these changes:

**a) Update the format matching** (around line 3237):

Add `"pfind_tsv"` to the format match arm:

```rust
"pfind_tsv" => ImportFormat::PFindTsv,
```

Update the error message to include `pfind_tsv`:

```rust
"unknown format: '{other}'. Supported: auto, custom_json, diann_parquet, pfind_spectra, pfind_tsv"
```

**b) Add the PFindTsv parsing branch** (around line 3253, alongside CustomJson and DiannParquet):

```rust
ImportFormat::PFindTsv => {
    protein_copilot_result_import::pfind_tsv::PFindTsvParser
        .parse(&result_path, &unimod)
        .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?
}
```

**c) Handle scan matching skip for pFind TSV**:

After the `psms` are parsed (around line 3275), the existing code calls `match_scans`. For pFind TSV, scans are already matched. Replace the scan matching section with a conditional:

```rust
let all_scans_present = psms.iter().all(|p| p.matched_scan.is_some());

let match_report = if all_scans_present {
    // pFind TSV already has scan numbers — build MatchReport directly
    use std::collections::HashSet;
    let raw_names: HashSet<&str> = psms.iter().map(|p| p.raw_name.as_str()).collect();
    let mut per_file = HashMap::new();
    for raw in &raw_names {
        let count = psms.iter().filter(|p| p.raw_name == *raw).count();
        per_file.insert(
            raw.to_string(),
            protein_copilot_result_import::FileMatchStats {
                total: count,
                matched: count,
                ms2_count: 0, // unknown without scanning mzML
            },
        );
    }
    protein_copilot_result_import::MatchReport {
        total_psms: psms.len(),
        matched: psms.len(),
        unmatched: 0,
        median_rt_delta_min: 0.0,
        max_rt_delta_min: 0.0,
        per_file,
    }
} else {
    // Normal path: RT-based scan matching
    let config = ScanMatcherConfig {
        rt_tolerance_min: input.rt_tolerance_min,
        mzml_dir: mzml_dir.clone(),
    };
    match_scans(&mut psms, &config, &|path| {
        protein_copilot_spectrum_io::create_indexed_reader(path).map_err(|e| {
            protein_copilot_result_import::ResultImportError::SpectrumIo(format!(
                "failed to open {}: {e}",
                path.display(),
            ))
        })
    })
    .map_err(|e| mcp_err(ErrorCode::INTERNAL_ERROR, e.to_string()))?
};
```

**d) Update format_name match** (around line 3335):

```rust
ImportFormat::PFindTsv => "pfind_tsv",
```

**e) Update the tool description string** to mention pFind TSV support.

**f) Update `ImportSearchResultsInput` doc comment** for `format` field:

```rust
/// Result file format. 'auto' detects from extension. Options: auto, custom_json, diann_parquet, pfind_spectra, pfind_tsv.
```

- [ ] **Step 2: Update `classify_entrapment_hits` tool**

In the format matching section (around line 3788):

Add `"pfind_tsv"` variant:

```rust
Some("pfind_tsv") => ResultFormat::PFindTsv,
```

- [ ] **Step 3: Add PFindTsvParser to imports**

At the top of `tools.rs`, update the `protein_copilot_result_import` use block:

```rust
use protein_copilot_result_import::{
    converter::build_search_result,
    custom_json::CustomJsonParser,
    diann::DiannParser,
    pfind_tsv::PFindTsvParser,
    scan_matcher::{match_scans, ScanMatcherConfig},
    unimod::UnimodDb,
    ImportFormat, ImportResult, ResultParser,
};
```

- [ ] **Step 4: Build to verify**

```bash
cargo build -p protein-copilot-mcp-server 2>&1
```

Expected: builds without errors.

- [ ] **Step 5: Commit**

```bash
git add crates/mcp-server/
git commit -m "feat(mcp-server): wire pFind TSV format in import and entrapment tools

- import_search_results: add pfind_tsv format, skip scan matching when scans present
- classify_entrapment_hits: add pfind_tsv format option
- Auto-detection works for .tsv files with pFind headers"
```

---

### Task 5: Full Build and Test Verification

**Files:** None (verification only)

- [ ] **Step 1: Run full test suite**

```bash
cargo test --workspace 2>&1
```

Expected: all tests pass (894+ tests).

- [ ] **Step 2: Run clippy**

```bash
cargo clippy --workspace 2>&1
```

Expected: no new warnings from our changes.

- [ ] **Step 3: Test with real file via MCP**

```bash
# Quick smoke test: use the fixture file
printf '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{},"clientInfo":{"name":"test","version":"1.0"}}}\n{"jsonrpc":"2.0","method":"notifications/initialized"}\n{"jsonrpc":"2.0","id":2,"method":"tools/call","params":{"name":"import_search_results","arguments":{"result_file":"tests/fixtures/pfind_sample.tsv","mzml_dir":"/tmp","format":"pfind_tsv"}}}\n' | RUST_LOG=info timeout 10 cargo run -p protein-copilot-mcp-server 2>/dev/null
```

Expected: error about mzML files not found (expected — no mzML for test data), but parsing should succeed (look for "parsed pFind TSV" in stderr if RUST_LOG enabled).

- [ ] **Step 4: Commit any fixes**

If any issues found, fix and commit.
