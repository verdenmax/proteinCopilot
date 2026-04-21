# Entrapment v2: Edit Distance + Substitution Type Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace Hamming-only similarity with Levenshtein edit distance + k-mer pre-filtering to reduce false L4 classifications for indel/dipeptide homologs, and annotate substitution types (Q/K, isobaric dipeptide) within L2.

**Architecture:** Add a new `levenshtein.rs` module for Wagner-Fischer algorithm with alignment backtracking. Extend `TargetDigestIndex` with a k-mer inverted index (`find_similar()`) that pre-filters candidates before edit distance computation. Upgrade `classify_single()` to use edit distance instead of Hamming-only, with a `categorize_substitution()` function. Add `SubstitutionType` enum and 3 new fields to `ClassifiedPsm`. All output layers (TSV, HTML, MCP) get updated columns.

**Tech Stack:** Rust (no new external crates needed), `residue_mass()` from search-engine/chemistry.rs

**Spec:** `docs/superpowers/specs/2025-07-18-entrapment-v2-design.md`

**Baseline:** 56 passing tests in `protein-copilot-entrapment-analysis`

---

### Task 1: SubstitutionType enum and ClassifiedPsm extension (types.rs)

**Files:**
- Modify: `crates/entrapment-analysis/src/types.rs`

- [ ] **Step 1: Write tests for SubstitutionType serialization and ClassifiedPsm new fields**

Add to the existing `mod tests` block in `types.rs`:

```rust
#[test]
fn test_substitution_type_serde() {
    let st = SubstitutionType::QKSubstitution;
    let json = serde_json::to_string(&st).unwrap();
    assert_eq!(json, r#""QKSubstitution"#);

    let st2: SubstitutionType = serde_json::from_str(&json).unwrap();
    assert_eq!(st2, SubstitutionType::QKSubstitution);
}

#[test]
fn test_substitution_type_isobaric_dipeptide_serde() {
    let st = SubstitutionType::IsobaricDipeptide {
        single_residue: 'N',
        dipeptide: "GG".to_string(),
    };
    let json = serde_json::to_string(&st).unwrap();
    assert!(json.contains("IsobaricDipeptide"));
    assert!(json.contains("GG"));
    let st2: SubstitutionType = serde_json::from_str(&json).unwrap();
    assert_eq!(st2, st);
}

#[test]
fn test_substitution_type_display() {
    assert_eq!(SubstitutionType::None.as_str(), "None");
    assert_eq!(SubstitutionType::LIIsomer.as_str(), "LIIsomer");
    assert_eq!(SubstitutionType::QKSubstitution.as_str(), "QKSubstitution");
    assert_eq!(SubstitutionType::NearIsobaric.as_str(), "NearIsobaric");
    assert_eq!(SubstitutionType::Distinguishable.as_str(), "Distinguishable");
    let idb = SubstitutionType::IsobaricDipeptide {
        single_residue: 'N',
        dipeptide: "GG".to_string(),
    };
    assert_eq!(idb.as_str(), "IsobaricDipeptide");
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p protein-copilot-entrapment-analysis types::tests 2>&1 | tail -10`
Expected: FAIL — `SubstitutionType` not defined

- [ ] **Step 3: Implement SubstitutionType enum and extend ClassifiedPsm**

Add `SubstitutionType` enum before the `ClassifiedPsm` struct:

```rust
/// Substitution type annotation for L2 classified PSMs (v2).
///
/// Informational only — does not affect L0-L4 level assignment.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum SubstitutionType {
    /// No substitution detected (L0/L4) or not applicable.
    None,
    /// I↔L isomer (L1).
    LIIsomer,
    /// Q↔K substitution (Δm ≈ 36.4 mDa).
    QKSubstitution,
    /// Isobaric dipeptide substitution (N↔GG or Q↔AG).
    IsobaricDipeptide {
        single_residue: char,
        dipeptide: String,
    },
    /// Other near-isobaric substitution (|Δm| < threshold).
    NearIsobaric,
    /// Distinguishable substitution (|Δm| ≥ threshold).
    Distinguishable,
}

impl SubstitutionType {
    /// Returns a short label for display.
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::None => "None",
            Self::LIIsomer => "LIIsomer",
            Self::QKSubstitution => "QKSubstitution",
            Self::IsobaricDipeptide { .. } => "IsobaricDipeptide",
            Self::NearIsobaric => "NearIsobaric",
            Self::Distinguishable => "Distinguishable",
        }
    }
}

impl fmt::Display for SubstitutionType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}
```

Add 3 new fields to `ClassifiedPsm`:

```rust
pub struct ClassifiedPsm {
    // ... existing fields unchanged ...

    /// Substitution type annotation (v2). Informational only.
    pub substitution_type: SubstitutionType,
    /// Edit distance to best target (v2). Equals Hamming distance for same-length matches.
    pub edit_distance: Option<u32>,
    /// Alignment detail string (v2), e.g. "D0→N" or "ins:G@5".
    pub alignment_detail: Option<String>,
}
```

- [ ] **Step 4: Fix all compilation errors from the new required fields**

Every site that constructs `ClassifiedPsm` must add the 3 new fields. Search all files:

```bash
cargo build -p protein-copilot-entrapment-analysis 2>&1 | grep "missing field"
```

In each constructor, add:

```rust
substitution_type: SubstitutionType::None,
edit_distance: None,
alignment_detail: None,
```

Files affected:
- `similarity.rs`: all `ClassifiedPsm { ... }` blocks (lines 81-91, 96-106, 114-124, 126-135, 179-188, 197-206) — set `SubstitutionType::None` for L0/L4, `SubstitutionType::LIIsomer` for L1
- `output.rs` test helpers: `make_classified_psm()` (line ~257)
- `report.rs` test helpers: `make_psm()` (line ~155)

- [ ] **Step 5: Run tests to verify they pass**

Run: `cargo test -p protein-copilot-entrapment-analysis 2>&1 | tail -5`
Expected: all 59 tests pass (56 existing + 3 new)

- [ ] **Step 6: Export SubstitutionType from lib.rs**

In `crates/entrapment-analysis/src/lib.rs`, add to the `pub use types::` line:

```rust
pub use types::{
    ClassifiedPsm, DiscriminabilityLevel, EntrapmentSummary, LevelCounts, PsmGroup,
    SubstitutionType, UnifiedPsm,
};
```

- [ ] **Step 7: Commit**

```bash
git add crates/entrapment-analysis/src/types.rs crates/entrapment-analysis/src/lib.rs \
       crates/entrapment-analysis/src/similarity.rs crates/entrapment-analysis/src/output.rs \
       crates/entrapment-analysis/src/report.rs
git commit -m "feat(entrapment): add SubstitutionType enum and extend ClassifiedPsm with v2 fields"
```

---

### Task 2: Levenshtein edit distance module (levenshtein.rs)

**Files:**
- Create: `crates/entrapment-analysis/src/levenshtein.rs`
- Modify: `crates/entrapment-analysis/src/lib.rs` (add `pub mod levenshtein;`)

- [ ] **Step 1: Create levenshtein.rs with test stubs**

Create `crates/entrapment-analysis/src/levenshtein.rs`:

```rust
//! Levenshtein edit distance with alignment backtracking and delta-mass calculation.
//!
//! Used to find similar target peptides that may differ by insertions, deletions,
//! or substitutions — not just same-length Hamming distance.

use protein_copilot_search_engine::digest::residue_mass;

/// Result of a Levenshtein alignment between two sequences.
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    /// Edit distance (number of substitutions + insertions + deletions).
    pub edit_distance: u32,
    /// Mass difference: sum of mass changes from edit operations (Da).
    /// Positive means target is heavier than query.
    pub delta_mass_da: f64,
    /// Human-readable alignment detail, e.g. "D0→N,ins:G@5".
    pub alignment_detail: String,
}

/// Compute the Levenshtein edit distance between two sequences.
///
/// Returns only the distance (no alignment backtracking).
/// Uses O(min(m,n)) space with single-row optimization.
pub fn edit_distance(a: &str, b: &str) -> u32 {
    let a_bytes = a.as_bytes();
    let b_bytes = b.as_bytes();
    let (short, long) = if a_bytes.len() <= b_bytes.len() {
        (a_bytes, b_bytes)
    } else {
        (b_bytes, a_bytes)
    };

    let n = short.len();
    let m = long.len();

    let mut prev: Vec<u32> = (0..=(n as u32)).collect();
    let mut curr = vec![0u32; n + 1];

    for i in 1..=m {
        curr[0] = i as u32;
        for j in 1..=n {
            let cost = if long[i - 1] == short[j - 1] { 0 } else { 1 };
            curr[j] = (prev[j] + 1)          // deletion
                .min(curr[j - 1] + 1)         // insertion
                .min(prev[j - 1] + cost);     // substitution
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    prev[n]
}

/// Compute Levenshtein alignment with full backtracking.
///
/// Returns edit distance, mass difference, and a human-readable detail string.
/// The mass difference is computed from the alignment: for each edit operation,
/// the residue mass change is accumulated.
///
/// - Substitution A→B: Δm += mass(B) - mass(A)
/// - Insertion of B (in target, not in query): Δm += mass(B)
/// - Deletion of A (in query, not in target): Δm -= mass(A)
///
/// `a` is the query (trap peptide), `b` is the target peptide.
pub fn align(a: &str, b: &str) -> AlignmentResult {
    let a_chars: Vec<char> = a.chars().collect();
    let b_chars: Vec<char> = b.chars().collect();
    let m = a_chars.len();
    let n = b_chars.len();

    // Build full DP matrix for backtracking
    let mut dp = vec![vec![0u32; n + 1]; m + 1];
    for i in 0..=m {
        dp[i][0] = i as u32;
    }
    for j in 0..=n {
        dp[0][j] = j as u32;
    }
    for i in 1..=m {
        for j in 1..=n {
            let cost = if a_chars[i - 1] == b_chars[j - 1] { 0 } else { 1 };
            dp[i][j] = (dp[i - 1][j] + 1)
                .min(dp[i][j - 1] + 1)
                .min(dp[i - 1][j - 1] + cost);
        }
    }

    // Backtrack to reconstruct alignment
    let mut ops: Vec<String> = Vec::new();
    let mut delta_mass: f64 = 0.0;
    let mut i = m;
    let mut j = n;

    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let cost = if a_chars[i - 1] == b_chars[j - 1] { 0 } else { 1 };
            if dp[i][j] == dp[i - 1][j - 1] + cost {
                if cost == 1 {
                    // Substitution: a[i-1] → b[j-1]
                    let ca = a_chars[i - 1];
                    let cb = b_chars[j - 1];
                    let dm = match (residue_mass(ca), residue_mass(cb)) {
                        (Some(ma), Some(mb)) => mb - ma,
                        _ => 0.0,
                    };
                    delta_mass += dm;
                    ops.push(format!("{ca}{}→{cb}", i - 1));
                }
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if j > 0 && dp[i][j] == dp[i][j - 1] + 1 {
            // Insertion: b[j-1] is extra in target (not in query)
            let cb = b_chars[j - 1];
            if let Some(mb) = residue_mass(cb) {
                delta_mass += mb;
            }
            ops.push(format!("ins:{cb}@{i}"));
            j -= 1;
        } else if i > 0 && dp[i][j] == dp[i - 1][j] + 1 {
            // Deletion: a[i-1] is in query but not in target
            let ca = a_chars[i - 1];
            if let Some(ma) = residue_mass(ca) {
                delta_mass -= ma;
            }
            ops.push(format!("del:{ca}@{}", i - 1));
            i -= 1;
        } else {
            // Should not reach here, but break to avoid infinite loop
            break;
        }
    }

    ops.reverse();
    let detail = ops.join(",");

    AlignmentResult {
        edit_distance: dp[m][n],
        delta_mass_da: delta_mass,
        alignment_detail: detail,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- edit_distance tests ---

    #[test]
    fn test_edit_distance_identical() {
        assert_eq!(edit_distance("PEPTIDEK", "PEPTIDEK"), 0);
    }

    #[test]
    fn test_edit_distance_one_substitution() {
        assert_eq!(edit_distance("PEPTIDEK", "PEPTLDEK"), 1);
    }

    #[test]
    fn test_edit_distance_one_insertion() {
        // "PEPTIDEK" vs "PEPGTIDEK" — G inserted at pos 3
        assert_eq!(edit_distance("PEPTIDEK", "PEPGTIDEK"), 1);
    }

    #[test]
    fn test_edit_distance_one_deletion() {
        // "PEPTIDEK" vs "PEPTDEK" — I deleted
        assert_eq!(edit_distance("PEPTIDEK", "PEPTDEK"), 1);
    }

    #[test]
    fn test_edit_distance_two_edits() {
        // Two substitutions
        assert_eq!(edit_distance("AGCDEK", "ANCDER"), 2);
    }

    #[test]
    fn test_edit_distance_empty() {
        assert_eq!(edit_distance("", ""), 0);
        assert_eq!(edit_distance("ABC", ""), 3);
        assert_eq!(edit_distance("", "ABC"), 3);
    }

    #[test]
    fn test_edit_distance_equals_hamming_for_same_length() {
        // For equal-length strings, edit distance <= hamming distance
        // (edit distance can be lower if substitution = del+ins is cheaper)
        // For single substitution, they're equal
        assert_eq!(edit_distance("DGFLLDGFPR", "NGFLLDGFPR"), 1);
    }

    // --- align tests ---

    #[test]
    fn test_align_identical() {
        let r = align("PEPTIDEK", "PEPTIDEK");
        assert_eq!(r.edit_distance, 0);
        assert!(r.delta_mass_da.abs() < 1e-6);
        assert_eq!(r.alignment_detail, "");
    }

    #[test]
    fn test_align_one_substitution_dn() {
        // D→N at position 0
        let r = align("DGFLLDGFPR", "NGFLLDGFPR");
        assert_eq!(r.edit_distance, 1);
        // N(114.042927) - D(115.026943) = -0.984016
        assert!((r.delta_mass_da - (-0.984016)).abs() < 0.001);
        assert_eq!(r.alignment_detail, "D0→N");
    }

    #[test]
    fn test_align_insertion() {
        // query: "PEPNDE", target: "PEPGGDE" — N→GG (one deletion + one insertion = 2 edits via Levenshtein)
        // Actually: "PEPNDE" (6) vs "PEPGGDE" (7), edit distance should be 2
        let r = align("PEPNDE", "PEPGGDE");
        assert_eq!(r.edit_distance, 2);
        assert!(r.alignment_detail.len() > 0);
    }

    #[test]
    fn test_align_qk_substitution() {
        // Q→K at position 4
        let r = align("PEPTQDEK", "PEPTKDEK");
        assert_eq!(r.edit_distance, 1);
        // K(128.094963) - Q(128.058578) = 0.036385
        assert!((r.delta_mass_da - 0.036385).abs() < 0.001);
        assert_eq!(r.alignment_detail, "Q4→K");
    }

    #[test]
    fn test_align_empty_strings() {
        let r = align("", "");
        assert_eq!(r.edit_distance, 0);
        assert_eq!(r.alignment_detail, "");

        let r2 = align("", "ABC");
        assert_eq!(r2.edit_distance, 3);

        let r3 = align("ABC", "");
        assert_eq!(r3.edit_distance, 3);
    }
}
```

- [ ] **Step 2: Register the module in lib.rs**

Add `pub mod levenshtein;` to `crates/entrapment-analysis/src/lib.rs`.

- [ ] **Step 3: Run tests to verify they pass**

Run: `cargo test -p protein-copilot-entrapment-analysis levenshtein::tests 2>&1 | tail -10`
Expected: all 11 levenshtein tests pass

- [ ] **Step 4: Commit**

```bash
git add crates/entrapment-analysis/src/levenshtein.rs crates/entrapment-analysis/src/lib.rs
git commit -m "feat(entrapment): add Levenshtein edit distance module with alignment backtracking"
```

---

### Task 3: SimilarityConfig v2 fields (config.rs)

**Files:**
- Modify: `crates/entrapment-analysis/src/config.rs`

- [ ] **Step 1: Write tests for new config fields and serde alias**

Add to the existing `mod tests` in `config.rs`:

```rust
#[test]
fn test_v2_defaults() {
    let cfg = EntrapmentConfig::from_yaml_str(MINIMAL_YAML).expect("parse");
    assert_eq!(cfg.similarity.len_tolerance, 2);
    assert!(cfg.similarity.enable_dipeptide_check);
    assert!(cfg.similarity.enable_qk_detection);
    // New field name
    assert!((cfg.similarity.delta_mass_threshold_da - 1.0).abs() < f64::EPSILON);
}

#[test]
fn test_delta_mz_alias_still_works() {
    let yaml = r#"
version: 1
target:
  rules:
    - type: AccessionContains
      any_of: ["_HUMAN"]
trap:
  rules:
    - type: AccessionContains
      any_of: ["_YEAST"]
similarity:
  delta_mz_threshold_da: 0.5
"#;
    let cfg = EntrapmentConfig::from_yaml_str(yaml).expect("parse with alias");
    assert!((cfg.similarity.delta_mass_threshold_da - 0.5).abs() < f64::EPSILON);
}

#[test]
fn test_v2_explicit_overrides() {
    let yaml = r#"
version: 1
target:
  rules:
    - type: AccessionContains
      any_of: ["_HUMAN"]
trap:
  rules:
    - type: AccessionContains
      any_of: ["_YEAST"]
similarity:
  max_mismatches: 3
  delta_mass_threshold_da: 0.8
  len_tolerance: 3
  enable_dipeptide_check: false
  enable_qk_detection: false
"#;
    let cfg = EntrapmentConfig::from_yaml_str(yaml).expect("parse");
    assert_eq!(cfg.similarity.max_mismatches, 3);
    assert!((cfg.similarity.delta_mass_threshold_da - 0.8).abs() < f64::EPSILON);
    assert_eq!(cfg.similarity.len_tolerance, 3);
    assert!(!cfg.similarity.enable_dipeptide_check);
    assert!(!cfg.similarity.enable_qk_detection);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p protein-copilot-entrapment-analysis config::tests 2>&1 | tail -10`
Expected: FAIL — `len_tolerance`, `delta_mass_threshold_da` not defined

- [ ] **Step 3: Implement SimilarityConfig changes**

Replace the `SimilarityConfig` struct and its defaults:

```rust
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimilarityConfig {
    /// Maximum edit distance (substitutions + insertions + deletions).
    /// Semantically replaces "mismatches" for v2 but field name kept for YAML compat.
    #[serde(default = "default_max_mismatches")]
    pub max_mismatches: u16,

    /// Mass-difference threshold (Da) separating L2 (near-isobaric) from L3.
    #[serde(
        default = "default_delta_mass_threshold_da",
        alias = "delta_mz_threshold_da"
    )]
    pub delta_mass_threshold_da: f64,

    /// Whether both ends of a peptide must be tryptic.
    #[serde(default = "default_require_tryptic_ends")]
    pub require_tryptic_ends: bool,

    /// Maximum number of missed cleavages to allow.
    #[serde(default = "default_max_missed_cleavages")]
    pub max_missed_cleavages: u32,

    /// Length tolerance: search target peptides within len ± len_tolerance.
    #[serde(default = "default_len_tolerance")]
    pub len_tolerance: usize,

    /// Enable isobaric dipeptide detection (N↔GG, Q↔AG).
    #[serde(default = "default_true")]
    pub enable_dipeptide_check: bool,

    /// Enable Q/K near-isobaric substitution detection.
    #[serde(default = "default_true")]
    pub enable_qk_detection: bool,
}

fn default_delta_mass_threshold_da() -> f64 {
    1.0
}
fn default_len_tolerance() -> usize {
    2
}
fn default_true() -> bool {
    true
}
```

Update the `Default` impl:

```rust
impl Default for SimilarityConfig {
    fn default() -> Self {
        Self {
            max_mismatches: default_max_mismatches(),
            delta_mass_threshold_da: default_delta_mass_threshold_da(),
            require_tryptic_ends: default_require_tryptic_ends(),
            max_missed_cleavages: default_max_missed_cleavages(),
            len_tolerance: default_len_tolerance(),
            enable_dipeptide_check: default_true(),
            enable_qk_detection: default_true(),
        }
    }
}
```

Remove the old `default_delta_mz_threshold_da()` function and rename `delta_mz_threshold_da` → `delta_mass_threshold_da` in the struct.

- [ ] **Step 4: Fix all references to the old field name**

Search and replace `delta_mz_threshold_da` → `delta_mass_threshold_da` across the crate:

```bash
grep -rn "delta_mz_threshold_da" crates/entrapment-analysis/
```

Files affected:
- `similarity.rs` line 191: `config.delta_mz_threshold_da` → `config.delta_mass_threshold_da`
- `config.rs` tests: update `FULL_YAML` constant and `test_parse_full_config` assertion

Also fix in `crates/mcp-server/src/tools.rs` if referenced.

- [ ] **Step 5: Run tests to verify they pass**

Run: `cargo test -p protein-copilot-entrapment-analysis 2>&1 | tail -5`
Expected: all tests pass (56 original + 3 type tests + 3 config tests = 62)

- [ ] **Step 6: Commit**

```bash
git add crates/entrapment-analysis/src/config.rs crates/entrapment-analysis/src/similarity.rs
git commit -m "feat(entrapment): add v2 config fields, rename delta_mz to delta_mass with serde alias"
```

---

### Task 4: k-mer index and find_similar() in TargetDigestIndex (digest.rs)

**Files:**
- Modify: `crates/entrapment-analysis/src/digest.rs`

- [ ] **Step 1: Write tests for k-mer index and find_similar()**

Add to the existing `mod tests` in `digest.rs`:

```rust
use crate::config::SimilarityConfig;

#[test]
fn test_kmer_index_built_on_from_fasta() {
    use std::io::Write;
    let mut f = tempfile::NamedTempFile::new().expect("create temp file");
    write!(f, ">sp|P001|TEST_HUMAN Test\nPEPTIDEKANSTHERPEPTIDERLASTPART\n").unwrap();
    let config = SimilarityConfig::default();
    let idx = TargetDigestIndex::from_fasta(f.path(), config.max_missed_cleavages)
        .expect("build index");
    // k-mer index should be populated
    assert!(!idx.kmer_index.is_empty());
    assert!(!idx.all_peptides.is_empty());
    assert!(idx.kmer_k >= 2);
}

#[test]
fn test_find_similar_exact_match_excluded() {
    use std::io::Write;
    let mut f = tempfile::NamedTempFile::new().expect("create temp file");
    write!(f, ">sp|P001|TEST_HUMAN Test\nPEPTIDEKANSTHERR\n").unwrap();
    let config = SimilarityConfig::default();
    let idx = TargetDigestIndex::from_fasta(f.path(), config.max_missed_cleavages).unwrap();
    // "PEPTIDEK" is in the index; find_similar should return matches but exact match filtered out
    let matches = idx.find_similar("PEPTIDEK", 2, 2, &config);
    // Exact match should NOT appear in find_similar results (handled by L0 path)
    for m in &matches {
        assert_ne!(m.target_peptide, "PEPTIDEK");
    }
}

#[test]
fn test_find_similar_one_substitution() {
    use std::io::Write;
    let mut f = tempfile::NamedTempFile::new().expect("create temp file");
    // Target has "NGFLLDGFPR", query is "DGFLLDGFPR" (D→N, edit=1)
    write!(f, ">sp|P001|TEST_HUMAN Test\nNGFLLDGFPR\n").unwrap();
    let config = SimilarityConfig::default();
    let idx = TargetDigestIndex::from_fasta(f.path(), config.max_missed_cleavages).unwrap();
    let matches = idx.find_similar("DGFLLDGFPR", 2, 2, &config);
    assert_eq!(matches.len(), 1);
    assert_eq!(matches[0].target_peptide, "NGFLLDGFPR");
    assert_eq!(matches[0].edit_distance, 1);
    assert!((matches[0].delta_mass_da.abs() - 0.984016).abs() < 0.001);
}

#[test]
fn test_find_similar_different_length_indel() {
    use std::io::Write;
    let mut f = tempfile::NamedTempFile::new().expect("create temp file");
    // Target: "PEPGGDEK" (8 chars), Query: "PEPNDEK" (7 chars)
    // N→GG is edit distance 2 (del N, ins G, ins G? Actually: depends on alignment)
    write!(f, ">sp|P001|TEST_HUMAN Test\nPEPGGDEKKAAAAAR\n").unwrap();
    let config = SimilarityConfig::default();
    let idx = TargetDigestIndex::from_fasta(f.path(), config.max_missed_cleavages).unwrap();
    let matches = idx.find_similar("PEPNDEK", 2, 2, &config);
    // Should find "PEPGGDEK" as a candidate (len 8 is within len_tolerance=2 of len 7)
    let found = matches.iter().any(|m| m.target_peptide == "PEPGGDEK");
    assert!(found, "should find PEPGGDEK as similar to PEPNDEK; found: {:?}", matches);
}

#[test]
fn test_find_similar_no_match_beyond_tolerance() {
    use std::io::Write;
    let mut f = tempfile::NamedTempFile::new().expect("create temp file");
    write!(f, ">sp|P001|TEST_HUMAN Test\nWWWWWWWWR\n").unwrap();
    let config = SimilarityConfig::default();
    let idx = TargetDigestIndex::from_fasta(f.path(), config.max_missed_cleavages).unwrap();
    // Completely different peptide
    let matches = idx.find_similar("PEPTIDEK", 2, 2, &config);
    assert!(matches.is_empty());
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p protein-copilot-entrapment-analysis digest::tests 2>&1 | tail -10`
Expected: FAIL — `kmer_index`, `find_similar` not defined

- [ ] **Step 3: Implement k-mer index in TargetDigestIndex**

Add the new struct `SimilarityMatch` and k-mer fields to `TargetDigestIndex`:

```rust
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use crate::config::SimilarityConfig;
use crate::levenshtein;
use crate::types::SubstitutionType;

/// A match found by the similarity search.
#[derive(Debug, Clone)]
pub struct SimilarityMatch {
    pub target_peptide: String,
    pub target_protein: String,
    pub edit_distance: u32,
    pub delta_mass_da: f64,
    pub alignment_detail: String,
    pub substitution_type: SubstitutionType,
}
```

Extend `TargetDigestIndex`:

```rust
pub struct TargetDigestIndex {
    // existing fields...
    pub by_length: HashMap<usize, Vec<TargetPeptide>>,
    pub exact_set: HashSet<String>,
    pub normalized_set: HashSet<String>,
    pub exact_to_protein: HashMap<String, String>,
    pub normalized_to_original: HashMap<String, (String, String)>,

    // v2: k-mer inverted index
    kmer_index: HashMap<u64, Vec<u32>>,
    all_peptides: Vec<TargetPeptide>,
    kmer_k: usize,
}
```

Add k-mer helper functions:

```rust
/// Hash a k-mer string to u64.
fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut hasher = DefaultHasher::new();
    kmer.hash(&mut hasher);
    hasher.finish()
}

/// Extract all k-mers from a sequence.
fn extract_kmers(seq: &[u8], k: usize) -> Vec<u64> {
    if seq.len() < k {
        return Vec::new();
    }
    (0..=(seq.len() - k))
        .map(|i| hash_kmer(&seq[i..i + k]))
        .collect()
}
```

In `from_fasta()`, after building existing structures, build the k-mer index:

```rust
// Determine k for pigeonhole guarantee: k = min_len / (max_edit + 1)
let min_peptide_len = 6usize; // our min digest length
let max_edit = 2u16; // default max_mismatches
let kmer_k = min_peptide_len / (max_edit as usize + 1);
let kmer_k = kmer_k.max(2); // at least 2

// Build flat peptide array and k-mer index
let mut all_peptides: Vec<TargetPeptide> = Vec::new();
let mut kmer_index: HashMap<u64, Vec<u32>> = HashMap::new();

for peptides in by_length.values() {
    for tp in peptides {
        let pid = all_peptides.len() as u32;
        let kmers = extract_kmers(tp.sequence.as_bytes(), kmer_k);
        for kh in kmers {
            kmer_index.entry(kh).or_default().push(pid);
        }
        all_peptides.push(tp.clone());
    }
}

// Deduplicate each posting list
for list in kmer_index.values_mut() {
    list.sort_unstable();
    list.dedup();
}
```

And return with the new fields:

```rust
Ok(Self {
    by_length,
    exact_set,
    normalized_set,
    exact_to_protein,
    normalized_to_original,
    kmer_index,
    all_peptides,
    kmer_k,
})
```

- [ ] **Step 4: Implement find_similar() method**

```rust
impl TargetDigestIndex {
    /// Find target peptides within `max_edit_dist` edit distance of `query`.
    ///
    /// Uses k-mer pre-filtering (pigeonhole principle) to reduce candidates,
    /// then computes full Levenshtein distance + alignment for survivors.
    pub fn find_similar(
        &self,
        query: &str,
        max_edit_dist: u16,
        len_tolerance: usize,
        config: &SimilarityConfig,
    ) -> Vec<SimilarityMatch> {
        let query_len = query.len();
        let min_len = query_len.saturating_sub(len_tolerance);
        let max_len = query_len + len_tolerance;

        // Extract k-mers from query
        let query_kmers = extract_kmers(query.as_bytes(), self.kmer_k);
        if query_kmers.is_empty() {
            return Vec::new();
        }

        // Collect candidate peptide IDs from k-mer hits
        let mut candidate_ids: HashSet<u32> = HashSet::new();
        for kh in &query_kmers {
            if let Some(ids) = self.kmer_index.get(kh) {
                for &id in ids {
                    candidate_ids.insert(id);
                }
            }
        }

        // Filter and compute edit distance
        let mut results = Vec::new();
        for &pid in &candidate_ids {
            let tp = &self.all_peptides[pid as usize];

            // Length filter
            if tp.sequence.len() < min_len || tp.sequence.len() > max_len {
                continue;
            }

            // Skip exact matches (handled by L0 path)
            if tp.sequence == query {
                continue;
            }

            // Quick edit distance check
            let dist = levenshtein::edit_distance(query, &tp.sequence);
            if dist > max_edit_dist as u32 {
                continue;
            }

            // Full alignment for survivors
            let alignment = levenshtein::align(query, &tp.sequence);

            results.push(SimilarityMatch {
                target_peptide: tp.sequence.clone(),
                target_protein: tp.protein_accession.clone(),
                edit_distance: alignment.edit_distance,
                delta_mass_da: alignment.delta_mass_da,
                alignment_detail: alignment.alignment_detail,
                substitution_type: SubstitutionType::None, // categorized later
            });
        }

        // Sort by edit distance, then by |delta_mass|
        results.sort_by(|a, b| {
            a.edit_distance.cmp(&b.edit_distance)
                .then_with(|| a.delta_mass_da.abs().partial_cmp(&b.delta_mass_da.abs())
                    .unwrap_or(std::cmp::Ordering::Equal))
        });

        results
    }
}
```

- [ ] **Step 5: Update empty_for_test() to include new fields**

```rust
#[cfg(test)]
pub fn empty_for_test() -> Self {
    Self {
        by_length: HashMap::new(),
        exact_set: HashSet::new(),
        normalized_set: HashSet::new(),
        exact_to_protein: HashMap::new(),
        normalized_to_original: HashMap::new(),
        kmer_index: HashMap::new(),
        all_peptides: Vec::new(),
        kmer_k: 2,
    }
}
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `cargo test -p protein-copilot-entrapment-analysis 2>&1 | tail -5`
Expected: all tests pass

- [ ] **Step 7: Commit**

```bash
git add crates/entrapment-analysis/src/digest.rs
git commit -m "feat(entrapment): add k-mer inverted index and find_similar() to TargetDigestIndex"
```

---

### Task 5: Upgrade classify_single() with edit distance + substitution type (similarity.rs)

**Files:**
- Modify: `crates/entrapment-analysis/src/similarity.rs`

- [ ] **Step 1: Write tests for the upgraded classification logic**

Add to the existing `mod tests` in `similarity.rs`:

```rust
use crate::types::SubstitutionType;

#[test]
fn test_classify_trap_l1_gets_li_isomer_type() {
    let psm = make_psm("PEPTIDEK");
    let mut index = TargetDigestIndex::empty_for_test();
    let norm = "PEPTLDEK";
    index.normalized_set.insert(norm.to_owned());
    index.normalized_to_original.insert(
        norm.to_owned(),
        ("PEPTLDEK".to_owned(), "P00002".to_owned()),
    );
    let config = SimilarityConfig::default();
    let result = classify_single(&psm, PsmGroup::Trap, &index, &config);
    assert_eq!(result.level, DiscriminabilityLevel::L1);
    assert_eq!(result.substitution_type, SubstitutionType::LIIsomer);
}

#[test]
fn test_classify_l2_near_isobaric_has_edit_distance() {
    let psm = make_psm("DGFLLDGFPR");
    let mut index = TargetDigestIndex::empty_for_test();
    index.by_length.insert(
        10,
        vec![TargetPeptide {
            sequence: "NGFLLDGFPR".to_owned(),
            protein_accession: "P00003".to_owned(),
            neutral_mass: 0.0,
        }],
    );
    let config = SimilarityConfig::default();
    let result = classify_single(&psm, PsmGroup::Trap, &index, &config);
    assert_eq!(result.level, DiscriminabilityLevel::L2);
    assert_eq!(result.edit_distance, Some(1));
    assert!(result.alignment_detail.is_some());
}

#[test]
fn test_classify_l0_has_none_substitution_type() {
    let psm = make_psm("PEPTIDEK");
    let mut index = TargetDigestIndex::empty_for_test();
    index.exact_set.insert("PEPTIDEK".to_owned());
    index.exact_to_protein.insert("PEPTIDEK".to_owned(), "P00001".to_owned());
    let config = SimilarityConfig::default();
    let result = classify_single(&psm, PsmGroup::Trap, &index, &config);
    assert_eq!(result.level, DiscriminabilityLevel::L0);
    assert_eq!(result.substitution_type, SubstitutionType::None);
}

#[test]
fn test_classify_l4_has_none_substitution_type() {
    let psm = make_psm("DGFLLDGFPR");
    let index = TargetDigestIndex::empty_for_test();
    let config = SimilarityConfig::default();
    let result = classify_single(&psm, PsmGroup::Trap, &index, &config);
    assert_eq!(result.level, DiscriminabilityLevel::L4);
    assert_eq!(result.substitution_type, SubstitutionType::None);
}
```

- [ ] **Step 2: Run tests to verify the new assertions fail (existing fields don't have substitution_type)**

Run: `cargo test -p protein-copilot-entrapment-analysis similarity::tests 2>&1 | tail -15`
Expected: Existing tests still pass, but new tests fail on substitution_type assertion.

- [ ] **Step 3: Add categorize_substitution() and helper functions**

Add to `similarity.rs`:

```rust
use crate::levenshtein;
use crate::types::SubstitutionType;

/// Known isobaric dipeptide pairs: (single_residue, dipeptide).
const ISOBARIC_DIPEPTIDES: &[(char, &str)] = &[
    ('N', "GG"), // 114.04293 Da
    ('Q', "AG"), // 128.05858 Da
];

/// Categorize the type of substitution between trap and best target match.
fn categorize_substitution(
    trap: &str,
    best_target: &str,
    edit_dist: u32,
    delta_mass: f64,
    alignment_detail: &str,
    config: &SimilarityConfig,
) -> SubstitutionType {
    let len_diff = (trap.len() as i64 - best_target.len() as i64).unsigned_abs() as usize;

    // Equal length, single substitution → check Q↔K
    if len_diff == 0 && edit_dist == 1 && config.enable_qk_detection {
        if is_qk_substitution(trap, best_target) {
            return SubstitutionType::QKSubstitution;
        }
    }

    // Length diff of 1, edit distance ≤ 2 → check isobaric dipeptide (N↔GG, Q↔AG)
    if len_diff == 1 && config.enable_dipeptide_check {
        if let Some((single, dipeptide)) = check_isobaric_dipeptide(trap, best_target) {
            return SubstitutionType::IsobaricDipeptide {
                single_residue: single,
                dipeptide,
            };
        }
    }

    // General categorization by delta mass
    if delta_mass.abs() < config.delta_mass_threshold_da {
        SubstitutionType::NearIsobaric
    } else {
        SubstitutionType::Distinguishable
    }
}

/// Check if the single differing position is a Q↔K swap.
fn is_qk_substitution(a: &str, b: &str) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let diffs: Vec<(char, char)> = a.chars().zip(b.chars())
        .filter(|(ca, cb)| ca != cb)
        .collect();
    if diffs.len() != 1 {
        return false;
    }
    let (ca, cb) = diffs[0];
    (ca == 'Q' && cb == 'K') || (ca == 'K' && cb == 'Q')
}

/// Check if the alignment represents an isobaric dipeptide substitution.
///
/// Looks for N↔GG or Q↔AG patterns where one sequence has the single residue
/// and the other has the dipeptide at the same position.
fn check_isobaric_dipeptide(shorter: &str, longer: &str) -> Option<(char, String)> {
    let (s, l) = if shorter.len() < longer.len() {
        (shorter, longer)
    } else if shorter.len() > longer.len() {
        (longer, shorter)
    } else {
        return None; // same length, not a dipeptide substitution
    };

    if l.len() != s.len() + 1 {
        return None;
    }

    // Find the position where they diverge
    let s_chars: Vec<char> = s.chars().collect();
    let l_chars: Vec<char> = l.chars().collect();

    for i in 0..s_chars.len() {
        if s_chars[i] != l_chars[i] {
            // Check if s[i] maps to l[i..i+2]
            if i + 1 < l_chars.len() {
                let single = s_chars[i];
                let dipeptide: String = l_chars[i..=i + 1].iter().collect();
                // Check if the rest matches (shifted by 1)
                let rest_matches = s_chars[i + 1..].iter()
                    .zip(l_chars[i + 2..].iter())
                    .all(|(a, b)| a == b);
                if rest_matches {
                    // Check against known isobaric pairs
                    for &(known_single, known_di) in ISOBARIC_DIPEPTIDES {
                        if (single == known_single && dipeptide == known_di) ||
                           // Check reverse: s has dipeptide part? No — s is shorter.
                           false
                        {
                            return Some((known_single, known_di.to_string()));
                        }
                    }
                }
            }
            return None;
        }
    }

    // Divergence at the very end: s is a prefix of l, extra char at end
    // This would be a simple insertion, not a dipeptide swap
    None
}
```

- [ ] **Step 4: Upgrade classify_single() to use edit distance for L2/L3**

Replace the L2/L3/L4 hamming scan section (after the L1 check) with:

```rust
    // --- L2/L3/L4: edit-distance scan (v2) --------------------------------

    // Phase A: same-length Hamming scan (fast path, backward compatible)
    let candidates = index.peptides_of_length(psm.peptide.len());
    let mut best_mm: u16 = u16::MAX;
    let mut best_dm: f64 = f64::MAX;
    let mut best_dp = String::new();
    let mut best_seq: Option<&str> = None;
    let mut best_prot: Option<&str> = None;

    for target in candidates {
        let (mm, dm, dp) = match hamming_diff(&psm.peptide, &target.sequence) {
            Some(v) => v,
            None => continue,
        };
        if mm == 0 { continue; }
        if mm > config.max_mismatches { continue; }
        if is_only_li_substitution(&psm.peptide, &target.sequence) { continue; }

        let abs_dm = dm.abs();
        if mm < best_mm || (mm == best_mm && abs_dm < best_dm) {
            best_mm = mm;
            best_dm = abs_dm;
            best_dp = dp;
            best_seq = Some(&target.sequence);
            best_prot = Some(&target.protein_accession);
        }
    }

    // Phase B: cross-length edit distance scan (v2 upgrade)
    // Search target peptides of different lengths via k-mer index
    let cross_matches = index.find_similar(
        &psm.peptide,
        config.max_mismatches,
        config.len_tolerance,
        config,
    );

    // Find best cross-length match (only consider matches with different length)
    let best_cross = cross_matches.iter()
        .filter(|m| m.target_peptide.len() != psm.peptide.len())
        .min_by(|a, b| {
            a.edit_distance.cmp(&b.edit_distance)
                .then_with(|| a.delta_mass_da.abs().partial_cmp(&b.delta_mass_da.abs())
                    .unwrap_or(std::cmp::Ordering::Equal))
        });

    // Determine overall best: compare Hamming best vs cross-length best
    enum BestMatch<'a> {
        Hamming { mm: u16, dm: f64, dp: String, seq: &'a str, prot: &'a str },
        CrossLength(SimilarityMatch),
        None,
    }

    let overall_best = match (best_mm < u16::MAX, best_cross) {
        (true, Some(cross)) => {
            // Compare: Hamming mm vs cross edit_distance
            if (best_mm as u32) <= cross.edit_distance
                || ((best_mm as u32) == cross.edit_distance && best_dm <= cross.delta_mass_da.abs())
            {
                BestMatch::Hamming { mm: best_mm, dm: best_dm, dp: best_dp, seq: best_seq.unwrap(), prot: best_prot.unwrap() }
            } else {
                BestMatch::CrossLength(cross.clone())
            }
        }
        (true, None) => BestMatch::Hamming { mm: best_mm, dm: best_dm, dp: best_dp, seq: best_seq.unwrap(), prot: best_prot.unwrap() },
        (false, Some(cross)) => BestMatch::CrossLength(cross.clone()),
        (false, None) => BestMatch::None,
    };

    match overall_best {
        BestMatch::None => ClassifiedPsm {
            psm: psm.clone(),
            group,
            level: DiscriminabilityLevel::L4,
            best_target_peptide: None,
            best_target_protein: None,
            mismatches: None,
            delta_mass_da: None,
            diff_positions: None,
            substitution_type: SubstitutionType::None,
            edit_distance: None,
            alignment_detail: None,
        },
        BestMatch::Hamming { mm, dm, dp, seq, prot } => {
            let sub_type = categorize_substitution(
                &psm.peptide, seq, mm as u32, dm, &dp, config,
            );
            let alignment = levenshtein::align(&psm.peptide, seq);
            let level = if dm.abs() < config.delta_mass_threshold_da {
                DiscriminabilityLevel::L2
            } else {
                DiscriminabilityLevel::L3
            };
            ClassifiedPsm {
                psm: psm.clone(),
                group,
                level,
                best_target_peptide: Some(seq.to_owned()),
                best_target_protein: Some(prot.to_owned()),
                mismatches: Some(mm),
                delta_mass_da: Some(dm),
                diff_positions: Some(dp),
                substitution_type: sub_type,
                edit_distance: Some(mm as u32),
                alignment_detail: Some(alignment.alignment_detail),
            }
        }
        BestMatch::CrossLength(cross) => {
            let sub_type = categorize_substitution(
                &psm.peptide, &cross.target_peptide,
                cross.edit_distance, cross.delta_mass_da,
                &cross.alignment_detail, config,
            );
            let level = if cross.delta_mass_da.abs() < config.delta_mass_threshold_da {
                DiscriminabilityLevel::L2
            } else {
                DiscriminabilityLevel::L3
            };
            ClassifiedPsm {
                psm: psm.clone(),
                group,
                level,
                best_target_peptide: Some(cross.target_peptide),
                best_target_protein: Some(cross.target_protein),
                mismatches: None, // Not applicable for cross-length
                delta_mass_da: Some(cross.delta_mass_da),
                diff_positions: None,
                substitution_type: sub_type,
                edit_distance: Some(cross.edit_distance),
                alignment_detail: Some(cross.alignment_detail),
            }
        }
    }
```

- [ ] **Step 5: Add import for SimilarityMatch**

At the top of `similarity.rs`:

```rust
use crate::digest::SimilarityMatch;
```

- [ ] **Step 6: Run all tests to verify correctness**

Run: `cargo test -p protein-copilot-entrapment-analysis 2>&1 | tail -10`
Expected: All tests pass (existing + new). The `best_mm == 1` L2 condition is now `dm < threshold`.

- [ ] **Step 7: Commit**

```bash
git add crates/entrapment-analysis/src/similarity.rs
git commit -m "feat(entrapment): upgrade classify_single with edit distance, substitution type, and cross-length matching"
```

---

### Task 6: TSV output — new columns (output.rs)

**Files:**
- Modify: `crates/entrapment-analysis/src/output.rs`

- [ ] **Step 1: Write test for new TSV columns**

Add to `output.rs` tests:

```rust
#[test]
fn test_write_classified_tsv_v2_columns() {
    use crate::types::SubstitutionType;
    let dir = tempfile::tempdir().expect("create temp dir");
    let path = dir.path().join("classified_v2.tsv");

    let psm = ClassifiedPsm {
        psm: UnifiedPsm {
            peptide: "PEPQDEK".to_owned(),
            charge: Some(2),
            precursor_mz: Some(500.0),
            retention_time: Some(10.0),
            scan_number: Some(100),
            spectrum_file: Some("test.raw".to_owned()),
            protein_ids: "sp|P001|TRAP_YEAST".to_owned(),
            q_value: Some(0.01),
        },
        group: PsmGroup::Trap,
        level: DiscriminabilityLevel::L2,
        best_target_peptide: Some("PEPKDEK".to_owned()),
        best_target_protein: Some("sp|P002|TARGET_HUMAN".to_owned()),
        mismatches: Some(1),
        delta_mass_da: Some(0.036385),
        diff_positions: Some("[3:Q->K]".to_owned()),
        substitution_type: SubstitutionType::QKSubstitution,
        edit_distance: Some(1),
        alignment_detail: Some("Q3→K".to_owned()),
    };

    write_classified_tsv(&[psm], &path).expect("write TSV");

    let content = std::fs::read_to_string(&path).expect("read TSV");
    let lines: Vec<&str> = content.lines().collect();
    // Header should contain new columns
    assert!(lines[0].contains("substitution_type"));
    assert!(lines[0].contains("edit_distance"));
    assert!(lines[0].contains("alignment_detail"));
    // Data row should contain values
    assert!(lines[1].contains("QKSubstitution"));
    assert!(lines[1].contains("Q3→K") || lines[1].contains("Q3"));
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p protein-copilot-entrapment-analysis output::tests::test_write_classified_tsv_v2_columns 2>&1 | tail -10`
Expected: FAIL — header doesn't contain new columns yet

- [ ] **Step 3: Add new columns to write_classified_tsv()**

In the `wtr.write_record` for the header, append:

```rust
"substitution_type",
"edit_distance",
"alignment_detail",
```

In the data row write, append:

```rust
&cp.substitution_type.to_string(),
&opt_to_string(&cp.edit_distance),
&opt_to_string(&cp.alignment_detail),
```

- [ ] **Step 4: Update make_classified_psm test helper**

Add the 3 new fields to the helper:

```rust
substitution_type: SubstitutionType::None,
edit_distance: None,
alignment_detail: None,
```

Add `use crate::types::SubstitutionType;` at the top of the test module.

- [ ] **Step 5: Run tests**

Run: `cargo test -p protein-copilot-entrapment-analysis output::tests 2>&1 | tail -10`
Expected: All output tests pass

- [ ] **Step 6: Commit**

```bash
git add crates/entrapment-analysis/src/output.rs
git commit -m "feat(entrapment): add substitution_type, edit_distance, alignment_detail columns to TSV output"
```

---

### Task 7: HTML report — new column and mDa display (report.rs + template)

**Files:**
- Modify: `crates/entrapment-analysis/src/report.rs`
- Modify: `crates/entrapment-analysis/templates/entrapment_report.html`

- [ ] **Step 1: Write test for new report column**

Add to `report.rs` tests:

```rust
#[test]
fn test_psm_row_includes_substitution_type() {
    use crate::types::SubstitutionType;
    let mut cp = make_psm("PEPQDEK", PsmGroup::Trap, DiscriminabilityLevel::L2);
    cp.substitution_type = SubstitutionType::QKSubstitution;
    cp.edit_distance = Some(1);
    cp.alignment_detail = Some("Q3→K".to_string());
    let row = PsmRow::from_classified(&cp);
    assert_eq!(row.substitution_type, "QKSubstitution");
    assert_eq!(row.edit_distance, "1");
    assert_eq!(row.alignment_detail, "Q3→K");
}
```

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL — `PsmRow` doesn't have `substitution_type` field

- [ ] **Step 3: Add new fields to PsmRow and from_classified()**

In `PsmRow`:

```rust
struct PsmRow {
    // ... existing fields ...
    substitution_type: String,
    edit_distance: String,
    alignment_detail: String,
}
```

In `PsmRow::from_classified()`:

```rust
substitution_type: cp.substitution_type.to_string(),
edit_distance: cp.edit_distance.map(|d| d.to_string()).unwrap_or_default(),
alignment_detail: cp.alignment_detail.clone().unwrap_or_default(),
```

- [ ] **Step 4: Update make_psm test helper**

Add to the test helper `make_psm()`:

```rust
substitution_type: SubstitutionType::None,
edit_distance: None,
alignment_detail: None,
```

Add `use crate::types::SubstitutionType;` at the top of the test module (if not already imported).

- [ ] **Step 5: Update HTML template**

In `entrapment_report.html`, find the table header row (the `<th>` elements for the PSM table) and add:

```html
<th>Substitution Type</th>
<th>Edit Dist</th>
<th>Alignment</th>
```

In the JavaScript function that builds table rows, add the new columns to match the PsmRow fields:

```javascript
// In the row building function, add cells for:
row.substitution_type
row.edit_distance
row.alignment_detail
```

Add mDa display for delta_mass values < 0.1 Da:

```javascript
// Format delta mass: show mDa for small values
function formatDeltaMass(val) {
    if (!val || val === '') return '';
    const num = parseFloat(val);
    if (isNaN(num)) return val;
    if (Math.abs(num) < 0.1) return (num * 1000).toFixed(1) + ' mDa';
    return num.toFixed(4) + ' Da';
}
```

- [ ] **Step 6: Run tests**

Run: `cargo test -p protein-copilot-entrapment-analysis report::tests 2>&1 | tail -10`
Expected: All report tests pass

- [ ] **Step 7: Commit**

```bash
git add crates/entrapment-analysis/src/report.rs crates/entrapment-analysis/templates/entrapment_report.html
git commit -m "feat(entrapment): add substitution_type column to HTML report, mDa display for small deltas"
```

---

### Task 8: MCP server output schema updates (tools.rs)

**Files:**
- Modify: `crates/mcp-server/src/tools.rs`

- [ ] **Step 1: Update FindSimilarTargetsOutput struct**

Add new fields:

```rust
struct FindSimilarTargetsOutput {
    // ... existing fields ...
    substitution_type: Option<String>,
    edit_distance: Option<u32>,
    alignment_detail: Option<String>,
}
```

- [ ] **Step 2: Update the find_similar_targets handler**

In the handler where `FindSimilarTargetsOutput` is constructed:

```rust
substitution_type: Some(result.substitution_type.to_string()),
edit_distance: result.edit_distance,
alignment_detail: result.alignment_detail,
```

- [ ] **Step 3: Update FindSimilarTargetsInput to accept len_tolerance**

Add optional parameter:

```rust
struct FindSimilarTargetsInput {
    // existing...
    /// Length tolerance for cross-length matching (default: 2)
    #[serde(default)]
    max_mismatches: Option<u16>,
}
```

No new field needed — the existing `max_mismatches` already covers this.

- [ ] **Step 4: Update the tool description to mention edit distance**

Change the description from:
```
"compares the query peptide against all same-length target peptides using Hamming distance"
```
to:
```
"compares the query peptide against target peptides using edit distance (Hamming for same-length, Levenshtein for cross-length). Returns closest matches with mass differences and substitution type annotations."
```

- [ ] **Step 5: Build to verify compilation**

Run: `cargo build -p protein-copilot-mcp-server 2>&1 | tail -5`
Expected: successful build

- [ ] **Step 6: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "feat(entrapment): update MCP tool outputs with substitution_type, edit_distance, alignment_detail"
```

---

### Task 9: Integration tests — cross-length matching and substitution detection

**Files:**
- Create: `crates/entrapment-analysis/tests/v2_edit_distance.rs`

- [ ] **Step 1: Create integration test file**

```rust
//! Integration tests for entrapment v2: edit distance + substitution type.

use std::io::Write;

use protein_copilot_entrapment_analysis::config::EntrapmentConfig;
use protein_copilot_entrapment_analysis::{
    ClassifiedPsm, DiscriminabilityLevel, EntrapmentAnalyzer, PsmGroup, SubstitutionType,
    UnifiedPsm,
};

fn make_psm(peptide: &str, protein: &str) -> UnifiedPsm {
    UnifiedPsm {
        peptide: peptide.to_string(),
        charge: Some(2),
        precursor_mz: None,
        retention_time: None,
        scan_number: None,
        spectrum_file: None,
        protein_ids: protein.to_string(),
        q_value: Some(0.01),
    }
}

fn make_config_yaml() -> String {
    r#"
version: 1
target:
  rules:
    - type: AccessionContains
      any_of: ["_HUMAN"]
trap:
  rules:
    - type: AccessionContains
      any_of: ["_YEAST", "_ECOLI"]
similarity:
  max_mismatches: 2
  delta_mass_threshold_da: 1.0
  len_tolerance: 2
"#
    .to_string()
}

#[test]
fn test_cross_length_indel_gets_l2_not_l4() {
    // Target has "PEPGGDEK" (8 aa), trap PSM has "PEPNDEK" (7 aa)
    // N↔GG isobaric dipeptide, delta_mass ≈ 0
    // v1 would classify as L4 (different length), v2 should classify as L2
    let mut fasta_file = tempfile::NamedTempFile::new().unwrap();
    write!(
        fasta_file,
        ">sp|P001|TEST_HUMAN Test protein\nPEPGGDEKLASTR\n"
    )
    .unwrap();

    let config = EntrapmentConfig::from_yaml_str(&make_config_yaml()).unwrap();
    let analyzer = EntrapmentAnalyzer::new(config, fasta_file.path()).unwrap();

    let psm = make_psm("PEPNDEK", "sp|Q001|TRAP_YEAST");
    let result = analyzer.classify(&psm).unwrap();

    assert_eq!(result.group, PsmGroup::Trap);
    // This is the key v2 assertion: should NOT be L4 anymore
    assert!(
        result.level == DiscriminabilityLevel::L2 || result.level == DiscriminabilityLevel::L3,
        "cross-length indel should be L2 or L3, got {:?}",
        result.level
    );
    assert!(result.edit_distance.is_some());
}

#[test]
fn test_qk_substitution_detected() {
    // Target: "PEPTKDEK", trap: "PEPTQDEK" — Q↔K substitution
    let mut fasta_file = tempfile::NamedTempFile::new().unwrap();
    write!(
        fasta_file,
        ">sp|P001|TEST_HUMAN Test\nPEPTKDEKLASTR\n"
    )
    .unwrap();

    let config = EntrapmentConfig::from_yaml_str(&make_config_yaml()).unwrap();
    let analyzer = EntrapmentAnalyzer::new(config, fasta_file.path()).unwrap();

    let psm = make_psm("PEPTQDEK", "sp|Q001|TRAP_YEAST");
    let result = analyzer.classify(&psm).unwrap();

    assert_eq!(result.level, DiscriminabilityLevel::L2);
    assert_eq!(result.substitution_type, SubstitutionType::QKSubstitution);
}

#[test]
fn test_backward_compatible_same_length_matching() {
    // Same-length D→N substitution should still work exactly as v1
    let mut fasta_file = tempfile::NamedTempFile::new().unwrap();
    write!(
        fasta_file,
        ">sp|P001|TEST_HUMAN Test\nNGFLLDGFPR\n"
    )
    .unwrap();

    let config = EntrapmentConfig::from_yaml_str(&make_config_yaml()).unwrap();
    let analyzer = EntrapmentAnalyzer::new(config, fasta_file.path()).unwrap();

    let psm = make_psm("DGFLLDGFPR", "sp|Q001|TRAP_YEAST");
    let result = analyzer.classify(&psm).unwrap();

    assert_eq!(result.level, DiscriminabilityLevel::L2);
    assert_eq!(result.mismatches, Some(1));
    assert!(result.delta_mass_da.unwrap().abs() < 1.0);
}

#[test]
fn test_true_trap_still_l4() {
    // Completely unrelated peptide should still be L4
    let mut fasta_file = tempfile::NamedTempFile::new().unwrap();
    write!(
        fasta_file,
        ">sp|P001|TEST_HUMAN Test\nAAAAAAAAR\n"
    )
    .unwrap();

    let config = EntrapmentConfig::from_yaml_str(&make_config_yaml()).unwrap();
    let analyzer = EntrapmentAnalyzer::new(config, fasta_file.path()).unwrap();

    let psm = make_psm("WWWWWWWWR", "sp|Q001|TRAP_YEAST");
    let result = analyzer.classify(&psm).unwrap();

    assert_eq!(result.level, DiscriminabilityLevel::L4);
    assert_eq!(result.substitution_type, SubstitutionType::None);
}
```

- [ ] **Step 2: Run integration tests**

Run: `cargo test -p protein-copilot-entrapment-analysis --test v2_edit_distance 2>&1 | tail -10`
Expected: All 4 integration tests pass

- [ ] **Step 3: Commit**

```bash
git add crates/entrapment-analysis/tests/v2_edit_distance.rs
git commit -m "test(entrapment): add v2 integration tests for cross-length matching and substitution detection"
```

---

### Task 10: Full regression test and cleanup

**Files:**
- All files from Tasks 1-9

- [ ] **Step 1: Run full test suite**

Run: `cargo test -p protein-copilot-entrapment-analysis 2>&1 | tail -10`
Expected: All tests pass (56 original + ~15 new unit + 4 integration)

- [ ] **Step 2: Run workspace build**

Run: `cargo build --release 2>&1 | tail -5`
Expected: No errors, no warnings

- [ ] **Step 3: Run clippy**

Run: `cargo clippy -p protein-copilot-entrapment-analysis -- -D warnings 2>&1 | tail -10`
Expected: No warnings

- [ ] **Step 4: Run full workspace tests**

Run: `cargo test 2>&1 | tail -10`
Expected: All workspace tests pass

- [ ] **Step 5: Final commit if any clippy/formatting fixes needed**

```bash
cargo fmt -p protein-copilot-entrapment-analysis
git add -A
git commit -m "chore(entrapment): clippy and formatting fixes for v2"
```

---

## File Change Summary

| File | Action | Description |
|------|--------|-------------|
| `crates/entrapment-analysis/src/types.rs` | Modify | Add `SubstitutionType` enum, 3 new fields on `ClassifiedPsm` |
| `crates/entrapment-analysis/src/levenshtein.rs` | **Create** | Wagner-Fischer algorithm + alignment backtracking + delta_mass |
| `crates/entrapment-analysis/src/config.rs` | Modify | Rename `delta_mz_threshold_da` → `delta_mass_threshold_da` with alias, add 3 v2 fields |
| `crates/entrapment-analysis/src/digest.rs` | Modify | Add k-mer index, `SimilarityMatch`, `find_similar()` method |
| `crates/entrapment-analysis/src/similarity.rs` | Modify | Upgrade `classify_single()` with edit distance + categorize_substitution |
| `crates/entrapment-analysis/src/output.rs` | Modify | Add 3 new TSV columns |
| `crates/entrapment-analysis/src/report.rs` | Modify | Add new PsmRow fields, mDa display |
| `crates/entrapment-analysis/templates/entrapment_report.html` | Modify | New table columns, mDa formatting |
| `crates/entrapment-analysis/src/lib.rs` | Modify | Add `pub mod levenshtein;`, export `SubstitutionType` |
| `crates/mcp-server/src/tools.rs` | Modify | Update output structs + tool description |
| `crates/entrapment-analysis/tests/v2_edit_distance.rs` | **Create** | 4 integration tests |
