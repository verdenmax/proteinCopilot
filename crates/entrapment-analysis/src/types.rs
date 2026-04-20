//! Core type definitions for entrapment analysis.
//!
//! Defines discriminability levels (L0–L4), PSM grouping, and result summary types
//! used throughout the entrapment classification pipeline.

use std::fmt;

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Enums
// ---------------------------------------------------------------------------

/// Discriminability level assigned to each trap-database PSM hit.
///
/// Levels range from L0 (razor attribution error – exact match in the target
/// database) to L4 (true trap hit with no close homolog).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DiscriminabilityLevel {
    /// Exact sequence match in the target database (razor error).
    L0,
    /// L/I (leucine / isoleucine) isomer match only.
    L1,
    /// One amino-acid mismatch with `delta_mass` below the near-isobaric threshold.
    L2,
    /// One or two mismatches but *not* near-isobaric (distinguishable homolog).
    L3,
    /// No close match in the target database (true trap hit).
    L4,
}

impl DiscriminabilityLevel {
    /// Returns the short string label for this level.
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::L0 => "L0",
            Self::L1 => "L1",
            Self::L2 => "L2",
            Self::L3 => "L3",
            Self::L4 => "L4",
        }
    }
}

impl fmt::Display for DiscriminabilityLevel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// High-level group a PSM is classified into.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PsmGroup {
    /// PSM matched a target-database protein.
    Target,
    /// PSM matched a trap-database protein.
    Trap,
    /// PSM cannot be unambiguously assigned to target or trap.
    Ambiguous,
}

impl fmt::Display for PsmGroup {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            Self::Target => "target",
            Self::Trap => "trap",
            Self::Ambiguous => "ambiguous",
        };
        f.write_str(s)
    }
}

// ---------------------------------------------------------------------------
// Structs
// ---------------------------------------------------------------------------

/// A unified PSM record normalised from any search-engine result format.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UnifiedPsm {
    /// Stripped amino-acid sequence (no modifications / flanking residues).
    pub peptide: String,
    /// Charge state of the precursor ion.
    pub charge: Option<i32>,
    /// Observed precursor *m/z*.
    pub precursor_mz: Option<f64>,
    /// Retention time in **minutes**.
    pub retention_time: Option<f64>,
    /// Scan number (1-based).
    pub scan_number: Option<u32>,
    /// Name of the spectrum / raw file.
    pub spectrum_file: Option<String>,
    /// Semicolon-separated protein accessions.
    pub protein_ids: String,
    /// False-discovery-rate q-value.
    pub q_value: Option<f64>,
}

/// A PSM that has been classified with group and discriminability information.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClassifiedPsm {
    /// The underlying unified PSM.
    pub psm: UnifiedPsm,
    /// High-level group assignment.
    pub group: PsmGroup,
    /// Discriminability level (L0–L4).
    pub level: DiscriminabilityLevel,
    /// Closest matching target-database peptide (if any).
    pub best_target_peptide: Option<String>,
    /// Protein accession of the best target match.
    pub best_target_protein: Option<String>,
    /// Hamming distance (number of mismatched positions) to the best target peptide.
    pub mismatches: Option<u8>,
    /// Mass difference to the best target peptide in Da.
    pub delta_mass_da: Option<f64>,
    /// Human-readable diff positions, e.g. `"[2:D->N,5:G->A]"`.
    pub diff_positions: Option<String>,
}

/// Per-level hit counts for the five discriminability levels.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct LevelCounts {
    /// Count of L0 (razor error) PSMs.
    pub l0: usize,
    /// Count of L1 (L/I isomer) PSMs.
    pub l1: usize,
    /// Count of L2 (near-isobaric) PSMs.
    pub l2: usize,
    /// Count of L3 (distinguishable homolog) PSMs.
    pub l3: usize,
    /// Count of L4 (true trap) PSMs.
    pub l4: usize,
}

impl LevelCounts {
    /// Returns the total number of classified PSMs across all levels.
    pub fn total(&self) -> usize {
        self.l0 + self.l1 + self.l2 + self.l3 + self.l4
    }

    /// Increments the counter for the given discriminability level by one.
    pub fn increment(&mut self, level: DiscriminabilityLevel) {
        match level {
            DiscriminabilityLevel::L0 => self.l0 += 1,
            DiscriminabilityLevel::L1 => self.l1 += 1,
            DiscriminabilityLevel::L2 => self.l2 += 1,
            DiscriminabilityLevel::L3 => self.l3 += 1,
            DiscriminabilityLevel::L4 => self.l4 += 1,
        }
    }
}

/// Summary statistics for an entrapment analysis run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EntrapmentSummary {
    /// Total number of PSMs processed.
    pub total_psms: usize,
    /// Number of PSMs assigned to the target group.
    pub target_psms: usize,
    /// Number of PSMs assigned to the trap group.
    pub trap_psms: usize,
    /// Number of ambiguous PSMs.
    pub ambiguous_psms: usize,
    /// Per-level hit counts.
    pub level_counts: LevelCounts,
    /// Top protein families contributing L0 (razor-error) hits.
    pub top_razor_families: Vec<RazorFamily>,
}

/// A protein family that contributes L0 (razor attribution error) hits.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RazorFamily {
    /// Protein family name.
    pub family: String,
    /// Number of L0 PSMs from this family.
    pub count: usize,
    /// An example peptide sequence from this family.
    pub example_peptide: String,
    /// An example trap-database protein accession.
    pub example_trap_protein: String,
    /// An example target-database protein accession.
    pub example_target_protein: String,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_level_display() {
        assert_eq!(DiscriminabilityLevel::L0.to_string(), "L0");
        assert_eq!(DiscriminabilityLevel::L1.to_string(), "L1");
        assert_eq!(DiscriminabilityLevel::L2.to_string(), "L2");
        assert_eq!(DiscriminabilityLevel::L3.to_string(), "L3");
        assert_eq!(DiscriminabilityLevel::L4.to_string(), "L4");
    }

    #[test]
    fn test_level_counts() {
        let mut counts = LevelCounts::default();
        assert_eq!(counts.total(), 0);

        counts.increment(DiscriminabilityLevel::L0);
        counts.increment(DiscriminabilityLevel::L0);
        counts.increment(DiscriminabilityLevel::L1);
        counts.increment(DiscriminabilityLevel::L2);
        counts.increment(DiscriminabilityLevel::L3);
        counts.increment(DiscriminabilityLevel::L3);
        counts.increment(DiscriminabilityLevel::L3);
        counts.increment(DiscriminabilityLevel::L4);

        assert_eq!(counts.l0, 2);
        assert_eq!(counts.l1, 1);
        assert_eq!(counts.l2, 1);
        assert_eq!(counts.l3, 3);
        assert_eq!(counts.l4, 1);
        assert_eq!(counts.total(), 8);
    }

    #[test]
    fn test_psm_group_display() {
        assert_eq!(PsmGroup::Target.to_string(), "target");
        assert_eq!(PsmGroup::Trap.to_string(), "trap");
        assert_eq!(PsmGroup::Ambiguous.to_string(), "ambiguous");
    }
}
