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
// SubstitutionType (v2)
// ---------------------------------------------------------------------------

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
    /// Retention time in **minutes** (apex / single value).
    pub retention_time: Option<f64>,
    /// Elution window start in **minutes** (e.g. DIA-NN `RT.Start`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rt_start: Option<f64>,
    /// Elution window end in **minutes** (e.g. DIA-NN `RT.Stop`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rt_stop: Option<f64>,
    /// Scan number (1-based).
    pub scan_number: Option<u32>,
    /// Name of the spectrum / raw file.
    pub spectrum_file: Option<String>,
    /// Semicolon-separated protein accessions.
    pub protein_ids: String,
    /// False-discovery-rate q-value.
    pub q_value: Option<f64>,
    /// Parsed modifications: (0-based position, delta_mass_da).
    #[serde(default)]
    pub modifications: Vec<(usize, f64)>,
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
    pub mismatches: Option<u16>,
    /// Mass difference to the best target peptide in Da.
    pub delta_mass_da: Option<f64>,
    /// Human-readable diff positions, e.g. `"[2:D->N,5:G->A]"`.
    pub diff_positions: Option<String>,
    /// Substitution type annotation (v2). Informational only.
    pub substitution_type: SubstitutionType,
    /// Edit distance to best target (v2). Equals Hamming distance for same-length matches.
    pub edit_distance: Option<u32>,
    /// Alignment detail string (v2), e.g. "D0→N" or "ins:G@5".
    pub alignment_detail: Option<String>,
    /// Fragment ion provenance analysis result (v3). None if not traced.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub provenance: Option<crate::provenance::FragmentProvenance>,
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
// v4 Multi-Target Provenance Types
// ---------------------------------------------------------------------------

/// Whether a co-eluting candidate is a light or heavy (SILAC) form.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum LabelForm {
    /// Light (unlabeled) form.
    Light,
    /// Heavy (SILAC-labeled) form with shifted precursor and residue deltas.
    Heavy {
        /// Heavy precursor m/z.
        precursor_mz_heavy: f64,
        /// (0-based position, delta_Da) for each labeled residue (K or R).
        residue_deltas: Vec<(usize, f64)>,
    },
}

/// A co-eluting target peptide candidate found within the same DIA window and RT range.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoElutingCandidate {
    /// Target peptide sequence.
    pub peptide: String,
    /// Protein accession(s).
    pub protein_ids: Vec<String>,
    /// Precursor m/z (light form; for Heavy, use `label_form` field).
    pub precursor_mz: f64,
    /// Charge state.
    pub charge: i32,
    /// Elution window start (minutes).
    pub rt_start: f64,
    /// Elution window stop (minutes).
    pub rt_stop: f64,
    /// Light or Heavy form.
    pub label_form: LabelForm,
    /// Modifications: (0-based position, delta_mass_Da).
    pub modifications: Vec<(usize, f64)>,
}

/// A single target ion match for a multi-target annotated peak.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TargetIonMatch {
    /// Index into `MultiTargetProvenance::candidates`.
    pub candidate_index: usize,
    /// Ion label, e.g. "b3+1", "y5+2".
    pub ion_label: String,
    /// Matching error in ppm.
    pub delta_ppm: f64,
}

/// An observed peak annotated with multi-target provenance information.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiAnnotatedPeak {
    /// Observed m/z.
    pub mz_observed: f64,
    /// Observed intensity.
    pub intensity: f64,
    /// Ion label from the trap peptide (if matched).
    pub trap_ion: Option<String>,
    /// All target ion matches (may be 0, 1, or many).
    pub target_matches: Vec<TargetIonMatch>,
}

/// Complete multi-target provenance result for one trap PSM.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiTargetProvenance {
    /// The trap PSM being analyzed.
    pub trap_peptide: String,
    /// Scan number of the analyzed spectrum.
    pub scan_number: u32,
    /// All co-eluting target candidates (light + heavy).
    pub candidates: Vec<CoElutingCandidate>,
    /// Per-peak multi-target annotation.
    pub annotated_peaks: Vec<MultiAnnotatedPeak>,
    /// Summary: count of peaks matching only trap ions.
    pub trap_only_count: u32,
    /// Summary: count of peaks matching at least one target (not trap).
    pub target_only_count: u32,
    /// Summary: count of peaks matching both trap and at least one target.
    pub shared_count: u32,
    /// Summary: count of peaks matching nothing.
    pub unassigned_count: u32,
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

    #[test]
    fn test_substitution_type_serde() {
        let st = SubstitutionType::QKSubstitution;
        let json = serde_json::to_string(&st).unwrap();
        assert_eq!(json, r#""QKSubstitution""#);
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
        assert_eq!(
            SubstitutionType::Distinguishable.as_str(),
            "Distinguishable"
        );
        let idb = SubstitutionType::IsobaricDipeptide {
            single_residue: 'N',
            dipeptide: "GG".to_string(),
        };
        assert_eq!(idb.as_str(), "IsobaricDipeptide");
    }

    #[test]
    fn classified_psm_provenance_default() {
        // Deserialize JSON without provenance field → None
        let json = r#"{"psm":{"peptide":"PEP","charge":2,"precursor_mz":300.0,"retention_time":null,"scan_number":null,"spectrum_file":null,"protein_ids":"P1","q_value":0.01},"group":"Target","level":"L4","best_target_peptide":null,"best_target_protein":null,"mismatches":null,"delta_mass_da":null,"diff_positions":null,"substitution_type":"None","edit_distance":null,"alignment_detail":null}"#;
        let cpsm: ClassifiedPsm = serde_json::from_str(json).unwrap();
        assert!(cpsm.provenance.is_none());
    }

    #[test]
    fn unified_psm_modifications_default() {
        let json = r#"{"peptide":"PEP","charge":2,"precursor_mz":300.0,"retention_time":null,"scan_number":null,"spectrum_file":null,"protein_ids":"P1","q_value":0.01}"#;
        let psm: UnifiedPsm = serde_json::from_str(json).unwrap();
        assert!(psm.modifications.is_empty());
    }

    #[test]
    fn unified_psm_modifications_roundtrip() {
        let psm = UnifiedPsm {
            peptide: "ACDFK".into(),
            charge: Some(2),
            precursor_mz: Some(300.0),
            retention_time: None,
            rt_start: None,
            rt_stop: None,
            scan_number: None,
            spectrum_file: None,
            protein_ids: "P1".into(),
            q_value: Some(0.01),
            modifications: vec![(1, 57.021464)],
        };
        let json = serde_json::to_string(&psm).unwrap();
        let deser: UnifiedPsm = serde_json::from_str(&json).unwrap();
        assert_eq!(deser.modifications.len(), 1);
        assert_eq!(deser.modifications[0].0, 1);
        assert!((deser.modifications[0].1 - 57.021464).abs() < 1e-6);
    }

    #[test]
    fn test_label_form_light() {
        let form = LabelForm::Light;
        assert!(matches!(form, LabelForm::Light));
    }

    #[test]
    fn test_label_form_heavy() {
        let form = LabelForm::Heavy {
            precursor_mz_heavy: 556.13,
            residue_deltas: vec![(9, 8.014199)],
        };
        if let LabelForm::Heavy { precursor_mz_heavy, residue_deltas } = form {
            assert!((precursor_mz_heavy - 556.13).abs() < 0.01);
            assert_eq!(residue_deltas.len(), 1);
        } else {
            panic!("expected Heavy");
        }
    }

    #[test]
    fn test_co_eluting_candidate() {
        let candidate = CoElutingCandidate {
            peptide: "STTSGHLVYK".to_string(),
            protein_ids: vec!["sp|P12345|EF1A_HUMAN".to_string()],
            precursor_mz: 548.12,
            charge: 2,
            rt_start: 34.5,
            rt_stop: 35.8,
            label_form: LabelForm::Light,
            modifications: vec![],
        };
        assert_eq!(candidate.peptide, "STTSGHLVYK");
        assert_eq!(candidate.charge, 2);
    }

    #[test]
    fn test_multi_annotated_peak() {
        let peak = MultiAnnotatedPeak {
            mz_observed: 285.155,
            intensity: 45230.0,
            trap_ion: Some("b3+1".to_string()),
            target_matches: vec![
                TargetIonMatch {
                    candidate_index: 0,
                    ion_label: "b3+1".to_string(),
                    delta_ppm: -2.1,
                },
            ],
        };
        assert!(peak.trap_ion.is_some());
        assert_eq!(peak.target_matches.len(), 1);
    }
}
