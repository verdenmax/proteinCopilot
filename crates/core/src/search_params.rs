//! Search parameter definitions for mass spectrometry database search.
//!
//! This module defines the types needed to configure a proteomics search:
//! - [`Enzyme`] — Digestion enzyme
//! - [`ModPosition`] — Modification position on the peptide
//! - [`Modification`] — Chemical modification (fixed or variable)
//! - [`ToleranceUnit`] / [`MassTolerance`] — Mass accuracy specification
//! - [`DecoyStrategy`] — Target-decoy approach
//! - [`SearchParams`] — Complete search configuration with validation

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::spectrum::AcquisitionMode;

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Validation errors for search parameters.
#[derive(Debug, Error)]
pub enum SearchParamsError {
    /// Database path is empty or whitespace-only.
    #[error("database_path must not be empty")]
    EmptyDatabasePath,

    /// Precursor mass tolerance is not a finite positive value.
    #[error("precursor tolerance must be a finite positive value, got {value}")]
    InvalidPrecursorTolerance {
        /// The invalid tolerance value.
        value: f64,
    },

    /// Fragment mass tolerance is not a finite positive value.
    #[error("fragment tolerance must be a finite positive value, got {value}")]
    InvalidFragmentTolerance {
        /// The invalid tolerance value.
        value: f64,
    },

    /// Too many missed cleavages.
    #[error("missed_cleavages must be <= {max}, got {actual}")]
    TooManyMissedCleavages {
        /// The value provided.
        actual: u32,
        /// The maximum allowed.
        max: u32,
    },

    /// A modification has a non-finite mass_delta (NaN or Infinity).
    #[error("modification \"{name}\" has non-finite mass_delta: {value}")]
    InvalidModificationMassDelta {
        /// Name of the modification.
        name: String,
        /// The invalid mass_delta value.
        value: f64,
    },

    /// max_variable_modifications exceeds limit.
    #[error("max_variable_modifications must be <= {max}, got {actual}")]
    TooManyVariableMods {
        /// The value provided.
        actual: u32,
        /// The maximum allowed.
        max: u32,
    },

    /// min_peptide_length is zero.
    #[error("min_peptide_length must be >= 1, got 0")]
    ZeroPeptideLength,

    /// min_peptide_length > max_peptide_length.
    #[error("min_peptide_length ({min}) must be <= max_peptide_length ({max})")]
    InvalidPeptideLengthRange {
        /// The minimum value provided.
        min: u32,
        /// The maximum value provided.
        max: u32,
    },
}

// ---------------------------------------------------------------------------
// Enzyme
// ---------------------------------------------------------------------------

/// Digestion enzyme used for protein cleavage.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub enum Enzyme {
    /// Trypsin — cleaves after K/R (not before P).
    Trypsin,
    /// Lys-C — cleaves after K.
    LysC,
    /// Glu-C — cleaves after D/E.
    GluC,
    /// Asp-N — cleaves before D.
    AspN,
    /// Chymotrypsin — cleaves after F/W/Y/L.
    Chymotrypsin,
    /// Trypsin/P — cleaves after K/R (including before P).
    TrypsinP,
    /// No specific cleavage rule.
    NonSpecific,
    /// User-defined enzyme with custom cleavage rule.
    Custom {
        /// Enzyme name.
        name: String,
        /// Cleavage rule (regex or engine-specific syntax).
        cleavage_rule: String,
    },
}

impl std::fmt::Display for Enzyme {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Enzyme::Trypsin => write!(f, "Trypsin"),
            Enzyme::LysC => write!(f, "Lys-C"),
            Enzyme::GluC => write!(f, "Glu-C"),
            Enzyme::AspN => write!(f, "Asp-N"),
            Enzyme::Chymotrypsin => write!(f, "Chymotrypsin"),
            Enzyme::TrypsinP => write!(f, "Trypsin/P"),
            Enzyme::NonSpecific => write!(f, "NonSpecific"),
            Enzyme::Custom { name, .. } => write!(f, "Custom({name})"),
        }
    }
}

// ---------------------------------------------------------------------------
// ModPosition
// ---------------------------------------------------------------------------

/// Position where a modification can occur.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, JsonSchema)]
pub enum ModPosition {
    /// Anywhere on the peptide.
    Anywhere,
    /// Any peptide N-terminus.
    AnyNTerm,
    /// Any peptide C-terminus.
    AnyCTerm,
    /// Protein N-terminus only.
    ProteinNTerm,
    /// Protein C-terminus only.
    ProteinCTerm,
}

// ---------------------------------------------------------------------------
// Modification
// ---------------------------------------------------------------------------

/// A chemical modification (fixed or variable) applied during search.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct Modification {
    /// Modification name (e.g. "Carbamidomethyl", "Oxidation").
    pub name: String,
    /// Mass shift in Daltons.
    /// Positive for mass-increasing modifications (e.g., +57.021 for Carbamidomethyl),
    /// negative for mass-decreasing (rare; e.g., -18.011 for dehydration).
    pub mass_delta: f64,
    /// Target residues (e.g. `['C']` for Carbamidomethyl, `['M']` for Oxidation).
    pub residues: Vec<char>,
    /// Where on the peptide/protein this modification can occur.
    pub position: ModPosition,
}

// ---------------------------------------------------------------------------
// ToleranceUnit & MassTolerance
// ---------------------------------------------------------------------------

/// Unit for mass tolerance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, JsonSchema)]
pub enum ToleranceUnit {
    /// Parts per million.
    Ppm,
    /// Daltons (absolute mass).
    Da,
}

impl std::fmt::Display for ToleranceUnit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ToleranceUnit::Ppm => write!(f, "ppm"),
            ToleranceUnit::Da => write!(f, "Da"),
        }
    }
}

/// Mass tolerance specification (value + unit).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct MassTolerance {
    /// Tolerance value (must be positive).
    pub value: f64,
    /// Tolerance unit.
    pub unit: ToleranceUnit,
}

impl std::fmt::Display for MassTolerance {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} {}", self.value, self.unit)
    }
}

// ---------------------------------------------------------------------------
// DecoyStrategy
// ---------------------------------------------------------------------------

/// Target-decoy strategy for FDR estimation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, JsonSchema)]
pub enum DecoyStrategy {
    /// Reverse protein sequences.
    Reverse,
    /// Shuffle protein sequences.
    Shuffle,
    /// No decoy database.
    None,
}

// ---------------------------------------------------------------------------
// SearchParams
// ---------------------------------------------------------------------------

/// Maximum allowed value for `missed_cleavages`.
const MAX_MISSED_CLEAVAGES: u32 = 5;

/// Maximum allowed value for `max_variable_modifications`.
const MAX_VARIABLE_MODS_LIMIT: u32 = 10;

fn default_max_variable_modifications() -> u32 {
    3
}
fn default_min_peptide_length() -> u32 {
    7
}
fn default_max_peptide_length() -> u32 {
    50
}

/// Complete search configuration for a proteomics database search.
///
/// Use [`SearchParams::validate()`] after construction or deserialization
/// to ensure all values are within acceptable ranges.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct SearchParams {
    /// Digestion enzyme.
    pub enzyme: Enzyme,
    /// Maximum number of missed cleavage sites (0–5).
    pub missed_cleavages: u32,
    /// Modifications always present (e.g. Carbamidomethyl on C).
    pub fixed_modifications: Vec<Modification>,
    /// Modifications that may or may not be present (e.g. Oxidation on M).
    pub variable_modifications: Vec<Modification>,
    /// Precursor ion mass tolerance.
    pub precursor_tolerance: MassTolerance,
    /// Fragment ion mass tolerance.
    pub fragment_tolerance: MassTolerance,
    /// Path to the FASTA protein database file.
    pub database_path: String,
    /// Target-decoy strategy for FDR estimation.
    pub decoy_strategy: DecoyStrategy,
    /// Data acquisition mode. `None` = auto-detect or not applicable.
    #[serde(default)]
    pub acquisition_mode: Option<AcquisitionMode>,

    /// Maximum number of variable modifications per peptide (default: 3).
    /// Limits combinatorial explosion during variable modification enumeration.
    #[serde(default = "default_max_variable_modifications")]
    pub max_variable_modifications: u32,

    /// Minimum peptide length in residues (default: 7).
    /// Peptides shorter than this are excluded from search results.
    #[serde(default = "default_min_peptide_length")]
    pub min_peptide_length: u32,

    /// Maximum peptide length in residues (default: 50).
    /// Peptides longer than this are excluded from search results.
    #[serde(default = "default_max_peptide_length")]
    pub max_peptide_length: u32,

    /// Search engine to use. `None` means default engine (SimpleSearch).
    #[serde(default)]
    pub engine: Option<String>,
}

impl SearchParams {
    /// Validates that all parameter values are within acceptable ranges.
    ///
    /// Checks:
    /// - `database_path` is not empty
    /// - `precursor_tolerance.value` is finite and > 0
    /// - `fragment_tolerance.value` is finite and > 0
    /// - `missed_cleavages` ≤ 5
    /// - All modification `mass_delta` values are finite
    pub fn validate(&self) -> Result<(), SearchParamsError> {
        if self.database_path.trim().is_empty() {
            return Err(SearchParamsError::EmptyDatabasePath);
        }
        if !self.precursor_tolerance.value.is_finite() || self.precursor_tolerance.value <= 0.0 {
            return Err(SearchParamsError::InvalidPrecursorTolerance {
                value: self.precursor_tolerance.value,
            });
        }
        if !self.fragment_tolerance.value.is_finite() || self.fragment_tolerance.value <= 0.0 {
            return Err(SearchParamsError::InvalidFragmentTolerance {
                value: self.fragment_tolerance.value,
            });
        }
        if self.missed_cleavages > MAX_MISSED_CLEAVAGES {
            return Err(SearchParamsError::TooManyMissedCleavages {
                actual: self.missed_cleavages,
                max: MAX_MISSED_CLEAVAGES,
            });
        }
        for m in self
            .fixed_modifications
            .iter()
            .chain(&self.variable_modifications)
        {
            if !m.mass_delta.is_finite() {
                return Err(SearchParamsError::InvalidModificationMassDelta {
                    name: m.name.clone(),
                    value: m.mass_delta,
                });
            }
        }
        if self.max_variable_modifications > MAX_VARIABLE_MODS_LIMIT {
            return Err(SearchParamsError::TooManyVariableMods {
                actual: self.max_variable_modifications,
                max: MAX_VARIABLE_MODS_LIMIT,
            });
        }
        if self.min_peptide_length == 0 {
            return Err(SearchParamsError::ZeroPeptideLength);
        }
        if self.min_peptide_length > self.max_peptide_length {
            return Err(SearchParamsError::InvalidPeptideLengthRange {
                min: self.min_peptide_length,
                max: self.max_peptide_length,
            });
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn carbamidomethyl() -> Modification {
        Modification {
            name: "Carbamidomethyl".to_string(),
            mass_delta: 57.021464,
            residues: vec!['C'],
            position: ModPosition::Anywhere,
        }
    }

    fn oxidation() -> Modification {
        Modification {
            name: "Oxidation".to_string(),
            mass_delta: 15.994915,
            residues: vec!['M'],
            position: ModPosition::Anywhere,
        }
    }

    fn valid_params() -> SearchParams {
        SearchParams {
            enzyme: Enzyme::Trypsin,
            missed_cleavages: 2,
            fixed_modifications: vec![carbamidomethyl()],
            variable_modifications: vec![oxidation()],
            precursor_tolerance: MassTolerance {
                value: 20.0,
                unit: ToleranceUnit::Ppm,
            },
            fragment_tolerance: MassTolerance {
                value: 0.02,
                unit: ToleranceUnit::Da,
            },
            database_path: "/data/uniprot_human.fasta".to_string(),
            decoy_strategy: DecoyStrategy::Reverse,
            acquisition_mode: None,
            max_variable_modifications: 3,
            min_peptide_length: 7,
            max_peptide_length: 50,
            engine: None,
        }
    }

    // -- Enzyme ---------------------------------------------------------

    #[test]
    fn enzyme_serde_roundtrip() {
        let enzymes = vec![
            Enzyme::Trypsin,
            Enzyme::LysC,
            Enzyme::NonSpecific,
            Enzyme::Custom {
                name: "MyEnzyme".to_string(),
                cleavage_rule: "KR".to_string(),
            },
        ];
        for enzyme in &enzymes {
            let json = serde_json::to_string(enzyme).unwrap();
            let back: Enzyme = serde_json::from_str(&json).unwrap();
            assert_eq!(*enzyme, back);
        }
    }

    #[test]
    fn enzyme_display() {
        assert_eq!(Enzyme::Trypsin.to_string(), "Trypsin");
        assert_eq!(Enzyme::LysC.to_string(), "Lys-C");
        assert_eq!(Enzyme::TrypsinP.to_string(), "Trypsin/P");
        assert_eq!(
            Enzyme::Custom {
                name: "X".to_string(),
                cleavage_rule: "R".to_string()
            }
            .to_string(),
            "Custom(X)"
        );
    }

    // -- Modification ---------------------------------------------------

    #[test]
    fn modification_serde_roundtrip() {
        let m = carbamidomethyl();
        let json = serde_json::to_string_pretty(&m).unwrap();
        let back: Modification = serde_json::from_str(&json).unwrap();
        assert_eq!(m, back);
    }

    #[test]
    fn mod_position_all_variants_roundtrip() {
        let positions = [
            ModPosition::Anywhere,
            ModPosition::AnyNTerm,
            ModPosition::AnyCTerm,
            ModPosition::ProteinNTerm,
            ModPosition::ProteinCTerm,
        ];
        for pos in &positions {
            let json = serde_json::to_string(pos).unwrap();
            let back: ModPosition = serde_json::from_str(&json).unwrap();
            assert_eq!(*pos, back);
        }
    }

    // -- MassTolerance --------------------------------------------------

    #[test]
    fn mass_tolerance_serde_roundtrip() {
        let t = MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        };
        let json = serde_json::to_string(&t).unwrap();
        let back: MassTolerance = serde_json::from_str(&json).unwrap();
        assert_eq!(t, back);
    }

    #[test]
    fn mass_tolerance_display() {
        let t = MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        };
        assert_eq!(t.to_string(), "20 ppm");

        let t2 = MassTolerance {
            value: 0.5,
            unit: ToleranceUnit::Da,
        };
        assert_eq!(t2.to_string(), "0.5 Da");
    }

    // -- DecoyStrategy --------------------------------------------------

    #[test]
    fn decoy_strategy_serde_roundtrip() {
        for ds in [
            DecoyStrategy::Reverse,
            DecoyStrategy::Shuffle,
            DecoyStrategy::None,
        ] {
            let json = serde_json::to_string(&ds).unwrap();
            let back: DecoyStrategy = serde_json::from_str(&json).unwrap();
            assert_eq!(ds, back);
        }
    }

    // -- SearchParams ---------------------------------------------------

    #[test]
    fn search_params_serde_roundtrip() {
        let params = valid_params();
        let json = serde_json::to_string_pretty(&params).unwrap();
        let back: SearchParams = serde_json::from_str(&json).unwrap();
        assert_eq!(params.enzyme, back.enzyme);
        assert_eq!(params.missed_cleavages, back.missed_cleavages);
        assert_eq!(
            params.fixed_modifications.len(),
            back.fixed_modifications.len()
        );
        assert_eq!(params.database_path, back.database_path);
        assert_eq!(params.decoy_strategy, back.decoy_strategy);
    }

    // -- Validation -----------------------------------------------------

    #[test]
    fn validate_passes_for_valid_params() {
        assert!(valid_params().validate().is_ok());
    }

    #[test]
    fn validate_rejects_empty_database_path() {
        let mut p = valid_params();
        p.database_path = "".to_string();
        assert!(matches!(
            p.validate(),
            Err(SearchParamsError::EmptyDatabasePath)
        ));
    }

    #[test]
    fn validate_rejects_whitespace_database_path() {
        let mut p = valid_params();
        p.database_path = "   ".to_string();
        assert!(matches!(
            p.validate(),
            Err(SearchParamsError::EmptyDatabasePath)
        ));
    }

    #[test]
    fn validate_rejects_zero_precursor_tolerance() {
        let mut p = valid_params();
        p.precursor_tolerance.value = 0.0;
        assert!(matches!(
            p.validate(),
            Err(SearchParamsError::InvalidPrecursorTolerance { .. })
        ));
    }

    #[test]
    fn validate_rejects_negative_fragment_tolerance() {
        let mut p = valid_params();
        p.fragment_tolerance.value = -1.0;
        assert!(matches!(
            p.validate(),
            Err(SearchParamsError::InvalidFragmentTolerance { .. })
        ));
    }

    #[test]
    fn validate_rejects_excessive_missed_cleavages() {
        let mut p = valid_params();
        p.missed_cleavages = 6;
        let err = p.validate().unwrap_err();
        assert!(matches!(
            err,
            SearchParamsError::TooManyMissedCleavages { actual: 6, max: 5 }
        ));
    }

    #[test]
    fn validate_allows_max_missed_cleavages() {
        let mut p = valid_params();
        p.missed_cleavages = 5;
        assert!(p.validate().is_ok());
    }

    // -- NaN/Infinity validation ----------------------------------------

    #[test]
    fn validate_rejects_nan_precursor_tolerance() {
        let mut p = valid_params();
        p.precursor_tolerance.value = f64::NAN;
        assert!(matches!(
            p.validate(),
            Err(SearchParamsError::InvalidPrecursorTolerance { .. })
        ));
    }

    #[test]
    fn validate_rejects_infinity_precursor_tolerance() {
        let mut p = valid_params();
        p.precursor_tolerance.value = f64::INFINITY;
        assert!(matches!(
            p.validate(),
            Err(SearchParamsError::InvalidPrecursorTolerance { .. })
        ));
    }

    #[test]
    fn validate_rejects_nan_fragment_tolerance() {
        let mut p = valid_params();
        p.fragment_tolerance.value = f64::NAN;
        assert!(matches!(
            p.validate(),
            Err(SearchParamsError::InvalidFragmentTolerance { .. })
        ));
    }

    #[test]
    fn validate_rejects_neg_infinity_fragment_tolerance() {
        let mut p = valid_params();
        p.fragment_tolerance.value = f64::NEG_INFINITY;
        assert!(matches!(
            p.validate(),
            Err(SearchParamsError::InvalidFragmentTolerance { .. })
        ));
    }

    #[test]
    fn validate_rejects_nan_fixed_modification_mass_delta() {
        let mut p = valid_params();
        p.fixed_modifications[0].mass_delta = f64::NAN;
        let err = p.validate().unwrap_err();
        assert!(matches!(
            err,
            SearchParamsError::InvalidModificationMassDelta { .. }
        ));
        assert!(err.to_string().contains("Carbamidomethyl"));
    }

    #[test]
    fn validate_rejects_infinity_variable_modification_mass_delta() {
        let mut p = valid_params();
        p.variable_modifications[0].mass_delta = f64::INFINITY;
        let err = p.validate().unwrap_err();
        assert!(matches!(
            err,
            SearchParamsError::InvalidModificationMassDelta { .. }
        ));
        assert!(err.to_string().contains("Oxidation"));
    }

    // -- New field validation ------------------------------------------

    #[test]
    fn new_fields_serde_default() {
        let json = r#"{
            "enzyme": "Trypsin",
            "missed_cleavages": 2,
            "fixed_modifications": [],
            "variable_modifications": [],
            "precursor_tolerance": {"value": 10.0, "unit": "Ppm"},
            "fragment_tolerance": {"value": 0.02, "unit": "Da"},
            "database_path": "test.fasta",
            "decoy_strategy": "Reverse"
        }"#;
        let params: SearchParams = serde_json::from_str(json).unwrap();
        assert_eq!(params.max_variable_modifications, 3);
        assert_eq!(params.min_peptide_length, 7);
        assert_eq!(params.max_peptide_length, 50);
    }

    #[test]
    fn validate_rejects_excessive_variable_mods() {
        let mut p = valid_params();
        p.max_variable_modifications = 11;
        assert!(matches!(
            p.validate(),
            Err(SearchParamsError::TooManyVariableMods {
                actual: 11,
                max: 10
            })
        ));
    }

    #[test]
    fn validate_allows_max_variable_mods() {
        let mut p = valid_params();
        p.max_variable_modifications = 10;
        assert!(p.validate().is_ok());
    }

    #[test]
    fn validate_rejects_zero_peptide_length() {
        let mut p = valid_params();
        p.min_peptide_length = 0;
        assert!(matches!(
            p.validate(),
            Err(SearchParamsError::ZeroPeptideLength)
        ));
    }

    #[test]
    fn validate_rejects_min_gt_max_peptide_length() {
        let mut p = valid_params();
        p.min_peptide_length = 20;
        p.max_peptide_length = 10;
        assert!(matches!(
            p.validate(),
            Err(SearchParamsError::InvalidPeptideLengthRange { min: 20, max: 10 })
        ));
    }

    #[test]
    fn validate_allows_equal_min_max_peptide_length() {
        let mut p = valid_params();
        p.min_peptide_length = 10;
        p.max_peptide_length = 10;
        assert!(p.validate().is_ok());
    }
}
