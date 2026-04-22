//! YAML configuration parsing for entrapment analysis.

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::EntrapmentError;

/// Top-level configuration for an entrapment analysis run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EntrapmentConfig {
    /// Schema version – must be `1`.
    pub version: u32,

    /// Rules that classify a protein accession as **target**.
    pub target: GroupConfig,

    /// Rules that classify a protein accession as **trap / entrapment**.
    pub trap: GroupConfig,

    /// How to resolve a protein that matches *both* target and trap rules.
    #[serde(default)]
    pub conflict_resolution: ConflictResolution,

    /// What to do with proteins that match *neither* target nor trap rules.
    #[serde(default)]
    pub unmatched: UnmatchedPolicy,

    /// Parameters governing sequence-similarity / homology scoring.
    #[serde(default)]
    pub similarity: SimilarityConfig,

    /// Parameters for fragment ion provenance analysis (v3).
    #[serde(default)]
    pub provenance: ProvenanceConfig,
}

/// Configuration for a single classification group (target **or** trap).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GroupConfig {
    /// Ordered list of classification rules – a protein matches the group if
    /// **any** rule matches.
    pub rules: Vec<Rule>,

    /// Optional list of FASTA files whose accessions belong to this group.
    #[serde(default)]
    pub fasta: Option<Vec<FastaRef>>,

    /// Optional path to a plain-text file (one accession per line).
    #[serde(default)]
    pub accession_list: Option<PathBuf>,
}

/// Reference to a FASTA file used for accession loading.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FastaRef {
    /// Path to the FASTA file.
    pub path: PathBuf,
}

/// A single classification rule.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum Rule {
    /// Match if the accession string contains **any** of the given substrings.
    AccessionContains {
        /// Substrings to search for (logical OR).
        any_of: Vec<String>,
    },

    /// Match if the accession string matches the given regex.
    AccessionRegex {
        /// Regex pattern (Rust `regex` crate syntax).
        pattern: String,
    },

    /// Match if the accession appears in the given FASTA file.
    Fasta {
        /// Path to the FASTA file.
        path: PathBuf,
    },

    /// Match if the accession appears in the given plain-text list
    /// (one accession per line).
    AccessionList {
        /// Path to the accession list file.
        path: PathBuf,
    },
}

/// Strategy for resolving proteins that match **both** target and trap rules.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ConflictResolution {
    /// Classify as target (default).
    #[default]
    PreferTarget,
    /// Classify as trap.
    PreferTrap,
    /// Mark the protein as ambiguous.
    MarkAmbiguous,
}

/// Policy for proteins that match **neither** target nor trap rules.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum UnmatchedPolicy {
    /// Silently ignore unmatched proteins (default).
    #[default]
    Ignore,
    /// Treat unmatched proteins as trap.
    Trap,
    /// Treat unmatched proteins as target.
    Target,
    /// Return an error if any protein is unmatched.
    Error,
}

/// Parameters controlling sequence-similarity / homology scoring.
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

/// Configuration for fragment ion provenance analysis (v3).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProvenanceConfig {
    /// Fragment-ion mass tolerance in ppm.
    #[serde(default = "default_fragment_tolerance_ppm")]
    pub fragment_tolerance_ppm: f64,

    /// Maximum charge state to consider for fragment ions.
    #[serde(default = "default_max_fragment_charge")]
    pub max_fragment_charge: i32,

    /// Fraction of explained intensity above which a spectrum is flagged as chimeric.
    #[serde(default = "default_chimera_threshold")]
    pub chimera_threshold: f64,

    /// Minimum number of peaks required to run provenance analysis on a spectrum.
    #[serde(default = "default_min_peaks_for_analysis")]
    pub min_peaks_for_analysis: u32,

    /// Which similarity levels to include in the provenance trace.
    #[serde(default = "default_levels_to_trace")]
    pub levels_to_trace: Vec<String>,

    /// Fallback RT tolerance in minutes when `scan_number` is unavailable.
    ///
    /// Used with `find_by_rt` to locate the matching MS2 scan.
    /// If the PSM carries `rt_start` / `rt_stop`, the tolerance is derived
    /// from those instead: `(rt_stop − rt_start) / 2`.
    #[serde(default = "default_rt_tolerance_min")]
    pub rt_tolerance_min: f64,
}

impl Default for ProvenanceConfig {
    fn default() -> Self {
        Self {
            fragment_tolerance_ppm: default_fragment_tolerance_ppm(),
            max_fragment_charge: default_max_fragment_charge(),
            chimera_threshold: default_chimera_threshold(),
            min_peaks_for_analysis: default_min_peaks_for_analysis(),
            levels_to_trace: default_levels_to_trace(),
            rt_tolerance_min: default_rt_tolerance_min(),
        }
    }
}

fn default_fragment_tolerance_ppm() -> f64 {
    20.0
}
fn default_max_fragment_charge() -> i32 {
    2
}
fn default_chimera_threshold() -> f64 {
    0.3
}
fn default_min_peaks_for_analysis() -> u32 {
    6
}
fn default_levels_to_trace() -> Vec<String> {
    vec!["L2".to_string(), "L3".to_string(), "L4".to_string()]
}

fn default_rt_tolerance_min() -> f64 {
    0.5 // 30 seconds, reasonable default for DIA elution windows
}

fn default_max_mismatches() -> u16 {
    2
}
fn default_delta_mass_threshold_da() -> f64 {
    1.0
}
fn default_require_tryptic_ends() -> bool {
    true
}
fn default_max_missed_cleavages() -> u32 {
    2
}
fn default_len_tolerance() -> usize {
    2
}
fn default_true() -> bool {
    true
}

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

impl EntrapmentConfig {
    /// Read and parse an [`EntrapmentConfig`] from a YAML file on disk.
    ///
    /// The file is read to a string, parsed, and then validated.
    pub fn from_yaml(path: &Path) -> Result<Self, EntrapmentError> {
        let contents =
            std::fs::read_to_string(path).map_err(|e| EntrapmentError::ConfigIoError {
                path: path.to_path_buf(),
                detail: e.to_string(),
            })?;
        let config: Self =
            serde_yaml::from_str(&contents).map_err(|e| EntrapmentError::ConfigError {
                detail: format!("YAML parse error in {}: {e}", path.display()),
            })?;
        config.validate()?;
        Ok(config)
    }

    /// Parse an [`EntrapmentConfig`] from a YAML string.
    pub fn from_yaml_str(yaml: &str) -> Result<Self, EntrapmentError> {
        let config: Self =
            serde_yaml::from_str(yaml).map_err(|e| EntrapmentError::ConfigError {
                detail: format!("YAML parse error: {e}"),
            })?;
        config.validate()?;
        Ok(config)
    }

    /// Validate the parsed configuration.
    ///
    /// - `version` must be `1`.
    /// - Both `target.rules` and `trap.rules` must be non-empty.
    pub fn validate(&self) -> Result<(), EntrapmentError> {
        if self.version != 1 {
            return Err(EntrapmentError::ConfigError {
                detail: format!("unsupported config version {}, expected 1", self.version),
            });
        }
        if self.target.rules.is_empty() {
            return Err(EntrapmentError::ConfigError {
                detail: "target group must have at least one rule".to_string(),
            });
        }
        if self.trap.rules.is_empty() {
            return Err(EntrapmentError::ConfigError {
                detail: "trap group must have at least one rule".to_string(),
            });
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const MINIMAL_YAML: &str = r#"
version: 1
target:
  rules:
    - type: AccessionContains
      any_of: ["_HUMAN"]
trap:
  rules:
    - type: AccessionContains
      any_of: ["_YEAST", "_ECOLI"]
"#;

    const FULL_YAML: &str = r#"
version: 1
target:
  rules:
    - type: AccessionContains
      any_of: ["_HUMAN"]
    - type: AccessionRegex
      pattern: "^sp\\|.*_HUMAN"
trap:
  rules:
    - type: AccessionContains
      any_of: ["_YEAST", "_ECOLI", "_DICDI"]
conflict_resolution: mark_ambiguous
unmatched: trap
similarity:
  max_mismatches: 3
  delta_mz_threshold_da: 0.5
  require_tryptic_ends: false
  max_missed_cleavages: 1
"#;

    #[test]
    fn test_parse_minimal_config() {
        let cfg =
            EntrapmentConfig::from_yaml_str(MINIMAL_YAML).expect("minimal config should parse");

        assert_eq!(cfg.version, 1);
        assert_eq!(cfg.target.rules.len(), 1);
        assert_eq!(cfg.trap.rules.len(), 1);

        // Defaults
        assert_eq!(cfg.conflict_resolution, ConflictResolution::PreferTarget);
        assert_eq!(cfg.unmatched, UnmatchedPolicy::Ignore);
        assert_eq!(cfg.similarity.max_mismatches, 2);
        assert!((cfg.similarity.delta_mass_threshold_da - 1.0).abs() < f64::EPSILON);
        assert!(cfg.similarity.require_tryptic_ends);
        assert_eq!(cfg.similarity.max_missed_cleavages, 2);
    }

    #[test]
    fn test_parse_full_config() {
        let cfg = EntrapmentConfig::from_yaml_str(FULL_YAML).expect("full config should parse");

        assert_eq!(cfg.version, 1);
        assert_eq!(cfg.target.rules.len(), 2);
        assert_eq!(cfg.trap.rules.len(), 1);

        // Check target rules
        match &cfg.target.rules[0] {
            Rule::AccessionContains { any_of } => {
                assert_eq!(any_of, &["_HUMAN"]);
            }
            other => panic!("expected AccessionContains, got {other:?}"),
        }
        match &cfg.target.rules[1] {
            Rule::AccessionRegex { pattern } => {
                assert_eq!(pattern, r"^sp\|.*_HUMAN");
            }
            other => panic!("expected AccessionRegex, got {other:?}"),
        }

        // Check trap rules
        match &cfg.trap.rules[0] {
            Rule::AccessionContains { any_of } => {
                assert_eq!(any_of, &["_YEAST", "_ECOLI", "_DICDI"]);
            }
            other => panic!("expected AccessionContains, got {other:?}"),
        }

        // Explicit overrides
        assert_eq!(cfg.conflict_resolution, ConflictResolution::MarkAmbiguous);
        assert_eq!(cfg.unmatched, UnmatchedPolicy::Trap);
        assert_eq!(cfg.similarity.max_mismatches, 3);
        assert!((cfg.similarity.delta_mass_threshold_da - 0.5).abs() < f64::EPSILON);
        assert!(!cfg.similarity.require_tryptic_ends);
        assert_eq!(cfg.similarity.max_missed_cleavages, 1);
    }

    #[test]
    fn test_reject_bad_version() {
        let yaml = r#"
version: 2
target:
  rules:
    - type: AccessionContains
      any_of: ["_HUMAN"]
trap:
  rules:
    - type: AccessionContains
      any_of: ["_YEAST"]
"#;
        let err = EntrapmentConfig::from_yaml_str(yaml).expect_err("version 2 should be rejected");
        let msg = err.to_string();
        assert!(
            msg.contains("unsupported config version 2"),
            "unexpected error message: {msg}"
        );
    }

    #[test]
    fn test_reject_empty_rules() {
        let yaml = r#"
version: 1
target:
  rules: []
trap:
  rules:
    - type: AccessionContains
      any_of: ["_YEAST"]
"#;
        let err = EntrapmentConfig::from_yaml_str(yaml)
            .expect_err("empty target rules should be rejected");
        let msg = err.to_string();
        assert!(
            msg.contains("target group must have at least one rule"),
            "unexpected error message: {msg}"
        );
    }

    #[test]
    fn test_v2_defaults() {
        let cfg = EntrapmentConfig::from_yaml_str(MINIMAL_YAML).expect("parse");
        assert_eq!(cfg.similarity.len_tolerance, 2);
        assert!(cfg.similarity.enable_dipeptide_check);
        assert!(cfg.similarity.enable_qk_detection);
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
    fn provenance_config_defaults() {
        let yaml = "version: 1\ntarget:\n  rules: []\ntrap:\n  rules: []\n";
        let config: EntrapmentConfig = serde_yaml::from_str(yaml).unwrap();
        assert!((config.provenance.fragment_tolerance_ppm - 20.0).abs() < 1e-6);
        assert_eq!(config.provenance.max_fragment_charge, 2);
        assert!((config.provenance.chimera_threshold - 0.3).abs() < 1e-6);
        assert_eq!(config.provenance.min_peaks_for_analysis, 6);
        assert_eq!(config.provenance.levels_to_trace, vec!["L2", "L3", "L4"]);
    }

    #[test]
    fn provenance_config_custom_values() {
        let yaml = "version: 1\ntarget:\n  rules: []\ntrap:\n  rules: []\nprovenance:\n  fragment_tolerance_ppm: 10.0\n  max_fragment_charge: 3\n  chimera_threshold: 0.5\n  min_peaks_for_analysis: 10\n  levels_to_trace: [\"L3\", \"L4\"]\n";
        let config: EntrapmentConfig = serde_yaml::from_str(yaml).unwrap();
        assert!((config.provenance.fragment_tolerance_ppm - 10.0).abs() < 1e-6);
        assert_eq!(config.provenance.max_fragment_charge, 3);
        assert!((config.provenance.chimera_threshold - 0.5).abs() < 1e-6);
        assert_eq!(config.provenance.min_peaks_for_analysis, 10);
        assert_eq!(config.provenance.levels_to_trace, vec!["L3", "L4"]);
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
}
