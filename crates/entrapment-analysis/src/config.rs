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
    /// Maximum number of amino-acid mismatches allowed before considering
    /// two peptides distinct.
    #[serde(default = "default_max_mismatches")]
    pub max_mismatches: u8,

    /// Precursor m/z tolerance (Da) for linking similar peptides.
    #[serde(default = "default_delta_mz_threshold_da")]
    pub delta_mz_threshold_da: f64,

    /// Whether both ends of a peptide must be tryptic.
    #[serde(default = "default_require_tryptic_ends")]
    pub require_tryptic_ends: bool,

    /// Maximum number of missed cleavages to allow.
    #[serde(default = "default_max_missed_cleavages")]
    pub max_missed_cleavages: u32,
}

fn default_max_mismatches() -> u8 {
    2
}
fn default_delta_mz_threshold_da() -> f64 {
    1.0
}
fn default_require_tryptic_ends() -> bool {
    true
}
fn default_max_missed_cleavages() -> u32 {
    2
}

impl Default for SimilarityConfig {
    fn default() -> Self {
        Self {
            max_mismatches: default_max_mismatches(),
            delta_mz_threshold_da: default_delta_mz_threshold_da(),
            require_tryptic_ends: default_require_tryptic_ends(),
            max_missed_cleavages: default_max_missed_cleavages(),
        }
    }
}

impl EntrapmentConfig {
    /// Read and parse an [`EntrapmentConfig`] from a YAML file on disk.
    ///
    /// The file is read to a string, parsed, and then validated.
    pub fn from_yaml(path: &Path) -> Result<Self, EntrapmentError> {
        let contents = std::fs::read_to_string(path).map_err(|e| {
            EntrapmentError::ConfigIoError {
                path: path.to_path_buf(),
                detail: e.to_string(),
            }
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
                detail: format!(
                    "unsupported config version {}, expected 1",
                    self.version
                ),
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
        let cfg = EntrapmentConfig::from_yaml_str(MINIMAL_YAML)
            .expect("minimal config should parse");

        assert_eq!(cfg.version, 1);
        assert_eq!(cfg.target.rules.len(), 1);
        assert_eq!(cfg.trap.rules.len(), 1);

        // Defaults
        assert_eq!(cfg.conflict_resolution, ConflictResolution::PreferTarget);
        assert_eq!(cfg.unmatched, UnmatchedPolicy::Ignore);
        assert_eq!(cfg.similarity.max_mismatches, 2);
        assert!((cfg.similarity.delta_mz_threshold_da - 1.0).abs() < f64::EPSILON);
        assert!(cfg.similarity.require_tryptic_ends);
        assert_eq!(cfg.similarity.max_missed_cleavages, 2);
    }

    #[test]
    fn test_parse_full_config() {
        let cfg = EntrapmentConfig::from_yaml_str(FULL_YAML)
            .expect("full config should parse");

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
        assert!((cfg.similarity.delta_mz_threshold_da - 0.5).abs() < f64::EPSILON);
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
        let err = EntrapmentConfig::from_yaml_str(yaml)
            .expect_err("version 2 should be rejected");
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
}
