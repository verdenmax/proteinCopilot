//! Target/trap classification tagger.
//!
//! Applies target and trap rules from the configuration to classify PSMs
//! by their protein accessions.  Each accession in a semicolon-separated
//! list is tested against the configured matchers and accession sets;
//! conflicts and unmatched accessions are resolved according to the
//! [`ConflictResolution`] and [`UnmatchedPolicy`] settings.

use std::collections::HashSet;
use std::io::BufRead;
use std::path::Path;

use regex::Regex;
use tracing::debug;

use protein_copilot_search_engine::fasta::parse_fasta;

use crate::config::{ConflictResolution, EntrapmentConfig, Rule, UnmatchedPolicy};
use crate::error::EntrapmentError;
use crate::types::PsmGroup;

// ---------------------------------------------------------------------------
// Internal matcher
// ---------------------------------------------------------------------------

/// A compiled classification rule used for substring or regex matching.
#[derive(Debug)]
enum Matcher {
    /// Matches if the accession string contains **any** of the given substrings.
    Contains(Vec<String>),
    /// Matches if the accession string matches the compiled regex.
    Regex(Regex),
}

// ---------------------------------------------------------------------------
// Tagger
// ---------------------------------------------------------------------------

/// Classifies protein accessions as *target*, *trap*, or *ambiguous*.
///
/// Built from an [`EntrapmentConfig`]; call [`Tagger::tag`] to classify a
/// semicolon-separated protein-ID string.
#[derive(Debug)]
pub struct Tagger {
    target_matchers: Vec<Matcher>,
    trap_matchers: Vec<Matcher>,
    target_accessions: HashSet<String>,
    trap_accessions: HashSet<String>,
    conflict_resolution: ConflictResolution,
    unmatched: UnmatchedPolicy,
}

impl Tagger {
    /// Build a new [`Tagger`] from the given configuration.
    ///
    /// This compiles all regex patterns, loads FASTA accessions, and reads
    /// accession-list files.  Returns [`EntrapmentError::ConfigError`] for
    /// invalid regex patterns, or I/O-related errors for inaccessible files.
    pub fn new(config: &EntrapmentConfig) -> Result<Self, EntrapmentError> {
        let mut target_matchers = Vec::new();
        let mut target_accessions = HashSet::new();
        let mut trap_matchers = Vec::new();
        let mut trap_accessions = HashSet::new();

        // --- target rules ---
        for rule in &config.target.rules {
            match rule {
                Rule::AccessionContains { any_of } => {
                    target_matchers.push(Matcher::Contains(any_of.clone()));
                }
                Rule::AccessionRegex { pattern } => {
                    let re = Regex::new(pattern).map_err(|e| EntrapmentError::ConfigError {
                        detail: format!("invalid target regex '{pattern}': {e}"),
                    })?;
                    target_matchers.push(Matcher::Regex(re));
                }
                Rule::Fasta { path } => {
                    load_fasta_accessions(path, &mut target_accessions)?;
                }
                Rule::AccessionList { path } => {
                    load_accession_list(path, &mut target_accessions)?;
                }
            }
        }

        // target GroupConfig.fasta
        if let Some(fasta_refs) = &config.target.fasta {
            for fasta_ref in fasta_refs {
                load_fasta_accessions(&fasta_ref.path, &mut target_accessions)?;
            }
        }
        // target GroupConfig.accession_list
        if let Some(list_path) = &config.target.accession_list {
            load_accession_list(list_path, &mut target_accessions)?;
        }

        // --- trap rules ---
        for rule in &config.trap.rules {
            match rule {
                Rule::AccessionContains { any_of } => {
                    trap_matchers.push(Matcher::Contains(any_of.clone()));
                }
                Rule::AccessionRegex { pattern } => {
                    let re = Regex::new(pattern).map_err(|e| EntrapmentError::ConfigError {
                        detail: format!("invalid trap regex '{pattern}': {e}"),
                    })?;
                    trap_matchers.push(Matcher::Regex(re));
                }
                Rule::Fasta { path } => {
                    load_fasta_accessions(path, &mut trap_accessions)?;
                }
                Rule::AccessionList { path } => {
                    load_accession_list(path, &mut trap_accessions)?;
                }
            }
        }

        // trap GroupConfig.fasta
        if let Some(fasta_refs) = &config.trap.fasta {
            for fasta_ref in fasta_refs {
                load_fasta_accessions(&fasta_ref.path, &mut trap_accessions)?;
            }
        }
        // trap GroupConfig.accession_list
        if let Some(list_path) = &config.trap.accession_list {
            load_accession_list(list_path, &mut trap_accessions)?;
        }

        debug!(
            target_matchers = target_matchers.len(),
            trap_matchers = trap_matchers.len(),
            target_accessions = target_accessions.len(),
            trap_accessions = trap_accessions.len(),
            "tagger initialised"
        );

        Ok(Self {
            target_matchers,
            trap_matchers,
            target_accessions,
            trap_accessions,
            conflict_resolution: config.conflict_resolution,
            unmatched: config.unmatched,
        })
    }

    /// Classify a semicolon-separated list of protein accessions.
    ///
    /// Returns [`PsmGroup::Target`], [`PsmGroup::Trap`], or
    /// [`PsmGroup::Ambiguous`] depending on the match results and the
    /// configured conflict-resolution / unmatched policies.
    pub fn tag(&self, protein_ids: &str) -> Result<PsmGroup, EntrapmentError> {
        let mut is_target = false;
        let mut is_trap = false;

        for accession in protein_ids.split(';') {
            let accession = accession.trim();
            if accession.is_empty() {
                continue;
            }

            if matches_any(accession, &self.target_matchers, &self.target_accessions) {
                is_target = true;
            }
            if matches_any(accession, &self.trap_matchers, &self.trap_accessions) {
                is_trap = true;
            }
        }

        match (is_target, is_trap) {
            (true, true) => match self.conflict_resolution {
                ConflictResolution::PreferTarget => Ok(PsmGroup::Target),
                ConflictResolution::PreferTrap => Ok(PsmGroup::Trap),
                ConflictResolution::MarkAmbiguous => Ok(PsmGroup::Ambiguous),
            },
            (true, false) => Ok(PsmGroup::Target),
            (false, true) => Ok(PsmGroup::Trap),
            (false, false) => match self.unmatched {
                UnmatchedPolicy::Ignore | UnmatchedPolicy::Target => Ok(PsmGroup::Target),
                UnmatchedPolicy::Trap => Ok(PsmGroup::Trap),
                UnmatchedPolicy::Error => Err(EntrapmentError::TaggingError {
                    detail: format!("unmatched protein accession(s): {protein_ids}"),
                }),
            },
        }
    }
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

/// Check whether `accession` is in the explicit accession set **or** matches
/// any of the compiled matchers.  Accession-set lookup (exact match) is
/// checked first for performance.
fn matches_any(accession: &str, matchers: &[Matcher], accession_set: &HashSet<String>) -> bool {
    // Exact accession set first
    if accession_set.contains(accession) {
        return true;
    }

    // Then pattern matchers
    for matcher in matchers {
        match matcher {
            Matcher::Contains(substrings) => {
                if substrings.iter().any(|s| accession.contains(s.as_str())) {
                    return true;
                }
            }
            Matcher::Regex(re) => {
                if re.is_match(accession) {
                    return true;
                }
            }
        }
    }
    false
}

/// Load accessions from a FASTA file via the search-engine crate parser.
///
/// For UniProt-format headers like `sp|P12345|ALBU_HUMAN`, inserts both
/// the full accession (`sp|P12345|ALBU_HUMAN`) and the bare UniProt ID
/// (`P12345`) so that either format in search results will match.
fn load_fasta_accessions(path: &Path, set: &mut HashSet<String>) -> Result<(), EntrapmentError> {
    let entries = parse_fasta(path).map_err(|e| EntrapmentError::FastaError {
        path: path.to_path_buf(),
        detail: e.to_string(),
    })?;
    for entry in entries {
        // Extract bare UniProt ID from "sp|P12345|NAME_SPECIES" or "tr|P12345|..."
        if let Some(bare_id) = extract_uniprot_id(&entry.accession) {
            set.insert(bare_id.to_string());
        }
        set.insert(entry.accession);
    }
    Ok(())
}

/// Extract the bare UniProt accession from a `db|ACCESSION|ENTRY_NAME` string.
///
/// Returns `Some("P12345")` for `"sp|P12345|ALBU_HUMAN"`, `None` for
/// accessions that don't follow this format.
fn extract_uniprot_id(accession: &str) -> Option<&str> {
    let parts: Vec<&str> = accession.split('|').collect();
    if parts.len() >= 2 && (parts[0] == "sp" || parts[0] == "tr") && !parts[1].is_empty() {
        Some(parts[1])
    } else {
        None
    }
}

/// Load accessions from a plain-text file (one per line).
///
/// Empty lines and lines starting with `#` are skipped.
fn load_accession_list(path: &Path, set: &mut HashSet<String>) -> Result<(), EntrapmentError> {
    let file = std::fs::File::open(path).map_err(|e| EntrapmentError::IoError {
        path: path.to_path_buf(),
        detail: e.to_string(),
    })?;
    let reader = std::io::BufReader::new(file);
    for line in reader.lines() {
        let line = line.map_err(|e| EntrapmentError::IoError {
            path: path.to_path_buf(),
            detail: e.to_string(),
        })?;
        let trimmed = line.trim().to_string();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        set.insert(trimmed);
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{
        ConflictResolution, EntrapmentConfig, GroupConfig, ProvenanceConfig, SimilarityConfig,
        UnmatchedPolicy,
    };

    /// Helper: build a minimal config from target/trap `AccessionContains` rules.
    fn make_config(
        target_substrings: Vec<String>,
        trap_substrings: Vec<String>,
        conflict: ConflictResolution,
        unmatched: UnmatchedPolicy,
    ) -> EntrapmentConfig {
        EntrapmentConfig {
            version: 1,
            target: GroupConfig {
                rules: vec![Rule::AccessionContains {
                    any_of: target_substrings,
                }],
                fasta: None,
                accession_list: None,
            },
            trap: GroupConfig {
                rules: vec![Rule::AccessionContains {
                    any_of: trap_substrings,
                }],
                fasta: None,
                accession_list: None,
            },
            conflict_resolution: conflict,
            unmatched,
            similarity: SimilarityConfig::default(),
            provenance: ProvenanceConfig::default(),
        }
    }

    #[test]
    fn test_simple_contains_rules() {
        let config = make_config(
            vec!["_HUMAN".into()],
            vec!["_YEAST".into()],
            ConflictResolution::PreferTarget,
            UnmatchedPolicy::Ignore,
        );
        let tagger = Tagger::new(&config).expect("tagger should build");

        assert_eq!(
            tagger.tag("sp|P12345|EF1A_HUMAN").expect("should classify"),
            PsmGroup::Target
        );
        assert_eq!(
            tagger.tag("sp|P99999|FAKE_YEAST").expect("should classify"),
            PsmGroup::Trap
        );
    }

    #[test]
    fn test_conflict_prefer_target() {
        // Protein matches both _HUMAN (target) and _HUMAN_YEAST (trap) via
        // substring — but we use a single accession that contains both.
        let config = make_config(
            vec!["_HUMAN".into()],
            vec!["_YEAST".into()],
            ConflictResolution::PreferTarget,
            UnmatchedPolicy::Ignore,
        );
        let tagger = Tagger::new(&config).expect("tagger should build");

        // Accession that matches both target ("_HUMAN") and trap ("_YEAST")
        assert_eq!(
            tagger
                .tag("sp|P00001|COMBO_HUMAN_YEAST")
                .expect("should classify"),
            PsmGroup::Target
        );
    }

    #[test]
    fn test_conflict_mark_ambiguous() {
        let config = make_config(
            vec!["_HUMAN".into()],
            vec!["_YEAST".into()],
            ConflictResolution::MarkAmbiguous,
            UnmatchedPolicy::Ignore,
        );
        let tagger = Tagger::new(&config).expect("tagger should build");

        assert_eq!(
            tagger
                .tag("sp|P00001|COMBO_HUMAN_YEAST")
                .expect("should classify"),
            PsmGroup::Ambiguous
        );
    }

    #[test]
    fn test_unmatched_error() {
        let config = make_config(
            vec!["_HUMAN".into()],
            vec!["_YEAST".into()],
            ConflictResolution::PreferTarget,
            UnmatchedPolicy::Error,
        );
        let tagger = Tagger::new(&config).expect("tagger should build");

        let result = tagger.tag("sp|P00001|UNKNOWN_MOUSE");
        assert!(result.is_err(), "unmatched should return error");
        let err_msg = result.unwrap_err().to_string();
        assert!(
            err_msg.contains("unmatched protein accession"),
            "unexpected error: {err_msg}"
        );
    }

    #[test]
    fn test_semicolon_separated_proteins() {
        let config = make_config(
            vec!["_HUMAN".into()],
            vec!["_YEAST".into()],
            ConflictResolution::PreferTarget,
            UnmatchedPolicy::Ignore,
        );
        let tagger = Tagger::new(&config).expect("tagger should build");

        // One target, one trap → conflict → PreferTarget
        assert_eq!(
            tagger
                .tag("sp|P12345|EF1A_HUMAN;sp|P99999|FAKE_YEAST")
                .expect("should classify"),
            PsmGroup::Target
        );

        // Same with MarkAmbiguous
        let config2 = make_config(
            vec!["_HUMAN".into()],
            vec!["_YEAST".into()],
            ConflictResolution::MarkAmbiguous,
            UnmatchedPolicy::Ignore,
        );
        let tagger2 = Tagger::new(&config2).expect("tagger should build");
        assert_eq!(
            tagger2
                .tag("sp|P12345|EF1A_HUMAN;sp|P99999|FAKE_YEAST")
                .expect("should classify"),
            PsmGroup::Ambiguous
        );
    }

    #[test]
    fn test_regex_rule() {
        let config = EntrapmentConfig {
            version: 1,
            target: GroupConfig {
                rules: vec![Rule::AccessionRegex {
                    pattern: r"^sp\|.*_HUMAN$".into(),
                }],
                fasta: None,
                accession_list: None,
            },
            trap: GroupConfig {
                rules: vec![Rule::AccessionRegex {
                    pattern: r"_YEAST$".into(),
                }],
                fasta: None,
                accession_list: None,
            },
            conflict_resolution: ConflictResolution::PreferTarget,
            unmatched: UnmatchedPolicy::Ignore,
            similarity: SimilarityConfig::default(),
            provenance: ProvenanceConfig::default(),
        };

        let tagger = Tagger::new(&config).expect("tagger should build");

        assert_eq!(
            tagger.tag("sp|P12345|EF1A_HUMAN").expect("should classify"),
            PsmGroup::Target
        );
        assert_eq!(
            tagger.tag("sp|P99999|FAKE_YEAST").expect("should classify"),
            PsmGroup::Trap
        );
        // "tr|..." does not match "^sp\|" → unmatched → Target (Ignore policy)
        assert_eq!(
            tagger
                .tag("tr|Q11111|THING_MOUSE")
                .expect("should classify"),
            PsmGroup::Target
        );
    }

    #[test]
    fn test_bad_regex_returns_config_error() {
        let config = EntrapmentConfig {
            version: 1,
            target: GroupConfig {
                rules: vec![Rule::AccessionRegex {
                    pattern: r"[invalid".into(),
                }],
                fasta: None,
                accession_list: None,
            },
            trap: GroupConfig {
                rules: vec![Rule::AccessionContains {
                    any_of: vec!["_YEAST".into()],
                }],
                fasta: None,
                accession_list: None,
            },
            conflict_resolution: ConflictResolution::PreferTarget,
            unmatched: UnmatchedPolicy::Ignore,
            similarity: SimilarityConfig::default(),
            provenance: ProvenanceConfig::default(),
        };

        let result = Tagger::new(&config);
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(
            err_msg.contains("invalid target regex"),
            "unexpected error: {err_msg}"
        );
    }

    #[test]
    fn test_unmatched_target_policy() {
        let config = make_config(
            vec!["_HUMAN".into()],
            vec!["_YEAST".into()],
            ConflictResolution::PreferTarget,
            UnmatchedPolicy::Target,
        );
        let tagger = Tagger::new(&config).expect("tagger should build");

        assert_eq!(
            tagger
                .tag("sp|P00001|UNKNOWN_MOUSE")
                .expect("should classify"),
            PsmGroup::Target
        );
    }

    #[test]
    fn test_unmatched_trap_policy() {
        let config = make_config(
            vec!["_HUMAN".into()],
            vec!["_YEAST".into()],
            ConflictResolution::PreferTarget,
            UnmatchedPolicy::Trap,
        );
        let tagger = Tagger::new(&config).expect("tagger should build");

        assert_eq!(
            tagger
                .tag("sp|P00001|UNKNOWN_MOUSE")
                .expect("should classify"),
            PsmGroup::Trap
        );
    }

    #[test]
    fn test_conflict_prefer_trap() {
        let config = make_config(
            vec!["_HUMAN".into()],
            vec!["_YEAST".into()],
            ConflictResolution::PreferTrap,
            UnmatchedPolicy::Ignore,
        );
        let tagger = Tagger::new(&config).expect("tagger should build");

        assert_eq!(
            tagger
                .tag("sp|P12345|EF1A_HUMAN;sp|P99999|FAKE_YEAST")
                .expect("should classify"),
            PsmGroup::Trap
        );
    }
}
