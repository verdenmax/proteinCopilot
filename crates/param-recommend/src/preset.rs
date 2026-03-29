//! Search presets — predefined parameter combinations for common scenarios.

use protein_copilot_core::search_params::{
    DecoyStrategy, Enzyme, MassTolerance, ModPosition, Modification, SearchParams, ToleranceUnit,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

/// A named preset of search parameters for a common experimental scenario.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct SearchPreset {
    /// Preset identifier (e.g., "standard", "phospho").
    pub name: String,
    /// Human-readable description.
    pub description: String,
    /// The recommended search parameters.
    pub params: SearchParams,
    /// Scenarios where this preset is applicable.
    pub applicable_scenarios: Vec<String>,
}

impl SearchPreset {
    /// Returns a copy of this preset with the database path set.
    ///
    /// Presets use a placeholder `"<database_path>"` by default.
    /// Call this method to produce a usable `SearchParams`.
    pub fn with_database(&self, database_path: &str) -> SearchParams {
        let mut params = self.params.clone();
        params.database_path = database_path.to_string();
        params
    }
}

// ---------------------------------------------------------------------------
// Common modifications (reused across presets)
// ---------------------------------------------------------------------------

fn carbamidomethyl_c() -> Modification {
    Modification {
        name: "Carbamidomethyl".to_string(),
        mass_delta: 57.021464,
        residues: vec!['C'],
        position: ModPosition::Anywhere,
    }
}

fn oxidation_m() -> Modification {
    Modification {
        name: "Oxidation".to_string(),
        mass_delta: 15.994915,
        residues: vec!['M'],
        position: ModPosition::Anywhere,
    }
}

fn phospho_sty() -> Modification {
    Modification {
        name: "Phospho".to_string(),
        mass_delta: 79.966331,
        residues: vec!['S', 'T', 'Y'],
        position: ModPosition::Anywhere,
    }
}

fn tmt6plex_k() -> Modification {
    Modification {
        name: "TMT6plex".to_string(),
        mass_delta: 229.162932,
        residues: vec!['K'],
        position: ModPosition::Anywhere,
    }
}

fn tmt6plex_nterm() -> Modification {
    Modification {
        name: "TMT6plex".to_string(),
        mass_delta: 229.162932,
        residues: vec![],
        position: ModPosition::AnyNTerm,
    }
}

/// Placeholder database path — must be replaced by the caller.
const PLACEHOLDER_DB: &str = "<database_path>";

// ---------------------------------------------------------------------------
// Built-in presets
// ---------------------------------------------------------------------------

/// Standard proteomics search preset.
pub fn standard_preset() -> SearchPreset {
    SearchPreset {
        name: "standard".to_string(),
        description: "Standard protein search: Trypsin, Carbamidomethyl(C) fixed, Oxidation(M) variable, 10ppm/0.02Da".to_string(),
        params: SearchParams {
            database_path: PLACEHOLDER_DB.to_string(),
            enzyme: Enzyme::Trypsin,
            missed_cleavages: 2,
            fixed_modifications: vec![carbamidomethyl_c()],
            variable_modifications: vec![oxidation_m()],
            precursor_tolerance: MassTolerance { value: 10.0, unit: ToleranceUnit::Ppm },
            fragment_tolerance: MassTolerance { value: 0.02, unit: ToleranceUnit::Da },
            decoy_strategy: DecoyStrategy::Reverse,
        },
        applicable_scenarios: vec![
            "General protein identification".to_string(),
            "HeLa cell line".to_string(),
            "Standard shotgun proteomics".to_string(),
        ],
    }
}

/// Phosphoproteomics search preset.
pub fn phospho_preset() -> SearchPreset {
    SearchPreset {
        name: "phospho".to_string(),
        description: "Phosphoproteomics: standard + Phospho(STY) variable modification".to_string(),
        params: SearchParams {
            database_path: PLACEHOLDER_DB.to_string(),
            enzyme: Enzyme::Trypsin,
            missed_cleavages: 2,
            fixed_modifications: vec![carbamidomethyl_c()],
            variable_modifications: vec![oxidation_m(), phospho_sty()],
            precursor_tolerance: MassTolerance {
                value: 10.0,
                unit: ToleranceUnit::Ppm,
            },
            fragment_tolerance: MassTolerance {
                value: 0.02,
                unit: ToleranceUnit::Da,
            },
            decoy_strategy: DecoyStrategy::Reverse,
        },
        applicable_scenarios: vec![
            "Phosphoproteomics".to_string(),
            "Enriched phosphopeptides".to_string(),
            "Signaling pathway analysis".to_string(),
        ],
    }
}

/// TMT-labeled search preset.
pub fn tmt_preset() -> SearchPreset {
    SearchPreset {
        name: "tmt".to_string(),
        description: "TMT6plex labeled search: TMT on K and N-term as fixed modifications"
            .to_string(),
        params: SearchParams {
            database_path: PLACEHOLDER_DB.to_string(),
            enzyme: Enzyme::Trypsin,
            missed_cleavages: 2,
            fixed_modifications: vec![carbamidomethyl_c(), tmt6plex_k(), tmt6plex_nterm()],
            variable_modifications: vec![oxidation_m()],
            precursor_tolerance: MassTolerance {
                value: 10.0,
                unit: ToleranceUnit::Ppm,
            },
            fragment_tolerance: MassTolerance {
                value: 0.02,
                unit: ToleranceUnit::Da,
            },
            decoy_strategy: DecoyStrategy::Reverse,
        },
        applicable_scenarios: vec![
            "TMT labeling quantification".to_string(),
            "Multiplexed proteomics".to_string(),
        ],
    }
}

/// Open search preset (wide precursor tolerance for PTM discovery).
pub fn open_search_preset() -> SearchPreset {
    SearchPreset {
        name: "open".to_string(),
        description: "Open search: 500 Da precursor tolerance for unknown modification discovery"
            .to_string(),
        params: SearchParams {
            database_path: PLACEHOLDER_DB.to_string(),
            enzyme: Enzyme::Trypsin,
            missed_cleavages: 2,
            fixed_modifications: vec![carbamidomethyl_c()],
            variable_modifications: vec![oxidation_m()],
            precursor_tolerance: MassTolerance {
                value: 500.0,
                unit: ToleranceUnit::Da,
            },
            fragment_tolerance: MassTolerance {
                value: 0.02,
                unit: ToleranceUnit::Da,
            },
            decoy_strategy: DecoyStrategy::Reverse,
        },
        applicable_scenarios: vec![
            "Unknown modification discovery".to_string(),
            "Mass shift analysis".to_string(),
            "PTM-centric proteomics".to_string(),
        ],
    }
}

/// Returns all built-in search presets.
pub fn all_presets() -> Vec<SearchPreset> {
    vec![
        standard_preset(),
        phospho_preset(),
        tmt_preset(),
        open_search_preset(),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_presets_have_4_entries() {
        assert_eq!(all_presets().len(), 4);
    }

    #[test]
    fn preset_names_unique() {
        let presets = all_presets();
        let mut names: Vec<&str> = presets.iter().map(|p| p.name.as_str()).collect();
        names.sort();
        names.dedup();
        assert_eq!(names.len(), 4);
    }

    #[test]
    fn presets_serde_roundtrip() {
        for preset in all_presets() {
            let json = serde_json::to_string(&preset).unwrap();
            let back: SearchPreset = serde_json::from_str(&json).unwrap();
            assert_eq!(preset.name, back.name);
            assert_eq!(preset.params.enzyme, back.params.enzyme);
        }
    }

    #[test]
    fn phospho_preset_has_phospho_mod() {
        let p = phospho_preset();
        assert!(p
            .params
            .variable_modifications
            .iter()
            .any(|m| m.name == "Phospho"));
    }

    #[test]
    fn tmt_preset_has_tmt_fixed() {
        let p = tmt_preset();
        assert!(p
            .params
            .fixed_modifications
            .iter()
            .any(|m| m.name == "TMT6plex"));
    }

    #[test]
    fn open_preset_has_wide_tolerance() {
        let p = open_search_preset();
        assert!(p.params.precursor_tolerance.value >= 500.0);
        assert_eq!(p.params.precursor_tolerance.unit, ToleranceUnit::Da);
    }
}
