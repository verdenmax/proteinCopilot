//! SearchParams → sage Parameters configuration builder.

use std::collections::HashMap;

use protein_copilot_core::search_params::{DecoyStrategy, Enzyme, SearchParams};
use sage_core::database::{Builder, EnzymeBuilder, Parameters};
use sage_core::ion_series::Kind;
use sage_core::modification::ModificationSpecificity;

use super::convert::{fixed_mod_to_sage, variable_mod_to_sage};

/// Build sage `Parameters` from our `SearchParams`.
///
/// The returned Parameters is ready to call `.build(fasta)` to create an IndexedDatabase.
pub fn build_sage_parameters(params: &SearchParams) -> Parameters {
    let enzyme = match &params.enzyme {
        Enzyme::Trypsin => EnzymeBuilder {
            cleave_at: Some("KR".into()),
            restrict: Some("P".into()),
            c_terminal: Some(true),
            missed_cleavages: Some(params.missed_cleavages as u8),
            min_len: Some(params.min_peptide_length as usize),
            max_len: Some(params.max_peptide_length as usize),
            semi_enzymatic: Some(false),
        },
        Enzyme::TrypsinP => EnzymeBuilder {
            cleave_at: Some("KR".into()),
            restrict: Some("".into()),
            c_terminal: Some(true),
            missed_cleavages: Some(params.missed_cleavages as u8),
            min_len: Some(params.min_peptide_length as usize),
            max_len: Some(params.max_peptide_length as usize),
            semi_enzymatic: Some(false),
        },
        Enzyme::LysC => EnzymeBuilder {
            cleave_at: Some("K".into()),
            restrict: Some("P".into()),
            c_terminal: Some(true),
            missed_cleavages: Some(params.missed_cleavages as u8),
            min_len: Some(params.min_peptide_length as usize),
            max_len: Some(params.max_peptide_length as usize),
            semi_enzymatic: Some(false),
        },
        Enzyme::GluC => EnzymeBuilder {
            cleave_at: Some("DE".into()),
            restrict: Some("".into()),
            c_terminal: Some(true),
            missed_cleavages: Some(params.missed_cleavages as u8),
            min_len: Some(params.min_peptide_length as usize),
            max_len: Some(params.max_peptide_length as usize),
            semi_enzymatic: Some(false),
        },
        Enzyme::AspN => EnzymeBuilder {
            cleave_at: Some("D".into()),
            restrict: Some("".into()),
            c_terminal: Some(false),
            missed_cleavages: Some(params.missed_cleavages as u8),
            min_len: Some(params.min_peptide_length as usize),
            max_len: Some(params.max_peptide_length as usize),
            semi_enzymatic: Some(false),
        },
        Enzyme::Chymotrypsin => EnzymeBuilder {
            cleave_at: Some("FWYL".into()),
            restrict: Some("".into()),
            c_terminal: Some(true),
            missed_cleavages: Some(params.missed_cleavages as u8),
            min_len: Some(params.min_peptide_length as usize),
            max_len: Some(params.max_peptide_length as usize),
            semi_enzymatic: Some(false),
        },
        Enzyme::NonSpecific => EnzymeBuilder {
            cleave_at: Some("".into()),
            restrict: Some("".into()),
            c_terminal: Some(true),
            missed_cleavages: Some(0),
            min_len: Some(params.min_peptide_length as usize),
            max_len: Some(params.max_peptide_length as usize),
            semi_enzymatic: Some(true),
        },
        Enzyme::Custom { cleavage_rule, .. } => EnzymeBuilder {
            cleave_at: Some(cleavage_rule.clone()),
            restrict: Some("".into()),
            c_terminal: Some(true),
            missed_cleavages: Some(params.missed_cleavages as u8),
            min_len: Some(params.min_peptide_length as usize),
            max_len: Some(params.max_peptide_length as usize),
            semi_enzymatic: Some(false),
        },
    };

    // Build static mods
    let mut static_mods: HashMap<ModificationSpecificity, f32> = HashMap::new();
    for m in &params.fixed_modifications {
        for (spec, mass) in fixed_mod_to_sage(m) {
            static_mods.insert(spec, mass);
        }
    }

    // Build variable mods
    let mut variable_mods: HashMap<ModificationSpecificity, Vec<f32>> = HashMap::new();
    for m in &params.variable_modifications {
        for (spec, masses) in variable_mod_to_sage(m) {
            variable_mods.entry(spec).or_default().extend(masses);
        }
    }

    let generate_decoys = params.decoy_strategy != DecoyStrategy::None;

    // Builder.static_mods takes Option<HashMap<String, f32>> — keys are string representations
    // of ModificationSpecificity. Builder::make_parameters() internally calls validate_mods()
    // which parses strings back into ModificationSpecificity.
    Builder {
        bucket_size: None,
        enzyme: Some(enzyme),
        peptide_min_mass: Some(500.0),
        peptide_max_mass: Some(5000.0),
        ion_kinds: Some(vec![Kind::B, Kind::Y]),
        min_ion_index: Some(2),
        static_mods: Some(
            static_mods
                .into_iter()
                .map(|(k, v)| (k.to_string(), v))
                .collect(),
        ),
        variable_mods: Some(
            variable_mods
                .into_iter()
                .map(|(k, v)| (k.to_string(), v))
                .collect(),
        ),
        max_variable_mods: Some(params.max_variable_modifications as usize),
        decoy_tag: Some("rev_".into()),
        generate_decoys: Some(generate_decoys),
        fasta: Some(params.database_path.clone()),
        prefilter_chunk_size: None,
        prefilter: None,
        prefilter_low_memory: None,
    }
    .make_parameters()
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::search_params::{
        DecoyStrategy, MassTolerance, ModPosition, Modification, ToleranceUnit,
    };

    fn test_params() -> SearchParams {
        SearchParams {
            enzyme: Enzyme::Trypsin,
            missed_cleavages: 2,
            fixed_modifications: vec![Modification {
                name: "Carbamidomethyl".into(),
                mass_delta: 57.021464,
                residues: vec!['C'],
                position: ModPosition::Anywhere,
            }],
            variable_modifications: vec![Modification {
                name: "Oxidation".into(),
                mass_delta: 15.994915,
                residues: vec!['M'],
                position: ModPosition::Anywhere,
            }],
            precursor_tolerance: MassTolerance {
                value: 20.0,
                unit: ToleranceUnit::Ppm,
            },
            fragment_tolerance: MassTolerance {
                value: 0.02,
                unit: ToleranceUnit::Da,
            },
            database_path: "/data/test.fasta".into(),
            decoy_strategy: DecoyStrategy::Reverse,
            acquisition_mode: None,
            max_variable_modifications: 3,
            min_peptide_length: 7,
            max_peptide_length: 50,
            engine: Some("Sage".into()),
        }
    }

    #[test]
    fn build_parameters_enzyme_trypsin() {
        let sage_params = build_sage_parameters(&test_params());
        let enzyme = sage_params.enzyme;
        assert_eq!(enzyme.cleave_at, Some("KR".into()));
        assert_eq!(enzyme.restrict, Some("P".into()));
        assert_eq!(enzyme.c_terminal, Some(true));
        assert_eq!(enzyme.missed_cleavages, Some(2));
    }

    #[test]
    fn build_parameters_peptide_length() {
        let sage_params = build_sage_parameters(&test_params());
        let enzyme = sage_params.enzyme;
        assert_eq!(enzyme.min_len, Some(7));
        assert_eq!(enzyme.max_len, Some(50));
    }

    #[test]
    fn build_parameters_static_mods() {
        let sage_params = build_sage_parameters(&test_params());
        assert!(sage_params
            .static_mods
            .contains_key(&ModificationSpecificity::Residue(b'C')));
        let mass = sage_params.static_mods[&ModificationSpecificity::Residue(b'C')];
        assert!((mass - 57.021464f32).abs() < 1e-4);
    }

    #[test]
    fn build_parameters_variable_mods() {
        let sage_params = build_sage_parameters(&test_params());
        assert!(sage_params
            .variable_mods
            .contains_key(&ModificationSpecificity::Residue(b'M')));
        let masses = &sage_params.variable_mods[&ModificationSpecificity::Residue(b'M')];
        assert_eq!(masses.len(), 1);
        assert!((masses[0] - 15.994915f32).abs() < 1e-4);
    }

    #[test]
    fn build_parameters_max_variable_mods() {
        let sage_params = build_sage_parameters(&test_params());
        assert_eq!(sage_params.max_variable_mods, 3);
    }

    #[test]
    fn build_parameters_decoy_generation() {
        let params = test_params();
        let sage_params = build_sage_parameters(&params);
        assert!(sage_params.generate_decoys);

        let mut no_decoy = params;
        no_decoy.decoy_strategy = DecoyStrategy::None;
        let sage_params2 = build_sage_parameters(&no_decoy);
        assert!(!sage_params2.generate_decoys);
    }

    #[test]
    fn build_parameters_fasta_path() {
        let sage_params = build_sage_parameters(&test_params());
        assert_eq!(sage_params.fasta, "/data/test.fasta");
    }
}
