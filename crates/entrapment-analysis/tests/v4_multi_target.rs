//! v4 multi-target provenance integration tests.

use protein_copilot_entrapment_analysis::coelution::{CoElutionIndex, DiaWindow};
use protein_copilot_entrapment_analysis::multi_provenance::trace_multi_target;
use protein_copilot_entrapment_analysis::types::*;

use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};

#[test]
fn test_full_pipeline_mock() {
    let targets: Vec<(UnifiedPsm, PsmGroup)> = vec![(
        UnifiedPsm {
            peptide: "STTSGHLVYK".to_string(),
            charge: Some(2),
            precursor_mz: Some(548.12),
            retention_time: Some(35.15),
            rt_start: Some(34.5),
            rt_stop: Some(35.8),
            scan_number: None,
            spectrum_file: Some("Rep1".to_string()),
            protein_ids: "sp|P12345|EF1A_HUMAN".to_string(),
            q_value: Some(0.001),
            modifications: vec![],
        },
        PsmGroup::Target,
    )];

    let psms: Vec<UnifiedPsm> = targets.iter().map(|(p, _)| p.clone()).collect();
    let groups: Vec<PsmGroup> = targets.iter().map(|(_, g)| *g).collect();

    let windows = vec![DiaWindow {
        center: 548.0,
        low: 546.0,
        high: 550.0,
    }];
    let index = CoElutionIndex::build(&psms, &groups, &windows, None, 20);

    let trap = UnifiedPsm {
        peptide: "STTTGHLIYK".to_string(),
        charge: Some(2),
        precursor_mz: Some(548.30),
        retention_time: Some(35.2),
        rt_start: Some(34.8),
        rt_stop: Some(35.6),
        scan_number: Some(12345),
        spectrum_file: Some("Rep1".to_string()),
        protein_ids: "sp|P99999|TRAP_YEAST".to_string(),
        q_value: Some(0.005),
        modifications: vec![],
    };

    let candidates = index.find_co_eluting(&trap, "Rep1");
    assert!(!candidates.is_empty(), "should find co-eluting targets");
    assert_eq!(candidates[0].peptide, "STTSGHLVYK");

    let tolerance = MassTolerance {
        value: 20.0,
        unit: ToleranceUnit::Ppm,
    };
    let observed_mz = vec![200.0, 300.0, 400.0, 500.0, 600.0, 700.0];
    let observed_int = vec![1000.0; 6];

    let result = trace_multi_target(
        &observed_mz,
        &observed_int,
        &trap.peptide,
        &trap.modifications,
        &candidates,
        &tolerance,
        2,
    );

    assert_eq!(result.annotated_peaks.len(), 6);
}
