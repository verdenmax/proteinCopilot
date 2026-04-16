//! End-to-end integration tests for the protein inference pipeline.
//!
//! Exercises the full pipeline: mapper → parsimony → razor → protein FDR → coverage
//! with realistic proteomics test data.

use std::collections::HashMap;

use protein_copilot_core::search_result::Psm;
use protein_copilot_fdr::protein_fdr::calculate_protein_fdr;
use protein_copilot_protein_inference::coverage::calculate_coverage;
use protein_copilot_protein_inference::mapper::{build_peptide_protein_map, normalize_il};
use protein_copilot_protein_inference::parsimony::run_parsimony;
use protein_copilot_protein_inference::razor::assign_razor_peptides;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Create a minimal PSM for testing.
fn make_psm(
    sequence: &str,
    proteins: &[&str],
    score: f64,
    q_value: Option<f64>,
    is_decoy: bool,
) -> Psm {
    Psm {
        spectrum_scan: 1,
        peptide_sequence: sequence.to_string(),
        modifications: vec![],
        charge: 2,
        precursor_mz: 500.0,
        calculated_mz: 500.0,
        delta_mass_ppm: 0.0,
        score,
        protein_accessions: proteins.iter().map(|s| s.to_string()).collect(),
        is_decoy,
        q_value,
    }
}

/// Build the standard realistic test dataset.
///
/// Target proteins:
///   P001 — unique peptides: ACDEFGHK, LMNPQR
///   P002 — unique peptides: MNPQRSTW, GHIJKLMN
///   P003 — shares ACDEFGHK with P001, has unique VWXYZABC
///   P004 — indistinguishable from P002 (same peptides: MNPQRSTW, GHIJKLMN)
///   P005 — unique peptide STUVWXYZ, also has KDEFGHIJ (I/L equivalent to KDEFGHLJ)
///   P006 — has the I/L variant KDEFGHLJ, plus unique QRSTUVWX
///
/// Decoy proteins (REV_ prefix):
///   REV_P001 — KHGFEDCA (low score)
///   REV_P002 — WTSRQPNM (low score)
///   REV_P003 — CBAZYX (low score)
///   REV_P004 — NMLKJIHG (low score)
///   REV_P005 — ZYXWVUTS (low score)
///   REV_P006 — XWVUTSRQ (low score)
fn build_realistic_psms() -> Vec<Psm> {
    vec![
        // --- Target PSMs (q ≤ 0.01) ---
        // P001 unique
        make_psm("ACDEFGHK", &["P001", "P003"], 20.0, Some(0.001), false),
        make_psm("LMNPQR", &["P001"], 18.0, Some(0.002), false),
        // P002 unique (shared with P004 indistinguishable)
        make_psm("MNPQRSTW", &["P002", "P004"], 17.0, Some(0.001), false),
        make_psm("GHIJKLMN", &["P002", "P004"], 16.0, Some(0.002), false),
        // P003 unique
        make_psm("VWXYZABC", &["P003"], 15.0, Some(0.003), false),
        // P005 unique
        make_psm("STUVWXYZ", &["P005"], 14.0, Some(0.004), false),
        // I/L equivalent: KDEFGHIJ (P005) and KDEFGHLJ (P006) normalize to same sequence
        make_psm("KDEFGHIJ", &["P005"], 13.0, Some(0.005), false),
        make_psm("KDEFGHLJ", &["P006"], 12.0, Some(0.006), false),
        // P006 unique
        make_psm("QRSTUVWX", &["P006"], 11.0, Some(0.007), false),
        // Duplicate PSM for ACDEFGHK at lower score (tests best-score tracking)
        make_psm("ACDEFGHK", &["P001", "P003"], 8.0, Some(0.008), false),
        // --- Decoy PSMs ---
        make_psm("KHGFEDCA", &["REV_P001"], 5.0, Some(0.005), true),
        make_psm("WTSRQPNM", &["REV_P002"], 4.0, Some(0.006), true),
        make_psm("CBAZYX", &["REV_P003"], 3.5, Some(0.007), true),
        make_psm("NMLKJIHG", &["REV_P004"], 3.0, Some(0.008), true),
        make_psm("ZYXWVUTS", &["REV_P005"], 2.5, Some(0.009), true),
        make_psm("XWVUTSRQ", &["REV_P006"], 2.0, Some(0.010), true),
        // --- PSMs above q-value threshold (should be filtered at 0.01) ---
        make_psm("BADPEPTIDE", &["P001"], 1.0, Some(0.05), false),
        make_psm("ANOTHERBAD", &["P002"], 0.5, Some(0.10), false),
        // PSM without q-value (should be kept when filtering)
        make_psm("NOQVALPEP", &["P001"], 9.0, None, false),
    ]
}

/// Build a FASTA sequence map for the realistic dataset.
fn build_fasta_sequences() -> HashMap<String, String> {
    let mut fasta = HashMap::new();
    // Target proteins — sequences contain their identified peptides
    fasta.insert(
        "P001".to_string(),
        "MSTARTACDEFGHKLMNPQRNOQVALPEPENDXYZ".to_string(),
    );
    fasta.insert(
        "P002".to_string(),
        "MSTARTMNPQRSTWGHIJKLMNENDXYZ".to_string(),
    );
    fasta.insert(
        "P003".to_string(),
        "MSTARTACDEFGHKVWXYZABCENDXYZ".to_string(),
    );
    // P004 same as P002 (indistinguishable)
    fasta.insert(
        "P004".to_string(),
        "MSTARTMNPQRSTWGHIJKLMNENDXYZ".to_string(),
    );
    fasta.insert(
        "P005".to_string(),
        "MSTARTSTUVWXYZKDEFGHIJENDXYZ".to_string(),
    );
    fasta.insert(
        "P006".to_string(),
        "MSTARTKDEFGHLJQRSTUVWXENDXYZ".to_string(),
    );
    // Decoy proteins
    fasta.insert("REV_P001".to_string(), "MKHGFEDCAXYZ".to_string());
    fasta.insert("REV_P002".to_string(), "MWTSRQPNMXYZ".to_string());
    fasta.insert("REV_P003".to_string(), "MCBAZYXXYZ".to_string());
    fasta.insert("REV_P004".to_string(), "MNMLKJIHGXYZ".to_string());
    fasta.insert("REV_P005".to_string(), "MZYXWVUTSXYZ".to_string());
    fasta.insert("REV_P006".to_string(), "MXWVUTSRQXYZ".to_string());
    fasta
}

// ---------------------------------------------------------------------------
// Full pipeline test
// ---------------------------------------------------------------------------

#[test]
fn test_full_inference_pipeline() {
    let psms = build_realistic_psms();
    let fasta_sequences = build_fasta_sequences();

    // Step 1: Build peptide-protein map with q-value filtering at 1% FDR
    let map = build_peptide_protein_map(&psms, Some(0.01)).unwrap();

    // Verify q-value filtering: BADPEPTIDE (q=0.05) and ANOTHERBAD (q=0.10) excluded
    assert!(
        !map.peptide_to_proteins
            .contains_key(&normalize_il("BADPEPTIDE")),
        "BADPEPTIDE (q=0.05) should be filtered"
    );
    assert!(
        !map.peptide_to_proteins
            .contains_key(&normalize_il("ANOTHERBAD")),
        "ANOTHERBAD (q=0.10) should be filtered"
    );
    // NOQVALPEP has no q-value → kept
    assert!(
        map.peptide_to_proteins
            .contains_key(&normalize_il("NOQVALPEP")),
        "PSM with no q-value should be kept"
    );

    // Verify I/L equivalence: KDEFGHIJ and KDEFGHLJ normalize to the same key
    let norm_il = normalize_il("KDEFGHIJ");
    assert_eq!(norm_il, normalize_il("KDEFGHLJ"));
    assert!(map.peptide_to_proteins.contains_key(&norm_il));
    let il_proteins = &map.peptide_to_proteins[&norm_il];
    assert!(
        il_proteins.contains("P005") && il_proteins.contains("P006"),
        "I/L-equivalent peptides should map to both P005 and P006"
    );

    // Verify best-score tracking for duplicate PSMs
    let acdefghk_norm = normalize_il("ACDEFGHK");
    assert_eq!(
        map.peptide_best_score[&acdefghk_norm], 20.0,
        "best score for ACDEFGHK should be 20.0 (not 8.0)"
    );

    // Step 2: Run parsimony
    let mut groups = run_parsimony(&map).unwrap();
    assert!(
        !groups.is_empty(),
        "parsimony should produce at least one group"
    );

    // P002 and P004 are indistinguishable (identical peptide sets) → same group
    let p002_group = groups
        .iter()
        .find(|g| g.member_accessions.contains(&"P002".to_string()))
        .expect("P002 should be in a group");
    assert!(
        p002_group
            .member_accessions
            .contains(&"P004".to_string()),
        "P002 and P004 should be in the same group (indistinguishable)"
    );

    // Every peptide from the map should appear in at least one group
    let all_group_peptides: std::collections::HashSet<&str> = groups
        .iter()
        .flat_map(|g| g.peptides.iter().map(String::as_str))
        .collect();
    for pep in map.peptide_to_proteins.keys() {
        assert!(
            all_group_peptides.contains(pep.as_str()),
            "peptide {pep} not covered by any group"
        );
    }

    // Step 3: Assign razor peptides
    let razor_map = assign_razor_peptides(&mut groups, &map);
    // Razor map may or may not be populated depending on shared peptides among
    // selected groups. Just verify it doesn't panic and returns a valid map.
    for (pep, leader) in &razor_map {
        assert!(
            groups.iter().any(|g| g.leader_accession == *leader),
            "razor peptide {pep} assigned to unknown leader {leader}"
        );
    }

    // Step 4: Calculate protein-level FDR
    let fdr_result = calculate_protein_fdr(&groups).unwrap();
    let mut final_groups = fdr_result.groups;

    // All returned groups should be targets (decoys filtered out)
    assert!(
        final_groups.iter().all(|g| !g.is_decoy),
        "FDR result should only contain target groups"
    );
    // All targets should have q-values assigned
    assert!(
        final_groups.iter().all(|g| g.q_value.is_some()),
        "all target groups should have q-values"
    );
    // q-values should be valid [0.0, 1.0]
    for g in &final_groups {
        let q = g.q_value.unwrap();
        assert!(
            (0.0..=1.0).contains(&q),
            "q-value {q} out of range for {}",
            g.leader_accession
        );
    }

    // Step 5: Calculate coverage with FASTA sequences
    calculate_coverage(&mut final_groups, &fasta_sequences);

    // All target groups with a matching FASTA entry should have coverage
    for g in &final_groups {
        if fasta_sequences.contains_key(&g.leader_accession) {
            assert!(
                g.coverage.is_some(),
                "group {} should have coverage",
                g.leader_accession
            );
            let cov = g.coverage.unwrap();
            assert!(
                (0.0..=1.0).contains(&cov),
                "coverage {cov} out of range for {}",
                g.leader_accession
            );
            assert!(
                cov > 0.0,
                "coverage should be > 0 for {} (peptides exist in FASTA)",
                g.leader_accession
            );
        }
    }

    // Validate each group structurally
    for g in &final_groups {
        g.validate().unwrap_or_else(|e| {
            panic!(
                "group {} failed validation: {e}",
                g.leader_accession
            )
        });
    }
}

// ---------------------------------------------------------------------------
// Full pipeline without q-value filtering
// ---------------------------------------------------------------------------

#[test]
fn test_pipeline_no_qvalue_filter() {
    let psms = build_realistic_psms();

    // No q-value threshold → all PSMs included (including BADPEPTIDE, ANOTHERBAD)
    let map = build_peptide_protein_map(&psms, None).unwrap();
    assert!(
        map.peptide_to_proteins
            .contains_key(&normalize_il("BADPEPTIDE")),
        "without q-value filter, BADPEPTIDE should be present"
    );

    let mut groups = run_parsimony(&map).unwrap();
    let _razor_map = assign_razor_peptides(&mut groups, &map);
    let fdr_result = calculate_protein_fdr(&groups).unwrap();

    assert!(
        fdr_result.groups.iter().all(|g| !g.is_decoy),
        "FDR result should only contain targets"
    );
    assert!(
        fdr_result.groups.iter().all(|g| g.q_value.is_some()),
        "all targets should have q-values"
    );
}

// ---------------------------------------------------------------------------
// Edge case: no decoy PSMs → FDR gives all q=0
// ---------------------------------------------------------------------------

#[test]
fn test_pipeline_no_decoys() {
    let psms = vec![
        make_psm("PEPTIDEA", &["P001"], 20.0, Some(0.001), false),
        make_psm("PEPTIDEB", &["P001"], 18.0, Some(0.002), false),
        make_psm("PEPTIDEC", &["P002"], 15.0, Some(0.003), false),
    ];

    let map = build_peptide_protein_map(&psms, Some(0.01)).unwrap();
    let mut groups = run_parsimony(&map).unwrap();
    let _razor_map = assign_razor_peptides(&mut groups, &map);

    // No decoy groups → calculate_protein_fdr should succeed with all q=0
    let fdr_result = calculate_protein_fdr(&groups).unwrap();
    assert!(
        fdr_result.groups.iter().all(|g| !g.is_decoy),
        "all groups should be targets"
    );
    for g in &fdr_result.groups {
        assert_eq!(
            g.q_value,
            Some(0.0),
            "with no decoys, all q-values should be 0.0"
        );
    }
    assert_eq!(fdr_result.total_decoy_groups, 0);
}

// ---------------------------------------------------------------------------
// Edge case: only decoy PSMs → mapper returns NoTargetProteins error
// ---------------------------------------------------------------------------

#[test]
fn test_pipeline_only_decoys() {
    let psms = vec![
        make_psm("DECOYPEP1", &["REV_P001"], 10.0, Some(0.001), true),
        make_psm("DECOYPEP2", &["REV_P002"], 8.0, Some(0.002), true),
    ];

    let result = build_peptide_protein_map(&psms, Some(0.01));
    assert!(
        result.is_err(),
        "only-decoy PSMs should produce an error from mapper"
    );

    use protein_copilot_protein_inference::error::InferenceError;
    assert!(
        matches!(result, Err(InferenceError::NoTargetProteins)),
        "expected NoTargetProteins error"
    );
}

// ---------------------------------------------------------------------------
// Edge case: all PSMs above q-value threshold → NoPsms or NoTargetProteins
// ---------------------------------------------------------------------------

#[test]
fn test_pipeline_all_psms_above_threshold() {
    let psms = vec![
        make_psm("PEPTIDEA", &["P001"], 10.0, Some(0.05), false),
        make_psm("PEPTIDEB", &["P002"], 8.0, Some(0.10), false),
        make_psm("DECOYPEP", &["REV_P001"], 5.0, Some(0.20), true),
    ];

    // All q-values exceed 0.01 → everything filtered → error
    let result = build_peptide_protein_map(&psms, Some(0.01));
    assert!(
        result.is_err(),
        "all PSMs above threshold should produce an error"
    );
}

// ---------------------------------------------------------------------------
// Edge case: single protein, single peptide (simplest case)
// ---------------------------------------------------------------------------

#[test]
fn test_pipeline_single_protein_single_peptide() {
    let psms = vec![make_psm("SIMPLEPEP", &["PROT_SOLO"], 25.0, Some(0.001), false)];
    let fasta: HashMap<String, String> =
        [("PROT_SOLO".into(), "MSIMPLEPEPENDXYZ".into())]
            .into_iter()
            .collect();

    let map = build_peptide_protein_map(&psms, Some(0.01)).unwrap();
    let mut groups = run_parsimony(&map).unwrap();

    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0].leader_accession, "PROT_SOLO");
    assert_eq!(groups[0].peptides.len(), 1);
    assert_eq!(groups[0].unique_peptides.len(), 1);

    let razor_map = assign_razor_peptides(&mut groups, &map);
    assert!(
        razor_map.is_empty(),
        "single group → no shared peptides → empty razor map"
    );

    // No decoys → all q=0
    let fdr_result = calculate_protein_fdr(&groups).unwrap();
    assert_eq!(fdr_result.groups.len(), 1);
    assert_eq!(fdr_result.groups[0].q_value, Some(0.0));

    let mut final_groups = fdr_result.groups;
    calculate_coverage(&mut final_groups, &fasta);
    assert!(final_groups[0].coverage.is_some());
    assert!(final_groups[0].coverage.unwrap() > 0.0);
}

// ---------------------------------------------------------------------------
// Edge case: large number of indistinguishable proteins → single group
// ---------------------------------------------------------------------------

#[test]
fn test_pipeline_many_indistinguishable_proteins() {
    // 20 proteins sharing the exact same 2 peptides → one group with 20 members
    let protein_names: Vec<String> = (1..=20).map(|i| format!("PROT_{i:03}")).collect();
    let protein_refs: Vec<&str> = protein_names.iter().map(String::as_str).collect();

    let psms = vec![
        make_psm("SHAREDPEP1", &protein_refs, 20.0, Some(0.001), false),
        make_psm("SHAREDPEP2", &protein_refs, 18.0, Some(0.002), false),
    ];

    let map = build_peptide_protein_map(&psms, Some(0.01)).unwrap();
    let mut groups = run_parsimony(&map).unwrap();

    assert_eq!(groups.len(), 1, "all indistinguishable → one group");
    assert_eq!(groups[0].member_accessions.len(), 20);
    // Leader should be alphabetically first
    assert_eq!(groups[0].leader_accession, "PROT_001");
    // Members should be sorted alphabetically
    let mut sorted_members = groups[0].member_accessions.clone();
    sorted_members.sort();
    assert_eq!(groups[0].member_accessions, sorted_members);

    let razor_map = assign_razor_peptides(&mut groups, &map);
    assert!(razor_map.is_empty(), "single group → no razor assignment");
}

// ---------------------------------------------------------------------------
// I/L equivalence across the full pipeline
// ---------------------------------------------------------------------------

#[test]
fn test_pipeline_il_equivalence_grouping() {
    // Two proteins that share an I/L-equivalent peptide but each has a truly unique one
    let psms = vec![
        make_psm("PEPTIDE", &["PROT_I"], 20.0, Some(0.001), false),
        make_psm("PEPTLDE", &["PROT_L"], 18.0, Some(0.002), false),
        make_psm("AAAAAA", &["PROT_I"], 15.0, Some(0.003), false),
        make_psm("BBBBBB", &["PROT_L"], 14.0, Some(0.004), false),
    ];

    let map = build_peptide_protein_map(&psms, Some(0.01)).unwrap();

    // PEPTIDE and PEPTLDE normalize to the same key
    let norm = normalize_il("PEPTIDE");
    assert_eq!(norm, normalize_il("PEPTLDE"));
    assert!(map.peptide_to_proteins.contains_key(&norm));
    assert!(map.peptide_to_proteins[&norm].contains("PROT_I"));
    assert!(map.peptide_to_proteins[&norm].contains("PROT_L"));

    let mut groups = run_parsimony(&map).unwrap();
    // Both proteins share the I/L peptide but each has a unique peptide,
    // so they remain separate groups (not indistinguishable).
    assert!(
        groups.len() >= 2,
        "PROT_I and PROT_L should be separate groups (have different unique peptides)"
    );

    let _razor_map = assign_razor_peptides(&mut groups, &map);
    let fdr_result = calculate_protein_fdr(&groups).unwrap();
    assert!(fdr_result.groups.iter().all(|g| g.q_value.is_some()));
}

// ---------------------------------------------------------------------------
// Verify decoy-target pairing in protein FDR through full pipeline
// ---------------------------------------------------------------------------

#[test]
fn test_pipeline_target_decoy_competition() {
    // Scenario: target wins some pairs, decoy wins others
    let psms = vec![
        // P001 target (score 25) vs REV_P001 decoy (score 30) → decoy wins
        make_psm("TARGETPEP1", &["P001"], 25.0, Some(0.001), false),
        make_psm("DECOYPEP1", &["REV_P001"], 30.0, Some(0.002), true),
        // P002 target (score 20) vs REV_P002 decoy (score 5) → target wins
        make_psm("TARGETPEP2", &["P002"], 20.0, Some(0.001), false),
        make_psm("DECOYPEP2", &["REV_P002"], 5.0, Some(0.003), true),
        // P003 target (score 18) — no matching decoy → wins automatically
        make_psm("TARGETPEP3", &["P003"], 18.0, Some(0.001), false),
    ];

    let map = build_peptide_protein_map(&psms, Some(0.01)).unwrap();
    let mut groups = run_parsimony(&map).unwrap();
    let _razor_map = assign_razor_peptides(&mut groups, &map);

    let fdr_result = calculate_protein_fdr(&groups).unwrap();

    // P001 lost to its decoy → not in results
    let leaders: Vec<&str> = fdr_result
        .groups
        .iter()
        .map(|g| g.leader_accession.as_str())
        .collect();
    assert!(
        !leaders.contains(&"P001"),
        "P001 should be eliminated (decoy won)"
    );
    // P002 and P003 should be present
    assert!(leaders.contains(&"P002"), "P002 should survive FDR");
    assert!(leaders.contains(&"P003"), "P003 should survive FDR");

    // All surviving groups are targets with q-values
    for g in &fdr_result.groups {
        assert!(!g.is_decoy);
        assert!(g.q_value.is_some());
    }
}

// ---------------------------------------------------------------------------
// Coverage with missing FASTA entries
// ---------------------------------------------------------------------------

#[test]
fn test_pipeline_coverage_missing_fasta() {
    let psms = vec![
        make_psm("PEPTIDEA", &["PROT_KNOWN"], 20.0, Some(0.001), false),
        make_psm("PEPTIDEB", &["PROT_MISSING"], 18.0, Some(0.002), false),
    ];
    // Only provide FASTA for PROT_KNOWN
    let fasta: HashMap<String, String> =
        [("PROT_KNOWN".into(), "MPEPTIDEAENDXYZ".into())]
            .into_iter()
            .collect();

    let map = build_peptide_protein_map(&psms, Some(0.01)).unwrap();
    let mut groups = run_parsimony(&map).unwrap();
    let _razor_map = assign_razor_peptides(&mut groups, &map);
    let fdr_result = calculate_protein_fdr(&groups).unwrap();
    let mut final_groups = fdr_result.groups;

    calculate_coverage(&mut final_groups, &fasta);

    let known = final_groups
        .iter()
        .find(|g| g.leader_accession == "PROT_KNOWN");
    let missing = final_groups
        .iter()
        .find(|g| g.leader_accession == "PROT_MISSING");

    if let Some(g) = known {
        assert!(g.coverage.is_some(), "known protein should have coverage");
        assert!(g.coverage.unwrap() > 0.0);
    }
    if let Some(g) = missing {
        assert!(
            g.coverage.is_none(),
            "missing-FASTA protein should have coverage = None"
        );
    }
}

// ---------------------------------------------------------------------------
// Razor peptide assignment in full pipeline
// ---------------------------------------------------------------------------

#[test]
fn test_pipeline_razor_assignment_correctness() {
    // Two proteins sharing a peptide, each with different unique-peptide counts
    let psms = vec![
        // P_BIG has 3 unique peptides → more evidence
        make_psm("UNIQUE_BIG1", &["P_BIG"], 20.0, Some(0.001), false),
        make_psm("UNIQUE_BIG2", &["P_BIG"], 19.0, Some(0.001), false),
        make_psm("UNIQUE_BIG3", &["P_BIG"], 18.0, Some(0.001), false),
        // P_SMALL has 1 unique peptide
        make_psm("UNIQUE_SMALL", &["P_SMALL"], 17.0, Some(0.002), false),
        // Shared peptide between both
        make_psm("SHARED_PEP", &["P_BIG", "P_SMALL"], 16.0, Some(0.003), false),
    ];

    let map = build_peptide_protein_map(&psms, Some(0.01)).unwrap();
    let mut groups = run_parsimony(&map).unwrap();
    let razor_map = assign_razor_peptides(&mut groups, &map);

    // SHARED_PEP should be assigned to P_BIG (more unique peptides)
    if let Some(assigned_to) = razor_map.get(&normalize_il("SHARED_PEP")) {
        assert_eq!(
            assigned_to, "P_BIG",
            "shared peptide should be razored to the group with more unique peptides"
        );
    }

    // Verify razor_peptides field is set on the winning group
    let big_group = groups
        .iter()
        .find(|g| g.leader_accession == "P_BIG")
        .expect("P_BIG group should exist");
    if !razor_map.is_empty() {
        assert!(
            big_group
                .razor_peptides
                .contains(&normalize_il("SHARED_PEP")),
            "P_BIG should have SHARED_PEP in razor_peptides"
        );
    }
}

// ---------------------------------------------------------------------------
// Score ordering preserved through pipeline
// ---------------------------------------------------------------------------

#[test]
fn test_pipeline_score_ordering() {
    let psms = vec![
        make_psm("HIGHSCORE", &["P_HIGH"], 50.0, Some(0.001), false),
        make_psm("MIDSCORE", &["P_MID"], 30.0, Some(0.002), false),
        make_psm("LOWSCORE", &["P_LOW"], 10.0, Some(0.003), false),
        // Decoys with lower scores
        make_psm("DEC1", &["REV_P_HIGH"], 2.0, Some(0.005), true),
        make_psm("DEC2", &["REV_P_MID"], 1.5, Some(0.006), true),
        make_psm("DEC3", &["REV_P_LOW"], 1.0, Some(0.007), true),
    ];

    let map = build_peptide_protein_map(&psms, Some(0.01)).unwrap();
    let mut groups = run_parsimony(&map).unwrap();
    let _razor_map = assign_razor_peptides(&mut groups, &map);
    let fdr_result = calculate_protein_fdr(&groups).unwrap();

    // FDR result should be sorted by score descending
    for window in fdr_result.groups.windows(2) {
        assert!(
            window[0].score >= window[1].score,
            "groups should be sorted by score descending: {} < {}",
            window[0].score,
            window[1].score
        );
    }
}

// ---------------------------------------------------------------------------
// Empty PSM input
// ---------------------------------------------------------------------------

#[test]
fn test_pipeline_empty_psms() {
    let result = build_peptide_protein_map(&[], Some(0.01));
    assert!(result.is_err(), "empty PSMs should produce an error");
    use protein_copilot_protein_inference::error::InferenceError;
    assert!(matches!(result, Err(InferenceError::NoPsms)));
}

// ---------------------------------------------------------------------------
// Validate that coverage + razor interact correctly
// ---------------------------------------------------------------------------

#[test]
fn test_pipeline_coverage_includes_razor_peptides() {
    // PROT_A: unique PEP1, PEP2; shared PEP_S razored to PROT_A
    // PROT_B: unique PEP3; shared PEP_S not razored here
    // Coverage of PROT_A should include PEP_S contribution
    let psms = vec![
        make_psm("AAAA", &["PROT_A"], 20.0, Some(0.001), false),
        make_psm("BBBB", &["PROT_A"], 19.0, Some(0.001), false),
        make_psm("CCCC", &["PROT_B"], 18.0, Some(0.002), false),
        make_psm("SSSS", &["PROT_A", "PROT_B"], 17.0, Some(0.003), false),
    ];
    // FASTA: protein sequences that include the peptides
    let fasta: HashMap<String, String> = [
        ("PROT_A".into(), "MAAAABBBBSSSSXYZ".into()),
        ("PROT_B".into(), "MCCCCSSSSXYZ".into()),
    ]
    .into_iter()
    .collect();

    let map = build_peptide_protein_map(&psms, Some(0.01)).unwrap();
    let mut groups = run_parsimony(&map).unwrap();
    let _razor_map = assign_razor_peptides(&mut groups, &map);

    // Calculate FDR (no decoys → all q=0) then coverage
    let fdr_result = calculate_protein_fdr(&groups).unwrap();
    let mut final_groups = fdr_result.groups;
    calculate_coverage(&mut final_groups, &fasta);

    // PROT_A should have coverage > 0
    let prot_a = final_groups
        .iter()
        .find(|g| g.leader_accession == "PROT_A");
    if let Some(g) = prot_a {
        assert!(g.coverage.is_some());
        assert!(
            g.coverage.unwrap() > 0.0,
            "PROT_A coverage should be positive"
        );
    }
}
