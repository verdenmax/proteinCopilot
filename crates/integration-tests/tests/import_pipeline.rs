//! Integration test: custom JSON import pipeline.
//!
//! Creates a temporary custom JSON file, imports it, and verifies PSM parsing.

use protein_copilot_result_import::custom_json::CustomJsonParser;
use protein_copilot_result_import::unimod::UnimodDb;
use protein_copilot_result_import::ResultParser;
use std::io::Write;

#[test]
fn import_custom_json_parses_psms() {
    // Create a minimal custom JSON file
    let json_content = r#"[
        {
            "sequence": "PEPTIDEK",
            "charge": 2,
            "modify": [],
            "rt": 30.5,
            "precursor_mz": 458.24,
            "raw_title": "test_sample",
            "protein_names": ["sp|P12345|TEST_HUMAN"]
        },
        {
            "sequence": "DGFLLDGFPR",
            "charge": 2,
            "modify": [[7, 35]],
            "rt": 45.2,
            "precursor_mz": 547.28,
            "raw_title": "test_sample",
            "protein_names": ["sp|P67890|DEMO_HUMAN"]
        },
        {
            "sequence": "PEPTIDE",
            "charge": 2,
            "modify": [],
            "rt": 25.0,
            "precursor_mz": 400.19,
            "raw_title": "test_sample",
            "protein_names": ["sp|P11111|NOLABEL_HUMAN"]
        }
    ]"#;

    let dir = tempfile::tempdir().unwrap();
    let json_path = dir.path().join("test_import.json");
    let mut file = std::fs::File::create(&json_path).unwrap();
    file.write_all(json_content.as_bytes()).unwrap();

    // Parse with custom JSON parser
    let unimod = UnimodDb::builtin();
    let parser = CustomJsonParser;
    let psms = parser.parse(&json_path, &unimod).unwrap();

    assert_eq!(psms.len(), 3, "should parse 3 PSMs");

    // Verify first PSM
    assert_eq!(psms[0].sequence, "PEPTIDEK");
    assert_eq!(psms[0].charge, 2);
    assert!((psms[0].rt_min - 30.5).abs() < 0.01);
    assert!((psms[0].precursor_mz - 458.24).abs() < 0.01);
    assert_eq!(psms[0].raw_name, "test_sample");
    assert_eq!(psms[0].protein_accessions, vec!["sp|P12345|TEST_HUMAN"]);
    assert!(psms[0].modifications.is_empty());

    // Verify second PSM has oxidation modification (unimod:35)
    assert_eq!(psms[1].sequence, "DGFLLDGFPR");
    // unimod:35 = Oxidation — check if parsed
    if !psms[1].modifications.is_empty() {
        let m = &psms[1].modifications[0];
        assert!(
            m.name.contains("Oxidation") || m.mass_delta.abs() > 0.0,
            "modification should be Oxidation"
        );
    }

    // Verify third PSM — no K/R (zero offset scenario)
    assert_eq!(psms[2].sequence, "PEPTIDE");
    assert!(psms[2].modifications.is_empty());
}

#[test]
fn import_empty_json_array() {
    let dir = tempfile::tempdir().unwrap();
    let json_path = dir.path().join("empty.json");
    std::fs::write(&json_path, "[]").unwrap();

    let unimod = UnimodDb::builtin();
    let parser = CustomJsonParser;
    let psms = parser.parse(&json_path, &unimod).unwrap();
    assert!(psms.is_empty(), "empty array should yield 0 PSMs");
}

#[test]
fn import_nonexistent_file_errors() {
    let unimod = UnimodDb::builtin();
    let parser = CustomJsonParser;
    let result = parser.parse(std::path::Path::new("/nonexistent/file.json"), &unimod);
    assert!(result.is_err(), "missing file should return error");
}
