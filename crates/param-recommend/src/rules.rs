//! Recommendation rules — the core logic of the parameter engine.
//!
//! This module contains the deterministic rules for:
//! - Instrument type inference (based on m/z range and peak density)
//! - Tolerance recommendation
//! - Enzyme recommendation
//! - Modification recommendation (based on experiment type)
//! - Confidence calculation
//! - Semantic conflict detection

use protein_copilot_core::ai_decision::AiDecision;
use protein_copilot_core::search_params::{SearchParams, ToleranceUnit};
use protein_copilot_core::spectrum::SpectrumSummary;

use crate::error::ParamRecommendError;
use crate::hints::UserHints;
use crate::preset;

/// Main recommendation entry point.
pub(crate) fn recommend(
    summary: &SpectrumSummary,
    hints: Option<&UserHints>,
) -> Result<AiDecision<SearchParams>, ParamRecommendError> {
    if summary.is_empty() {
        return Err(ParamRecommendError::EmptySummary);
    }

    // Validate summary fields (reject NaN/Inf values)
    summary
        .validate()
        .map_err(|e| ParamRecommendError::InvalidSummary {
            field: "summary",
            detail: e.to_string(),
        })?;

    // Step 1: Select base preset
    let experiment_type = hints
        .and_then(|h| h.experiment_type.as_deref())
        .unwrap_or("standard");

    let lower_type = experiment_type.to_lowercase();
    let is_open_search = lower_type.contains("open");

    let mut base = select_preset(&lower_type, is_open_search);

    // Step 2: Adjust tolerance based on instrument inference
    let instrument = infer_instrument(summary, hints);
    // Don't override tolerance for open search (uses Da-based tolerance)
    if !is_open_search {
        apply_tolerance(&mut base, &instrument);
    }

    // Step 3: Apply enzyme hint override
    if let Some(enzyme) = hints.and_then(|h| h.enzyme.clone()) {
        base.enzyme = enzyme;
    }

    // Step 3.5: Detect DIA acquisition mode from isolation window width
    let dia_detected = detect_dia(summary);
    if dia_detected {
        base.acquisition_mode = Some(protein_copilot_core::spectrum::AcquisitionMode::DIA);
    }

    // Step 4: Build explanation, evidence, alternatives
    let explanation = build_explanation(summary, &instrument, experiment_type, dia_detected);
    let evidence = build_evidence(summary, &instrument, dia_detected);
    let alternatives = build_alternatives(experiment_type);
    let confidence = compute_confidence(hints, &instrument);

    // Step 5: Detect semantic conflicts and append warnings
    let warnings = detect_conflicts(&base, is_open_search);
    let final_explanation = if warnings.is_empty() {
        explanation
    } else {
        format!("{explanation}\n\n⚠ Warnings:\n{}", warnings.join("\n"))
    };

    // Step 6: Build input_summary
    let input_summary = format!(
        "{} spectra, m/z range [{:.0}-{:.0}], median {} peaks/spectrum, RT [{:.0}-{:.0}] sec",
        summary.total_spectra,
        summary.mz_range[0],
        summary.mz_range[1],
        summary.median_peaks_per_spectrum,
        summary.rt_range_sec[0],
        summary.rt_range_sec[1],
    );

    Ok(AiDecision {
        decision: base,
        confidence,
        explanation: final_explanation,
        input_summary,
        alternatives,
        evidence,
    })
}

// ---------------------------------------------------------------------------
// DIA detection
// ---------------------------------------------------------------------------

/// Detect DIA acquisition mode based on median isolation window width.
///
/// Returns `true` if `median_isolation_window_da > 5.0 Da`.
fn detect_dia(summary: &SpectrumSummary) -> bool {
    summary
        .median_isolation_window_da
        .map(|w| w > 5.0)
        .unwrap_or(false)
}

// ---------------------------------------------------------------------------
// Instrument inference
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
enum InstrumentClass {
    HighResolution, // Orbitrap, Q-Exactive
    LowResolution,  // TOF, ion trap
    General,        // cannot determine
}

fn infer_instrument(summary: &SpectrumSummary, hints: Option<&UserHints>) -> InstrumentClass {
    // User hint overrides auto-inference
    if let Some(hint) = hints.and_then(|h| h.instrument_type.as_deref()) {
        let lower = hint.to_lowercase();
        if lower.contains("orbitrap") || lower.contains("exactive") || lower.contains("hires") {
            return InstrumentClass::HighResolution;
        }
        if lower.contains("tof") || lower.contains("trap") || lower.contains("lowres") {
            return InstrumentClass::LowResolution;
        }
    }

    // Scoring-based inference from data characteristics
    let mz_upper = summary.mz_range[1];
    let median_peaks = summary.median_peaks_per_spectrum;
    let mut hi_score: i32 = 0;
    let mut lo_score: i32 = 0;

    if mz_upper > 1800.0 {
        hi_score += 2;
    } else if mz_upper > 1500.0 {
        hi_score += 1;
    } else if mz_upper < 1200.0 {
        lo_score += 1;
    }

    if median_peaks > 300 {
        hi_score += 2;
    } else if median_peaks > 200 {
        hi_score += 1;
    } else if median_peaks < 100 {
        lo_score += 1;
    }

    if hi_score >= 2 {
        InstrumentClass::HighResolution
    } else if lo_score >= 2 {
        InstrumentClass::LowResolution
    } else {
        InstrumentClass::General
    }
}

// ---------------------------------------------------------------------------
// Preset selection
// ---------------------------------------------------------------------------

fn select_preset(lower_type: &str, is_open_search: bool) -> SearchParams {
    // Open search takes priority: it defines the tolerance strategy,
    // then we layer experiment-specific modifications on top.
    if is_open_search {
        let mut params = preset::open_search_preset().params;
        // Layer experiment-specific modifications onto open search base
        if lower_type.contains("phospho") {
            merge_modifications(&mut params, &preset::phospho_preset().params);
        } else if lower_type.contains("tmt") {
            merge_modifications(&mut params, &preset::tmt_preset().params);
        } else if lower_type.contains("silac") {
            merge_modifications(&mut params, &preset::silac_preset().params);
        }
        params
    } else if lower_type.contains("phospho") {
        preset::phospho_preset().params
    } else if lower_type.contains("tmt") {
        preset::tmt_preset().params
    } else if lower_type.contains("silac") {
        preset::silac_preset().params
    } else {
        preset::standard_preset().params
    }
}

/// Merge fixed and variable modifications from `source` into `target`,
/// avoiding duplicates (by modification name).
fn merge_modifications(target: &mut SearchParams, source: &SearchParams) {
    for m in &source.fixed_modifications {
        if !target.fixed_modifications.iter().any(|t| t.name == m.name) {
            target.fixed_modifications.push(m.clone());
        }
    }
    for m in &source.variable_modifications {
        if !target
            .variable_modifications
            .iter()
            .any(|t| t.name == m.name)
        {
            target.variable_modifications.push(m.clone());
        }
    }
}

// ---------------------------------------------------------------------------
// Tolerance adjustment
// ---------------------------------------------------------------------------

fn apply_tolerance(params: &mut SearchParams, instrument: &InstrumentClass) {
    use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};

    match instrument {
        InstrumentClass::HighResolution => {
            params.precursor_tolerance = MassTolerance {
                value: 10.0,
                unit: ToleranceUnit::Ppm,
            };
            params.fragment_tolerance = MassTolerance {
                value: 20.0,
                unit: ToleranceUnit::Ppm,
            };
        }
        InstrumentClass::LowResolution => {
            params.precursor_tolerance = MassTolerance {
                value: 20.0,
                unit: ToleranceUnit::Ppm,
            };
            params.fragment_tolerance = MassTolerance {
                value: 0.1,
                unit: ToleranceUnit::Da,
            };
        }
        InstrumentClass::General => {
            params.precursor_tolerance = MassTolerance {
                value: 15.0,
                unit: ToleranceUnit::Ppm,
            };
            params.fragment_tolerance = MassTolerance {
                value: 20.0,
                unit: ToleranceUnit::Ppm,
            };
        }
    }
}

// ---------------------------------------------------------------------------
// Explanation, evidence, alternatives
// ---------------------------------------------------------------------------

fn build_explanation(
    summary: &SpectrumSummary,
    instrument: &InstrumentClass,
    experiment_type: &str,
    dia_detected: bool,
) -> String {
    let instrument_desc = match instrument {
        InstrumentClass::HighResolution => "high-resolution instrument (e.g., Orbitrap)",
        InstrumentClass::LowResolution => "low-resolution instrument (e.g., TOF/ion trap)",
        InstrumentClass::General => "instrument with moderate resolution",
    };

    let tolerance_desc = match instrument {
        InstrumentClass::HighResolution => "precursor tolerance 10 ppm, fragment tolerance 0.02 Da",
        InstrumentClass::LowResolution => "precursor tolerance 20 ppm, fragment tolerance 0.1 Da",
        InstrumentClass::General => "precursor tolerance 15 ppm, fragment tolerance 0.05 Da",
    };

    let mut text = format!(
        "Based on m/z range [{:.0}-{:.0}] and median {:.0} peaks/spectrum, \
         inferred {instrument_desc}. Recommending {tolerance_desc}. \
         Experiment type: \"{experiment_type}\", using Trypsin digestion with \
         appropriate modification set.",
        summary.mz_range[0], summary.mz_range[1], summary.median_peaks_per_spectrum,
    );

    if dia_detected {
        text.push_str(
            " DIA acquisition mode detected (median isolation window > 5 Da). \
             Use extract_dia_precursors tool before running search.",
        );
    }

    text
}

fn build_evidence(summary: &SpectrumSummary, instrument: &InstrumentClass, dia_detected: bool) -> Vec<String> {
    let mut evidence = vec![
        format!(
            "m/z range: {:.0}-{:.0}",
            summary.mz_range[0], summary.mz_range[1]
        ),
        format!(
            "Median peaks per spectrum: {}",
            summary.median_peaks_per_spectrum
        ),
        format!(
            "Total spectra: {} ({} MS1, {} MS2)",
            summary.total_spectra, summary.ms1_count, summary.ms2_count
        ),
    ];

    if !summary.precursor_charge_distribution.is_empty() {
        let mut charges: Vec<_> = summary.precursor_charge_distribution.iter().collect();
        charges.sort_by_key(|(c, _)| *c);
        let dist: Vec<String> = charges.iter().map(|(c, n)| format!("{c}+: {n}")).collect();
        evidence.push(format!("Charge distribution: {}", dist.join(", ")));
    }

    evidence.push(format!("Instrument inference: {:?}", instrument));

    if dia_detected {
        if let Some(w) = summary.median_isolation_window_da {
            evidence.push(format!("DIA detected: median isolation window {w:.1} Da"));
        }
    }

    evidence
}

fn build_alternatives(experiment_type: &str) -> Vec<String> {
    let lower = experiment_type.to_lowercase();
    let mut alts = Vec::new();

    if !lower.contains("open") {
        alts.push("Open search with 500 Da tolerance for PTM discovery".to_string());
    }
    if !lower.contains("phospho") {
        alts.push("Phosphoproteomics search with Phospho(STY) modification".to_string());
    }
    if !lower.contains("tmt") {
        alts.push("TMT-labeled search with TMT6plex modifications".to_string());
    }
    if !lower.contains("silac") {
        alts.push("SILAC labeled search with heavy K/R modifications".to_string());
    }
    if lower != "standard" && !lower.contains("standard") {
        alts.push("Standard search without special modifications".to_string());
    }

    alts
}

// ---------------------------------------------------------------------------
// Confidence
// ---------------------------------------------------------------------------

fn compute_confidence(hints: Option<&UserHints>, instrument: &InstrumentClass) -> f64 {
    let mut confidence: f64 = match instrument {
        InstrumentClass::HighResolution | InstrumentClass::LowResolution => 0.80,
        InstrumentClass::General => 0.70,
    };

    if let Some(h) = hints {
        if h.experiment_type.is_some() {
            confidence += 0.10;
        }
        if h.instrument_type.is_some() {
            confidence += 0.10;
        }
        if h.enzyme.is_some() {
            confidence += 0.05;
        }
    }

    confidence.min(0.95)
}

// ---------------------------------------------------------------------------
// Semantic conflict detection
// ---------------------------------------------------------------------------

fn detect_conflicts(params: &SearchParams, is_open_search: bool) -> Vec<String> {
    use protein_copilot_core::search_params::Enzyme;

    let mut warnings = Vec::new();

    // NonSpecific enzyme + missed_cleavages > 0 is contradictory
    if params.enzyme == Enzyme::NonSpecific && params.missed_cleavages > 0 {
        warnings.push(format!(
            "- NonSpecific enzyme with missed_cleavages={} is contradictory \
             (NonSpecific has no cleavage sites to miss)",
            params.missed_cleavages
        ));
    }

    // Open search + narrow fragment tolerance (skip if using open preset intentionally)
    if !is_open_search && params.precursor_tolerance.value > 100.0 {
        let narrow_frag = match params.fragment_tolerance.unit {
            ToleranceUnit::Da => params.fragment_tolerance.value < 0.05,
            ToleranceUnit::Ppm => params.fragment_tolerance.value < 10.0,
        };
        if narrow_frag {
            warnings.push(
                "- Wide precursor tolerance with very narrow fragment tolerance \
                 may produce few matches; consider widening fragment tolerance"
                    .to_string(),
            );
        }
    }

    // Open search + high missed cleavages → huge search space
    if params.precursor_tolerance.value > 100.0 && params.missed_cleavages > 3 {
        warnings.push(format!(
            "- Open search with missed_cleavages={} will dramatically increase search \
             space and runtime; consider reducing to ≤ 2",
            params.missed_cleavages
        ));
    }

    warnings
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::search_params::ToleranceUnit;
    use protein_copilot_core::spectrum::{SpectrumFormat, SpectrumSummary};
    use std::collections::HashMap;

    fn high_res_summary() -> SpectrumSummary {
        let mut charge_dist = HashMap::new();
        charge_dist.insert(2, 5000);
        charge_dist.insert(3, 3000);

        SpectrumSummary {
            file_path: "/data/sample.mzML".to_string(),
            format: SpectrumFormat::MzML,
            total_spectra: 10000,
            ms1_count: 500,
            ms2_count: 9500,
            mz_range: [100.0, 2000.0],
            rt_range_sec: [0.0, 3600.0],
            precursor_charge_distribution: charge_dist,
            median_peaks_per_spectrum: 256,
            median_isolation_window_da: None,
        }
    }

    fn low_res_summary() -> SpectrumSummary {
        let mut charge_dist = HashMap::new();
        charge_dist.insert(2, 3000);

        SpectrumSummary {
            file_path: "/data/low_res.mgf".to_string(),
            format: SpectrumFormat::Mgf,
            total_spectra: 5000,
            ms1_count: 0,
            ms2_count: 5000,
            mz_range: [100.0, 1100.0],
            rt_range_sec: [0.0, 1800.0],
            precursor_charge_distribution: charge_dist,
            median_peaks_per_spectrum: 50,
            median_isolation_window_da: None,
        }
    }

    fn empty_summary() -> SpectrumSummary {
        SpectrumSummary {
            file_path: "/data/empty.mgf".to_string(),
            format: SpectrumFormat::Mgf,
            total_spectra: 0,
            ms1_count: 0,
            ms2_count: 0,
            mz_range: [0.0, 0.0],
            rt_range_sec: [0.0, 0.0],
            precursor_charge_distribution: HashMap::new(),
            median_peaks_per_spectrum: 0,
            median_isolation_window_da: None,
        }
    }

    // -- Instrument inference -------------------------------------------

    #[test]
    fn high_res_data_gets_10ppm() {
        let result = recommend(&high_res_summary(), None).unwrap();
        assert_eq!(result.decision.precursor_tolerance.value, 10.0);
        assert_eq!(result.decision.fragment_tolerance.value, 20.0);
        assert_eq!(result.decision.fragment_tolerance.unit, ToleranceUnit::Ppm);
    }

    #[test]
    fn low_res_data_gets_20ppm() {
        let result = recommend(&low_res_summary(), None).unwrap();
        assert_eq!(result.decision.precursor_tolerance.value, 20.0);
        assert_eq!(result.decision.fragment_tolerance.value, 0.1);
        assert_eq!(result.decision.fragment_tolerance.unit, ToleranceUnit::Da);
    }

    #[test]
    fn instrument_hint_overrides_inference() {
        let hints = UserHints {
            instrument_type: Some("Orbitrap".to_string()),
            ..Default::default()
        };
        let result = recommend(&low_res_summary(), Some(&hints)).unwrap();
        // Low-res data, but Orbitrap hint → high-res tolerance
        assert_eq!(result.decision.precursor_tolerance.value, 10.0);
    }

    // -- Experiment type ------------------------------------------------

    #[test]
    fn phospho_hint_adds_phospho_mod() {
        let hints = UserHints {
            experiment_type: Some("phosphorylation".to_string()),
            ..Default::default()
        };
        let result = recommend(&high_res_summary(), Some(&hints)).unwrap();
        assert!(result
            .decision
            .variable_modifications
            .iter()
            .any(|m| m.name == "Phospho"));
    }

    #[test]
    fn tmt_hint_adds_tmt_fixed() {
        let hints = UserHints {
            experiment_type: Some("TMT".to_string()),
            ..Default::default()
        };
        let result = recommend(&high_res_summary(), Some(&hints)).unwrap();
        assert!(result
            .decision
            .fixed_modifications
            .iter()
            .any(|m| m.name == "TMT6plex"));
    }

    #[test]
    fn default_uses_trypsin() {
        let result = recommend(&high_res_summary(), None).unwrap();
        assert_eq!(
            result.decision.enzyme,
            protein_copilot_core::search_params::Enzyme::Trypsin
        );
    }

    // -- Confidence -----------------------------------------------------

    #[test]
    fn confidence_higher_with_hints() {
        let no_hints = recommend(&high_res_summary(), None).unwrap();
        let with_hints = recommend(
            &high_res_summary(),
            Some(&UserHints {
                experiment_type: Some("standard".to_string()),
                instrument_type: Some("Orbitrap".to_string()),
                ..Default::default()
            }),
        )
        .unwrap();

        assert!(with_hints.confidence > no_hints.confidence);
        assert!(with_hints.confidence <= 0.95);
    }

    #[test]
    fn confidence_in_valid_range() {
        let result = recommend(&high_res_summary(), None).unwrap();
        assert!(result.confidence >= 0.0 && result.confidence <= 1.0);
    }

    // -- AiDecision fields ----------------------------------------------

    #[test]
    fn result_has_non_empty_explanation() {
        let result = recommend(&high_res_summary(), None).unwrap();
        assert!(!result.explanation.is_empty());
        assert!(result.explanation.contains("m/z range"));
    }

    #[test]
    fn result_has_evidence() {
        let result = recommend(&high_res_summary(), None).unwrap();
        assert!(!result.evidence.is_empty());
    }

    #[test]
    fn result_has_alternatives() {
        let result = recommend(&high_res_summary(), None).unwrap();
        assert!(!result.alternatives.is_empty());
    }

    #[test]
    fn result_has_input_summary() {
        let result = recommend(&high_res_summary(), None).unwrap();
        assert!(result.input_summary.contains("10000 spectra"));
    }

    // -- Determinism ----------------------------------------------------

    #[test]
    fn deterministic_same_input_same_output() {
        let s = high_res_summary();
        let r1 = recommend(&s, None).unwrap();
        let r2 = recommend(&s, None).unwrap();
        assert_eq!(r1.decision, r2.decision);
        assert_eq!(r1.confidence, r2.confidence);
        assert_eq!(r1.explanation, r2.explanation);
    }

    // -- Empty file -----------------------------------------------------

    #[test]
    fn empty_summary_returns_error() {
        let result = recommend(&empty_summary(), None);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("empty"));
    }

    // -- Semantic conflicts ---------------------------------------------

    #[test]
    fn nonspecific_with_missed_cleavages_warns() {
        use protein_copilot_core::search_params::Enzyme;
        let mut params = preset::standard_preset().params;
        params.enzyme = Enzyme::NonSpecific;
        params.missed_cleavages = 2;
        let warnings = detect_conflicts(&params, false);
        assert!(!warnings.is_empty());
        assert!(warnings[0].contains("NonSpecific"));
    }

    #[test]
    fn open_search_preset_no_self_conflict() {
        let result = recommend(
            &high_res_summary(),
            Some(&UserHints {
                experiment_type: Some("open".to_string()),
                ..Default::default()
            }),
        )
        .unwrap();
        // Open search preset should NOT warn about its own narrow fragment tolerance
        assert!(
            !result.explanation.contains("⚠ Warnings"),
            "Open search preset should not trigger self-conflict"
        );
    }

    // -- Compound experiment types --------------------------------------

    #[test]
    fn tmt_open_gets_open_tolerance_with_tmt_mods() {
        let hints = UserHints {
            experiment_type: Some("tmt open".to_string()),
            ..Default::default()
        };
        let result = recommend(&high_res_summary(), Some(&hints)).unwrap();
        // Must use open search tolerance (500 Da), not TMT's 10 ppm
        assert_eq!(
            result.decision.precursor_tolerance.value, 500.0,
            "Compound 'tmt open' should use 500 Da precursor tolerance"
        );
        // Must also include TMT modifications
        assert!(
            result
                .decision
                .fixed_modifications
                .iter()
                .any(|m| m.name.contains("TMT")),
            "Compound 'tmt open' should include TMT modifications"
        );
    }

    // -- Enzyme override ------------------------------------------------

    #[test]
    fn enzyme_hint_overrides_default() {
        use protein_copilot_core::search_params::Enzyme;
        let hints = UserHints {
            enzyme: Some(Enzyme::LysC),
            ..Default::default()
        };
        let result = recommend(&high_res_summary(), Some(&hints)).unwrap();
        assert_eq!(result.decision.enzyme, Enzyme::LysC);
    }

    // -- Middle resolution (General) ------------------------------------

    #[test]
    fn mid_res_data_gets_15ppm() {
        let mut summary = high_res_summary();
        summary.mz_range = [100.0, 1600.0]; // borderline
        summary.median_peaks_per_spectrum = 150; // between thresholds
        let result = recommend(&summary, None).unwrap();
        assert_eq!(result.decision.precursor_tolerance.value, 15.0);
        assert_eq!(result.decision.fragment_tolerance.value, 20.0);
        assert_eq!(result.decision.fragment_tolerance.unit, ToleranceUnit::Ppm);
    }

    // -- Open search selection ------------------------------------------

    #[test]
    fn open_hint_selects_open_preset() {
        let hints = UserHints {
            experiment_type: Some("open".to_string()),
            ..Default::default()
        };
        let result = recommend(&high_res_summary(), Some(&hints)).unwrap();
        assert_eq!(
            result.decision.precursor_tolerance.value, 500.0,
            "Open search should keep 500 Da tolerance"
        );
    }

    // -- Confidence meets spec ------------------------------------------

    #[test]
    fn single_hint_reaches_090() {
        let hints = UserHints {
            experiment_type: Some("standard".to_string()),
            ..Default::default()
        };
        let result = recommend(&high_res_summary(), Some(&hints)).unwrap();
        assert!(
            result.confidence >= 0.90,
            "Single hint should reach 0.90+, got {}",
            result.confidence
        );
    }

    // -- Validation passthrough ----------------------------------------

    #[test]
    fn recommended_params_validate() {
        let result = recommend(&high_res_summary(), None).unwrap();
        let mut params = result.decision;
        params.database_path = "/data/test.fasta".to_string();
        assert!(params.validate().is_ok());
    }

    // -- with_database helper ------------------------------------------

    #[test]
    fn preset_with_database_validates() {
        let params = preset::standard_preset().with_database("/data/human.fasta");
        assert!(params.validate().is_ok());
    }

    // -- DIA detection --------------------------------------------------

    #[test]
    fn dia_detected_from_wide_isolation_window() {
        let mut summary = high_res_summary();
        summary.median_isolation_window_da = Some(25.0);
        let result = recommend(&summary, None).unwrap();
        assert_eq!(
            result.decision.acquisition_mode,
            Some(protein_copilot_core::spectrum::AcquisitionMode::DIA),
        );
        assert!(result.explanation.contains("DIA acquisition mode detected"));
        assert!(result.evidence.iter().any(|e| e.contains("DIA detected")));
    }

    #[test]
    fn dda_not_flagged_with_narrow_window() {
        let mut summary = high_res_summary();
        summary.median_isolation_window_da = Some(2.0);
        let result = recommend(&summary, None).unwrap();
        assert_eq!(result.decision.acquisition_mode, None);
        assert!(!result.explanation.contains("DIA"));
    }

    #[test]
    fn no_isolation_window_stays_none() {
        let summary = high_res_summary(); // median_isolation_window_da is None
        let result = recommend(&summary, None).unwrap();
        assert_eq!(result.decision.acquisition_mode, None);
        assert!(!result.explanation.contains("DIA"));
    }
}
