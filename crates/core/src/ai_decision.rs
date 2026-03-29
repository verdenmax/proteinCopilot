//! Structured AI decision output wrapper.
//!
//! This module defines [`AiDecision<T>`], a generic wrapper that enforces
//! structured output for all AI/LLM-assisted decisions in ProteinCopilot.
//!
//! Every recommendation, interpretation, or diagnostic produced by the AI
//! layer must include confidence, explanation, evidence, and alternatives.
//! This ensures auditability and reproducibility per §2.4 of
//! copilot-instructions.md.
//!
//! The `decision` field is generic — it can wrap a [`SearchParams`] for
//! parameter recommendations, a `String` for free-form interpretations,
//! or any domain-specific type.
//!
//! [`SearchParams`]: crate::search_params::SearchParams

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use thiserror::Error;

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Errors related to AI decision validation.
#[derive(Debug, Error)]
pub enum AiDecisionError {
    /// Confidence value is outside the valid range [0.0, 1.0].
    #[error("confidence must be in [0.0, 1.0], got {value}")]
    InvalidConfidence {
        /// The actual confidence value.
        value: f64,
    },

    /// A required text field is empty or whitespace-only.
    #[error("{field} must not be empty")]
    EmptyField {
        /// Name of the field.
        field: &'static str,
    },
}

// ---------------------------------------------------------------------------
// AiDecision<T>
// ---------------------------------------------------------------------------

/// Structured wrapper for all AI-assisted decisions.
///
/// Every time the LLM (via Agent/Skill) makes a recommendation or
/// interpretation, the result is wrapped in this struct. This ensures:
///
/// - **Auditability**: `explanation` + `evidence` record *why* a decision
///   was made.
/// - **Calibration**: `confidence` gives a quantitative self-assessment.
/// - **Exploration**: `alternatives` lists other options considered.
/// - **Context**: `input_summary` captures the data the decision was based on.
///
/// # Example JSON output (per copilot-instructions.md §2.4)
///
/// ```json
/// {
///   "decision": "推荐使用 Trypsin 作为消化酶",
///   "confidence": 0.92,
///   "explanation": "输入数据的末端碎裂模式符合 Trypsin 消化特征...",
///   "input_summary": "检测到 12,345 张谱图，平均母离子质量 1,200 Da...",
///   "alternatives": ["Lys-C", "Chymotrypsin"],
///   "evidence": ["末端碎裂模式分析", "母离子质量分布"]
/// }
/// ```
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct AiDecision<T> {
    /// The concrete decision or recommendation value.
    pub decision: T,
    /// Confidence score in \[0.0, 1.0\] (1.0 = fully confident).
    pub confidence: f64,
    /// Human-readable explanation of the reasoning behind this decision.
    pub explanation: String,
    /// Summary of the input data on which this decision was based.
    pub input_summary: String,
    /// Other options that were considered but not chosen.
    pub alternatives: Vec<String>,
    /// Evidence or observations supporting this decision.
    pub evidence: Vec<String>,
}

impl<T> AiDecision<T> {
    /// Validates the decision wrapper fields.
    ///
    /// Checks:
    /// - `confidence` is finite and in \[0.0, 1.0\]
    /// - `explanation` is not empty
    /// - `input_summary` is not empty
    pub fn validate(&self) -> Result<(), AiDecisionError> {
        if !self.confidence.is_finite() || !(0.0..=1.0).contains(&self.confidence) {
            return Err(AiDecisionError::InvalidConfidence {
                value: self.confidence,
            });
        }
        if self.explanation.trim().is_empty() {
            return Err(AiDecisionError::EmptyField {
                field: "explanation",
            });
        }
        if self.input_summary.trim().is_empty() {
            return Err(AiDecisionError::EmptyField {
                field: "input_summary",
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
    use crate::search_params::{
        DecoyStrategy, Enzyme, MassTolerance, Modification, SearchParams, ToleranceUnit,
    };

    fn sample_string_decision() -> AiDecision<String> {
        AiDecision {
            decision: "推荐使用 Trypsin 作为消化酶".to_string(),
            confidence: 0.92,
            explanation: "输入数据的末端碎裂模式符合 Trypsin 消化特征".to_string(),
            input_summary: "检测到 12,345 张谱图，平均母离子质量 1,200 Da".to_string(),
            alternatives: vec!["Lys-C".to_string(), "Chymotrypsin".to_string()],
            evidence: vec!["末端碎裂模式分析".to_string(), "母离子质量分布".to_string()],
        }
    }

    fn sample_params_decision() -> AiDecision<SearchParams> {
        AiDecision {
            decision: SearchParams {
                database_path: "/data/uniprot_human.fasta".to_string(),
                enzyme: Enzyme::Trypsin,
                missed_cleavages: 2,
                fixed_modifications: vec![Modification {
                    name: "Carbamidomethyl".to_string(),
                    mass_delta: 57.021464,
                    residues: vec!['C'],
                    position: crate::search_params::ModPosition::Anywhere,
                }],
                variable_modifications: vec![],
                precursor_tolerance: MassTolerance {
                    value: 20.0,
                    unit: ToleranceUnit::Ppm,
                },
                fragment_tolerance: MassTolerance {
                    value: 0.02,
                    unit: ToleranceUnit::Da,
                },
                decoy_strategy: DecoyStrategy::Reverse,
            },
            confidence: 0.85,
            explanation: "Based on HeLa cell line data characteristics".to_string(),
            input_summary: "12,345 MS2 spectra, median precursor m/z 650".to_string(),
            alternatives: vec!["Open search with wider tolerance".to_string()],
            evidence: vec![
                "Precursor mass distribution analysis".to_string(),
                "Fragmentation pattern consistent with HCD".to_string(),
            ],
        }
    }

    // -- Serde roundtrip ------------------------------------------------

    #[test]
    fn string_decision_serde_roundtrip() {
        let d = sample_string_decision();
        let json = serde_json::to_string_pretty(&d).unwrap();
        let back: AiDecision<String> = serde_json::from_str(&json).unwrap();
        assert_eq!(d.decision, back.decision);
        assert_eq!(d.confidence, back.confidence);
        assert_eq!(d.alternatives.len(), back.alternatives.len());
        assert_eq!(d.evidence.len(), back.evidence.len());
    }

    #[test]
    fn params_decision_serde_roundtrip() {
        let d = sample_params_decision();
        let json = serde_json::to_string_pretty(&d).unwrap();
        let back: AiDecision<SearchParams> = serde_json::from_str(&json).unwrap();
        assert_eq!(d.decision.database_path, back.decision.database_path);
        assert_eq!(d.confidence, back.confidence);
        assert_eq!(d.explanation, back.explanation);
    }

    #[test]
    fn json_output_has_required_fields_string() {
        // Verify compliance with copilot-instructions.md §2.4
        let d = sample_string_decision();
        let json = serde_json::to_string(&d).unwrap();
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v.get("decision").is_some());
        assert!(v.get("confidence").is_some());
        assert!(v.get("explanation").is_some());
        assert!(v.get("input_summary").is_some());
        assert!(v.get("alternatives").is_some());
        assert!(v.get("evidence").is_some());
    }

    #[test]
    fn json_output_has_required_fields_search_params() {
        // Sub-task 1.1.5.2: AiDecision<SearchParams> JSON format per §2.4
        let d = sample_params_decision();
        let json = serde_json::to_string(&d).unwrap();
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v.get("decision").is_some());
        assert!(
            v["decision"].is_object(),
            "decision should be a SearchParams object"
        );
        assert!(v.get("confidence").is_some());
        assert!(v["confidence"].is_f64());
        assert!(v.get("explanation").is_some());
        assert!(v.get("input_summary").is_some());
        assert!(v.get("alternatives").is_some());
        assert!(v["alternatives"].is_array());
        assert!(v.get("evidence").is_some());
        assert!(v["evidence"].is_array());
    }

    // -- Validation -----------------------------------------------------

    #[test]
    fn validate_passes_for_valid_data() {
        assert!(sample_string_decision().validate().is_ok());
        assert!(sample_params_decision().validate().is_ok());
    }

    #[test]
    fn validate_accepts_boundary_confidence() {
        let mut d = sample_string_decision();
        d.confidence = 0.0;
        assert!(d.validate().is_ok());
        d.confidence = 1.0;
        assert!(d.validate().is_ok());
    }

    #[test]
    fn validate_rejects_confidence_above_one() {
        let mut d = sample_string_decision();
        d.confidence = 1.01;
        let err = d.validate().unwrap_err();
        assert!(err.to_string().contains("confidence"));
    }

    #[test]
    fn validate_rejects_negative_confidence() {
        let mut d = sample_string_decision();
        d.confidence = -0.1;
        let err = d.validate().unwrap_err();
        assert!(err.to_string().contains("confidence"));
    }

    #[test]
    fn validate_rejects_nan_confidence() {
        let mut d = sample_string_decision();
        d.confidence = f64::NAN;
        assert!(d.validate().is_err());
    }

    #[test]
    fn validate_rejects_empty_explanation() {
        let mut d = sample_string_decision();
        d.explanation = "".to_string();
        let err = d.validate().unwrap_err();
        assert!(err.to_string().contains("explanation"));
    }

    #[test]
    fn validate_rejects_whitespace_only_explanation() {
        let mut d = sample_string_decision();
        d.explanation = "   ".to_string();
        let err = d.validate().unwrap_err();
        assert!(err.to_string().contains("explanation"));
    }

    #[test]
    fn validate_rejects_empty_input_summary() {
        let mut d = sample_string_decision();
        d.input_summary = "".to_string();
        let err = d.validate().unwrap_err();
        assert!(err.to_string().contains("input_summary"));
    }

    #[test]
    fn validate_accepts_empty_alternatives_and_evidence() {
        let mut d = sample_string_decision();
        d.alternatives = vec![];
        d.evidence = vec![];
        assert!(d.validate().is_ok());
    }
}
