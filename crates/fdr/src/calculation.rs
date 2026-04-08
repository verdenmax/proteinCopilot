//! Target-decoy FDR calculation and q-value assignment.
//!
//! Implements the standard target-decoy approach:
//! 1. Sort PSMs by score (descending)
//! 2. At each score threshold: FDR = #decoys / #targets
//! 3. q-value = minimum FDR at this score or better
//! 4. Enforce monotonicity (q-values never increase with better scores)

use crate::error::FdrError;

/// A PSM with score and target/decoy label, for FDR calculation.
#[derive(Debug, Clone)]
pub struct ScoredPsm {
    /// Index into the original PSM vector.
    pub index: usize,
    /// Search engine score (higher = better).
    pub score: f64,
    /// Whether this is a decoy hit.
    pub is_decoy: bool,
}

/// Calculates q-values for a list of scored PSMs using target-decoy approach.
///
/// Returns: `Vec<(usize, f64)>` — (original_index, q_value) for each PSM.
///
/// Algorithm:
/// 1. Sort by score descending
/// 2. Walk down: at each position, FDR = decoys_so_far / targets_so_far
/// 3. Walk back up: enforce monotonicity (q = min(current_fdr, next_q))
pub fn calculate_fdr(psms: &[ScoredPsm]) -> Result<Vec<(usize, f64)>, FdrError> {
    if psms.is_empty() {
        return Err(FdrError::NoPsms);
    }

    if psms.iter().any(|p| !p.score.is_finite()) {
        return Err(FdrError::InvalidScore);
    }

    // Sort by score descending
    let mut sorted: Vec<&ScoredPsm> = psms.iter().collect();
    sorted.sort_by(|a, b| b.score.total_cmp(&a.score));

    // Calculate FDR at each position
    let mut targets: u64 = 0;
    let mut decoys: u64 = 0;
    let mut raw_fdrs: Vec<f64> = Vec::with_capacity(sorted.len());

    for psm in &sorted {
        if psm.is_decoy {
            decoys += 1;
        } else {
            targets += 1;
        }
        let fdr = if targets > 0 {
            decoys as f64 / targets as f64
        } else {
            1.0
        };
        raw_fdrs.push(fdr.min(1.0));
    }

    // Enforce monotonicity: walk backward
    let mut q_values = raw_fdrs;
    for i in (0..q_values.len().saturating_sub(1)).rev() {
        q_values[i] = q_values[i].min(q_values[i + 1]);
    }

    let result: Vec<(usize, f64)> = sorted
        .iter()
        .zip(q_values.iter())
        .map(|(psm, &q)| (psm.index, q))
        .collect();

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fdr_basic_target_decoy() {
        let psms = vec![
            ScoredPsm { index: 0, score: 0.9, is_decoy: false },
            ScoredPsm { index: 1, score: 0.8, is_decoy: false },
            ScoredPsm { index: 2, score: 0.7, is_decoy: false },
            ScoredPsm { index: 3, score: 0.3, is_decoy: true },
            ScoredPsm { index: 4, score: 0.2, is_decoy: true },
        ];
        let result = calculate_fdr(&psms).unwrap();

        let q0 = result.iter().find(|(i, _)| *i == 0).unwrap().1;
        assert!(q0 < 0.01, "top target should have very low q-value: {q0}");

        let q3 = result.iter().find(|(i, _)| *i == 3).unwrap().1;
        assert!(q3 > 0.3, "decoy should have high q-value: {q3}");
    }

    #[test]
    fn fdr_monotonicity_enforced() {
        let psms = vec![
            ScoredPsm { index: 0, score: 0.9, is_decoy: false },
            ScoredPsm { index: 1, score: 0.8, is_decoy: true },
            ScoredPsm { index: 2, score: 0.7, is_decoy: false },
            ScoredPsm { index: 3, score: 0.6, is_decoy: false },
        ];
        let result = calculate_fdr(&psms).unwrap();

        // Verify monotonicity: q-values should not decrease when going to worse scores
        let mut sorted_by_score: Vec<(f64, f64)> = result
            .iter()
            .map(|(idx, q)| {
                let score = psms[*idx].score;
                (score, *q)
            })
            .collect();
        sorted_by_score.sort_by(|a, b| b.0.total_cmp(&a.0));

        for window in sorted_by_score.windows(2) {
            assert!(
                window[0].1 <= window[1].1 + 1e-12,
                "q-values should be monotonically non-decreasing: {:?}",
                sorted_by_score
            );
        }
    }

    #[test]
    fn fdr_empty_returns_error() {
        assert!(matches!(calculate_fdr(&[]), Err(FdrError::NoPsms)));
    }

    #[test]
    fn fdr_nan_score_returns_error() {
        let psms = vec![ScoredPsm {
            index: 0,
            score: f64::NAN,
            is_decoy: false,
        }];
        assert!(matches!(calculate_fdr(&psms), Err(FdrError::InvalidScore)));
    }

    #[test]
    fn fdr_all_targets_zero_fdr() {
        let psms = vec![
            ScoredPsm { index: 0, score: 0.9, is_decoy: false },
            ScoredPsm { index: 1, score: 0.8, is_decoy: false },
        ];
        let result = calculate_fdr(&psms).unwrap();
        for (_, q) in &result {
            assert!(q.abs() < 1e-9, "all targets, no decoys: FDR should be 0");
        }
    }

    #[test]
    fn fdr_q_values_bounded() {
        let psms = vec![
            ScoredPsm { index: 0, score: 0.1, is_decoy: true },
            ScoredPsm { index: 1, score: 0.2, is_decoy: true },
            ScoredPsm { index: 2, score: 0.3, is_decoy: false },
        ];
        let result = calculate_fdr(&psms).unwrap();
        for (_, q) in &result {
            assert!(*q >= 0.0 && *q <= 1.0, "q-value should be in [0,1]: {q}");
        }
    }
}
