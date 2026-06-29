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
    let _span = tracing::info_span!("calculate_fdr", psm_count = psms.len()).entered();

    if psms.is_empty() {
        return Err(FdrError::NoPsms);
    }

    if psms.iter().any(|p| !p.score.is_finite()) {
        return Err(FdrError::InvalidScore);
    }

    // Check that at least one decoy hit exists for target-decoy FDR
    if !psms.iter().any(|p| p.is_decoy) {
        return Err(FdrError::NoDecoyHits);
    }

    // Sort by score descending
    let mut sorted: Vec<&ScoredPsm> = psms.iter().collect();
    sorted.sort_by(|a, b| b.score.total_cmp(&a.score));

    // Calculate FDR at each position
    let mut targets: u64 = 0;
    let mut decoys: u64 = 0;
    let mut raw_fdrs: Vec<f64> = Vec::with_capacity(sorted.len());

    let total = sorted.len();
    let progress_interval: usize = 5000;
    let loop_start = std::time::Instant::now();

    for (i, psm) in sorted.iter().enumerate() {
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

        if (i + 1) % progress_interval == 0 || i + 1 == total {
            let elapsed = loop_start.elapsed().as_secs_f64();
            let rate = if elapsed > 0.0 {
                (i + 1) as f64 / elapsed
            } else {
                0.0
            };
            let eta = if rate > 0.0 {
                (total - i - 1) as f64 / rate
            } else {
                0.0
            };
            tracing::info!(
                progress = i + 1,
                total = total,
                rate = format!("{:.0}/s", rate),
                eta_sec = format!("{:.1}", eta),
                "calculating FDR"
            );
        }
    }

    // Flatten tie groups: PSMs with identical scores must receive the SAME FDR.
    // For each maximal run of equal scores [lo, hi) in the descending-sorted array,
    // set every position to the FDR at the LAST position of the run (after counting
    // ALL targets+decoys at that score). This removes input-order dependence for ties.
    let mut lo = 0;
    while lo < sorted.len() {
        let mut hi = lo + 1;
        while hi < sorted.len()
            && sorted[hi].score.total_cmp(&sorted[lo].score) == std::cmp::Ordering::Equal
        {
            hi += 1;
        }
        let tie_fdr = raw_fdrs[hi - 1];
        for fdr in &mut raw_fdrs[lo..hi] {
            *fdr = tie_fdr;
        }
        lo = hi;
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

    let passing = result.iter().filter(|(_, q)| *q <= 0.01).count();
    tracing::info!(
        target = targets,
        decoy = decoys,
        psms_at_1pct = passing,
        "FDR calculated"
    );

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fdr_basic_target_decoy() {
        let psms = vec![
            ScoredPsm {
                index: 0,
                score: 0.9,
                is_decoy: false,
            },
            ScoredPsm {
                index: 1,
                score: 0.8,
                is_decoy: false,
            },
            ScoredPsm {
                index: 2,
                score: 0.7,
                is_decoy: false,
            },
            ScoredPsm {
                index: 3,
                score: 0.3,
                is_decoy: true,
            },
            ScoredPsm {
                index: 4,
                score: 0.2,
                is_decoy: true,
            },
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
            ScoredPsm {
                index: 0,
                score: 0.9,
                is_decoy: false,
            },
            ScoredPsm {
                index: 1,
                score: 0.8,
                is_decoy: true,
            },
            ScoredPsm {
                index: 2,
                score: 0.7,
                is_decoy: false,
            },
            ScoredPsm {
                index: 3,
                score: 0.6,
                is_decoy: false,
            },
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
            ScoredPsm {
                index: 0,
                score: 0.9,
                is_decoy: false,
            },
            ScoredPsm {
                index: 1,
                score: 0.8,
                is_decoy: false,
            },
        ];
        let result = calculate_fdr(&psms);
        assert!(
            matches!(result, Err(FdrError::NoDecoyHits)),
            "all targets, no decoys: should return NoDecoyHits error"
        );
    }

    #[test]
    fn fdr_tied_score_target_decoy_same_qvalue() {
        // Two PSMs with identical score, one target one decoy.
        // They MUST receive the same q-value regardless of which sorts first.
        let psms = vec![
            ScoredPsm {
                index: 0,
                score: 5.0,
                is_decoy: false,
            },
            ScoredPsm {
                index: 1,
                score: 5.0,
                is_decoy: true,
            },
        ];
        let result = calculate_fdr(&psms).unwrap();
        let q0 = result.iter().find(|(i, _)| *i == 0).unwrap().1;
        let q1 = result.iter().find(|(i, _)| *i == 1).unwrap().1;
        assert!(
            (q0 - q1).abs() < f64::EPSILON,
            "tied PSMs must share one q-value: target q={q0}, decoy q={q1}"
        );
    }

    #[test]
    fn fdr_order_independence_with_ties() {
        // Build a set containing a tied (target, decoy) pair at score 8.0.
        // q-values per original index must be identical regardless of input order.
        let psms = vec![
            ScoredPsm {
                index: 0,
                score: 10.0,
                is_decoy: false,
            },
            ScoredPsm {
                index: 1,
                score: 8.0,
                is_decoy: false,
            },
            ScoredPsm {
                index: 2,
                score: 8.0,
                is_decoy: true,
            },
            ScoredPsm {
                index: 3,
                score: 6.0,
                is_decoy: false,
            },
            ScoredPsm {
                index: 4,
                score: 5.0,
                is_decoy: true,
            },
        ];

        let q_of = |res: &[(usize, f64)]| -> std::collections::HashMap<usize, f64> {
            res.iter().copied().collect()
        };

        let forward = q_of(&calculate_fdr(&psms).unwrap());

        let mut reversed = psms.clone();
        reversed.reverse();
        let backward = q_of(&calculate_fdr(&reversed).unwrap());

        for idx in 0..psms.len() {
            assert!(
                (forward[&idx] - backward[&idx]).abs() < f64::EPSILON,
                "q-value for index {idx} must be order-independent: {} vs {}",
                forward[&idx],
                backward[&idx]
            );
        }
    }

    #[test]
    fn fdr_q_values_bounded() {
        let psms = vec![
            ScoredPsm {
                index: 0,
                score: 0.1,
                is_decoy: true,
            },
            ScoredPsm {
                index: 1,
                score: 0.2,
                is_decoy: true,
            },
            ScoredPsm {
                index: 2,
                score: 0.3,
                is_decoy: false,
            },
        ];
        let result = calculate_fdr(&psms).unwrap();
        for (_, q) in &result {
            assert!(*q >= 0.0 && *q <= 1.0, "q-value should be in [0,1]: {q}");
        }
    }
}
