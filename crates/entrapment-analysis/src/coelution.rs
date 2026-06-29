//! Co-elution index for efficient lookup of target PSMs co-eluting with a
//! given trap PSM within the same DIA isolation window.
//!
//! [`CoElutionIndex`] groups target PSMs by run name (spectrum file) and
//! sorts them by `rt_start` so that binary search + forward scan yields all
//! overlapping entries in O(log n + k) time.

use std::collections::HashMap;

use crate::config::SilacConfig;
use crate::types::{CoElutingCandidate, LabelForm, PsmGroup, UnifiedPsm};

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// A DIA isolation window described by its center, lower, and upper m/z bounds.
#[derive(Debug, Clone, Copy)]
pub struct DiaWindow {
    /// Center m/z of the window.
    pub center: f64,
    /// Lower m/z bound (inclusive).
    pub low: f64,
    /// Upper m/z bound (inclusive).
    pub high: f64,
}

// ---------------------------------------------------------------------------
// Internal types
// ---------------------------------------------------------------------------

/// A pre-processed target entry stored in the per-run index.
struct TargetEntry {
    peptide: String,
    protein_ids: Vec<String>,
    precursor_mz: f64,
    charge: i32,
    rt_start: f64,
    rt_stop: f64,
    modifications: Vec<(usize, f64)>,
}

// ---------------------------------------------------------------------------
// CoElutionIndex
// ---------------------------------------------------------------------------

/// Index of target PSMs organised by run name for fast co-elution queries.
///
/// Each run's entries are sorted by `rt_start` so that
/// [`find_co_eluting`](Self::find_co_eluting) can use binary search
/// (`partition_point`) followed by a bounded forward scan.
pub struct CoElutionIndex {
    by_run: HashMap<String, Vec<TargetEntry>>,
    dia_windows: Vec<DiaWindow>,
    silac: Option<SilacConfig>,
    max_candidates: usize,
}

impl CoElutionIndex {
    /// Build a new [`CoElutionIndex`] from classified PSMs.
    ///
    /// Only `Target` PSMs with valid `rt_start`, `rt_stop`, `precursor_mz`,
    /// and `spectrum_file` are included.  Each run's entries are sorted by
    /// `rt_start` to enable binary search in queries.
    pub fn build(
        psms: &[UnifiedPsm],
        groups: &[PsmGroup],
        dia_windows: &[DiaWindow],
        silac: Option<&SilacConfig>,
        max_candidates: usize,
    ) -> Self {
        let mut by_run: HashMap<String, Vec<TargetEntry>> = HashMap::new();

        for (psm, group) in psms.iter().zip(groups.iter()) {
            if *group != PsmGroup::Target {
                continue;
            }

            // Require all essential fields.
            let (rt_start, rt_stop, precursor_mz, spectrum_file) = match (
                psm.rt_start,
                psm.rt_stop,
                psm.precursor_mz,
                psm.spectrum_file.as_deref(),
            ) {
                (Some(s), Some(e), Some(mz), Some(f)) => (s, e, mz, f),
                _ => continue,
            };

            let charge = psm.charge.unwrap_or(2);

            let entry = TargetEntry {
                peptide: psm.peptide.clone(),
                protein_ids: psm
                    .protein_ids
                    .split(';')
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty())
                    .collect(),
                precursor_mz,
                charge,
                rt_start,
                rt_stop,
                modifications: psm.modifications.clone(),
            };

            by_run
                .entry(spectrum_file.to_string())
                .or_default()
                .push(entry);
        }

        // Sort each run's entries by rt_start for binary search.
        for entries in by_run.values_mut() {
            entries.sort_by(|a, b| {
                a.rt_start
                    .partial_cmp(&b.rt_start)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        }

        Self {
            by_run,
            dia_windows: dia_windows.to_vec(),
            silac: silac.cloned(),
            max_candidates,
        }
    }

    /// Find all target PSMs co-eluting with a given trap PSM in the specified
    /// run.
    ///
    /// A target co-elutes when:
    /// 1. Its RT window overlaps with the trap's RT window.
    /// 2. Both the trap and target precursor m/z fall within the same DIA
    ///    isolation window.
    /// 3. The target peptide is not identical to the trap peptide.
    ///
    /// If SILAC is configured and `enable_heavy_search` is `true`, a heavy
    /// candidate is also emitted for each qualifying light target.
    ///
    /// Results are capped at `max_candidates`.
    pub fn find_co_eluting(&self, trap: &UnifiedPsm, run: &str) -> Vec<CoElutingCandidate> {
        let entries = match self.by_run.get(run) {
            Some(e) => e,
            None => return Vec::new(),
        };

        let trap_rt_start = match trap.rt_start {
            Some(v) => v,
            None => return Vec::new(),
        };
        let trap_rt_stop = match trap.rt_stop {
            Some(v) => v,
            None => return Vec::new(),
        };
        let trap_mz = match trap.precursor_mz {
            Some(v) => v,
            None => return Vec::new(),
        };

        let trap_window = self.find_dia_window(trap_mz);

        let mut candidates = Vec::new();

        // Entries are sorted by `rt_start`, so only that key may be binary
        // searched. A target overlaps the trap iff
        //   entry.rt_start <= trap_rt_stop   AND   entry.rt_stop >= trap_rt_start.
        //
        // `partition_point` requires a predicate that is monotonic over the
        // slice. Only `rt_start` is sorted, so we use it for the UPPER bound
        // (first entry whose rt_start > trap_rt_stop). The lower bound depends on
        // `rt_stop`, which is NOT sorted and therefore must be filtered linearly.
        let end_idx = entries.partition_point(|e| e.rt_start <= trap_rt_stop);

        for entry in entries[..end_idx].iter() {
            if entry.rt_stop < trap_rt_start {
                continue; // ends before the trap starts → no RT overlap
            }

            // RT overlap confirmed: entry.rt_start <= trap_rt_stop (by end_idx)
            // and entry.rt_stop >= trap_rt_start (by the guard above).

            // Skip same peptide.
            if entry.peptide == trap.peptide {
                continue;
            }

            // Check same DIA window.
            let entry_window = self.find_dia_window(entry.precursor_mz);
            match (trap_window, entry_window) {
                (Some(tw), Some(ew)) => {
                    if (tw - ew).abs() > f64::EPSILON {
                        continue;
                    }
                }
                _ => continue,
            }

            // Add light candidate.
            candidates.push(CoElutingCandidate {
                peptide: entry.peptide.clone(),
                protein_ids: entry.protein_ids.clone(),
                precursor_mz: entry.precursor_mz,
                charge: entry.charge,
                rt_start: entry.rt_start,
                rt_stop: entry.rt_stop,
                label_form: LabelForm::Light,
                modifications: entry.modifications.clone(),
            });

            if candidates.len() >= self.max_candidates {
                return candidates;
            }

            // If SILAC configured and heavy search enabled, add heavy candidate.
            if let Some(ref silac) = self.silac {
                if silac.enable_heavy_search {
                    let heavy_delta = compute_heavy_delta(&entry.peptide, silac);
                    if heavy_delta.abs() > f64::EPSILON {
                        let charge_f = entry.charge as f64;
                        let heavy_mz = entry.precursor_mz + heavy_delta / charge_f;
                        let residue_deltas = compute_residue_deltas(&entry.peptide, silac);

                        candidates.push(CoElutingCandidate {
                            peptide: entry.peptide.clone(),
                            protein_ids: entry.protein_ids.clone(),
                            precursor_mz: entry.precursor_mz,
                            charge: entry.charge,
                            rt_start: entry.rt_start,
                            rt_stop: entry.rt_stop,
                            label_form: LabelForm::Heavy {
                                precursor_mz_heavy: heavy_mz,
                                residue_deltas,
                            },
                            modifications: entry.modifications.clone(),
                        });

                        if candidates.len() >= self.max_candidates {
                            return candidates;
                        }
                    }
                }
            }
        }

        candidates
    }

    /// Find which DIA window a given precursor m/z falls into.
    ///
    /// Returns the window center if found, or `None` if no window contains
    /// the given m/z.
    fn find_dia_window(&self, mz: f64) -> Option<f64> {
        for w in &self.dia_windows {
            if mz >= w.low && mz <= w.high {
                return Some(w.center);
            }
        }
        None
    }
}

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

/// Compute the total heavy-label mass delta for a peptide sequence.
///
/// Sums `heavy_k_delta` for each Lysine (K) and `heavy_r_delta` for each
/// Arginine (R) in the sequence.
fn compute_heavy_delta(sequence: &str, silac: &SilacConfig) -> f64 {
    let mut delta = 0.0;
    for ch in sequence.chars() {
        match ch {
            'K' => delta += silac.heavy_k_delta,
            'R' => delta += silac.heavy_r_delta,
            _ => {}
        }
    }
    delta
}

/// Compute per-residue heavy-label deltas for a peptide sequence.
///
/// Returns a vec of `(0-based position, delta_Da)` for each K and R residue.
fn compute_residue_deltas(sequence: &str, silac: &SilacConfig) -> Vec<(usize, f64)> {
    let mut deltas = Vec::new();
    for (i, ch) in sequence.chars().enumerate() {
        match ch {
            'K' => deltas.push((i, silac.heavy_k_delta)),
            'R' => deltas.push((i, silac.heavy_r_delta)),
            _ => {}
        }
    }
    deltas
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::PsmGroup;

    fn make_target_psm(
        peptide: &str,
        mz: f64,
        rt_start: f64,
        rt_stop: f64,
        run: &str,
    ) -> (UnifiedPsm, PsmGroup) {
        let psm = UnifiedPsm {
            peptide: peptide.to_string(),
            charge: Some(2),
            precursor_mz: Some(mz),
            retention_time: Some((rt_start + rt_stop) / 2.0),
            rt_start: Some(rt_start),
            rt_stop: Some(rt_stop),
            scan_number: None,
            spectrum_file: Some(run.to_string()),
            protein_ids: "sp|P12345|TEST_HUMAN".to_string(),
            q_value: Some(0.001),
            modifications: vec![],
        };
        (psm, PsmGroup::Target)
    }

    fn make_trap_psm(peptide: &str, mz: f64, rt_start: f64, rt_stop: f64, run: &str) -> UnifiedPsm {
        UnifiedPsm {
            peptide: peptide.to_string(),
            charge: Some(2),
            precursor_mz: Some(mz),
            retention_time: Some((rt_start + rt_stop) / 2.0),
            rt_start: Some(rt_start),
            rt_stop: Some(rt_stop),
            scan_number: None,
            spectrum_file: Some(run.to_string()),
            protein_ids: "sp|P99999|TRAP_YEAST".to_string(),
            q_value: Some(0.005),
            modifications: vec![],
        }
    }

    fn default_windows() -> Vec<DiaWindow> {
        vec![
            DiaWindow {
                center: 500.0,
                low: 487.5,
                high: 512.5,
            },
            DiaWindow {
                center: 525.0,
                low: 512.5,
                high: 537.5,
            },
            DiaWindow {
                center: 550.0,
                low: 537.5,
                high: 562.5,
            },
        ]
    }

    #[test]
    fn test_build_index_and_query_basic() {
        // Two overlapping targets found, one far-away target not found.
        let (psm1, g1) = make_target_psm("AAAAAK", 500.0, 30.0, 35.0, "run1");
        let (psm2, g2) = make_target_psm("BBBBCK", 505.0, 32.0, 37.0, "run1");
        let (psm3, g3) = make_target_psm("CCCCDK", 500.0, 60.0, 65.0, "run1"); // far away

        let psms = vec![psm1, psm2, psm3];
        let groups = vec![g1, g2, g3];
        let windows = default_windows();

        let index = CoElutionIndex::build(&psms, &groups, &windows, None, 20);

        // Trap elutes 31.0–34.0, mz=502.0 => same DIA window as AAAAAK and BBBBCK.
        let trap = make_trap_psm("TRAPPEP", 502.0, 31.0, 34.0, "run1");
        let result = index.find_co_eluting(&trap, "run1");

        assert_eq!(
            result.len(),
            2,
            "expected 2 co-eluting targets, got {}",
            result.len()
        );
        let peptides: Vec<&str> = result.iter().map(|c| c.peptide.as_str()).collect();
        assert!(peptides.contains(&"AAAAAK"), "AAAAAK should co-elute");
        assert!(peptides.contains(&"BBBBCK"), "BBBBCK should co-elute");
        // CCCCDK should NOT be found (RT 60-65 doesn't overlap 31-34).
        assert!(!peptides.contains(&"CCCCDK"), "CCCCDK should NOT co-elute");
    }

    #[test]
    fn test_no_rt_overlap() {
        let (psm1, g1) = make_target_psm("AAAAAK", 500.0, 10.0, 15.0, "run1");
        let psms = vec![psm1];
        let groups = vec![g1];
        let windows = default_windows();

        let index = CoElutionIndex::build(&psms, &groups, &windows, None, 20);

        // Trap at RT 20-25, no overlap with target 10-15.
        let trap = make_trap_psm("TRAPPEP", 502.0, 20.0, 25.0, "run1");
        let result = index.find_co_eluting(&trap, "run1");

        assert!(result.is_empty(), "expected no co-eluting targets");
    }

    #[test]
    fn test_different_dia_window() {
        let (psm1, g1) = make_target_psm("AAAAAK", 550.0, 30.0, 35.0, "run1");
        let psms = vec![psm1];
        let groups = vec![g1];
        let windows = default_windows();

        let index = CoElutionIndex::build(&psms, &groups, &windows, None, 20);

        // Trap at mz=500 => DIA window center 500, but target at mz=550 => window center 550.
        let trap = make_trap_psm("TRAPPEP", 500.0, 30.0, 35.0, "run1");
        let result = index.find_co_eluting(&trap, "run1");

        assert!(
            result.is_empty(),
            "expected no co-eluting targets (different DIA window)"
        );
    }

    #[test]
    fn test_silac_heavy_pairing() {
        // Target peptide ending in K should produce both light + heavy candidates.
        let (psm1, g1) = make_target_psm("STTSGHLVYK", 548.0, 30.0, 35.0, "run1");
        let psms = vec![psm1];
        let groups = vec![g1];
        let windows = default_windows();

        let silac = SilacConfig::default();
        let index = CoElutionIndex::build(&psms, &groups, &windows, Some(&silac), 20);

        let trap = make_trap_psm("TRAPPEP", 540.0, 31.0, 34.0, "run1");
        let result = index.find_co_eluting(&trap, "run1");

        assert_eq!(result.len(), 2, "expected light + heavy candidates");

        let light = result
            .iter()
            .find(|c| matches!(c.label_form, LabelForm::Light));
        assert!(light.is_some(), "should have a Light candidate");

        let heavy = result
            .iter()
            .find(|c| matches!(c.label_form, LabelForm::Heavy { .. }));
        assert!(heavy.is_some(), "should have a Heavy candidate");

        if let Some(h) = heavy {
            if let LabelForm::Heavy {
                precursor_mz_heavy,
                residue_deltas,
            } = &h.label_form
            {
                // STTSGHLVYK has 1 K => heavy_delta = 8.014199, charge=2 => +4.0071
                let expected_heavy_mz = 548.0 + 8.014199 / 2.0;
                assert!(
                    (precursor_mz_heavy - expected_heavy_mz).abs() < 0.001,
                    "heavy mz mismatch: got {precursor_mz_heavy}, expected {expected_heavy_mz}"
                );
                // One K residue at position 9 (0-based).
                assert_eq!(residue_deltas.len(), 1);
                assert_eq!(residue_deltas[0].0, 9);
                assert!((residue_deltas[0].1 - 8.014199).abs() < 1e-6);
            } else {
                panic!("expected Heavy label form");
            }
        }
    }

    #[test]
    fn test_overlap_with_smaller_rt_start_is_found() {
        // Regression for the `partition_point` bug: a genuinely-overlapping
        // target whose `rt_start` is SMALLER than the trap's must still be
        // returned. The previous code called `partition_point` on the unsorted
        // `rt_stop` key, which skipped such entries (false negatives).
        //
        // Entries (all in the same DIA window, mz ~500), sorted by rt_start:
        //   E0 rt 10–100  -> overlaps trap 50–60 (small rt_start, large rt_stop)
        //   E1 rt 20–25   -> ends before trap starts (no overlap)
        //   E2 rt 30–35   -> ends before trap starts (no overlap)
        //   E3 rt 55–70   -> overlaps trap 50–60
        //   E4 rt 80–90   -> starts after trap ends (no overlap)
        let (p0, g0) = make_target_psm("OVERLAPBIGK", 500.0, 10.0, 100.0, "run1");
        let (p1, g1) = make_target_psm("NOOVERLAPAK", 501.0, 20.0, 25.0, "run1");
        let (p2, g2) = make_target_psm("NOOVERLAPBK", 502.0, 30.0, 35.0, "run1");
        let (p3, g3) = make_target_psm("OVERLAPLATEK", 503.0, 55.0, 70.0, "run1");
        let (p4, g4) = make_target_psm("AFTERTRAPXK", 504.0, 80.0, 90.0, "run1");

        let psms = vec![p0, p1, p2, p3, p4];
        let groups = vec![g0, g1, g2, g3, g4];
        let windows = default_windows();
        let index = CoElutionIndex::build(&psms, &groups, &windows, None, 20);

        // Trap elutes 50–60, mz=505 => same DIA window (center 500) as all targets.
        let trap = make_trap_psm("TRAPPEP", 505.0, 50.0, 60.0, "run1");
        let result = index.find_co_eluting(&trap, "run1");
        let peptides: Vec<&str> = result.iter().map(|c| c.peptide.as_str()).collect();

        assert!(
            peptides.contains(&"OVERLAPBIGK"),
            "overlapping target with a smaller rt_start than the trap must be found; got {peptides:?}"
        );
        assert!(
            peptides.contains(&"OVERLAPLATEK"),
            "later-eluting overlapping target must be found; got {peptides:?}"
        );
        assert!(
            !peptides.contains(&"NOOVERLAPAK") && !peptides.contains(&"NOOVERLAPBK"),
            "targets ending before the trap starts must NOT be found; got {peptides:?}"
        );
        assert!(
            !peptides.contains(&"AFTERTRAPXK"),
            "target starting after the trap ends must NOT be found; got {peptides:?}"
        );
        assert_eq!(
            result.len(),
            2,
            "exactly two overlapping targets expected; got {peptides:?}"
        );
    }

    #[test]
    fn test_max_candidates_cap() {
        // Create 25 overlapping targets, cap at 20.
        let mut psms = Vec::new();
        let mut groups = Vec::new();
        for i in 0..25 {
            let peptide = format!("PEP{i:03}X");
            let (psm, g) = make_target_psm(&peptide, 500.0 + (i as f64) * 0.1, 30.0, 35.0, "run1");
            psms.push(psm);
            groups.push(g);
        }
        let windows = default_windows();
        let index = CoElutionIndex::build(&psms, &groups, &windows, None, 20);

        let trap = make_trap_psm("TRAPPEP", 502.0, 31.0, 34.0, "run1");
        let result = index.find_co_eluting(&trap, "run1");

        assert!(
            result.len() <= 20,
            "expected at most 20 candidates, got {}",
            result.len()
        );
    }
}
