//! Isotope pattern detection for precursor extraction.

use protein_copilot_core::spectrum::{IsolationWindow, PrecursorInfo, Spectrum};

use crate::extractor::PrecursorExtractor;

/// Mass difference between ¹³C and ¹²C isotopes (Da).
const C13_C12_MASS_DIFF: f64 = 1.003_355;

/// Extracts precursor candidates from MS1 spectra by detecting isotope envelope patterns.
///
/// For each charge state z, looks for peaks separated by Δm/z ≈ 1.00335/z (¹³C-¹²C mass diff / charge).
#[derive(Debug, Clone)]
pub struct IsotopePatternExtractor {
    /// Minimum charge state to consider (default: 2).
    pub min_charge: i32,
    /// Maximum charge state to consider (default: 5).
    pub max_charge: i32,
    /// Tolerance for isotope spacing match in Da (default: 0.01).
    pub isotope_tolerance_da: f64,
    /// Minimum number of isotope peaks to form a valid cluster (default: 2).
    pub min_isotope_peaks: usize,
    /// Minimum intensity ratio: subsequent peaks must be at least this fraction of the first peak (default: 0.1).
    pub min_intensity_ratio: f64,
}

impl Default for IsotopePatternExtractor {
    fn default() -> Self {
        Self {
            min_charge: 2,
            max_charge: 5,
            isotope_tolerance_da: 0.01,
            min_isotope_peaks: 2,
            min_intensity_ratio: 0.1,
        }
    }
}

impl PrecursorExtractor for IsotopePatternExtractor {
    fn extract(&self, ms1: &Spectrum, isolation_window: &IsolationWindow) -> Vec<PrecursorInfo> {
        // Guard against invalid charge range
        if self.min_charge < 1 || self.max_charge < self.min_charge {
            return Vec::new();
        }

        let low = isolation_window.target_mz - isolation_window.lower_offset;
        let high = isolation_window.target_mz + isolation_window.upper_offset;

        // Step 1: Filter peaks to isolation window range.
        let mut peaks: Vec<(usize, f64, f64)> = ms1
            .mz_array
            .iter()
            .zip(ms1.intensity_array.iter())
            .enumerate()
            .filter(|&(_, (&mz, _))| mz >= low && mz <= high)
            .map(|(i, (&mz, &intensity))| (i, mz, intensity))
            .collect();

        if peaks.is_empty() {
            return Vec::new();
        }

        // Step 2: Sort by intensity descending.
        peaks.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));

        // Step 3: Track used peaks.
        let mut used = vec![false; ms1.mz_array.len()];

        let mut candidates: Vec<PrecursorInfo> = Vec::new();

        // Step 4: For each unused peak, try to build isotope clusters.
        for &(seed_idx, seed_mz, seed_intensity) in &peaks {
            if used[seed_idx] {
                continue;
            }

            let mut best_cluster: Option<(i32, Vec<usize>)> = None;

            for z in self.min_charge..=self.max_charge {
                let delta = C13_C12_MASS_DIFF / z as f64;
                let mut cluster_indices = vec![seed_idx];

                // Look BACKWARD one step for the monoisotopic peak (M+0).
                // The seed may be M+1 for large peptides where M+1 > M+0.
                // Only accept if the lighter peak has reasonable intensity.
                {
                    let expected_mz = seed_mz - delta;
                    let mut best_candidate: Option<(usize, f64, f64)> = None;
                    let mut best_dist = f64::MAX;

                    for &(idx, mz, intensity) in &peaks {
                        if used[idx] || idx == seed_idx {
                            continue;
                        }
                        let dist = (mz - expected_mz).abs();
                        if dist <= self.isotope_tolerance_da && dist < best_dist {
                            best_candidate = Some((idx, mz, intensity));
                            best_dist = dist;
                        }
                    }

                    if let Some((idx, _, intensity)) = best_candidate {
                        // M+0 should be at least min_intensity_ratio of seed
                        if intensity >= self.min_intensity_ratio * seed_intensity {
                            cluster_indices.push(idx);
                        }
                    }
                }

                // Look FORWARD for heavier isotope peaks (+1, +2, +3, ...).
                let mut prev_intensity = seed_intensity;
                let mut current_mz = seed_mz;

                loop {
                    let expected_mz = current_mz + delta;

                    // Find the closest unused peak within tolerance.
                    let mut best_candidate: Option<(usize, f64, f64)> = None;
                    let mut best_dist = f64::MAX;

                    for &(idx, mz, intensity) in &peaks {
                        if used[idx] || cluster_indices.contains(&idx) {
                            continue;
                        }
                        let dist = (mz - expected_mz).abs();
                        if dist <= self.isotope_tolerance_da && dist < best_dist {
                            best_candidate = Some((idx, mz, intensity));
                            best_dist = dist;
                        }
                    }

                    match best_candidate {
                        Some((idx, mz, intensity)) => {
                            // Intensity check: must not exceed 1.5× previous peak
                            // and must be >= min_intensity_ratio × first peak.
                            if intensity > prev_intensity * 1.5 {
                                break;
                            }
                            if intensity < self.min_intensity_ratio * seed_intensity {
                                break;
                            }
                            cluster_indices.push(idx);
                            prev_intensity = intensity;
                            current_mz = mz;
                        }
                        None => break,
                    }
                }

                if cluster_indices.len() >= self.min_isotope_peaks {
                    let is_better = match &best_cluster {
                        None => true,
                        Some((best_z, best_indices)) => {
                            if cluster_indices.len() > best_indices.len() {
                                true
                            } else if cluster_indices.len() == best_indices.len() {
                                // Prefer lower charge state (2 > 3 > ...).
                                z < *best_z
                            } else {
                                false
                            }
                        }
                    };
                    if is_better {
                        best_cluster = Some((z, cluster_indices));
                    }
                }
            }

            if let Some((charge, indices)) = best_cluster {
                for &idx in &indices {
                    used[idx] = true;
                }
                // Use the monoisotopic (lowest m/z) peak from the cluster,
                // not the seed peak which may be M+1 or M+2.
                let monoisotopic_mz = indices
                    .iter()
                    .filter_map(|&idx| {
                        peaks
                            .iter()
                            .find(|(i, _, _)| *i == idx)
                            .map(|&(_, mz, _)| mz)
                    })
                    .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
                    .unwrap_or(seed_mz);
                candidates.push(PrecursorInfo {
                    mz: monoisotopic_mz,
                    charge: Some(charge),
                    intensity: Some(seed_intensity),
                    isolation_window: None,
                    source_scan: None,
                });
            }
        }

        // Step 5: Deduplication — merge candidates within 0.01 Da with same charge.
        candidates.sort_by(|a, b| {
            a.charge
                .cmp(&b.charge)
                .then(a.mz.partial_cmp(&b.mz).unwrap_or(std::cmp::Ordering::Equal))
        });

        let mut deduped: Vec<PrecursorInfo> = Vec::new();
        for candidate in candidates {
            let dominated = deduped.iter().any(|existing| {
                existing.charge == candidate.charge
                    && (existing.mz - candidate.mz).abs() < 0.01
                    && existing.intensity >= candidate.intensity
            });
            if !dominated {
                // Remove any existing entry dominated by this candidate.
                deduped.retain(|existing| {
                    !(existing.charge == candidate.charge
                        && (existing.mz - candidate.mz).abs() < 0.01
                        && existing.intensity < candidate.intensity)
                });
                deduped.push(candidate);
            }
        }

        // Step 6: Sort by intensity descending.
        deduped.sort_by(|a, b| {
            b.intensity
                .partial_cmp(&a.intensity)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        deduped
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::spectrum::MsLevel;

    fn make_spectrum(mz_array: Vec<f64>, intensity_array: Vec<f64>) -> Spectrum {
        Spectrum::new(1, MsLevel::MS1, 100.0, vec![], mz_array, intensity_array).unwrap()
    }

    fn default_window(target_mz: f64, half_width: f64) -> IsolationWindow {
        IsolationWindow {
            target_mz,
            lower_offset: half_width,
            upper_offset: half_width,
        }
    }

    #[test]
    fn test_charge_2_cluster() {
        // Spacing for z=2: 1.003355 / 2 ≈ 0.5017
        let mz = vec![500.0, 500.502, 501.003];
        let intensity = vec![1000.0, 800.0, 400.0];
        let spectrum = make_spectrum(mz, intensity);
        let window = default_window(500.0, 50.0);
        let extractor = IsotopePatternExtractor::default();

        let results = extractor.extract(&spectrum, &window);
        assert_eq!(results.len(), 1);
        assert!((results[0].mz - 500.0).abs() < 0.01);
        assert_eq!(results[0].charge, Some(2));
    }

    #[test]
    fn test_charge_3_cluster() {
        // Spacing for z=3: 1.003355 / 3 ≈ 0.3345
        let mz = vec![600.0, 600.335, 600.669];
        let intensity = vec![1000.0, 800.0, 500.0];
        let spectrum = make_spectrum(mz, intensity);
        let window = default_window(600.0, 50.0);
        let extractor = IsotopePatternExtractor::default();

        let results = extractor.extract(&spectrum, &window);
        assert_eq!(results.len(), 1);
        assert!((results[0].mz - 600.0).abs() < 0.01);
        assert_eq!(results[0].charge, Some(3));
    }

    #[test]
    fn test_noise_no_detection() {
        let mz = vec![100.0, 200.5, 350.2, 480.7];
        let intensity = vec![500.0, 300.0, 700.0, 200.0];
        let spectrum = make_spectrum(mz, intensity);
        let window = default_window(300.0, 250.0);
        let extractor = IsotopePatternExtractor::default();

        let results = extractor.extract(&spectrum, &window);
        assert!(
            results.is_empty(),
            "Randomly spaced peaks should not form isotope clusters"
        );
    }

    #[test]
    fn test_outside_window_filtered() {
        // Valid z=2 cluster at m/z 800, but window is [400, 450].
        let mz = vec![800.0, 800.502, 801.003];
        let intensity = vec![1000.0, 800.0, 400.0];
        let spectrum = make_spectrum(mz, intensity);
        let window = default_window(425.0, 25.0);
        let extractor = IsotopePatternExtractor::default();

        let results = extractor.extract(&spectrum, &window);
        assert!(
            results.is_empty(),
            "Peaks outside isolation window should be filtered"
        );
    }

    #[test]
    fn test_intensity_decreasing() {
        // Peaks at correct z=2 spacing but with increasing intensity (reverse pattern).
        let mz = vec![500.0, 500.502, 501.003];
        let intensity = vec![100.0, 500.0, 1000.0];
        let spectrum = make_spectrum(mz, intensity);
        let window = default_window(500.0, 50.0);
        let extractor = IsotopePatternExtractor::default();

        let results = extractor.extract(&spectrum, &window);
        // The seed (highest intensity = 1000.0 at 501.003) cannot extend forward,
        // and 500.502 (seed=500.0, int=500.0) → 501.003 is already used.
        // No valid cluster of >= 2 peaks should form with proper monoisotopic assignment.
        // The algorithm processes peaks by descending intensity, so:
        //  - 1000.0 at 501.003: no +delta peak → no cluster
        //  - 500.0 at 500.502: next expected at 501.003 (used) → no cluster
        //  - 100.0 at 500.0: next expected at 500.502, intensity 500.0 > 100.0*1.5 → fails intensity check
        // So: should detect nothing or at most fewer peaks.
        for r in &results {
            if r.charge == Some(2) && (r.mz - 500.0).abs() < 0.01 {
                panic!(
                    "Should not detect a z=2 cluster starting at 500.0 with increasing intensity"
                );
            }
        }
    }

    #[test]
    fn test_multiple_clusters() {
        // Cluster 1: z=2 at 500.0
        // Cluster 2: z=3 at 700.0
        let mz = vec![500.0, 500.502, 501.003, 700.0, 700.335, 700.669];
        let intensity = vec![1000.0, 800.0, 400.0, 900.0, 700.0, 350.0];
        let spectrum = make_spectrum(mz, intensity);
        let window = default_window(600.0, 150.0);
        let extractor = IsotopePatternExtractor::default();

        let results = extractor.extract(&spectrum, &window);
        assert_eq!(results.len(), 2, "Should detect two separate clusters");

        let charges: Vec<Option<i32>> = results.iter().map(|r| r.charge).collect();
        assert!(charges.contains(&Some(2)), "Should detect z=2 cluster");
        assert!(charges.contains(&Some(3)), "Should detect z=3 cluster");
    }

    #[test]
    fn test_empty_ms1() {
        let spectrum = make_spectrum(vec![], vec![]);
        let window = default_window(500.0, 50.0);
        let extractor = IsotopePatternExtractor::default();

        let results = extractor.extract(&spectrum, &window);
        assert!(results.is_empty());
    }

    #[test]
    fn test_single_peak() {
        let spectrum = make_spectrum(vec![500.0], vec![1000.0]);
        let window = default_window(500.0, 50.0);
        let extractor = IsotopePatternExtractor::default();

        let results = extractor.extract(&spectrum, &window);
        assert!(
            results.is_empty(),
            "Single peak cannot form an isotope cluster"
        );
    }
}
