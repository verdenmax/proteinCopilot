//! Scan matching: associates imported PSMs with mzML MS2 scans.
//!
//! Algorithm:
//! 1. Pre-scan mzML to collect (scan_number, rt_min, isolation_window) for all MS2 spectra
//! 2. Sort by RT
//! 3. For each PSM: binary search for RT-proximate candidates, filter by isolation window,
//!    pick the closest RT match

use std::collections::HashMap;
use std::ffi::OsStr;
use std::path::{Path, PathBuf};

use protein_copilot_spectrum_io::reader::SpectrumReader;

use crate::{FileMatchStats, ImportedPsm, MatchReport, ResultImportError};

/// MS2 spectrum info extracted from mzML for scan matching.
#[derive(Debug, Clone)]
pub struct Ms2Info {
    pub scan_number: u32,
    pub rt_min: f64,
    /// (target_mz, lower_offset, upper_offset)
    pub isolation_window: Option<(f64, f64, f64)>,
}

/// Scan matcher configuration.
pub struct ScanMatcherConfig {
    pub rt_tolerance_min: f64,
    pub mzml_dir: PathBuf,
}

/// Match imported PSMs to mzML MS2 scans.
///
/// PSMs are matched by:
/// 1. `raw_name` → mzML file (raw_name + ".mzML" in `mzml_dir`)
/// 2. RT proximity (within `rt_tolerance_min`)
/// 3. precursor_mz falls within the MS2's isolation window
///
/// Returns the mutated PSMs (with `matched_scan` and `rt_delta_min` filled)
/// and a `MatchReport` with quality statistics.
pub type ReaderFactory = dyn Fn(&Path) -> Result<Box<dyn SpectrumReader>, ResultImportError>;

pub fn match_scans(
    psms: &mut [ImportedPsm],
    config: &ScanMatcherConfig,
    reader_factory: &ReaderFactory,
) -> Result<MatchReport, ResultImportError> {
    let _span = tracing::info_span!("match_scans",
        psm_count = psms.len(),
    ).entered();

    // Group PSMs by raw_name
    let mut groups: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, psm) in psms.iter().enumerate() {
        groups.entry(psm.raw_name.clone()).or_default().push(i);
    }

    let mut per_file = HashMap::new();
    let mut all_rt_deltas: Vec<f64> = Vec::new();
    let mut total_matched = 0usize;
    let mut total_unmatched = 0usize;
    let total_psm_count = psms.len();
    let progress_interval: usize = 1000;
    let loop_start = std::time::Instant::now();
    let mut processed_count: usize = 0;

    for (raw_name, indices) in &groups {
        // `raw_name` originates from an external search-result file. Confine it
        // to `mzml_dir`: it must be a bare file-name component so it cannot
        // escape the directory via path separators or `..` (path traversal).
        validate_raw_name(raw_name)?;

        let mzml_path = config.mzml_dir.join(format!("{raw_name}.mzML"));
        let mzml_path_lower = config.mzml_dir.join(format!("{raw_name}.mzml"));

        let actual_path = if mzml_path.exists() {
            mzml_path
        } else if mzml_path_lower.exists() {
            mzml_path_lower
        } else {
            let available = list_mzml_files(&config.mzml_dir);
            return Err(ResultImportError::MzmlNotFound {
                raw_name: raw_name.clone(),
                dir: config.mzml_dir.clone(),
                available,
            });
        };

        let reader = reader_factory(&actual_path)?;

        // Pre-scan all MS2 spectra
        let ms2_infos = collect_ms2_info(&*reader, &actual_path)?;
        let ms2_count = ms2_infos.len();

        // Sort by RT for binary search
        let mut sorted_ms2 = ms2_infos;
        sorted_ms2.sort_by(|a, b| {
            a.rt_min
                .partial_cmp(&b.rt_min)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let mut file_matched = 0usize;
        let mut file_unmatched = 0usize;

        for &idx in indices {
            let psm = &psms[idx];
            if let Some((scan, delta)) = find_best_match(
                &sorted_ms2,
                psm.rt_min,
                psm.precursor_mz,
                config.rt_tolerance_min,
            ) {
                let psm_mut = &mut psms[idx];
                psm_mut.matched_scan = Some(scan);
                psm_mut.rt_delta_min = Some(delta);
                all_rt_deltas.push(delta.abs());
                file_matched += 1;
            } else {
                file_unmatched += 1;
            }
            processed_count += 1;
            if processed_count % progress_interval == 0 || processed_count == total_psm_count {
                let elapsed = loop_start.elapsed().as_secs_f64();
                let rate = if elapsed > 0.0 { processed_count as f64 / elapsed } else { 0.0 };
                let eta = if rate > 0.0 { (total_psm_count - processed_count) as f64 / rate } else { 0.0 };
                tracing::info!(
                    progress = processed_count,
                    total = total_psm_count,
                    rate = format!("{:.0}/s", rate),
                    eta_sec = format!("{:.1}", eta),
                    "matching scans"
                );
            }
        }

        total_matched += file_matched;
        total_unmatched += file_unmatched;

        per_file.insert(
            raw_name.clone(),
            FileMatchStats {
                total: indices.len(),
                matched: file_matched,
                ms2_count,
            },
        );
    }

    all_rt_deltas.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median_rt_delta = if all_rt_deltas.is_empty() {
        0.0
    } else {
        all_rt_deltas[all_rt_deltas.len() / 2]
    };
    let max_rt_delta = all_rt_deltas.last().copied().unwrap_or(0.0);

    tracing::info!(matched = total_matched, unmatched = total_unmatched, "scan matching complete");

    Ok(MatchReport {
        total_psms: psms.len(),
        matched: total_matched,
        unmatched: total_unmatched,
        median_rt_delta_min: median_rt_delta,
        max_rt_delta_min: max_rt_delta,
        per_file,
    })
}

/// Find the best MS2 scan matching a given RT and precursor m/z.
///
/// Single-lookup convenience wrapper around `collect_ms2_info` + `find_best_match`.
/// Used by `annotate_spectrum` and `extract_xic` when the user provides RT instead
/// of scan_number.
///
/// # Arguments
/// * `file` — path to the mzML file
/// * `rt_min` — target retention time in minutes
/// * `precursor_mz` — precursor m/z to match against isolation windows
/// * `rt_tolerance_min` — RT tolerance in minutes (default: 0.5)
/// * `reader` — spectrum reader for the file
pub fn find_scan_by_rt(
    file: &Path,
    rt_min: f64,
    precursor_mz: f64,
    rt_tolerance_min: f64,
    reader: &dyn SpectrumReader,
) -> Result<u32, ResultImportError> {
    // Delegate to SpectrumReader::find_by_rt() which uses O(log N) binary
    // search on IndexedMzMLReader (falls back to read_all on other readers).
    reader
        .find_by_rt(file, rt_min, precursor_mz, rt_tolerance_min)
        .map_err(|e| ResultImportError::SpectrumIo(e.to_string()))?
        .map(|(scan, _delta)| scan)
        .ok_or(ResultImportError::NoMatchingScan {
            rt_min,
            tolerance_min: rt_tolerance_min,
            precursor_mz,
        })
}

/// Collect (scan, rt, isolation_window) for all MS2 spectra from a reader.
pub fn collect_ms2_info(
    reader: &dyn SpectrumReader,
    path: &Path,
) -> Result<Vec<Ms2Info>, ResultImportError> {
    let meta = reader
        .list_ms2_meta(path)
        .map_err(|e| ResultImportError::SpectrumIo(e.to_string()))?;
    Ok(meta
        .into_iter()
        .map(|m| Ms2Info {
            scan_number: m.scan_number,
            rt_min: m.rt_min,
            isolation_window: m.isolation_window,
        })
        .collect())
}

/// Find the best matching MS2 for a given RT and precursor m/z.
///
/// Returns `(scan_number, rt_delta_min)` or `None` if no match found.
pub fn find_best_match(
    sorted_ms2: &[Ms2Info],
    psm_rt_min: f64,
    psm_precursor_mz: f64,
    rt_tolerance_min: f64,
) -> Option<(u32, f64)> {
    if sorted_ms2.is_empty() {
        return None;
    }

    // Binary search for the start of the RT window
    let insert_pos = sorted_ms2.partition_point(|m| m.rt_min < psm_rt_min - rt_tolerance_min);

    let mut best: Option<(u32, f64)> = None;

    for ms2 in sorted_ms2[insert_pos..].iter() {
        let delta = ms2.rt_min - psm_rt_min;
        if delta > rt_tolerance_min {
            break; // past the RT window
        }
        if delta.abs() > rt_tolerance_min {
            continue;
        }

        // Check isolation window
        if let Some((target, lower, upper)) = ms2.isolation_window {
            let low = target - lower;
            let high = target + upper;
            if psm_precursor_mz < low || psm_precursor_mz > high {
                continue; // precursor_mz outside isolation window
            }
        }
        // else: no isolation window info → accept based on RT only (DDA fallback)

        match &best {
            None => best = Some((ms2.scan_number, delta)),
            Some((_, best_delta)) => {
                if delta.abs() < best_delta.abs() {
                    best = Some((ms2.scan_number, delta));
                }
            }
        }
    }

    best
}

/// Validates that an externally supplied `raw_name` is a bare file-name
/// component, so that `mzml_dir.join(format!("{raw_name}.mzML"))` cannot escape
/// `mzml_dir` via path separators or `..` (path-traversal confinement).
fn validate_raw_name(raw_name: &str) -> Result<(), ResultImportError> {
    let invalid = raw_name.is_empty()
        || raw_name == "."
        || raw_name == ".."
        || raw_name.contains('/')
        || raw_name.contains('\\')
        || raw_name.contains(std::path::MAIN_SEPARATOR)
        || Path::new(raw_name).file_name() != Some(OsStr::new(raw_name));
    if invalid {
        return Err(ResultImportError::InvalidRawName {
            raw_name: raw_name.to_string(),
        });
    }
    Ok(())
}

/// List .mzML files in a directory for error messages.
fn list_mzml_files(dir: &Path) -> String {
    match std::fs::read_dir(dir) {
        Ok(entries) => {
            let files: Vec<String> = entries
                .filter_map(|e| e.ok())
                .filter(|e| {
                    e.path()
                        .extension()
                        .is_some_and(|ext| ext.eq_ignore_ascii_case("mzml"))
                })
                .map(|e| e.file_name().to_string_lossy().to_string())
                .collect();
            if files.is_empty() {
                "none".to_string()
            } else {
                files.join(", ")
            }
        }
        Err(_) => "directory not readable".to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn make_psm(raw_name: &str) -> ImportedPsm {
        ImportedPsm {
            sequence: "PEPTIDEK".to_string(),
            charge: 2,
            precursor_mz: 500.0,
            rt_min: 100.0,
            modifications: Vec::new(),
            score: None,
            q_value: None,
            protein_accessions: Vec::new(),
            raw_name: raw_name.to_string(),
            matched_scan: None,
            rt_delta_min: None,
        }
    }

    #[test]
    fn match_scans_rejects_path_traversal_raw_name() {
        let dir = TempDir::new().unwrap();
        let config = ScanMatcherConfig {
            rt_tolerance_min: 0.5,
            mzml_dir: dir.path().to_path_buf(),
        };
        let called = std::rc::Rc::new(std::cell::RefCell::new(Vec::<PathBuf>::new()));
        let called_in = std::rc::Rc::clone(&called);
        let factory = move |p: &Path| -> Result<Box<dyn SpectrumReader>, ResultImportError> {
            called_in.borrow_mut().push(p.to_path_buf());
            Err(ResultImportError::Other("factory-reached".to_string()))
        };

        let mut psms = vec![make_psm("../evil")];
        let err = match_scans(&mut psms, &config, &factory).unwrap_err();

        assert!(
            matches!(err, ResultImportError::InvalidRawName { .. }),
            "expected InvalidRawName, got {err:?}"
        );
        assert!(
            called.borrow().is_empty(),
            "reader factory must not be invoked for a traversal raw_name"
        );
        // The PSM must not have been matched/escaped.
        assert!(psms[0].matched_scan.is_none());
    }

    #[test]
    fn match_scans_accepts_bare_raw_name() {
        let dir = TempDir::new().unwrap();
        // The legitimate file must exist so matching reaches the read path.
        std::fs::write(dir.path().join("sample1.mzML"), b"").unwrap();
        let config = ScanMatcherConfig {
            rt_tolerance_min: 0.5,
            mzml_dir: dir.path().to_path_buf(),
        };
        let called = std::rc::Rc::new(std::cell::RefCell::new(Vec::<PathBuf>::new()));
        let called_in = std::rc::Rc::clone(&called);
        let factory = move |p: &Path| -> Result<Box<dyn SpectrumReader>, ResultImportError> {
            called_in.borrow_mut().push(p.to_path_buf());
            Err(ResultImportError::Other("factory-reached".to_string()))
        };

        let mut psms = vec![make_psm("sample1")];
        let err = match_scans(&mut psms, &config, &factory).unwrap_err();

        // Confinement passed; we reached the (factory-driven) read path, not
        // an InvalidRawName rejection.
        assert!(
            matches!(err, ResultImportError::Other(_)),
            "expected to reach the factory read path, got {err:?}"
        );
        assert_eq!(called.borrow().len(), 1);
        assert!(called.borrow()[0].ends_with("sample1.mzML"));
    }

    fn make_ms2_infos() -> Vec<Ms2Info> {
        vec![
            Ms2Info {
                scan_number: 10,
                rt_min: 100.0,
                isolation_window: Some((500.0, 12.5, 12.5)),
            },
            Ms2Info {
                scan_number: 20,
                rt_min: 200.0,
                isolation_window: Some((600.0, 12.5, 12.5)),
            },
            Ms2Info {
                scan_number: 30,
                rt_min: 300.0,
                isolation_window: Some((500.0, 12.5, 12.5)),
            },
            Ms2Info {
                scan_number: 40,
                rt_min: 400.0,
                isolation_window: Some((700.0, 12.5, 12.5)),
            },
            Ms2Info {
                scan_number: 50,
                rt_min: 500.0,
                isolation_window: Some((500.0, 25.0, 25.0)),
            },
        ]
    }

    #[test]
    fn find_best_match_exact_rt_and_mz() {
        let ms2s = make_ms2_infos();
        let result = find_best_match(&ms2s, 100.0, 500.0, 30.0);
        assert_eq!(result, Some((10, 0.0)));
    }

    #[test]
    fn find_best_match_within_tolerance() {
        let ms2s = make_ms2_infos();
        // RT=198, mz=605 → should match scan 20 (RT=200, window 587.5–612.5)
        let result = find_best_match(&ms2s, 198.0, 605.0, 30.0);
        assert_eq!(result.unwrap().0, 20);
        assert!((result.unwrap().1 - 2.0).abs() < 0.01);
    }

    #[test]
    fn find_best_match_mz_outside_window() {
        let ms2s = make_ms2_infos();
        // RT=100, mz=550 → scan 10 has window 487.5–512.5, mz=550 is outside
        let result = find_best_match(&ms2s, 100.0, 550.0, 30.0);
        assert!(result.is_none());
    }

    #[test]
    fn find_best_match_rt_outside_tolerance() {
        let ms2s = make_ms2_infos();
        // RT=150, tolerance=30 → nearest scan 10 (RT=100) is 50 min away
        let result = find_best_match(&ms2s, 150.0, 500.0, 30.0);
        assert!(result.is_none());
    }

    #[test]
    fn find_best_match_wide_dia_window() {
        let ms2s = make_ms2_infos();
        // scan 50: RT=500, window 475–525 (wide DIA)
        let result = find_best_match(&ms2s, 502.0, 520.0, 30.0);
        assert_eq!(result.unwrap().0, 50);
    }

    #[test]
    fn find_best_match_picks_closest_rt() {
        let ms2s = vec![
            Ms2Info {
                scan_number: 1,
                rt_min: 100.0,
                isolation_window: Some((500.0, 25.0, 25.0)),
            },
            Ms2Info {
                scan_number: 2,
                rt_min: 105.0,
                isolation_window: Some((500.0, 25.0, 25.0)),
            },
        ];
        // PSM at RT=103 → closer to scan 2 (105)
        let result = find_best_match(&ms2s, 103.0, 500.0, 30.0);
        assert_eq!(result.unwrap().0, 2);
    }

    #[test]
    fn find_best_match_no_isolation_window_fallback() {
        let ms2s = vec![Ms2Info {
            scan_number: 1,
            rt_min: 100.0,
            isolation_window: None,
        }];
        let result = find_best_match(&ms2s, 105.0, 999.0, 30.0);
        assert_eq!(result.unwrap().0, 1);
    }

    #[test]
    fn find_best_match_empty_ms2_list() {
        let result = find_best_match(&[], 100.0, 500.0, 30.0);
        assert!(result.is_none());
    }
}
