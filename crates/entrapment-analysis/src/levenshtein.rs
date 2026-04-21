//! Levenshtein edit distance with alignment backtracking and delta-mass calculation.
//!
//! Used to find similar target peptides that may differ by insertions, deletions,
//! or substitutions — not just same-length Hamming distance.

use protein_copilot_search_engine::digest::residue_mass;

/// Result of a Levenshtein alignment between two sequences.
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    /// Edit distance (number of substitutions + insertions + deletions).
    pub edit_distance: u32,
    /// Mass difference: sum of mass changes from edit operations (Da).
    /// Positive means target is heavier than query.
    pub delta_mass_da: f64,
    /// Human-readable alignment detail, e.g. "D0→N,ins:G@5".
    pub alignment_detail: String,
}

/// Compute the Levenshtein edit distance between two sequences.
///
/// Returns only the distance (no alignment backtracking).
/// Uses O(min(m,n)) space with single-row optimization.
pub fn edit_distance(a: &str, b: &str) -> u32 {
    let a_bytes = a.as_bytes();
    let b_bytes = b.as_bytes();
    let (short, long) = if a_bytes.len() <= b_bytes.len() {
        (a_bytes, b_bytes)
    } else {
        (b_bytes, a_bytes)
    };

    let n = short.len();
    let m = long.len();

    let mut prev: Vec<u32> = (0..=(n as u32)).collect();
    let mut curr = vec![0u32; n + 1];

    for i in 1..=m {
        curr[0] = i as u32;
        for j in 1..=n {
            let cost = if long[i - 1] == short[j - 1] { 0 } else { 1 };
            curr[j] = (prev[j] + 1) // deletion
                .min(curr[j - 1] + 1) // insertion
                .min(prev[j - 1] + cost); // substitution
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    prev[n]
}

/// Compute Levenshtein alignment with full backtracking.
///
/// Returns edit distance, mass difference, and a human-readable detail string.
/// The mass difference is computed from the alignment: for each edit operation,
/// the residue mass change is accumulated.
///
/// - Substitution A→B: Δm += mass(B) - mass(A)
/// - Insertion of B (in target, not in query): Δm += mass(B)
/// - Deletion of A (in query, not in target): Δm -= mass(A)
///
/// `a` is the query (trap peptide), `b` is the target peptide.
pub fn align(a: &str, b: &str) -> AlignmentResult {
    let a_chars: Vec<char> = a.chars().collect();
    let b_chars: Vec<char> = b.chars().collect();
    let m = a_chars.len();
    let n = b_chars.len();

    // Build full DP matrix for backtracking
    let mut dp = vec![vec![0u32; n + 1]; m + 1];
    for i in 0..=m {
        dp[i][0] = i as u32;
    }
    for j in 0..=n {
        dp[0][j] = j as u32;
    }
    for i in 1..=m {
        for j in 1..=n {
            let cost = if a_chars[i - 1] == b_chars[j - 1] {
                0
            } else {
                1
            };
            dp[i][j] = (dp[i - 1][j] + 1)
                .min(dp[i][j - 1] + 1)
                .min(dp[i - 1][j - 1] + cost);
        }
    }

    // Backtrack to reconstruct alignment
    let mut ops: Vec<String> = Vec::new();
    let mut delta_mass: f64 = 0.0;
    let mut i = m;
    let mut j = n;

    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let cost = if a_chars[i - 1] == b_chars[j - 1] {
                0
            } else {
                1
            };
            if dp[i][j] == dp[i - 1][j - 1] + cost {
                if cost == 1 {
                    // Substitution: a[i-1] → b[j-1]
                    let ca = a_chars[i - 1];
                    let cb = b_chars[j - 1];
                    let dm = match (residue_mass(ca), residue_mass(cb)) {
                        (Some(ma), Some(mb)) => mb - ma,
                        _ => 0.0,
                    };
                    delta_mass += dm;
                    ops.push(format!("{ca}{}→{cb}", i - 1));
                }
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if j > 0 && dp[i][j] == dp[i][j - 1] + 1 {
            // Insertion: b[j-1] is extra in target (not in query)
            let cb = b_chars[j - 1];
            if let Some(mb) = residue_mass(cb) {
                delta_mass += mb;
            }
            ops.push(format!("ins:{cb}@{i}"));
            j -= 1;
        } else if i > 0 && dp[i][j] == dp[i - 1][j] + 1 {
            // Deletion: a[i-1] is in query but not in target
            let ca = a_chars[i - 1];
            if let Some(ma) = residue_mass(ca) {
                delta_mass -= ma;
            }
            ops.push(format!("del:{ca}@{}", i - 1));
            i -= 1;
        } else {
            // Should not reach here, but break to avoid infinite loop
            break;
        }
    }

    ops.reverse();
    let detail = ops.join(",");

    AlignmentResult {
        edit_distance: dp[m][n],
        delta_mass_da: delta_mass,
        alignment_detail: detail,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- edit_distance tests ---

    #[test]
    fn test_edit_distance_identical() {
        assert_eq!(edit_distance("PEPTIDEK", "PEPTIDEK"), 0);
    }

    #[test]
    fn test_edit_distance_one_substitution() {
        assert_eq!(edit_distance("PEPTIDEK", "PEPTLDEK"), 1);
    }

    #[test]
    fn test_edit_distance_one_insertion() {
        // "PEPTIDEK" vs "PEPGTIDEK" — G inserted at pos 3
        assert_eq!(edit_distance("PEPTIDEK", "PEPGTIDEK"), 1);
    }

    #[test]
    fn test_edit_distance_one_deletion() {
        // "PEPTIDEK" vs "PEPTDEK" — I deleted
        assert_eq!(edit_distance("PEPTIDEK", "PEPTDEK"), 1);
    }

    #[test]
    fn test_edit_distance_two_edits() {
        // Two substitutions
        assert_eq!(edit_distance("AGCDEK", "ANCDER"), 2);
    }

    #[test]
    fn test_edit_distance_empty() {
        assert_eq!(edit_distance("", ""), 0);
        assert_eq!(edit_distance("ABC", ""), 3);
        assert_eq!(edit_distance("", "ABC"), 3);
    }

    #[test]
    fn test_edit_distance_equals_hamming_for_same_length() {
        // For equal-length strings with single substitution, edit distance equals hamming
        assert_eq!(edit_distance("DGFLLDGFPR", "NGFLLDGFPR"), 1);
    }

    // --- align tests ---

    #[test]
    fn test_align_identical() {
        let r = align("PEPTIDEK", "PEPTIDEK");
        assert_eq!(r.edit_distance, 0);
        assert!(r.delta_mass_da.abs() < 1e-6);
        assert_eq!(r.alignment_detail, "");
    }

    #[test]
    fn test_align_one_substitution_dn() {
        // D→N at position 0
        let r = align("DGFLLDGFPR", "NGFLLDGFPR");
        assert_eq!(r.edit_distance, 1);
        // N(114.042927) - D(115.026943) = -0.984016
        assert!((r.delta_mass_da - (-0.984016)).abs() < 0.001);
        assert_eq!(r.alignment_detail, "D0→N");
    }

    #[test]
    fn test_align_insertion() {
        // query: "PEPNDE", target: "PEPGGDE" — N→GG, edit distance 2
        let r = align("PEPNDE", "PEPGGDE");
        assert_eq!(r.edit_distance, 2);
        assert!(!r.alignment_detail.is_empty());
    }

    #[test]
    fn test_align_qk_substitution() {
        // Q→K at position 4
        let r = align("PEPTQDEK", "PEPTKDEK");
        assert_eq!(r.edit_distance, 1);
        // K(128.094963) - Q(128.058578) = 0.036385
        assert!((r.delta_mass_da - 0.036385).abs() < 0.001);
        assert_eq!(r.alignment_detail, "Q4→K");
    }

    #[test]
    fn test_align_empty_strings() {
        let r = align("", "");
        assert_eq!(r.edit_distance, 0);
        assert_eq!(r.alignment_detail, "");

        let r2 = align("", "ABC");
        assert_eq!(r2.edit_distance, 3);

        let r3 = align("ABC", "");
        assert_eq!(r3.edit_distance, 3);
    }
}
