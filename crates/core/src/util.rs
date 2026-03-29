//! Shared utility functions used across ProteinCopilot crates.

/// Computes the median of a sorted f64 slice.
///
/// - Returns `0.0` for empty slices.
/// - For odd-length slices, returns the middle element.
/// - For even-length slices, returns the average of the two middle elements.
///
/// The input slice **must be sorted** before calling this function.
pub fn compute_median(sorted: &[f64]) -> f64 {
    let len = sorted.len();
    if len == 0 {
        0.0
    } else if len % 2 == 0 {
        (sorted[len / 2 - 1] + sorted[len / 2]) / 2.0
    } else {
        sorted[len / 2]
    }
}

/// Computes the median of a sorted u32 slice.
///
/// Same semantics as [`compute_median`] but for integer arrays.
/// Returns the lower-middle value for even-length slices (integer division).
pub fn compute_median_u32(sorted: &[u32]) -> u32 {
    let len = sorted.len();
    if len == 0 {
        0
    } else {
        sorted[len / 2]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_empty() {
        assert_eq!(compute_median(&[]), 0.0);
    }

    #[test]
    fn median_single() {
        assert_eq!(compute_median(&[5.0]), 5.0);
    }

    #[test]
    fn median_odd() {
        assert_eq!(compute_median(&[1.0, 3.0, 5.0]), 3.0);
    }

    #[test]
    fn median_even() {
        assert!((compute_median(&[1.0, 3.0]) - 2.0).abs() < 0.001);
        assert!((compute_median(&[0.6, 0.8]) - 0.7).abs() < 0.001);
    }

    #[test]
    fn median_u32_empty() {
        assert_eq!(compute_median_u32(&[]), 0);
    }

    #[test]
    fn median_u32_odd() {
        assert_eq!(compute_median_u32(&[1, 3, 5]), 3);
    }
}
