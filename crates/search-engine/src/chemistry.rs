//! Amino acid chemistry constants shared across the search engine.
//!
//! Contains monoisotopic residue masses and mass calculation utilities.

/// Monoisotopic mass of a single amino acid residue (Da).
///
/// Returns `None` for unknown/non-standard residues (e.g. `*`, `X`, `B`, `Z`).
pub fn residue_mass(aa: char) -> Option<f64> {
    match aa {
        'G' => Some(57.021464),
        'A' => Some(71.037114),
        'V' => Some(99.068414),
        'L' => Some(113.084064),
        'I' => Some(113.084064),
        'P' => Some(97.052764),
        'F' => Some(147.068414),
        'W' => Some(186.079313),
        'M' => Some(131.040485),
        'S' => Some(87.032028),
        'T' => Some(101.047679),
        'C' => Some(103.009185),
        'Y' => Some(163.063329),
        'H' => Some(137.058912),
        'D' => Some(115.026943),
        'E' => Some(129.042593),
        'N' => Some(114.042927),
        'Q' => Some(128.058578),
        'K' => Some(128.094963),
        'R' => Some(156.101111),
        _ => None,
    }
}

/// Returns `true` if every character in the sequence is a standard amino acid.
pub fn is_standard_sequence(sequence: &str) -> bool {
    sequence.chars().all(|c| residue_mass(c).is_some())
}

/// Water mass: H₂O added for intact peptide (N-term H + C-term OH).
pub const WATER_MASS: f64 = 18.010565;

/// Proton mass (Da).
pub const PROTON_MASS: f64 = 1.007276;

/// Mass loss for H₂O neutral loss from fragment ions (Da).
/// Loss of water from side-chain hydroxyl or carboxyl groups (-18.010565 Da).
pub const H2O_LOSS: f64 = 18.010565;

/// Mass loss for NH₃ neutral loss from fragment ions (Da).
/// Loss of ammonia from side-chain amine groups (-17.026549 Da).
pub const NH3_LOSS: f64 = 17.026549;

/// Calculates the monoisotopic neutral mass of a peptide sequence.
///
/// Returns `None` if the sequence contains any non-standard amino acid.
pub fn peptide_mass(sequence: &str) -> Option<f64> {
    let mut sum = 0.0;
    for aa in sequence.chars() {
        sum += residue_mass(aa)?;
    }
    Some(sum + WATER_MASS)
}

/// Calculates the m/z value for a peptide at a given charge state.
///
/// Formula: m/z = (neutral_mass + charge × proton_mass) / charge
///
/// # Panics
///
/// Panics if `charge` is zero or negative (physically impossible in MS).
pub fn peptide_mz(neutral_mass: f64, charge: i32) -> f64 {
    assert!(charge > 0, "charge must be > 0, got {charge}");
    (neutral_mass + charge as f64 * PROTON_MASS) / charge as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_residues_return_some() {
        for aa in "GAVLIPFWMSTCYHDENKQR".chars() {
            assert!(
                residue_mass(aa).is_some(),
                "Standard amino acid '{aa}' should return Some"
            );
        }
    }

    #[test]
    fn nonstandard_residues_return_none() {
        for aa in "OBJUXZ*".chars() {
            assert!(
                residue_mass(aa).is_none(),
                "Non-standard character '{aa}' should return None"
            );
        }
    }

    #[test]
    fn peptide_mass_with_nonstandard_returns_none() {
        assert!(peptide_mass("PEPTIDE").is_some());
        assert!(peptide_mass("PEPT*DE").is_none(), "stop codon");
        assert!(peptide_mass("PEPTXDE").is_none(), "unknown X");
        assert!(peptide_mass("PEPTODE").is_none(), "non-standard O");
    }

    #[test]
    fn is_standard_sequence_works() {
        assert!(is_standard_sequence("PEPTIDEK"));
        assert!(!is_standard_sequence("PEPT*DE"));
        assert!(!is_standard_sequence("PEPTODE"));
    }
}
