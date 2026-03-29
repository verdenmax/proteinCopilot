//! Amino acid chemistry constants shared across the search engine.
//!
//! Contains monoisotopic residue masses and mass calculation utilities.

/// Monoisotopic mass of a single amino acid residue (Da).
///
/// Returns 0.0 for unknown residues.
pub fn residue_mass(aa: char) -> f64 {
    match aa {
        'G' => 57.021464,
        'A' => 71.037114,
        'V' => 99.068414,
        'L' => 113.084064,
        'I' => 113.084064,
        'P' => 97.052764,
        'F' => 147.068414,
        'W' => 186.079313,
        'M' => 131.040485,
        'S' => 87.032028,
        'T' => 101.047679,
        'C' => 103.009185,
        'Y' => 163.063329,
        'H' => 137.058912,
        'D' => 115.026943,
        'E' => 129.042593,
        'N' => 114.042927,
        'Q' => 128.058578,
        'K' => 128.094963,
        'R' => 156.101111,
        _ => 0.0,
    }
}

/// Water mass: H₂O added for intact peptide (N-term H + C-term OH).
pub const WATER_MASS: f64 = 18.010565;

/// Proton mass (Da).
pub const PROTON_MASS: f64 = 1.007276;

/// Calculates the monoisotopic neutral mass of a peptide sequence.
pub fn peptide_mass(sequence: &str) -> f64 {
    let sum: f64 = sequence.chars().map(residue_mass).sum();
    sum + WATER_MASS
}

/// Calculates the m/z value for a peptide at a given charge state.
///
/// Formula: m/z = (neutral_mass + charge × proton_mass) / charge
pub fn peptide_mz(neutral_mass: f64, charge: i32) -> f64 {
    (neutral_mass + charge as f64 * PROTON_MASS) / charge as f64
}
