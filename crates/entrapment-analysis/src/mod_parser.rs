//! Parser for DIA-NN `Modified.Sequence` format.
//!
//! Strips UniMod annotations from peptide sequences and returns a list of
//! [`ParsedModification`] entries with 0-based positions, delta masses, and
//! UniMod accession numbers.

/// A single post-translational modification parsed from a modified sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct ParsedModification {
    /// 0-based residue index in the stripped (unmodified) sequence.
    pub position: usize,
    /// Mass delta in Da.
    pub delta_mass: f64,
    /// UniMod accession number.
    pub unimod_id: u32,
}

/// Look up the mono-isotopic delta mass for a UniMod accession ID.
///
/// Returns `None` for unknown IDs. Covers the 13 most common accessions
/// encountered in DIA-NN / pFind / MaxQuant output.
pub fn unimod_delta_mass(id: u32) -> Option<f64> {
    match id {
        1 => Some(42.010565),    // Acetyl (N-term)
        4 => Some(57.021464),    // Carbamidomethyl (C)
        5 => Some(43.005814),    // Carbamyl
        7 => Some(0.984016),     // Deamidated (N, Q)
        21 => Some(79.966331),   // Phospho (S, T, Y)
        27 => Some(-18.010565),  // Glu->pyro-Glu (loss of H₂O)
        28 => Some(-17.026549),  // Gln->pyro-Glu (loss of NH₃)
        34 => Some(14.015650),   // Methyl
        35 => Some(15.994915),   // Oxidation (M)
        121 => Some(114.042927), // GG (ubiquitin remnant, K)
        214 => Some(229.162932), // TMT6plex / TMTpro
        259 => Some(8.014199),   // Label:13C(6)15N(2) heavy K (SILAC)
        267 => Some(10.008269),  // Label:13C(6)15N(4) heavy R (SILAC)
        _ => None,
    }
}

/// Parse a DIA-NN `Modified.Sequence` string into a stripped sequence and
/// a vector of modifications.
///
/// The format consists of amino-acid residues optionally followed by
/// `(UniMod:N)` annotations.  N-terminal modifications appear as
/// `(UniMod:N)` at the very start, *before* any residue character,
/// and are assigned position 0 (the first residue).
///
/// # Examples
///
/// ```
/// use protein_copilot_entrapment_analysis::mod_parser::parse_modified_sequence;
///
/// let (seq, mods) = parse_modified_sequence("AAAC(UniMod:4)DFK");
/// assert_eq!(seq, "AAACDFK");
/// assert_eq!(mods.len(), 1);
/// assert_eq!(mods[0].position, 3);
/// ```
pub fn parse_modified_sequence(modified_seq: &str) -> (String, Vec<ParsedModification>) {
    if modified_seq.is_empty() {
        return (String::new(), Vec::new());
    }

    let mut stripped = String::with_capacity(modified_seq.len());
    let mut mods = Vec::new();
    let bytes = modified_seq.as_bytes();
    let len = bytes.len();
    let mut i = 0;

    while i < len {
        if bytes[i] == b'(' {
            // Try to parse a "(UniMod:N)" annotation.
            if let Some(end) = find_closing_paren(bytes, i) {
                let inner = &modified_seq[i + 1..end]; // content between ( and )
                if let Some(id) = parse_unimod_tag(inner) {
                    let delta_mass = unimod_delta_mass(id).unwrap_or(0.0);
                    // Position is the index of the *preceding* residue.
                    // If the annotation appears at the start (N-terminal mod),
                    // position is 0.
                    let position = if stripped.is_empty() { 0 } else { stripped.len() - 1 };
                    mods.push(ParsedModification {
                        position,
                        delta_mass,
                        unimod_id: id,
                    });
                }
                i = end + 1; // skip past ')'
            } else {
                // Malformed — no closing paren; treat '(' as a regular char.
                stripped.push('(');
                i += 1;
            }
        } else {
            stripped.push(bytes[i] as char);
            i += 1;
        }
    }

    (stripped, mods)
}

/// Find the index of the closing ')' that matches the '(' at `start`.
fn find_closing_paren(bytes: &[u8], start: usize) -> Option<usize> {
    let mut depth = 0u32;
    for (j, &b) in bytes.iter().enumerate().skip(start) {
        match b {
            b'(' => depth += 1,
            b')' => {
                depth -= 1;
                if depth == 0 {
                    return Some(j);
                }
            }
            _ => {}
        }
    }
    None
}

/// Parse `"UniMod:N"` → `Some(N)`.
fn parse_unimod_tag(inner: &str) -> Option<u32> {
    let stripped = inner.trim();
    let id_str = stripped.strip_prefix("UniMod:")?;
    id_str.parse::<u32>().ok()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_unmodified_sequence() {
        let (seq, mods) = parse_modified_sequence("PEPTIDE");
        assert_eq!(seq, "PEPTIDE");
        assert!(mods.is_empty());
    }

    #[test]
    fn parse_single_carbamidomethyl() {
        let (seq, mods) = parse_modified_sequence("AAAC(UniMod:4)DFK");
        assert_eq!(seq, "AAACDFK");
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].position, 3);
        assert!((mods[0].delta_mass - 57.021464).abs() < 1e-6);
        assert_eq!(mods[0].unimod_id, 4);
    }

    #[test]
    fn parse_oxidation() {
        let (seq, mods) = parse_modified_sequence("PEPTM(UniMod:35)DE");
        assert_eq!(seq, "PEPTMDE");
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].position, 4);
        assert!((mods[0].delta_mass - 15.994915).abs() < 1e-6);
        assert_eq!(mods[0].unimod_id, 35);
    }

    #[test]
    fn parse_multiple_modifications() {
        let (seq, mods) = parse_modified_sequence("AC(UniMod:4)DEFM(UniMod:35)GK");
        assert_eq!(seq, "ACDEFMGK");
        assert_eq!(mods.len(), 2);
        assert_eq!(mods[0].position, 1);
        assert_eq!(mods[0].unimod_id, 4);
        assert_eq!(mods[1].position, 5);
        assert_eq!(mods[1].unimod_id, 35);
    }

    #[test]
    fn parse_nterm_modification() {
        let (seq, mods) = parse_modified_sequence("(UniMod:1)PEPTIDE");
        assert_eq!(seq, "PEPTIDE");
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].position, 0);
        assert!((mods[0].delta_mass - 42.010565).abs() < 1e-6);
        assert_eq!(mods[0].unimod_id, 1);
    }

    #[test]
    fn parse_empty_string() {
        let (seq, mods) = parse_modified_sequence("");
        assert_eq!(seq, "");
        assert!(mods.is_empty());
    }

    #[test]
    fn parse_unknown_unimod_id() {
        let (seq, mods) = parse_modified_sequence("PEP(UniMod:99999)TIDE");
        assert_eq!(seq, "PEPTIDE");
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].position, 2);
        assert!((mods[0].delta_mass - 0.0).abs() < 1e-6);
        assert_eq!(mods[0].unimod_id, 99999);
    }
}
