//! Decoy database generation.
//!
//! Generates decoy protein sequences for target-decoy FDR estimation.
//! Supports reverse and shuffle strategies.

use protein_copilot_core::search_params::DecoyStrategy;

/// A decoy protein entry.
#[derive(Debug, Clone)]
pub struct DecoyProtein {
    /// Accession prefixed with "REV_" or "SHUF_".
    pub accession: String,
    /// Original protein description.
    pub description: String,
    /// Decoy sequence (reversed or shuffled).
    pub sequence: String,
}

/// Generates decoy proteins from target sequences.
///
/// - `Reverse`: reverses each protein sequence but keeps the last amino acid
///   (C-terminal K/R) in place for tryptic digestion compatibility.
/// - `Shuffle`: randomly shuffles each protein sequence (deterministic with seed).
/// - `None`: returns empty vec.
pub fn generate_decoys(
    proteins: &[(String, String, String)], // (accession, description, sequence)
    strategy: DecoyStrategy,
) -> Vec<DecoyProtein> {
    match strategy {
        DecoyStrategy::None => Vec::new(),
        DecoyStrategy::Reverse => proteins
            .iter()
            .map(|(acc, desc, seq)| {
                let decoy_seq = reverse_sequence(seq);
                DecoyProtein {
                    accession: format!("REV_{acc}"),
                    description: desc.clone(),
                    sequence: decoy_seq,
                }
            })
            .collect(),
        DecoyStrategy::Shuffle => {
            use rand::seq::SliceRandom;
            use rand::SeedableRng;
            let mut rng = rand::rngs::StdRng::seed_from_u64(42);
            proteins
                .iter()
                .map(|(acc, desc, seq)| {
                    let mut chars: Vec<char> = seq.chars().collect();
                    if chars.len() > 1 {
                        let last = chars.len() - 1;
                        chars[..last].shuffle(&mut rng);
                    }
                    DecoyProtein {
                        accession: format!("SHUF_{acc}"),
                        description: desc.clone(),
                        sequence: chars.into_iter().collect(),
                    }
                })
                .collect()
        }
    }
}

/// Reverses a protein sequence, keeping the last amino acid in place.
///
/// This preserves C-terminal residues (K/R for trypsin) so decoy
/// peptides have similar enzymatic properties to target peptides.
fn reverse_sequence(seq: &str) -> String {
    let chars: Vec<char> = seq.chars().collect();
    if chars.len() <= 1 {
        return seq.to_string();
    }
    let last = chars[chars.len() - 1];
    let mut middle: Vec<char> = chars[..chars.len() - 1].to_vec();
    middle.reverse();
    middle.push(last);
    middle.into_iter().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_keeps_last_residue() {
        assert_eq!(reverse_sequence("PEPTIDEK"), "EDITPEPK");
    }

    #[test]
    fn reverse_single_char() {
        assert_eq!(reverse_sequence("K"), "K");
    }

    #[test]
    fn reverse_empty() {
        assert_eq!(reverse_sequence(""), "");
    }

    #[test]
    fn generate_decoys_reverse() {
        let proteins = vec![(
            "P001".to_string(),
            "Test".to_string(),
            "PEPTIDEK".to_string(),
        )];
        let decoys = generate_decoys(&proteins, DecoyStrategy::Reverse);
        assert_eq!(decoys.len(), 1);
        assert_eq!(decoys[0].accession, "REV_P001");
        assert_eq!(decoys[0].sequence, "EDITPEPK");
    }

    #[test]
    fn generate_decoys_none() {
        let proteins = vec![(
            "P001".to_string(),
            "Test".to_string(),
            "PEPTIDEK".to_string(),
        )];
        let decoys = generate_decoys(&proteins, DecoyStrategy::None);
        assert!(decoys.is_empty());
    }

    #[test]
    fn generate_decoys_shuffle_deterministic() {
        let proteins = vec![(
            "P001".to_string(),
            "Test".to_string(),
            "PEPTIDEK".to_string(),
        )];
        let d1 = generate_decoys(&proteins, DecoyStrategy::Shuffle);
        let d2 = generate_decoys(&proteins, DecoyStrategy::Shuffle);
        assert_eq!(d1[0].sequence, d2[0].sequence);
        assert_eq!(d1[0].sequence.chars().last(), Some('K'));
    }

    #[test]
    fn generate_decoys_shuffle_prefix() {
        let proteins = vec![(
            "P001".to_string(),
            "Test".to_string(),
            "PEPTIDEK".to_string(),
        )];
        let decoys = generate_decoys(&proteins, DecoyStrategy::Shuffle);
        assert!(decoys[0].accession.starts_with("SHUF_"));
    }

    #[test]
    fn reverse_preserves_length() {
        let seq = "MKWVTFISLLFLFSSAYSRGVFRR";
        let rev = reverse_sequence(seq);
        assert_eq!(rev.len(), seq.len());
        assert_eq!(rev.chars().last(), Some('R'));
    }
}
