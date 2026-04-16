//! Parsimony (greedy set cover) algorithm for protein inference.
//!
//! Given a peptide-to-protein mapping, finds the minimal set of proteins
//! that explains all observed peptides. Groups indistinguishable proteins
//! (those with identical peptide sets) into [`ProteinGroup`]s.
