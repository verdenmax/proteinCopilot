//! Razor peptide assignment.
//!
//! Shared peptides (those mapping to multiple protein groups) are assigned
//! to the group with the most unique peptide evidence ("razor" logic).
//! This maximizes quantitative accuracy by avoiding double-counting.
