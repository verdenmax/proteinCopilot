//! # ProteinCopilot Core
//!
//! Shared domain types and traits for the ProteinCopilot platform.
//!
//! This crate defines the data structures used across all MCP modules:
//! - Spectrum data (mass spectrometry spectra and summaries)
//! - Search parameters (enzyme, modifications, tolerances)
//! - Search results (PSM, peptide, protein level)
//! - AI decision wrapper (structured LLM output)
//! - Search engine adapter trait
//! - Run metadata and error types

pub mod ai_decision;
pub mod engine;
pub mod error;
pub mod label;
pub mod progress;
pub mod run_metadata;
pub mod search_params;
pub mod search_result;
pub mod spectrum;
pub mod protein_group;
pub mod util;
