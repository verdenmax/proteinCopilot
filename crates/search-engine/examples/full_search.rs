//! Full pipeline demo: spectrum file → recommend params → search → report
//!
//! Usage:
//!   cargo run -p protein-copilot-search-engine --example full_search -- <spectrum.mgf|mzML> <database.fasta> [output_dir]
//!
//! Example:
//!   cargo run -p protein-copilot-search-engine --example full_search -- data/sample.mgf data/human.fasta ./output

use std::env;
use std::path::Path;

use protein_copilot_core::engine::SearchEngineAdapter;
use protein_copilot_core::progress::noop_progress;
use protein_copilot_param_recommend::ParamRecommender;
use protein_copilot_report::ReportGenerator;
use protein_copilot_search_engine::SimpleSearchEngine;
use protein_copilot_spectrum_io::{create_reader, detect_format};

#[tokio::main]
async fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: cargo run -p protein-copilot-search-engine --example full_search -- <spectrum_file> <fasta_db>"
        );
        eprintln!("  spectrum_file: .mgf or .mzML file");
        eprintln!("  fasta_db:      .fasta protein database");
        std::process::exit(1);
    }

    let spectrum_path = Path::new(&args[1]);
    let fasta_path = &args[2];

    println!("╔══════════════════════════════════════════════╗");
    println!("║   ProteinCopilot — Full Search Pipeline      ║");
    println!("╚══════════════════════════════════════════════╝\n");

    // ──────────────────────────────────────────────────────
    // Step 1: Read spectrum file → SpectrumSummary
    // ──────────────────────────────────────────────────────
    println!("▶ Step 1: Reading spectrum file...");
    let file_info = match detect_format(spectrum_path) {
        Ok(info) => {
            println!(
                "  File:   {} ({}, {:.2} MB)",
                info.path,
                info.format,
                info.file_size_bytes as f64 / 1_048_576.0
            );
            info
        }
        Err(e) => {
            eprintln!("  ✗ Error: {e}");
            std::process::exit(1);
        }
    };

    let reader = create_reader(&file_info);
    let summary = match reader.read_summary(spectrum_path) {
        Ok(s) => {
            println!("  Total spectra:  {}", s.total_spectra);
            println!("  MS1 / MS2:      {} / {}", s.ms1_count, s.ms2_count);
            println!(
                "  m/z range:      {:.1} - {:.1}",
                s.mz_range[0], s.mz_range[1]
            );
            println!(
                "  RT range:       {:.1} - {:.1} sec",
                s.rt_range_min[0], s.rt_range_min[1]
            );
            println!("  Median peaks:   {}", s.median_peaks_per_spectrum);
            s
        }
        Err(e) => {
            eprintln!("  ✗ Error: {e}");
            std::process::exit(1);
        }
    };

    // ──────────────────────────────────────────────────────
    // Step 2: Recommend search parameters
    // ──────────────────────────────────────────────────────
    println!("\n▶ Step 2: Recommending search parameters...");
    let recommender = ParamRecommender;
    let recommendation = match recommender.recommend(&summary, None) {
        Ok(r) => {
            println!("  Enzyme:         {:?}", r.decision.enzyme);
            println!(
                "  Precursor tol:  {} {:?}",
                r.decision.precursor_tolerance.value, r.decision.precursor_tolerance.unit
            );
            println!(
                "  Fragment tol:   {} {:?}",
                r.decision.fragment_tolerance.value, r.decision.fragment_tolerance.unit
            );
            println!(
                "  Fixed mods:     {:?}",
                r.decision
                    .fixed_modifications
                    .iter()
                    .map(|m| m.name.as_str())
                    .collect::<Vec<_>>()
            );
            println!(
                "  Variable mods:  {:?}",
                r.decision
                    .variable_modifications
                    .iter()
                    .map(|m| m.name.as_str())
                    .collect::<Vec<_>>()
            );
            println!("  Confidence:     {:.2}", r.confidence);
            println!(
                "  Explanation:    {}",
                &r.explanation[..r.explanation.len().min(100)]
            );
            r
        }
        Err(e) => {
            eprintln!("  ✗ Error: {e}");
            std::process::exit(1);
        }
    };

    // ──────────────────────────────────────────────────────
    // Step 3: Set database path and run search
    // ──────────────────────────────────────────────────────
    println!("\n▶ Step 3: Running search...");
    let mut params = recommendation.decision;
    params.database_path = fasta_path.to_string();

    let engine = SimpleSearchEngine::new();
    println!("  Engine:         {}", engine.engine_info().name);
    println!("  Database:       {fasta_path}");

    let mut diag = protein_copilot_core::diagnostics::SearchDiagnostics::new();
    let result = match engine
        .search(
            &params,
            &[spectrum_path.to_path_buf()],
            noop_progress(),
            &mut diag,
        )
        .await
    {
        Ok(r) => r,
        Err(e) => {
            eprintln!("  ✗ Search failed: {e}");
            std::process::exit(1);
        }
    };

    // ──────────────────────────────────────────────────────
    // Step 4: Display results
    // ──────────────────────────────────────────────────────
    println!("\n▶ Step 4: Search Results");
    println!("  ─────────────────────────────────────");
    println!("  Run ID:              {}", result.run_id);
    println!(
        "  Duration:            {:.2} sec",
        result.summary.search_duration_sec
    );
    println!(
        "  Spectra searched:    {}",
        result.summary.total_spectra_searched
    );
    println!("  Total PSMs:          {}", result.summary.total_psms);
    println!(
        "  Identification rate: {:.1}%",
        result.summary.identification_rate * 100.0
    );
    println!(
        "  Unique peptides:     {}",
        result.summary.unique_peptides_at_1pct_fdr
    );
    println!(
        "  Protein groups:      {}",
        result.summary.protein_groups_at_1pct_fdr
    );
    println!("  Median score:        {:.4}", result.summary.median_score);
    println!(
        "  Median Δppm:         {:.2}",
        result.summary.median_delta_mass_ppm
    );

    // Top PSMs
    if !result.psms.is_empty() {
        println!("\n  Top PSMs (up to 5):");
        for psm in result.psms.iter().take(5) {
            println!(
                "    Scan {:>4} | {} | charge={} | score={:.4} | Δ={:.1}ppm | {}",
                psm.spectrum_scan,
                psm.peptide_sequence,
                psm.charge,
                psm.score,
                psm.delta_mass_ppm,
                psm.protein_accessions.join(",")
            );
        }
    }

    // Proteins
    if !result.proteins.is_empty() {
        println!("\n  Proteins identified:");
        for prot in &result.proteins {
            println!(
                "    {} | {:.0}% coverage | {} peptides | {}",
                prot.accession,
                prot.coverage * 100.0,
                prot.peptide_count,
                &prot.description[..prot.description.len().min(50)]
            );
        }
    }

    // ──────────────────────────────────────────────────────
    // Step 5: Generate report & export
    // ──────────────────────────────────────────────────────
    let output_dir = if args.len() >= 4 {
        Path::new(&args[3]).to_path_buf()
    } else {
        Path::new("./output").to_path_buf()
    };

    println!("\n▶ Step 5: Generating report...");
    println!("  Output dir:  {}", output_dir.display());

    // Re-generate summary with proper FDR filtering
    let summary = ReportGenerator::generate_summary(&result);
    println!(
        "  Summary:     {} PSMs at 1% FDR, {:.1}% identification rate",
        summary.psms_at_1pct_fdr,
        summary.identification_rate * 100.0
    );

    // Export TSV files
    match ReportGenerator::export_tsv(&result, &output_dir) {
        Ok(()) => {
            println!("  ✓ Exported:  psm.tsv, peptide.tsv, protein.tsv");
        }
        Err(e) => eprintln!("  ✗ TSV export error: {e}"),
    }

    // Export JSON
    let json_path = output_dir.join("result.json");
    match ReportGenerator::export_json(&result, &json_path) {
        Ok(()) => println!("  ✓ Exported:  result.json"),
        Err(e) => eprintln!("  ✗ JSON export error: {e}"),
    }

    // Export metadata
    let meta_path = output_dir.join("run_metadata.json");
    match ReportGenerator::export_metadata(&result.metadata, &meta_path) {
        Ok(()) => println!("  ✓ Exported:  run_metadata.json"),
        Err(e) => eprintln!("  ✗ Metadata export error: {e}"),
    }

    println!(
        "\n✓ Pipeline complete! Results in: {}",
        output_dir.display()
    );
}
