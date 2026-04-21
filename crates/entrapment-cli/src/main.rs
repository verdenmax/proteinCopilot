//! CLI for entrapment (trap-database) hit classification.
//!
//! Provides three subcommands:
//!
//! - **analyze** – run the full entrapment analysis pipeline
//! - **report** – regenerate an HTML report from a classified TSV
//! - **inspect** – inspect a single peptide's similarity to the target database

use std::path::{Path, PathBuf};
use std::process;

use clap::{Parser, Subcommand, ValueEnum};
use tracing::info;

use protein_copilot_entrapment_analysis::config::{EntrapmentConfig, SimilarityConfig};
use protein_copilot_entrapment_analysis::digest::TargetDigestIndex;
use protein_copilot_entrapment_analysis::loader::{self, ResultFormat};
use protein_copilot_entrapment_analysis::output::{
    file_sha256, write_classified_tsv, write_razor_errors_tsv, write_run_metadata, RunMetadata,
};
use protein_copilot_entrapment_analysis::report;
use protein_copilot_entrapment_analysis::similarity::classify_single;
use protein_copilot_entrapment_analysis::{EntrapmentAnalyzer, PsmGroup, UnifiedPsm};

// ---------------------------------------------------------------------------
// CLI definition
// ---------------------------------------------------------------------------

/// Entrapment analysis – classify trap PSMs by homology.
#[derive(Parser)]
#[command(
    name = "entrapment",
    about = "Entrapment analysis - classify trap PSMs by homology"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Run full entrapment analysis pipeline
    Analyze {
        /// Path to search results file (.parquet or .tsv)
        #[arg(short, long)]
        results: String,
        /// Path to YAML config file
        #[arg(short, long)]
        config: String,
        /// Path to target FASTA database
        #[arg(short, long)]
        target_fasta: String,
        /// Result format (auto-detect from extension if omitted)
        #[arg(short, long)]
        format: Option<FormatArg>,
        /// Output directory (default: ./output/entrapment)
        #[arg(short, long, default_value = "output/entrapment")]
        out: String,
        /// Directory containing mzML files for provenance tracing (optional)
        #[arg(long)]
        mzml_dir: Option<PathBuf>,
    },
    /// Regenerate HTML report from classified TSV
    Report {
        /// Path to classified TSV file
        #[arg(short = 'l', long)]
        classified: String,
        /// Output HTML path (default: entrapment_report.html in same dir)
        #[arg(short, long)]
        out: Option<String>,
    },
    /// Inspect a single peptide's similarity to target database
    Inspect {
        /// Peptide sequence to inspect
        #[arg(short, long)]
        peptide: String,
        /// Path to target FASTA database
        #[arg(short, long)]
        target_fasta: String,
        /// Optional YAML config for similarity settings
        #[arg(short, long)]
        config: Option<String>,
    },
}

/// Supported result-file formats for the CLI `--format` flag.
#[derive(Clone, ValueEnum)]
enum FormatArg {
    /// DIA-NN parquet format
    DiannParquet,
    /// Generic tab-separated values
    GenericTsv,
}

impl FormatArg {
    /// Convert to the library's [`ResultFormat`].
    fn to_result_format(&self) -> ResultFormat {
        match self {
            Self::DiannParquet => ResultFormat::DiannParquet,
            Self::GenericTsv => ResultFormat::GenericTsv,
        }
    }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

fn main() {
    // Initialise tracing subscriber (env-filter honours RUST_LOG)
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .init();

    let cli = Cli::parse();

    let result = match cli.command {
        Commands::Analyze {
            results,
            config,
            target_fasta,
            format,
            out,
            mzml_dir,
        } => run_analyze(&results, &config, &target_fasta, format.as_ref(), &out, mzml_dir.as_deref()),
        Commands::Report { classified, out } => run_report(&classified, out.as_deref()),
        Commands::Inspect {
            peptide,
            target_fasta,
            config,
        } => run_inspect(&peptide, &target_fasta, config.as_deref()),
    };

    if let Err(e) = result {
        eprintln!("Error: {e}");
        process::exit(1);
    }
}

// ---------------------------------------------------------------------------
// Subcommand implementations
// ---------------------------------------------------------------------------

/// Run the full entrapment analysis pipeline.
fn run_analyze(
    results_path: &str,
    config_path: &str,
    fasta_path: &str,
    format_arg: Option<&FormatArg>,
    out_dir: &str,
    mzml_dir: Option<&Path>,
) -> Result<(), Box<dyn std::error::Error>> {
    let results_path = Path::new(results_path);
    let config_path = Path::new(config_path);
    let fasta_path = Path::new(fasta_path);
    let out_dir = PathBuf::from(out_dir);

    // 1. Load config from YAML
    info!(path = %config_path.display(), "loading config");
    let config = EntrapmentConfig::from_yaml(config_path)?;

    // 2. Detect or use specified format
    let format = match format_arg {
        Some(fa) => fa.to_result_format(),
        None => ResultFormat::from_path(results_path)?,
    };

    // 3. Load PSMs
    info!(path = %results_path.display(), "loading PSMs");
    let psms = loader::load_psms(results_path, &format, None)?;
    info!(count = psms.len(), "loaded PSMs");

    // 4. Build analyser
    info!(fasta = %fasta_path.display(), "building entrapment analyser");
    let analyser = EntrapmentAnalyzer::new(config.clone(), fasta_path)?;

    // 5. Classify all PSMs
    info!("classifying PSMs");
    let mut classified = analyser.classify_all(&psms)?;

    // 5b. Provenance tracing (optional)
    if let Some(mzml_dir) = mzml_dir {
        use protein_copilot_entrapment_analysis::trace_provenance_batch;

        println!(
            "Running provenance tracing with mzML files from: {}",
            mzml_dir.display()
        );
        match trace_provenance_batch(&mut classified, mzml_dir, &config) {
            Ok(count) => println!("Provenance traced for {} PSMs", count),
            Err(e) => eprintln!("Warning: provenance tracing failed: {}", e),
        }
    }

    // 6. Create output directory
    std::fs::create_dir_all(&out_dir)?;

    // 7. Write outputs
    let classified_tsv_path = out_dir.join("classified.tsv");
    write_classified_tsv(&classified, &classified_tsv_path)?;

    let razor_tsv_path = out_dir.join("razor_errors.tsv");
    write_razor_errors_tsv(&classified, &razor_tsv_path)?;

    // Build run metadata
    let summary = analyser.summary(&classified);

    let config_snapshot = serde_json::to_value(&config).unwrap_or(serde_json::Value::Null);

    let metadata = RunMetadata {
        tool_version: env!("CARGO_PKG_VERSION").to_owned(),
        run_timestamp: chrono::Utc::now().to_rfc3339(),
        input_file: results_path.display().to_string(),
        input_sha256: file_sha256(results_path)?,
        fasta_file: fasta_path.display().to_string(),
        fasta_sha256: file_sha256(fasta_path)?,
        config_snapshot,
        total_psms: summary.total_psms,
        trap_psms: summary.trap_psms,
        level_counts: summary.level_counts.clone(),
    };

    let metadata_path = out_dir.join("run_metadata.json");
    write_run_metadata(&metadata, &metadata_path)?;

    // 8. Generate HTML report
    let report_path = out_dir.join("entrapment_report.html");
    report::render_report(&summary, &classified, &report_path)?;
    info!(path = %report_path.display(), "HTML report generated");

    // 9. Print summary to stdout
    println!("=== Entrapment Analysis Summary ===");
    println!("Total PSMs:     {}", summary.total_psms);
    println!("  Target:       {}", summary.target_psms);
    println!("  Trap:         {}", summary.trap_psms);
    println!("  Ambiguous:    {}", summary.ambiguous_psms);
    println!();
    println!("Trap PSM breakdown by discriminability level:");
    println!("  L0 (razor error):         {}", summary.level_counts.l0);
    println!("  L1 (L/I isomer):          {}", summary.level_counts.l1);
    println!("  L2 (near-isobaric):       {}", summary.level_counts.l2);
    println!("  L3 (distinguishable):     {}", summary.level_counts.l3);
    println!("  L4 (true trap):           {}", summary.level_counts.l4);
    println!();
    println!("Output written to: {}", out_dir.display());

    Ok(())
}

/// Regenerate an HTML report from a classified TSV.
///
/// Not yet implemented – placeholder for Task 10.
fn run_report(
    _classified_path: &str,
    _out_path: Option<&str>,
) -> Result<(), Box<dyn std::error::Error>> {
    Err("Report regeneration not yet implemented (see Task 10). Use 'analyze' subcommand which generates the HTML report automatically.".into())
}

/// Inspect a single peptide's similarity to the target database.
fn run_inspect(
    peptide: &str,
    fasta_path: &str,
    config_path: Option<&str>,
) -> Result<(), Box<dyn std::error::Error>> {
    let fasta_path = Path::new(fasta_path);

    // 1. Build default or loaded SimilarityConfig
    let similarity_config = match config_path {
        Some(path) => {
            let cfg = EntrapmentConfig::from_yaml(Path::new(path))?;
            cfg.similarity
        }
        None => SimilarityConfig::default(),
    };

    // 2. Build TargetDigestIndex from FASTA
    info!(fasta = %fasta_path.display(), "building target digest index");
    let index = TargetDigestIndex::from_fasta(
        fasta_path,
        similarity_config.max_missed_cleavages,
        similarity_config.max_mismatches,
    )?;
    info!(peptides = index.len(), "target digest index ready");

    // 3. Create dummy UnifiedPsm with the peptide
    let psm = UnifiedPsm {
        peptide: peptide.to_owned(),
        charge: None,
        precursor_mz: None,
        retention_time: None,
        scan_number: None,
        spectrum_file: None,
        protein_ids: String::new(),
        q_value: None,
        modifications: Vec::new(),
    };

    // 4. Classify the single peptide (always as Trap to trigger similarity checks)
    let result = classify_single(&psm, PsmGroup::Trap, &index, &similarity_config);

    // 5. Print results
    println!("=== Peptide Inspection ===");
    println!("Peptide:           {peptide}");
    println!("Level:             {}", result.level);
    println!(
        "Best target match: {}",
        result.best_target_peptide.as_deref().unwrap_or("(none)")
    );
    println!(
        "Best target prot:  {}",
        result.best_target_protein.as_deref().unwrap_or("(none)")
    );
    println!(
        "Mismatches:        {}",
        result
            .mismatches
            .map(|m| m.to_string())
            .unwrap_or_else(|| "N/A".to_owned())
    );
    println!(
        "Delta mass (Da):   {}",
        result
            .delta_mass_da
            .map(|d| format!("{d:.6}"))
            .unwrap_or_else(|| "N/A".to_owned())
    );
    println!(
        "Diff positions:    {}",
        result.diff_positions.as_deref().unwrap_or("N/A")
    );

    Ok(())
}
