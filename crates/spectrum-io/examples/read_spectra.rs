//! Quick CLI tool to test spectrum-io with real data files.
//!
//! Usage: cargo run -p protein-copilot-spectrum-io --example read_spectra -- <file>

use std::env;
use std::path::Path;

use protein_copilot_spectrum_io::{create_reader, detect_format};

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: cargo run -p protein-copilot-spectrum-io --example read_spectra -- <file.mgf|file.mzML>");
        std::process::exit(1);
    }

    let path = Path::new(&args[1]);
    println!("=== ProteinCopilot Spectrum Reader ===\n");

    // Step 1: detect format
    let info = match detect_format(path) {
        Ok(info) => {
            println!("File:   {}", info.path);
            println!("Format: {}", info.format);
            println!(
                "Size:   {:.2} MB",
                info.file_size_bytes as f64 / 1_048_576.0
            );
            info
        }
        Err(e) => {
            eprintln!("Error: {e}");
            std::process::exit(1);
        }
    };

    let reader = create_reader(&info);

    // Step 2: streaming summary
    println!("\n--- Summary (streaming) ---");
    match reader.read_summary(path) {
        Ok(s) => {
            println!("Total spectra:    {}", s.total_spectra);
            println!("MS1 count:        {}", s.ms1_count);
            println!("MS2 count:        {}", s.ms2_count);
            println!(
                "m/z range:        {:.2} - {:.2} Da",
                s.mz_range.0, s.mz_range.1
            );
            println!(
                "RT range:         {:.1} - {:.1} sec",
                s.rt_range_sec.0, s.rt_range_sec.1
            );
            println!("Median peaks:     {}", s.median_peaks_per_spectrum);

            let mut charges: Vec<_> = s.precursor_charge_distribution.iter().collect();
            charges.sort_by_key(|(c, _)| *c);
            println!("Charge distribution:");
            for (charge, count) in charges {
                println!("  charge {charge:+}: {count}");
            }

            match s.validate() {
                Ok(()) => println!("✓ Summary validation passed"),
                Err(e) => eprintln!("⚠ Validation failed: {e}"),
            }
        }
        Err(e) => eprintln!("Error: {e}"),
    }

    // Step 3: first 3 spectra
    println!("\n--- Sample spectra ---");
    for scan in 1..=3 {
        match reader.read_spectrum(path, scan) {
            Ok(sp) => {
                print!(
                    "Scan {}: {} peaks, RT={:.1}s",
                    sp.scan_number,
                    sp.num_peaks(),
                    sp.retention_time_sec
                );
                if let Some(p) = sp.precursors.first() {
                    print!(", precursor m/z={:.4}, charge={:?}", p.mz, p.charge);
                }
                println!();
            }
            Err(e) => eprintln!("Scan {scan}: {e}"),
        }
    }

    println!("\nDone!");
}
