use protein_copilot_core::spectrum::*;

fn main() {
    let result = Spectrum::new(
        1,
        MsLevel::MS2,
        10.0,
        vec![PrecursorInfo {
            mz: 500.0,
            charge: None,
            intensity: None,
            isolation_window: Some(IsolationWindow {
                target_mz: 500.0,
                lower_offset: -1.0,
                upper_offset: 12.5,
            }),
            source_scan: None,
        }],
        vec![100.0],
        vec![1000.0],
    );
    
    match result {
        Err(e) => println!("Error message: {}", e),
        Ok(_) => println!("Unexpectedly succeeded"),
    }
}
