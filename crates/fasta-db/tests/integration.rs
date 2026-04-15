//! Integration tests for fasta-db (requires network access).
//! Run with: cargo test -p protein-copilot-fasta-db -- --ignored

use protein_copilot_fasta_db::{download_database, get_database_info, list_databases, DownloadStatus};
use tempfile::TempDir;

#[tokio::test]
#[ignore = "requires network access"]
async fn download_ecoli_end_to_end() {
    let dir = TempDir::new().unwrap();
    let cache_dir = dir.path();

    // Step 1: List — ecoli should be Available
    let list = list_databases(cache_dir).unwrap();
    let ecoli = list.iter().find(|d| d.id == "ecoli_swissprot").unwrap();
    assert!(
        matches!(ecoli.status, DownloadStatus::Available),
        "ecoli should not be cached yet"
    );

    // Step 2: Download
    let result = download_database("ecoli_swissprot", cache_dir, false)
        .await
        .unwrap();
    assert!(result.protein_count > 100, "E.coli should have >100 proteins");
    assert!(result.file_size_bytes > 0);

    // Step 3: List — ecoli should now be Downloaded
    let list = list_databases(cache_dir).unwrap();
    let ecoli = list.iter().find(|d| d.id == "ecoli_swissprot").unwrap();
    assert!(
        matches!(ecoli.status, DownloadStatus::Downloaded { .. }),
        "ecoli should be cached now"
    );

    // Step 4: Info
    let info = get_database_info("ecoli_swissprot", cache_dir).unwrap();
    assert_eq!(info.protein_count, result.protein_count);
    assert!(!info.first_accessions.is_empty());
    assert!(info.sha256.len() == 64);

    // Step 5: Re-download with force=false — should return cached
    let result2 = download_database("ecoli_swissprot", cache_dir, false)
        .await
        .unwrap();
    assert_eq!(result2.protein_count, result.protein_count);

    // Step 6: Re-download with force=true — should re-download
    let result3 = download_database("ecoli_swissprot", cache_dir, true)
        .await
        .unwrap();
    assert!(result3.protein_count > 0);
}

#[tokio::test]
#[ignore = "requires network access"]
async fn unknown_database_returns_error() {
    let dir = TempDir::new().unwrap();
    let result = download_database("nonexistent_db", dir.path(), false).await;
    assert!(result.is_err());
    let err = result.unwrap_err().to_string();
    assert!(err.contains("unknown database"));
    assert!(err.contains("human_swissprot")); // should list available
}
