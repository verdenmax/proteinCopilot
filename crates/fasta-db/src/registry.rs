//! Built-in database registry — static definitions of common proteomics FASTA databases.

/// A built-in database entry in the registry.
#[derive(Debug, Clone)]
pub struct DatabaseEntry {
    pub id: &'static str,
    pub species: &'static str,
    pub taxonomy_id: u32,
    pub db_type: &'static str,
    pub description: &'static str,
    pub url: &'static str,
}

const BUILTIN_DATABASES: &[DatabaseEntry] = &[
    DatabaseEntry {
        id: "human_swissprot",
        species: "Homo sapiens",
        taxonomy_id: 9606,
        db_type: "Swiss-Prot",
        description: "Human reviewed proteome (UniProt Swiss-Prot)",
        url: "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(reviewed:true)+AND+(model_organism:9606)",
    },
    DatabaseEntry {
        id: "mouse_swissprot",
        species: "Mus musculus",
        taxonomy_id: 10090,
        db_type: "Swiss-Prot",
        description: "Mouse reviewed proteome (UniProt Swiss-Prot)",
        url: "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(reviewed:true)+AND+(model_organism:10090)",
    },
    DatabaseEntry {
        id: "ecoli_swissprot",
        species: "Escherichia coli (K12)",
        taxonomy_id: 83333,
        db_type: "Swiss-Prot",
        description: "E. coli K12 reviewed proteome (UniProt Swiss-Prot)",
        url: "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(reviewed:true)+AND+(model_organism:83333)",
    },
    DatabaseEntry {
        id: "yeast_swissprot",
        species: "Saccharomyces cerevisiae",
        taxonomy_id: 559292,
        db_type: "Swiss-Prot",
        description: "Yeast reviewed proteome (UniProt Swiss-Prot)",
        url: "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(reviewed:true)+AND+(model_organism:559292)",
    },
    DatabaseEntry {
        id: "arabidopsis_swissprot",
        species: "Arabidopsis thaliana",
        taxonomy_id: 3702,
        db_type: "Swiss-Prot",
        description: "Arabidopsis reviewed proteome (UniProt Swiss-Prot)",
        url: "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(reviewed:true)+AND+(model_organism:3702)",
    },
    DatabaseEntry {
        id: "crap",
        species: "Contaminants",
        taxonomy_id: 0,
        db_type: "cRAP",
        description: "Common Repository of Adventitious Proteins (contaminants)",
        url: "https://ftp.thegpm.org/fasta/cRAP/crap.fasta",
    },
];

/// Returns all built-in database entries.
pub fn all_databases() -> &'static [DatabaseEntry] {
    BUILTIN_DATABASES
}

/// Returns all built-in database IDs.
pub fn all_database_ids() -> Vec<&'static str> {
    BUILTIN_DATABASES.iter().map(|e| e.id).collect()
}

/// Looks up a database entry by ID. Returns `None` if not found.
pub fn get_database(id: &str) -> Option<&'static DatabaseEntry> {
    BUILTIN_DATABASES.iter().find(|e| e.id == id)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registry_has_expected_entries() {
        let entries = all_databases();
        assert!(
            entries.len() >= 6,
            "should have at least 6 built-in databases"
        );

        let human = get_database("human_swissprot");
        assert!(human.is_some(), "human_swissprot must exist");
        let human = human.unwrap();
        assert_eq!(human.taxonomy_id, 9606);
        assert!(human.url.contains("rest.uniprot.org"));

        let crap = get_database("crap");
        assert!(crap.is_some(), "crap must exist");
    }

    #[test]
    fn unknown_database_returns_none() {
        assert!(get_database("nonexistent_db").is_none());
    }

    #[test]
    fn all_ids_returns_complete_list() {
        let ids = all_database_ids();
        assert!(ids.contains(&"human_swissprot"));
        assert!(ids.contains(&"crap"));
        assert_eq!(ids.len(), all_databases().len());
    }
}
