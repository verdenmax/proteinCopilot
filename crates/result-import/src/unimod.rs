//! Unimod modification database: maps record_id → name + mono_mass.
//!
//! Provides a builtin table of ~20 common modifications and an XML parser
//! for the full unimod.xml database.

use std::collections::HashMap;
use std::path::Path;

use protein_copilot_core::search_params::{ModPosition, Modification};

use crate::ResultImportError;

/// A single Unimod modification entry.
#[derive(Debug, Clone)]
pub struct UnimodEntry {
    pub record_id: u32,
    pub title: String,
    pub mono_mass: f64,
    /// Residues this modification can occur on (empty = any).
    pub residues: Vec<char>,
}

/// Unimod modification database.
pub struct UnimodDb {
    entries: HashMap<u32, UnimodEntry>,
}

impl UnimodDb {
    /// Create a database with ~20 common builtin modifications.
    pub fn builtin() -> Self {
        let mut entries = HashMap::new();
        let common = vec![
            (1, "Acetyl", 42.010565, vec![]),
            (4, "Carbamidomethyl", 57.021464, vec!['C']),
            (5, "Carbamyl", 43.005814, vec![]),
            (7, "Deamidated", 0.984016, vec!['N', 'Q']),
            (21, "Phospho", 79.966331, vec!['S', 'T', 'Y']),
            (28, "Gln->pyro-Glu", -17.026549, vec!['Q']),
            (27, "Glu->pyro-Glu", -18.010565, vec!['E']),
            (34, "Methyl", 14.015650, vec![]),
            (35, "Oxidation", 15.994915, vec!['M', 'W', 'H']),
            (36, "Dimethyl", 28.031300, vec![]),
            (37, "Trimethyl", 42.046950, vec![]),
            (39, "Dehydrated", -18.010565, vec![]),
            (40, "Formyl", 27.994915, vec![]),
            (121, "GlyGly", 114.042927, vec!['K']),
            (188, "Label:13C(6)", 6.020129, vec!['K', 'R']),
            (199, "Label:13C(6)15N(2)", 8.014199, vec!['K']),
            (259, "Label:13C(6)15N(4)", 10.008269, vec!['R']),
            (214, "Label:2H(4)", 4.025107, vec![]),
            (267, "Label:13C(6)15N(1)", 7.017165, vec![]),
            (268, "iTRAQ4plex", 144.102063, vec![]),
            (737, "TMT6plex", 229.162932, vec![]),
            (738, "TMTpro", 304.207146, vec![]),
        ];
        for (id, title, mass, residues) in common {
            entries.insert(
                id,
                UnimodEntry {
                    record_id: id,
                    title: title.to_string(),
                    mono_mass: mass,
                    residues,
                },
            );
        }
        Self { entries }
    }

    /// Parse the full Unimod XML database.
    pub fn from_xml(path: &Path) -> Result<Self, ResultImportError> {
        use quick_xml::events::Event;
        use quick_xml::Reader;

        if !path.exists() {
            return Err(ResultImportError::FileNotFound {
                path: path.to_path_buf(),
            });
        }

        let xml_bytes = std::fs::read(path)?;
        let mut reader = Reader::from_reader(xml_bytes.as_slice());
        reader.config_mut().trim_text(true);

        let mut entries = HashMap::new();
        let mut buf = Vec::new();

        // State for current <umod:mod> element
        let mut current_id: Option<u32> = None;
        let mut current_title: Option<String> = None;
        let mut current_mass: Option<f64> = None;
        let mut current_residues: Vec<char> = Vec::new();

        loop {
            match reader.read_event_into(&mut buf) {
                Ok(Event::Start(ref e) | Event::Empty(ref e)) => {
                    let local = e.local_name();
                    let local_str = std::str::from_utf8(local.as_ref()).unwrap_or("");

                    if local_str == "mod" {
                        // Reset state
                        current_id = None;
                        current_title = None;
                        current_mass = None;
                        current_residues.clear();

                        for attr in e.attributes().flatten() {
                            let key = std::str::from_utf8(attr.key.as_ref()).unwrap_or("");
                            let val = std::str::from_utf8(&attr.value).unwrap_or("");
                            match key {
                                "record_id" => current_id = val.parse().ok(),
                                "title" => current_title = Some(val.to_string()),
                                _ => {}
                            }
                        }
                    } else if local_str == "delta" {
                        for attr in e.attributes().flatten() {
                            let key = std::str::from_utf8(attr.key.as_ref()).unwrap_or("");
                            let val = std::str::from_utf8(&attr.value).unwrap_or("");
                            if key == "mono_mass" {
                                current_mass = val.parse().ok();
                            }
                        }
                    } else if local_str == "specificity" {
                        for attr in e.attributes().flatten() {
                            let key = std::str::from_utf8(attr.key.as_ref()).unwrap_or("");
                            let val = std::str::from_utf8(&attr.value).unwrap_or("");
                            if key == "site" && val.len() == 1 {
                                let Some(ch) = val.chars().next() else {
                                    continue;
                                };
                                if ch.is_ascii_uppercase() && !current_residues.contains(&ch) {
                                    current_residues.push(ch);
                                }
                            }
                        }
                    }
                }
                Ok(Event::End(ref e)) => {
                    let local_name = e.local_name();
                    let local = std::str::from_utf8(local_name.as_ref()).unwrap_or("");
                    if local == "mod" {
                        if let (Some(id), Some(title), Some(mass)) =
                            (current_id, current_title.take(), current_mass)
                        {
                            entries.insert(
                                id,
                                UnimodEntry {
                                    record_id: id,
                                    title,
                                    mono_mass: mass,
                                    residues: current_residues.clone(),
                                },
                            );
                        }
                    }
                }
                Ok(Event::Eof) => break,
                Err(e) => return Err(ResultImportError::XmlError(e)),
                _ => {}
            }
            buf.clear();
        }

        tracing::info!(
            "loaded {} Unimod modifications from {}",
            entries.len(),
            path.display()
        );
        Ok(Self { entries })
    }

    /// Look up a modification by Unimod record_id.
    pub fn get(&self, record_id: u32) -> Option<&UnimodEntry> {
        self.entries.get(&record_id)
    }

    /// Convert a Unimod record_id at a given position to a core::Modification.
    ///
    /// `position` is 1-based. The residue at that position in `sequence` is
    /// used to set the `residues` field of the resulting Modification.
    pub fn to_modification(
        &self,
        record_id: u32,
        position: usize,
        sequence: &str,
    ) -> Result<Modification, ResultImportError> {
        let entry = self
            .entries
            .get(&record_id)
            .ok_or(ResultImportError::UnknownUnimodId(record_id))?;

        let residue = if position >= 1 && position <= sequence.len() {
            sequence.chars().nth(position - 1).unwrap_or('X')
        } else {
            return Err(ResultImportError::InvalidModPosition {
                position,
                seq_len: sequence.len(),
            });
        };

        Ok(Modification {
            name: entry.title.clone(),
            mass_delta: entry.mono_mass,
            residues: vec![residue],
            position: ModPosition::Anywhere,
        })
    }

    /// Number of entries in the database.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the database is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builtin_has_common_mods() {
        let db = UnimodDb::builtin();
        assert!(db.len() >= 20);

        // Oxidation
        let ox = db.get(35).expect("Oxidation should be builtin");
        assert_eq!(ox.title, "Oxidation");
        assert!((ox.mono_mass - 15.994915).abs() < 0.001);

        // Carbamidomethyl
        let cam = db.get(4).expect("Carbamidomethyl should be builtin");
        assert_eq!(cam.title, "Carbamidomethyl");
        assert!((cam.mono_mass - 57.021464).abs() < 0.001);

        // Phospho
        let phos = db.get(21).expect("Phospho should be builtin");
        assert_eq!(phos.title, "Phospho");
        assert!((phos.mono_mass - 79.966331).abs() < 0.001);
    }

    #[test]
    fn to_modification_oxidation_on_m() {
        let db = UnimodDb::builtin();
        let m = db.to_modification(35, 5, "PEPTMIDE").unwrap();
        assert_eq!(m.name, "Oxidation");
        assert!((m.mass_delta - 15.994915).abs() < 0.001);
        assert_eq!(m.residues, vec!['M']);
        assert_eq!(m.position, ModPosition::Anywhere);
    }

    #[test]
    fn to_modification_unknown_id_errors() {
        let db = UnimodDb::builtin();
        assert!(db.to_modification(99999, 1, "PEPTIDE").is_err());
    }

    #[test]
    fn to_modification_invalid_position_errors() {
        let db = UnimodDb::builtin();
        assert!(db.to_modification(35, 0, "PEPTIDE").is_err());
        assert!(db.to_modification(35, 100, "PEPTIDE").is_err());
    }

    #[test]
    fn from_xml_loads_real_unimod() {
        let xml_path = Path::new("/home/verden/pfind/2025-fall/code/ms2-met/unimod.xml");
        if !xml_path.exists() {
            eprintln!("skipping XML test: unimod.xml not found");
            return;
        }
        let db = UnimodDb::from_xml(xml_path).unwrap();
        assert!(db.len() > 1000, "expected >1000 mods, got {}", db.len());

        let ox = db.get(35).expect("Oxidation");
        assert_eq!(ox.title, "Oxidation");
        assert!((ox.mono_mass - 15.994915).abs() < 0.001);
    }
}
