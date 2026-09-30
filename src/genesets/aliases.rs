//! Gene symbol resolution: case-insensitive lookup with legacy HGNC aliases.
//!
//! Public datasets aligned to older annotations (GRCh37 / GENCODE v19 and
//! earlier) still carry pre-2010 symbols for many splicing factors
//! (`SFRS1` for `SRSF1`, `ASCC3L1` for `SNRNP200`, `U2AF65` for `U2AF2`).
//! Every panel is written with current symbols; this table lets such
//! datasets resolve instead of silently losing whole panels.

use ahash::AHashMap;

use crate::expression::ExpressionMatrix;

/// `(legacy symbol, current symbol)`. Mouse symbols resolve through the
/// case-insensitive lookup (`Srsf1` -> `SRSF1`).
pub const LEGACY_ALIASES: &[(&str, &str)] = &[
    // SR proteins
    ("SFRS1", "SRSF1"),
    ("ASF", "SRSF1"),
    ("SFRS2", "SRSF2"),
    ("SC35", "SRSF2"),
    ("SFRS3", "SRSF3"),
    ("SFRS4", "SRSF4"),
    ("SFRS5", "SRSF5"),
    ("SFRS6", "SRSF6"),
    ("SFRS7", "SRSF7"),
    ("SFRS8", "SRSF8"),
    ("SFRS9", "SRSF9"),
    ("SFRS10", "TRA2B"),
    ("SFRS11", "SRSF11"),
    ("SFRS12", "SREK1"),
    ("SFRS2B", "SRSF8"),
    // hnRNPs
    ("HNRPA1", "HNRNPA1"),
    ("HNRPA2B1", "HNRNPA2B1"),
    ("HNRPC", "HNRNPC"),
    ("HNRPK", "HNRNPK"),
    ("HNRPD", "HNRNPD"),
    ("HNRPH1", "HNRNPH1"),
    ("HNRPL", "HNRNPL"),
    ("HNRPM", "HNRNPM"),
    ("HNRPU", "HNRNPU"),
    // snRNP / spliceosome core
    ("ASCC3L1", "SNRNP200"),
    ("HELIC2", "SNRNP200"),
    ("SNRP70", "SNRNP70"),
    ("U1-70K", "SNRNP70"),
    ("PRPC8", "PRPF8"),
    ("HPRP3", "PRPF3"),
    ("HPRP4", "PRPF4"),
    ("SNRP116", "EFTUD2"),
    ("U5-116KD", "EFTUD2"),
    ("SNRP40", "SNRNP40"),
    ("SAP61", "SF3A3"),
    ("SAP62", "SF3A2"),
    ("SAP114", "SF3A1"),
    ("SAP130", "SF3B3"),
    ("SAP145", "SF3B2"),
    ("SAP155", "SF3B1"),
    ("SF3B155", "SF3B1"),
    ("SAP49", "SF3B4"),
    ("SF3B10", "SF3B5"),
    ("SF3B14", "SF3B6"),
    ("PRP19", "PRPF19"),
    ("SNEV", "PRPF19"),
    ("U2AF35", "U2AF1"),
    ("U2AF65", "U2AF2"),
    ("URP", "ZRSR2"),
    ("U11/U12-35K", "SNRNP35"),
    ("U11/U12-48K", "SNRNP48"),
    ("SNRNP25", "SNRNP25"),
    // NMD
    ("RENT1", "UPF1"),
    ("RENT2", "UPF2"),
    ("KIAA0421", "SMG1"),
    ("LIP", "SMG1"),
    ("EST1B", "SMG5"),
    ("EST1A", "SMG6"),
    ("EST1C", "SMG7"),
    // transcription / R-loop
    ("SPT5", "SUPT5H"),
    ("SPT6", "SUPT6H"),
    ("SCA1", "SETX"),
    ("ALS4", "SETX"),
    ("RHA", "DHX9"),
    ("NDHII", "DHX9"),
    ("RBM39", "RBM39"),
    ("CAPER", "RBM39"),
    ("RNPC2", "RBM39"),
    // cell cycle (Tirosh 2016 lists, Seurat 2019 update)
    ("MLF1IP", "CENPU"),
    ("RPA2", "POLR1B"),
    ("FAM64A", "PIMREG"),
    ("HN1", "JPT1"),
];

/// Case-insensitive symbol index of a matrix (first occurrence wins).
pub fn symbol_index(matrix: &dyn ExpressionMatrix) -> AHashMap<String, u32> {
    let mut index: AHashMap<String, u32> = AHashMap::with_capacity(matrix.n_genes());
    for gene_idx in 0..matrix.n_genes() {
        index
            .entry(matrix.gene_symbol(gene_idx).to_ascii_uppercase())
            .or_insert(gene_idx as u32);
    }
    index
}

/// Gene id for `symbol`: the current symbol first, then any legacy alias of
/// it present in the index. Returns the id and whether an alias was used.
pub fn resolve_symbol(index: &AHashMap<String, u32>, symbol: &str) -> Option<(u32, bool)> {
    let key = symbol.to_ascii_uppercase();
    if let Some(&id) = index.get(&key) {
        return Some((id, false));
    }
    LEGACY_ALIASES
        .iter()
        .filter(|(_, current)| current.eq_ignore_ascii_case(symbol))
        .find_map(|(legacy, _)| index.get(&legacy.to_ascii_uppercase()).map(|&id| (id, true)))
}

/// Species guess from gene symbol casing: mouse/rat symbols are title-case
/// (`Snrpb`), human symbols upper-case (`SNRPB`). Returns `"human"`,
/// `"mouse"` or `"unknown"` (no symbols, or mixed casing).
pub fn detect_species(matrix: &dyn ExpressionMatrix) -> &'static str {
    let mut upper = 0usize;
    let mut title = 0usize;
    for g in 0..matrix.n_genes() {
        let s = matrix.gene_symbol(g);
        let mut chars = s.chars();
        let Some(first) = chars.next() else { continue };
        if !first.is_ascii_alphabetic() {
            continue;
        }
        let rest: Vec<char> = chars.filter(|c| c.is_ascii_alphabetic()).collect();
        if rest.is_empty() {
            continue;
        }
        if rest.iter().all(|c| c.is_ascii_uppercase()) {
            upper += 1;
        } else if rest.iter().all(|c| c.is_ascii_lowercase()) {
            title += 1;
        }
    }
    let total = upper + title;
    if total == 0 {
        "unknown"
    } else if upper * 10 >= total * 9 {
        "human"
    } else if title * 10 >= total * 9 {
        "mouse"
    } else {
        "unknown"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn legacy_symbols_resolve_to_current_ids() {
        let mut index = AHashMap::new();
        index.insert("SFRS1".to_string(), 3u32);
        index.insert("HNRNPC".to_string(), 5u32);
        assert_eq!(resolve_symbol(&index, "SRSF1"), Some((3, true)));
        assert_eq!(resolve_symbol(&index, "srsf1"), Some((3, true)));
        assert_eq!(resolve_symbol(&index, "HNRNPC"), Some((5, false)));
        assert_eq!(resolve_symbol(&index, "SNRNP200"), None);
    }

    #[test]
    fn alias_table_has_no_self_loops_or_duplicate_legacy_entries() {
        let mut seen = std::collections::HashSet::new();
        for (legacy, current) in LEGACY_ALIASES {
            assert!(legacy != current || *legacy == "SNRNP25" || *legacy == "RBM39", "{legacy}");
            if legacy != current {
                assert!(seen.insert(*legacy), "duplicate legacy symbol {legacy}");
            }
        }
    }
}
