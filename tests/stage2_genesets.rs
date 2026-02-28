use std::fs;
use std::path::Path;

use kira_spliceqc::expression::ExpressionMatrix;
use kira_spliceqc::genesets::loader::load_catalog;
use kira_spliceqc::pipeline::stage2_genesets::aggregate;
use tempfile::tempdir;

struct MockMatrix {
    genes: Vec<String>,
    cells: Vec<String>,
    counts: Vec<Vec<u32>>,
    libsizes: Vec<u64>,
}

impl MockMatrix {
    fn new(genes: Vec<&str>, cells: Vec<&str>, counts: Vec<Vec<u32>>) -> Self {
        let mut libsizes = vec![0u64; cells.len()];
        for g in 0..genes.len() {
            for c in 0..cells.len() {
                libsizes[c] += counts[g][c] as u64;
            }
        }
        Self {
            genes: genes.into_iter().map(|s| s.to_string()).collect(),
            cells: cells.into_iter().map(|s| s.to_string()).collect(),
            counts,
            libsizes,
        }
    }
}

impl ExpressionMatrix for MockMatrix {
    fn n_genes(&self) -> usize {
        self.genes.len()
    }

    fn n_cells(&self) -> usize {
        self.cells.len()
    }

    fn gene_symbol(&self, gene_id: usize) -> &str {
        &self.genes[gene_id]
    }

    fn cell_name(&self, cell_id: usize) -> &str {
        &self.cells[cell_id]
    }

    fn libsize(&self, cell_id: usize) -> u64 {
        self.libsizes[cell_id]
    }

    fn count(&self, gene_id: usize, cell_id: usize) -> u32 {
        self.counts[gene_id][cell_id]
    }
}

fn write_geneset_tsv(path: &Path, contents: &str) {
    fs::write(path, contents).unwrap();
}

#[test]
fn tiny_aggregation_mean() {
    let matrix = MockMatrix::new(
        vec!["A", "B", "C"],
        vec!["c1", "c2"],
        vec![vec![10, 0], vec![0, 0], vec![10, 10]],
    );

    let temp = tempdir().unwrap();
    let tsv = temp.path().join("genesets.tsv");
    write_geneset_tsv(
        &tsv,
        "geneset_id\taxis\tgene_symbol\nGS1\tAX1\tA\nGS1\tAX1\tC\n",
    );

    let catalog = load_catalog(&tsv, &matrix).unwrap();
    let activity = aggregate(&matrix, &catalog).unwrap();

    assert_eq!(activity.genesets, vec!["GS1".to_string()]);
    assert_eq!(activity.axes, vec!["AX1".to_string()]);
    assert_eq!(activity.n_cells, 2);

    let v0 = activity.value(0, 0);
    let v1 = activity.value(0, 1);

    let expected0 = (matrix.log_cp10k(0, 0) + matrix.log_cp10k(2, 0)) / 2.0;
    let expected1 = (matrix.log_cp10k(0, 1) + matrix.log_cp10k(2, 1)) / 2.0;

    assert!((v0 - expected0).abs() < 1e-6);
    assert!((v1 - expected1).abs() < 1e-6);
}

#[test]
fn missing_gene_nan() {
    let matrix = MockMatrix::new(vec!["A"], vec!["c1"], vec![vec![1]]);
    let temp = tempdir().unwrap();
    let tsv = temp.path().join("genesets.tsv");
    write_geneset_tsv(&tsv, "geneset_id\taxis\tgene_symbol\nGS1\tAX1\tZ\n");

    let catalog = load_catalog(&tsv, &matrix).unwrap();
    let activity = aggregate(&matrix, &catalog).unwrap();

    assert!(activity.value(0, 0).is_nan());
}

#[test]
fn deterministic_values() {
    let matrix = MockMatrix::new(
        vec!["A", "B"],
        vec!["c1", "c2"],
        vec![vec![1, 2], vec![3, 4]],
    );
    let temp = tempdir().unwrap();
    let tsv = temp.path().join("genesets.tsv");
    write_geneset_tsv(
        &tsv,
        "geneset_id\taxis\tgene_symbol\nGS1\tAX1\tB\nGS1\tAX1\tA\n",
    );

    let catalog = load_catalog(&tsv, &matrix).unwrap();
    let first = aggregate(&matrix, &catalog).unwrap();
    let second = aggregate(&matrix, &catalog).unwrap();

    assert_eq!(first.values, second.values);
}
