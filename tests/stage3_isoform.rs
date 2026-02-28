use kira_spliceqc::expression::ExpressionMatrix;
use kira_spliceqc::genesets::{Geneset, GenesetCatalog};
use kira_spliceqc::pipeline::stage3_isoform::compute;

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

fn catalog_with_ids(ids: Vec<u32>) -> GenesetCatalog {
    GenesetCatalog {
        genesets: vec![Geneset {
            id: "SRSF_SR".to_string(),
            axis: "REG".to_string(),
            gene_ids: ids,
            missing: Vec::new(),
        }],
    }
}

#[test]
fn uniform_expression_high_entropy() {
    let matrix = MockMatrix::new(vec!["G1", "G2"], vec!["c1"], vec![vec![10], vec![10]]);
    let catalog = catalog_with_ids(vec![0, 1]);
    let metrics = compute(&matrix, &catalog).unwrap();

    assert!((metrics.entropy[0] - 1.0).abs() < 1e-5);
    assert!((metrics.dispersion[0] - 1.0).abs() < 1e-5);
}

#[test]
fn skewed_expression_low_entropy() {
    let matrix = MockMatrix::new(vec!["G1", "G2"], vec!["c1"], vec![vec![100], vec![1]]);
    let catalog = catalog_with_ids(vec![0, 1]);
    let metrics = compute(&matrix, &catalog).unwrap();

    assert!(metrics.entropy[0] < 0.5);
    assert!(metrics.dispersion[0] < 0.8);
}

#[test]
fn deterministic_outputs() {
    let matrix = MockMatrix::new(
        vec!["G1", "G2"],
        vec!["c1", "c2"],
        vec![vec![1, 2], vec![3, 4]],
    );
    let catalog = catalog_with_ids(vec![0, 1]);
    let first = compute(&matrix, &catalog).unwrap();
    let second = compute(&matrix, &catalog).unwrap();

    assert_eq!(first.entropy, second.entropy);
    assert_eq!(first.dispersion, second.dispersion);
    assert_eq!(first.z_entropy, second.z_entropy);
}
