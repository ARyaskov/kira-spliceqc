//! Spliced/unspliced layer loading (input level L1): detection next to a 10x
//! directory, the STARsolo `Gene/` <-> `Velocyto/` sibling layout with barcode
//! matching, AnnData `layers/` groups, and the `--layers` override.

use std::fs;
use std::path::Path;

use hdf5::types::VarLenUnicode;
use kira_spliceqc::cli::config::RunMode;
use kira_spliceqc::expression::ExpressionMatrix;
use kira_spliceqc::io::layers::LayerLocation;
use kira_spliceqc::pipeline::stage0_input::run_stage0_with_layers;
use kira_spliceqc::pipeline::stage1_expression::run_stage1_full;
use tempfile::tempdir;

/// 3 genes x 3 cells; cell/gene names deliberately unsorted so reindexing is
/// exercised (stage 1 sorts genes and cells lexicographically).
fn write_main(dir: &Path) {
    fs::create_dir_all(dir).unwrap();
    fs::write(
        dir.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n3 3 4\n1 1 10\n2 1 20\n3 2 30\n1 3 40\n",
    )
    .unwrap();
    fs::write(dir.join("features.tsv"), "g1\tGeneB\ng2\tGeneA\ng3\tGeneC\n").unwrap();
    fs::write(dir.join("barcodes.tsv"), "cellB\ncellA\ncellC\n").unwrap();
}

fn write_layer(path: &Path, body: &str) {
    fs::write(
        path,
        format!("%%MatrixMarket matrix coordinate integer general\n{body}"),
    )
    .unwrap();
}

#[test]
fn layers_next_to_matrix_are_detected_and_reindexed() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("data");
    write_main(&input);
    // Raw index space: gene 1 = GeneB, cell 1 = cellB.
    write_layer(&input.join("spliced.mtx"), "3 3 2\n1 1 7\n2 3 2\n");
    write_layer(&input.join("unspliced.mtx"), "3 3 2\n1 1 3\n3 2 5\n");

    let stage0 = run_stage0_with_layers(&input, RunMode::Standalone, None, None).unwrap();
    assert_eq!(stage0.layers, Some(LayerLocation::MtxDir(input.clone())));

    let out = tempdir().unwrap();
    let stage1 = run_stage1_full(&stage0, out.path()).unwrap();
    let layers = stage1.layers.expect("layers loaded");
    let m = &stage1.matrix;

    // After sorting: cells = [cellA, cellB, cellC], genes = [GeneA, GeneB, GeneC].
    let cell = |name: &str| (0..m.n_cells()).find(|&c| m.cell_name(c) == name).unwrap();
    let gene = |name: &str| (0..m.n_genes()).find(|&g| m.gene_symbol(g) == name).unwrap();

    assert_eq!(layers.spliced.count(gene("GeneB"), cell("cellB")), 7);
    assert_eq!(layers.unspliced.count(gene("GeneB"), cell("cellB")), 3);
    assert_eq!(layers.spliced.count(gene("GeneA"), cell("cellC")), 2);
    assert_eq!(layers.unspliced.count(gene("GeneC"), cell("cellA")), 5);
    assert_eq!(layers.spliced.cell_total(cell("cellB")), 7);
    assert_eq!(layers.unspliced.cell_total(cell("cellA")), 5);
    assert!(layers.ambiguous.is_none());
    assert_eq!(layers.cells_without_layers, 0);
    assert!(layers.source.starts_with("mtx-dir:"));
}

#[test]
fn starsolo_sibling_layout_matches_cells_by_barcode() {
    let dir = tempdir().unwrap();
    let solo = dir.path().join("Solo.out");
    let gene_dir = solo.join("Gene").join("filtered");
    write_main(&gene_dir);

    // Velocyto/filtered carries only two of the three cells, in another order.
    let velo = solo.join("Velocyto").join("filtered");
    fs::create_dir_all(&velo).unwrap();
    fs::write(velo.join("barcodes.tsv"), "cellC\ncellB\n").unwrap();
    fs::write(velo.join("features.tsv"), "g1\tGeneB\ng2\tGeneA\ng3\tGeneC\n").unwrap();
    // layer column 1 = cellC, column 2 = cellB
    write_layer(&velo.join("spliced.mtx"), "3 2 2\n1 1 4\n1 2 6\n");
    write_layer(&velo.join("unspliced.mtx"), "3 2 1\n2 2 9\n");
    write_layer(&velo.join("ambiguous.mtx"), "3 2 1\n3 1 1\n");

    let stage0 = run_stage0_with_layers(&gene_dir, RunMode::Standalone, None, None).unwrap();
    assert_eq!(stage0.layers, Some(LayerLocation::MtxDir(velo.clone())));

    let out = tempdir().unwrap();
    let stage1 = run_stage1_full(&stage0, out.path()).unwrap();
    let layers = stage1.layers.unwrap();
    let m = &stage1.matrix;
    let cell = |name: &str| (0..m.n_cells()).find(|&c| m.cell_name(c) == name).unwrap();
    let gene = |name: &str| (0..m.n_genes()).find(|&g| m.gene_symbol(g) == name).unwrap();

    assert_eq!(layers.spliced.count(gene("GeneB"), cell("cellC")), 4);
    assert_eq!(layers.spliced.count(gene("GeneB"), cell("cellB")), 6);
    assert_eq!(layers.unspliced.count(gene("GeneA"), cell("cellB")), 9);
    assert_eq!(layers.ambiguous.as_ref().unwrap().count(gene("GeneC"), cell("cellC")), 1);
    // cellA has no column in the layer files.
    assert_eq!(layers.spliced.cell_total(cell("cellA")), 0);
    assert_eq!(layers.cells_without_layers, 1);
}

#[test]
fn layers_override_wins_and_dimension_mismatch_is_an_error() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("data");
    write_main(&input);
    let elsewhere = dir.path().join("layers");
    fs::create_dir_all(&elsewhere).unwrap();
    write_layer(&elsewhere.join("spliced.mtx"), "3 3 1\n1 1 1\n");
    write_layer(&elsewhere.join("unspliced.mtx"), "3 3 1\n1 1 1\n");

    let stage0 =
        run_stage0_with_layers(&input, RunMode::Standalone, None, Some(&elsewhere)).unwrap();
    assert_eq!(stage0.layers, Some(LayerLocation::MtxDir(elsewhere.clone())));
    let out = tempdir().unwrap();
    assert!(run_stage1_full(&stage0, out.path()).unwrap().layers.is_some());

    // Wrong gene count in the layer -> LayerMismatch.
    write_layer(&elsewhere.join("unspliced.mtx"), "2 3 1\n1 1 1\n");
    let out = tempdir().unwrap();
    let err = run_stage1_full(&stage0, out.path()).unwrap_err().to_string();
    assert!(err.contains("layers do not match"), "{err}");
}

#[test]
fn no_layers_means_level_zero_only() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("data");
    write_main(&input);
    let stage0 = run_stage0_with_layers(&input, RunMode::Standalone, None, None).unwrap();
    assert!(stage0.layers.is_none());
    let out = tempdir().unwrap();
    assert!(run_stage1_full(&stage0, out.path()).unwrap().layers.is_none());
}

fn write_sparse_group(file: &hdf5::File, path: &str, encoding: &str, indptr: &[i32], indices: &[i32], data: &[f32], shape: [u64; 2]) {
    let g = file.create_group(path).unwrap();
    g.new_attr::<VarLenUnicode>()
        .create("encoding-type")
        .unwrap()
        .write_scalar(&unsafe { VarLenUnicode::from_str_unchecked(encoding) })
        .unwrap();
    g.new_attr::<u64>().shape(2).create("shape").unwrap().write(&shape).unwrap();
    g.new_dataset::<i32>().shape(indptr.len()).create("indptr").unwrap().write(indptr).unwrap();
    g.new_dataset::<i32>().shape(indices.len()).create("indices").unwrap().write(indices).unwrap();
    g.new_dataset::<f32>().shape(data.len()).create("data").unwrap().write(data).unwrap();
}

fn write_strings(group: &hdf5::Group, name: &str, values: &[&str]) {
    let v: Vec<VarLenUnicode> = values
        .iter()
        .map(|s| unsafe { VarLenUnicode::from_str_unchecked(*s) })
        .collect();
    group
        .new_dataset::<VarLenUnicode>()
        .shape(values.len())
        .create(name)
        .unwrap()
        .write(&v)
        .unwrap();
}

#[test]
fn h5ad_layers_are_read_in_x_order() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("data.h5ad");
    let file = hdf5::File::create(&path).unwrap();
    // X: 2 cells x 3 genes, CSR. cell0: gene0=5, gene2=1; cell1: gene1=2.
    write_sparse_group(&file, "X", "csr_matrix", &[0, 2, 3], &[0, 2, 1], &[5.0, 1.0, 2.0], [2, 3]);
    // spliced (CSR): cell0 gene0=4; cell1 gene1=2.
    write_sparse_group(&file, "layers/spliced", "csr_matrix", &[0, 1, 2], &[0, 1], &[4.0, 2.0], [2, 3]);
    // unspliced (CSC to exercise the other encoding): gene0: cell0=1; gene2: cell0=1.
    write_sparse_group(&file, "layers/unspliced", "csc_matrix", &[0, 1, 1, 2], &[0, 0], &[1.0, 1.0], [2, 3]);
    let var = file.create_group("var").unwrap();
    write_strings(&var, "_index", &["GeneZ", "GeneY", "GeneX"]);
    let obs = file.create_group("obs").unwrap();
    write_strings(&obs, "_index", &["cell2", "cell1"]);
    drop(file);

    let stage0 = run_stage0_with_layers(&path, RunMode::Standalone, None, None).unwrap();
    assert_eq!(stage0.layers, Some(LayerLocation::H5ad(path.clone())));
    let out = tempdir().unwrap();
    let stage1 = run_stage1_full(&stage0, out.path()).unwrap();
    let layers = stage1.layers.unwrap();
    let m = &stage1.matrix;
    let cell = |name: &str| (0..m.n_cells()).find(|&c| m.cell_name(c) == name).unwrap();
    let gene = |name: &str| (0..m.n_genes()).find(|&g| m.gene_symbol(g) == name).unwrap();

    assert_eq!(m.count(gene("GeneZ"), cell("cell2")), 5);
    assert_eq!(layers.spliced.count(gene("GeneZ"), cell("cell2")), 4);
    assert_eq!(layers.spliced.count(gene("GeneY"), cell("cell1")), 2);
    assert_eq!(layers.unspliced.count(gene("GeneZ"), cell("cell2")), 1);
    assert_eq!(layers.unspliced.count(gene("GeneX"), cell("cell2")), 1);
    assert_eq!(layers.unspliced.cell_total(cell("cell1")), 0);
    assert!(layers.source.starts_with("h5ad-layers:"));
}
