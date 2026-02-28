use std::fs;

use std::path::Path;

use hdf5::types::VarLenUnicode;
use kira_spliceqc::cli::config::RunMode;
use kira_spliceqc::input::error::InputError;
use kira_spliceqc::input::{GeneSymbolSource, InputKind};
use kira_spliceqc::pipeline::stage0_input::run_stage0;
use tempfile::tempdir;

fn write_tenx(dir: &Path) {
    let matrix = "%%MatrixMarket matrix coordinate integer general\n%\n2 3 2\n1 1 1\n2 3 1\n";
    fs::write(dir.join("matrix.mtx"), matrix).unwrap();

    let features = "g1\tGeneA\ng2\tGeneB\n";
    fs::write(dir.join("features.tsv"), features).unwrap();

    let barcodes = "cell1\ncell2\ncell3\n";
    fs::write(dir.join("barcodes.tsv"), barcodes).unwrap();
}

fn write_h5ad_csr(path: &Path, n_cells: i64, n_genes: i64) {
    let file = hdf5::File::create(path).unwrap();

    let x_group = file.create_group("X").unwrap();
    let attr = x_group
        .new_attr::<VarLenUnicode>()
        .create("encoding-type")
        .unwrap();
    attr.write_scalar(&unsafe { VarLenUnicode::from_str_unchecked("csr_matrix") })
        .unwrap();

    x_group
        .new_dataset::<i32>()
        .shape(2)
        .create("data")
        .unwrap()
        .write(&[1, 2])
        .unwrap();
    x_group
        .new_dataset::<i32>()
        .shape(2)
        .create("indices")
        .unwrap()
        .write(&[0, 1])
        .unwrap();
    x_group
        .new_dataset::<i32>()
        .shape(3)
        .create("indptr")
        .unwrap()
        .write(&[0, 1, 2])
        .unwrap();
    x_group
        .new_dataset::<i64>()
        .shape(2)
        .create("shape")
        .unwrap()
        .write(&[n_cells, n_genes])
        .unwrap();

    let var_group = file.create_group("var").unwrap();
    var_group
        .new_dataset::<i32>()
        .shape(n_genes as usize)
        .create("gene_symbols")
        .unwrap()
        .write(&vec![1; n_genes as usize])
        .unwrap();
    var_group
        .new_dataset::<i32>()
        .shape(n_genes as usize)
        .create("_index")
        .unwrap()
        .write(&vec![1; n_genes as usize])
        .unwrap();

    let obs_group = file.create_group("obs").unwrap();
    obs_group
        .new_dataset::<i32>()
        .shape(n_cells as usize)
        .create("_index")
        .unwrap()
        .write(&vec![1; n_cells as usize])
        .unwrap();
    obs_group
        .new_dataset::<i32>()
        .shape(n_cells as usize)
        .create("qc")
        .unwrap()
        .write(&vec![1; n_cells as usize])
        .unwrap();
}

fn write_h5ad_dense(path: &Path, n_cells: usize, n_genes: usize) {
    let file = hdf5::File::create(path).unwrap();
    file.new_dataset::<i32>()
        .shape(n_cells * n_genes)
        .create("X")
        .unwrap()
        .write(&vec![0; n_cells * n_genes])
        .unwrap();
    let var_group = file.create_group("var").unwrap();
    var_group
        .new_dataset::<i32>()
        .shape(n_genes)
        .create("gene_symbols")
        .unwrap()
        .write(&vec![1; n_genes])
        .unwrap();
    let obs_group = file.create_group("obs").unwrap();
    obs_group
        .new_dataset::<i32>()
        .shape(n_cells)
        .create("_index")
        .unwrap()
        .write(&vec![1; n_cells])
        .unwrap();
}

#[test]
fn valid_tenx_directory() {
    let temp = tempdir().unwrap();
    write_tenx(temp.path());

    let descriptor = run_stage0(temp.path(), RunMode::Standalone).unwrap();
    assert_eq!(descriptor.n_genes, 2);
    assert_eq!(descriptor.n_cells, 3);
    assert!(!descriptor.has_metadata);
    assert!(!descriptor.has_multiple_samples);

    match descriptor.kind {
        InputKind::TenX(tenx) => {
            assert_eq!(tenx.root, temp.path());
            assert!(!tenx.compressed);
        }
        _ => panic!("expected TenX input"),
    }
}

#[test]
fn missing_file_errors() {
    let temp = tempdir().unwrap();
    write_tenx(temp.path());
    fs::remove_file(temp.path().join("barcodes.tsv")).unwrap();

    let err = run_stage0(temp.path(), RunMode::Standalone).unwrap_err();
    match err {
        InputError::MissingFile(name) => assert_eq!(name, "barcodes.tsv"),
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn valid_h5ad_csr() {
    let temp = tempdir().unwrap();
    let path = temp.path().join("sample.h5ad");
    write_h5ad_csr(&path, 3, 2);

    let descriptor = run_stage0(&path, RunMode::Standalone, None).unwrap();
    assert_eq!(descriptor.n_genes, 2);
    assert_eq!(descriptor.n_cells, 3);

    match descriptor.kind {
        InputKind::H5AD(h5ad) => {
            assert!(h5ad.x_is_csr);
            assert_eq!(h5ad.gene_symbol_source, GeneSymbolSource::VarGeneSymbols);
        }
        _ => panic!("expected H5AD input"),
    }
}

#[test]
fn h5ad_dense_rejected() {
    let temp = tempdir().unwrap();
    let path = temp.path().join("dense.h5ad");
    write_h5ad_dense(&path, 2, 2);

    let err = run_stage0(&path, RunMode::Standalone, None).unwrap_err();
    match err {
        InputError::UnsupportedH5ADLayout => {}
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn deterministic_detection() {
    let temp = tempdir().unwrap();
    write_tenx(temp.path());

    let first = run_stage0(temp.path(), RunMode::Standalone).unwrap();
    let second = run_stage0(temp.path(), RunMode::Standalone).unwrap();

    assert_eq!(first, second);
}
