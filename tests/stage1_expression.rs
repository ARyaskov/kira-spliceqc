use std::fs;
use std::path::Path;

use hdf5::types::VarLenUnicode;
use kira_spliceqc::cli::config::RunMode;
use kira_spliceqc::expression::ExpressionMatrix;
use kira_spliceqc::pipeline::stage0_input::run_stage0;
use kira_spliceqc::pipeline::stage1_expression::run_stage1;
use tempfile::tempdir;

fn write_tenx(dir: &Path) {
    let matrix =
        "%%MatrixMarket matrix coordinate integer general\n%\n3 3 3\n1 2 1\n2 1 2\n3 3 3\n";
    fs::write(dir.join("matrix.mtx"), matrix).unwrap();

    let features = "g1\tGeneB\ng2\tGeneA\ng3\tGeneC\n";
    fs::write(dir.join("features.tsv"), features).unwrap();

    let barcodes = "cell2\ncell1\ncell3\n";
    fs::write(dir.join("barcodes.tsv"), barcodes).unwrap();
}

fn write_h5ad_csr(path: &Path) {
    let file = hdf5::File::create(path).unwrap();
    let x_group = file.create_group("X").unwrap();
    let attr = x_group
        .new_attr::<VarLenUnicode>()
        .create("encoding-type")
        .unwrap();
    attr.write_scalar(&unsafe { VarLenUnicode::from_str_unchecked("csr_matrix") })
        .unwrap();

    x_group
        .new_dataset::<i64>()
        .shape(2)
        .create("data")
        .unwrap()
        .write(&[1_i64, 2])
        .unwrap();
    x_group
        .new_dataset::<i64>()
        .shape(2)
        .create("indices")
        .unwrap()
        .write(&[0_i64, 1])
        .unwrap();
    x_group
        .new_dataset::<i64>()
        .shape(3)
        .create("indptr")
        .unwrap()
        .write(&[0_i64, 1, 2])
        .unwrap();
    x_group
        .new_dataset::<i64>()
        .shape(2)
        .create("shape")
        .unwrap()
        .write(&[2_i64, 2])
        .unwrap();

    let var_group = file.create_group("var").unwrap();
    let gene_symbols = vec![
        unsafe { VarLenUnicode::from_str_unchecked("GeneB") },
        unsafe { VarLenUnicode::from_str_unchecked("GeneA") },
    ];
    var_group
        .new_dataset::<VarLenUnicode>()
        .shape(2)
        .create("gene_symbols")
        .unwrap()
        .write(&gene_symbols)
        .unwrap();
    var_group
        .new_dataset::<VarLenUnicode>()
        .shape(2)
        .create("_index")
        .unwrap()
        .write(&gene_symbols)
        .unwrap();

    let obs_group = file.create_group("obs").unwrap();
    let cells = vec![
        unsafe { VarLenUnicode::from_str_unchecked("cell2") },
        unsafe { VarLenUnicode::from_str_unchecked("cell1") },
    ];
    obs_group
        .new_dataset::<VarLenUnicode>()
        .shape(2)
        .create("_index")
        .unwrap()
        .write(&cells)
        .unwrap();
}

#[test]
fn tiny_mtx_counts_and_norm() {
    let temp = tempdir().unwrap();
    write_tenx(temp.path());

    let descriptor = run_stage0(temp.path(), RunMode::Standalone).unwrap();
    let matrix = run_stage1(&descriptor, temp.path()).unwrap();

    assert_eq!(matrix.n_genes(), 3);
    assert_eq!(matrix.n_cells(), 3);

    assert_eq!(matrix.gene_symbol(0), "GeneA");
    assert_eq!(matrix.gene_symbol(1), "GeneB");
    assert_eq!(matrix.gene_symbol(2), "GeneC");

    assert_eq!(matrix.cell_name(0), "cell1");
    assert_eq!(matrix.cell_name(1), "cell2");
    assert_eq!(matrix.cell_name(2), "cell3");

    assert_eq!(matrix.count(0, 1), 2);
    assert_eq!(matrix.count(1, 0), 1);
    assert_eq!(matrix.count(2, 2), 3);
    assert_eq!(matrix.count(0, 0), 0);

    assert_eq!(matrix.libsize(0), 1);
    assert_eq!(matrix.libsize(1), 2);
    assert_eq!(matrix.libsize(2), 3);

    let cp10k = matrix.cp10k(0, 1);
    let log_cp10k = matrix.log_cp10k(0, 1);
    assert!((cp10k - 10000.0).abs() < 1e-6);
    assert!((log_cp10k - (10001.0_f32).ln()).abs() < 1e-6);
}

#[test]
fn deterministic_expr_bin() {
    let temp = tempdir().unwrap();
    write_tenx(temp.path());

    let descriptor = run_stage0(temp.path(), RunMode::Standalone).unwrap();

    let out1 = temp.path().join("out1");
    let out2 = temp.path().join("out2");
    fs::create_dir_all(&out1).unwrap();
    fs::create_dir_all(&out2).unwrap();

    run_stage1(&descriptor, &out1).unwrap();
    run_stage1(&descriptor, &out2).unwrap();

    let bytes1 = fs::read(out1.join("expr.bin")).unwrap();
    let bytes2 = fs::read(out2.join("expr.bin")).unwrap();
    assert_eq!(bytes1, bytes2);
}

#[test]
fn h5ad_smoke() {
    let temp = tempdir().unwrap();
    let path = temp.path().join("sample.h5ad");
    write_h5ad_csr(&path);

    let descriptor = run_stage0(&path, RunMode::Standalone, None).unwrap();
    let matrix = run_stage1(&descriptor, temp.path()).unwrap();

    assert_eq!(matrix.n_genes(), 2);
    assert_eq!(matrix.n_cells(), 2);
    assert_eq!(matrix.count(0, 0), 2);
}
