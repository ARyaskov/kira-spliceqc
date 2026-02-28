use std::fs;
use std::io::Write;
use std::path::Path;
use std::sync::{Arc, Mutex};

use kira_spliceqc::cli::config::{AnalysisMode, RunConfig, RunMode};
use kira_spliceqc::cli::run::run_pipeline;
use kira_spliceqc::expression::ExpressionMatrix;
use kira_spliceqc::expression::MmapExpressionMatrix;
use kira_spliceqc::input::detect::{detect_prefix, resolve_shared_cache_filename};
use kira_spliceqc::input::error::InputError;
use kira_spliceqc::input::shared_cache::{crc64_ecma, validate_and_open};
use tempfile::tempdir;

#[derive(Clone)]
struct BufferWriter {
    buffer: Arc<Mutex<Vec<u8>>>,
}

struct BufferGuard {
    buffer: Arc<Mutex<Vec<u8>>>,
}

impl<'a> tracing_subscriber::fmt::MakeWriter<'a> for BufferWriter {
    type Writer = BufferGuard;

    fn make_writer(&'a self) -> Self::Writer {
        BufferGuard {
            buffer: self.buffer.clone(),
        }
    }
}

impl Write for BufferGuard {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.buffer.lock().unwrap().extend_from_slice(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

fn write_tenx(dir: &Path) {
    let matrix = "%%MatrixMarket matrix coordinate integer general\n%\n5 2 5\n1 1 1\n2 1 1\n3 2 1\n4 2 1\n5 1 1\n";
    fs::write(dir.join("matrix.mtx"), matrix).unwrap();

    let features = "g1\tSNRPC\ng2\tSF3A1\ng3\tSF3B1\ng4\tSRSF1\ng5\tHNRNPA1\n";
    fs::write(dir.join("features.tsv"), features).unwrap();

    let barcodes = "cell1\ncell2\n";
    fs::write(dir.join("barcodes.tsv"), barcodes).unwrap();
}

fn write_prefixed_tenx(dir: &Path, prefix: &str) {
    let matrix = "%%MatrixMarket matrix coordinate integer general\n%\n5 2 5\n1 1 1\n2 1 1\n3 2 1\n4 2 1\n5 1 1\n";
    fs::write(dir.join(format!("{prefix}_matrix.mtx")), matrix).unwrap();

    let features = "g1\tSNRPC\ng2\tSF3A1\ng3\tSF3B1\ng4\tSRSF1\ng5\tHNRNPA1\n";
    fs::write(dir.join(format!("{prefix}_features.tsv")), features).unwrap();

    let barcodes = "cell1\ncell2\n";
    fs::write(dir.join(format!("{prefix}_barcodes.tsv")), barcodes).unwrap();
}

fn write_valid_cache(path: &Path) {
    let genes = vec!["SNRPC", "SF3A1", "SF3B1", "SRSF1", "HNRNPA1"];
    let barcodes = vec!["cell1", "cell2"];
    let col_ptr = vec![0u64, 3, 5];
    let row_idx = vec![0u32, 1, 4, 2, 3];
    let values = vec![1u32, 1, 1, 1, 1];
    write_cache_file(path, &genes, &barcodes, &col_ptr, &row_idx, &values);
}

fn write_cache_file(
    path: &Path,
    genes: &[&str],
    barcodes: &[&str],
    col_ptr: &[u64],
    row_idx: &[u32],
    values: &[u32],
) {
    assert_eq!(col_ptr.len(), barcodes.len() + 1);
    assert_eq!(row_idx.len(), values.len());
    assert_eq!(col_ptr.last().copied().unwrap(), row_idx.len() as u64);

    let genes_table = encode_string_table(genes);
    let barcodes_table = encode_string_table(barcodes);
    let nnz = row_idx.len();

    let mut file = vec![0u8; 256];

    let genes_offset = align64(file.len());
    file.resize(genes_offset, 0);
    file.extend_from_slice(&genes_table);
    let genes_bytes = genes_table.len();

    let barcodes_offset = align64(file.len());
    file.resize(barcodes_offset, 0);
    file.extend_from_slice(&barcodes_table);
    let barcodes_bytes = barcodes_table.len();

    let col_ptr_offset = align64(file.len());
    file.resize(col_ptr_offset, 0);
    for value in col_ptr {
        file.extend_from_slice(&value.to_le_bytes());
    }

    let row_idx_offset = align64(file.len());
    file.resize(row_idx_offset, 0);
    for value in row_idx {
        file.extend_from_slice(&value.to_le_bytes());
    }

    let values_offset = align64(file.len());
    file.resize(values_offset, 0);
    for value in values {
        file.extend_from_slice(&value.to_le_bytes());
    }

    let file_bytes = file.len();
    file[0..4].copy_from_slice(b"KORG");
    file[4..6].copy_from_slice(&1u16.to_le_bytes());
    file[6..8].copy_from_slice(&0u16.to_le_bytes());
    file[8..12].copy_from_slice(&0x1234_5678u32.to_le_bytes());
    file[12..16].copy_from_slice(&256u32.to_le_bytes());
    file[16..24].copy_from_slice(&(genes.len() as u64).to_le_bytes());
    file[24..32].copy_from_slice(&(barcodes.len() as u64).to_le_bytes());
    file[32..40].copy_from_slice(&(nnz as u64).to_le_bytes());
    file[40..48].copy_from_slice(&(genes_offset as u64).to_le_bytes());
    file[48..56].copy_from_slice(&(genes_bytes as u64).to_le_bytes());
    file[56..64].copy_from_slice(&(barcodes_offset as u64).to_le_bytes());
    file[64..72].copy_from_slice(&(barcodes_bytes as u64).to_le_bytes());
    file[72..80].copy_from_slice(&(col_ptr_offset as u64).to_le_bytes());
    file[80..88].copy_from_slice(&(row_idx_offset as u64).to_le_bytes());
    file[88..96].copy_from_slice(&(values_offset as u64).to_le_bytes());
    file[96..104].copy_from_slice(&0u64.to_le_bytes());
    file[104..112].copy_from_slice(&0u64.to_le_bytes());
    file[112..120].copy_from_slice(&(file_bytes as u64).to_le_bytes());
    file[120..128].copy_from_slice(&0u64.to_le_bytes());
    file[128..136].copy_from_slice(&0u64.to_le_bytes());

    let mut header = file[0..256].to_vec();
    header[120..128].fill(0);
    let crc = crc64_ecma(&header);
    file[120..128].copy_from_slice(&crc.to_le_bytes());

    fs::write(path, file).unwrap();
}

fn encode_string_table(items: &[&str]) -> Vec<u8> {
    let mut offsets = Vec::with_capacity(items.len() + 1);
    offsets.push(0u32);
    let mut blob = Vec::<u8>::new();
    for item in items {
        blob.extend_from_slice(item.as_bytes());
        offsets.push(blob.len() as u32);
    }

    let mut out = Vec::new();
    out.extend_from_slice(&(items.len() as u32).to_le_bytes());
    for offset in offsets {
        out.extend_from_slice(&offset.to_le_bytes());
    }
    out.extend_from_slice(&blob);
    out
}

fn align64(value: usize) -> usize {
    let rem = value % 64;
    if rem == 0 { value } else { value + (64 - rem) }
}

fn read_u32_le(bytes: &[u8], offset: usize) -> u32 {
    let raw: [u8; 4] = bytes[offset..offset + 4].try_into().unwrap();
    u32::from_le_bytes(raw)
}

#[test]
fn prefix_detection_prefixed_and_non_prefixed() {
    let plain = tempdir().unwrap();
    write_tenx(plain.path());
    assert_eq!(detect_prefix(plain.path()).unwrap(), None);

    let prefixed = tempdir().unwrap();
    write_prefixed_tenx(prefixed.path(), "XYZ");
    assert_eq!(
        detect_prefix(prefixed.path()).unwrap(),
        Some("XYZ".to_string())
    );
}

#[test]
fn shared_cache_filename_resolution() {
    assert_eq!(resolve_shared_cache_filename(None), "kira-organelle.bin");
    assert_eq!(
        resolve_shared_cache_filename(Some("XYZ")),
        "XYZ.kira-organelle.bin"
    );
}

#[test]
fn shared_cache_read_mmap_validation_and_traversal() {
    let temp = tempdir().unwrap();
    let cache_path = temp.path().join("kira-organelle.bin");
    write_valid_cache(&cache_path);

    let validated = validate_and_open(&cache_path).unwrap();
    assert_eq!(validated.n_genes, 5);
    assert_eq!(validated.n_cells, 2);
    assert_eq!(validated.nnz, 5);
    assert_eq!(
        validated.genes,
        vec![
            "SNRPC".to_string(),
            "SF3A1".to_string(),
            "SF3B1".to_string(),
            "SRSF1".to_string(),
            "HNRNPA1".to_string()
        ]
    );
    assert_eq!(
        validated.barcodes,
        vec!["cell1".to_string(), "cell2".to_string()]
    );
    assert_eq!(validated.col_ptr, vec![0, 3, 5]);

    let row_idx = (0..validated.nnz)
        .map(|i| read_u32_le(&validated.mmap, validated.row_idx_offset + i * 4))
        .collect::<Vec<_>>();
    let values = (0..validated.nnz)
        .map(|i| read_u32_le(&validated.mmap, validated.values_u32_offset + i * 4))
        .collect::<Vec<_>>();
    assert_eq!(row_idx, vec![0, 1, 4, 2, 3]);
    assert_eq!(values, vec![1, 1, 1, 1, 1]);

    let matrix = MmapExpressionMatrix::open_shared_cache(&cache_path).unwrap();
    assert_eq!(matrix.n_genes(), 5);
    assert_eq!(matrix.n_cells(), 2);
    assert_eq!(matrix.gene_symbol(0), "SNRPC");
    assert_eq!(matrix.cell_name(1), "cell2");
    assert_eq!(matrix.count(0, 0), 1);
    assert_eq!(matrix.count(3, 1), 1);
    assert_eq!(matrix.count(2, 0), 0);
}

#[test]
fn shared_cache_crc_tamper_is_rejected() {
    let temp = tempdir().unwrap();
    let cache_path = temp.path().join("kira-organelle.bin");
    write_valid_cache(&cache_path);

    let mut bytes = fs::read(&cache_path).unwrap();
    bytes[16] ^= 0x01;
    fs::write(&cache_path, bytes).unwrap();

    let err = validate_and_open(&cache_path).unwrap_err();
    match err {
        InputError::InvalidSharedCache(msg) => assert!(msg.contains("CRC64")),
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn pipeline_mode_uses_cache_when_present() {
    let input = tempdir().unwrap();
    let cache_path = input.path().join("kira-organelle.bin");
    write_valid_cache(&cache_path);

    let out = tempdir().unwrap();
    let config = RunConfig {
        input: input.path().to_path_buf(),
        out_dir: out.path().to_path_buf(),
        cache_path: None,
        mode: AnalysisMode::Cell,
        run_mode: RunMode::Pipeline,
        output_json: true,
        output_tsv: false,
        extended: false,
        threads: None,
    };

    run_pipeline(config).unwrap();
    assert!(
        out.path()
            .join("kira-spliceqc")
            .join("spliceqc.tsv")
            .exists()
    );
    assert!(
        out.path()
            .join("kira-spliceqc")
            .join("summary.json")
            .exists()
    );
    assert!(
        out.path()
            .join("kira-spliceqc")
            .join("pipeline_step.json")
            .exists()
    );
}

#[test]
fn pipeline_mode_missing_cache_falls_back_with_warn() {
    let input = tempdir().unwrap();
    write_prefixed_tenx(input.path(), "XYZ");
    let expected_cache = input.path().join("XYZ.kira-organelle.bin");

    let out = tempdir().unwrap();
    let config = RunConfig {
        input: input.path().to_path_buf(),
        out_dir: out.path().to_path_buf(),
        cache_path: None,
        mode: AnalysisMode::Cell,
        run_mode: RunMode::Pipeline,
        output_json: true,
        output_tsv: false,
        extended: false,
        threads: None,
    };

    let buffer = Arc::new(Mutex::new(Vec::new()));
    let writer = BufferWriter {
        buffer: buffer.clone(),
    };
    let subscriber = tracing_subscriber::fmt()
        .with_writer(writer)
        .with_ansi(false)
        .with_env_filter(tracing_subscriber::EnvFilter::new("warn"))
        .finish();
    let dispatch = tracing::Dispatch::new(subscriber);
    tracing::dispatcher::with_default(&dispatch, || {
        run_pipeline(config).unwrap();
    });

    let logs = String::from_utf8(buffer.lock().unwrap().clone()).unwrap();
    assert_eq!(
        logs.matches("shared cache not found in pipeline mode")
            .count(),
        1
    );
    assert!(logs.contains(&expected_cache.display().to_string()));
}

#[test]
fn pipeline_mode_invalid_cache_is_hard_error() {
    let input = tempdir().unwrap();
    write_tenx(input.path());
    fs::write(input.path().join("kira-organelle.bin"), b"corrupt").unwrap();

    let out = tempdir().unwrap();
    let config = RunConfig {
        input: input.path().to_path_buf(),
        out_dir: out.path().to_path_buf(),
        cache_path: None,
        mode: AnalysisMode::Cell,
        run_mode: RunMode::Pipeline,
        output_json: true,
        output_tsv: false,
        extended: false,
        threads: None,
    };

    let err = run_pipeline(config).unwrap_err();
    assert!(err.to_string().contains("invalid shared cache"));
}
