use std::fs;
use std::io::Write;
use std::path::Path;
use std::sync::{Arc, Mutex};

use kira_spliceqc::cli::config::{AnalysisMode, RunConfig, RunMode};
use kira_spliceqc::cli::run::run_pipeline;

/// Both tests run the pipeline, but only one captures its logs through a
/// thread-local dispatcher: run concurrently, the other test's uncaptured
/// `info!` calls can cache the callsite interest as "never" before the
/// capturing dispatcher sees them, and stage lines go missing. Serialize.
static PIPELINE_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
use sha2::{Digest, Sha256};
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

fn write_tenx_extended(dir: &Path) {
    let genes = [
        "SNRPC", "SF3A1", "SF3B1", "SRSF1", "HNRNPA1", "SNRNP35", "UPF1", "POLR2A", "U2AF1",
        "PRPF3", "PRPF8",
    ];

    let mut matrix = String::new();
    matrix.push_str("%%MatrixMarket matrix coordinate integer general\n%\n");
    matrix.push_str(&format!("{} 2 {}\n", genes.len(), genes.len() * 2));

    // Cell 2 is not a scaled copy of cell 1: identical cp10k profiles would
    // give a zero MAD and undefined z-scores for every panel.
    for (idx, _gene) in genes.iter().enumerate() {
        let gene_id = idx + 1;
        matrix.push_str(&format!("{gene_id} 1 {}\n", 1 + idx % 2));
        matrix.push_str(&format!("{gene_id} 2 {}\n", 1 + (idx + 1) % 3));
    }

    fs::write(dir.join("matrix.mtx"), matrix).unwrap();

    let mut features = String::new();
    for (idx, gene) in genes.iter().enumerate() {
        features.push_str(&format!("g{}\t{}\n", idx + 1, gene));
    }
    fs::write(dir.join("features.tsv"), features).unwrap();

    let barcodes = "cell1\ncell2\n";
    fs::write(dir.join("barcodes.tsv"), barcodes).unwrap();
}

#[test]
fn extended_pipeline_logs_and_json() {
    let _guard = PIPELINE_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let input_dir = tempdir().unwrap();
    write_tenx_extended(input_dir.path());

    let out_dir = tempdir().unwrap();
    let config = RunConfig {
        input: input_dir.path().to_path_buf(),
        out_dir: out_dir.path().to_path_buf(),
        cache_path: None,
        layers: None,
        junctions: None,
        metadata: None,
        stratify_by: None,
        reference: None,
        catalog: None,
        min_counts: 0,
        min_genes: 0,
        mode: AnalysisMode::Cell,
        run_mode: RunMode::Standalone,
        output_json: true,
        output_tsv: false,
        extended: true,
        threads: None,
        experimental_signatures: true,
    };

    let buffer = Arc::new(Mutex::new(Vec::new()));
    let writer = BufferWriter {
        buffer: buffer.clone(),
    };
    let subscriber = tracing_subscriber::fmt()
        .with_writer(writer)
        .with_ansi(false)
        .with_env_filter(tracing_subscriber::EnvFilter::new("info"))
        .finish();
    let dispatch = tracing::Dispatch::new(subscriber);

    tracing::dispatcher::with_default(&dispatch, || {
        tracing::callsite::rebuild_interest_cache();
        run_pipeline(config).unwrap();
    });

    let logs = String::from_utf8(buffer.lock().unwrap().clone()).unwrap();
    for stage in [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 7] {
        assert!(
            logs.contains(&format!("starting Stage {stage}")),
            "missing stage {stage} start"
        );
    }

    let json_path = out_dir.path().join("cells.json");
    let data = fs::read_to_string(json_path).unwrap();
    let v: serde_json::Value = serde_json::from_str(&data).unwrap();
    assert!(v.get("splicing_noise").is_some());
    assert!(v.get("cryptic_risk").is_some());
    assert!(v.get("collapse").is_some());
    assert!(v.get("timecourse").is_some());
    assert!(v.get("cell_cycle").is_some());
    assert!(v["cell_cycle"]["phase_counts"].is_object());
}

#[test]
fn extended_pipeline_deterministic_hash() {
    let _guard = PIPELINE_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let input_dir = tempdir().unwrap();
    write_tenx_extended(input_dir.path());

    let out_dir1 = tempdir().unwrap();
    let out_dir2 = tempdir().unwrap();

    let config1 = RunConfig {
        input: input_dir.path().to_path_buf(),
        out_dir: out_dir1.path().to_path_buf(),
        cache_path: None,
        layers: None,
        junctions: None,
        metadata: None,
        stratify_by: None,
        reference: None,
        catalog: None,
        min_counts: 0,
        min_genes: 0,
        mode: AnalysisMode::Cell,
        run_mode: RunMode::Standalone,
        output_json: true,
        output_tsv: false,
        extended: true,
        threads: None,
        experimental_signatures: true,
    };
    let config2 = RunConfig {
        input: input_dir.path().to_path_buf(),
        out_dir: out_dir2.path().to_path_buf(),
        cache_path: None,
        layers: None,
        junctions: None,
        metadata: None,
        stratify_by: None,
        reference: None,
        catalog: None,
        min_counts: 0,
        min_genes: 0,
        mode: AnalysisMode::Cell,
        run_mode: RunMode::Standalone,
        output_json: true,
        output_tsv: false,
        extended: true,
        threads: None,
        experimental_signatures: true,
    };

    run_pipeline(config1).unwrap();
    run_pipeline(config2).unwrap();

    let bytes1 = fs::read(out_dir1.path().join("cells.json")).unwrap();
    let bytes2 = fs::read(out_dir2.path().join("cells.json")).unwrap();

    let mut hasher1 = Sha256::new();
    hasher1.update(&bytes1);
    let hash1 = hasher1.finalize();

    let mut hasher2 = Sha256::new();
    hasher2.update(&bytes2);
    let hash2 = hasher2.finalize();

    assert_eq!(hash1[..], hash2[..]);
}
