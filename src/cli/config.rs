use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct RunConfig {
    pub input: PathBuf,
    pub out_dir: PathBuf,
    pub cache_path: Option<PathBuf>,
    /// Explicit spliced/unspliced layer source (directory with
    /// `spliced.mtx`/`unspliced.mtx`, or an AnnData file with `layers/`).
    /// Auto-detected when absent.
    pub layers: Option<PathBuf>,
    pub mode: AnalysisMode,
    pub run_mode: RunMode,
    pub output_json: bool,
    pub output_tsv: bool,
    pub extended: bool,
    pub threads: Option<usize>,
    /// Write experimental composite signatures (SIS/class, SOS/RLR/SII and
    /// their flags, cryptic risk, collapse) to the per-cell outputs. Pipeline
    /// mode implies this because the pipeline contract is built on them.
    pub experimental_signatures: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnalysisMode {
    Cell,
    Sample,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunMode {
    Standalone,
    Pipeline,
}
