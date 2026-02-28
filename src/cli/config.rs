use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct RunConfig {
    pub input: PathBuf,
    pub out_dir: PathBuf,
    pub cache_path: Option<PathBuf>,
    pub mode: AnalysisMode,
    pub run_mode: RunMode,
    pub output_json: bool,
    pub output_tsv: bool,
    pub extended: bool,
    pub threads: Option<usize>,
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
