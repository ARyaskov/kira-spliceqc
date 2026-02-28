use std::path::{Path, PathBuf};

use crate::input::error::InputError;

pub const SHARED_CACHE_BASENAME: &str = "kira-organelle.bin";

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DetectedInput {
    TenX(TenXPaths),
    H5AD(PathBuf),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TenXPaths {
    pub root: PathBuf,
    pub matrix_path: PathBuf,
    pub features_path: PathBuf,
    pub barcodes_path: PathBuf,
    pub compressed: bool,
}

pub fn detect_input(path: &Path) -> Result<DetectedInput, InputError> {
    let meta = std::fs::metadata(path).map_err(|e| InputError::io(path, e))?;

    if meta.is_dir() {
        let ds = kira_scio::discover(path).map_err(|e| InputError::UnsupportedInput(e.message))?;
        let features_path = ds
            .features
            .or(ds.genes)
            .ok_or_else(|| InputError::MissingFile("features.tsv".to_string()))?;
        let barcodes_path = ds
            .barcodes
            .ok_or_else(|| InputError::MissingFile("barcodes.tsv".to_string()))?;

        let compressed = ds.matrix.extension().is_some_and(|ext| ext == "gz")
            || features_path.extension().is_some_and(|ext| ext == "gz")
            || barcodes_path.extension().is_some_and(|ext| ext == "gz");

        return Ok(DetectedInput::TenX(TenXPaths {
            root: path.to_path_buf(),
            matrix_path: ds.matrix,
            features_path,
            barcodes_path,
            compressed,
        }));
    }

    if meta.is_file() {
        if path.extension().is_some_and(|ext| ext == "h5ad") {
            return Ok(DetectedInput::H5AD(path.to_path_buf()));
        }
        return Err(InputError::UnsupportedInput(path.display().to_string()));
    }

    Err(InputError::UnsupportedInput(path.display().to_string()))
}

pub fn detect_prefix(dir: &Path) -> Result<Option<String>, InputError> {
    kira_scio::detect_prefix(dir).map_err(|e| InputError::UnsupportedInput(e.to_string()))
}

pub fn resolve_shared_cache_filename(prefix: Option<&str>) -> String {
    kira_shared_sc_cache::resolve_shared_cache_filename(prefix)
}
