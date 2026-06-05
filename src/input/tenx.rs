use std::path::Path;

use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

use crate::input::TenXInput;
use crate::input::detect::TenXPaths;
use crate::input::error::InputError;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TenXValidation {
    pub input: TenXInput,
    pub n_genes: usize,
    pub n_cells: usize,
    pub has_multiple_samples: bool,
    pub has_metadata: bool,
}

pub fn validate(paths: TenXPaths) -> Result<TenXValidation, InputError> {
    // Metadata-only probe — stage 1 will do the full read once.
    let reader = Reader::with_options(
        &paths.root,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    );
    let metadata = reader
        .read_metadata()
        .map_err(|e| InputError::UnsupportedInput(e.message))?;

    let n_genes = metadata.n_genes;
    let n_cells = metadata.n_cells;
    let n_features = metadata.gene_symbols.len();
    let n_barcodes = metadata.barcodes.len();

    if n_genes == 0 || n_cells == 0 || n_features == 0 || n_barcodes == 0 {
        return Err(InputError::InvalidMatrixDimensions);
    }
    if n_genes != n_features {
        return Err(InputError::DimensionMismatch {
            expected: format!("matrix rows = {n_features}"),
            found: format!("{n_genes}"),
        });
    }
    if n_cells != n_barcodes {
        return Err(InputError::DimensionMismatch {
            expected: format!("matrix cols = {n_barcodes}"),
            found: format!("{n_cells}"),
        });
    }

    check_dims(n_genes, n_cells)?;

    let has_multiple_samples = detect_multiple_samples(&paths.root)?;
    let has_metadata =
        paths.root.join("metadata.tsv").exists() || paths.root.join("metadata.tsv.gz").exists();

    Ok(TenXValidation {
        input: TenXInput {
            root: paths.root,
            matrix_path: paths.matrix_path,
            features_path: paths.features_path,
            barcodes_path: paths.barcodes_path,
            compressed: paths.compressed,
        },
        n_genes,
        n_cells,
        has_multiple_samples,
        has_metadata,
    })
}

fn detect_multiple_samples(root: &Path) -> Result<bool, InputError> {
    let mut count = 0usize;
    let entries = std::fs::read_dir(root).map_err(|e| InputError::io(root, e))?;
    for entry in entries {
        let entry = entry.map_err(|e| InputError::io(root, e))?;
        let path = entry.path();
        if !path.is_dir() {
            continue;
        }
        if kira_scio::discover(&path).is_ok() {
            count += 1;
            if count > 1 {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn check_dims(n_genes: usize, n_cells: usize) -> Result<(), InputError> {
    if n_genes == 0 || n_cells == 0 {
        return Err(InputError::InvalidMatrixDimensions);
    }
    if n_genes > u32::MAX as usize || n_cells > u32::MAX as usize {
        return Err(InputError::DimensionTooLarge { n_genes, n_cells });
    }
    Ok(())
}
