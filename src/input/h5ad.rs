use std::path::Path;

use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

use crate::input::error::InputError;
use crate::input::{GeneSymbolSource, H5ADInput};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct H5ADValidation {
    pub input: H5ADInput,
    pub n_genes: usize,
    pub n_cells: usize,
    pub has_multiple_samples: bool,
    pub has_metadata: bool,
}

pub fn validate(path: &Path) -> Result<H5ADValidation, InputError> {
    // Cheap metadata-only probe; stage 1 will do the full read once.
    let reader = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::H5ad),
            strict: true,
        },
    );
    let metadata = reader
        .read_metadata()
        .map_err(|e| InputError::UnsupportedInput(e.message))?;

    let n_genes = metadata.n_genes;
    let n_cells = metadata.n_cells;
    check_dims(n_genes, n_cells)?;

    Ok(H5ADValidation {
        input: H5ADInput {
            path: path.to_path_buf(),
            x_is_csr: true,
            gene_symbol_source: GeneSymbolSource::VarNamesFallback,
        },
        n_genes,
        n_cells,
        has_multiple_samples: false,
        has_metadata: true,
    })
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
