use std::path::Path;

use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

use crate::input::error::InputError;

#[derive(Debug, Clone)]
pub struct RawMatrix {
    pub genes: Vec<String>,
    pub cells: Vec<String>,
    pub triplets: Vec<(u32, u32, u32)>,
}

pub mod h5ad;
pub mod mtx;

/// Shared CSC → triplets loader for stage 1 ingestion.
///
/// Validates that values are non-negative integers (with a permissive tolerance
/// for rounding artifacts from upstream normalization tools) and that no count
/// overflows u32.
pub(crate) fn read_via_scio(path: &Path, format: DetectedFormat) -> Result<RawMatrix, InputError> {
    const FRAC_TOL: f32 = 1e-4;
    const U32_MAX_F32: f32 = u32::MAX as f32;

    let reader = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(format),
            strict: true,
        },
    );
    let canonical = reader
        .read_all()
        .map_err(|e| InputError::UnsupportedInput(e.message))?;

    let nnz = canonical.matrix.values.len();
    let mut triplets = Vec::with_capacity(nnz);
    let row_idx = canonical.matrix.row_idx.as_slice();
    let values = canonical.matrix.values.as_slice();
    for (col, window) in canonical.matrix.col_ptr.windows(2).enumerate() {
        let s = window[0] as usize;
        let e = window[1] as usize;
        for idx in s..e {
            let row = row_idx[idx];
            let value = values[idx];
            if !value.is_finite() || value < 0.0 || value > U32_MAX_F32 {
                return Err(InputError::InvalidSparseMatrix);
            }
            if value.fract().abs() > FRAC_TOL {
                return Err(InputError::InvalidSparseMatrix);
            }
            triplets.push((row, col as u32, value as u32));
        }
    }

    Ok(RawMatrix {
        genes: canonical.metadata.gene_symbols,
        cells: canonical.metadata.barcodes,
        triplets,
    })
}
