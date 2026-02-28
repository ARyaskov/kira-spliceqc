use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

use crate::input::TenXInput;
use crate::input::error::InputError;
use crate::io::RawMatrix;

pub fn read_tenx(input: &TenXInput) -> Result<RawMatrix, InputError> {
    let reader = Reader::with_options(
        &input.root,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    );
    let canonical = reader
        .read_all()
        .map_err(|e| InputError::UnsupportedInput(e.message))?;

    let mut triplets = Vec::with_capacity(canonical.matrix.values.len());
    for (col, window) in canonical.matrix.col_ptr.windows(2).enumerate() {
        for idx in window[0]..window[1] {
            let row = canonical.matrix.row_idx[idx];
            let value = canonical.matrix.values[idx];
            if value < 0.0 || value.fract().abs() > 1e-6 {
                return Err(InputError::InvalidSparseMatrix);
            }
            triplets.push((row as u32, col as u32, value as u32));
        }
    }

    Ok(RawMatrix {
        genes: canonical.metadata.gene_symbols,
        cells: canonical.metadata.barcodes,
        triplets,
    })
}
