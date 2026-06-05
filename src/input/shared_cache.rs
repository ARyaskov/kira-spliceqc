use std::path::Path;

use crate::input::error::InputError;

pub fn validate_dimensions(path: &Path) -> Result<(usize, usize), InputError> {
    let (n_genes, n_cells) = kira_shared_sc_cache::validate_dimensions(path)
        .map_err(|e| InputError::InvalidSharedCache(e.to_string()))?;
    Ok((n_genes as usize, n_cells as usize))
}

pub fn crc64_ecma(data: &[u8]) -> u64 {
    kira_shared_sc_cache::crc64_ecma(data)
}
