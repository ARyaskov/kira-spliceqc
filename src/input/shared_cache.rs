use std::path::Path;

use crate::input::error::InputError;

#[derive(Debug)]
pub struct SharedCacheValidation {
    pub n_genes: usize,
    pub n_cells: usize,
    pub nnz: usize,
    pub genes: Vec<String>,
    pub barcodes: Vec<String>,
    pub col_ptr: Vec<u64>,
    pub row_idx: Vec<u32>,
    pub values_u32: Vec<u32>,
    pub libsizes: Vec<u64>,
}

pub fn validate_and_open(path: &Path) -> Result<SharedCacheValidation, InputError> {
    let cache = kira_shared_sc_cache::read_shared_cache_owned(path)
        .map_err(|e| InputError::InvalidSharedCache(e.to_string()))?;

    let mut libsizes = vec![0u64; cache.n_cells as usize];
    for cell in 0..cache.n_cells as usize {
        let start = cache.col_ptr[cell] as usize;
        let end = cache.col_ptr[cell + 1] as usize;
        let mut sum = 0u64;
        for v in &cache.values_u32[start..end] {
            sum += *v as u64;
        }
        libsizes[cell] = sum;
    }

    Ok(SharedCacheValidation {
        n_genes: cache.n_genes as usize,
        n_cells: cache.n_cells as usize,
        nnz: cache.nnz as usize,
        genes: cache.genes,
        barcodes: cache.barcodes,
        col_ptr: cache.col_ptr,
        row_idx: cache.row_idx,
        values_u32: cache.values_u32,
        libsizes,
    })
}

pub fn validate_dimensions(path: &Path) -> Result<(usize, usize), InputError> {
    let (n_genes, n_cells) = kira_shared_sc_cache::validate_dimensions(path)
        .map_err(|e| InputError::InvalidSharedCache(e.to_string()))?;
    Ok((n_genes as usize, n_cells as usize))
}

pub fn crc64_ecma(data: &[u8]) -> u64 {
    kira_shared_sc_cache::crc64_ecma(data)
}
