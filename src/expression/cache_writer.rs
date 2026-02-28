use std::fs::File;
use std::io::{Seek, SeekFrom, Write};
use std::path::Path;

use crate::input::error::InputError;

const MAGIC: &[u8; 8] = b"KIRAEXP1";
const VERSION: u32 = 1;
const HEADER_SIZE: u64 = 52;

pub struct CacheData {
    pub n_genes: usize,
    pub n_cells: usize,
    pub gene_symbols: Vec<String>,
    pub cell_names: Vec<String>,
    pub triplets: Vec<(u32, u32, u32)>,
    pub libsizes: Vec<u64>,
}

pub fn write_expr_bin(path: &Path, data: CacheData) -> Result<(), InputError> {
    let mut file = File::create(path).map_err(|e| InputError::io(path, e))?;

    file.seek(SeekFrom::Start(HEADER_SIZE))
        .map_err(|e| InputError::io(path, e))?;

    let counts_offset = HEADER_SIZE;
    let mut triplet_idx = 0usize;
    for gene_id in 0..data.n_genes {
        let row_start = triplet_idx;
        while triplet_idx < data.triplets.len() && data.triplets[triplet_idx].0 as usize == gene_id
        {
            triplet_idx += 1;
        }
        let nnz = (triplet_idx - row_start) as u32;
        let row_offset = file
            .stream_position()
            .map_err(|e| InputError::io(path, e))?
            + 8
            + 4;
        write_u64(&mut file, row_offset, path)?;
        write_u32(&mut file, nnz, path)?;
        for (gene, cell, count) in &data.triplets[row_start..triplet_idx] {
            debug_assert_eq!(*gene as usize, gene_id);
            write_u32(&mut file, *cell, path)?;
            write_u32(&mut file, *count, path)?;
        }
    }

    let libsize_offset = file
        .stream_position()
        .map_err(|e| InputError::io(path, e))?;
    for value in &data.libsizes {
        write_u64(&mut file, *value, path)?;
    }

    let gene_index_offset = file
        .stream_position()
        .map_err(|e| InputError::io(path, e))?;
    for name in &data.gene_symbols {
        file.write_all(name.as_bytes())
            .map_err(|e| InputError::io(path, e))?;
        file.write_all(&[0]).map_err(|e| InputError::io(path, e))?;
    }

    let cell_index_offset = file
        .stream_position()
        .map_err(|e| InputError::io(path, e))?;
    for name in &data.cell_names {
        file.write_all(name.as_bytes())
            .map_err(|e| InputError::io(path, e))?;
        file.write_all(&[0]).map_err(|e| InputError::io(path, e))?;
    }

    file.seek(SeekFrom::Start(0))
        .map_err(|e| InputError::io(path, e))?;

    file.write_all(MAGIC).map_err(|e| InputError::io(path, e))?;
    write_u32(&mut file, VERSION, path)?;
    write_u32(&mut file, data.n_genes as u32, path)?;
    write_u32(&mut file, data.n_cells as u32, path)?;
    write_u64(&mut file, counts_offset, path)?;
    write_u64(&mut file, libsize_offset, path)?;
    write_u64(&mut file, gene_index_offset, path)?;
    write_u64(&mut file, cell_index_offset, path)?;

    Ok(())
}

fn write_u32<W: Write>(writer: &mut W, value: u32, path: &Path) -> Result<(), InputError> {
    writer
        .write_all(&value.to_le_bytes())
        .map_err(|e| InputError::io(path, e))
}

fn write_u64<W: Write>(writer: &mut W, value: u64, path: &Path) -> Result<(), InputError> {
    writer
        .write_all(&value.to_le_bytes())
        .map_err(|e| InputError::io(path, e))
}
