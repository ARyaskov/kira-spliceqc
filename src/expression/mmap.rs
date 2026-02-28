use std::fs::File;
use std::path::Path;

use memmap2::Mmap;

use crate::expression::ExpressionMatrix;
use crate::input::error::InputError;
use crate::input::shared_cache::{SharedCacheValidation, validate_and_open};

const EXPR_MAGIC: &[u8; 8] = b"KIRAEXP1";
const EXPR_VERSION: u32 = 1;
const EXPR_HEADER_SIZE: usize = 52;

#[derive(Debug)]
enum MmapBackend {
    ExprBin {
        row_offsets: Vec<u64>,
        row_nnz: Vec<u32>,
    },
    SharedCache {
        col_ptr: Vec<u64>,
        row_idx: Vec<u32>,
        values_u32: Vec<u32>,
        nnz: usize,
    },
}

#[derive(Debug)]
pub struct MmapExpressionMatrix {
    mmap: Mmap,
    n_genes: usize,
    n_cells: usize,
    gene_symbols: Vec<String>,
    cell_names: Vec<String>,
    libsizes: Vec<u64>,
    cell_nnz: Vec<u64>,
    backend: MmapBackend,
}

impl MmapExpressionMatrix {
    pub fn open(path: &Path) -> Result<Self, InputError> {
        let file = File::open(path).map_err(|e| InputError::io(path, e))?;
        let mmap = unsafe { Mmap::map(&file).map_err(|e| InputError::io(path, e))? };
        if mmap.len() < EXPR_HEADER_SIZE {
            return Err(InputError::InvalidSparseMatrix);
        }

        if &mmap[0..8] != EXPR_MAGIC {
            return Err(InputError::InvalidSparseMatrix);
        }

        let version = read_u32(&mmap, 8)?;
        if version != EXPR_VERSION {
            return Err(InputError::InvalidSparseMatrix);
        }

        let n_genes = read_u32(&mmap, 12)? as usize;
        let n_cells = read_u32(&mmap, 16)? as usize;
        let counts_offset = read_u64(&mmap, 20)?;
        let libsize_offset = read_u64(&mmap, 28)?;
        let gene_index_offset = read_u64(&mmap, 36)?;
        let cell_index_offset = read_u64(&mmap, 44)?;

        let file_len = mmap.len() as u64;
        if !(counts_offset < libsize_offset
            && libsize_offset < gene_index_offset
            && gene_index_offset < cell_index_offset
            && cell_index_offset <= file_len)
        {
            return Err(InputError::InvalidSparseMatrix);
        }

        let (row_offsets, row_nnz) =
            parse_row_index(&mmap, n_genes, counts_offset, libsize_offset)?;
        let libsizes = parse_libsizes(&mmap, n_cells, libsize_offset, gene_index_offset)?;
        let gene_symbols =
            parse_expr_strings(&mmap, gene_index_offset, cell_index_offset, n_genes)?;
        let cell_names = parse_expr_strings(&mmap, cell_index_offset, file_len, n_cells)?;

        let cell_nnz = compute_cell_nnz_expr(&mmap, n_cells, &row_offsets, &row_nnz)?;

        Ok(Self {
            mmap,
            n_genes,
            n_cells,
            gene_symbols,
            cell_names,
            libsizes,
            cell_nnz,
            backend: MmapBackend::ExprBin {
                row_offsets,
                row_nnz,
            },
        })
    }

    pub fn open_shared_cache(path: &Path) -> Result<Self, InputError> {
        let file = File::open(path).map_err(|e| InputError::io(path, e))?;
        let mmap = unsafe { Mmap::map(&file).map_err(|e| InputError::io(path, e))? };
        let SharedCacheValidation {
            n_genes,
            n_cells,
            nnz,
            genes,
            barcodes,
            col_ptr,
            row_idx,
            values_u32,
            libsizes,
        } = validate_and_open(path)?;

        let cell_nnz = col_ptr
            .windows(2)
            .map(|w| w[1].saturating_sub(w[0]))
            .collect::<Vec<_>>();

        Ok(Self {
            mmap,
            n_genes,
            n_cells,
            gene_symbols: genes,
            cell_names: barcodes,
            libsizes,
            cell_nnz,
            backend: MmapBackend::SharedCache {
                col_ptr,
                row_idx,
                values_u32,
                nnz,
            },
        })
    }
}

impl ExpressionMatrix for MmapExpressionMatrix {
    fn n_genes(&self) -> usize {
        self.n_genes
    }

    fn n_cells(&self) -> usize {
        self.n_cells
    }

    fn gene_symbol(&self, gene_id: usize) -> &str {
        &self.gene_symbols[gene_id]
    }

    fn cell_name(&self, cell_id: usize) -> &str {
        &self.cell_names[cell_id]
    }

    fn libsize(&self, cell_id: usize) -> u64 {
        self.libsizes[cell_id]
    }

    fn nnz_cell(&self, cell_id: usize) -> u64 {
        self.cell_nnz[cell_id]
    }

    fn count(&self, gene_id: usize, cell_id: usize) -> u32 {
        match &self.backend {
            MmapBackend::ExprBin {
                row_offsets,
                row_nnz,
            } => {
                let nnz = row_nnz[gene_id] as usize;
                if nnz == 0 {
                    return 0;
                }
                let start = row_offsets[gene_id] as usize;
                let mut low = 0usize;
                let mut high = nnz;
                while low < high {
                    let mid = (low + high) / 2;
                    let offset = start + mid * 8;
                    let cell = read_u32_unchecked(&self.mmap, offset);
                    if cell == cell_id as u32 {
                        return read_u32_unchecked(&self.mmap, offset + 4);
                    }
                    if cell < cell_id as u32 {
                        low = mid + 1;
                    } else {
                        high = mid;
                    }
                }
                0
            }
            MmapBackend::SharedCache {
                col_ptr,
                row_idx,
                values_u32,
                nnz,
            } => {
                let start = col_ptr[cell_id] as usize;
                let end = col_ptr[cell_id + 1] as usize;
                if end > *nnz || start >= end {
                    return 0;
                }
                let mut low = start;
                let mut high = end;
                while low < high {
                    let mid = (low + high) / 2;
                    let row = row_idx[mid];
                    if row == gene_id as u32 {
                        return values_u32[mid];
                    }
                    if row < gene_id as u32 {
                        low = mid + 1;
                    } else {
                        high = mid;
                    }
                }
                0
            }
        }
    }
}

fn compute_cell_nnz_expr(
    mmap: &Mmap,
    n_cells: usize,
    row_offsets: &[u64],
    row_nnz: &[u32],
) -> Result<Vec<u64>, InputError> {
    let mut cell_nnz = vec![0u64; n_cells];
    for (row_offset, &nnz) in row_offsets.iter().zip(row_nnz.iter()) {
        let start = *row_offset as usize;
        for i in 0..(nnz as usize) {
            let cell = read_u32(mmap, start + i * 8)? as usize;
            if cell >= n_cells {
                return Err(InputError::InvalidSparseMatrix);
            }
            cell_nnz[cell] += 1;
        }
    }
    Ok(cell_nnz)
}

fn parse_row_index(
    mmap: &Mmap,
    n_genes: usize,
    counts_offset: u64,
    libsize_offset: u64,
) -> Result<(Vec<u64>, Vec<u32>), InputError> {
    let mut row_offsets = Vec::with_capacity(n_genes);
    let mut row_nnz = Vec::with_capacity(n_genes);
    let mut cursor = counts_offset as usize;
    let end = libsize_offset as usize;

    for _ in 0..n_genes {
        if cursor + 12 > end {
            return Err(InputError::InvalidSparseMatrix);
        }
        let row_offset = read_u64(mmap, cursor)?;
        let nnz = read_u32(mmap, cursor + 8)?;
        let row_end = row_offset
            .checked_add((nnz as u64) * 8)
            .ok_or(InputError::InvalidSparseMatrix)?;
        if row_offset < counts_offset || row_end > libsize_offset {
            return Err(InputError::InvalidSparseMatrix);
        }
        row_offsets.push(row_offset);
        row_nnz.push(nnz);
        let row_bytes = (nnz as usize).saturating_mul(8);
        cursor = cursor + 12 + row_bytes;
    }

    Ok((row_offsets, row_nnz))
}

fn parse_libsizes(
    mmap: &Mmap,
    n_cells: usize,
    libsize_offset: u64,
    gene_index_offset: u64,
) -> Result<Vec<u64>, InputError> {
    let mut libsizes = Vec::with_capacity(n_cells);
    let mut cursor = libsize_offset as usize;
    let end = gene_index_offset as usize;
    for _ in 0..n_cells {
        if cursor + 8 > end {
            return Err(InputError::InvalidSparseMatrix);
        }
        libsizes.push(read_u64(mmap, cursor)?);
        cursor += 8;
    }
    Ok(libsizes)
}

fn parse_expr_strings(
    mmap: &Mmap,
    start: u64,
    end: u64,
    expected: usize,
) -> Result<Vec<String>, InputError> {
    if expected == 0 {
        return Ok(Vec::new());
    }
    let bytes = &mmap[start as usize..end as usize];
    let mut names = Vec::with_capacity(expected);
    let mut current = Vec::new();

    for &b in bytes {
        if b == 0 {
            let name = String::from_utf8(current).map_err(|_| InputError::InvalidSparseMatrix)?;
            names.push(name);
            current = Vec::new();
            if names.len() == expected {
                break;
            }
        } else {
            current.push(b);
        }
    }

    if names.len() != expected {
        return Err(InputError::InvalidSparseMatrix);
    }
    Ok(names)
}

fn read_u32(mmap: &Mmap, offset: usize) -> Result<u32, InputError> {
    if offset + 4 > mmap.len() {
        return Err(InputError::InvalidSparseMatrix);
    }
    let bytes: [u8; 4] = mmap[offset..offset + 4].try_into().unwrap();
    Ok(u32::from_le_bytes(bytes))
}

fn read_u64(mmap: &Mmap, offset: usize) -> Result<u64, InputError> {
    if offset + 8 > mmap.len() {
        return Err(InputError::InvalidSparseMatrix);
    }
    let bytes: [u8; 8] = mmap[offset..offset + 8].try_into().unwrap();
    Ok(u64::from_le_bytes(bytes))
}

fn read_u32_unchecked(mmap: &Mmap, offset: usize) -> u32 {
    let bytes: [u8; 4] = mmap[offset..offset + 4].try_into().unwrap();
    u32::from_le_bytes(bytes)
}
