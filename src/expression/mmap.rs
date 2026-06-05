use std::fs::File;
use std::path::Path;

use kira_shared_sc_cache::{SharedCacheMmap, mmap_shared_cache};
use memmap2::Mmap;

use crate::expression::ExpressionMatrix;
use crate::input::error::InputError;

const EXPR_MAGIC: &[u8; 8] = b"KIRAEXP1";
const EXPR_VERSION: u32 = 1;
const EXPR_HEADER_SIZE: usize = 52;

/// CSR layout — one row per gene, `(cell, count)` pairs sorted by cell.
struct ExprBinBackend {
    mmap: Mmap,
    row_offsets: Vec<u64>,
    row_nnz: Vec<u32>,
}

/// CSC layout — borrows the validated `kira-organelle.bin` shared cache,
/// zero-copy slices for `col_ptr`/`row_idx`/`values_u32`.
struct SharedCacheBackend {
    cache: SharedCacheMmap,
}

enum Backend {
    ExprBin(ExprBinBackend),
    SharedCache(SharedCacheBackend),
}

pub struct MmapExpressionMatrix {
    n_genes: usize,
    n_cells: usize,
    gene_symbols: Vec<String>,
    cell_names: Vec<String>,
    libsizes: Vec<u64>,
    /// Pre-computed for ExprBin; for SharedCache we derive O(1) from col_ptr.
    cell_nnz: Option<Vec<u64>>,
    backend: Backend,
}

impl std::fmt::Debug for MmapExpressionMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MmapExpressionMatrix")
            .field("n_genes", &self.n_genes)
            .field("n_cells", &self.n_cells)
            .field(
                "backend",
                &match self.backend {
                    Backend::ExprBin(_) => "ExprBin",
                    Backend::SharedCache(_) => "SharedCache",
                },
            )
            .finish()
    }
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
            n_genes,
            n_cells,
            gene_symbols,
            cell_names,
            libsizes,
            cell_nnz: Some(cell_nnz),
            backend: Backend::ExprBin(ExprBinBackend {
                mmap,
                row_offsets,
                row_nnz,
            }),
        })
    }

    pub fn open_shared_cache(path: &Path) -> Result<Self, InputError> {
        let cache = mmap_shared_cache(path)
            .map_err(|e| InputError::InvalidSharedCache(e.to_string()))?;

        let n_genes = cache.n_genes;
        let n_cells = cache.n_cells;
        let gene_symbols = cache.genes.clone();
        let cell_names = cache.barcodes.clone();

        // Compute libsizes once from zero-copy values slice.
        let col_ptr = cache.col_ptr();
        let values = cache.values_u32();
        let mut libsizes = vec![0u64; n_cells];
        for cell in 0..n_cells {
            let s = col_ptr[cell] as usize;
            let e = col_ptr[cell + 1] as usize;
            let mut sum = 0u64;
            for &v in &values[s..e] {
                sum += v as u64;
            }
            libsizes[cell] = sum;
        }

        Ok(Self {
            n_genes,
            n_cells,
            gene_symbols,
            cell_names,
            libsizes,
            cell_nnz: None, // computed O(1) via col_ptr deltas.
            backend: Backend::SharedCache(SharedCacheBackend { cache }),
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
        match (&self.backend, &self.cell_nnz) {
            (_, Some(v)) => v[cell_id],
            (Backend::SharedCache(b), None) => {
                let cp = b.cache.col_ptr();
                cp[cell_id + 1] - cp[cell_id]
            }
            // Unreachable: SharedCache always has col_ptr; ExprBin always pre-computes cell_nnz.
            _ => 0,
        }
    }

    fn count(&self, gene_id: usize, cell_id: usize) -> u32 {
        match &self.backend {
            Backend::ExprBin(b) => exprbin_count(b, gene_id, cell_id),
            Backend::SharedCache(b) => sharedcache_count(b, gene_id, cell_id),
        }
    }

    fn gather_panel_log_cp10k(
        &self,
        panel_sorted: &[u32],
        cell: usize,
        out: &mut Vec<f32>,
    ) {
        debug_assert!(is_sorted_unique(panel_sorted));
        out.reserve(panel_sorted.len());
        let scale = 1e4_f32 / self.libsizes[cell].max(1) as f32;

        match &self.backend {
            Backend::SharedCache(b) => {
                let cp = b.cache.col_ptr();
                let rows = b.cache.row_idx();
                let vals = b.cache.values_u32();
                let s = cp[cell] as usize;
                let e = cp[cell + 1] as usize;
                let col_rows = &rows[s..e];
                let col_vals = &vals[s..e];
                // Inline sparse merge; missing genes contribute log(1+0)=0.
                let mut pi = 0usize;
                let mut ri = 0usize;
                while pi < panel_sorted.len() {
                    if ri >= col_rows.len() {
                        out.push(0.0);
                        pi += 1;
                        continue;
                    }
                    let pg = panel_sorted[pi];
                    let rg = col_rows[ri];
                    if pg == rg {
                        let x = col_vals[ri] as f32;
                        out.push((1.0 + scale * x).ln());
                        pi += 1;
                        ri += 1;
                    } else if pg < rg {
                        out.push(0.0);
                        pi += 1;
                    } else {
                        ri += 1;
                    }
                }
            }
            Backend::ExprBin(b) => {
                for &g in panel_sorted {
                    let c = exprbin_count(b, g as usize, cell) as f32;
                    out.push((1.0 + scale * c).ln());
                }
            }
        }
    }

    fn panel_ln1p_scaled_sum(&self, panel_sorted: &[u32], cell: usize, scale: f32) -> f32 {
        debug_assert!(is_sorted_unique(panel_sorted));
        let mut sum = 0.0_f32;
        match &self.backend {
            Backend::SharedCache(b) => {
                let cp = b.cache.col_ptr();
                let rows = b.cache.row_idx();
                let vals = b.cache.values_u32();
                let s = cp[cell] as usize;
                let e = cp[cell + 1] as usize;
                // Only matched genes contribute; zeros add ln(1)=0.
                let col_rows = &rows[s..e];
                let col_vals = &vals[s..e];
                let mut pi = 0usize;
                let mut ri = 0usize;
                while pi < panel_sorted.len() && ri < col_rows.len() {
                    let pg = panel_sorted[pi];
                    let rg = col_rows[ri];
                    if pg == rg {
                        let v = col_vals[ri] as f32;
                        sum += (1.0 + scale * v).ln();
                        pi += 1;
                        ri += 1;
                    } else if pg < rg {
                        pi += 1;
                    } else {
                        ri += 1;
                    }
                }
            }
            Backend::ExprBin(b) => {
                for &g in panel_sorted {
                    let c = exprbin_count(b, g as usize, cell) as f32;
                    if c > 0.0 {
                        sum += (1.0 + scale * c).ln();
                    }
                }
            }
        }
        sum
    }

    fn panel_count_sum_and_detected(
        &self,
        panel_sorted: &[u32],
        cell: usize,
    ) -> (u64, usize) {
        debug_assert!(is_sorted_unique(panel_sorted));
        let mut sum = 0u64;
        let mut detected = 0usize;
        match &self.backend {
            Backend::SharedCache(b) => {
                let cp = b.cache.col_ptr();
                let rows = b.cache.row_idx();
                let vals = b.cache.values_u32();
                let s = cp[cell] as usize;
                let e = cp[cell + 1] as usize;
                let col_rows = &rows[s..e];
                let col_vals = &vals[s..e];
                let mut pi = 0usize;
                let mut ri = 0usize;
                while pi < panel_sorted.len() && ri < col_rows.len() {
                    let pg = panel_sorted[pi];
                    let rg = col_rows[ri];
                    if pg == rg {
                        let v = col_vals[ri];
                        if v > 0 {
                            detected += 1;
                            sum += v as u64;
                        }
                        pi += 1;
                        ri += 1;
                    } else if pg < rg {
                        pi += 1;
                    } else {
                        ri += 1;
                    }
                }
            }
            Backend::ExprBin(b) => {
                for &g in panel_sorted {
                    let c = exprbin_count(b, g as usize, cell);
                    if c > 0 {
                        detected += 1;
                        sum += c as u64;
                    }
                }
            }
        }
        (sum, detected)
    }
}

#[inline]
fn is_sorted_unique(panel: &[u32]) -> bool {
    panel.windows(2).all(|w| w[0] < w[1])
}

fn exprbin_count(b: &ExprBinBackend, gene_id: usize, cell_id: usize) -> u32 {
    let nnz = b.row_nnz[gene_id] as usize;
    if nnz == 0 {
        return 0;
    }
    let start = b.row_offsets[gene_id] as usize;
    let target = cell_id as u32;
    let mut low = 0usize;
    let mut high = nnz;
    while low < high {
        let mid = (low + high) / 2;
        let offset = start + mid * 8;
        let cell = read_u32_unchecked(&b.mmap, offset);
        if cell == target {
            return read_u32_unchecked(&b.mmap, offset + 4);
        }
        if cell < target {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    0
}

fn sharedcache_count(b: &SharedCacheBackend, gene_id: usize, cell_id: usize) -> u32 {
    let cp = b.cache.col_ptr();
    let start = cp[cell_id] as usize;
    let end = cp[cell_id + 1] as usize;
    if start >= end {
        return 0;
    }
    let rows = &b.cache.row_idx()[start..end];
    let vals = &b.cache.values_u32()[start..end];
    let target = gene_id as u32;
    match rows.binary_search(&target) {
        Ok(idx) => vals[idx],
        Err(_) => 0,
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
