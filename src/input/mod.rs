use std::path::PathBuf;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InputKind {
    TenX(TenXInput),
    H5AD(H5ADInput),
    OrganelleCache(OrganelleCacheInput),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InputDescriptor {
    pub kind: InputKind,
    pub n_genes: usize,
    pub n_cells: usize,
    pub has_multiple_samples: bool,
    pub has_metadata: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TenXInput {
    pub root: PathBuf,
    pub matrix_path: PathBuf,
    pub features_path: PathBuf,
    pub barcodes_path: PathBuf,
    pub compressed: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct H5ADInput {
    pub path: PathBuf,
    pub x_is_csr: bool,
    pub gene_symbol_source: GeneSymbolSource,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OrganelleCacheInput {
    pub root: PathBuf,
    pub cache_path: PathBuf,
    pub prefix: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeneSymbolSource {
    VarGeneSymbols,
    VarNamesFallback,
}

pub mod detect;
pub mod error;
pub mod h5ad;
pub mod shared_cache;
pub mod tenx;
