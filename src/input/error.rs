use std::path::PathBuf;

#[derive(Debug, thiserror::Error)]
pub enum InputError {
    #[error("unsupported input: {0}")]
    UnsupportedInput(String),
    #[error("missing file: {0}")]
    MissingFile(String),
    #[error("invalid MatrixMarket header")]
    InvalidMatrixMarket,
    #[error("expected coordinate sparse MatrixMarket, found {0}")]
    InvalidMatrixMarketType(String),
    #[error("features.tsv must have at least 2 columns")]
    InvalidFeaturesFormat,
    #[error("features.tsv contains empty gene symbol on line {0}")]
    EmptyGeneSymbol(usize),
    #[error("barcodes.tsv must have one column")]
    InvalidBarcodesFormat,
    #[error("barcodes.tsv contains empty entry on line {0}")]
    EmptyBarcode(usize),
    #[error("dimension mismatch: expected {expected}, found {found}")]
    DimensionMismatch { expected: String, found: String },
    #[error("matrix dimensions exceed u32::MAX: genes={n_genes}, cells={n_cells}")]
    DimensionTooLarge { n_genes: usize, n_cells: usize },
    #[error("invalid matrix dimensions")]
    InvalidMatrixDimensions,
    #[error("missing gene symbols")]
    MissingGeneSymbols,
    #[error("geneset catalog missing: {0}")]
    GenesetCatalogMissing(String),
    #[error("invalid geneset catalog: {0}")]
    InvalidGenesetCatalog(String),
    #[error("empty regulator gene set")]
    EmptyRegulatorSet,
    #[error("insufficient core spliceosome genesets")]
    InsufficientCoreGenesets,
    #[error("insufficient core genesets for imbalance")]
    InsufficientCoreGenesetsImbalance,
    #[error("length mismatch: {0}")]
    LengthMismatch(String),
    #[error("output serialization error: {0}")]
    OutputSerialization(String),
    #[error("empty transcription coupling geneset")]
    EmptyCouplingGeneset,
    #[error("insufficient genesets for exon/intron bias")]
    InsufficientGenesetsExonIntron,
    #[error("insufficient assembly phase genesets")]
    InsufficientAssemblyPhaseGenesets,
    #[error("missing cryptic splicing risk signals")]
    MissingCrypticRiskSignals,
    #[error("missing spliceosome collapse signals")]
    MissingCollapseSignals,
    #[error("missing timecourse signals")]
    MissingTimecourseSignals,
    #[error("insufficient splicing noise genesets")]
    InsufficientSplicingNoiseGenesets,
    #[error("unsupported h5ad layout")]
    UnsupportedH5ADLayout,
    #[error("invalid sparse matrix")]
    InvalidSparseMatrix,
    #[error("invalid shared cache: {0}")]
    InvalidSharedCache(String),
    #[error("gene index overflow")]
    GeneIndexOverflow,
    #[error("cell index overflow")]
    CellIndexOverflow,
    #[error("missing dataset: {0}")]
    MissingDataset(String),
    #[error("io error at {path}: {source}")]
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    #[error("hdf5 error: {0}")]
    Hdf5(#[from] hdf5::Error),
}

impl InputError {
    pub fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        Self::Io {
            path: path.into(),
            source,
        }
    }
}
