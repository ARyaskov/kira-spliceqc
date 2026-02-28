use std::path::Path;

use crate::expression::ExpressionMatrix;
use crate::input::error::InputError;

#[derive(Debug, Clone)]
pub struct Geneset {
    pub id: String,
    pub axis: String,
    pub gene_ids: Vec<u32>,
    pub missing: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct GenesetCatalog {
    pub genesets: Vec<Geneset>,
}

pub mod catalog;
pub mod loader;

pub fn load_catalog(
    path: &Path,
    matrix: &dyn ExpressionMatrix,
) -> Result<GenesetCatalog, InputError> {
    loader::load_catalog(path, matrix)
}
