use std::path::{Path, PathBuf};

pub fn default_catalog_path() -> PathBuf {
    let relative = Path::new("resources")
        .join("genesets")
        .join("splicing_genesets.tsv");
    if relative.is_file() {
        return relative;
    }

    let manifest = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("resources")
        .join("genesets")
        .join("splicing_genesets.tsv");
    if manifest.is_file() {
        return manifest;
    }

    if let Ok(exe) = std::env::current_exe()
        && let Some(dir) = exe.parent()
    {
        let sibling = dir
            .join("resources")
            .join("genesets")
            .join("splicing_genesets.tsv");
        if sibling.is_file() {
            return sibling;
        }
        let parent = dir
            .join("..")
            .join("resources")
            .join("genesets")
            .join("splicing_genesets.tsv");
        if parent.is_file() {
            return parent;
        }
    }

    relative
}
