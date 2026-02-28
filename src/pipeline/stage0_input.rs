use std::path::{Path, PathBuf};

use tracing::{error, info, warn};

use crate::cli::config::RunMode;
use crate::input::detect::{DetectedInput, detect_input};
use crate::input::error::InputError;
use crate::input::shared_cache::validate_dimensions;
use crate::input::{InputDescriptor, InputKind, OrganelleCacheInput};
use crate::input::{h5ad, tenx};

pub fn run_stage0(
    path: &Path,
    run_mode: RunMode,
    cache_override: Option<&Path>,
) -> Result<InputDescriptor, InputError> {
    if run_mode == RunMode::Pipeline {
        if let Some(cache_path) = cache_override {
            return descriptor_from_cache(path, cache_path.to_path_buf(), None);
        }
        if let Some(descriptor) = try_pipeline_cache(path)? {
            return Ok(descriptor);
        }
    }

    let result = match detect_input(path)? {
        DetectedInput::TenX(paths) => {
            info!("detected input format: 10x");
            let validation = tenx::validate(paths)?;
            info!(
                n_genes = validation.n_genes,
                n_cells = validation.n_cells,
                "input dimensions"
            );
            Ok(InputDescriptor {
                kind: InputKind::TenX(validation.input),
                n_genes: validation.n_genes,
                n_cells: validation.n_cells,
                has_multiple_samples: validation.has_multiple_samples,
                has_metadata: validation.has_metadata,
            })
        }
        DetectedInput::H5AD(path) => {
            info!("detected input format: h5ad");
            let validation = h5ad::validate(&path)?;
            info!(
                n_genes = validation.n_genes,
                n_cells = validation.n_cells,
                "input dimensions"
            );
            Ok(InputDescriptor {
                kind: InputKind::H5AD(validation.input),
                n_genes: validation.n_genes,
                n_cells: validation.n_cells,
                has_multiple_samples: validation.has_multiple_samples,
                has_metadata: validation.has_metadata,
            })
        }
    };

    if let Err(ref err) = result {
        error!(error = ?err, "stage0 input validation failed");
    }

    result
}

fn try_pipeline_cache(path: &Path) -> Result<Option<InputDescriptor>, InputError> {
    let meta = std::fs::metadata(path).map_err(|e| InputError::io(path, e))?;
    if !meta.is_dir() {
        return Ok(None);
    }

    let prefix = crate::input::detect::detect_prefix(path)?;
    let cache_name = crate::input::detect::resolve_shared_cache_filename(prefix.as_deref());
    let cache_path = path.join(cache_name);

    if !cache_path.exists() {
        warn!(
            expected_cache_path = %cache_path.display(),
            "shared cache not found in pipeline mode, falling back to MatrixMarket input"
        );
        return Ok(None);
    }

    Ok(Some(descriptor_from_cache(path, cache_path, prefix)?))
}

fn descriptor_from_cache(
    root: &Path,
    cache_path: PathBuf,
    prefix: Option<String>,
) -> Result<InputDescriptor, InputError> {
    let (n_genes, n_cells) = validate_dimensions(&cache_path)?;
    info!("detected input format: shared-cache");
    info!(n_genes = n_genes, n_cells = n_cells, "input dimensions");

    Ok(InputDescriptor {
        kind: InputKind::OrganelleCache(OrganelleCacheInput {
            root: root.to_path_buf(),
            cache_path,
            prefix,
        }),
        n_genes,
        n_cells,
        has_multiple_samples: false,
        has_metadata: root.join("metadata.tsv").exists() || root.join("metadata.tsv.gz").exists(),
    })
}
