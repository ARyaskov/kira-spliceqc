use std::path::{Path, PathBuf};

use tracing::{error, info, warn};

use crate::cli::config::RunMode;
use crate::input::detect::{DetectedInput, detect_input};
use crate::input::error::InputError;
use crate::input::shared_cache::validate_dimensions;
use crate::input::{InputDescriptor, InputKind, OrganelleCacheInput};
use crate::input::{h5ad, tenx};
use crate::io::junctions::{JunctionLocation, detect_junctions};
use crate::io::layers::{LayerLocation, detect_mtx_layers, h5ad_has_layers};

pub fn run_stage0(
    path: &Path,
    run_mode: RunMode,
    cache_override: Option<&Path>,
) -> Result<InputDescriptor, InputError> {
    run_stage0_with_layers(path, run_mode, cache_override, None)
}

/// Stage 0 with an explicit spliced/unspliced layer source (`--layers`).
/// Without an override, layers are auto-detected next to the input (see
/// `io::layers::detect_mtx_layers`) or inside the AnnData file.
pub fn run_stage0_with_layers(
    path: &Path,
    run_mode: RunMode,
    cache_override: Option<&Path>,
    layers_override: Option<&Path>,
) -> Result<InputDescriptor, InputError> {
    run_stage0_full(path, run_mode, cache_override, layers_override, None)
}

/// Stage 0 with explicit layer (`--layers`) and junction (`--junctions`)
/// sources; both are auto-detected when absent.
pub fn run_stage0_full(
    path: &Path,
    run_mode: RunMode,
    cache_override: Option<&Path>,
    layers_override: Option<&Path>,
    junctions_override: Option<&Path>,
) -> Result<InputDescriptor, InputError> {
    let mut descriptor = run_stage0_inner(path, run_mode, cache_override)?;
    descriptor.layers = detect_layers(&descriptor, layers_override)?;
    match &descriptor.layers {
        Some(location) => {
            info!(layers = %location.describe(), "spliced/unspliced layers detected (input level L1)")
        }
        None => info!("no spliced/unspliced layers found (input level L0 only)"),
    }
    if let Some(p) = junctions_override
        && !p.is_dir()
    {
        return Err(InputError::MissingFile(p.display().to_string()));
    }
    descriptor.junctions = match &descriptor.kind {
        InputKind::TenX(tenx) => detect_junctions(&tenx.root, junctions_override),
        InputKind::OrganelleCache(cache) => detect_junctions(&cache.root, junctions_override),
        InputKind::H5AD(_) => junctions_override.map(|p| JunctionLocation::MtxDir(p.to_path_buf())),
    };
    match &descriptor.junctions {
        Some(location) => {
            info!(junctions = %location.describe(), "junction matrix detected (input level L2)")
        }
        None => info!("no junction matrix found (Tier B unavailable)"),
    }
    Ok(descriptor)
}

fn detect_layers(
    descriptor: &InputDescriptor,
    layers_override: Option<&Path>,
) -> Result<Option<LayerLocation>, InputError> {
    if let Some(p) = layers_override {
        if !p.exists() {
            return Err(InputError::MissingFile(p.display().to_string()));
        }
        return Ok(Some(if p.is_dir() {
            LayerLocation::MtxDir(p.to_path_buf())
        } else {
            LayerLocation::H5ad(p.to_path_buf())
        }));
    }
    Ok(match &descriptor.kind {
        InputKind::TenX(tenx) => detect_mtx_layers(&tenx.root, None),
        InputKind::H5AD(h5) => {
            h5ad_has_layers(&h5.path).then(|| LayerLocation::H5ad(h5.path.clone()))
        }
        InputKind::OrganelleCache(cache) => detect_mtx_layers(&cache.root, None),
    })
}

fn run_stage0_inner(
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
                layers: None,
                junctions: None,
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
                layers: None,
                junctions: None,
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
        layers: None,
        junctions: None,
    })
}
