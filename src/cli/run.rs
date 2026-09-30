use std::time::Instant;

use tracing::{info, warn};

use crate::cli::PipelineContext;
use crate::cli::config::{AnalysisMode, RunConfig, RunMode};
use crate::expression::ExpressionMatrix;
use crate::genesets::catalog::default_catalog_path;
use crate::genesets::controls::ControlPool;
use crate::genesets::load_catalog;
use crate::metrics::splicing_instability::panel_gene_ids as instability_panel_gene_ids;
use crate::input::error::InputError;
use crate::output::pipeline_contract;
use crate::pipeline::stage0_input::run_stage0_with_layers;
use crate::pipeline::stage1_expression::{Stage1Output, run_stage1_full};
use crate::pipeline::stage2_genesets::{aggregate_with_controls as run_stage2, standardize_activity};
use crate::pipeline::stage3_isoform::compute as run_stage3;
use crate::pipeline::stage4_missplicing::compute as compute_missplicing;
use crate::pipeline::stage5_imbalance::compute as compute_imbalance;
use crate::pipeline::stage6_sis::run_stage6;
use crate::pipeline::stage7_output::{OutputOptions, run_stage7};
use crate::pipeline::stage8_coupling::compute as compute_coupling;
use crate::pipeline::stage9_exon_intron_bias::compute as compute_exon_intron;
use crate::pipeline::stage10_assembly_phase::compute as compute_assembly;
use crate::pipeline::stage11_splicing_noise::compute as compute_splicing_noise;
use crate::pipeline::stage12_cryptic_risk::compute as compute_cryptic_risk;
use crate::pipeline::stage13_collapse::compute as compute_collapse;
use crate::pipeline::stage15_splicing_instability::compute as compute_splicing_instability;
use crate::pipeline::stage16_unspliced::run_stage16;
use crate::pipeline::stage17_intron_retention::run_stage17;
use crate::pipeline::stage18_cell_cycle::run_stage18;
use crate::metrics::cell_cycle::cell_cycle_gene_ids;
use crate::reference::Strata;

/// Scratch directory (inside the output directory) for the stage-1 cache.
const EXPR_CACHE_DIR: &str = ".kira-spliceqc-cache";

#[derive(Debug, thiserror::Error)]
pub enum SpliceQcError {
    #[error("invalid input: {0}")]
    InvalidInput(String),
    #[error("pipeline failure: {0}")]
    PipelineFailure(String),
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

pub fn run_pipeline(config: RunConfig) -> Result<(), SpliceQcError> {
    if config.mode != AnalysisMode::Cell {
        return Err(SpliceQcError::InvalidInput(
            "sample mode is not implemented".to_string(),
        ));
    }

    if let Some(threads) = config.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .map_err(|e| SpliceQcError::InvalidInput(format!("failed to set threads: {e}")))?;
    }

    info!(
        target: "kira_spliceqc::cli::run",
        simd_backend = crate::simd::backend(),
        "compute backend selected"
    );

    // The pipeline contract (regimes, fidelity, confidence) is built entirely
    // on the composite signatures, so pipeline mode implies them.
    let experimental = config.experimental_signatures || config.run_mode == RunMode::Pipeline;
    if experimental {
        info!(
            target: "kira_spliceqc::cli::run",
            implied_by_pipeline_mode = !config.experimental_signatures,
            "experimental composite signatures enabled (unvalidated; see METRICS.md)"
        );
    }

    let effective_out_dir = if config.run_mode == RunMode::Pipeline {
        pipeline_contract::pipeline_out_dir(&config.out_dir)
    } else {
        config.out_dir.clone()
    };
    std::fs::create_dir_all(&effective_out_dir)?;

    let stage0 = run_logged(0, || {
        run_stage0_with_layers(
            &config.input,
            config.run_mode,
            config.cache_path.as_deref(),
            config.layers.as_deref(),
        )
    })?;
    // Stage 1 materialises an internal expression cache (expr.bin). It is an
    // implementation detail, so it lives in a hidden scratch directory that is
    // removed once the run completes instead of polluting the output directory.
    let cache_dir = effective_out_dir.join(EXPR_CACHE_DIR);
    let Stage1Output {
        matrix: stage1,
        layers: stage1_layers,
        metadata: stage1_metadata,
    } = run_logged(1, || run_stage1_full(&stage0, &cache_dir, config.metadata.as_deref()))?;
    let strata = Strata::from_metadata(&stage1_metadata, stage1.n_cells(), config.stratify_by.as_deref());
    info!(
        target: "kira_spliceqc::cli::run",
        reference = strata.mode.as_str(),
        column = strata.column.as_deref().unwrap_or("-"),
        strata = strata.n_strata(),
        metadata_columns = stage1_metadata.columns.len(),
        "reference strata"
    );
    if let Some(layers) = &stage1_layers {
        info!(
            target: "kira_spliceqc::cli::run",
            source = layers.source.as_str(),
            "input level L1 available: spliced/unspliced layers loaded"
        );
    }

    // Load catalog ONCE and pass &GenesetCatalog to all downstream consumers
    // (previously stages 2, 3 and pipeline contract each reloaded it).
    let catalog_path = default_catalog_path();
    let catalog = load_catalog(&catalog_path, &stage1)?;

    // Control-gene pool for depth correction: every gene of any splicing panel
    // is excluded from the pool (see genesets::controls).
    let controls = {
        let mut exclude: Vec<u32> = catalog.genesets.iter().flat_map(|g| g.gene_ids.iter().copied()).collect();
        exclude.extend(instability_panel_gene_ids(&stage1));
        exclude.extend(cell_cycle_gene_ids(&stage1));
        let pool = ControlPool::new(stage1.gene_mean_log_cp10k(), &exclude);
        info!(
            target: "kira_spliceqc::cli::run",
            eligible_control_genes = pool.n_eligible(),
            controls_per_gene = pool.controls_per_gene,
            "control-gene pool built"
        );
        pool
    };

    let libsizes: Vec<u64> = (0..stage1.n_cells()).map(|c| stage1.libsize(c)).collect();
    let stage2_raw = run_logged(2, || run_stage2(&stage1, &catalog, Some(&controls)))?;
    let stage2 = standardize_activity(&stage2_raw, &strata, Some(&libsizes));
    let stage3 = run_logged(3, || run_stage3(&stage1, &catalog, &strata))?;
    let stage4 = run_logged(4, || compute_missplicing(&stage2))?;
    let stage5 = run_logged(5, || compute_imbalance(&stage2))?;
    let stage6 = run_logged(6, || run_stage6(&stage3, &stage4, &stage5))?;
    let stage15 = run_logged(15, || {
        Ok(compute_splicing_instability(&stage1, Some(&controls), &strata))
    })?;
    let stage18 = run_logged(18, || Ok(run_stage18(&stage1, Some(&controls))))?;
    let (stage16, stage17) = match &stage1_layers {
        Some(layers) => (
            Some(run_logged(16, || Ok(run_stage16(layers, &strata)))?),
            Some(run_logged(17, || Ok(run_stage17(layers, &strata)))?),
        ),
        None => {
            info!(target: "kira_spliceqc::cli::run", "skipping Stages 16-17 (no spliced/unspliced layers)");
            (None, None)
        }
    };

    let (stage8, stage9, stage10, stage11, stage12, stage13, stage14) = if config.extended {
        let stage8 = run_logged(8, || compute_coupling(&stage2))?;
        let stage9 = run_logged(9, || compute_exon_intron(&stage2))?;
        let stage10 = run_logged(10, || compute_assembly(&stage2))?;
        let stage11 = run_logged(11, || compute_splicing_noise(&stage2_raw))?;
        let stage12 = run_logged(12, || compute_cryptic_risk(&stage3, &stage5))?;
        let stage13 = run_logged(13, || compute_collapse(&stage5, &stage6))?;
        (
            Some(stage8),
            Some(stage9),
            Some(stage10),
            Some(stage11),
            Some(stage12),
            Some(stage13),
            None,
        )
    } else {
        (None, None, None, None, None, None, None)
    };

    let cell_names = (0..stage1.n_cells())
        .map(|i| stage1.cell_name(i).to_string())
        .collect::<Vec<_>>();

    let context = PipelineContext {
        stage0,
        stage1,
        stage1_layers,
        stage2,
        stage3,
        stage4,
        stage5,
        stage6,
        stage8,
        stage9,
        stage10,
        stage11,
        stage12,
        stage13,
        stage14,
        stage15,
        stage16,
        stage17,
        stage18,
    };

    if config.extended && context.stage14.is_none() {
        info!(target: "kira_spliceqc::cli::run", "skipping Stage 14 (no timecourse metadata)");
    }

    let summary = run_logged(7, || {
        run_stage7(
            &effective_out_dir,
            &cell_names,
            &context.stage3,
            &context.stage4,
            &context.stage5,
            &context.stage6,
            context.stage8.as_ref(),
            context.stage9.as_ref(),
            context.stage10.as_ref(),
            context.stage11.as_ref(),
            context.stage12.as_ref(),
            context.stage13.as_ref(),
            context.stage14.as_ref(),
            &context.stage15,
            context.stage16.as_ref(),
            context.stage17.as_ref(),
            &context.stage18,
            &strata,
            OutputOptions {
                json: config.output_json,
                tsv: config.output_tsv,
                experimental,
            },
        )
    })?;

    if config.run_mode == RunMode::Pipeline {
        info!(
            target: "kira_spliceqc::cli::run",
            "starting pipeline contract artifact generation"
        );
        pipeline_contract::write_pipeline_contract(
            &effective_out_dir,
            &context.stage0,
            &context.stage1,
            &context.stage4,
            &context.stage5,
            &context.stage6,
            context.stage8.as_ref(),
            &context.stage15,
            context.stage16.as_ref(),
            context.stage17.as_ref(),
            &context.stage18,
            &strata,
            &catalog,
        )?;
        info!(
            target: "kira_spliceqc::cli::run",
            "finished pipeline contract artifact generation"
        );
    }

    print!("{summary}");

    // Drop the mmap before removing the scratch directory it points at.
    drop(context);
    if cache_dir.is_dir()
        && let Err(err) = std::fs::remove_dir_all(&cache_dir)
    {
        warn!(
            target: "kira_spliceqc::cli::run",
            path = %cache_dir.display(),
            error = %err,
            "failed to remove expression cache directory"
        );
    }
    Ok(())
}

fn run_logged<T, F>(stage: usize, f: F) -> Result<T, SpliceQcError>
where
    F: FnOnce() -> Result<T, InputError>,
{
    info!(target: "kira_spliceqc::cli::run", "starting Stage {stage}");
    let start = Instant::now();
    let result = f().map_err(SpliceQcError::from)?;
    let elapsed = start.elapsed();
    let elapsed_fmt = format!("{elapsed:.2?}");
    info!(
        target: "kira_spliceqc::cli::run",
        "finished Stage {stage} (elapsed: {elapsed_fmt})"
    );
    Ok(result)
}

impl From<InputError> for SpliceQcError {
    fn from(err: InputError) -> Self {
        SpliceQcError::PipelineFailure(err.to_string())
    }
}
