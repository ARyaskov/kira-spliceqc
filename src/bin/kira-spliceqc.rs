use std::path::PathBuf;

use clap::{Args, Parser, Subcommand, ValueEnum};
use kira_spliceqc::cli::config::{AnalysisMode, RunConfig, RunMode};
use kira_spliceqc::cli::run::{SpliceQcError, build_reference_file, run_pipeline};
use kira_spliceqc::validation::evaluate::{Pair, default_pairs, evaluate};
use kira_spliceqc::validation::simulate::{Effect, SimulationConfig, simulate};
use tracing_subscriber::EnvFilter;

#[derive(Parser)]
#[command(name = "kira-spliceqc", version, author)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Commands>,
    #[command(flatten)]
    pub run: RunArgs,
}

#[derive(Subcommand)]
pub enum Commands {
    Run(RunArgs),
    /// Reference-file commands.
    Reference(ReferenceCommand),
    /// Write a synthetic dataset with known truth (tier-1 validation).
    Simulate(SimulateArgs),
    /// Score a run's cells.tsv against a truth table (AUROC, AUPRC, flag
    /// precision/recall).
    Validate(ValidateArgs),
}

#[derive(Args, Clone)]
pub struct SimulateArgs {
    /// Output directory (10x matrix + layers + sj/ junctions + metadata.tsv + truth.tsv).
    #[arg(long)]
    pub out: PathBuf,
    #[arg(long, default_value_t = 2000)]
    pub n_cells: usize,
    #[arg(long, default_value_t = 800)]
    pub n_filler_genes: usize,
    #[arg(long, default_value_t = 80)]
    pub n_junction_genes: usize,
    #[arg(long, default_value_t = 0x5EED)]
    pub seed: u64,
    /// Fraction of cells with cryptic 3' splice-site usage (SF3B1-like).
    #[arg(long, default_value_t = 0.05)]
    pub cryptic_fraction: f64,
    /// Share of donor reads at the cryptic acceptor in those cells.
    #[arg(long, default_value_t = 0.15)]
    pub cryptic_ratio: f64,
    /// Fraction of cells with global intron retention.
    #[arg(long, default_value_t = 0.05)]
    pub ir_fraction: f64,
    #[arg(long, default_value_t = 2.0)]
    pub ir_fold: f64,
    /// Fraction of damaged cells (unspliced counts / 10).
    #[arg(long, default_value_t = 0.03)]
    pub damaged_fraction: f64,
    /// Fraction of cells with exon skipping.
    #[arg(long, default_value_t = 0.03)]
    pub skip_fraction: f64,
    #[arg(long, default_value_t = 0.25)]
    pub skip_ratio: f64,
}

#[derive(Args, Clone)]
pub struct ValidateArgs {
    /// Run output directory (standalone: contains cells.tsv; pipeline: <out>/kira-spliceqc).
    #[arg(long)]
    pub run: PathBuf,
    /// Truth table: barcode column + boolean truth columns (+ optional cell_type).
    #[arg(long)]
    pub truth: PathBuf,
    /// Scoring pairs `truth:metric[:flag[:sign]]`; defaults cover the
    /// simulate truth columns.
    #[arg(long)]
    pub pair: Vec<String>,
    /// Output JSON (a .md sibling is written too).
    #[arg(long)]
    pub out: PathBuf,
}

#[derive(Args, Clone)]
pub struct ReferenceCommand {
    #[command(subcommand)]
    pub action: ReferenceAction,
}

#[derive(Subcommand, Clone)]
pub enum ReferenceAction {
    /// Build a reference (ref.json) from a control dataset with
    /// spliced/unspliced layers; use it in later runs with --reference.
    Build(ReferenceBuildArgs),
}

#[derive(Args, Clone)]
pub struct ReferenceBuildArgs {
    /// Control dataset (10x directory or .h5ad) with spliced/unspliced layers.
    #[arg(long)]
    pub input: PathBuf,
    /// Output reference file (ref.json).
    #[arg(long)]
    pub out: PathBuf,
    #[arg(long)]
    pub layers: Option<PathBuf>,
    #[arg(long)]
    pub metadata: Option<PathBuf>,
    #[arg(long)]
    pub stratify_by: Option<String>,
    #[arg(long, default_value_t = 500)]
    pub min_counts: u64,
    #[arg(long, default_value_t = 200)]
    pub min_genes: u64,
    #[arg(long)]
    pub threads: Option<usize>,
}

#[derive(ValueEnum, Clone, Copy)]
pub enum ModeArg {
    Cell,
    Sample,
}

#[derive(ValueEnum, Clone, Copy)]
pub enum RunModeArg {
    Standalone,
    Pipeline,
}

#[derive(Args, Clone)]
pub struct RunArgs {
    #[arg(long)]
    pub input: Option<PathBuf>,
    #[arg(long)]
    pub out: Option<PathBuf>,
    #[arg(long)]
    pub cache: Option<PathBuf>,
    /// Spliced/unspliced layer source: a directory with spliced.mtx and
    /// unspliced.mtx (STARsolo Velocyto, kb-python) or an .h5ad with layers/.
    /// Auto-detected next to the input when omitted.
    #[arg(long)]
    pub layers: Option<PathBuf>,
    /// Junction count matrix directory (STARsolo Solo.out/SJ/<subset>:
    /// matrix.mtx + features.tsv + barcodes.tsv). Auto-detected as the SJ/
    /// sibling of a STARsolo Gene/ directory.
    #[arg(long)]
    pub junctions: Option<PathBuf>,
    /// Cell metadata table (barcode + columns, tab-separated, header line).
    /// Auto-detected as metadata.tsv[.gz] next to a 10x directory; .h5ad
    /// inputs use obs.
    #[arg(long)]
    pub metadata: Option<PathBuf>,
    /// Metadata column that defines reference strata (default: cell_type,
    /// then cluster aliases, else one global stratum).
    #[arg(long)]
    pub stratify_by: Option<String>,
    /// External reference (ref.json from `reference build`); Tier A
    /// deviations and flags are then relative to the reference strata.
    #[arg(long)]
    pub reference: Option<PathBuf>,
    /// Geneset catalog TSV (geneset_id, axis, gene_symbol[, ensembl_id]).
    /// Default: resources/genesets/splicing_genesets.tsv or the embedded copy.
    #[arg(long)]
    pub catalog: Option<PathBuf>,
    /// Cells with fewer UMIs are flagged LOW_DEPTH and excluded from
    /// reference norms (0 disables).
    #[arg(long, default_value_t = 500)]
    pub min_counts: u64,
    /// Cells with fewer detected genes are flagged LOW_DEPTH (0 disables).
    #[arg(long, default_value_t = 200)]
    pub min_genes: u64,
    #[arg(long, value_enum, default_value = "cell")]
    pub mode: ModeArg,
    #[arg(long)]
    pub json: bool,
    #[arg(long)]
    pub tsv: bool,
    #[arg(long)]
    pub extended: bool,
    #[arg(long)]
    pub threads: Option<usize>,
    /// Write experimental composite signatures (sis/class, SOS/RLR/SII, flags,
    /// cryptic risk, collapse) to the per-cell outputs. Implied by
    /// `--run-mode pipeline`.
    #[arg(long)]
    pub experimental_signatures: bool,
    #[arg(long, value_enum, default_value = "standalone")]
    pub run_mode: RunModeArg,
}

fn main() {
    let cli = Cli::parse();
    init_tracing();

    let result = match cli.command {
        Some(Commands::Run(args)) => execute_run(args),
        Some(Commands::Reference(cmd)) => match cmd.action {
            ReferenceAction::Build(args) => execute_reference_build(args),
        },
        Some(Commands::Simulate(args)) => execute_simulate(args),
        Some(Commands::Validate(args)) => execute_validate(args),
        None => execute_run(cli.run),
    };

    if let Err(err) = result {
        handle_error(err);
    }
}

fn init_tracing() {
    let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info"));
    tracing_subscriber::fmt().with_env_filter(filter).init();
}

fn handle_error(err: SpliceQcError) -> ! {
    match err {
        SpliceQcError::InvalidInput(msg) => {
            eprintln!("error: {msg}");
            std::process::exit(1);
        }
        SpliceQcError::Io(io) => {
            eprintln!("error: {io}");
            std::process::exit(1);
        }
        SpliceQcError::PipelineFailure(msg) => {
            eprintln!("error: {msg}");
            std::process::exit(2);
        }
    }
}

fn execute_simulate(args: SimulateArgs) -> Result<(), SpliceQcError> {
    let config = SimulationConfig {
        n_cells: args.n_cells,
        n_filler_genes: args.n_filler_genes,
        n_junction_genes: args.n_junction_genes,
        seed: args.seed,
        cryptic_fraction: args.cryptic_fraction,
        cryptic_ratio: args.cryptic_ratio,
        ir_fraction: args.ir_fraction,
        ir_fold: args.ir_fold,
        damaged_fraction: args.damaged_fraction,
        skip_fraction: args.skip_fraction,
        skip_ratio: args.skip_ratio,
    };
    let effects =
        simulate(&config, &args.out).map_err(|e| SpliceQcError::PipelineFailure(e.to_string()))?;
    let count = |e: Effect| effects.iter().filter(|x| **x == e).count();
    println!(
        "simulated {} cells into {}: cryptic {}, ir {}, damaged {}, skip {} (junctions in {}/sj, run with --junctions)",
        effects.len(),
        args.out.display(),
        count(Effect::Cryptic),
        count(Effect::Ir),
        count(Effect::Damaged),
        count(Effect::Skip),
        args.out.display()
    );
    Ok(())
}

fn execute_validate(args: ValidateArgs) -> Result<(), SpliceQcError> {
    let cells = if args.run.join("cells.tsv").is_file() {
        args.run.join("cells.tsv")
    } else {
        args.run.join("kira-spliceqc").join("cells.tsv")
    };
    let pairs = if args.pair.is_empty() {
        default_pairs()
    } else {
        args.pair
            .iter()
            .map(|p| Pair::parse(p))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| SpliceQcError::InvalidInput(e.to_string()))?
    };
    let report = evaluate(&cells, &args.truth, &pairs)
        .map_err(|e| SpliceQcError::PipelineFailure(e.to_string()))?;
    let json = serde_json::to_string_pretty(&report)
        .map_err(|e| SpliceQcError::PipelineFailure(e.to_string()))?;
    std::fs::write(&args.out, json)?;
    let md = args.out.with_extension("md");
    std::fs::write(&md, report.to_markdown())?;
    print!("{}", report.to_markdown());
    println!(
        "
written: {} and {}",
        args.out.display(),
        md.display()
    );
    Ok(())
}

fn execute_reference_build(args: ReferenceBuildArgs) -> Result<(), SpliceQcError> {
    let config = RunConfig {
        input: args.input,
        out_dir: args
            .out
            .parent()
            .map(|p| p.to_path_buf())
            .unwrap_or_default(),
        cache_path: None,
        layers: args.layers,
        junctions: None,
        metadata: args.metadata,
        stratify_by: args.stratify_by,
        reference: None,
        catalog: None,
        min_counts: args.min_counts,
        min_genes: args.min_genes,
        mode: AnalysisMode::Cell,
        run_mode: RunMode::Standalone,
        output_json: false,
        output_tsv: false,
        extended: false,
        threads: args.threads,
        experimental_signatures: false,
    };
    build_reference_file(config, &args.out)
}

fn execute_run(args: RunArgs) -> Result<(), SpliceQcError> {
    let input = args
        .input
        .ok_or_else(|| SpliceQcError::InvalidInput("--input is required".to_string()))?;
    let out = args
        .out
        .ok_or_else(|| SpliceQcError::InvalidInput("--out is required".to_string()))?;

    let config = RunConfig {
        input,
        out_dir: out,
        cache_path: args.cache,
        layers: args.layers,
        junctions: args.junctions,
        metadata: args.metadata,
        stratify_by: args.stratify_by,
        reference: args.reference,
        catalog: args.catalog,
        min_counts: args.min_counts,
        min_genes: args.min_genes,
        mode: match args.mode {
            ModeArg::Cell => AnalysisMode::Cell,
            ModeArg::Sample => AnalysisMode::Sample,
        },
        run_mode: match args.run_mode {
            RunModeArg::Standalone => RunMode::Standalone,
            RunModeArg::Pipeline => RunMode::Pipeline,
        },
        output_json: args.json,
        output_tsv: args.tsv,
        extended: args.extended,
        threads: args.threads,
        experimental_signatures: args.experimental_signatures,
    };
    run_pipeline(config)
}
