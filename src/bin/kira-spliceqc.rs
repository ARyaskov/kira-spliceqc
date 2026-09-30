use std::path::PathBuf;

use clap::{Args, Parser, Subcommand, ValueEnum};
use kira_spliceqc::cli::config::{AnalysisMode, RunConfig, RunMode};
use kira_spliceqc::cli::run::{SpliceQcError, build_reference_file, run_pipeline};
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

fn execute_reference_build(args: ReferenceBuildArgs) -> Result<(), SpliceQcError> {
    let config = RunConfig {
        input: args.input,
        out_dir: args.out.parent().map(|p| p.to_path_buf()).unwrap_or_default(),
        cache_path: None,
        layers: args.layers,
        metadata: args.metadata,
        stratify_by: args.stratify_by,
        reference: None,
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
        metadata: args.metadata,
        stratify_by: args.stratify_by,
        reference: args.reference,
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
