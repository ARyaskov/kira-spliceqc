use std::path::PathBuf;

use clap::{Args, Parser, Subcommand, ValueEnum};
use kira_spliceqc::cli::config::{AnalysisMode, RunConfig, RunMode};
use kira_spliceqc::cli::run::{SpliceQcError, run_pipeline};
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
    #[arg(long, value_enum, default_value = "standalone")]
    pub run_mode: RunModeArg,
}

fn main() {
    let cli = Cli::parse();
    init_tracing();

    let run_args = match cli.command {
        Some(Commands::Run(args)) => args,
        None => cli.run,
    };

    let result = execute_run(run_args);

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
    };
    run_pipeline(config)
}
