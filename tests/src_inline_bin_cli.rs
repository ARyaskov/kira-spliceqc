use clap::Parser;

#[path = "../src/bin/kira-spliceqc.rs"]
mod bin_main;

#[test]
fn run_mode_defaults_to_standalone() {
    let cli =
        bin_main::Cli::try_parse_from(["kira-spliceqc", "--input", "in", "--out", "out"]).unwrap();
    let run = match cli.command {
        Some(bin_main::Commands::Run(args)) => args,
        None => cli.run,
    };
    assert!(matches!(run.run_mode, bin_main::RunModeArg::Standalone));
}

#[test]
fn run_mode_accepts_pipeline() {
    let cli = bin_main::Cli::try_parse_from([
        "kira-spliceqc",
        "--input",
        "in",
        "--out",
        "out",
        "--run-mode",
        "pipeline",
    ])
    .unwrap();
    let run = match cli.command {
        Some(bin_main::Commands::Run(args)) => args,
        None => cli.run,
    };
    assert!(matches!(run.run_mode, bin_main::RunModeArg::Pipeline));
}
