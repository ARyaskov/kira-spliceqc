pub const SPLICEQC_INSTABILITY_PANEL_V1: &str = "SPLICEQC_INSTABILITY_PANEL_V1";

pub const SPLICEOSOME_PANEL: &[&str] = &[
    "SNRPB", "SNRPD1", "SNRPD2", "SNRPD3", "SNRPE", "SNRPF", "SNRPG", "SF3A1", "SF3A2", "SF3A3",
    "SF3B1", "SF3B2", "SF3B3", "SF3B4", "SF3B5", "PRPF3", "PRPF4", "PRPF6", "PRPF8", "PRPF19",
    "U2AF1", "U2AF2",
];

pub const SPLICING_RBP_PANEL: &[&str] = &[
    "HNRNPA1",
    "HNRNPA2B1",
    "HNRNPC",
    "HNRNPK",
    "SRSF1",
    "SRSF2",
    "SRSF3",
    "SRSF6",
    "SRSF7",
    "RBM39",
    "RBM10",
    "RBM17",
];

pub const RLOOP_RESOLUTION_PANEL: &[&str] = &[
    "SETX", "DDX5", "DDX21", "DHX9", "RNASEH1", "RNASEH2A", "RNASEH2B", "RNASEH2C", "BRCA1",
    "BRCA2",
];

pub const CONFLICT_RISK_PANEL: &[&str] = &["TOP1", "TOP2A", "TOP2B", "POLR2A", "SUPT5H", "SUPT6H"];

pub const NMD_PANEL: &[&str] = &["UPF1", "UPF2", "UPF3B", "SMG1", "SMG5", "SMG6", "SMG7"];

pub const PANEL_TRIM_FRAC: f32 = 0.10;
pub const MIN_GENES_PER_PANEL_CELL: usize = 3;
pub const EPS_ROBUST: f32 = 1e-6;
