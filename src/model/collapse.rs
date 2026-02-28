#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpliceosomeCollapseStatus {
    NoCollapse,
    Collapse,
    Inconclusive,
}

#[derive(Debug, Clone)]
pub struct SpliceosomeCollapseMetrics {
    pub collapse_status: Vec<SpliceosomeCollapseStatus>,
    pub core_suppression: Vec<bool>,
    pub high_imbalance: Vec<bool>,
    pub low_sis: Vec<bool>,
}
