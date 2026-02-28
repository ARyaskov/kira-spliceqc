#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpliceIntegrityClass {
    Intact,
    Stressed,
    Impaired,
    Broken,
}

#[derive(Debug, Clone)]
pub struct SpliceIntegrityMetrics {
    pub sis: Vec<f32>,
    pub class: Vec<SpliceIntegrityClass>,
    pub p_missplicing: Vec<f32>,
    pub p_imbalance: Vec<f32>,
    pub p_entropy_z: Vec<f32>,
    pub p_entropy_abs: Vec<f32>,
}
