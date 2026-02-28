#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SplicingTrajectoryClass {
    Adaptive,
    Degenerative,
    Oscillatory,
    Inconclusive,
}

#[derive(Debug, Clone)]
pub struct TimecourseSplicingMetrics {
    pub trajectory: SplicingTrajectoryClass,
    pub delta_sis: Vec<f32>,
    pub delta_entropy: Vec<f32>,
    pub delta_imbalance: Vec<f32>,
}
