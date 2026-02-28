#[derive(Debug, Clone)]
pub struct MissplicingMetrics {
    pub b_core: Vec<f32>,
    pub b_u12: Vec<f32>,
    pub b_nmd: Vec<f32>,
    pub b_srhn: Vec<f32>,
    pub burden: Vec<f32>,
    pub burden_star: Vec<f32>,
}
