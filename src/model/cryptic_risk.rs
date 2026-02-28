#[derive(Debug, Clone)]
pub struct CrypticSplicingRiskMetrics {
    pub cryptic_risk: Vec<f32>,
    pub x_sr_hnrnp: Vec<f32>,
    pub x_entropy: Vec<f32>,
    pub x_nmd: Vec<f32>,
}
