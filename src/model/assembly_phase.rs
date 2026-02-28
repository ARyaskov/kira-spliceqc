#[derive(Debug, Clone)]
pub struct AssemblyPhaseImbalanceMetrics {
    pub z_ea: Vec<f32>,
    pub z_b: Vec<f32>,
    pub z_cat: Vec<f32>,
    pub ea_imbalance: Vec<f32>,
    pub b_imbalance: Vec<f32>,
    pub cat_imbalance: Vec<f32>,
}
