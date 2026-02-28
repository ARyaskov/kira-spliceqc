#[derive(Debug, Clone)]
pub struct SpliceosomeImbalanceMetrics {
    pub z_u1: Vec<f32>,
    pub z_u2: Vec<f32>,
    pub z_sf3b: Vec<f32>,
    pub z_srsf: Vec<f32>,
    pub z_hnrnp: Vec<f32>,
    pub z_u12: Vec<f32>,
    pub z_nmd: Vec<f32>,

    pub axis_sr_hnrnp: Vec<f32>,
    pub axis_u2_u1: Vec<f32>,
    pub axis_u12_major: Vec<f32>,
    pub axis_nmd: Vec<f32>,

    pub imbalance: Vec<f32>,
}
