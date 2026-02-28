#[derive(Debug, Clone)]
pub struct IsoformDispersionMetrics {
    pub entropy: Vec<f32>,
    pub dispersion: Vec<f32>,
    pub z_entropy: Vec<f32>,
}
