#[derive(Debug, Clone)]
pub struct SplicingNoiseMetrics {
    pub noise_index: f32,
    pub per_geneset_noise: Vec<(String, f32)>,
}
