#[derive(Debug, Clone)]
pub struct ExonIntronDefinitionMetrics {
    pub exon_definition_bias: Vec<f32>,
    pub z_srsf: Vec<f32>,
    pub z_u2af: Vec<f32>,
    pub z_hnrnp: Vec<f32>,
}
