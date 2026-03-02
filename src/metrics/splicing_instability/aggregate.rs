use crate::model::splicing_instability::SplicingInstabilityClusterStat;

pub fn aggregate_cluster_stats(
    _cluster_ids: Option<&[String]>,
    _sos: &[f32],
    _rlr: &[f32],
    _sii: &[f32],
    _splice_overload_high: &[bool],
    _rloop_risk_high: &[bool],
    _splicing_instability_high: &[bool],
    _genome_instability_splicing_flag: &[bool],
) -> Vec<SplicingInstabilityClusterStat> {
    // Junction/cluster-aware mode is optional and currently unavailable in kira-spliceqc runtime.
    // Keep deterministic empty output until cluster labels are plumbed into pipeline context.
    Vec::new()
}
