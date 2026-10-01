//! Tier B: per-cell splicing metrics from junction counts (input level L2).
//!
//! The annotation comes from the junctions themselves: every junction the
//! aligner marked `annotated` contributes an annotated donor and acceptor,
//! so no GTF is needed. Unannotated junctions are then classified:
//!
//! - **cryptic 3' splice site**: the donor is annotated and the acceptor lies
//!   `CRYPTIC_MIN..=CRYPTIC_MAX` nt upstream (toward the donor) of an
//!   annotated acceptor of the same strand. This is the SF3B1-mutant
//!   phenotype (Darman et al. 2015 Cell Reports; Alsafadi et al. 2016 Nature
//!   Communications). The canonical partner is the junction donor -> that
//!   annotated acceptor; the per-cell fraction is cryptic / (cryptic +
//!   canonical), i.e. cryptic usage at affected donors.
//! - **exon skipping**: annotated junctions `D -> A1` and `D2 -> A` exist
//!   with `A1` before `D2`, i.e. the junction `D -> A` skips at least one
//!   annotated exon (annotated skip junctions count too). Inclusion UMIs are
//!   those of all such partner junctions; the per-cell fraction is
//!   skip / (skip + inclusion), the usual PSI complement (Shen et al. 2014
//!   PNAS, rMATS event definition).
//! - everything else unannotated is **novel**.
//!
//! `splice_site_shift` follows the idea of SpliZ (Olivieri et al. 2022
//! Nature Methods): for every donor with >= 2 acceptors (and every acceptor
//! with >= 2 donors) the acceptor rank of each UMI is a position variable;
//! a cell's mean rank at a site is compared with the reference stratum's
//! per-UMI mean and variance, `z = (r_c - r_s) / sqrt(v_s / n_c)`, and the
//! cell's score is `median_sites |z| / 0.6745`, ~1 under the null.
//!
//! All per-cell proportions and their flags use the reference machinery of
//! Tier A (logit deviations with overdispersion, BH-adjusted flags).

use std::collections::BTreeMap;

use ahash::AHashMap;
use rayon::prelude::*;

use crate::expression::junctions::JunctionSet;
use crate::io::junctions::Strand;
use crate::model::junctions::JunctionMetrics;
use crate::reference::{
    Strata, StratumStat, apply_proportion_norms, flag_outliers, logit_deviation_by_stratum,
    proportion_norms, robust_z_by_stratum, robust_z_by_stratum_and_depth,
};
use crate::stats::robust::median;

/// Cells with fewer junction UMIs have every Tier B value undefined.
pub const MIN_JUNCTION_UMIS: u64 = 200;
/// A ratio (cryptic, skip) is defined when its denominator has this many UMIs.
pub const MIN_RATIO_UMIS: u64 = 20;
/// Cryptic acceptor window upstream of the annotated acceptor, in nt.
pub const CRYPTIC_MIN: u64 = 10;
pub const CRYPTIC_MAX: u64 = 50;
/// A site group counts for a cell when it has this many UMIs there.
pub const MIN_SITE_UMIS: u32 = 3;
/// A cell's shift score needs this many site groups.
pub const MIN_SITE_GROUPS: usize = 5;
/// Median of |N(0,1)|.
const MEDIAN_ABS_NORMAL: f32 = 0.6745;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Class {
    Annotated,
    CrypticAcceptor,
    Novel,
}

/// Static classification of every junction.
pub struct JunctionAnnotation {
    class: Vec<Class>,
    /// Junction is the canonical partner of at least one cryptic junction.
    canonical_partner: Vec<bool>,
    /// Junction skips at least one annotated exon.
    skip: Vec<bool>,
    /// Junction is an inclusion partner of at least one skip junction.
    inclusion_partner: Vec<bool>,
    /// `(group id, rank)` memberships per junction (donor-anchored and
    /// acceptor-anchored alternative-site groups).
    site_membership: Vec<Vec<(u32, u32)>>,
    n_groups: usize,
}

impl JunctionAnnotation {
    pub fn n_cryptic(&self) -> usize {
        self.class
            .iter()
            .filter(|c| **c == Class::CrypticAcceptor)
            .count()
    }
    pub fn n_skip(&self) -> usize {
        self.skip.iter().filter(|s| **s).count()
    }
}

pub fn annotate(set: &JunctionSet) -> JunctionAnnotation {
    let junctions = &set.junctions;
    let n = junctions.len();
    let key = |chrom: &str, strand: Strand| (chrom.to_string(), strand.as_char());

    // Annotated sites per (chrom, strand): donors and sorted acceptors, plus
    // the junction index for (donor, acceptor) lookups.
    let mut donors: AHashMap<(String, char), BTreeMap<u64, Vec<u32>>> = AHashMap::new();
    let mut acceptors: AHashMap<(String, char), BTreeMap<u64, Vec<u32>>> = AHashMap::new();
    for (i, j) in junctions.iter().enumerate() {
        if !j.annotated || j.strand == Strand::Unknown {
            continue;
        }
        let k = key(&j.chrom, j.strand);
        donors
            .entry(k.clone())
            .or_default()
            .entry(j.donor())
            .or_default()
            .push(i as u32);
        acceptors
            .entry(k)
            .or_default()
            .entry(j.acceptor())
            .or_default()
            .push(i as u32);
    }

    let mut class = vec![Class::Novel; n];
    let mut canonical_partner = vec![false; n];
    for (i, j) in junctions.iter().enumerate() {
        if j.annotated {
            class[i] = Class::Annotated;
            continue;
        }
        if j.strand == Strand::Unknown {
            continue;
        }
        let k = key(&j.chrom, j.strand);
        let Some(donor_map) = donors.get(&k) else {
            continue;
        };
        let Some(by_donor) = donor_map.get(&j.donor()) else {
            continue;
        };
        // Annotated acceptor downstream of this acceptor by 10..=50 nt.
        let acc = j.acceptor();
        let (lo, hi) = match j.strand {
            Strand::Plus => (acc + CRYPTIC_MIN, acc + CRYPTIC_MAX),
            Strand::Minus => (
                acc.saturating_sub(CRYPTIC_MAX),
                acc.saturating_sub(CRYPTIC_MIN),
            ),
            Strand::Unknown => unreachable!(),
        };
        // Canonical partner: an annotated junction from the same donor whose
        // acceptor is in the window.
        let partner = by_donor
            .iter()
            .copied()
            .find(|&p| (lo..=hi).contains(&junctions[p as usize].acceptor()));
        if let Some(p) = partner {
            class[i] = Class::CrypticAcceptor;
            canonical_partner[p as usize] = true;
        }
    }

    // Exon skipping: D -> A skips an exon when annotated D -> A1 and D2 -> A
    // exist with A1 before D2 (transcript direction).
    let mut skip = vec![false; n];
    let mut inclusion_partner = vec![false; n];
    for (i, j) in junctions.iter().enumerate() {
        if j.strand == Strand::Unknown {
            continue;
        }
        let k = key(&j.chrom, j.strand);
        let (Some(donor_map), Some(acceptor_map)) = (donors.get(&k), acceptors.get(&k)) else {
            continue;
        };
        let (Some(from_donor), Some(to_acceptor)) =
            (donor_map.get(&j.donor()), acceptor_map.get(&j.acceptor()))
        else {
            continue;
        };
        let d = j.donor();
        let a = j.acceptor();
        let mut partners: Vec<u32> = Vec::new();
        for &p in from_donor {
            let a1 = junctions[p as usize].acceptor();
            // A1 strictly between D and A
            let inside = match j.strand {
                Strand::Plus => a1 > d && a1 < a,
                _ => a1 < d && a1 > a,
            };
            if inside {
                for &q in to_acceptor {
                    let d2 = junctions[q as usize].donor();
                    let after_a1 = match j.strand {
                        Strand::Plus => d2 > a1 && d2 < a,
                        _ => d2 < a1 && d2 > a,
                    };
                    if after_a1 {
                        partners.push(p);
                        partners.push(q);
                    }
                }
            }
        }
        if !partners.is_empty() {
            skip[i] = true;
            for p in partners {
                inclusion_partner[p as usize] = true;
            }
        }
    }

    // Alternative-site groups: all junctions (any annotation) sharing a
    // donor, ranked by acceptor position in transcript direction; likewise
    // sharing an acceptor, ranked by donor.
    let mut by_donor_all: AHashMap<(String, char, u64), Vec<u32>> = AHashMap::new();
    let mut by_acceptor_all: AHashMap<(String, char, u64), Vec<u32>> = AHashMap::new();
    for (i, j) in junctions.iter().enumerate() {
        if j.strand == Strand::Unknown {
            continue;
        }
        by_donor_all
            .entry((j.chrom.clone(), j.strand.as_char(), j.donor()))
            .or_default()
            .push(i as u32);
        by_acceptor_all
            .entry((j.chrom.clone(), j.strand.as_char(), j.acceptor()))
            .or_default()
            .push(i as u32);
    }
    let mut site_membership: Vec<Vec<(u32, u32)>> = vec![Vec::new(); n];
    let mut n_groups = 0u32;
    let mut groups: Vec<(Vec<u32>, bool)> = Vec::new();
    for (_, members) in by_donor_all {
        if members.len() >= 2 {
            groups.push((members, true));
        }
    }
    for (_, members) in by_acceptor_all {
        if members.len() >= 2 {
            groups.push((members, false));
        }
    }
    // Deterministic group order.
    groups.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    for (mut members, donor_anchored) in groups {
        let strand = junctions[members[0] as usize].strand;
        members.sort_by_key(|&m| {
            let j = &junctions[m as usize];
            let pos = if donor_anchored {
                j.acceptor()
            } else {
                j.donor()
            };
            // transcript direction: ascending on +, descending on -
            match strand {
                Strand::Minus => u64::MAX - pos,
                _ => pos,
            }
        });
        for (rank, m) in members.iter().enumerate() {
            site_membership[*m as usize].push((n_groups, rank as u32));
        }
        n_groups += 1;
    }

    JunctionAnnotation {
        class,
        canonical_partner,
        skip,
        inclusion_partner,
        site_membership,
        n_groups: n_groups as usize,
    }
}

struct CellSums {
    total: u64,
    annotated: u64,
    cryptic: u64,
    canonical: u64,
    skip: u64,
    inclusion: u64,
    /// group id -> (umis, sum rank * umis)
    sites: Vec<(u32, u32, f64)>,
}

fn cell_sums(set: &JunctionSet, ann: &JunctionAnnotation, cell: usize) -> CellSums {
    let (rows, vals) = set.counts.cell(cell);
    let mut s = CellSums {
        total: 0,
        annotated: 0,
        cryptic: 0,
        canonical: 0,
        skip: 0,
        inclusion: 0,
        sites: Vec::new(),
    };
    let mut site_acc: AHashMap<u32, (u32, f64)> = AHashMap::new();
    for (&j, &k) in rows.iter().zip(vals) {
        let j = j as usize;
        let k64 = k as u64;
        s.total += k64;
        match ann.class[j] {
            Class::Annotated => s.annotated += k64,
            Class::CrypticAcceptor => s.cryptic += k64,
            Class::Novel => {}
        }
        if ann.canonical_partner[j] {
            s.canonical += k64;
        }
        if ann.skip[j] {
            s.skip += k64;
        }
        if ann.inclusion_partner[j] {
            s.inclusion += k64;
        }
        for &(g, rank) in &ann.site_membership[j] {
            let e = site_acc.entry(g).or_insert((0, 0.0));
            e.0 += k;
            e.1 += rank as f64 * k as f64;
        }
    }
    s.sites = site_acc
        .into_iter()
        .map(|(g, (n, sr))| (g, n, sr))
        .collect();
    s.sites.sort_unstable_by_key(|t| t.0);
    s
}

fn ratio(num: u64, denom: u64) -> f32 {
    if denom >= MIN_RATIO_UMIS {
        num as f32 / denom as f32
    } else {
        f32::NAN
    }
}

pub fn compute(set: &JunctionSet, strata: &Strata) -> JunctionMetrics {
    let n_cells = set.n_cells();
    debug_assert_eq!(strata.n_cells(), n_cells);
    let ann = annotate(set);

    let sums: Vec<CellSums> = (0..n_cells)
        .into_par_iter()
        .map(|c| cell_sums(set, &ann, c))
        .collect();

    let mut junction_umis = Vec::with_capacity(n_cells);
    let mut annotated_umis = Vec::with_capacity(n_cells);
    let mut unannotated_fraction = Vec::with_capacity(n_cells);
    let mut cryptic_umis = Vec::with_capacity(n_cells);
    let mut cryptic_fraction = Vec::with_capacity(n_cells);
    let mut cryptic_trials = Vec::with_capacity(n_cells);
    let mut skip_umis = Vec::with_capacity(n_cells);
    let mut skip_fraction = Vec::with_capacity(n_cells);
    let mut skip_trials = Vec::with_capacity(n_cells);
    let mut undefined_cells = 0usize;
    for s in &sums {
        let defined = s.total >= MIN_JUNCTION_UMIS;
        if !defined {
            undefined_cells += 1;
        }
        junction_umis.push(s.total);
        annotated_umis.push(s.annotated);
        unannotated_fraction.push(if defined {
            (s.total - s.annotated) as f32 / s.total as f32
        } else {
            f32::NAN
        });
        cryptic_umis.push(s.cryptic);
        let ct = s.cryptic + s.canonical;
        cryptic_fraction.push(if defined {
            ratio(s.cryptic, ct)
        } else {
            f32::NAN
        });
        cryptic_trials.push(ct);
        skip_umis.push(s.skip);
        let st = s.skip + s.inclusion;
        skip_fraction.push(if defined { ratio(s.skip, st) } else { f32::NAN });
        skip_trials.push(st);
    }

    // Site-usage shift: reference per stratum from the stratum's non-excluded cells.
    let members = strata.members();
    let n_groups = ann.n_groups;
    let mut stratum_stats: Vec<Vec<(f64, f64, f64)>> = Vec::with_capacity(members.len()); // (n, sum r, sum r^2) per group
    for cells in &members {
        let mut acc = vec![(0f64, 0f64, 0f64); n_groups];
        for &c in cells {
            if sums[c].total < MIN_JUNCTION_UMIS {
                continue;
            }
            let (rows, vals) = set.counts.cell(c);
            for (&j, &k) in rows.iter().zip(vals) {
                for &(g, rank) in &ann.site_membership[j as usize] {
                    let e = &mut acc[g as usize];
                    e.0 += k as f64;
                    e.1 += rank as f64 * k as f64;
                    e.2 += (rank as f64).powi(2) * k as f64;
                }
            }
        }
        stratum_stats.push(acc);
    }
    let shift_rows: Vec<(f32, u32)> = (0..n_cells)
        .into_par_iter()
        .map(|c| {
            let s = &sums[c];
            if s.total < MIN_JUNCTION_UMIS {
                return (f32::NAN, 0);
            }
            let stats = &stratum_stats[strata.labels[c] as usize];
            let mut zs: Vec<f32> = Vec::new();
            for &(g, n, sr) in &s.sites {
                if n < MIN_SITE_UMIS {
                    continue;
                }
                let (tn, tr, tr2) = stats[g as usize];
                if tn < 2.0 * MIN_SITE_UMIS as f64 {
                    continue;
                }
                let mean = tr / tn;
                let var = (tr2 / tn - mean * mean).max(0.0);
                if var <= 0.0 {
                    continue;
                }
                let r_c = sr / n as f64;
                let z = (r_c - mean) / (var / n as f64).sqrt();
                zs.push(z.abs() as f32);
            }
            let used = zs.len() as u32;
            if zs.len() < MIN_SITE_GROUPS {
                (f32::NAN, used)
            } else {
                (median(&zs) / MEDIAN_ABS_NORMAL, used)
            }
        })
        .collect();
    let mut splice_site_shift = Vec::with_capacity(n_cells);
    let mut site_groups_used = Vec::with_capacity(n_cells);
    for (v, u) in shift_rows {
        splice_site_shift.push(v);
        site_groups_used.push(u);
    }

    // Deviations and flags.
    let (_, cryptic_reference) =
        logit_deviation_by_stratum(&cryptic_fraction, &cryptic_trials, strata);
    let cryptic_norms: Vec<_> = proportion_norms(&cryptic_fraction, &cryptic_trials, strata)
        .into_iter()
        .map(Some)
        .collect();
    let cryptic_dev = apply_proportion_norms(
        &cryptic_fraction,
        &cryptic_trials,
        &strata.labels,
        &cryptic_norms,
    );
    let cryptic_high = flag_outliers(&cryptic_dev, strata, 1.0);

    let (_, skip_reference) = logit_deviation_by_stratum(&skip_fraction, &skip_trials, strata);
    let skip_norms: Vec<_> = proportion_norms(&skip_fraction, &skip_trials, strata)
        .into_iter()
        .map(Some)
        .collect();
    let skip_dev =
        apply_proportion_norms(&skip_fraction, &skip_trials, &strata.labels, &skip_norms);
    let skip_high = flag_outliers(&skip_dev, strata, 1.0);

    // The shift score rises with junction depth (more sites, more UMIs per
    // site), so its deviation is standardized within stratum and
    // junction-depth bin like every depth-sensitive metric.
    let (_, shift_reference) = robust_z_by_stratum(&splice_site_shift, strata);
    let (shift_dev, _) = robust_z_by_stratum_and_depth(&splice_site_shift, strata, &junction_umis);
    let shift_high = flag_outliers(&shift_dev, strata, 1.0);

    JunctionMetrics {
        source: set.source.clone(),
        n_junctions: set.n_junctions(),
        n_annotated: set.junctions.iter().filter(|j| j.annotated).count(),
        n_cryptic_acceptor_junctions: ann.n_cryptic(),
        n_skip_junctions: ann.n_skip(),
        n_site_groups: n_groups,
        min_junction_umis: MIN_JUNCTION_UMIS,
        min_ratio_umis: MIN_RATIO_UMIS,
        cells_without_junctions: set.cells_without_junctions,
        undefined_cells,
        junction_umis,
        annotated_umis,
        unannotated_junction_fraction: unannotated_fraction,
        cryptic_3ss_umis: cryptic_umis,
        cryptic_3ss_fraction: cryptic_fraction,
        cryptic_3ss_fraction_dev: cryptic_dev,
        cryptic_3ss_high: cryptic_high,
        exon_skip_umis: skip_umis,
        exon_skip_fraction: skip_fraction,
        exon_skip_fraction_dev: skip_dev,
        exon_skip_high: skip_high,
        splice_site_shift,
        splice_site_shift_dev: shift_dev,
        splice_site_shift_high: shift_high,
        site_groups_used,
        cryptic_reference: as_stats(cryptic_reference),
        skip_reference: as_stats(skip_reference),
        shift_reference: as_stats(shift_reference),
    }
}

fn as_stats(v: Vec<StratumStat>) -> Vec<StratumStat> {
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expression::LayerMatrix;
    use crate::io::junctions::Junction;

    fn j(chrom: &str, start: u64, end: u64, strand: Strand, annotated: bool) -> Junction {
        Junction {
            chrom: chrom.to_string(),
            start,
            end,
            strand,
            annotated,
        }
    }

    /// Three-exon gene on +: exon1 1..100, exon2 201..300, exon3 401..500.
    /// Junctions: 0 = 101-200 (annotated), 1 = 301-400 (annotated),
    /// 2 = 101-400 skip (unannotated), 3 = 101-180 cryptic acceptor (-20 nt),
    /// 4 = 101-195 (-5 nt, outside the window -> novel), 5 = 150-200 novel donor.
    fn junctions() -> Vec<Junction> {
        vec![
            j("chr1", 101, 200, Strand::Plus, true),
            j("chr1", 301, 400, Strand::Plus, true),
            j("chr1", 101, 400, Strand::Plus, false),
            j("chr1", 101, 180, Strand::Plus, false),
            j("chr1", 101, 195, Strand::Plus, false),
            j("chr1", 150, 200, Strand::Plus, false),
        ]
    }

    fn set(cells: Vec<Vec<(u32, u32)>>) -> JunctionSet {
        let n_cells = cells.len();
        let mut triplets = Vec::new();
        for (c, entries) in cells.into_iter().enumerate() {
            for (jid, k) in entries {
                triplets.push((jid, c as u32, k));
            }
        }
        JunctionSet {
            junctions: junctions(),
            counts: LayerMatrix::from_triplets(6, n_cells, triplets),
            source: "test".to_string(),
            cells_without_junctions: 0,
        }
    }

    #[test]
    fn classification_from_annotated_junctions() {
        let s = set(vec![vec![]]);
        let ann = annotate(&s);
        assert_eq!(ann.class[0], Class::Annotated);
        assert_eq!(ann.class[2], Class::Novel); // skip junction: not cryptic
        assert!(ann.skip[2]);
        assert!(ann.inclusion_partner[0] && ann.inclusion_partner[1]);
        assert_eq!(ann.class[3], Class::CrypticAcceptor);
        assert!(ann.canonical_partner[0]);
        assert_eq!(ann.class[4], Class::Novel); // 5 nt: outside 10..50
        assert_eq!(ann.class[5], Class::Novel);
        assert!(!ann.skip[5]);
        // Donor 101 has acceptors 180, 195, 200, 400 -> one donor group;
        // acceptor 200 has donors 101, 150 and acceptor 400 has donors 101
        // (skip) and 301 -> two acceptor groups.
        assert_eq!(ann.n_groups, 3);
    }

    #[test]
    fn per_cell_fractions_and_flags() {
        // 100 normal cells: 150 UMIs on each annotated junction, 2 cryptic, 3 skip;
        // 3 mutant-like cells: 40 cryptic UMIs; 2 skip-heavy cells: 60 skip UMIs.
        let mut cells = Vec::new();
        for c in 0..105 {
            let jitter = (c % 4) as u32;
            let cryptic = if c < 3 { 40 } else { 2 + jitter % 2 };
            let skip = if (3..5).contains(&c) {
                60
            } else {
                3 + jitter % 3
            };
            cells.push(vec![
                (0, 150 + jitter),
                (1, 150 + jitter),
                (2, skip),
                (3, cryptic),
                (5, 10),
            ]);
        }
        let s = set(cells);
        let m = compute(&s, &Strata::global(105));
        assert_eq!(m.undefined_cells, 0);
        assert_eq!(m.n_cryptic_acceptor_junctions, 1);
        assert_eq!(m.n_skip_junctions, 1);
        assert!((m.cryptic_3ss_fraction[0] - 40.0 / 190.0).abs() < 1e-5);
        assert!(m.cryptic_3ss_fraction[10] < 0.03);
        for c in 0..3 {
            assert!(
                m.cryptic_3ss_high[c],
                "cell {c} dev {}",
                m.cryptic_3ss_fraction_dev[c]
            );
        }
        assert_eq!(m.cryptic_3ss_high.iter().filter(|f| **f).count(), 3);
        for c in 3..5 {
            assert!(
                m.exon_skip_high[c],
                "cell {c} dev {}",
                m.exon_skip_fraction_dev[c]
            );
        }
        assert_eq!(m.exon_skip_high.iter().filter(|f| **f).count(), 2);
        assert!(m.exon_skip_fraction[3] > 0.15);
        assert!(m.unannotated_junction_fraction[10] < 0.1);
        // Too few site groups for a shift score in this toy (2 groups < 5).
        assert!(m.splice_site_shift[0].is_nan());
    }

    #[test]
    fn shallow_cells_are_undefined() {
        let s = set(vec![vec![(0, 10), (1, 10)], vec![(0, 150), (1, 150)]]);
        let m = compute(&s, &Strata::global(2));
        assert!(m.unannotated_junction_fraction[0].is_nan());
        assert!(m.cryptic_3ss_fraction[0].is_nan());
        assert_eq!(m.undefined_cells, 1);
        // Defined cell but no cryptic/canonical UMIs beyond the floor -> ratio defined (0 / 300 >= 20).
        assert_eq!(m.cryptic_3ss_fraction[1], 0.0);
    }
}
