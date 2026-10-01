# Metric cards

One card per production metric: what it measures, the formula, the range,
when it is undefined, the confounders it is known to track and the
literature. Deviations (`_dev`) and flags follow the rules in
[Reference strata](../reference.md). The full specification with every
constant is in [Full specification](specification.md).

## Tier A: direct measurements from spliced / unspliced layers

### `unspliced_fraction`

- **Measures**: share of a cell's layer UMIs that are unspliced (intronic).
- **Formula**: \\(U_c / (S_c + U_c)\\); ambiguous UMIs excluded. Wilson 95 %
  interval in `unspliced_fraction_ci_low` / `_ci_high`.
- **Range**: 0-1. Whole cells 0.1-0.3, nuclei 0.5-0.7 (10x 3'); protocol- and
  cell-type-dependent.
- **Undefined** when \\(S_c + U_c < 100\\).
- **Deviation**: logit scale against the stratum (`unspliced_fraction_dev`),
  with binomial plus overdispersion variance.
- **Flag**: `nuclear_fraction_flag`, deviation \\(\le -3\\) (fraction far
  *below* the stratum): cytoplasmic-only droplets, damaged cells, cell-free RNA.
- **Confounders**: protocol (nuclei vs cells), cell type (neurons high),
  quantification tool (velocyto vs STARsolo vs kb assign introns differently),
  library size through the interval width only.
- **Literature**: La Manno et al. 2018 Nature; Muskovic & Powell 2021 Genome
  Biology (DropletQC); Soneson et al. 2021 PLoS Comput Biol (quantification
  differences).

### `intron_retention_index`

- **Measures**: whether the cell retains introns more than its stratum, gene
  by gene.
- **Formula**: per gene \\(g\\) with \\(\ge 10\\) layer UMIs, the log2 ratio of
  the cell's beta-shrunk unspliced ratio to the stratum's pooled ratio
  (shrinkage \\(K = 10\\) pseudo-counts); precision-weighted mean over genes
  (weights capped at 50 UMIs) with a delta-method standard error.
  `ir_genes_used` and `ir_gene_dispersion` describe the support.
- **Range**: unbounded, centred on 0; \\(\pm 1\\) is a two-fold change.
- **Undefined** with fewer than 5 usable genes.
- **Deviation**: scaled deviation using the per-cell SE and the stratum's
  overdispersion, within library-size bins.
- **Flag**: `intron_retention_high`, deviation \\(\ge 3\\).
- **Confounders**: library size (handled by depth bins), nuclear fraction
  (a damaged cell is globally low, not gene-specifically high), 3' bias of the
  protocol.
- **Literature**: Middleton et al. 2017 Genome Biology (IRFinder);
  Pellagatti et al. 2018 Blood (IR in splicing-factor-mutant MDS); Dvinge &
  Bradley 2015 Genome Med.

## Tier B: direct measurements from junction counts

Defined for cells with at least 200 junction UMIs; ratios need at least 20
UMIs in the denominator.

### `cryptic_3ss_fraction`

- **Measures**: usage of unannotated acceptors 10-50 nt upstream of an
  annotated acceptor of the same donor: the SF3B1-mutant phenotype.
- **Formula**: cryptic UMIs / (cryptic + canonical-partner UMIs) over the
  affected donors.
- **Range**: 0-1; typically below 0.02 in wild-type cells.
- **Flag**: `cryptic_3ss_high`, logit deviation \\(\ge 3\\).
- **Confounders**: annotation completeness (an unannotated but genuine
  acceptor counts as cryptic everywhere, which the stratum reference absorbs),
  mapping artefacts near repeats.
- **Literature**: Darman et al. 2015 Cell Reports; Alsafadi et al. 2016
  Nature Communications; Nam et al. 2019 Nature (single-cell SF3B1 K700E).

### `exon_skip_fraction`

- **Measures**: skipping of annotated exons.
- **Formula**: skip UMIs / (skip + inclusion UMIs / 2), inclusion counted on
  both flanking junctions (rMATS junction-count PSI complement).
- **Range**: 0-1.
- **Flag**: `exon_skip_high`.
- **Confounders**: 3' coverage bias (fewer internal junctions), gene-length
  mix of the stratum.
- **Literature**: Shen et al. 2014 PNAS (rMATS); Wang et al. 2008 Nature.

### `splice_site_shift`

- **Measures**: whether the cell's choice among alternative donors /
  acceptors differs from the stratum, regardless of annotation (SpliZ idea).
- **Formula**: for every site with \\(\ge 2\\) partners, each UMI carries the
  rank of its partner in transcript direction; the cell's mean rank is
  standardized by the site's per-UMI mean and variance; the cell score is the
  median over its \\(\ge 5\\) sites with \\(\ge 3\\) UMIs.
- **Range**: unbounded, centred on 0; positive = downstream shift.
- **Flag**: `splice_site_shift_high`, deviation standardized within
  junction-depth bins (the raw score rises with depth).
- **Confounders**: junction depth, number of expressed multi-partner sites.
- **Literature**: Olivieri et al. 2022 Nature Methods (SpliZ).

### `unannotated_junction_fraction`

- **Measures**: share of junction UMIs on unannotated junctions; a global
  novelty / noise indicator, reported without a flag.

## Expression signatures (`*_expr`, any input)

Robust z-scores of geneset activity (trimmed mean of log CP10K minus a
control-gene background of 50 genes of matching mean expression per panel
gene, Tirosh et al. 2016) within stratum and library-size bin. They say that
the splicing machinery is expressed differently, not that splicing is
different: use them as covariates and hypotheses, not as evidence.

| column | geneset / panel |
| --- | --- |
| `spliceosome_core_expr` | snRNP / SF3 / U2AF panel |
| `splicing_rbp_expr` | hnRNPs, SR proteins, RBM10 / RBM39 |
| `rloop_resolution_expr` | SETX, RNASEH1 / 2A, DHX9, DDX5, TOP1, ... |
| `conflict_risk_expr` | TOP1, TOP2B, POLR2A, SUPT5H, SUPT6H |
| `nmd_factor_expr` | UPF1-3B, SMG1 / 5 / 6 / 7 |
| `regulator_entropy_expr`, `regulator_dispersion_expr` | entropy of SR / hnRNP regulator distribution |
| `missplicing_burden_expr`, `spliceosome_imbalance_expr`, `coupling_stress_expr`, `exon_definition_bias_expr`, `*_phase_imbalance_expr` | derived contrasts of catalog genesets |
| `s_score_expr`, `g2m_score_expr` | Tirosh 2016 cell-cycle lists |

Known confounders: cell cycle (spliceosome and R-loop panels rise in S / G2M;
`cycling` is reported next to them), proliferation and ribosome biogenesis
programmes, stress responses.

## Experimental composites

`sis`, `class`, `SOS`, `RLR`, `SII` and their flags combine the expression
signatures and are gated behind `--experimental-signatures` (implied in
pipeline mode) until validated on tiers 2-3; their flags are calibrated on the
null model like every other flag.
