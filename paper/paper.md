---
title: 'kira-spliceqc: deterministic, explainable splicing quality control for single-cell RNA-seq'
tags:
  - Rust
  - single-cell RNA-seq
  - splicing
  - quality control
  - intron retention
authors:
  - name: Andrei Riaskov
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Independent researcher
    index: 1
date: 1 October 2026
bibliography: paper.bib
---

# Summary

`kira-spliceqc` is a command-line tool that scores every cell of a single-cell
RNA-seq dataset for splicing abnormalities and reports calibrated flags next to
the raw measurements they rest on. It consumes the three levels of evidence a
typical processing run already produces: gene counts (10x MatrixMarket, `.h5ad`
or a shared binary cache), spliced / unspliced count layers (kb-python, STARsolo
Velocyto, `.h5ad` layers) and STARsolo junction count matrices. From these it
derives, per cell, an unspliced fraction with a Wilson interval, a beta-shrunk
intron retention index, cryptic 3' splice-site usage of the kind caused by
hotspot mutations in *SF3B1* [@darman2015; @alsafadi2016], exon skipping in the
rMATS event definition [@shen2014], a SpliZ-like splice-site shift score
[@olivieri2022], and a set of expression signatures of the splicing machinery
(spliceosome, SR proteins, hnRNPs, minor spliceosome, NMD) with a control-gene
depth correction [@tirosh2016]. Every deviation is computed against reference
strata (cell type, cluster or an external control dataset) and library-size
bins; flags are false-discovery-rate controlled and calibrated on a Poisson null
model to at most one per cent of cells. The implementation is deterministic
(byte-identical outputs across runs and thread counts), single-binary, and
writes a MultiQC custom-content table, a pipeline-contract table and full
provenance (parameters, catalog and reference hashes, constants).

# Statement of need

Splicing defects are a hallmark of myelodysplastic syndromes, chronic
lymphocytic leukaemia, uveal melanoma and several solid tumours, and nuclear
RNA leakage from damaged cells is a dominant technical artefact of droplet
single-cell protocols. Both appear in routine single-cell data, yet standard QC
(mitochondrial fraction, library size, doublet calls) does not look at splicing.
Existing tools cover single aspects: velocyto and scVelo [@lamanno2018;
@bergen2020] estimate unspliced fractions for dynamics rather than QC, DropletQC
[@muskovic2021] flags damaged cells from the nuclear fraction only, SpliZ
[@olivieri2022] and scQuint [@benegas2022] quantify splice-site usage but need
alignments and per-gene modelling, and IRFinder [@middleton2017] works on bulk
pseudobulks. None of them combines the three levels of evidence per cell,
exposes calibrated flags with reference strata, or runs as a dependency-free
step of a pipeline.

`kira-spliceqc` fills that gap for pipeline authors and analysts: a single
command adds a splicing QC layer to any dataset with the evidence it has
(expression only, layers, or junctions), labels what is measured versus what is
an expression proxy (`_expr` suffix; composite indices are gated behind an
experimental flag until validated), and ships a validation package. A built-in
simulator spikes cryptic splice sites, intron retention, damaged cells and exon
skipping into a Poisson dataset with known truth; a `validate` subcommand
reports AUROC and flag precision / recall / FPR; continuous integration enforces
AUROC $\geq$ 0.95 and FPR $\leq$ 1 % for every effect and the absence of
library-size correlation on the null model. Positive-control benchmarks on public
data (Genotyping of Transcriptomes CLL and MDS cohorts [@nam2019], genome-wide
Perturb-seq [@replogle2022]) and agreement benchmarks against velocyto, DropletQC,
SpliZ and IRFinder are specified in the repository and will be published with a
DOI once run.

# Acknowledgements

The author thanks the maintainers of STARsolo, kb-python and the scverse
ecosystem, whose output formats this tool reads.

# References
