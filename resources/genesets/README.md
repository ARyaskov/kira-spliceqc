# Geneset catalog

`splicing_genesets.tsv` is the default catalog embedded in the binary
(`--catalog PATH` substitutes another file; the hash of the catalog in use is
recorded in `summary.json.provenance`). It is versioned separately from the
tool: **catalog version 1.0.0** (the `# Splicing Genesets (v1)` header). Panel
membership changes only with a catalog major version; the pipeline contract
records `panel_version` for the stage-15 panels.

Format: tab-separated `geneset_id`, `axis`, `gene_symbol` and an optional
fourth `ensembl_id` column; `#` lines are comments; symbols are current HGNC
symbols (legacy aliases and case-insensitive matching are applied at load
time, see `src/genesets/aliases.rs`).

## Genesets

| geneset_id | axis | genes | what it represents | sources |
| --- | --- | --- | --- | --- |
| `U1_CORE` | CORE_SPLICEOSOME | SNRPC, SNRPA, SNRNP70 | U1 snRNP, 5' splice-site recognition | Wahl, Will & Lührmann 2009 Cell; Kondo et al. 2015 eLife |
| `U2_CORE` | CORE_SPLICEOSOME | SF3A1, SF3A2, SF3B1 | U2 snRNP SF3a / SF3b, branch-point recognition | Wahl et al. 2009; Cretu et al. 2016 Mol Cell |
| `U4_U6_CORE` | CORE_SPLICEOSOME | PRPF3, PRPF4, SNRNP200 | U4/U6 di-snRNP and Brr2 helicase | Wahl et al. 2009; Agafonov et al. 2016 Science |
| `U5_CORE` | CORE_SPLICEOSOME | PRPF8, EFTUD2, SNRNP40 | U5 snRNP core | Wahl et al. 2009; Grainger & Beggs 2005 RNA |
| `SF3B_AXIS` | CORE_SPLICEOSOME | SF3B1, SF3B2, SF3B3 | SF3b complex, target of SF3b inhibitors and hotspot mutations | Cretu et al. 2016; Darman et al. 2015 Cell Rep; Alsafadi et al. 2016 Nat Commun |
| `SRSF_SR` | REGULATORS | SRSF1, SRSF2, SRSF3 | SR proteins, exon definition | Long & Caceres 2009 Biochem J |
| `HNRNP` | REGULATORS | HNRNPA1, HNRNPA2B1, HNRNPC | hnRNPs, silencing / antagonists of SR proteins | Geuens, Bouhy & Timmerman 2016 Hum Genet |
| `MINOR_U12` | MINOR_SPLICEOSOME | SNRNP35, SNRNP48, ZRSR2 | U12-type minor spliceosome | Turunen et al. 2013 WIREs RNA; Madan et al. 2015 Nat Commun |
| `NMD_SURVEILLANCE` | SURVEILLANCE | UPF1, UPF2, SMG1 | nonsense-mediated decay core | Kurosaki, Popp & Maquat 2019 Nat Rev Mol Cell Biol |
| `CPA_3P_END` | SURVEILLANCE | CPSF1, CSTF1, CPSF3 | cleavage and polyadenylation | Shi & Manley 2015 Genes Dev |
| `TRANSCRIPTION_COUPLING` | TRANSCRIPTION | POLR2A, CDK9, SUPT5H, SUPT6H, AFF4, ELL | Pol II elongation machinery (co-transcriptional splicing) | Saldi et al. 2016 J Mol Biol; Herzel et al. 2017 Nat Rev Mol Cell Biol |
| `U2AF_AXIS` | U2AF | U2AF1, U2AF2, SF1, ZRSR2 | 3' splice-site recognition (U2AF heterodimer, SF1) | Wahl et al. 2009; Yoshida et al. 2011 Nature |
| `SPLICE_EA_PHASE` | ASSEMBLY_EA | SNRPC, SNRPA, SNRNP70, U2AF1, U2AF2, SF1 | E / A complex assembly | Wahl et al. 2009 |
| `SPLICE_B_PHASE` | ASSEMBLY_B | PRPF3, PRPF4, PRPF6, PRPF8, SNRNP200, EFTUD2 | tri-snRNP recruitment, B complex | Agafonov et al. 2016; Bertram et al. 2017 Nature |
| `SPLICE_CATALYTIC_PHASE` | ASSEMBLY_CAT | PRPF8, CDC5L, SLU7, DHX15 | catalytic activation, exon ligation, disassembly | Fica & Nagai 2017 Nat Struct Mol Biol |

Stage-15 panels (`SPLICEQC_INSTABILITY_PANEL_V1`, defined in
`src/metrics/splicing_instability/panels.rs`): `SPLICEOSOME_PANEL` (22 snRNP /
SF3 / U2AF genes), `SPLICING_RBP_PANEL` (hnRNPs, SR proteins, RBM10/RBM39),
`RLOOP_RESOLUTION_PANEL` (SETX, RNASEH1/2A, DHX9, DDX5, TOP1, ...),
`CONFLICT_RISK_PANEL` (TOP1, TOP2B, POLR2A, SUPT5H, SUPT6H), `NMD_PANEL` (UPF1-3B,
SMG1/5/6/7). Cell-cycle lists (`S_GENES`, `G2M_GENES`) follow Tirosh et al.
2016 Science with the Seurat 2019 symbol update.

## Proposing a panel

Open a "New panel" issue (template in `.github/ISSUE_TEMPLATE`). A panel is
accepted when it has a stated biological hypothesis, at least three genes
with literature support, a positive-control dataset in `benchmarks/datasets.tsv`
and a null-model result (`tests/null_model.rs`) showing no library-size
correlation. Panels are added in a catalog minor version; removals or
renamings wait for a catalog major version.
