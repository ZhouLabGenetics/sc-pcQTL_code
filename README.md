# sc-pcQTL

This repository contains the analysis code for sc-pcQTL, from single-cell
gene-pair screening through cluster-PC QTL mapping and downstream trait
analyses. Generated results and controlled-access inputs are not included.

## Workflow

1. `modules/01_pairwise_hurdle_clustering`: hurdle-based gene-pair screening
   and local cluster calling.
2. `modules/02_cluster_pc_qtl`: cluster principal components and SAIGE-QTL
   cis-pcQTL mapping.
3. `modules/03_pcqtl_summary_and_followup`: pcQTL summaries, ACAT, and targeted
   GIMAP SMR follow-up.
4. `modules/04_cluster_annotation_enrichment`: neighboring-null annotation
   enrichment.
5. `modules/05_simulations`: model-based and permutation-based simulations.
6. `modules/06_finngen_susie_coloc`: OneK1K QTL fine-mapping and FinnGen
   colocalization, including the observed-variant `coloc.abf` sensitivity run.
7. `modules/07_strict_signal_grouping`: colocalized-signal grouping.
8. `modules/08_mechanistic_gene_effects`: PIP-weighted gene effects and
   regulatory annotations.
9. `modules/09_sldsc_heritability`: pcQTL/eQTL S-LDSC analysis.
10. `modules/10_representative_loci`: GIMAP and supplementary locus analyses.
11. `modules/11_manuscript_outputs`: final figure and supplementary-table
    generation.
12. `modules/12_joint_score_clustering_sensitivity`: joint-score clustering
    sensitivity analysis.

Run modules 01-04 and 06-11 in numerical order. Module 05 is independent, and
module 12 is a parallel sensitivity analysis. Each module README lists its
inputs, commands, and outputs.

## Setup

The workflow requires R, Python, `fasthurdle`, SAIGE-QTL, PLINK, SMR,
TensorQTL, LDSC, and liftOver. Package and input requirements are listed in
`docs/dependencies.md` and `docs/input_inventory.md`.

```bash
cp config/example.env config/local.env
# Edit paths in config/local.env, then:
source config/local.env
pip install -r requirements.txt
```

Primary analyses use `config/celltype_eligibility.tsv`, which retains OneK1K
cell types with at least 10,000 cells. External OneK1K, FinnGen, annotation,
and LD-reference files must be obtained from their original providers. The
workflow is designed for Slurm-based Linux systems, although the underlying R
and Python scripts can be run directly after inputs are prepared.

## Citation

Citation metadata for sc-pcQTL are provided in `CITATION.cff`.

## License

Source code is released under the MIT License. External software and data
remain subject to their original licenses and access terms.
