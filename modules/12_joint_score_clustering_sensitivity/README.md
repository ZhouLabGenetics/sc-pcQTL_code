# Joint-score clustering sensitivity

This module is an isolated alternative to the primary component-union pair
screen in module 01. It does not modify or overwrite the primary workflow.
For each ordered response-predictor direction, it combines the one-degree-of-
freedom count and detection score statistics as

```text
T_joint = T_count + T_detection ~ chi-square(2)
```

If `M` is the number of all unordered within-chromosome pairs among the
1%-filtered autosomal genes, each ordered-direction joint test uses the Bonferroni
threshold `0.05/(2M)`. An unordered pair is retained when either direction
passes. The downstream greedy cluster caller is unchanged: filtered genes are
ordered by chromosome and position, windows are evaluated from 50 to 2 genes,
and a window is accepted when at least 70% of its unordered pairs are
significant.

Only pairs separated by at most 49 positions in filtered-gene order are fitted.
This is exact for the maximum 50-gene clustering window; the multiplicity
denominator remains the complete `2M` directional-test universe.

## Separation from the primary workflow

Use a dedicated `SC_PCQTL_JOINT_SCORE_WORK_DIR` for each cell type. The module
writes the same downstream interface as module 01 beneath
`results/method2_sc_hurdle`, but never reads or replaces primary association or
cluster outputs. The primary `filtered_genes.tsv` is required as an immutable
input so that the sensitivity analysis changes the pair test rather than the
gene set.

## Requirements

- R with `data.table`.
- `fasthurdle` 1.2.0 with the count and detection score
  test interfaces.
- The same count matrix, gene coordinates, covariates, and library-size values
  used by module 01.
- The central `config/celltype_eligibility.tsv` manifest.

## Per-cell-type run order

1. Copy `config.example.env`, set all required paths, and source it.
2. Run `01_stage_inputs.R`. This validates and stages the primary filtered-gene
   set, covariates, library size, and chromosome count matrices.
3. Run `02_plan_chunks.R [responses_per_task]`.
4. Run `03_run_joint_score_chunk.R <task_id>` for every task in
   `joint_score_task_manifest.tsv`.
5. Run `04_merge_joint_score_chunks.R`.
6. Run `05_call_clusters_chr.R <chromosome>` for chromosomes 1-22.
7. Run `06_merge_clusters.R`.

After steps 1-3, `submit_joint_score_array.sh` submits steps 4-7 as a Slurm
array plus a dependent merge/cluster job. Resource settings are controlled by
the `SC_PCQTL_SLURM_*` environment variables documented in
`config.example.env`.

## Primary-versus-sensitivity comparison

After all 10 eligible cell types have completed, run:

```bash
Rscript 07_compare_primary_joint_score.R \
  /path/to/primary/celltypes \
  /path/to/joint_score/celltypes \
  /path/to/comparison_output
```

Each root must contain one directory per cell type with the standard
`cluster_identification/results/method2_sc_hurdle` layout. The script writes:

- `table_s2_joint_score_cluster_sensitivity.tsv`, the per-cell-type and
  aggregate pair/cluster comparison displayed as Supplementary Table S2;
- `representative_cluster_retention.tsv`, an exact-gene-set retention summary
  for the four manuscript loci.

The comparison is limited to the pairwise-screen and cluster levels; downstream
QTL mapping and colocalization are not rerun in this sensitivity branch.
