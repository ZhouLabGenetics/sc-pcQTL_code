# Module 04 — Cluster annotation enrichment (neighboring null)

Tests whether local sc-pcQTL clusters are enriched for shared genomic
annotations relative to length- and size-matched neighboring null gene sets.
Donor-level Spearman signs are computed only to classify retained sc-pcQTL
clusters as positive-only, negative-only, or mixed; pseudobulk clusters are not
constructed or tested for enrichment.

- **Environment variables:** `SC_PCQTL_CLUSTER_ENRICHMENT_ROOT`, `SC_PCQTL_UPSTREAM_ROOT`, `SC_PCQTL_UPSTREAM_DATA_ROOT`, optional `SC_PCQTL_CLUSTER_ENRICHMENT_FIGURE_OUTPUT_DIR`.
- **Upstream inputs:** cluster definitions from module 01, cell-level counts and
  donor labels used to calculate sign strata, and downloaded annotation
  resources (GO, ABC enhancer links, CTCF, TAD boundaries, paralogs, gene
  models).

## Scripts (run order)

Numbered `00`–`13`:

| Stage | Scripts |
|---|---|
| Stage inputs / collect clusters | `00_stage_inputs.sh`, `01_collect_add_cov_clusters.R` |
| Correlation-type signs | `10_stage_pb_slim_input.sh`, `02_pb_spearman_celltype.R` |
| Build matched neighboring null sets | `04_build_null_sets.R` |
| Annotations | `05_download_annotations.sh`, `06_prepare_annotations.R`, `07_annotate_cluster_sets.R` |
| Enrichment + final plot | `08_run_enrichment.R`, `12_plot_publication_all_strata.R` |
| Enrichment post-processing (Slurm) | `13_submit_pcqtl_style_enrichment.slurm` |

`run_enrichment_after_pair_signs.sh` runs the enrichment stages after the
per-cell-type sign files have been generated.

## Running

Run `00_stage_inputs.sh` and `01_collect_add_cov_clusters.R`, then run
`10_stage_pb_slim_input.sh` and `02_pb_spearman_celltype.R` once per primary
cell type. Finally, run `run_enrichment_after_pair_signs.sh`. Annotation
downloads require network access to the resource providers. Input staging reads
the central cell-type eligibility manifest and copies only the 10 primary cell
types.

## Outputs

`enrichment_by_method.tsv` and the final figure
`cluster_annotation_enrichment_add_cov.pdf` (plus its PNG preview; written to
`SC_PCQTL_CLUSTER_ENRICHMENT_FIGURE_OUTPUT_DIR` when set) and the supporting
supplementary table. Tests with a minimum expected contingency-
table cell count below 1 are skipped. Script
`12_plot_publication_all_strata.R` is the sole figure-reproduction entry point.
