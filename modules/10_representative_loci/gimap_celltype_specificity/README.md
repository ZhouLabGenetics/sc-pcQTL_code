# GIMAP Cell-Type Specificity Analysis

This workflow compares the `CD8_NC` GIMAP pcQTL across OneK1K cell types.

## Run

Source `config/example.env`, then run the wrapper from this code checkout:

```bash
bash "${SC_PCQTL_CODE_ROOT}/modules/10_representative_loci/gimap_celltype_specificity/scripts/refresh_gimap_specificity_analysis.sh"
```

The wrapper writes tables and logs beneath `SC_PCQTL_GIMAP_SPECIFICITY_ROOT`.

## Outputs

- `data/sc_gimap_best_pcqtl_by_celltype.tsv`
- `data/sc_gimap_module_completeness_by_celltype.tsv`

## Inputs

Single-cell summaries are read from:

- `${SC_PCQTL_UPSTREAM_CELLTYPES_ROOT}/<celltype>/cluster_identification/results/method2_sc_hurdle/merged_clusters/cluster_gene_assignments.tsv` for the largest GIMAP-family cluster in the 10 primary cell types meeting the 10,000-cell eligibility threshold, including zero categories
- `${SC_PCQTL_PCQTL_SUMMARY_ROOT}/data/all_pcqtl_results.tsv`

Set the paths in `config/example.env` (`SC_PCQTL_UPSTREAM_CELLTYPES_ROOT`,
`SC_PCQTL_PCQTL_SUMMARY_ROOT`, and `SC_PCQTL_GIMAP_SPECIFICITY_ROOT`) and
source it before running. Runtime
outputs are written outside the code checkout.
