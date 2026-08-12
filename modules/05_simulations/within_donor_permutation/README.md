# Within-donor permutation diagnostic

This sensitivity workflow reproduces Supplementary Figure S2. It is parallel
to the primary real-data global gene-wise permutation null and does not
replace that analysis.

## Analysis contract

The global and within-donor analyses use the same 200 genes, 124,968 cells from
982 donors, 1% nonzero-cell eligibility cutoff, unadjusted Poisson count and
binomial detection components, two ordered response--predictor directions as
separate test records, and donor-mean pseudobulk Spearman comparator. The only
intended change is the permutation unit:

- the primary null independently permutes each gene across all cells;
- this diagnostic independently permutes each gene among cells from the same
  donor.

Within-donor permutation preserves every donor-gene sum and mean exactly. It
therefore retains donor-level co-expression and is not a strict null for the
pseudobulk comparator. Its threshold-crossing proportions are reported as
apparent positive fractions, not empirical false-positive rates.

## Inputs

Source `config/example.env` from the repository root. The workflow consumes:

- `COQTL_RAW_COUNTS_FILE`, the unpermuted OneK1K count matrix;
- `SC_PCQTL_R_LIBS`, an R library containing the published `fasthurdle` 1.1.1
  and `fastglm` 0.0.4 versions;
- `$COQTL_SIM_DATA_ROOT/shuffle_full/full_realcounts/sc_counts_shuffle_null.rds`,
  which fixes the primary real-data permutation genes, cells, and donors; and
- `$COQTL_SIM_NULL_ROOT/results/full_realcounts`, the primary global-permutation
  directional p-value output.

Runtime data and results are written beneath `COQTL_WITHIN_DONOR_ROOT`, outside
the source checkout.

## Run

```bash
source config/local.env
modules/05_simulations/within_donor_permutation/submit_pipeline.sh
```

The submitted job requests 12 CPUs and runs 12 one-core forked workers. The
published run is pinned to `fasthurdle` 1.1.1 and `fastglm` 0.0.4, matching the
primary real-data permutation null. The workflow stops on a version mismatch and
writes the loaded versions to `analysis_versions.tsv`.

Key outputs are:

- `threshold_fraction_comparison_global_vs_within.tsv`;
- `qq_panel_summary.tsv`;
- `permutation_validation.tsv`;
- `plots/qq_global_vs_within_donor.{pdf,png,tiff}`.
