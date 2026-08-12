# Real-data permutation-based null calibration

This workflow builds the `fpr_comparison.tsv` source table used by the
real-data permutation calibration panel of Figure 2. Code is read from this directory;
runtime data and results are written to `COQTL_SIM_DATA_ROOT` and
`COQTL_SIM_NULL_ROOT`.

## Current implementation

1. `generate_shuffle_null_full.R` samples 200 genes and preserves the selected
   OneK1K marginal count distributions while independently permuting each gene
   across cells.
2. `test_sc_hurdle_null.R` fits unadjusted Poisson/binomial hurdle models in
   both ordered response--predictor directions. The two directions remain
   separate test records; count--detection component-union status is evaluated within each
   direction.
3. `test_pb_spearman_null.R` recomputes donor means from the permuted matrix and
   applies Spearman tests.
4. `compare_fpr_results.R` applies the 1% nonzero filter and writes the FPR
   source table at four nominal thresholds.

No age, sex, genotype-PC, PEER-factor, or library-size term is fitted in this
permutation simulation because independent gene-wise permutation has already
broken expression alignment with the measured covariates. Pair-mean matrices
are retained only as diagnostic artifacts; the final comparator requires the
directional p-value files and will not fall back to pair means.

## Run

After sourcing `config/example.env`:

```bash
./slurm/submit_full_pipeline.sh full_realcounts
```

The final source file is
`$COQTL_SIM_NULL_ROOT/results/full_realcounts/fpr_comparison.tsv`. Plot styling
is implemented only in module 11.
