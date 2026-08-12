# Count-model calibration diagnostic

This supplementary diagnostic compares Poisson and Negative-Binomial
positive-count models under the same real-data permutation null. Wald and
score-test settings are shown in the original 3-by-4 QQ grid used as
Supplementary Figure S3. The count matrix, eligible genes, ordered gene-pair
tests, and binomial detection component are held fixed.

A separate controlled simulation evaluates the scope of this comparison under
three positive-count regimes: zero-truncated Poisson, shifted Poisson, and
overdispersed zero-truncated Negative-Binomial. The shifted-Poisson regime uses
`1 + Poisson(1.0)` for the predictor and `1 + Poisson(0.6)` for the response,
matching the boundary-case diagnostic described in the Supplementary
Information.

The manuscript analysis used and this script requires `fasthurdle` 1.2.0.

The input RDS is produced by
`../real_data_permutation_null/scripts/generate_shuffle_null_full.R` and must
contain a cells-by-genes `counts` matrix and matching `genes` vector. Runtime
data and results remain outside this repository.

```bash
Rscript 01_compare_count_models.R \
  "$COQTL_SIM_DATA_ROOT/shuffle_full/full_realcounts/sc_counts_shuffle_null.rds" \
  "$COQTL_SIM_NULL_ROOT/count_model_calibration"

Rscript 02_plot_count_model_qq.R \
  "$COQTL_SIM_NULL_ROOT/count_model_calibration" \
  "$COQTL_SIM_NULL_ROOT/count_model_calibration/qq_distribution_assumptions_score_test_v12.pdf"

Rscript 03_simulate_count_regimes.R \
  "$COQTL_SIM_NULL_ROOT/count_model_calibration/controlled_regimes" \
  20000 20000 20 20260812

Rscript 04_plot_count_regimes.R \
  "$COQTL_SIM_NULL_ROOT/count_model_calibration/controlled_regimes" \
  "$COQTL_SIM_NULL_ROOT/count_model_calibration/controlled_regimes/count_regime_calibration.pdf"
```

The first script writes directional count- and detection-component p-values
for all four settings and compact rejection-rate summaries. The second script
reproduces the original Supplementary Figure S3 panel. The third and fourth
scripts generate the controlled-regime calibration and dispersion panels,
including threshold-specific rejection rates and Negative-Binomial dispersion
summaries. Runtime p-values and source tables are not stored in this repository.
