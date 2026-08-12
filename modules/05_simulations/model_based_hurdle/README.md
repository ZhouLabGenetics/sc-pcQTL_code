# Model-based hurdle simulation

This is the model-based calibration/power workflow used to build the source
table for Figure 2a,c. Runtime data are written beneath
`SC_PCQTL_MODEL_SIM_WORK_ROOT`; the code checkout remains read-only.

## Final design

- 300 donors, 100 genes and 50 replicates per scenario.
- Five percent of all unordered gene pairs are signal pairs, sampled
  preferentially from genes with reference detection rates and positive-count
  intensities at or above their respective medians.
- Seven scenarios: null plus low/high detection-only, count-only and joint
  signals. Detection-only settings use the internal `zero_only` label.
- Detection is Bernoulli/logit and positive counts are zero-truncated Poisson.
- The count fit is Poisson with `log_total_read_counts` as an offset; the
  detection fit is binomial/logit with library size as an ordinary covariate.
- Simulated expression has no donor-level covariate effect. Randomized donor
  sex is retained as a null adjustment covariate; age, genotype PCs and PEER
  factors are constant and are omitted automatically by the hurdle worker.
- Pseudobulk Pearson and Spearman, computed across donor means, are the
  aggregation baselines.
- Hurdle count and detection directions are retained as separate tests in this
  simulation. The component-union decision is evaluated within each direction;
  no directional minimum is taken for simulation calibration or power.

## Scripts

1. `01_fit_reference_from_real_data.R` estimates empirical cell-count,
   library-size and gene-level parameters.
2. `02_simulate_hurdle_counts.R` generates one scenario/replicate.
3. `03_run_hurdle_association.R` runs the released Poisson hurdle worker.
4. `06_run_donor_aggregate_cor.R` runs pseudobulk Pearson/Spearman.
5. `10_make_target_figures.R` writes only the final source tables:
   `plots_target/alpha_curve_replicate_metrics.tsv` and
   `plots_target/alpha_curve_summary.tsv`. Module 11 draws the manuscript
   figure from the summary table.

## Run

Source `config/example.env`, set `COQTL_RAW_COUNTS_FILE`, and run:

```bash
./submit_large_test_slurm.sh
```

The non-Slurm full runner is `run_large_test.sh`.
