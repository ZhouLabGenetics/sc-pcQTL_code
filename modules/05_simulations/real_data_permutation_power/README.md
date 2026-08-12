# Real-data permutation-based signal-injection power

This workflow builds the `power_comparison.tsv` source tables used by Figure
2d. Runtime data are written beneath `COQTL_SIM_DATA_ROOT/power_full`; analysis
outputs are written beneath `COQTL_SIM_POWER_ROOT`.

## Current implementation

- 200 sampled genes.
- Four injected modules of 20 genes.
- Effect strengths `0.02`, `0.04` and `0.06`.
- Per-gene permutation supplies the null matrix.
- Hurdle models are unadjusted after permutation; no measured covariates or
  library-size terms are included.
- Both ordered response--predictor directions are retained as separate tests,
  and count--detection component-union status is evaluated within direction.
- Single-cell hurdle and pseudobulk Spearman are evaluated on the same
  SC-derived eligible-pair mask.
- Output is summarized at nominal thresholds `0.05`, `0.01`, `0.005` and
  `0.001`; module 11 performs all manuscript plotting.

This real-data permutation simulation intentionally differs from the covariate-adjusted
real-data screen: independent gene-wise permutation removes alignment with the
measured covariates, so the simulated hurdle fits are unadjusted.

## Run

After sourcing `config/example.env`:

```bash
./slurm/submit_full_pipeline.sh full_realcounts
```

The runner refuses to mix a new run with an existing data/result directory;
use a new config label or archive the previous run explicitly. Final source
tables are under
`$COQTL_SIM_POWER_ROOT/results/<label>/compare/strength_s*/power_comparison.tsv`.
