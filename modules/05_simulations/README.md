# Module 05 — Simulations (model-based and real-data permutation-based)

Calibration and performance simulations reported in the manuscript. Five
self-contained sub-pipelines; `shared/` holds the common environment helper
(`module_env.sh`) and path resolver (`sim_paths.R`).

- **Environment variables:** `COQTL_RAW_COUNTS_FILE`, `COQTL_GENE_INFO_FILE`, `COQTL_SIM_DATA_ROOT`, `COQTL_SIM_ROOT`, `COQTL_SIM_NULL_ROOT`, `COQTL_SIM_POWER_ROOT`, `COQTL_WITHIN_DONOR_ROOT`, `COQTL_RSCRIPT`, `SC_PCQTL_HURDLE_SCRIPT_DIR`. Each sub-pipeline also reads a `config/paper_defaults.env`.
- **Requires:** the custom R package `fasthurdle` (see `../../docs/dependencies.md`).

## Sub-pipelines

| Subdirectory | Purpose |
|---|---|
| `model_based_hurdle/` | Fit a reference, simulate hurdle counts, run SC/donor tests, and build the Figure 2a,c source table. |
| `real_data_permutation_null/` | Per-gene permutation null and directional hurdle source-table workflow for Figure 2b. |
| `real_data_permutation_power/` | Signal-injection power comparing the single-cell hurdle test against pseudobulk Spearman (`scripts/`, `slurm/`; see its `README.md`). |
| `within_donor_permutation/` | Matched permutation-unit diagnostic for Supplementary Figure S2; it preserves donor-gene means and is not used as the calibration null. |
| `count_model_calibration/` | Poisson versus negative-binomial calibration under the OneK1K permutation null and controlled positive-count regimes. |

## Running

Use `model_based_hurdle/submit_large_test_slurm.sh` for the model-based run and
the `real_data_permutation_*/slurm/submit_full_pipeline.sh` Slurm chains. Source
`config/example.env`; all runtime roots are external to this code checkout.
Run `within_donor_permutation/submit_pipeline.sh` only after the primary
real-data permutation null output is available.

## Outputs

Source TSVs consumed by the module-11 Figure 2 compositor.
