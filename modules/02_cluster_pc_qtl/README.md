# Module 02 — Cluster-PC construction and SAIGE-QTL cis mapping

Turns the local gene clusters from module 01 into principal-component (pcQTL)
phenotypes and maps them in cis. Cluster PCs are computed from SCTransform-corrected
counts (the staged per-cell-type matrices in `COQTL_DATA_ROOT`, or the explicit
`COQTL_COUNT_FILE`); SAIGE-QTL then runs
its null-model → single-variant → region steps for each cluster PC.

- **Environment variables:** `SC_PCQTL_CELL_TYPE`, `SC_PCQTL_UPSTREAM_ROOT` (the pipeline directory containing `celltypes/`), `SC_PCQTL_UPSTREAM_CELLTYPES_ROOT`, `SC_PCQTL_CLUSTER_FILE`, `SC_PCQTL_CLUSTER_GENE_FILE`, `SC_PCQTL_PCQTL_DIR`, `COQTL_DATA_ROOT` (or `COQTL_COUNT_FILE`), `ONEK1K_RAW_GENOTYPE_DIR`, `ONEK1K_VARIANCE_RATIO_PLINK_PREFIX`, `PLINK`, `SAIGEQTL_SIF`, `SAIGEQTL_RUNTIME`.
- **Requires:** SAIGE-QTL (Singularity image) and PLINK (see `../../docs/dependencies.md`).
- **Upstream inputs:** cluster definitions from module 01; OneK1K genotypes and per-cell-type SCTransform-corrected count matrices.

## Scripts (run order)

| Script | Purpose |
|---|---|
| `config.R` | Per-cell-type paths (or set `SC_PCQTL_CELL_TYPE` / `SC_PCQTL_UPSTREAM_CELLTYPES_ROOT`). |
| `step2_cluster_pca.R` | Single-cell PCA per gene cluster → cluster PCs. |
| `step3_prepare_saige_inputs.R` | Region files and cluster-PC phenotype map for SAIGE-QTL. |
| `step3_run_saige_cluster.sh` | Runs SAIGE-QTL (step1 → step2 → step3) for one cluster. |

## Running

Per cell type: `submit_step2.sh` → `submit_step3_prepare.sh <step2_job_id>` →
`submit_step3_saige.sh <prep_job_id>` (each prints its Slurm job id).
`submit_all.sh` chains Step 2 PCA -> Step 3 prepare -> Step 3 SAIGE-QTL. Set
`SC_PCQTL_CELL_TYPE` explicitly for every cell-type run.

## Outputs

Cluster-PC cis-QTL (pcQTL) summary statistics per cell type, consumed by
module 03 (summaries/follow-up) and module 06 (fine-mapping / colocalization).
