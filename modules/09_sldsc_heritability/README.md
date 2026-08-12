# Module 09 - pcQTL/eQTL S-LDSC

This module tests disease-heritability enrichment of pcQTL and single-gene
eQTL annotations across FinnGen R12 traits.

## Analysis settings

- Pool per-phenotype BH-significant cis variants across the 10 eligible cell
  types and restrict both annotations to MAF >= 0.05.
- Match variants to the 1000 Genomes Phase 3 European LDSC panel.
- Use baseline-LD v2.2 after removing its molecular-QTL MaxCPP columns.
- Start from 247 prespecified FinnGen traits, retain traits with heritability
  Z score greater than 4, and prune correlated traits at `|rg| >= 0.7` and
  FDR below 0.05.
- Fit marginal pcQTL/eQTL and joint pcQTL-plus-eQTL models, followed by
  random-effects meta-analysis across the retained traits.

## Inputs

Set `SC_PCQTL_SLDSC_BASE_ROOT`, `SC_PCQTL_SLDSC_ROOT`,
`FINNGEN_SLDSC_ENDPOINT_MANIFEST`, `FINNGEN_R12_RG_MATRIX`,
`COQTL_UPSTREAM_CELLTYPES_DIR`, and `ONEK1K_EQTL_ROOT`. The workflow also
requires LDSC, baseline-LD v2.2, 1000 Genomes European reference files, and
FinnGen summary statistics.

## Run order

1. `74_build_prespecified_finngen_manifest.py`
2. `73_initialize_10cell_sldsc_run.sh`
3. `72_select_independent_finngen_traits.py`
4. `75_download_finngen_sumstats.py`
5. `76_munge_finngen_r12_for_ldsc.py`
6. `PREPARE_INPUTS=1 bash scripts/71_submit_sldsc_model.sh main` for a clean
   first run; omit `PREPARE_INPUTS=1` only when reusing validated inputs.
7. `70_collect_sldsc_models.sh`
8. `68_make_pcqtl_eqtl_heritability_figure.R`

The `functional` and `ngenes` modes of `71_submit_sldsc_model.sh` generate the
secondary annotation-stratified results. Outputs are written beneath
`SC_PCQTL_SLDSC_ROOT`.
