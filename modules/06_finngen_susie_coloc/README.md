# Module 06 - QTL fine-mapping and FinnGen colocalization

This module fine-maps OneK1K pcQTLs and single-gene eQTLs, imports FinnGen R12
SuSiE results, and runs `coloc.susie`. It also contains the observed-variant
`coloc.abf` sensitivity analysis.

## Inputs

- pcQTL summary statistics from module 02;
- OneK1K single-gene eQTL summary statistics;
- OneK1K genotypes for QTL LD;
- FinnGen R12 SuSiE fine-mapping files;
- hg19-to-hg38 liftOver resources.

Set paths in `config/config.R` or through the corresponding environment
variables, including `SC_PCQTL_FORMAL_COLOC_ROOT`, `SC_PCQTL_UPSTREAM_ROOT`,
`FINNGEN_ROOT`, `ONEK1K_EQTL_ROOT`, and the OneK1K genotype prefixes. QTL
fine-mapping uses the donor sample size in the SAIGE-QTL `N` field.

## Primary run order

1. `00_prepare_onk1k_liftover_maps.sh`
2. `01_build_qtl_manifests.R`
3. `04_make_qtl_ld.sh`
4. `05_finemap_qtl.R`
5. `20_build_finngen_coloc_gwas_phenotype_standard.R`
6. `21_download_finngen_official_susie_all_finemapped.R`
7. `21_import_finngen_official_susie_all_finemapped.R`
8. `21_merge_finngen_official_susie_all_finemapped_import_chunks.R`
9. `19_coloc_qtl_gwas_susie_official_finngen.R`
10. `19_merge_qtl_gwas_susie_official_finngen_chunks.R`

The Slurm entry points are `submit_qtl_all.sh` and
`submit_qtl_gwas_susie_official_finngen_all_finemapped_pipeline.sh`.

Post-colocalization scripts under
`05_post_coloc_susie_official_finngen_all_finemapped/scripts/` group signals
and compute PIP-weighted nominal effects.

## coloc.abf sensitivity

`02_build_finngen_endpoint_manifest.R` defines the prespecified 247-trait
universe used by this analysis. The `12_*` scripts use the corresponding raw
FinnGen summary statistics to run `coloc.abf` in 250-kb windows centered on QTL
lead variants, with allele-compatible observed variants, at least 50 shared
SNPs, and an MHC exclusion. `14_plot_abf_susie_sensitivity_comparison.R`
generates the ABF-versus-SuSiE comparison.

## Outputs

The module writes QTL credible sets, QTL-GWAS colocalization results, and ABF
sensitivity summaries beneath `SC_PCQTL_FORMAL_COLOC_ROOT`. These outputs feed
modules 07, 08, and 10.
