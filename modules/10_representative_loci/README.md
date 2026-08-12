# Module 10 — Final representative loci and GIMAP follow-up

This module produces the main GIMAP locus source data and the three retained
supplementary locus examples. All regional association panels use hg38,
the full OneK1K cis-QTL window, and one 1000 Genomes Phase 3 EUR display-LD
reference anchored to the pcQTL lead. FinnGen in-sample LD remains confined to
its official fine-mapping results and is not reused as display LD.

## Final locus set

- `CD8_ET` ASCL2–C11orf21–TSPAN32 PC2 / hematocrit.
- `CD8_ET` EIF3K–ACTN4 PC2 / platelet count.
- `B_IN` EVI2B–EVI2A PC2 / type 2 diabetes.

Each supplementary block contains aligned FinnGen GWAS, pcQTL, GENCODE gene
track, and every same-cluster gene eQTL panel, followed by coloc and
PIP-weighted loading/effect summaries.

## Scripts

| Script | Purpose |
|---|---|
| `00_prepare_gencode_gene_coords.R` | Prepare GENCODE v38 gene and exon tracks. |
| `01_plot_finngen_style_gimap.R` | Prepare the main GIMAP regional plot-data tables. |
| `02_build_publication_locus_manifest.R` | Build the high-confidence locus manifest. |
| `05_plot_single_panels_one_locus.R` | Export complete source tables for one retained locus. |
| `06_compose_supplement_locus_blocks_with_eqtls.R` | Render the three final all-eQTL supplementary blocks. |
| `eur_ld_display_helpers.R` | Compute unified 1000G EUR display LD. |
| `run_final_locus_figures.sh` | Prepare the main GIMAP data and render the retained supplementary set end to end. |

Set `SC_PCQTL_FORMAL_COLOC_ROOT`, `SC_PCQTL_LOCUSZOOM_WORK_ROOT`,
`EUR_PLINK_DIR`, and the R library variables described in the top-level input
inventory. `SC_PCQTL_LOCUS_OUTPUT_DIR` controls where the three final blocks are
written. `00_prepare_gencode_gene_coords.R` writes both the gene-coordinate and
collapsed protein-coding exon tracks beneath
`SC_PCQTL_LOCUSZOOM_WORK_ROOT/data/gencode`. The final main-text GIMAP composite
is assembled only by module 11.

## GIMAP specificity

`gimap_celltype_specificity/` reproduces the cross-cell-type summaries and the
donor-level pseudobulk TensorQTL mixing analysis.
