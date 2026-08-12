# Input Inventory

The release is source code only. Large, restricted, or provider-governed data are not redistributed.

Required input classes:

- OneK1K per-cell or per-cell-type SCTransform-corrected count tables with donor covariates and library-size columns.
- OneK1K gene coordinate table in the coordinate convention used for cluster construction.
- OneK1K donor genotypes in PLINK format for SAIGE-QTL and QTL fine-mapping LD.
- OneK1K original single-gene SAIGE-QTL tarballs for eQTL fine-mapping,
  S-LDSC annotation construction, locus panels, and targeted SMR follow-up.
- FinnGen R12 SuSiE fine-mapping files and trait metadata for colocalization;
  raw summary statistics, official SNP-heritability estimates, and the Finnish
  genetic-correlation matrix for the 247 prespecified S-LDSC traits.
- Genomic annotation resources for GO, ABC enhancer links, CTCF peaks, TAD boundaries, paralogs, and gene models.
- LDSC baselineLD v2.2, 1000 Genomes Phase 3 EUR LD scores/weights/frequencies, and HapMap3 regression SNPs excluding the MHC.
- Outputs from modules 01-10 for manuscript tables and display panels.
- Harmonized GIMAP PC3 and FinnGen R12 endpoint 3019198 SMR `.ma` inputs, the seven-gene CD8-naive-T eQTL BESD, and OneK1K chromosome 7 LD.
- Staged per-cell raw-count/covariate tables for the seven GIMAP genes in the
  10 primary cell types, used only for the donor-level pseudobulk mixing analysis.

The main environment variables are listed in `config/example.env`.
