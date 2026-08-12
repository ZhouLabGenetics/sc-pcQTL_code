# GIMAP Pseudobulk TensorQTL Analysis

This workflow constructs donor-level pseudobulk mixtures from OneK1K single-cell data and
analyzes GIMAP single-gene and cluster-PC phenotypes with TensorQTL, matching the
bulk pcQTL logic more closely than the single-cell SAIGE-QTL workflow.
The TensorQTL models use donor-level age, sex, and genotype principal components
as covariates. Cell-type-specific OneK1K PEER factors are deliberately excluded
because they are not defined after synthetic cell-type mixing.

The analysis uses two phenotype sets (`single_gene_only` and `pc_only`) and the
following scenarios: `CD8_NC` only, observed
cell-count weighting, equal cell-type weighting, and a 0%, 5%, 10%, 25%, 50%,
and 75% CD8_NC gradient. `cd8_nc_only` is the 100% endpoint.

## Prepare inputs

Source `config/example.env`. In particular, set
`SC_PCQTL_GIMAP_MIXING_ROOT`, `SC_PCQTL_GIMAP_SELECTED_READCOUNT_DIR`,
`ONEK1K_MAF005_GENOTYPE_PREFIX`, and `SC_PCQTL_FORMAL_COLOC_ROOT`.

```bash
bash "${SC_PCQTL_CODE_ROOT}/modules/10_representative_loci/gimap_celltype_specificity/mixing/pseudobulk_tensorqtl/scripts/refresh_prepare.sh"
```

## TensorQTL environment

The default system Python currently lacks TensorQTL. Create the local environment:

```bash
bash "${SC_PCQTL_CODE_ROOT}/modules/10_representative_loci/gimap_celltype_specificity/mixing/pseudobulk_tensorqtl/scripts/00_setup_tensorqtl_env.sh"
```

Alternatively set `TENSORQTL_PYTHON` to an existing Python with TensorQTL installed.

## Run TensorQTL

Run one scenario:

```bash
bash "${SC_PCQTL_CODE_ROOT}/modules/10_representative_loci/gimap_celltype_specificity/mixing/pseudobulk_tensorqtl/scripts/03_run_tensorqtl_job.sh" cd8_nc_only both
```

Run all jobs:

```bash
bash "${SC_PCQTL_CODE_ROOT}/modules/10_representative_loci/gimap_celltype_specificity/mixing/pseudobulk_tensorqtl/scripts/03_run_tensorqtl_job.sh" all both
```

After TensorQTL completes, export the full cis-nominal input, run the GIMAP
SuSiE sensitivity, and build the Figure S11c source table:

```bash
"${TENSORQTL_PYTHON}" "${SC_PCQTL_CODE_ROOT}/modules/10_representative_loci/gimap_celltype_specificity/mixing/pseudobulk_tensorqtl/scripts/04_summarize_tensorqtl_outputs.py"
sbatch --output="${SC_PCQTL_GIMAP_MIXING_ROOT}/logs/gimap_coloc_susie_%j.out" \
  --error="${SC_PCQTL_GIMAP_MIXING_ROOT}/logs/gimap_coloc_susie_%j.err" \
  "${SC_PCQTL_CODE_ROOT}/modules/10_representative_loci/gimap_celltype_specificity/mixing/pseudobulk_tensorqtl/scripts/06_submit_coloc_susie.sbatch"
Rscript "${SC_PCQTL_CODE_ROOT}/modules/10_representative_loci/gimap_celltype_specificity/mixing/pseudobulk_tensorqtl/scripts/07_make_coloc_susie_publication_summary.R"
```

## Outputs

- `data/gimap_pseudobulk_tensorqtl_phenotype_level_results.tsv`
- `coloc_susie/mix_all_sumstats.tsv`
- `coloc_susie/mix_coloc_susie_summary.tsv`
- `coloc_susie/gimap_mixing_coloc_susie_publication_summary.tsv`

## Interpretation boundary

This is a sensitivity analysis of how donor-level cell-type mixing changes
colocalization with FinnGen lymphocyte count.

QTL fine-mapping uses the cohort-level OneK1K chr7 donor LD matrix for every
mixture. FinnGen is not re-fine-mapped: colocalization uses its official R12
in-sample SuSiE object. OneK1K QTL coordinates are lifted from hg19 to hg38 and
alleles are harmonized before `coloc.susie`.
