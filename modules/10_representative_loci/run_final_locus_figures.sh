#!/usr/bin/env bash
set -euo pipefail

MODULE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT="${COQTL_RSCRIPT:-${R_SCRIPT:-Rscript}}"
: "${SC_PCQTL_FORMAL_COLOC_ROOT:?Set SC_PCQTL_FORMAL_COLOC_ROOT}"
: "${SC_PCQTL_LOCUSZOOM_WORK_ROOT:=${SC_PCQTL_FORMAL_COLOC_ROOT}/10_publication_locuszoom_redesign}"
export SC_PCQTL_LOCUSZOOM_WORK_ROOT

"${RSCRIPT}" "${MODULE_DIR}/00_prepare_gencode_gene_coords.R"
"${RSCRIPT}" "${MODULE_DIR}/01_plot_finngen_style_gimap.R"
"${RSCRIPT}" "${MODULE_DIR}/02_build_publication_locus_manifest.R"

figure_ids=(
  high_confidence_shared_mass_pass_cd8_et_SC_chr11_cluster_002_PC2_3009542
  high_confidence_shared_mass_pass_cd8_et_SC_chr19_cluster_024_PC2_3007461
  high_confidence_shared_mass_pass_b_in_SC_chr17_cluster_004_PC2_T2D
)
for figure_id in "${figure_ids[@]}"; do
  FIGURE_ID="${figure_id}" "${RSCRIPT}" "${MODULE_DIR}/05_plot_single_panels_one_locus.R"
done

"${RSCRIPT}" "${MODULE_DIR}/06_compose_supplement_locus_blocks_with_eqtls.R"
