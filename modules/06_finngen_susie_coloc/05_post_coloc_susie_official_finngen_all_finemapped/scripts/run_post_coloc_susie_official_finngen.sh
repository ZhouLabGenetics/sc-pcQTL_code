#!/usr/bin/env bash
set -euo pipefail

CODE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ROOT_DIR="${SC_PCQTL_FORMAL_COLOC_ROOT:?Set SC_PCQTL_FORMAL_COLOC_ROOT to the external formal-colocalization work directory.}"
mkdir -p "${ROOT_DIR}/05_post_coloc_susie_official_finngen_all_finemapped"

bash "${CODE_ROOT}/05_post_coloc_susie_official_finngen_all_finemapped/scripts/00_stage_upstream_inputs.sh"
Rscript "${CODE_ROOT}/05_post_coloc_susie_official_finngen_all_finemapped/scripts/01_group_susie_signals.R"
Rscript "${CODE_ROOT}/05_post_coloc_susie_official_finngen_all_finemapped/scripts/02_compute_pip_weighted_nominal_effects.R"
