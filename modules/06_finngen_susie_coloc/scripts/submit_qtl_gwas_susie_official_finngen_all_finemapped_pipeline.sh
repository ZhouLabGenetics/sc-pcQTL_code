#!/usr/bin/env bash
set -euo pipefail

CODE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ROOT="${SC_PCQTL_FORMAL_COLOC_ROOT:?Set SC_PCQTL_FORMAL_COLOC_ROOT to the external formal-colocalization work directory.}"
mkdir -p "$ROOT"
cd "$ROOT"

DOWNLOAD_CHUNKS="${1:-192}"
IMPORT_CHUNKS="${2:-192}"
COLOC_CHUNKS="${3:-240}"
DOWNLOAD_MAX="${4:-32}"
IMPORT_MAX="${5:-12}"
COLOC_MAX="${6:-48}"

PHENO_MANIFEST="$ROOT/manifests/finngen_coloc_gwas_phenotype_standard_all_sources.tsv"
IMPORT_DIR="$ROOT/results/fine_mapping/gwas_all_finemapped"
GWAS_MANIFEST="$IMPORT_DIR/finngen_official_cs_by_cluster.tsv"
PAIR_MANIFEST="$ROOT/manifests/qtl_gwas_susie_official_finngen_all_finemapped_pairs.tsv"
COLOC_DIR="$ROOT/results/coloc/qtl_gwas_susie_official_finngen_all_finemapped"
POST_DIR="$ROOT/05_post_coloc_susie_official_finngen_all_finemapped"

mkdir -p logs "$IMPORT_DIR/official_import_chunks" "$COLOC_DIR/chunks" "$POST_DIR/logs"

Rscript "$CODE_ROOT/scripts/20_build_finngen_coloc_gwas_phenotype_standard.R"

download_job=$(sbatch --parsable \
  --partition=normal \
  --job-name=fg_susie_all_dl \
  --array=1-${DOWNLOAD_CHUNKS}%${DOWNLOAD_MAX} \
  --time=12:00:00 \
  --mem=4G \
  --cpus-per-task=1 \
  --output=logs/fg_susie_all_dl_%A_%a.out \
  --error=logs/fg_susie_all_dl_%A_%a.err \
  --wrap="Rscript $CODE_ROOT/scripts/21_download_finngen_official_susie_all_finemapped.R --manifest $PHENO_MANIFEST --chunk-index \${SLURM_ARRAY_TASK_ID} --n-chunks ${DOWNLOAD_CHUNKS}")

import_job=$(sbatch --parsable \
  --dependency=afterok:${download_job} \
  --partition=normal \
  --job-name=fg_susie_all_import \
  --array=1-${IMPORT_CHUNKS}%${IMPORT_MAX} \
  --time=18:00:00 \
  --mem=64G \
  --cpus-per-task=1 \
  --output=logs/fg_susie_all_import_%A_%a.out \
  --error=logs/fg_susie_all_import_%A_%a.err \
  --wrap="Rscript $CODE_ROOT/scripts/21_import_finngen_official_susie_all_finemapped.R --phenotype-manifest $PHENO_MANIFEST --out-dir $IMPORT_DIR --chunk-index \${SLURM_ARRAY_TASK_ID} --n-chunks ${IMPORT_CHUNKS}")

import_merge_job=$(sbatch --parsable \
  --dependency=afterok:${import_job} \
  --partition=normal \
  --job-name=fg_susie_all_merge \
  --time=01:00:00 \
  --mem=16G \
  --cpus-per-task=1 \
  --output=logs/fg_susie_all_merge_%j.out \
  --error=logs/fg_susie_all_merge_%j.err \
  --wrap="Rscript $CODE_ROOT/scripts/21_merge_finngen_official_susie_all_finemapped_import_chunks.R --out-dir $IMPORT_DIR")

pair_job=$(sbatch --parsable \
  --dependency=afterok:${import_merge_job} \
  --partition=normal \
  --job-name=formal_susie_all_pairs \
  --time=01:00:00 \
  --mem=16G \
  --cpus-per-task=1 \
  --output=logs/formal_susie_all_pairs_%j.out \
  --error=logs/formal_susie_all_pairs_%j.err \
  --wrap="rm -f $PAIR_MANIFEST && Rscript $CODE_ROOT/scripts/19_coloc_qtl_gwas_susie_official_finngen.R --gwas-manifest $GWAS_MANIFEST --pair-manifest $PAIR_MANIFEST --out-dir $COLOC_DIR --build-pair-manifest-only")

coloc_job=$(sbatch --parsable \
  --dependency=afterok:${pair_job} \
  --partition=normal \
  --job-name=formal_susie_all \
  --array=1-${COLOC_CHUNKS}%${COLOC_MAX} \
  --time=06:00:00 \
  --mem=24G \
  --cpus-per-task=1 \
  --output=logs/formal_susie_all_%A_%a.out \
  --error=logs/formal_susie_all_%A_%a.err \
  --wrap="Rscript $CODE_ROOT/scripts/19_coloc_qtl_gwas_susie_official_finngen.R --gwas-manifest $GWAS_MANIFEST --pair-manifest $PAIR_MANIFEST --out-dir $COLOC_DIR --chunk-index \${SLURM_ARRAY_TASK_ID} --n-chunks ${COLOC_CHUNKS}")

coloc_merge_job=$(sbatch --parsable \
  --dependency=afterok:${coloc_job} \
  --partition=normal \
  --job-name=formal_susie_all_merge \
  --time=02:00:00 \
  --mem=32G \
  --cpus-per-task=1 \
  --output=logs/formal_susie_all_merge_%j.out \
  --error=logs/formal_susie_all_merge_%j.err \
  --wrap="Rscript $CODE_ROOT/scripts/19_merge_qtl_gwas_susie_official_finngen_chunks.R --out-dir $COLOC_DIR")

post_job=$(sbatch --parsable \
  --dependency=afterok:${coloc_merge_job} \
  --partition=normal \
  --job-name=post_susie_all \
  --time=08:00:00 \
  --mem=48G \
  --cpus-per-task=4 \
  --output=$POST_DIR/logs/post_susie_all_%j.out \
  --error=$POST_DIR/logs/post_susie_all_%j.err \
  --wrap="bash $CODE_ROOT/05_post_coloc_susie_official_finngen_all_finemapped/scripts/run_post_coloc_susie_official_finngen.sh")

echo "FinnGen all-finemapped download job: $download_job"
echo "FinnGen all-finemapped import job: $import_job"
echo "FinnGen all-finemapped import merge job: $import_merge_job"
echo "All-finemapped pair manifest job: $pair_job"
echo "All-finemapped SuSiE coloc job: $coloc_job"
echo "All-finemapped SuSiE coloc merge job: $coloc_merge_job"
echo "All-finemapped post-coloc job: $post_job"
echo "Main output directory: $COLOC_DIR"
