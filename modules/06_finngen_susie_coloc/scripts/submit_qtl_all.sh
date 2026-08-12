#!/usr/bin/env bash
set -euo pipefail

CODE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ROOT="${SC_PCQTL_FORMAL_COLOC_ROOT:?Set SC_PCQTL_FORMAL_COLOC_ROOT to the external formal-colocalization work directory.}"
mkdir -p "$ROOT"
cd "$ROOT"
mkdir -p logs slurm results/fine_mapping

Rscript "$CODE_ROOT/scripts/01_build_qtl_manifests.R"

n_ld=$(($(wc -l < manifests/qtl_ld_tasks.tsv) - 1))
n_qtl=$(awk -F'\t' 'NR>1 && $NF=="ready" {n++} END{print n+0}' manifests/qtl_finemap_tasks.tsv)

ld_job=$(sbatch --parsable \
  --array=1-${n_ld}%80 \
  --job-name=formal_coloc_ld \
  --output=logs/formal_coloc_ld_%A_%a.out \
  --error=logs/formal_coloc_ld_%A_%a.err \
  --time=12:00:00 \
  --mem=12G \
  --cpus-per-task=1 \
  --wrap="bash $CODE_ROOT/scripts/04_make_qtl_ld.sh $ROOT/manifests/qtl_ld_tasks.tsv \${SLURM_ARRAY_TASK_ID}")

qtl_job=$(sbatch --parsable \
  --dependency=afterok:${ld_job} \
  --array=1-${n_qtl}%80 \
  --job-name=formal_coloc_qtl \
  --output=logs/formal_coloc_qtl_%A_%a.out \
  --error=logs/formal_coloc_qtl_%A_%a.err \
  --time=04:00:00 \
  --mem=8G \
  --cpus-per-task=1 \
  --wrap="Rscript $CODE_ROOT/scripts/05_finemap_qtl.R --tasks $ROOT/manifests/qtl_finemap_tasks.tsv --task-index \${SLURM_ARRAY_TASK_ID}")

merge_job=$(sbatch --parsable \
  --dependency=afterok:${qtl_job} \
  --job-name=formal_coloc_qtl_merge \
  --output=logs/formal_coloc_qtl_merge_%A.out \
  --error=logs/formal_coloc_qtl_merge_%A.err \
  --time=04:00:00 \
  --mem=16G \
  --cpus-per-task=1 \
  --wrap="cd $ROOT && awk 'FNR==1 && NR!=1 {next} {print}' results/fine_mapping/qtl_finemap_status_task_*.tsv > results/fine_mapping/qtl_finemap_status.tsv")

echo "QTL LD array job: $ld_job"
echo "QTL fine-map array job: $qtl_job"
echo "QTL merge/coloc job: $merge_job"
