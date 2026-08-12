#!/usr/bin/env bash
set -euo pipefail

ROOT="${SC_PCQTL_FORMAL_COLOC_ROOT:?Set SC_PCQTL_FORMAL_COLOC_ROOT to the external formal-colocalization work directory.}"
TASKS="${1:-$ROOT/manifests/qtl_ld_tasks.tsv}"
PLINK="${PLINK:-plink}"
TASK_INDEX="${2:-${SLURM_ARRAY_TASK_ID:-}}"

if [[ ! -s "$TASKS" ]]; then
  echo "Missing LD task manifest: $TASKS" >&2
  exit 1
fi

if [[ -n "$TASK_INDEX" ]]; then
  TASK_STREAM=$(mktemp /tmp/formal_coloc_ld_task.XXXXXX.tsv)
  { head -n 1 "$TASKS"; awk -v n="$TASK_INDEX" 'NR==n+1 {print}' "$TASKS"; } > "$TASK_STREAM"
else
  TASK_STREAM="$TASKS"
fi

tail -n +2 "$TASK_STREAM" | while IFS=$'\t' read -r celltype cluster_id chr start end region_file genotype_prefix ld_prefix ld_matrix ld_variants; do
  [[ -z "$celltype" || "$celltype" == "NA" ]] && continue
  mkdir -p "$(dirname "$ld_prefix")"
  if [[ -s "$ld_matrix" && -s "$ld_variants" ]]; then
    echo "exists $celltype $cluster_id"
    continue
  fi

  tmp="${ld_prefix}.filtered"
  echo "LD $celltype $cluster_id chr$chr:$start-$end"
  "$PLINK" \
    --bfile "$genotype_prefix" \
    --chr "$chr" --from-bp "$start" --to-bp "$end" \
    --maf 0.05 \
    --make-bed \
    --out "$tmp"

  awk 'BEGIN{OFS="\t"; print "chr","variant_id","cm","pos","other_allele","effect_allele"} {print $1,$2,$3,$4,$5,$6}' \
    "${tmp}.bim" > "$ld_variants"

  "$PLINK" \
    --bfile "$tmp" \
    --threads 1 \
    --r square gz \
    --out "$ld_prefix"

  if [[ "${ld_prefix}.ld.gz" != "$ld_matrix" ]]; then
    mv "${ld_prefix}.ld.gz" "$ld_matrix"
  fi
  rm -f "${tmp}.bed" "${tmp}.bim" "${tmp}.fam" "${tmp}.log" "${tmp}.nosex" "${ld_prefix}.log" "${ld_prefix}.nosex"
done

if [[ -n "$TASK_INDEX" ]]; then
  rm -f "$TASK_STREAM"
fi
