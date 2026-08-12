#!/usr/bin/env bash
set -euo pipefail

ROOT="${SC_PCQTL_FORMAL_COLOC_ROOT:?Set SC_PCQTL_FORMAL_COLOC_ROOT to the external formal-colocalization work directory.}"
SOURCE_DIR="${SOURCE_DIR:-${ONEK1K_SUBSET_GENOTYPE_DIR:-}}"
TARGET_DIR="${TARGET_DIR:-$ROOT/resources/onek1k_liftover}"

mkdir -p "$TARGET_DIR"

for chr in $(seq 1 22); do
  src="$SOURCE_DIR/chr${chr}.posmap"
  dst="$TARGET_DIR/chr${chr}.hg19_to_hg38.posmap"
  if [[ ! -s "$src" ]]; then
    echo "Missing required OneK1K liftOver posmap: $src" >&2
    exit 1
  fi
  cp -f "$src" "$dst"

  unmapped="$SOURCE_DIR/chr${chr}.unmapped.snps"
  if [[ -e "$unmapped" ]]; then
    cp -f "$unmapped" "$TARGET_DIR/chr${chr}.unmapped.snps"
  fi
done

echo "Prepared OneK1K liftOver posmaps under $TARGET_DIR"
