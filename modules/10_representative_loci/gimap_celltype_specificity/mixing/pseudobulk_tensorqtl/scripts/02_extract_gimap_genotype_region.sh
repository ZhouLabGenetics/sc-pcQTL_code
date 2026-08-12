#!/usr/bin/env bash
set -euo pipefail
BASE="${SC_PCQTL_GIMAP_MIXING_ROOT:?Set SC_PCQTL_GIMAP_MIXING_ROOT.}"
PLINK="${PLINK:-$(command -v plink || echo plink)}"
SRC_TEMPLATE="${ONEK1K_MAF005_GENOTYPE_PREFIX:?Set ONEK1K_MAF005_GENOTYPE_PREFIX.}"
if [[ "${SRC_TEMPLATE}" == *'%s'* ]]; then
  printf -v SRC "${SRC_TEMPLATE}" 7
else
  SRC="${SRC_TEMPLATE}"
fi
OUT="${BASE}/genotypes/gimap_chr7_union_1mb_maf005"
LD_OUT="${BASE}/coloc_susie/gimap_union"
mkdir -p "${BASE}/genotypes" "${BASE}/coloc_susie" "${BASE}/logs"
for ext in bed bim fam; do
  [[ -s "${SRC}.${ext}" ]] || { echo "ERROR: missing genotype input ${SRC}.${ext}" >&2; exit 1; }
done
"${PLINK}" --bfile "${SRC}" --chr 7 --from-bp 149450629 --to-bp 151737348 --make-bed --out "${OUT}" --silent
printf '%s\n' "${OUT}" > "${BASE}/genotypes/gimap_tensorqtl_plink_prefix.txt"
"${PLINK}" --bfile "${OUT}" --r square gz --out "${LD_OUT}" --silent
awk 'BEGIN {OFS="\t"; print "chr","variant_id","cm","pos","other_allele","effect_allele"} {print $1,$2,$3,$4,$5,$6}' \
  "${OUT}.bim" > "${BASE}/coloc_susie/gimap_union.ld_variants.tsv"
echo "Extracted genotype prefix: ${OUT}"
echo "Wrote OneK1K donor LD: ${LD_OUT}.ld.gz"
