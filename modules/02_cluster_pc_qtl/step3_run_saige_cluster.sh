#!/bin/bash
# ============================================================================
# Run SAIGE-QTL (step1 → step2 → step3) for ONE cluster
# Invoked via SLURM array: SLURM_ARRAY_TASK_ID == 1-based cluster index
#   in the de-duplicated cluster list extracted from cluster_pc_map.tsv
# ============================================================================
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
PCQTL_DIR="${SC_PCQTL_PCQTL_DIR:?Set SC_PCQTL_PCQTL_DIR before running SAIGE-QTL.}"

CLUSTER_IDX=${1:?Usage: $0 <cluster_index>}

# ── paths ───────────────────────────────────────────────────────────────────
SIF="${SAIGEQTL_SIF:?Set SAIGEQTL_SIF to the SAIGE-QTL Singularity image.}"
CONTAINER_RUNTIME="${SAIGEQTL_RUNTIME:-}"
if [[ -z "${CONTAINER_RUNTIME}" ]]; then
    CONTAINER_RUNTIME="$(command -v apptainer || command -v singularity || true)"
elif [[ "${CONTAINER_RUNTIME}" != */* ]]; then
    CONTAINER_RUNTIME="$(command -v "${CONTAINER_RUNTIME}" || true)"
fi
if [[ -z "${CONTAINER_RUNTIME}" || ! -x "${CONTAINER_RUNTIME}" ]]; then
    echo "ERROR: set SAIGEQTL_RUNTIME to an executable apptainer or singularity binary." >&2
    exit 1
fi
SAIGE_DIR="${PCQTL_DIR}/step3_saige"
STEP2_DIR="${PCQTL_DIR}/step2_pca"
GENO_DIR="${ONEK1K_RAW_GENOTYPE_DIR:?Set ONEK1K_RAW_GENOTYPE_DIR to the PLINK genotype directory.}"
VR_PLINK="${ONEK1K_VARIANCE_RATIO_PLINK_PREFIX:?Set ONEK1K_VARIANCE_RATIO_PLINK_PREFIX to the pruned PLINK prefix for variance-ratio estimation.}"
PC_MAP="${SAIGE_DIR}/cluster_pc_map.tsv"
SAIGE_DATA_BIND="${SAIGE_DATA_BIND:-$(dirname "$(dirname "${GENO_DIR}")"):$(dirname "$(dirname "${GENO_DIR}")")}"

# ── resolve cluster_id from index ───────────────────────────────────────────
CLUSTER_ID=$(awk -F'\t' 'NR>1{print $1}' "${PC_MAP}" | awk '!seen[$0]++' | sed -n "${CLUSTER_IDX}p")
if [[ -z "${CLUSTER_ID:-}" ]]; then
    echo "No cluster at index ${CLUSTER_IDX} — nothing to do."
    exit 0
fi

# ── cluster metadata from pc_map ────────────────────────────────────────────
CHR=$(awk -F'\t' -v cid="${CLUSTER_ID}" 'NR>1 && $1==cid {print $4; exit}' "${PC_MAP}")
# PC names (one per line)
PCS=$(awk -F'\t' -v cid="${CLUSTER_ID}" 'NR>1 && $1==cid {print $2}' "${PC_MAP}")

PHENO_FILE="${STEP2_DIR}/${CLUSTER_ID}/pheno_with_pcs.tsv"
REGION_FILE="${SAIGE_DIR}/regions/${CLUSTER_ID}_region.txt"

echo "============================================================"
echo " SAIGE-QTL: ${CLUSTER_ID}  chr${CHR}"
echo " PCs      : $(echo ${PCS} | tr '\n' ' ')"
echo " Start    : $(date)"
echo "============================================================"

if [[ ! -f "${PHENO_FILE}" ]]; then
    echo "ERROR: Phenotype file not found: ${PHENO_FILE}"
    exit 1
fi
if [[ ! -f "${REGION_FILE}" ]]; then
    echo "ERROR: Region file not found: ${REGION_FILE}"
    exit 1
fi

# ── output directories ──────────────────────────────────────────────────────
mkdir -p "${SAIGE_DIR}/step1/${CLUSTER_ID}"
mkdir -p "${SAIGE_DIR}/step2/${CLUSTER_ID}"
mkdir -p "${SAIGE_DIR}/step3/${CLUSTER_ID}"

# ── iterate over PCs ────────────────────────────────────────────────────────
while IFS= read -r PC; do
    [[ -z "${PC}" ]] && continue

    STEP1_PFX="${SAIGE_DIR}/step1/${CLUSTER_ID}/${PC}"
    STEP2_OUT="${SAIGE_DIR}/step2/${CLUSTER_ID}/${PC}"
    STEP3_OUT="${SAIGE_DIR}/step3/${CLUSTER_ID}/${PC}_genePval"

    # ─── Step 1: fitNULLGLMM ─────────────────────────────────────────────
    echo "--- Step 1: ${PC} ---"
    if "${CONTAINER_RUNTIME}" exec \
        --bind "${SCRIPT_DIR}:${SCRIPT_DIR}" \
        --bind "${PCQTL_DIR}:${PCQTL_DIR}" \
        --bind "${SAIGE_DATA_BIND}" \
        --cleanenv "${SIF}" \
        step1_fitNULLGLMM_qtl.R \
        --useSparseGRMtoFitNULL=FALSE \
        --useGRMtoFitNULL=FALSE \
        --phenoFile="${PHENO_FILE}" \
        --phenoCol="${PC}" \
        --covarColList=age,sex,pc1,pc2,pc3,pc4,pc5,pc6,pf1,pf2 \
        --sampleCovarColList=age,sex,pc1,pc2,pc3,pc4,pc5,pc6,pf1,pf2 \
        --sampleIDColinphenoFile=individual \
        --traitType=quantitative \
        --outputPrefix="${STEP1_PFX}" \
        --invNormalize=TRUE \
        --skipVarianceRatioEstimation=FALSE \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE \
        --isCovariateTransform=TRUE \
        --skipModelFitting=FALSE \
        --tol=0.00001 \
        --plinkFile="${VR_PLINK}" \
        --IsOverwriteVarianceRatioFile=TRUE
    then
        echo "  Step 1 OK"
    else
        echo "  ERROR: Step 1 failed for ${PC} — skipping step 2 + 3"
        continue
    fi

    # verify outputs
    if [[ ! -f "${STEP1_PFX}.rda" ]] || [[ ! -f "${STEP1_PFX}.varianceRatio.txt" ]]; then
        echo "  ERROR: Step 1 outputs missing for ${PC}"
        continue
    fi

    # ─── Step 2: SPAGMMATtest ────────────────────────────────────────────
    echo "--- Step 2: ${PC} ---"
    if "${CONTAINER_RUNTIME}" exec \
        --bind "${SCRIPT_DIR}:${SCRIPT_DIR}" \
        --bind "${PCQTL_DIR}:${PCQTL_DIR}" \
        --bind "${SAIGE_DATA_BIND}" \
        --cleanenv "${SIF}" \
        step2_tests_qtl.R \
        --bedFile="${GENO_DIR}/full_genome_chr${CHR}.bed" \
        --bimFile="${GENO_DIR}/full_genome_chr${CHR}.bim" \
        --famFile="${GENO_DIR}/full_genome_chr${CHR}.fam" \
        --SAIGEOutputFile="${STEP2_OUT}" \
        --chrom="${CHR}" \
        --minMAF=0.05 \
        --LOCO=FALSE \
        --GMMATmodelFile="${STEP1_PFX}.rda" \
        --varianceRatioFile="${STEP1_PFX}.varianceRatio.txt" \
        --rangestoIncludeFile="${REGION_FILE}" \
        --markers_per_chunk=10000
    then
        echo "  Step 2 OK"
    else
        echo "  ERROR: Step 2 failed for ${PC} — skipping step 3"
        continue
    fi

    # ─── Step 3: ACAT gene-level p-value ─────────────────────────────────
    echo "--- Step 3: ${PC} ---"
    if "${CONTAINER_RUNTIME}" exec \
        --bind "${SCRIPT_DIR}:${SCRIPT_DIR}" \
        --bind "${PCQTL_DIR}:${PCQTL_DIR}" \
        --bind "${SAIGE_DATA_BIND}" \
        --cleanenv "${SIF}" \
        step3_gene_pvalue_qtl.R \
        --assocFile="${STEP2_OUT}" \
        --geneName="${PC}" \
        --genePval_outputFile="${STEP3_OUT}"
    then
        echo "  Step 3 OK"
    else
        echo "  WARNING: Step 3 (ACAT) failed for ${PC} — may have no valid p-values"
    fi

done <<< "${PCS}"

echo "============================================================"
echo " SAIGE-QTL ${CLUSTER_ID}: DONE  $(date)"
echo "============================================================"
