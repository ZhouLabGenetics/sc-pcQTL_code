#!/usr/bin/env bash
set -euo pipefail

ROOT="${SC_PCQTL_CLUSTER_ENRICHMENT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
RAW="${ROOT}/annotations/raw"
mkdir -p "${RAW}"
META="${ROOT}/annotations/metadata.tsv"

download_if_missing() {
  local url="$1"
  local out="$2"
  if [[ ! -s "${out}" ]]; then
    echo "[download] ${url}"
    curl -L --retry 3 --fail -o "${out}" "${url}"
  fi
}

echo -e "resource\turl\tlocal_file\tgenome_build\tnote" > "${META}"

GENCODE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
download_if_missing "${GENCODE_URL}" "${RAW}/gencode.v19.annotation.gtf.gz"
echo -e "gencode_v19\t${GENCODE_URL}\t${RAW}/gencode.v19.annotation.gtf.gz\tGRCh37/hg19\tGene coordinates and strand" >> "${META}"

GOA_URL="https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz"
download_if_missing "${GOA_URL}" "${RAW}/goa_human.gaf.gz"
echo -e "goa_human\t${GOA_URL}\t${RAW}/goa_human.gaf.gz\tgene_symbols\tGO biological process annotations" >> "${META}"

HOMOLOGY_URL="https://ftp.ensembl.org/pub/grch37/release-75/tsv/homo_sapiens/Homo_sapiens.GRCh37.75.homologies.tsv.gz"
if curl -L --retry 2 --fail -o "${RAW}/Homo_sapiens.GRCh37.75.homologies.tsv.gz.tmp" "${HOMOLOGY_URL}"; then
  mv "${RAW}/Homo_sapiens.GRCh37.75.homologies.tsv.gz.tmp" "${RAW}/Homo_sapiens.GRCh37.75.homologies.tsv.gz"
  echo -e "ensembl_homologies\t${HOMOLOGY_URL}\t${RAW}/Homo_sapiens.GRCh37.75.homologies.tsv.gz\tGRCh37/hg19\tParalog source" >> "${META}"
else
  rm -f "${RAW}/Homo_sapiens.GRCh37.75.homologies.tsv.gz.tmp"
  echo -e "ensembl_homologies\t${HOMOLOGY_URL}\tMISSING\tGRCh37/hg19\tRelease-75 FTP download failed; prepare script falls back to BioMart GRCh37 human paralog export when available" >> "${META}"
fi

CTCF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz"
download_if_missing "${CTCF_URL}" "${RAW}/wgEncodeRegTfbsClusteredV3.bed.gz"
echo -e "ucsc_encode_tfbs_clustered_v3_ctcf\t${CTCF_URL}\t${RAW}/wgEncodeRegTfbsClusteredV3.bed.gz\tGRCh37/hg19\tFallback broad UCSC/ENCODE TFBS clustered V3; prepare script only uses it if GM12878 CTCF is missing" >> "${META}"

GM12878_CTCF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Ctcfsc15914c20StdPk.narrowPeak.gz"
download_if_missing "${GM12878_CTCF_URL}" "${RAW}/wgEncodeSydhTfbsGm12878Ctcfsc15914c20StdPk.narrowPeak.gz"
echo -e "encode_gm12878_ctcf_narrowpeak\t${GM12878_CTCF_URL}\t${RAW}/wgEncodeSydhTfbsGm12878Ctcfsc15914c20StdPk.narrowPeak.gz\tGRCh37/hg19\tGM12878 CTCF narrowPeak; used as immune/blood-derived CTCF reference for OneK1K cell types" >> "${META}"

TAD_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_primary%2Breplicate_Arrowhead_domainlist.txt.gz"
download_if_missing "${TAD_URL}" "${RAW}/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz"
echo -e "rao2014_gm12878_arrowhead_tads\t${TAD_URL}\t${RAW}/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz\tGRCh37/hg19\tGM12878 primary+replicate Arrowhead domains from GSE63525; prepare script converts domain starts/ends to boundaries" >> "${META}"

NASSER_ABC_URL="https://mitra.stanford.edu/engreitz/oak/public/Nasser2021/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"
download_if_missing "${NASSER_ABC_URL}" "${RAW}/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"
echo -e "nasser2021_avg_hic_abc\t${NASSER_ABC_URL}\t${RAW}/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz\tGRCh37/hg19\tPrimary ABC source matching pcQTL Figure 1E code; prepare script keeps ABC.Score > 0.1 and class != promoter" >> "${META}"

ABC_URL="https://www.encodeproject.org/files/ENCFF811DOR/@@download/ENCFF811DOR.bed.gz"
download_if_missing "${ABC_URL}" "${RAW}/ENCFF811DOR_GM12878_ABC_thresholded_GRCh38.bed.gz"
echo -e "encode_gm12878_abc_thresholded\t${ABC_URL}\t${RAW}/ENCFF811DOR_GM12878_ABC_thresholded_GRCh38.bed.gz\tGRCh38_lifted_to_hg19\tFallback ENCODE GM12878 ABC thresholded element-gene links; used only if Nasser2021 ABC is missing" >> "${META}"

LIFTOVER_URL="https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver"
if [[ ! -x "${RAW}/liftOver" ]]; then
  echo "[download] ${LIFTOVER_URL}"
  curl -L --retry 3 --fail -o "${RAW}/liftOver" "${LIFTOVER_URL}"
  chmod +x "${RAW}/liftOver"
fi
echo -e "ucsc_liftover_binary\t${LIFTOVER_URL}\t${RAW}/liftOver\ttool\tUCSC liftOver binary used for ABC GRCh38 to hg19 conversion" >> "${META}"

CHAIN_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
download_if_missing "${CHAIN_URL}" "${RAW}/hg38ToHg19.over.chain.gz"
echo -e "ucsc_hg38_to_hg19_chain\t${CHAIN_URL}\t${RAW}/hg38ToHg19.over.chain.gz\thg38_to_hg19\tUCSC chain used for ABC coordinate conversion" >> "${META}"

cat > "${RAW}/README_optional_regulatory_annotations.txt" <<'EOF'
Regulatory annotations downloaded by this workflow:

1. ABC enhancer-gene links:
   Primary: Nasser2021 AllPredictions AvgHiC ABC full file, matching the pcQTL Figure 1E
   code path. The prepare script keeps ABC.Score > 0.1 and class != promoter.
   Fallback: ENCODE GM12878 thresholded ABC element-gene links (ENCFF811DOR), downloaded
   in GRCh38 and lifted to hg19.

2. CTCF peaks:
   Primary: ENCODE/SYDH GM12878 CTCF narrowPeak in hg19, selected as an immune/blood-derived
   reference for OneK1K. Fallback: UCSC/ENCODE hg19 TFBS clustered V3 filtered to CTCF.

3. TAD boundaries:
   Rao et al. GM12878 primary+replicate Arrowhead domains from GSE63525.
   Domain starts and ends are treated as boundaries and expanded by config BOUNDARY_WINDOW_BP.
EOF
echo "[download] core annotation download complete"
