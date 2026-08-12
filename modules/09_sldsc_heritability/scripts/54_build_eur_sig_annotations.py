#!/usr/bin/env python3
"""Build EUR-reference DENSE significant-cis-QTL binary annotations for LDSC.

Consumes the dense FDR-significant SNP sets from 53_build_sig_cisqtl_snpsets.py
and emits TWO single-column binary annotation files per chromosome (one per QTL
type), matched in definition, over the 1000G EUR Phase3 BIM:

  single/QTLsig_pcQTL_gt0/qtlsig_pcQTL.{chr}.annot.gz   col QTLsig_pcQTL_gt0
  single/QTLsig_eQTL_gt0/qtlsig_eQTL.{chr}.annot.gz     col QTLsig_eQTL_gt0

Single-column output keeps it compatible with 39_align_single_annotation_*.py
and the two-prefix joint runner (mirrors the existing window-model layout).
SNP matching is chromosome + position + UNORDERED allele pair (same key as 53).

Usage: 54_build_eur_sig_annotations.py <chr>
Env: EUR_PLINK_DIR, QTLSIG_PCQTL_SOURCE, QTLSIG_EQTL_SOURCE, EUR_QTLSIG_SINGLE_DIR
"""
import csv
import gzip
import os
import sys
from pathlib import Path

from work_root import work_root


def module_dir_from_script() -> Path:
    return work_root()


def load_sig_keys(path: Path) -> set:
    keys = set()
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chrom = str(row["chr"]).replace("chr", "")
            pos = str(row["pos"]).strip()
            a1, a2 = str(row["a1"]).upper(), str(row["a2"]).upper()
            lo, hi = (a1, a2) if a1 <= a2 else (a2, a1)
            keys.add((chrom, pos, lo, hi))
    return keys


def bim_rows(path: Path):
    with path.open() as handle:
        for line in handle:
            chrom, snp, cm, bp, a1, a2 = line.rstrip("\n").split()[:6]
            yield chrom.replace("chr", ""), snp, cm, bp, a1.upper(), a2.upper()


def write_annot(bim, keys, out_path: Path, col: str):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = nsig = 0
    with gzip.open(out_path, "wt", newline="") as handle:
        w = csv.DictWriter(handle, delimiter="\t", fieldnames=["CHR", "BP", "SNP", "CM", col])
        w.writeheader()
        for chrom_bim, snp, cm, bp, a1, a2 in bim:
            lo, hi = (a1, a2) if a1 <= a2 else (a2, a1)
            v = int((chrom_bim, bp, lo, hi) in keys)
            n += 1
            nsig += v
            w.writerow({"CHR": chrom_bim, "BP": bp, "SNP": snp, "CM": cm, col: v})
    return n, nsig


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: 54_build_eur_sig_annotations.py <chr>")
    chrom = str(sys.argv[1]).replace("chr", "")
    module_dir = module_dir_from_script()
    plink_dir = Path(os.environ.get(
        "EUR_PLINK_DIR",
        module_dir / "resources" / "eur_sldsc_ref" / "1000G_EUR_Phase3_plink"))
    pc_src = Path(os.environ.get(
        "QTLSIG_PCQTL_SOURCE",
        module_dir / "annotations" / "sig_cisqtl" / "sig_pcQTL_snps.tsv"))
    eq_src = Path(os.environ.get(
        "QTLSIG_EQTL_SOURCE",
        module_dir / "annotations" / "sig_cisqtl" / "sig_eQTL_snps.tsv"))
    single_dir = Path(os.environ.get(
        "EUR_QTLSIG_SINGLE_DIR",
        module_dir / "resources" / "eur_sldsc_custom_annotations" / "single"))
    qc_dir = module_dir / "qc" / "eur_sldsc_annotations"
    qc_dir.mkdir(parents=True, exist_ok=True)

    bim = plink_dir / f"1000G.EUR.QC.{chrom}.bim"
    if not bim.exists():
        raise FileNotFoundError(bim)
    pc_keys = load_sig_keys(pc_src)
    eq_keys = load_sig_keys(eq_src)

    n, n_pc = write_annot(bim_rows(bim), pc_keys,
                          single_dir / "QTLsig_pcQTL_gt0" / f"qtlsig_pcQTL.{chrom}.annot.gz",
                          "QTLsig_pcQTL_gt0")
    _, n_eq = write_annot(bim_rows(bim), eq_keys,
                          single_dir / "QTLsig_eQTL_gt0" / f"qtlsig_eQTL.{chrom}.annot.gz",
                          "QTLsig_eQTL_gt0")
    with (qc_dir / f"qtlsig.{chrom}.qc.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t")
        w.writerow(["chr", "n_bim_snps", "n_pcqtl_sig", "n_eqtl_sig"])
        w.writerow([chrom, n, n_pc, n_eq])
    print(f"chr{chrom}: bim={n} pcQTL_sig={n_pc} eQTL_sig={n_eq}")


if __name__ == "__main__":
    main()
