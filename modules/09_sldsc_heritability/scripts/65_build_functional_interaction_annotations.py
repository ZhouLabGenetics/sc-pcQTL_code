#!/usr/bin/env python3
"""Build functional-class x QTL-type interaction annotations to test whether
pcQTL's promoter-proximal localization translates into higher disease per-SNP
heritability (tau*) than eQTL, separately in PROMOTER vs ENHANCER elements.

Four disjoint-by-function binary annotations on the 1000G EUR BIM (same SNP order
as the existing single_maf05 qtlsig annots, which we reuse for QTL membership):
  pcQTL_prom : pcQTL-sig AND promoter        eQTL_prom : eQTL-sig AND promoter
  pcQTL_enh  : pcQTL-sig AND enhancer(¬prom)  eQTL_enh  : eQTL-sig AND enhancer(¬prom)
promoter = Promoter_UCSC|TSS_Hoffman|H3K4me3_Trynka ; enhancer = (Enhancer_Andersson|
Enhancer_Hoffman|H3K27ac_Hnisz|SuperEnhancer_Hnisz) AND NOT promoter. Masks from
baselineLD_noQTL, matched to the BIM by SNP rsid. Usage: 65_...py <chr>
"""
import csv, gzip, os, sys
from pathlib import Path

from work_root import work_root


def md(): return work_root()
PROM = ["Promoter_UCSC", "TSS_Hoffman", "H3K4me3_Trynka"]
ENH = ["Enhancer_Andersson", "Enhancer_Hoffman", "H3K27ac_Hnisz", "SuperEnhancer_Hnisz"]


def snp_masks(bl_path):
    """SNP -> (promoter, enhancer) from baselineLD annot."""
    prom, enh = {}, {}
    with gzip.open(bl_path, "rt") as h:
        r = csv.reader(h, delimiter="\t"); hdr = next(r)
        si = hdr.index("SNP")
        pi = [hdr.index(c) for c in PROM if c in hdr]
        ei = [hdr.index(c) for c in ENH if c in hdr]
        for row in r:
            p = any(row[i] == "1" for i in pi)
            e = any(row[i] == "1" for i in ei)
            snp = row[si]
            prom[snp] = p; enh[snp] = e and not p
    return prom, enh


def qtl_annot_rows(path):
    """yield (CHR,BP,SNP,CM,val) from a single-col qtlsig annot."""
    with gzip.open(path, "rt") as h:
        r = csv.reader(h, delimiter="\t"); hdr = next(r)
        ci, bi, si, mi = hdr.index("CHR"), hdr.index("BP"), hdr.index("SNP"), hdr.index("CM")
        vi = len(hdr) - 1
        for row in r:
            yield row[ci], row[bi], row[si], row[mi], int(row[vi])


def write_annot(path, rows, col):
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with gzip.open(path, "wt", newline="") as h:
        w = csv.writer(h, delimiter="\t"); w.writerow(["CHR", "BP", "SNP", "CM", col])
        for c, bp, snp, cm, v in rows:
            w.writerow([c, bp, snp, cm, v]); n += v
    return n


def main():
    chrom = str(sys.argv[1]).replace("chr", "")
    m = md()
    rd = m / "resources/eur_sldsc_ref"
    bl = rd / "baselineLD_noQTL" / f"baselineLD_noQTL.{chrom}.annot.gz"
    pc_annot = m / f"resources/eur_sldsc_custom_annotations/single_maf05/QTLsig_pcQTL_gt0/qtlsig_pcQTL.{chrom}.annot.gz"
    eq_annot = m / f"resources/eur_sldsc_custom_annotations/single_maf05/QTLsig_eQTL_gt0/qtlsig_eQTL.{chrom}.annot.gz"
    outdir = m / "resources/eur_sldsc_custom_annotations/single_maf05_funcinteract"
    prom, enh = snp_masks(bl)
    # eQTL membership by SNP (pcQTL annot defines the row order = BIM order)
    eqmem = {snp: v for _, _, snp, _, v in qtl_annot_rows(eq_annot)}

    rows_pcp, rows_pce, rows_eqp, rows_eqe = [], [], [], []
    for c, bp, snp, cm, pv in qtl_annot_rows(pc_annot):
        ev = eqmem.get(snp, 0)
        p = 1 if prom.get(snp, False) else 0
        e = 1 if enh.get(snp, False) else 0
        rows_pcp.append((c, bp, snp, cm, pv & p)); rows_pce.append((c, bp, snp, cm, pv & e))
        rows_eqp.append((c, bp, snp, cm, ev & p)); rows_eqe.append((c, bp, snp, cm, ev & e))
    n1 = write_annot(outdir / "pcQTL_prom" / f"fi_pcQTL_prom.{chrom}.annot.gz", rows_pcp, "pcQTL_prom")
    n2 = write_annot(outdir / "pcQTL_enh" / f"fi_pcQTL_enh.{chrom}.annot.gz", rows_pce, "pcQTL_enh")
    n3 = write_annot(outdir / "eQTL_prom" / f"fi_eQTL_prom.{chrom}.annot.gz", rows_eqp, "eQTL_prom")
    n4 = write_annot(outdir / "eQTL_enh" / f"fi_eQTL_enh.{chrom}.annot.gz", rows_eqe, "eQTL_enh")
    qc = m / "qc/eur_sldsc_annotations"; qc.mkdir(parents=True, exist_ok=True)
    with (qc / f"funcinteract.{chrom}.qc.tsv").open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t"); w.writerow(["chr", "pcQTL_prom", "pcQTL_enh", "eQTL_prom", "eQTL_enh"])
        w.writerow([chrom, n1, n2, n3, n4])
    print(f"chr{chrom}: pcQTL_prom={n1} pcQTL_enh={n2} eQTL_prom={n3} eQTL_enh={n4}")


if __name__ == "__main__":
    main()
