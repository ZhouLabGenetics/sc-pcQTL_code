#!/usr/bin/env python3
"""Build disjoint co-regulation-degree annotations from the n_genes-binned
pcQTL sig-SNP sets (script 66). A SNP significant across clusters of different
sizes is assigned to its MAX cluster-size bin (tests 'tags the highest co-regulation
it can'). Three mutually-exclusive binary annotations on the 1000G EUR BIM:
  ng_small : max cluster size 2         ng_med : max size 3-4        ng_large : max size >=5
Usage: 67_build_pcqtl_ngenes_annotations.py <chr>
"""
import csv, gzip, os, sys
from pathlib import Path

from work_root import work_root


def md(): return work_root()


def load_keys(path):
    keys = set()
    with Path(path).open(newline="") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            a1, a2 = str(r["a1"]).upper(), str(r["a2"]).upper()
            lo, hi = (a1, a2) if a1 <= a2 else (a2, a1)
            keys.add((str(r["chr"]).replace("chr", ""), str(r["pos"]).strip(), lo, hi))
    return keys


def bim_rows(path):
    with Path(path).open() as h:
        for line in h:
            c, snp, cm, bp, a1, a2 = line.rstrip("\n").split()[:6]
            yield c.replace("chr", ""), snp, cm, bp, a1.upper(), a2.upper()


def write_one(rows, out_path, col):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = ns = 0
    with gzip.open(out_path, "wt", newline="") as h:
        w = csv.writer(h, delimiter="\t"); w.writerow(["CHR", "BP", "SNP", "CM", col])
        for c, snp, cm, bp, v in rows:
            w.writerow([c, bp, snp, cm, v]); n += 1; ns += v
    return n, ns


def main():
    chrom = str(sys.argv[1]).replace("chr", "")
    m = md(); sig = m / "annotations/sig_cisqtl"
    plink = Path(os.environ.get("EUR_PLINK_DIR", m / "resources/eur_sldsc_ref/1000G_EUR_Phase3_plink"))
    small = load_keys(sig / "pcQTL_ngenes_small_maf05.tsv")
    med = load_keys(sig / "pcQTL_ngenes_med_maf05.tsv")
    large = load_keys(sig / "pcQTL_ngenes_large_maf05.tsv")
    outdir = m / "resources/eur_sldsc_custom_annotations/single_maf05_ngenes"
    bim = plink / f"1000G.EUR.QC.{chrom}.bim"
    if not bim.exists(): raise FileNotFoundError(bim)

    # exclude the (extended) MHC region: large co-regulated clusters are heavily
    # MHC-concentrated (chr6), which has atypical LD + extreme immune enrichment and
    # is standardly excluded from S-LDSC. Removing it tests co-regulation degree
    # cleanly (not an MHC artifact). Extended MHC hg19: chr6:25-35 Mb.
    def in_mhc(c, bp):
        return c == "6" and 25_000_000 <= int(bp) <= 35_000_000

    def classify(which):
        for c, snp, cm, bp, a1, a2 in bim_rows(bim):
            if in_mhc(c, bp):
                yield c, snp, cm, bp, 0; continue
            lo, hi = (a1, a2) if a1 <= a2 else (a2, a1)
            k = (c, bp, lo, hi)
            is_large = k in large
            is_med = (k in med) and not is_large
            is_small = (k in small) and not is_med and not is_large  # disjoint by MAX size
            v = {"large": is_large, "med": is_med, "small": is_small}[which]
            yield c, snp, cm, bp, int(v)

    res = {}
    for b in ["small", "med", "large"]:
        n, ns = write_one(classify(b), outdir / f"ng_{b}" / f"ng_{b}.{chrom}.annot.gz", f"ng_{b}")
        res[b] = ns
    qc = m / "qc/eur_sldsc_annotations"; qc.mkdir(parents=True, exist_ok=True)
    with (qc / f"ngenes.{chrom}.qc.tsv").open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t"); w.writerow(["chr", "ng_small", "ng_med", "ng_large"]); w.writerow([chrom, res["small"], res["med"], res["large"]])
    print(f"chr{chrom}: ng_small={res['small']} ng_med={res['med']} ng_large={res['large']}")


if __name__ == "__main__":
    main()
