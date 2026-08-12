#!/usr/bin/env python3
"""Stratify pcQTL FDR-significant SNPs by the co-regulation degree of their
cluster (number of genes in the cluster x PC), to test whether variants driving
MORE co-regulated genes carry more disease per-SNP heritability -- the core
multi-gene premise of pcQTL. Reuses script 53's per-feature BH extraction; adds an
n_genes lookup from each cell type's PCA master summary, keyed by cell type and
cluster, and pools significant SNPs into three bins: 2 genes, 3-4 genes, and
at least 5 genes. MAF >= 0.05 is applied to match the main annotation design.

Output: annotations/sig_cisqtl/pcQTL_ngenes_{small,med,large}_maf05.tsv
"""
import csv, importlib.util, re, sys
from pathlib import Path
import numpy as np
import pandas as pd

CODE_ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("s53", CODE_ROOT/"scripts/53_build_sig_cisqtl_snpsets.py")
s53 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s53)
WORK_ROOT = s53.MODULE_DIR

FDR = 0.05

def load_ngenes():
    """Complete per-cluster n_genes from each cell type's PCA master summary
    ({ct}/pcQTL/step2_pca/all_clusters_summary.tsv) -> (celltype, cluster_id)->n_genes.
    n_genes is a per-cluster property (constant across that cluster's PCs)."""
    m = {}
    for ct in s53.PCQTL_CELLTYPES:
        f = s53.UPSTREAM_CELLTYPES_DIR / ct / "pcQTL" / "step2_pca" / "all_clusters_summary.tsv"
        if not f.exists():
            sys.stderr.write(f"missing cluster summary: {f}\n"); continue
        with f.open() as h:
            for r in csv.DictReader(h, delimiter="\t"):
                try: m[(ct, r["cluster_id"])] = int(float(r["n_genes"]))
                except Exception: pass
    return m

def bin_of(n):
    if n <= 2: return "small"
    if n <= 4: return "med"
    return "large"

def main():
    ng = load_ngenes()
    print(f"n_genes map: {len(ng)} (celltype,cluster) entries", flush=True)
    bins = {"small": [], "med": [], "large": []}
    n_feat = {"small":0,"med":0,"large":0}; n_nolookup = 0
    for ct in s53.PCQTL_CELLTYPES:
        step2 = s53.UPSTREAM_CELLTYPES_DIR / ct / "pcQTL" / "step3_saige" / "step2"
        if not step2.is_dir():
            sys.stderr.write(f"missing {step2}\n"); continue
        files = [p for p in step2.glob("*/PC*") if p.is_file() and p.suffix != ".index"]
        for f in files:
            cluster_id = f.parent.name
            n = ng.get((ct, cluster_id))
            if n is None: n_nolookup += 1; continue
            df = s53.read_feature(str(f))
            if df is None: continue
            sig = s53.collect_significant(df, FDR)
            if len(sig) == 0: continue
            b = bin_of(n); bins[b].append(sig); n_feat[b]+=1
        print(f"  {ct}: done", flush=True)
    outdir = WORK_ROOT/"annotations/sig_cisqtl"; outdir.mkdir(parents=True, exist_ok=True)
    qc = []
    for b, frames in bins.items():
        ded = s53.dedupe(frames)
        # MAF>=0.05 filter (match the maf05 comparison)
        af = pd.to_numeric(ded["af"], errors="coerce")
        ded = ded[(af >= 0.05) & (af <= 0.95)].reset_index(drop=True)
        out = outdir/f"pcQTL_ngenes_{b}_maf05.tsv"
        ded.to_csv(out, sep="\t", index=False)
        qc.append((b, n_feat[b], len(ded)))
        print(f"bin {b}: {n_feat[b]} cluster×PC features -> {len(ded)} maf05 sig SNPs -> {out.name}", flush=True)
    print(f"features with no n_genes lookup (skipped): {n_nolookup}")
    with (outdir/"pcQTL_ngenes_qc.tsv").open("w",newline="") as h:
        w=csv.writer(h,delimiter="\t"); w.writerow(["bin","n_features","n_snps_maf05"]); w.writerows(qc)

if __name__ == "__main__":
    main()
