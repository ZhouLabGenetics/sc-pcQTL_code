#!/usr/bin/env python3
"""Build dense FDR-significant cis-QTL SNP sets for pcQTL and eQTL.

The credible-set MaxPIP annotation
is restricted to 95%-credible-set SNPs (union ~41k) and is too sparse for
well-powered S-LDSC. This builds the standard "all significant molecular-QTL"
annotation instead: per molecular feature, BH-adjust the cis-SNP p-values and
keep SNPs with BH-FDR < threshold, then union across features and cell types.
The SAME procedure is applied to pcQTL and eQTL so the two annotations are
defined identically (fair comparison).

Significance definition matches the project standard
(`pcqtl_compare/01_collect_results.R::summarise_step2`): per-feature BH<0.05.
Feature = cluster x PC for pcQTL, gene for eQTL.

Sources (identical SAIGE-QTL singleVar schema:
CHR POS MarkerID Allele1 Allele2 [AC_Allele2] AF_Allele2 ... p.value [...] N):
  pcQTL : {UPSTREAM_CELLTYPES_DIR}/{ct}/pcQTL/step3_saige/step2/{cluster_id}/{PC}
  eQTL  : {ONEK1K_EQTL_ROOT}/cis_{CellType}.tar.gz  (per-gene *.singleVar.txt)

Outputs (under <out-dir>, default module 09 annotations/sig_cisqtl):
  sig_{qtl}_snps.tsv        chr pos a1 a2 af  (union, deduped by chr:pos:allele-pair)
  sig_cisqtl_snpset_qc.tsv  per-(qtl,celltype) feature/SNP counts + union totals

Run examples:
  # Reduced validation run: one pcQTL cell type and one eQTL archive
  python3 53_build_sig_cisqtl_snpsets.py --qtl pcQTL --celltype cd8_nc
  python3 53_build_sig_cisqtl_snpsets.py --qtl eQTL  --celltype CD8_NC
  # full run (all cell types, both QTL types) then merge partials
  python3 53_build_sig_cisqtl_snpsets.py --qtl pcQTL
  python3 53_build_sig_cisqtl_snpsets.py --qtl eQTL
  python3 53_build_sig_cisqtl_snpsets.py --merge
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import sys
import tarfile
from pathlib import Path

sys.dont_write_bytecode = True

import numpy as np
import pandas as pd

from work_root import work_root

RELEASE_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(RELEASE_ROOT / "bin"))
from celltype_eligibility import primary_celltypes, primary_eqtl_celltypes

MODULE_DIR = work_root()                                  # .../04_formal_colocalization/08_heritability_enrichment_maxpip
WORKFLOW_ROOT = MODULE_DIR.parents[1]                      # .../coQTL_workflow
UPSTREAM_CELLTYPES_DIR = Path(os.environ.get(
    "COQTL_UPSTREAM_CELLTYPES_DIR",
    str(WORKFLOW_ROOT / "03_analysis_celltypes" / "01_upstream_main_pipeline_add_cov" / "celltypes")))
ONEK1K_EQTL_ROOT = Path(os.environ.get(
    "ONEK1K_EQTL_ROOT", "/path/to/OneK1K/original_single_gene_eqtl_tarballs"))

PCQTL_CELLTYPES = primary_celltypes()
EQTL_CELLTYPES = primary_eqtl_celltypes()

USECOLS = ["CHR", "POS", "Allele1", "Allele2", "AF_Allele2", "p.value"]


def bh_significant(pvals: np.ndarray, fdr: float) -> np.ndarray:
    """Return a boolean mask of Benjamini-Hochberg q < fdr (per feature)."""
    p = np.asarray(pvals, dtype=float)
    ok = np.isfinite(p)
    out = np.zeros(p.shape, dtype=bool)
    if ok.sum() == 0:
        return out
    pv = p[ok]
    n = pv.size
    order = np.argsort(pv)
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1)
    q = pv * n / ranks
    # enforce monotonicity over increasing p
    qs = q[order]
    qs = np.minimum.accumulate(qs[::-1])[::-1]
    q[order] = qs
    out_ok = q < fdr
    idx = np.flatnonzero(ok)
    out[idx[out_ok]] = True
    return out


def read_feature(buf) -> pd.DataFrame | None:
    try:
        df = pd.read_csv(buf, sep="\t", usecols=lambda c: c in USECOLS,
                         dtype={"CHR": str, "POS": "Int64", "Allele1": str,
                                "Allele2": str})
    except Exception:
        return None
    if df is None or df.empty or "p.value" not in df.columns:
        return None
    return df


def collect_significant(df: pd.DataFrame, fdr: float) -> pd.DataFrame:
    p = pd.to_numeric(df["p.value"], errors="coerce").to_numpy()
    mask = bh_significant(p, fdr)
    if not mask.any():
        return pd.DataFrame(columns=["chr", "pos", "a1", "a2", "af"])
    sub = df.loc[mask, ["CHR", "POS", "Allele1", "Allele2", "AF_Allele2"]].copy()
    sub.columns = ["chr", "pos", "a1", "a2", "af"]
    sub["chr"] = sub["chr"].astype(str).str.replace("^chr", "", regex=True)
    sub = sub.dropna(subset=["pos"])
    return sub


def dedupe(frames: list[pd.DataFrame]) -> pd.DataFrame:
    frames = [f for f in frames if f is not None and len(f) > 0]
    if not frames:
        return pd.DataFrame(columns=["chr", "pos", "a1", "a2", "af"])
    allsnps = pd.concat(frames, ignore_index=True)
    allsnps["pos"] = allsnps["pos"].astype("int64")
    # unordered allele pair key so the two QTL types match the EUR BIM regardless of orientation
    a = allsnps[["a1", "a2"]].to_numpy()
    lo = np.minimum(a[:, 0], a[:, 1])
    hi = np.maximum(a[:, 0], a[:, 1])
    allsnps["key"] = allsnps["chr"].astype(str) + ":" + allsnps["pos"].astype(str) + ":" + lo + ":" + hi
    allsnps = allsnps.drop_duplicates("key").drop(columns="key")
    return allsnps.reset_index(drop=True)


def filter_maf(df: pd.DataFrame, maf_min: float = 0.05) -> pd.DataFrame:
    """Restrict an allele-frequency table to MAF >= maf_min."""
    out = df.copy()
    af = pd.to_numeric(out["af"], errors="coerce")
    maf = np.minimum(af, 1.0 - af)
    return out.loc[np.isfinite(maf) & (maf >= maf_min)].reset_index(drop=True)


def run_pcqtl(celltypes, fdr, max_files):
    qc = []
    frames = []
    for ct in celltypes:
        step2 = UPSTREAM_CELLTYPES_DIR / ct / "pcQTL" / "step3_saige" / "step2"
        if not step2.is_dir():
            raise FileNotFoundError(f"Missing final pcQTL step2 directory for {ct}: {step2}")
        files = sorted(p for p in step2.glob("*/PC*") if p.is_file() and p.suffix != ".index")
        if max_files:
            files = files[:max_files]
        ct_frames = []
        nfeat = 0
        for f in files:
            df = read_feature(str(f))
            if df is None:
                continue
            nfeat += 1
            ct_frames.append(collect_significant(df, fdr))
        ct_union = dedupe(ct_frames)
        qc.append({"qtl": "pcQTL", "celltype": ct, "n_features": nfeat,
                   "n_sig_snps_ct_union": len(ct_union)})
        frames.append(ct_union)
        sys.stderr.write(f"[pcQTL] {ct}: {nfeat} features -> {len(ct_union)} sig SNPs (ct union)\n")
    return dedupe(frames), qc


def run_eqtl(celltypes, fdr, max_files):
    qc = []
    frames = []
    tarballs = [ONEK1K_EQTL_ROOT / f"cis_{c}.tar.gz" for c in celltypes]
    for tb in tarballs:
        if not tb.is_file():
            raise FileNotFoundError(f"Missing final OneK1K eQTL tarball: {tb}")
        ct = tb.name[len("cis_"):-len(".tar.gz")]
        ct_frames = []
        nfeat = 0
        with tarfile.open(tb, "r:gz") as tar:
            for member in tar:
                if not (member.isfile() and member.name.endswith(".singleVar.txt")):
                    continue
                if max_files and nfeat >= max_files:
                    break
                fobj = tar.extractfile(member)
                if fobj is None:
                    continue
                df = read_feature(io.TextIOWrapper(fobj, encoding="utf-8"))
                if df is None:
                    continue
                nfeat += 1
                ct_frames.append(collect_significant(df, fdr))
        ct_union = dedupe(ct_frames)
        qc.append({"qtl": "eQTL", "celltype": ct, "n_features": nfeat,
                   "n_sig_snps_ct_union": len(ct_union)})
        frames.append(ct_union)
        sys.stderr.write(f"[eQTL] {ct}: {nfeat} genes -> {len(ct_union)} sig SNPs (ct union)\n")
    return dedupe(frames), qc


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--qtl", choices=["pcQTL", "eQTL"])
    ap.add_argument("--celltype", help="single cell type (pcQTL token e.g. cd8_nc; eQTL token e.g. CD8_NC)")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--max-files", type=int, default=0, help="cap features per cell type (testing)")
    ap.add_argument("--out-dir", default=str(MODULE_DIR / "annotations" / "sig_cisqtl"))
    ap.add_argument("--merge", action="store_true", help="merge per-qtl partial files into final union + QC")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    part_dir = out_dir / "partials"
    part_dir.mkdir(exist_ok=True)

    if args.merge:
        qc_rows = []
        celltype_qc_rows = []
        for qtl in ("pcQTL", "eQTL"):
            expected_tags = PCQTL_CELLTYPES if qtl == "pcQTL" else EQTL_CELLTYPES
            per_celltype = [part_dir / f"sig_{qtl}_{tag}.tsv" for tag in expected_tags]
            per_celltype_qc = [part_dir / f"qc_{qtl}_{tag}.tsv" for tag in expected_tags]
            present = [path for path in per_celltype if path.is_file()]
            serial_part = part_dir / f"sig_{qtl}_all.tsv"

            if len(present) == len(per_celltype):
                selected = per_celltype
                source_mode = "per-cell-type array partials"
                unexpected = sorted(set(part_dir.glob(f"sig_{qtl}_*.tsv")) - set(per_celltype))
                if unexpected:
                    raise ValueError(
                        f"Unexpected {qtl} partials would make provenance ambiguous: "
                        + ", ".join(str(path) for path in unexpected)
                    )
                missing_qc = [path for path in per_celltype_qc if not path.is_file()]
                if missing_qc:
                    raise FileNotFoundError(
                        f"Missing {len(missing_qc)} current {qtl} cell-type QC files: "
                        + ", ".join(str(path) for path in missing_qc)
                    )
                for path in per_celltype_qc:
                    celltype_qc_rows.extend(pd.read_csv(path, sep="\t").to_dict("records"))
            elif not present and serial_part.is_file():
                selected = [serial_part]
                source_mode = "all-cell-type serial partial"
            else:
                missing = [str(path) for path in per_celltype if not path.is_file()]
                raise FileNotFoundError(
                    f"Incomplete current {qtl} partial set; missing {len(missing)} files: "
                    + ", ".join(missing)
                )

            union = dedupe([
                pd.read_csv(path, sep="\t", dtype={"chr": str}) for path in selected
            ])
            union_maf05 = filter_maf(union, maf_min=0.05)
            union.to_csv(out_dir / f"sig_{qtl}_snps.tsv", sep="\t", index=False)
            union_maf05.to_csv(out_dir / f"sig_{qtl}_snps_maf05.tsv", sep="\t", index=False)
            qc_rows.append({
                "qtl": qtl,
                "n_partials": len(selected),
                "n_sig_snps_union": len(union),
                "n_sig_snps_union_maf_ge_0.05": len(union_maf05),
            })
            sys.stderr.write(
                f"[merge] {qtl}: {source_mode} -> {len(union)} union SNPs; "
                f"{len(union_maf05)} with MAF >= 0.05\n"
            )
        pd.DataFrame(qc_rows).to_csv(out_dir / "sig_cisqtl_snpset_qc.tsv", sep="\t", index=False)
        if celltype_qc_rows:
            pd.DataFrame(celltype_qc_rows).to_csv(
                out_dir / "sig_cisqtl_celltype_qc.tsv", sep="\t", index=False
            )
        return 0

    if not args.qtl:
        ap.error("provide --qtl pcQTL|eQTL (or --merge)")

    cts = [args.celltype] if args.celltype else (
        PCQTL_CELLTYPES if args.qtl == "pcQTL" else EQTL_CELLTYPES
    )
    if args.qtl == "pcQTL":
        union, qc = run_pcqtl(cts, args.fdr, args.max_files)
    else:
        union, qc = run_eqtl(cts, args.fdr, args.max_files)

    tag = args.celltype if args.celltype else "all"
    union.to_csv(part_dir / f"sig_{args.qtl}_{tag}.tsv", sep="\t", index=False)
    pd.DataFrame(qc).to_csv(part_dir / f"qc_{args.qtl}_{tag}.tsv", sep="\t", index=False)
    sys.stderr.write(f"[{args.qtl}] wrote {len(union)} union SNPs for '{tag}' to {part_dir}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
