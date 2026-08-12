#!/usr/bin/env python3
"""Meta-analyse the per-endpoint Delta-tau* (pcQTL minus eQTL) across FinnGen
endpoints to a single headline estimate. Fixed-effect inverse-variance and
DerSimonian-Laird random-effects, with a symmetric 95% CI (reported whatever the
sign) and an I^2 heterogeneity summary.

Input : {summary_dir}/{prefix}_delta.tsv  (from 58)
Output: {summary_dir}/{prefix}_delta_meta.tsv  + a short stdout report
"""
import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.stats import t as student_t

from work_root import work_root


def module_dir_from_script() -> Path:
    return work_root()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary-dir", default="results/heritability_enrichment/eur_sldsc_qtlsig_joint_noqtl_baseline_summary")
    ap.add_argument("--prefix", default="finngen_r12_eur_sldsc_qtlsig_joint_noqtl")
    ap.add_argument("--trait-manifest", default="",
                    help="optional TSV whose phenocode column defines the meta-analysis set")
    ap.add_argument("--analysis-label", default="all_available_traits")
    ap.add_argument("--output-suffix", default="")
    args = ap.parse_args()

    module_dir = module_dir_from_script()
    sdir = Path(args.summary_dir)
    if not sdir.is_absolute():
        sdir = module_dir / sdir
    delta_path = sdir / f"{args.prefix}_delta.tsv"
    if not delta_path.exists():
        raise FileNotFoundError(delta_path)

    with delta_path.open(newline="") as h:
        rows = list(csv.DictReader(h, delimiter="\t"))
    selected = None
    if args.trait_manifest:
        trait_path = Path(args.trait_manifest)
        if not trait_path.is_absolute():
            trait_path = module_dir / trait_path
        with trait_path.open(newline="") as h:
            selected = {row["phenocode"] for row in csv.DictReader(h, delimiter="\t")}
        if not selected:
            raise ValueError(f"Trait manifest selects no endpoints: {trait_path}")
        rows = [row for row in rows if row.get("phenocode") in selected]
        if not rows:
            raise ValueError(f"No Delta-tau* rows match {trait_path}")
    d = np.array([float(r["delta_tau_star"]) for r in rows], dtype=float)
    se = np.array([float(r["delta_se"]) for r in rows], dtype=float)
    ok = np.isfinite(d) & np.isfinite(se) & (se > 0)
    d, se = d[ok], se[ok]
    k = d.size
    if k == 0:
        raise SystemExit("no usable per-endpoint Delta-tau* rows")

    w = 1.0 / se ** 2
    d_fe = float(np.sum(w * d) / np.sum(w))
    se_fe = float(math.sqrt(1.0 / np.sum(w)))
    Q = float(np.sum(w * (d - d_fe) ** 2))
    df = k - 1
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (Q - df) / c) if c > 0 else 0.0
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    wr = 1.0 / (se ** 2 + tau2)
    d_re = float(np.sum(wr * d) / np.sum(wr))
    se_re = float(math.sqrt(1.0 / np.sum(wr)))

    def zp(est, s):
        z = est / s
        p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
        return z, p

    z_fe, p_fe = zp(d_fe, se_fe)
    z_re, p_re = zp(d_re, se_re)

    # Knapp-Hartung small-sample correction for the random-effects model
    # (recommended over uncorrected DL at small k; t reference, df=k-1).
    se_kh = float(math.sqrt(np.sum(wr * (d - d_re) ** 2) / ((k - 1) * np.sum(wr)))) if k > 1 else float("nan")
    t_kh = d_re / se_kh if (math.isfinite(se_kh) and se_kh > 0) else float("nan")
    p_kh = 2 * (1 - float(student_t.cdf(abs(t_kh), k - 1))) if math.isfinite(t_kh) else float("nan")
    tc = float(student_t.ppf(0.975, k - 1))

    # assumption-light sign test for DIRECTION: one-sided P(X >= n_pos | Binom(k, 0.5))
    n_pos = int(np.sum(d > 0))
    sign_p_one_sided = sum(math.comb(k, i) for i in range(n_pos, k + 1)) / (2 ** k)

    out = [
        {"analysis_set": args.analysis_label, "n_traits_selected": len(selected) if selected is not None else k,
         "model": "fixed_effect", "k_endpoints": k, "delta_tau_star": d_fe, "se": se_fe,
         "ci95_lo": d_fe - 1.96 * se_fe, "ci95_hi": d_fe + 1.96 * se_fe, "z": z_fe, "p": p_fe,
         "Q": Q, "df": df, "I2_pct": I2, "tau2": tau2},
        {"analysis_set": args.analysis_label, "n_traits_selected": len(selected) if selected is not None else k,
         "model": "random_effect_DL", "k_endpoints": k, "delta_tau_star": d_re, "se": se_re,
         "ci95_lo": d_re - 1.96 * se_re, "ci95_hi": d_re + 1.96 * se_re, "z": z_re, "p": p_re,
         "Q": Q, "df": df, "I2_pct": I2, "tau2": tau2},
        {"analysis_set": args.analysis_label, "n_traits_selected": len(selected) if selected is not None else k,
         "model": "random_effect_KH", "k_endpoints": k, "delta_tau_star": d_re, "se": se_kh,
         "ci95_lo": d_re - tc * se_kh, "ci95_hi": d_re + tc * se_kh, "z": t_kh, "p": p_kh,
         "Q": Q, "df": df, "I2_pct": I2, "tau2": tau2},
        {"analysis_set": args.analysis_label, "n_traits_selected": len(selected) if selected is not None else k,
         "model": "sign_test_n_positive", "k_endpoints": k, "delta_tau_star": float(n_pos), "se": float("nan"),
         "ci95_lo": float("nan"), "ci95_hi": float("nan"), "z": float("nan"), "p": sign_p_one_sided,
         "Q": float("nan"), "df": float("nan"), "I2_pct": float("nan"), "tau2": float("nan")},
    ]
    suffix = f"_{args.output_suffix}" if args.output_suffix else ""
    out_path = sdir / f"{args.prefix}_delta_meta{suffix}.tsv"
    with out_path.open("w", newline="") as h:
        wcsv = csv.DictWriter(h, delimiter="\t", fieldnames=list(out[0].keys()))
        wcsv.writeheader(); wcsv.writerows(out)

    n_pos = int(np.sum(d > 0))
    print(f"Delta-tau* (pcQTL - eQTL) meta across k={k} endpoints "
          f"({n_pos} positive); positive => pcQTL > eQTL.")
    for r in out:
        print(f"  {r['model']:18s} Delta={r['delta_tau_star']:+.4g} "
              f"[95% CI {r['ci95_lo']:+.4g}, {r['ci95_hi']:+.4g}] p={r['p']:.3g}")
    print(f"  heterogeneity: Q={Q:.2f} df={df} I2={I2:.0f}% tau2={tau2:.3g}")
    print(f"  wrote {out_path}")


if __name__ == "__main__":
    main()
