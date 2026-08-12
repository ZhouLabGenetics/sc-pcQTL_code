#!/usr/bin/env python3
"""Paper-style (Hujoel et al. 2019 AJHG / Finucane et al. 2018) random-effects
meta-analysis of S-LDSC heritability ENRICHMENT and standardized effect size
TAU* across traits, reported per annotation.

Method match:
  - Enrichment = %h2(c)/%SNP(c); meta = DerSimonian-Laird random-effects across
    traits (= R rmeta::meta.summaries, as in Finucane/Hujoel). Reported as the
    headline descriptive "Nx enriched" number with a 95% CI.
  - tau* = M * sd_c * tau_c / h2 (standardized effect size; conditional, unique to
    the focal annotation). DL random-effects meta across traits + z-test p-value.
  - We HONESTLY report both for each annotation (no requirement that pcQTL > eQTL).

Input : a collector .tsv from script 63 with per-annotation rows
        (cols: phenocode, phenotype, annotation, Enrichment, Enrichment_std_error,
         Enrichment_p, tau_star, tau_star_se, ...).
Output: {out_dir}/{prefix}_annotation_meta.tsv  (one row per annotation)
        {out_dir}/{prefix}_top_traits.tsv        (top traits per annotation)
        + a stdout report in paper style.
"""
import argparse, csv, math
from pathlib import Path
import numpy as np

from work_root import work_root


def md(): return work_root()
def load(p):
    with open(p, newline="") as h: return list(csv.DictReader(h, delimiter="\t"))
def f(x):
    try: return float("nan") if x in ("", "NA", "nan", "NaN", None) else float(x)
    except Exception: return float("nan")
def Phi(x): return 0.5 * (1 + math.erf(x / math.sqrt(2)))
def zp(z): return 2 * (1 - Phi(abs(z)))


def dl_meta(est, se):
    """DerSimonian-Laird random-effects meta (matches rmeta::meta.summaries)."""
    est = np.asarray(est, float); se = np.asarray(se, float)
    ok = np.isfinite(est) & np.isfinite(se) & (se > 0)
    est, se = est[ok], se[ok]; k = est.size
    if k == 0: return None
    w = 1.0 / se ** 2
    fe = float(np.sum(w * est) / np.sum(w))
    Q = float(np.sum(w * (est - fe) ** 2)); dfree = k - 1
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (Q - dfree) / c) if c > 0 else 0.0
    I2 = max(0.0, (Q - dfree) / Q) * 100 if Q > 0 else 0.0
    wr = 1.0 / (se ** 2 + tau2)
    re = float(np.sum(wr * est) / np.sum(wr)); re_se = float(math.sqrt(1.0 / np.sum(wr)))
    return dict(k=k, re=re, re_se=re_se, lo=re - 1.96 * re_se, hi=re + 1.96 * re_se,
                Q=Q, I2=I2, tau2=tau2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--collector", required=True, help="per-annotation .tsv from script 63")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--label-map", default="QTLsig_pcQTL_gt0=pcQTL;QTLsig_eQTL_gt0=eQTL",
                    help="semicolon list of annotation=display")
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--trait-manifest", default="",
                    help="optional TSV whose phenocode column defines the meta-analysis set")
    ap.add_argument("--analysis-label", default="all_available_traits")
    ap.add_argument("--output-suffix", default="",
                    help="optional filename suffix, without the leading underscore")
    a = ap.parse_args()
    coll = Path(a.collector); coll = coll if coll.is_absolute() else md() / coll
    rows = load(coll)
    selected = None
    if a.trait_manifest:
        trait_path = Path(a.trait_manifest)
        trait_path = trait_path if trait_path.is_absolute() else md() / trait_path
        selected_rows = load(trait_path)
        selected = {row["phenocode"] for row in selected_rows}
        if not selected:
            raise ValueError(f"Trait manifest selects no endpoints: {trait_path}")
        rows = [row for row in rows if row.get("phenocode") in selected]
        if not rows:
            raise ValueError(f"No collector rows match the selected endpoints in {trait_path}")
    lab = dict(kv.split("=") for kv in a.label_map.split(";") if "=" in kv)
    out = Path(a.out_dir); out = out if out.is_absolute() else md() / out
    out.mkdir(parents=True, exist_ok=True)

    anns = []
    for r in rows:
        if r["annotation"] not in anns: anns.append(r["annotation"])

    meta_rows = []; top_rows = []
    print(f"Paper-style random-effects (DerSimonian-Laird) meta across traits — {coll.name}")
    print(f"{'annotation':14s} {'k':>4s} {'Enrichment (95% CI)':>26s} {'tau* (95% CI)':>24s} {'tau*_p':>9s}")
    for ann in anns:
        ar = [r for r in rows if r["annotation"] == ann]
        em = dl_meta([f(r["Enrichment"]) for r in ar], [f(r.get("Enrichment_std_error")) for r in ar])
        tm = dl_meta([f(r["tau_star"]) for r in ar], [f(r["tau_star_se"]) for r in ar])
        disp = lab.get(ann, ann)
        if em and tm:
            tz = tm["re"] / tm["re_se"]; tp = zp(tz)
            n_enr_sig = sum(1 for r in ar if math.isfinite(f(r.get("Enrichment_p"))) and f(r["Enrichment_p"]) < 0.05 and f(r["Enrichment"]) > 1)
            meta_rows.append({"analysis_set": a.analysis_label,
                              "n_traits_selected": len(selected) if selected is not None else len({r["phenocode"] for r in rows}),
                              "annotation": disp, "k_traits": em["k"],
                              "enrichment_re": em["re"], "enrichment_ci_lo": em["lo"], "enrichment_ci_hi": em["hi"],
                              "enrichment_I2_pct": em["I2"], "n_traits_enrich_sig_p05": n_enr_sig,
                              "tau_star_re": tm["re"], "tau_star_se": tm["re_se"], "tau_star_ci_lo": tm["lo"],
                              "tau_star_ci_hi": tm["hi"], "tau_star_z": tz, "tau_star_p": tp, "tau_star_I2_pct": tm["I2"]})
            print(f"{disp:14s} {em['k']:>4d} {em['re']:>9.2f} [{em['lo']:>5.2f},{em['hi']:>5.2f}]   "
                  f"{tm['re']:>+8.4f} [{tm['lo']:>+6.3f},{tm['hi']:>+6.3f}] {tp:>9.2g}")
        # top traits for this annotation (by enrichment among adequately-powered, h2g>0.01)
        cand = [r for r in ar if math.isfinite(f(r.get("h2g"))) and f(r["h2g"]) > 0.01
                and math.isfinite(f(r["Enrichment"])) and math.isfinite(f(r.get("Enrichment_p")))]
        cand.sort(key=lambda r: f(r["Enrichment_p"]))
        for r in cand[:a.top]:
            top_rows.append({"analysis_set": a.analysis_label,
                             "annotation": disp, "phenocode": r["phenocode"], "phenotype": r["phenotype"],
                             "h2g": f(r["h2g"]), "Enrichment": f(r["Enrichment"]),
                             "Enrichment_se": f(r.get("Enrichment_std_error")), "Enrichment_p": f(r["Enrichment_p"]),
                             "tau_star": f(r["tau_star"]), "tau_star_z": f(r.get("Coefficient_z-score"))})

    suffix = f"_{a.output_suffix}" if a.output_suffix else ""
    if meta_rows:
        with (out / f"{a.prefix}_annotation_meta{suffix}.tsv").open("w", newline="") as h:
            w = csv.DictWriter(h, delimiter="\t", fieldnames=list(meta_rows[0].keys())); w.writeheader(); w.writerows(meta_rows)
    if top_rows:
        with (out / f"{a.prefix}_top_traits{suffix}.tsv").open("w", newline="") as h:
            w = csv.DictWriter(h, delimiter="\t", fieldnames=list(top_rows[0].keys())); w.writeheader(); w.writerows(top_rows)
    print(f"\nwrote {a.prefix}_annotation_meta{suffix}.tsv ({len(meta_rows)} annotations) + "
          f"{a.prefix}_top_traits{suffix}.tsv ({len(top_rows)} rows); analysis={a.analysis_label}")


if __name__ == "__main__":
    main()
