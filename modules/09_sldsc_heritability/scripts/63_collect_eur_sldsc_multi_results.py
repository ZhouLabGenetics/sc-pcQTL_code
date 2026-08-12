#!/usr/bin/env python3
"""Collect corrected S-LDSC models with a noQTL baseline and MAF >= 0.05.

Computes tau* for each named custom annotation and, optionally, the
per-endpoint Delta-tau* for a designated pair with a block-jackknife SE.

tau* = Coefficient * sqrt(p(1-p)) * M / h2g (binary annotations); M from
baselineLD_noQTL M_5_50 col0. .part_delete custom columns are the LAST
len(names) columns in --ref-ld-chr / manifest 'prefixes' order = `names` order.
"""
import argparse, csv, math, re
from pathlib import Path
import numpy as np

from work_root import work_root


def md(): return work_root()
def read_tsv(p):
    with Path(p).open(newline="") as h: return list(csv.DictReader(h, delimiter="\t"))
def pf(v):
    try:
        return float("nan") if v in ("","NA","nan","NaN",None) else float(v)
    except Exception: return float("nan")
def parse_log(p):
    o={"h2g":float("nan"),"intercept":"","ratio":"","fin":False}
    if not Path(p).exists(): return o
    t=Path(p).read_text(errors="replace"); o["fin"]="Analysis finished" in t
    m=re.search(r"Intercept:\s+([-+0-9.eE]+)",t);  o["intercept"]=m.group(1) if m else ""
    m=re.search(r"Ratio:\s+([-+0-9.eE]+)",t);      o["ratio"]=m.group(1) if m else ""
    m=re.search(r"Total Observed scale h2:\s+([-+0-9.eE]+)",t); o["h2g"]=float(m.group(1)) if m else float("nan")
    return o
def total_M(bp):
    return sum(float(Path(f"{bp}{c}.l2.M_5_50").read_text().split()[0]) for c in range(1,23))
def sd_bin(p):
    return math.sqrt(p*(1-p)) if (math.isfinite(p) and 0<p<1) else float("nan")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--manifest",required=True)
    ap.add_argument("--names",required=True,help="semicolon annotation names in prefix order")
    ap.add_argument("--out-dir",required=True); ap.add_argument("--prefix",required=True)
    ap.add_argument("--delta-pair",default="",help="nameA,nameB -> per-endpoint Delta-tau*")
    ap.add_argument("--baseline-prefix",default=None)
    a=ap.parse_args()
    M=total_M(a.baseline_prefix or (md()/"resources/eur_sldsc_ref/baselineLD_noQTL/baselineLD_noQTL."))
    names=[x.strip() for x in a.names.split(";") if x.strip()]; nC=len(names)
    out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    man=Path(a.manifest);  man=man if man.is_absolute() else md()/man
    rows=[]; deltas=[]
    pair=[x.strip() for x in a.delta_pair.split(",")] if a.delta_pair else None
    for task in read_tsv(man):
        op=Path(task["out_prefix"]); res=Path(str(op)+".results"); log=parse_log(Path(str(op)+".log"))
        pdl=Path(str(op)+".part_delete")
        if not res.exists() or not log["fin"]: continue
        h2g=log["h2g"]; recs=read_tsv(res); per={}
        for nm in names:
            r=next((x for x in recs if nm in x["Category"]),None)
            if not r: continue
            p=pf(r.get("Prop._SNPs")); coef=pf(r.get("Coefficient")); cse=pf(r.get("Coefficient_std_error"))
            sc=(sd_bin(p)*M/h2g) if (math.isfinite(sd_bin(p)) and math.isfinite(h2g) and h2g) else float("nan")
            per[nm]={"tau_star":coef*sc,"tau_star_se":cse*sc,"scale":sc,"prop":p}
            rows.append({"phenocode":task["phenocode"],"phenotype":task["phenotype"],
                         "category_group":task["category"],"annotation":nm,"M":M,"h2g":h2g,
                         "intercept":log["intercept"],"ratio":log["ratio"],
                         "tau_star":coef*sc,"tau_star_se":cse*sc,**{k:r.get(k) for k in
                         ("Prop._SNPs","Enrichment","Enrichment_std_error","Enrichment_p","Coefficient","Coefficient_std_error","Coefficient_z-score")}})
        if pair and all(n in per for n in pair):
            A,B=per[pair[0]],per[pair[1]]; d=A["tau_star"]-B["tau_star"]
            se=float("nan"); meth="independent_se"
            if pdl.exists():
                try:
                    arr=np.loadtxt(pdl)
                    if arr.ndim==2 and arr.shape[1]>=nC:
                        ia,ib=names.index(pair[0]),names.index(pair[1])
                        db=arr[:,-nC+ia]*A["scale"]-arr[:,-nC+ib]*B["scale"]
                        n=db.size; se=math.sqrt((n-1)/n*np.sum((db-db.mean())**2)); meth="block_jackknife"
                except Exception: pass
            if not math.isfinite(se): se=math.sqrt(A["tau_star_se"]**2+B["tau_star_se"]**2)
            z=d/se if se>0 else float("nan")
            deltas.append({"phenocode":task["phenocode"],"phenotype":task["phenotype"],
                           "category_group":task["category"],f"tau_star_{pair[0]}":A["tau_star"],
                           f"tau_star_{pair[1]}":B["tau_star"],"delta_tau_star":d,"delta_se":se,
                           "delta_z":z,"se_method":meth})
    if rows:
        with (out/f"{a.prefix}.tsv").open("w",newline="") as h:
            w=csv.DictWriter(h,delimiter="\t",fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    if deltas:
        with (out/f"{a.prefix}_delta.tsv").open("w",newline="") as h:
            w=csv.DictWriter(h,delimiter="\t",fieldnames=list(deltas[0].keys())); w.writeheader(); w.writerows(deltas)
    print(f"M={M:.0f}; {len(rows)} annotation rows, {len(deltas)} delta rows; names={names}")


if __name__=="__main__":
    main()
