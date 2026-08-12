#!/usr/bin/env python3
"""Generalized S-LDSC runner: condition on baselineLD_noQTL + an arbitrary list
of custom annotation LD-score prefixes (semicolon-separated 'prefixes' column in
the manifest). Handles the corrected joint (pcQTL+eQTL), promoter/enhancer
interaction, cluster-size, and marginal single-annotation models.

Forces EUR_BASELINE_PREFIX to baselineLD_noQTL by default (overridable via
QTLSIG_BASELINE_PREFIX) so it never silently inherits env.sh's full baseline.
"""
import argparse
import csv
import os
import subprocess
from pathlib import Path

from work_root import work_root


def module_dir_from_script() -> Path:
    return work_root()


def read_tsv(path: Path):
    with path.open(newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def select_task(rows, task_id):
    for r in rows:
        if r["task_id"] == task_id:
            return r
    raise SystemExit(f"No task_id={task_id}")


def has(p, s):  # log finished + no traceback
    return p.exists() and s in p.read_text(errors="replace")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("task_id")
    ap.add_argument("--manifest", required=True)
    args = ap.parse_args()
    md = module_dir_from_script()
    man = Path(args.manifest)
    if not man.is_absolute():
        man = md / man
    task = select_task(read_tsv(man), args.task_id)

    rd = Path(os.environ.get("EUR_SLDSC_RESOURCE_DIR", md / "resources" / "eur_sldsc_ref"))
    baseline = os.environ.get("QTLSIG_BASELINE_PREFIX",
                              os.environ.get("EUR_BASELINE_PREFIX",
                                             str(md / "resources" / "eur_sldsc_ref" / "baselineLD_noQTL" / "baselineLD_noQTL.")))
    # safety: if EUR_BASELINE_PREFIX was inherited as the FULL baseline, force noQTL
    if "baselineLD_noQTL" not in baseline:
        baseline = str(md / "resources" / "eur_sldsc_ref" / "baselineLD_noQTL" / "baselineLD_noQTL.")
    weight = os.environ.get("EUR_WEIGHT_PREFIX", str(rd / "1000G_Phase3_weights_hm3_no_MHC" / "weights.hm3_noMHC."))
    frq = os.environ.get("EUR_FRQ_PREFIX", str(rd / "1000G_Phase3_frq" / "1000G.EUR.QC."))
    ldsc_py = os.environ.get("LDSC_PY", str(md / "tools" / "ldsc" / "ldsc.py"))
    ldsc_python = os.environ.get("LDSC_PYTHON", str(md / "tools" / "ldsc_py3_venv" / "bin" / "python"))

    prefixes = [p.strip() for p in task["prefixes"].split(";") if p.strip()]
    sumstats = Path(task["sumstats_path"])
    out_prefix = Path(task["out_prefix"])
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    if Path(str(out_prefix) + ".results").exists() and has(Path(str(out_prefix) + ".log"), "Analysis finished"):
        print(f"{task['phenocode']} exists; skip"); return
    if not sumstats.exists():
        raise FileNotFoundError(sumstats)

    cmd = [ldsc_python, ldsc_py, "--h2", str(sumstats),
           "--ref-ld-chr", ",".join([baseline] + prefixes),
           "--w-ld-chr", weight, "--frqfile-chr", frq,
           "--overlap-annot", "--print-coefficients", "--print-delete-vals",
           "--out", str(out_prefix)]
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
