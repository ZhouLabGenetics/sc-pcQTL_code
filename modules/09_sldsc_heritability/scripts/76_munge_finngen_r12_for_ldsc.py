#!/usr/bin/env python3
"""Pre-clean one FinnGen R12 summary-stat file and run LDSC munging."""

import argparse
import csv
import gzip
import os
import subprocess
from pathlib import Path

from work_root import work_root


VALID_ALLELES = {"A", "C", "G", "T"}


def read_rows(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def valid_existing(prefix: Path):
    output = Path(str(prefix) + ".sumstats.gz")
    log = Path(str(prefix) + ".log")
    if not output.is_file() or not log.is_file() or "Conversion finished" not in log.read_text(errors="replace"):
        return False
    with gzip.open(output, "rb") as handle:
        handle.read(1)
    return True


def clean_sumstats(row, out_path: Path, qc_path: Path):
    total = kept = bad = multi = no_rsid = 0
    with gzip.open(row["summary_stats_path"], "rt", newline="") as in_handle, gzip.open(out_path, "wt", newline="") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        writer = csv.DictWriter(out_handle, delimiter="\t", fieldnames=["SNP", "A1", "A2", "P", "BETA", "FRQ", "N"])
        writer.writeheader()
        for record in reader:
            total += 1
            rsids = (record.get("rsids") or "").strip()
            if not rsids:
                no_rsid += 1
                continue
            if "," in rsids:
                multi += 1
                continue
            a1 = (record.get("alt") or "").upper()
            a2 = (record.get("ref") or "").upper()
            pvalue = record.get("pval") or ""
            beta = record.get("beta") or ""
            frequency = record.get("af_alt") or ""
            if a1 not in VALID_ALLELES or a2 not in VALID_ALLELES or not pvalue or not beta or not frequency:
                bad += 1
                continue
            writer.writerow({"SNP": rsids, "A1": a1, "A2": a2, "P": pvalue, "BETA": beta, "FRQ": frequency, "N": row["total_n"]})
            kept += 1
    with qc_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["phenocode", "input_rows", "kept_rows", "no_rsid", "multi_rsid", "bad_rows"])
        writer.writerow([row["phenocode"], total, kept, no_rsid, multi, bad])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("task_id")
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()
    module_dir = work_root()
    rows = read_rows(Path(args.manifest))
    matches = [row for row in rows if row["task_id"] == args.task_id or row["phenocode"] == args.task_id]
    if len(matches) != 1:
        raise ValueError(f"Expected one manifest row for {args.task_id}; found {len(matches)}")
    row = matches[0]

    out_root = module_dir / "resources" / "finngen_r12_sldsc_sumstats" / row["phenocode"]
    qc_root = module_dir / "qc" / "eur_sldsc_munge"
    out_root.mkdir(parents=True, exist_ok=True)
    qc_root.mkdir(parents=True, exist_ok=True)
    out_prefix = out_root / f"finngen_R12_{row['phenocode']}"
    if valid_existing(out_prefix):
        print(f"{row['phenocode']} already munged; skip")
        return

    raw = Path(row["summary_stats_path"])
    if not raw.is_file():
        raise FileNotFoundError(raw)
    hm3 = Path(os.environ.get("HM3_SNPLIST", module_dir / "resources" / "eur_sldsc_ref" / "w_hm3.snplist"))
    munge = Path(os.environ.get("MUNGE_SUMSTATS", module_dir / "tools" / "ldsc" / "munge_sumstats.py"))
    ldsc_python = os.environ.get("LDSC_PYTHON", str(module_dir / "tools" / "ldsc_py3_venv" / "bin" / "python"))
    clean_path = out_root / f"finngen_R12_{row['phenocode']}.preclean.tsv.gz"
    clean_sumstats(row, clean_path, qc_root / f"{row['phenocode']}.preclean_qc.tsv")
    subprocess.run([
        ldsc_python, str(munge), "--sumstats", str(clean_path), "--snp", "SNP",
        "--a1", "A1", "--a2", "A2", "--p", "P", "--signed-sumstats", "BETA,0",
        "--frq", "FRQ", "--N", row["total_n"], "--merge-alleles", str(hm3),
        "--chunksize", os.environ.get("LDSC_MUNGE_CHUNKSIZE", "500000"), "--out", str(out_prefix),
    ], check=True)
    print(Path(str(out_prefix) + ".sumstats.gz"))


if __name__ == "__main__":
    main()
