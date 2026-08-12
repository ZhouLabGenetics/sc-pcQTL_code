#!/usr/bin/env python3
"""Align a single-annotation LDSC resource to the baseline LD-score SNP order."""

import argparse
import csv
import gzip
import shutil
from pathlib import Path

from work_root import work_root


def module_dir_from_script() -> Path:
    return work_root()


def read_snps(path: Path):
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row["SNP"] for row in reader]


def align_one(chrom: int, baseline_prefix: Path, custom_prefix: Path, annotation_name: str, backup_dir: Path):
    baseline_ld = Path(f"{baseline_prefix}{chrom}.l2.ldscore.gz")
    custom_ld = Path(f"{custom_prefix}{chrom}.l2.ldscore.gz")
    if not baseline_ld.exists():
        raise FileNotFoundError(baseline_ld)
    if not custom_ld.exists():
        raise FileNotFoundError(custom_ld)

    target_snps = read_snps(baseline_ld)
    target = set(target_snps)
    rows = {}
    scanned = 0
    with gzip.open(custom_ld, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "SNP" not in reader.fieldnames:
            raise ValueError(f"Missing SNP column in {custom_ld}")
        ld_cols = [col for col in reader.fieldnames if col not in ("CHR", "SNP", "BP")]
        if len(ld_cols) != 1:
            raise ValueError(f"Expected exactly one LD-score column in {custom_ld}; found {ld_cols}")
        ld_col = ld_cols[0]
        for row in reader:
            scanned += 1
            snp = row["SNP"]
            if snp in target:
                rows[snp] = row

    missing = [snp for snp in target_snps if snp not in rows]
    if missing:
        raise ValueError(f"{custom_ld} missing {len(missing)} baseline SNPs; first missing {missing[:5]}")

    out_header = ["CHR", "SNP", "BP", f"{annotation_name}L2"]
    tmp_path = custom_ld.with_suffix(custom_ld.suffix + ".aligned_tmp")
    with gzip.open(tmp_path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=out_header, lineterminator="\n")
        writer.writeheader()
        for snp in target_snps:
            row = rows[snp]
            writer.writerow({
                "CHR": row["CHR"],
                "SNP": row["SNP"],
                "BP": row["BP"],
                f"{annotation_name}L2": row[ld_col],
            })

    backup_dir.mkdir(parents=True, exist_ok=True)
    backup_path = backup_dir / custom_ld.name
    if not backup_path.exists():
        shutil.move(str(custom_ld), str(backup_path))
    else:
        custom_ld.unlink()
    shutil.move(str(tmp_path), str(custom_ld))
    return len(target_snps), scanned


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotation-name", required=True)
    parser.add_argument("--custom-prefix", required=True)
    parser.add_argument("--baseline-prefix", default=None)
    parser.add_argument("--chrom", default="all")
    args = parser.parse_args()

    module_dir = module_dir_from_script()
    baseline_prefix = Path(args.baseline_prefix) if args.baseline_prefix else (
        module_dir / "resources" / "eur_sldsc_ref" / "baselineLD_noQTL" / "baselineLD_noQTL."
    )
    custom_prefix = Path(args.custom_prefix)
    if not custom_prefix.is_absolute():
        custom_prefix = module_dir / custom_prefix
    backup_dir = custom_prefix.parent / "unaligned_full"

    chroms = range(1, 23) if args.chrom == "all" else [int(args.chrom)]
    for chrom in chroms:
        matched, scanned = align_one(chrom, baseline_prefix, custom_prefix, args.annotation_name, backup_dir)
        print(f"aligned chr{chrom}: matched={matched} scanned={scanned}")


if __name__ == "__main__":
    main()
