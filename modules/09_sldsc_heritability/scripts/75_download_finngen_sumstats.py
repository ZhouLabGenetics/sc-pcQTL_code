#!/usr/bin/env python3
"""Download one missing FinnGen summary-stat file with atomic validation."""

import argparse
import csv
import gzip
import os
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("task_id")
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()
    with Path(args.manifest).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    matches = [row for row in rows if row["task_id"] == args.task_id]
    if len(matches) != 1:
        raise ValueError(f"Expected one row for task {args.task_id}; found {len(matches)}")
    row = matches[0]
    target = Path(row["summary_stats_path"])
    target.parent.mkdir(parents=True, exist_ok=True)

    if target.is_file() and target.stat().st_size > 0:
        with gzip.open(target, "rb") as handle:
            handle.read(1)
        print(f"{row['phenocode']} already exists; skip")
        return
    if not row.get("source_url"):
        raise ValueError(f"No source URL for {row['phenocode']}")

    partial = Path(str(target) + ".part")
    cmd = [
        "curl", "--fail", "--location", "--retry", "8", "--retry-delay", "10",
        "--continue-at", "-", "--output", str(partial), row["source_url"],
    ]
    subprocess.run(cmd, check=True)
    with gzip.open(partial, "rb") as handle:
        while handle.read(8 * 1024 * 1024):
            pass
    os.replace(partial, target)
    print(f"downloaded {row['phenocode']} bytes={target.stat().st_size}")


if __name__ == "__main__":
    main()
