#!/usr/bin/env python3
"""Read and validate the primary sc-pcQTL cell-type analysis set."""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path


MIN_CELLTYPE_CELLS = 10_000
DEFAULT_MANIFEST = Path(__file__).resolve().parents[1] / "config/celltype_eligibility.tsv"


def load_eligibility(path: str | Path | None = None) -> list[dict[str, str]]:
    manifest = Path(path or os.environ.get("SC_PCQTL_CELLTYPE_MANIFEST", DEFAULT_MANIFEST))
    with manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"celltype", "eqtl_celltype", "display_label", "n_cells", "include_primary"}
    if not rows or required.difference(rows[0]):
        raise ValueError(f"Invalid cell-type eligibility manifest: {manifest}")
    for row in rows:
        n_cells = int(row["n_cells"])
        include = row["include_primary"].upper() == "TRUE"
        if include != (n_cells >= MIN_CELLTYPE_CELLS):
            raise ValueError(
                f"Eligibility mismatch for {row['celltype']}: n_cells={n_cells}, include_primary={include}"
            )
    return rows


def primary_celltypes(path: str | Path | None = None) -> list[str]:
    return [row["celltype"] for row in load_eligibility(path) if row["include_primary"].upper() == "TRUE"]


def primary_eqtl_celltypes(path: str | Path | None = None) -> list[str]:
    return [row["eqtl_celltype"] for row in load_eligibility(path) if row["include_primary"].upper() == "TRUE"]


def excluded_celltypes(path: str | Path | None = None) -> list[str]:
    return [row["celltype"] for row in load_eligibility(path) if row["include_primary"].upper() != "TRUE"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=None)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--list-primary", action="store_true")
    action.add_argument("--list-excluded", action="store_true")
    action.add_argument("--list-primary-mapping", action="store_true")
    action.add_argument("--require-primary", metavar="CELLTYPE")
    args = parser.parse_args()

    rows = load_eligibility(args.manifest)
    primary = [row for row in rows if row["include_primary"].upper() == "TRUE"]
    if args.list_primary:
        print("\n".join(row["celltype"] for row in primary))
        return 0
    if args.list_excluded:
        print("\n".join(row["celltype"] for row in rows if row["include_primary"].upper() != "TRUE"))
        return 0
    if args.list_primary_mapping:
        print("\n".join(f"{row['celltype']}\t{row['eqtl_celltype']}" for row in primary))
        return 0

    match = next((row for row in rows if row["celltype"] == args.require_primary), None)
    if match is None:
        print(f"Cell type is absent from the eligibility manifest: {args.require_primary}", file=sys.stderr)
        return 2
    if match["include_primary"].upper() != "TRUE":
        reason = match.get("exclusion_reason", "excluded from primary analysis")
        print(f"Cell type {args.require_primary} is excluded from primary analysis: {reason}", file=sys.stderr)
        return 2
    print(args.require_primary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
