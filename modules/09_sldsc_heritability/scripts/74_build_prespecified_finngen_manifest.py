#!/usr/bin/env python3
"""Convert the prespecified FinnGen phenotype table into an S-LDSC manifest.

No FinnGen fine-mapping or colocalization field is consulted when defining the
trait universe.
"""

import argparse
import csv
import hashlib
import os
from pathlib import Path


def as_bool(value):
    return str(value).strip().lower() in {"1", "true", "t", "yes"}


def checksum(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_rows(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-manifest", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--raw-root", default=os.environ.get("FINNGEN_ROOT", "work/finngen"))
    parser.add_argument("--expected-count", type=int, default=247)
    args = parser.parse_args()

    source = Path(args.source_manifest).resolve()
    out = Path(args.out).resolve()
    raw_root = Path(args.raw_root).resolve()
    with source.open(newline="") as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    if "in_finngen_main_endpoints" in (source_rows[0] if source_rows else {}):
        source_rows = [row for row in source_rows if as_bool(row["in_finngen_main_endpoints"])]
    if len(source_rows) != args.expected_count:
        raise ValueError(f"Expected {args.expected_count} prespecified traits; found {len(source_rows)}")
    if len({row["phenocode"] for row in source_rows}) != len(source_rows):
        raise ValueError("Prespecified phenotype table contains duplicate phenocodes")

    rows = []
    missing = []
    for task_id, source_row in enumerate(source_rows, 1):
        code = source_row["phenocode"]
        cases = int(float(source_row.get("num_cases") or 0))
        controls = int(float(source_row.get("num_controls") or 0))
        raw = raw_root / f"finngen_R12_{code}.gz"
        row = {
            "task_id": task_id,
            "phenocode": code,
            "phenotype": source_row.get("phenotype", ""),
            "category": source_row.get("category", ""),
            "num_cases": cases,
            "num_controls": controls,
            "total_n": cases + controls,
            "summary_stats_path": str(raw),
            "source_url": source_row.get("path_https", ""),
            "raw_summary_exists": str(raw.is_file() and raw.stat().st_size > 0).upper(),
        }
        rows.append(row)
        if row["raw_summary_exists"] != "TRUE":
            missing.append(row.copy())

    fields = list(rows[0])
    write_rows(out, rows, fields)
    for index, row in enumerate(missing, 1):
        row["task_id"] = index
    write_rows(out.parent / "missing_raw_finngen_endpoints.tsv", missing, fields)
    provenance = [{
        "source_manifest": str(source),
        "source_manifest_sha256": checksum(source),
        "selection_rule": "in_finngen_main_endpoints_TRUE",
        "gwas_finemapping_used": "FALSE",
        "n_prespecified": len(rows),
        "n_raw_present": len(rows) - len(missing),
        "n_raw_missing": len(missing),
    }]
    write_rows(out.parent / "prespecified247_manifest_provenance.tsv", provenance, list(provenance[0]))
    print(f"prespecified={len(rows)} raw_present={len(rows) - len(missing)} raw_missing={len(missing)}")


if __name__ == "__main__":
    main()
