#!/usr/bin/env python3
"""Run the two manuscript-targeted GIMAP SMR/HEIDI analyses.

This runner deliberately accepts already harmonized SMR inputs. It never infers
or substitutes a GWAS source: ``--finngen-gwas-ma`` must be the staged FinnGen
R12 endpoint 3019198 summary used by the manuscript.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
from pathlib import Path
from typing import Any


GENES = ("GIMAP8", "GIMAP7", "GIMAP4", "GIMAP6", "GIMAP2", "GIMAP1", "GIMAP5")
MA_COLUMNS = ("SNP", "A1", "A2", "freq", "b", "se", "p", "n")
HEIDI_EQTL_P_CUTOFF = 0.01
HEIDI_MIN_M = 3
HEIDI_P_CUTOFF = 0.01
SMR_P_CUTOFF = 0.05 / len(GENES)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cluster-pc3-ma", required=True, type=Path)
    parser.add_argument("--finngen-gwas-ma", required=True, type=Path)
    parser.add_argument("--finngen-ma-qc", required=True, type=Path,
                        help="QC JSON emitted by prepare_finngen3019198_ma.py.")
    parser.add_argument("--eqtl-besd", required=True, type=Path,
                        help="Prefix of the seven-gene CD8-naive-T eQTL BESD.")
    parser.add_argument("--ld-prefix", required=True, type=Path,
                        help="Prefix of the OneK1K chromosome 7 PLINK LD files.")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--smr-bin", type=Path,
                        default=Path(os.environ.get("SMR_BIN", "smr")))
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def require_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty {label}: {path}")


def validate_ma(path: Path, label: str) -> None:
    require_file(path, label)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header in {label}: {path}")
        missing = set(MA_COLUMNS) - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing {label} columns {sorted(missing)}: {path}")
        first = next(reader, None)
        if first is None:
            raise ValueError(f"No variants in {label}: {path}")
        for field in ("freq", "b", "se", "p", "n"):
            value = float(first[field])
            if not math.isfinite(value):
                raise ValueError(f"Non-finite {field} in {label}: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_finngen_ma_qc(args: argparse.Namespace) -> dict[str, Any]:
    require_file(args.finngen_ma_qc, "FinnGen 3019198 .ma QC JSON")
    qc = json.loads(args.finngen_ma_qc.read_text(encoding="utf-8"))
    expected = {
        "endpoint": "3019198",
        "sample_size": 183481,
        "source_build": "GRCh38",
        "target_build": "GRCh37_hg19",
    }
    for field, value in expected.items():
        if qc.get(field) != value:
            raise ValueError(
                f"FinnGen .ma QC field {field!r} is {qc.get(field)!r}, expected {value!r}"
            )
    recorded_hash = qc.get("output", {}).get("ma_sha256")
    observed_hash = sha256(args.finngen_gwas_ma)
    if recorded_hash != observed_hash:
        raise ValueError(
            "FinnGen .ma does not match its QC record: "
            f"recorded sha256={recorded_hash!r}, observed={observed_hash!r}"
        )
    if int(qc.get("counts", {}).get("variants_written", 0)) <= 0:
        raise ValueError("FinnGen .ma QC reports no written variants")
    return qc


def validate_inputs(args: argparse.Namespace) -> dict[str, Any]:
    validate_ma(args.cluster_pc3_ma, "cluster-PC3 SMR .ma")
    validate_ma(args.finngen_gwas_ma, "FinnGen R12 endpoint 3019198 SMR .ma")
    finngen_qc = validate_finngen_ma_qc(args)
    for suffix in (".besd", ".epi", ".esi"):
        require_file(Path(f"{args.eqtl_besd}{suffix}"), f"eQTL BESD {suffix}")
    for suffix in (".bed", ".bim", ".fam"):
        require_file(Path(f"{args.ld_prefix}{suffix}"), f"OneK1K chr7 LD {suffix}")
    if not args.dry_run:
        if args.smr_bin.is_absolute():
            require_file(args.smr_bin, "SMR executable")
        elif not any((directory / args.smr_bin).is_file() for directory in map(Path, os.environ.get("PATH", "").split(os.pathsep))):
            raise FileNotFoundError(f"SMR executable not found on PATH: {args.smr_bin}")
    return finngen_qc


def as_float(value: str | None) -> float | None:
    if value is None or value.strip() in {"", "NA", "NaN", "nan", "None"}:
        return None
    parsed = float(value)
    return parsed if math.isfinite(parsed) else None


def bh_adjust(values: list[float | None]) -> list[float | None]:
    ranked = sorted(((value, index) for index, value in enumerate(values) if value is not None))
    adjusted: list[float | None] = [None] * len(values)
    running = 1.0
    for rank in range(len(ranked), 0, -1):
        value, index = ranked[rank - 1]
        running = min(running, value * len(ranked) / rank)
        adjusted[index] = running
    return adjusted


def parse_best_smr(path: Path) -> dict[str, str]:
    if not path.is_file() or path.stat().st_size == 0:
        return {}
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    usable = [(as_float(row.get("p_SMR")), row) for row in rows]
    usable = [(pvalue, row) for pvalue, row in usable if pvalue is not None]
    return min(usable, key=lambda item: item[0])[1] if usable else {}


def status_for(p_smr: float | None, p_heidi: float | None) -> str:
    if p_smr is None:
        return "no_result"
    if p_smr < SMR_P_CUTOFF and p_heidi is not None and p_heidi > HEIDI_P_CUTOFF:
        return "smr_bonferroni_heidi_pass"
    if p_smr < SMR_P_CUTOFF and p_heidi is not None and p_heidi <= HEIDI_P_CUTOFF:
        return "smr_bonferroni_heidi_inconsistent"
    return "weak_or_unclear"


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def run_outcome(
    args: argparse.Namespace,
    outcome: str,
    ma_path: Path,
    output_name: str,
) -> list[dict[str, Any]]:
    run_dir = args.out_dir / "smr_runs" / outcome
    probe_dir = args.out_dir / "probe_lists"
    log_dir = args.out_dir / "logs" / outcome
    run_dir.mkdir(parents=True, exist_ok=True)
    probe_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    for gene in GENES:
        probe = probe_dir / f"{gene}.txt"
        probe.write_text(f"{gene}\n", encoding="utf-8")
        out_prefix = run_dir / gene
        command = [
            str(args.smr_bin),
            "--bfile", str(args.ld_prefix),
            "--gwas-summary", str(ma_path),
            "--beqtl-summary", str(args.eqtl_besd),
            "--extract-probe", str(probe),
            "--peqtl-smr", "0.05",
            "--peqtl-heidi", str(HEIDI_EQTL_P_CUTOFF),
            "--heidi-min-m", str(HEIDI_MIN_M),
            "--disable-freq-ck",
            "--smr",
            "--out", str(out_prefix),
        ]
        if args.dry_run:
            print("DRY-RUN", " ".join(command))
            result: dict[str, str] = {}
        else:
            with (log_dir / f"{gene}.log").open("w", encoding="utf-8") as log:
                subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
            result = parse_best_smr(Path(f"{out_prefix}.smr"))

        p_smr = as_float(result.get("p_SMR"))
        p_heidi = as_float(result.get("p_HEIDI"))
        rows.append({
            "outcome": outcome,
            "celltype": "cd8_nc",
            "cluster": "GIMAP seven-gene cluster",
            "phenotype": "PC3" if outcome == "cluster_pc3" else "FinnGen R12 endpoint 3019198",
            "gene": gene,
            "top_snp": result.get("topSNP", ""),
            "b_smr": result.get("b_SMR", ""),
            "se_smr": result.get("se_SMR", ""),
            "p_smr": result.get("p_SMR", ""),
            "p_heidi": result.get("p_HEIDI", ""),
            "nsnp_heidi": result.get("nsnp_HEIDI", ""),
            "p_gwas_top_snp": result.get("p_GWAS", ""),
            "p_eqtl_top_snp": result.get("p_eQTL", ""),
            "heidi_eqtl_p_cutoff": HEIDI_EQTL_P_CUTOFF,
            "heidi_min_m": HEIDI_MIN_M,
            "bonferroni_cutoff": SMR_P_CUTOFF,
            "heidi_p_cutoff": HEIDI_P_CUTOFF,
            "status": status_for(p_smr, p_heidi),
        })

    q_values = bh_adjust([as_float(row["p_smr"]) for row in rows])
    for row, q_value in zip(rows, q_values):
        row["fdr_within_7_genes"] = "" if q_value is None else q_value

    fields = [
        "outcome", "celltype", "cluster", "phenotype", "gene", "top_snp",
        "b_smr", "se_smr", "p_smr", "fdr_within_7_genes", "p_heidi",
        "nsnp_heidi", "p_gwas_top_snp", "p_eqtl_top_snp",
        "heidi_eqtl_p_cutoff", "heidi_min_m",
        "bonferroni_cutoff", "heidi_p_cutoff", "status",
    ]
    write_tsv(args.out_dir / output_name, rows, fields)
    return rows


def file_record(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": sha256(path),
    }


def smr_binary_banner(path: Path) -> dict[str, Any]:
    """Record the executable's own version banner separately from archive naming."""
    result = subprocess.run(
        [str(path), "--help"], capture_output=True, text=True, check=False
    )
    banner = "\n".join(part for part in (result.stdout, result.stderr) if part)
    version_match = re.search(r"\* Version\s+([^\s]+)", banner)
    build_match = re.search(r"\* Build at\s+(.+)", banner)
    return {
        "version": version_match.group(1) if version_match else "unparsed",
        "build": build_match.group(1).strip() if build_match else "unparsed",
        "help_exit_code": result.returncode,
    }


def main() -> int:
    args = parse_args()
    finngen_qc = validate_inputs(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    cluster_rows = run_outcome(
        args, "cluster_pc3", args.cluster_pc3_ma, "gimap_cluster_pcqtl_smr_summary.tsv"
    )
    gwas_rows = run_outcome(
        args, "finngen_3019198", args.finngen_gwas_ma, "gimap_gwas_smr_summary.tsv"
    )
    write_tsv(
        args.out_dir / "gimap_targeted_smr_summary.tsv",
        cluster_rows + gwas_rows,
        list((cluster_rows or gwas_rows)[0].keys()),
    )

    provenance = {
        "analysis": "targeted GIMAP SMR/HEIDI",
        "smr_release_archive_label": "1.4.0",
        "smr_binary_banner": smr_binary_banner(args.smr_bin) if args.smr_bin.is_file() else None,
        "genes": list(GENES),
        "peqtl_smr": 0.05,
        "peqtl_heidi": HEIDI_EQTL_P_CUTOFF,
        "heidi_min_m": HEIDI_MIN_M,
        "smr_bonferroni_cutoff": SMR_P_CUTOFF,
        "heidi_p_cutoff": HEIDI_P_CUTOFF,
        "frequency_check_disabled": True,
        "dry_run": args.dry_run,
        "finngen_ma_harmonization": finngen_qc,
        "inputs": {
            "cluster_pc3_ma": file_record(args.cluster_pc3_ma),
            "finngen_r12_endpoint_3019198_ma": file_record(args.finngen_gwas_ma),
            "finngen_r12_endpoint_3019198_ma_qc": file_record(args.finngen_ma_qc),
            "eqtl_besd": {
                suffix: file_record(Path(f"{args.eqtl_besd}{suffix}"))
                for suffix in (".besd", ".epi", ".esi")
            },
            "onek1k_chr7_ld": {
                suffix: file_record(Path(f"{args.ld_prefix}{suffix}"))
                for suffix in (".bed", ".bim", ".fam")
            },
            "smr_binary": file_record(args.smr_bin) if args.smr_bin.is_file() else str(args.smr_bin),
        },
    }
    (args.out_dir / "run_provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
