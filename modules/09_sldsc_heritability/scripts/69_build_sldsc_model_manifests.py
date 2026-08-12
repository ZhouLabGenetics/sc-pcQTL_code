#!/usr/bin/env python3
"""Build endpoint-level S-LDSC task manifests from a QC-approved trait set."""

import argparse
import csv
import re
from pathlib import Path

from work_root import work_root


def read_rows(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def safe_slug(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_]+", "_", value).strip("_")
    if not slug:
        raise ValueError("Analysis prefix is empty after sanitization")
    return slug


def write_manifest(module_dir: Path, endpoints, analysis_prefix: str, name: str,
                   prefixes, output_subdir: str):
    out = module_dir / "manifests" / f"{analysis_prefix}_{name}_tasks.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["task_id", "phenocode", "phenotype", "category", "sumstats_path", "prefixes", "out_prefix"]
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for task_id, endpoint in enumerate(endpoints, 1):
            phenocode = endpoint["phenocode"]
            writer.writerow({
                "task_id": task_id,
                "phenocode": phenocode,
                "phenotype": endpoint.get("phenotype", ""),
                "category": endpoint.get("category", ""),
                "sumstats_path": module_dir / "resources" / "finngen_r12_sldsc_sumstats" / phenocode / f"finngen_R12_{phenocode}.sumstats.gz",
                "prefixes": ";".join(str(module_dir / prefix) for prefix in prefixes),
                "out_prefix": module_dir / "results" / "heritability_enrichment" / output_subdir / phenocode / f"finngen_R12_{phenocode}_{output_subdir}",
            })
    print(f"{out}\t{len(endpoints)}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("main", "functional", "ngenes"))
    parser.add_argument("--endpoint-manifest", required=True,
                        help="Traits that passed the prespecified univariate h2 QC")
    parser.add_argument("--analysis-prefix", default="prespecified247_h2qc")
    args = parser.parse_args()

    module_dir = work_root()
    endpoints = read_rows(Path(args.endpoint_manifest).resolve())
    if not endpoints:
        raise ValueError("Endpoint manifest is empty")
    if len({row["phenocode"] for row in endpoints}) != len(endpoints):
        raise ValueError("Endpoint manifest contains duplicate phenocodes")
    analysis_prefix = safe_slug(args.analysis_prefix)

    if args.mode == "main":
        qtl = "resources/eur_sldsc_custom_annotations/single_maf05"
        write_manifest(module_dir, endpoints, analysis_prefix, "joint", [f"{qtl}/QTLsig_pcQTL_gt0/qtlsig_pcQTL.", f"{qtl}/QTLsig_eQTL_gt0/qtlsig_eQTL."], f"eur_sldsc_qtlsig_{analysis_prefix}_joint")
        write_manifest(module_dir, endpoints, analysis_prefix, "marg_pc", [f"{qtl}/QTLsig_pcQTL_gt0/qtlsig_pcQTL."], f"eur_sldsc_qtlsig_{analysis_prefix}_marg_pc")
        write_manifest(module_dir, endpoints, analysis_prefix, "marg_eq", [f"{qtl}/QTLsig_eQTL_gt0/qtlsig_eQTL."], f"eur_sldsc_qtlsig_{analysis_prefix}_marg_eq")
    elif args.mode == "functional":
        base = "resources/eur_sldsc_custom_annotations/single_maf05_funcinteract"
        write_manifest(module_dir, endpoints, analysis_prefix, "funcinteract", [f"{base}/pcQTL_prom/fi_pcQTL_prom.", f"{base}/pcQTL_enh/fi_pcQTL_enh.", f"{base}/eQTL_prom/fi_eQTL_prom.", f"{base}/eQTL_enh/fi_eQTL_enh."], f"eur_sldsc_qtlsig_{analysis_prefix}_funcinteract")
    else:
        base = "resources/eur_sldsc_custom_annotations/single_maf05_ngenes"
        write_manifest(module_dir, endpoints, analysis_prefix, "ngenes", [f"{base}/ng_small/ng_small.", f"{base}/ng_med/ng_med.", f"{base}/ng_large/ng_large."], f"eur_sldsc_qtlsig_{analysis_prefix}_ngenes")


if __name__ == "__main__":
    main()
