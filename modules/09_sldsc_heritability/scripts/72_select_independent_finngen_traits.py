#!/usr/bin/env python3
"""Select approximately independent FinnGen traits from the official R12 FIN rg matrix.

The focal S-LDSC results are never used for selection. Traits first pass an
official univariate LDSC h2 Z threshold. Significant genetic-correlation edges
are then defined by BH-adjusted pairwise rg P values and an absolute-rg cutoff.
Traits are greedily clumped in descending h2 Z and effective sample size so
that no retained pair is connected by a qualifying high-rg edge.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
from pathlib import Path


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def as_bool(value):
    return str(value).strip().lower() in {"1", "true", "t", "yes"}


def bh_adjust(values):
    """Benjamini-Hochberg adjustment, preserving the input order."""
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i])
    adjusted = [float("nan")] * n
    running = 1.0
    for rank_index in range(n - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[original_index] * n / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def effective_n(row):
    cases = as_float(row.get("num_cases"))
    controls = as_float(row.get("num_controls"))
    if cases > 0 and controls > 0:
        return 4.0 / (1.0 / cases + 1.0 / controls)
    return as_float(row.get("total_n"))


def threshold_slug(value):
    return f"rg{value:.2f}".rstrip("0").rstrip(".").replace(".", "p")


def sha256(path: Path):
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
    parser.add_argument("--endpoint-manifest", required=True)
    parser.add_argument("--rg-matrix", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--h2-z-min", type=float, default=4.0)
    parser.add_argument("--rg-thresholds", default="0.6,0.7,0.8")
    parser.add_argument("--primary-rg-threshold", type=float, default=0.7)
    parser.add_argument("--rg-fdr", type=float, default=0.05)
    args = parser.parse_args()

    endpoint_path = Path(args.endpoint_manifest).resolve()
    rg_path = Path(args.rg_matrix).resolve()
    out_dir = Path(args.out_dir).resolve()
    endpoints = read_tsv(endpoint_path)
    if not endpoints:
        raise ValueError("Endpoint manifest is empty")
    endpoint_by_code = {row["phenocode"]: row for row in endpoints}
    if len(endpoint_by_code) != len(endpoints):
        raise ValueError("Endpoint manifest contains duplicate phenocodes")
    target = set(endpoint_by_code)

    diagonal = {}
    pairs = {}
    total_rows = 0
    with gzip.open(rg_path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"p1", "p2", "rg", "se", "p", "h2_obs", "h2_obs_se", "CONVERGED"}
        missing_columns = required - set(reader.fieldnames or [])
        if missing_columns:
            raise ValueError(f"Official rg matrix lacks columns: {sorted(missing_columns)}")
        for row in reader:
            total_rows += 1
            p1, p2 = row["p1"], row["p2"]
            if p1 not in target and p2 not in target:
                continue
            if p1 == p2 and p1 in target:
                diagonal[p1] = row
            if p1 in target and p2 in target and p1 != p2:
                key = tuple(sorted((p1, p2)))
                current = pairs.get(key)
                record = {
                    "rg": as_float(row["rg"]),
                    "se": as_float(row["se"]),
                    "p": as_float(row["p"]),
                    "converged": as_bool(row["CONVERGED"]),
                }
                if current is not None:
                    comparable = (current["rg"], current["se"], current["p"], current["converged"])
                    incoming = (record["rg"], record["se"], record["p"], record["converged"])
                    if comparable != incoming:
                        raise ValueError(f"Conflicting duplicate rg rows for {key}")
                pairs[key] = record

    missing_diagonal = sorted(target - set(diagonal))
    if missing_diagonal:
        raise ValueError(f"Official rg matrix is incomplete for {len(missing_diagonal)} endpoints: {missing_diagonal}")

    h2 = {}
    eligible = []
    for code in sorted(target):
        row = diagonal[code]
        estimate = as_float(row["h2_obs"])
        se = as_float(row["h2_obs_se"])
        z = estimate / se if math.isfinite(estimate) and math.isfinite(se) and se > 0 else float("nan")
        h2[code] = {"h2_obs": estimate, "h2_obs_se": se, "h2_z": z,
                    "converged": as_bool(row["CONVERGED"])}
        if h2[code]["converged"] and math.isfinite(z) and z > args.h2_z_min:
            eligible.append(code)

    expected_pairs = {tuple(sorted((eligible[i], eligible[j])))
                      for i in range(len(eligible)) for j in range(i + 1, len(eligible))}
    absent_pairs = sorted(expected_pairs - set(pairs))
    if absent_pairs:
        raise ValueError(f"Official rg matrix lacks {len(absent_pairs)} h2-eligible trait pairs; first: {absent_pairs[:5]}")

    tested_keys = []
    tested_p = []
    for key in sorted(expected_pairs):
        record = pairs[key]
        if record["converged"] and math.isfinite(record["p"]):
            tested_keys.append(key)
            tested_p.append(record["p"])
    adjusted = bh_adjust(tested_p)
    for key, qvalue in zip(tested_keys, adjusted):
        pairs[key]["q"] = qvalue
    for key in expected_pairs - set(tested_keys):
        pairs[key]["q"] = float("nan")

    thresholds = sorted({float(x) for x in args.rg_thresholds.split(",") if x.strip()})
    if args.primary_rg_threshold not in thresholds:
        thresholds.append(args.primary_rg_threshold)
        thresholds.sort()
    summary = []
    all_fields = list(endpoints[0]) + [
        "official_h2_obs", "official_h2_obs_se", "official_h2_z", "passes_h2_qc",
        "rg_threshold", "rg_fdr", "rg_component_id", "component_size",
        "representative_phenocode", "rg_with_representative", "rg_q_with_representative",
        "selected_for_meta", "selection_reason", "effective_n",
    ]

    h2_eligible_rows = []
    for endpoint in endpoints:
        code = endpoint["phenocode"]
        if code not in eligible:
            continue
        h2_eligible_rows.append({
            **endpoint,
            "official_h2_obs": h2[code]["h2_obs"],
            "official_h2_obs_se": h2[code]["h2_obs_se"],
            "official_h2_z": h2[code]["h2_z"],
            "passes_h2_qc": "TRUE",
        })
    h2_fields = list(endpoints[0]) + [
        "official_h2_obs", "official_h2_obs_se", "official_h2_z", "passes_h2_qc",
    ]
    if "task_id" in h2_fields:
        for task_id, row in enumerate(h2_eligible_rows, 1):
            row["task_id"] = task_id
    write_rows(out_dir / "h2_eligible_endpoints.tsv", h2_eligible_rows, h2_fields)

    for threshold in thresholds:
        adjacency = {code: set() for code in eligible}
        edge_count = 0
        for key in sorted(expected_pairs):
            record = pairs[key]
            qvalue = record.get("q", float("nan"))
            if (record["converged"] and math.isfinite(record["rg"]) and
                    math.isfinite(qvalue) and qvalue < args.rg_fdr and abs(record["rg"]) >= threshold):
                left, right = key
                adjacency[left].add(right)
                adjacency[right].add(left)
                edge_count += 1

        priority = sorted(
            eligible,
            key=lambda code: (-h2[code]["h2_z"], -effective_n(endpoint_by_code[code]), code),
        )
        selected_order = []
        pruned_by = {}
        for code in priority:
            if code in pruned_by:
                continue
            selected_order.append(code)
            pruned_by[code] = code
            for neighbor in sorted(adjacency[code]):
                if neighbor not in pruned_by:
                    pruned_by[neighbor] = code
        selected = set(selected_order)
        violating_retained_edges = [
            (left, right) for left in selected for right in adjacency[left]
            if left < right and right in selected
        ]
        if violating_retained_edges:
            raise AssertionError(f"rg clumping retained correlated pairs: {violating_retained_edges[:5]}")

        edge_rows = []
        for left, right in sorted(expected_pairs):
            record = pairs[(left, right)]
            qvalue = record.get("q", float("nan"))
            if not (record["converged"] and math.isfinite(record["rg"]) and
                    math.isfinite(qvalue) and qvalue < args.rg_fdr and abs(record["rg"]) >= threshold):
                continue
            edge_rows.append({
                "phenocode1": left,
                "phenotype1": endpoint_by_code[left].get("phenotype", ""),
                "phenocode2": right,
                "phenotype2": endpoint_by_code[right].get("phenotype", ""),
                "rg": record["rg"],
                "rg_se": record["se"],
                "rg_p": record["p"],
                "rg_bh_q": qvalue,
                "abs_rg_threshold": threshold,
                "rg_fdr": args.rg_fdr,
            })
        edge_fields = list(edge_rows[0]) if edge_rows else [
            "phenocode1", "phenotype1", "phenocode2", "phenotype2", "rg",
            "rg_se", "rg_p", "rg_bh_q", "abs_rg_threshold", "rg_fdr",
        ]

        clumps = {representative: [] for representative in selected_order}
        for code, representative in pruned_by.items():
            clumps[representative].append(code)
        component_info = {}
        for index, representative in enumerate(selected_order, 1):
            members = sorted(clumps[representative])
            component_id = f"{threshold_slug(threshold)}_c{index:03d}"
            for code in members:
                component_info[code] = (component_id, len(members), representative)

        selection_rows = []
        selected_rows = []
        for endpoint in endpoints:
            code = endpoint["phenocode"]
            passes = code in eligible
            if passes:
                component_id, component_size, representative = component_info[code]
                keep = code in selected
                reason = "rg_clump_representative" if keep else f"directly_correlated_with_representative:{representative}"
                if keep:
                    representative_rg = representative_q = ""
                else:
                    record = pairs[tuple(sorted((code, representative)))]
                    representative_rg = record["rg"]
                    representative_q = record.get("q", float("nan"))
            else:
                component_id, component_size, representative, keep = "", 0, "", False
                reason = "failed_official_h2_z_qc"
                representative_rg = representative_q = ""
            extra = {
                "official_h2_obs": h2[code]["h2_obs"],
                "official_h2_obs_se": h2[code]["h2_obs_se"],
                "official_h2_z": h2[code]["h2_z"],
                "passes_h2_qc": str(passes).upper(),
                "rg_threshold": threshold,
                "rg_fdr": args.rg_fdr,
                "rg_component_id": component_id,
                "component_size": component_size,
                "representative_phenocode": representative,
                "rg_with_representative": representative_rg,
                "rg_q_with_representative": representative_q,
                "selected_for_meta": str(keep).upper(),
                "selection_reason": reason,
                "effective_n": effective_n(endpoint),
            }
            selection_rows.append({**endpoint, **extra})
            if keep:
                selected_rows.append({**endpoint, **extra})

        slug = threshold_slug(threshold)
        write_rows(out_dir / f"rg_pruning_edges_{slug}.tsv", edge_rows, edge_fields)
        write_rows(out_dir / f"trait_selection_{slug}.tsv", selection_rows, all_fields)
        write_rows(out_dir / f"selected_endpoints_{slug}.tsv", selected_rows, all_fields)
        summary.append({
            "analysis": "primary" if threshold == args.primary_rg_threshold else "sensitivity",
            "rg_threshold": threshold,
            "rg_fdr": args.rg_fdr,
            "h2_z_min_strict": args.h2_z_min,
            "n_manifest_endpoints": len(endpoints),
            "n_h2_eligible": len(eligible),
            "n_rg_tests": len(tested_keys),
            "n_rg_edges": edge_count,
            "n_rg_clumps": len(selected_order),
            "n_selected_for_meta": len(selected_rows),
        })

    write_rows(out_dir / "trait_selection_summary.tsv", summary, list(summary[0]))
    provenance = [{
        "rg_matrix": str(rg_path),
        "rg_matrix_sha256": sha256(rg_path),
        "rg_matrix_size_bytes": rg_path.stat().st_size,
        "rg_matrix_rows": total_rows,
        "endpoint_manifest": str(endpoint_path),
        "endpoint_manifest_sha256": sha256(endpoint_path),
        "h2_z_min_strict": args.h2_z_min,
        "rg_fdr": args.rg_fdr,
        "primary_rg_threshold": args.primary_rg_threshold,
        "selection_order": "official_h2_z_desc,effective_n_desc,phenocode_asc",
        "edge_definition": "CONVERGED and rg_BH_q<rg_fdr and abs(rg)>=threshold",
        "clumping_definition": "greedy_direct_edge_clumping;retained_traits_have_no_qualifying_pairwise_rg_edge",
    }]
    write_rows(out_dir / "trait_selection_provenance.tsv", provenance, list(provenance[0]))
    for row in summary:
        print("\t".join(f"{key}={value}" for key, value in row.items()))


if __name__ == "__main__":
    main()
