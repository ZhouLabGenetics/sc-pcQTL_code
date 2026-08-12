#!/usr/bin/env python3
"""Build a OneK1K-aligned SMR .ma file for FinnGen R12 endpoint 3019198.

FinnGen laboratory-value summary statistics are reported on GRCh38, whereas
the retained OneK1K GIMAP BESD and PLINK LD reference use GRCh37/hg19. This
script reverse-maps FinnGen chromosome 7 positions through the OneK1K
hg19-to-hg38 position map, then aligns every retained variant to the allele
order in the OneK1K PLINK .bim file.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import TextIO


REQUIRED_FINNGEN_COLUMNS = {
    "#chrom", "pos", "ref", "alt", "pval", "beta", "sebeta", "af_alt"
}
COMPLEMENT = str.maketrans("ACGT", "TGCA")
PALINDROMIC = {frozenset(("A", "T")), frozenset(("C", "G"))}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--finngen-gwas", required=True, type=Path,
                        help="FinnGen R12 endpoint 3019198 raw summary .gz.")
    parser.add_argument("--posmap", required=True, type=Path,
                        help="OneK1K chr7 hg19-to-hg38 two-column position map.")
    parser.add_argument("--ld-bim", required=True, type=Path,
                        help="OneK1K chromosome 7 PLINK .bim file.")
    parser.add_argument("--output-ma", required=True, type=Path)
    parser.add_argument("--qc-json", required=True, type=Path)
    parser.add_argument("--chromosome", default="7")
    parser.add_argument("--sample-size", type=int, default=183481,
                        help="FinnGen endpoint 3019198 analyzed sample size.")
    parser.add_argument("--endpoint", default="3019198")
    return parser.parse_args()


def require_nonempty(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty {label}: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_chromosome(value: str) -> str:
    return value.strip().removeprefix("chr")


def load_reverse_posmap(path: Path, chromosome: str) -> tuple[dict[int, int], set[int], Counter[str]]:
    candidates: dict[int, set[int]] = defaultdict(set)
    qc: Counter[str] = Counter()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                qc["posmap_malformed_rows"] += 1
                continue
            old_id, new_pos = fields[:2]
            old_parts = old_id.removeprefix("chr").split(":", 1)
            if len(old_parts) != 2 or normalize_chromosome(old_parts[0]) != chromosome:
                continue
            try:
                candidates[int(new_pos)].add(int(old_parts[1]))
            except ValueError:
                qc["posmap_malformed_rows"] += 1
                continue
            qc["posmap_rows_for_chromosome"] += 1

    unique = {new: next(iter(old)) for new, old in candidates.items() if len(old) == 1}
    ambiguous = {new for new, old in candidates.items() if len(old) > 1}
    qc["posmap_unique_hg38_positions"] = len(unique)
    qc["posmap_ambiguous_hg38_positions"] = len(ambiguous)
    return unique, ambiguous, qc


def load_bim(path: Path, chromosome: str) -> tuple[dict[int, list[tuple[str, str, str]]], Counter[str]]:
    by_position: dict[int, list[tuple[str, str, str]]] = defaultdict(list)
    qc: Counter[str] = Counter()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 6:
                qc["bim_malformed_rows"] += 1
                continue
            chrom, snp, _, pos, allele1, allele2 = fields[:6]
            if normalize_chromosome(chrom) != chromosome:
                continue
            try:
                position = int(pos)
            except ValueError:
                qc["bim_malformed_rows"] += 1
                continue
            by_position[position].append((snp, allele1.upper(), allele2.upper()))
            qc["bim_rows_for_chromosome"] += 1
    qc["bim_unique_positions"] = len(by_position)
    qc["bim_multiallelic_positions"] = sum(len(records) > 1 for records in by_position.values())
    return by_position, qc


def complement(allele: str) -> str | None:
    if len(allele) != 1 or allele not in "ACGT":
        return None
    return allele.translate(COMPLEMENT)


def align_to_bim(
    ref: str,
    alt: str,
    records: list[tuple[str, str, str]],
) -> tuple[tuple[str, str, str, bool] | None, str]:
    """Return (SNP, BIM A1, BIM A2, beta_flip) and an alignment label."""
    ref = ref.upper()
    alt = alt.upper()
    if frozenset((ref, alt)) in PALINDROMIC:
        return None, "palindromic_removed"

    alt_comp = complement(alt)
    ref_comp = complement(ref)
    matches: list[tuple[tuple[str, str, str, bool], str]] = []
    for snp, a1, a2 in records:
        if alt == a1 and ref == a2:
            matches.append(((snp, a1, a2, False), "direct"))
        elif alt == a2 and ref == a1:
            matches.append(((snp, a1, a2, True), "swapped"))
        elif alt_comp == a1 and ref_comp == a2:
            matches.append(((snp, a1, a2, False), "complement_direct"))
        elif alt_comp == a2 and ref_comp == a1:
            matches.append(((snp, a1, a2, True), "complement_swapped"))

    unique = {(record, label) for record, label in matches}
    if not unique:
        return None, "allele_mismatch"
    if len(unique) > 1:
        return None, "multiple_bim_matches"
    return next(iter(unique))


def open_gwas(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def finite_float(value: str) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def build_ma(args: argparse.Namespace) -> dict[str, object]:
    chromosome = normalize_chromosome(args.chromosome)
    reverse_map, ambiguous_map, map_qc = load_reverse_posmap(args.posmap, chromosome)
    bim, bim_qc = load_bim(args.ld_bim, chromosome)
    counts: Counter[str] = Counter()
    seen_snps: set[str] = set()

    args.output_ma.parent.mkdir(parents=True, exist_ok=True)
    with open_gwas(args.finngen_gwas) as source, args.output_ma.open("w", newline="", encoding="utf-8") as target:
        reader = csv.DictReader(source, delimiter="\t")
        missing = REQUIRED_FINNGEN_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"FinnGen summary is missing columns: {sorted(missing)}")
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(("SNP", "A1", "A2", "freq", "b", "se", "p", "n"))

        chromosome_started = False
        for row in reader:
            counts["finngen_rows_scanned"] += 1
            row_chrom = normalize_chromosome(row["#chrom"])
            if row_chrom != chromosome:
                if chromosome_started:
                    break
                continue
            chromosome_started = True
            counts["finngen_rows_on_chromosome"] += 1

            try:
                hg38_pos = int(row["pos"])
            except ValueError:
                counts["invalid_position"] += 1
                continue
            if hg38_pos in ambiguous_map:
                counts["ambiguous_liftover_position"] += 1
                continue
            hg19_pos = reverse_map.get(hg38_pos)
            if hg19_pos is None:
                counts["not_in_liftover_map"] += 1
                continue
            counts["liftover_mapped"] += 1

            bim_records = bim.get(hg19_pos)
            if not bim_records:
                counts["not_in_ld_bim"] += 1
                continue
            counts["ld_position_matched"] += 1

            aligned, label = align_to_bim(row["ref"], row["alt"], bim_records)
            counts[label] += 1
            if aligned is None:
                continue
            snp, allele1, allele2, beta_flip = aligned
            if snp in seen_snps:
                counts["duplicate_snp_removed"] += 1
                continue

            beta = finite_float(row["beta"])
            se = finite_float(row["sebeta"])
            pvalue = finite_float(row["pval"])
            af_alt = finite_float(row["af_alt"])
            if beta is None or se is None or pvalue is None or af_alt is None:
                counts["nonfinite_statistic"] += 1
                continue
            if se <= 0 or not 0 < af_alt < 1 or not 0 <= pvalue <= 1:
                counts["invalid_statistic_range"] += 1
                continue
            if pvalue == 0:
                pvalue = 1e-300
                counts["zero_pvalue_clamped"] += 1

            if beta_flip:
                beta = -beta
                af_alt = 1 - af_alt
                counts["effect_flipped_to_bim_a1"] += 1
            else:
                counts["effect_already_bim_a1"] += 1
            writer.writerow((
                snp, allele1, allele2, f"{af_alt:.10g}", f"{beta:.10g}",
                f"{se:.10g}", f"{pvalue:.10g}", args.sample_size,
            ))
            seen_snps.add(snp)
            counts["variants_written"] += 1

    if counts["variants_written"] == 0:
        raise RuntimeError("No allele-compatible FinnGen variants were written")

    return {
        "analysis": "FinnGen R12 endpoint 3019198 to OneK1K-aligned SMR .ma",
        "endpoint": args.endpoint,
        "trait_type": "quantitative_inverse_rank_normalized",
        "sample_size": args.sample_size,
        "source_build": "GRCh38",
        "target_build": "GRCh37_hg19",
        "effect_allele_source": "FinnGen alt",
        "output_effect_allele": "OneK1K PLINK BIM A1",
        "palindromic_policy": "remove A/T and C/G SNPs",
        "inputs": {
            "finngen_gwas": str(args.finngen_gwas.resolve()),
            "finngen_gwas_sha256": sha256(args.finngen_gwas),
            "posmap": str(args.posmap.resolve()),
            "posmap_sha256": sha256(args.posmap),
            "ld_bim": str(args.ld_bim.resolve()),
            "ld_bim_sha256": sha256(args.ld_bim),
        },
        "output": {
            "ma": str(args.output_ma.resolve()),
            "ma_sha256": sha256(args.output_ma),
        },
        "counts": dict(sorted((map_qc + bim_qc + counts).items())),
    }


def main() -> int:
    args = parse_args()
    require_nonempty(args.finngen_gwas, "FinnGen GWAS summary")
    require_nonempty(args.posmap, "OneK1K liftOver position map")
    require_nonempty(args.ld_bim, "OneK1K LD BIM")
    if args.sample_size <= 0:
        raise ValueError("--sample-size must be positive")
    qc = build_ma(args)
    args.qc_json.parent.mkdir(parents=True, exist_ok=True)
    args.qc_json.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(qc["counts"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
