#!/usr/bin/env python3
"""Build a baselineLD v2.2 resource with molecular QTL MaxCPP covariates removed."""

import csv
import gzip
from pathlib import Path

from work_root import work_root


DROP_LDSCORE_COLUMNS = {
    "GTEx_eQTL_MaxCPPL2",
    "BLUEPRINT_H3K27acQTL_MaxCPPL2",
    "BLUEPRINT_H3K4me1QTL_MaxCPPL2",
    "BLUEPRINT_DNA_methylation_MaxCPPL2",
}

DROP_ANNOT_COLUMNS = {name.removesuffix("L2") for name in DROP_LDSCORE_COLUMNS}


def module_dir_from_script() -> Path:
    return work_root()


def filter_table(in_path: Path, out_path: Path, drop_columns: set[str]):
    with gzip.open(in_path, "rt", newline="") as in_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing header in {in_path}")
        keep = [name for name in reader.fieldnames if name not in drop_columns]
        missing = sorted(drop_columns - set(reader.fieldnames))
        if missing:
            raise ValueError(f"{in_path} missing expected columns: {missing}")
        with gzip.open(out_path, "wt", compresslevel=1, newline="") as out_handle:
            writer = csv.DictWriter(out_handle, delimiter="\t", fieldnames=keep, lineterminator="\n")
            writer.writeheader()
            for row in reader:
                writer.writerow({name: row[name] for name in keep})


def filter_m_file(in_path: Path, out_path: Path, annotation_names: list[str], keep_names: list[str]):
    values = in_path.read_text().strip().split()
    if len(values) != len(annotation_names):
        raise ValueError(f"{in_path} has {len(values)} values but {len(annotation_names)} annotations")
    lookup = dict(zip(annotation_names, values))
    out_path.write_text("\t".join(lookup[name] for name in keep_names) + "\n")


def main():
    module_dir = module_dir_from_script()
    source_prefix = module_dir / "resources" / "eur_sldsc_ref" / "baselineLD."
    out_dir = module_dir / "resources" / "eur_sldsc_ref" / "baselineLD_noQTL"
    out_dir.mkdir(parents=True, exist_ok=True)

    for chrom in range(1, 23):
        in_ld = Path(f"{source_prefix}{chrom}.l2.ldscore.gz")
        out_ld = out_dir / f"baselineLD_noQTL.{chrom}.l2.ldscore.gz"
        filter_table(in_ld, out_ld, DROP_LDSCORE_COLUMNS)

        with gzip.open(in_ld, "rt") as handle:
            header = handle.readline().strip().split("\t")
        annotation_names = header[3:]
        keep_names = [name for name in annotation_names if name not in DROP_LDSCORE_COLUMNS]

        in_annot = Path(f"{source_prefix}{chrom}.annot.gz")
        out_annot = out_dir / f"baselineLD_noQTL.{chrom}.annot.gz"
        filter_table(in_annot, out_annot, DROP_ANNOT_COLUMNS)

        filter_m_file(Path(f"{source_prefix}{chrom}.l2.M"), out_dir / f"baselineLD_noQTL.{chrom}.l2.M", annotation_names, keep_names)
        filter_m_file(
            Path(f"{source_prefix}{chrom}.l2.M_5_50"),
            out_dir / f"baselineLD_noQTL.{chrom}.l2.M_5_50",
            annotation_names,
            keep_names,
        )
        source_log = Path(f"{source_prefix}{chrom}.log")
        out_log = out_dir / f"baselineLD_noQTL.{chrom}.log"
        if source_log.exists() and not out_log.exists():
            out_log.write_text(
                source_log.read_text(errors="replace")
                + "\n\nConstructed by removing molecular QTL MaxCPP covariates from baselineLD v2.2.\n"
            )
        print(out_ld)


if __name__ == "__main__":
    main()
