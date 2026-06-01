#!/usr/bin/env python3
"""Validate the sample sheet for the BSF population-genomics pipeline.

Checks performed:
  1. Required columns are present: sample, fastq_1, fastq_2
  2. Sample names are non-empty and unique
  3. Every FASTQ path listed actually exists on disk

Usage:
  python validate_samplesheet.py samplesheet.csv

Exit code 0 = valid, 1 = at least one problem found
(so a pipeline / Nextflow can stop early on a bad sheet).
"""
import csv
import sys
from pathlib import Path

REQUIRED = ["sample", "fastq_1", "fastq_2"]


def main(path):
    errors = []
    seen = set()

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)

        # 1. required columns present?
        missing = [c for c in REQUIRED if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit(f"ERROR: missing column(s): {', '.join(missing)}")

        # 2 & 3. check every row (line 2 = first data row, after the header)
        for i, row in enumerate(reader, start=2):
            name = (row["sample"] or "").strip()
            if not name:
                errors.append(f"line {i}: empty sample name")
            elif name in seen:
                errors.append(f"line {i}: duplicate sample '{name}'")
            else:
                seen.add(name)

            for col in ("fastq_1", "fastq_2"):
                f = (row[col] or "").strip()
                if not f:
                    errors.append(f"line {i}: empty {col}")
                elif not Path(f).is_file():
                    errors.append(f"line {i}: file not found -> {f}")

    if errors:
        print("\n".join(errors))
        sys.exit(f"Sample sheet INVALID ({len(errors)} problem(s)).")

    print(f"Sample sheet OK: {len(seen)} samples, all FASTQ files found.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python validate_samplesheet.py samplesheet.csv")
    main(sys.argv[1])
