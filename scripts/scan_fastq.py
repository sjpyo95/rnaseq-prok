#!/usr/bin/env python3
"""
Scan a directory for paired-end FASTQ files and print a draft samples.csv.

Usage: python3 scripts/scan_fastq.py <data_dir>

Handles common naming conventions:
  sample_R1.fastq.gz / sample_R2.fastq.gz
  sample_R1_001.fastq.gz / sample_R2_001.fastq.gz  (Illumina bcl2fastq)
  sample_1.fastq.gz / sample_2.fastq.gz
"""
import sys
import re
from pathlib import Path


def find_r1_files(root: Path) -> list[Path]:
    patterns = [
        "*_R1.fastq.gz", "*_R1.fq.gz",
        "*_R1_001.fastq.gz", "*_R1_001.fq.gz",
        "*_1.fastq.gz", "*_1.fq.gz",
    ]
    found = []
    for pat in patterns:
        found.extend(root.rglob(pat))
    return sorted(set(found))


def infer_r2(r1: Path) -> Path:
    name = r1.name
    for src, dst in [("_R1_001.", "_R2_001."), ("_R1.", "_R2."), ("_1.", "_2.")]:
        if src in name:
            return r1.parent / name.replace(src, dst, 1)
    return r1.parent / (name + ".r2_not_found")


def infer_sample_id(r1: Path) -> str:
    name = r1.name
    name = re.sub(r"_R1(_001)?\.(fastq|fq)\.gz$", "", name)
    name = re.sub(r"_1\.(fastq|fq)\.gz$", "", name)
    return name


def main():
    if len(sys.argv) < 2:
        print("Usage: scan_fastq.py <data_dir>", file=sys.stderr)
        sys.exit(1)

    root = Path(sys.argv[1]).expanduser().resolve()
    if not root.exists():
        print(f"ERROR: directory not found: {root}", file=sys.stderr)
        sys.exit(1)

    r1_files = find_r1_files(root)
    if not r1_files:
        print(f"No paired FASTQ files found in {root}", file=sys.stderr)
        print("Expected names like: sample_R1.fastq.gz / sample_R2.fastq.gz", file=sys.stderr)
        sys.exit(1)

    rows = []
    missing_r2 = []
    for r1 in r1_files:
        r2 = infer_r2(r1)
        sample_id = infer_sample_id(r1)
        r2_path = str(r2) if r2.exists() else f"MISSING ({r2.name})"
        rows.append((sample_id, str(r1), r2_path, r2.exists()))
        if not r2.exists():
            missing_r2.append(r2)

    print(f"# Found {len(rows)} sample(s) in {root}")
    print("sample_id\tcondition\treplicate\tstrandedness\tfq1\tfq2")
    for sample_id, fq1, fq2, _ in rows:
        print(f"{sample_id}\tUNKNOWN\t1\tunknown\t{fq1}\t{fq2}")

    if missing_r2:
        print(f"\n# WARNING: R2 not found for {len(missing_r2)} sample(s):", file=sys.stderr)
        for p in missing_r2:
            print(f"#   {p}", file=sys.stderr)


if __name__ == "__main__":
    main()
