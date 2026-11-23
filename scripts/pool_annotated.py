#!/usr/bin/env python3
"""
pool_annotated.py — pool per-run output TSVs of annotate_junctions.py into one file.

Use this when each sample/run was annotated separately (separate --outdir per
run). If you annotate a whole folder at once (`--junction <folder>`), the
outputs are already pooled and you don't need this.

Usage:
    python3 scripts/pool_annotated.py [--which WHICH] [--out pooled.tsv] PATH [PATH ...]

Each PATH is either an output TSV file or an output directory containing it.
WHICH selects the file: annotated_reads (default), junction_summary, disagreements.
Writes to stdout unless --out is given. A single header row is kept; '#' comment
lines and repeated header rows from later inputs are dropped. .gz inputs OK.
"""

import argparse
import gzip
import os
import sys

WHICH_CHOICES = ("annotated_reads", "junction_summary", "disagreements")


def open_text(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, encoding="utf-8", errors="replace")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help="annotated TSV files or annotate_junctions.py output dirs")
    ap.add_argument("--which", default="annotated_reads", choices=WHICH_CHOICES,
                    help="which output to pool (default: annotated_reads)")
    ap.add_argument("--out", help="output file (default: stdout)")
    args = ap.parse_args()

    files = []
    for p in args.paths:
        if os.path.isdir(p):
            p = os.path.join(p, args.which + ".tsv")
        if not os.path.isfile(p):
            sys.exit(f"ERROR: not found: {p}")
        files.append(p)

    out = open(args.out, "w", encoding="utf-8") if args.out else sys.stdout
    header = None
    n_rows = 0
    try:
        for f in files:
            with open_text(f) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    if header is None:
                        header = line
                        out.write(line)
                        continue
                    cols = line.split("\t")
                    hcols = header.split("\t")
                    if cols[:2] == hcols[:2]:  # repeated header row of another input
                        continue
                    out.write(line)
                    n_rows += 1
    finally:
        if args.out:
            out.close()
    print(f"pooled {n_rows:,} data rows from {len(files)} file(s) "
          f"({args.which}) -> {args.out or 'stdout'}", file=sys.stderr)


if __name__ == "__main__":
    main()
