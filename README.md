# bioinfo-random-scripts

> Lot's of little scripts that I use that don't fit neatly elsewhere.

A home for small bioinformatics helper scripts — each one self-contained, with
its own doc and a tiny example so it can be tried immediately.

## Repository layout

| Path | Purpose |
|---|---|
| `scripts/` | runnable scripts (bash / python), one per task |
| `docs/` | usage notes and workflow docs, one per script |
| `examples/` | tiny example inputs (+ expected command) per script |
| `tests/` | tests, where a script has them |
| `data/` | local data — **everything in here is gitignored** except its README |

Conventions:

- Python scripts are **standard-library only** where possible, with `--help`
  via `argparse`; bash scripts use `set -euo pipefail`.
- Scripts read/write paths given as arguments — no hardcoded locations.
- Raw data and results are **never committed**: `.gitignore` broadly excludes
  reads/alignments/calls/junction files, reference & annotation files, tabular
  exports, archives, big binaries, and the `data/`, `results*/`, `output/`,
  `out/`, `work/`, `tmp/` directories. Put local data in `data/`; only small
  synthetic files under `docs/` and `examples/` are tracked.
- Every script gets a row in the table below and a doc in `docs/`.

## Scripts

| Script | What it does | Docs |
|---|---|---|
| `scripts/annotate_junctions.py` | Annotate DRAGEN RNA `Chimeric.out.junction` split-read breakpoints — donor, acceptor, segment 1 and segment 2 — with gene names (both strands captured, same-strand view kept alongside), plus a per-breakpoint summary and a disagreements report | `docs/annotate_junctions.md` |
| `scripts/add_junction_header.sh` | Prepend a descriptive `#` header to a headerless `Chimeric.out.junction` | `docs/annotate_junctions.md` |
| `scripts/list_junction_files.sh` | List all `*.junction` / `*.junction.gz` files under a folder (recursive), or echo a single file | `docs/annotate_junctions.md` |

## Quick start (junction annotation)

```bash
# 1. optionally add a header to the headerless DRAGEN junction file
scripts/add_junction_header.sh Chimeric.out.junction Chimeric.out.junction.with_header

# 2. annotate all four coordinate sets with gene names
python3 scripts/annotate_junctions.py \
  --junction Chimeric.out.junction \
  --genes   genes.bed \        # BED, or a GTF/GFF with --feature gene
  --outdir  results

# try it on the bundled example first:
python3 scripts/annotate_junctions.py \
  --junction examples/chimeric_junction/Chimeric.out.junction \
  --genes    examples/chimeric_junction/genes.bed \
  --outdir   /tmp/junction_results

# or point --junction at a folder of DRAGEN runs (collects *.junction recursively):
python3 scripts/annotate_junctions.py --junction data/ --genes data/genes.bed --outdir data/results
```

## License

MIT — see [LICENSE](LICENSE). Feel free to use, copy and modify; attribution appreciated.
