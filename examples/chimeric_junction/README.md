# Example: junction annotation

Tiny synthetic inputs (modeled on the Illumina DRAGEN doc example) to try
`scripts/annotate_junctions.py` on. `genes.bed` deliberately contains an
opposite-strand overlapping gene (`GENE_ANTISENSE` inside `GENE_A`) so you can
see the both-strands capture and the `gene_*_aware` columns differ.

```bash
python3 scripts/annotate_junctions.py \
  --junction examples/chimeric_junction/Chimeric.out.junction \
  --genes    examples/chimeric_junction/genes.bed \
  --outdir   /tmp/junction_results

# same-strand-only as the primary output:
python3 scripts/annotate_junctions.py \
  --junction examples/chimeric_junction/Chimeric.out.junction \
  --genes    examples/chimeric_junction/genes.bed \
  --outdir   /tmp/junction_results --strand aware

# optional: add a header to the junction file first
scripts/add_junction_header.sh examples/chimeric_junction/Chimeric.out.junction /tmp/with_header.junction
```

Files:

- `Chimeric.out.junction` — 4 synthetic split-read lines (14 cols, no header)
- `genes.bed` — 4 genes, 6-col BED (0-based half-open)
