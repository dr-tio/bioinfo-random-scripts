# data/

Local data goes here. **Everything in this folder is gitignored** (the `data/*`
rule in the root `.gitignore`) — only this README is tracked.

Use it for raw inputs and run outputs, e.g.:

- `data/Chimeric.out.junction` (DRAGEN RNA output, one per run — subfolders are
  fine, the annotator collects `*.junction` recursively from a folder)
- `data/genes.bed` / `data/gencode.gtf` (gene annotation)
- `data/results/` (annotator output)

```bash
python3 scripts/annotate_junctions.py \
  --junction data/            \
  --genes    data/genes.bed   \
  --outdir   data/results
```

Never commit real sample data — if a file must be shared, put a small synthetic
version in `examples/` instead.
