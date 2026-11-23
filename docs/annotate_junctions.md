# Annotating DRAGEN `Chimeric.out.junction` with gene names

Goal: for every split read in the junction file, annotate **all four** coordinate
sets — donor, acceptor, segment 1, segment 2 — with gene names, capturing genes on
**both strands** (strand-agnostic primary output; same-strand view kept alongside),
then compare them to see where they differ.

Files:

- `scripts/annotate_junctions.py` — the annotator (pure Python stdlib, no installs)
- `scripts/add_junction_header.sh` — adds a `#` header to the headerless junction file
- `examples/chimeric_junction/` — tiny synthetic inputs to try it on

## 0. Inputs

| File | What it is |
|---|---|
| `Chimeric.out.junction` | DRAGEN RNA pipeline output (14 cols, no header — see below). **Or a folder**: all `*.junction` / `*.junction.gz` files under it are collected recursively (same convention as OGM-toolkit `process_smap.collect_smap_files`) |
| gene file | **BED** (`chrom start end [name] [score] [strand]`, 0-based half-open) **or** GTF/GFF (use `--feature gene`) |

Junction columns (1-based in the file):
`1 donor_chr  2 donor_pos  3 donor_strand  4 acceptor_chr  5 acceptor_pos  6 acceptor_strand  7-9 placeholders  10 read  11 seg1_start  12 seg1_cigar  13 seg2_start  14 seg2_cigar`

### Adding a header to the junction file (optional)

`Chimeric.out.junction` ships headerless. To make it self-describing:

```bash
scripts/add_junction_header.sh Chimeric.out.junction Chimeric.out.junction.with_header
```

The annotator skips `#` lines, so original and headed files both work.

### Gene file format

BED example (start is 0-based, end 1-based/half-open; name and strand optional):

```
chr19  580000  581000  GENE_A  0  +
```

GTF/GFF works directly (reads the `gene_name` attribute) — e.g. GENCODE or NCBI
annotation GTFs: `--genes gencode.v44.annotation.gtf --feature gene`.

**NCBI gene table** (the `Accession / Begin / End / Chromosome / Orientation /
Name / Symbol / Gene ID / Gene Type / ...` layout — tab- or pipe-delimited, header
optional) is supported natively: `Symbol` is used as the gene name, `Begin`
(1-based) is converted to 0-based, `Orientation` (plus/minus) becomes the strand,
and the plain `Chromosome` numbering (`1`) is auto-reconciled with `chr1` or
`NC_000001.11` in the junction file. Try it: `--genes examples/chimeric_junction/genes_ncbi.tsv`.

If your table is some other custom layout (e.g. 1-based `gene  chrom  start  end  strand`),
convert to BED first, e.g.:

```bash
awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ {print $2, $3-1, $4, $1, 0, $5}' raw.tsv > genes.bed
```

Gotchas:

- The **assembly must match the DRAGEN run's reference** (GRCh38 vs hg19). The
  script prints the chromosome-match % at startup, so a mismatch is obvious.
- Chromosome naming (`chr19` vs `19` vs `NC_000019.10`) is auto-reconciled
  (`--chrom-mode` to force). Extra scaffolds in the gene file (HLA, decoys,
  `chrUn_*`, alts) are harmless — they simply never match.
- `.gz` inputs work.

## 1. Quick sanity check on inputs (30 s)

```bash
cut -f1,4 Chimeric.out.junction | tr '\t' '\n' | sort | uniq -c | sort -rn | head   # junction chroms
cut -f1 genes.bed | sort -u | head                                                     # gene-file chroms
wc -l Chimeric.out.junction
```

## 2. Run

```bash
python3 scripts/annotate_junctions.py \
  --junction Chimeric.out.junction \
  --genes   genes.bed \
  --outdir  results

# GTF instead of BED:            --genes gencode.v44.annotation.gtf --feature gene
# keep only supported junctions: --min-reads 3
# same-strand-only primary out:  --strand aware   (default: any = both strands captured,
#                                        every overlapping gene is listed, comma-joined)
# skip the large per-read file:  --no-reads-file   (writes only the two summary TSVs)
```

Try the bundled example first:

```bash
python3 scripts/annotate_junctions.py \
  --junction examples/chimeric_junction/Chimeric.out.junction \
  --genes    examples/chimeric_junction/genes.bed \
  --outdir   /tmp/junction_results

# folder of DRAGEN runs instead of a single file — every output row gets
# `source_file` + `sample` columns (sample = part of the file name before
# '.Chimeric.out.junction', e.g. SampleA.Chimeric.out.junction -> SampleA;
# bare Chimeric.out.junction falls back to the parent folder):
python3 scripts/annotate_junctions.py --junction path/to/runs_folder --genes genes.bed --outdir results

# just list the junction files a folder input would pick up (no annotation):
python3 scripts/annotate_junctions.py --junction path/to/runs_folder --list-files
scripts/list_junction_files.sh path/to/runs_folder    # shell-only equivalent
```

## 2b. Full pipeline runthrough

```bash
git pull
# inputs (all gitignored, nothing data-like can be committed):
#   data/SampleA.Chimeric.out.junction, data/SampleB.Chimeric.out.junction, ... , data/genes.xlsx

# 1. gene table: Excel "Save As -> Text (Tab delimited) .txt" (or xlsx_to_tsv.py) -> data/genes.tsv
# 2. optional cosmetic headers on the junction files (folder runs skip the copies):
for f in $(scripts/list_junction_files.sh data/); do scripts/add_junction_header.sh "$f"; done
# 3. annotate the whole folder at once — outputs are ALREADY pooled, one row per
#    split read with source_file + sample columns:
python3 scripts/annotate_junctions.py --junction data/ --genes data/genes.tsv --outdir data/results --min-reads 3
# 4. outputs in data/results/:
#      annotated_reads.tsv    pooled per-read annotations (all samples in one file)
#      junction_summary.tsv   per sample + breakpoint, sorted by n_reads
#      disagreements.tsv      where the 4 annotations / strand modes differ
# 5. only needed if you annotated runs into SEPARATE outdirs — pool them after:
python3 scripts/pool_annotated.py --which annotated_reads   --out data/results/pooled_reads.tsv   outA outB outC
python3 scripts/pool_annotated.py --which junction_summary  --out data/results/pooled_summary.tsv outA outB outC
python3 scripts/pool_annotated.py --which disagreements     --out data/results/pooled_disagreements.tsv outA outB outC
```

## 3. Outputs (in `--outdir`)

| File | Use it for |
|---|---|
| `annotated_reads.tsv` | one row per split read; leading `source_file` + `sample` columns (folder input); `gene_donor/acceptor/seg1/seg2` = primary annotation (both strands by default), `gene_*_aware` = same-strand-only view |
| `junction_summary.tsv` | one row per **sample +** donor/acceptor breakpoint, sorted by `n_reads`, with `event_class` |
| `disagreements.tsv` | the "see the differences" file: rows where the annotations differ (strand-aware vs any, donor vs acceptor gene, breakpoint vs segment gene), with a `reasons` column |

## 4. How to read the results

What gets annotated per junction line:

- **donor** = point (col1, col2) — first base of the donor intron (1-based → 0-based point)
- **acceptor** = point (col4, col5)
- **segment 1** = span on the donor chrom, col11 → col11 + reference length of the col12 CIGAR (stops at DRAGEN's `p` chimera-gap op; `M/D/N/=/X` consume reference, `S/H/I` don't)
- **segment 2** = span on the acceptor chrom from col13/col14

`event_class` in `junction_summary.tsv`:

- `interchromosomal` — donor/acceptor on different chromosomes (translocation / classic fusion candidate)
- `intragenic_same_gene` — same chromosome, same gene both sides (readthrough, tandem dup, circRNA)
- `same_chrom_diff_genes` — same chromosome, different genes
- `unannotated` — not inside any gene (check `gene_*_any` / intergenic)

Other columns:

- `strand_match = no` → donor and acceptor strands (cols 3/6) differ: usually artifact or non-coding; an in-frame fusion needs `strand_match = yes`.
- Both-strands capture means loci with overlapping sense/antisense genes show both names comma-joined (e.g. `GENE_A,GENE_ANTISENSE`). `disagreements.tsv` flags those as `strand_matters` so you can cross-check with your own tools.
- Donor/acceptor positions are the **first intronic base** (exon boundary); segment start/end is the actual aligned read span. Where `seg1_gene != donor_gene`, the breakpoint sits at an exon edge and the segment tells you which exon the read supports.

## 5. Optional: bedtools alternative

```bash
bedtools intersect -a <(awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2-1,$2,"donor",0,$3,$10}' Chimeric.out.junction) \
                  -b genes.bed -wa -wb        # add -s for same-strand only
# repeat for acceptor with $4,$5,$6; for segments derive the end from the CIGAR first
```
