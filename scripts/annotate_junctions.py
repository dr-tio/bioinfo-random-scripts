#!/usr/bin/env python3
"""
annotate_junctions.py — annotate DRAGEN RNA `Chimeric.out.junction` split-read
breakpoints with gene names.

--junction accepts a single file OR a folder; for a folder, all
*.junction / *.junction.gz files are collected recursively (same convention as
OGM-toolkit's process_smap.collect_smap_files). --list-files just prints them.

Each line of Chimeric.out.junction = one chimerically aligned read (14 columns):
  1  donor chromosome
  2  first base of the donor intron (1-based)
  3  donor strand
  4  acceptor chromosome
  5  first base of the acceptor intron (1-based)
  6  acceptor strand
  7-9 placeholders (1, *, *)
  10 read name
  11 first base of the first segment (on the + strand; donor chromosome)
  12 CIGAR of the first segment
  13 first base of the second segment (acceptor chromosome)
  14 CIGAR of the second segment

For every line we annotate ALL FOUR coordinate sets:
  donor      point  (col1, col2)                 strand = col3
  acceptor   point  (col4, col5)                 strand = col6
  segment 1  span   (col1, col11 .. +reflen(col12)) strand = col3
  segment 2  span   (col4, col13 .. +reflen(col14)) strand = col6

Each annotation is done twice: strand-agnostic (any strand — the DEFAULT primary
output, so overlapping genes on both strands are all captured) and strand-aware
(same strand only, kept in *_aware / *_any columns), so you can see where strand
would have changed the call. Switch primary with --strand aware.

Gene file accepted formats (auto-detected):
  - BED:  chrom start end [name] [score] [strand]   (0-based half-open)
  - GTF/GFF: use --feature gene (default) and gene_name/Name/gene_id is used
  - NCBI gene table (Accession/Begin/End/Chromosome/Orientation/Name/Symbol/
    Gene ID/...; tab- or pipe-delimited, header optional): Symbol is used as
    the gene name, Begin (1-based) is converted to 0-based, Orientation
    (plus/minus) becomes the strand, Chromosome (plain "1") is used as contig.
Chromosome naming is auto-reconciled (chr1 vs 1 vs NC_000001.11).

Outputs in --outdir:
  annotated_reads.tsv   one row per junction line (per split read); source_file +
                        sample columns identify the origin (sample = part of the
                        file name before '.Chimeric.out.junction')
  junction_summary.tsv  one row per sample + unique donor/acceptor breakpoint, read counts, event class
  disagreements.tsv     junctions where the annotations differ (strand-aware vs any,
                        donor vs acceptor gene, breakpoint vs segment gene)

Pure standard library — no dependencies. Python 3.8+.
"""

import argparse
import csv
import gzip
import os
import re
import sys
from bisect import bisect_right
from collections import defaultdict

# ---------------------------------------------------------------- CIGAR ----
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=Xp])")  # includes DRAGEN's lowercase 'p' gap op
REF_CONSUMING = set("MDN=X")  # ops that consume reference bases


def ref_len(cigar: str) -> int:
    """Reference span of one segment's CIGAR.

    M/D/N/=/X consume reference; S/H/I do not. The DRAGEN-specific 'p' op marks
    the chimera gap - everything after it belongs to the OTHER segment, so we
    stop there (otherwise segment 1 would wrongly extend past the breakpoint).
    """
    total = 0
    for m in CIGAR_RE.finditer(cigar):
        n, op = int(m.group(1)), m.group(2)
        if op == "p":
            break
        if op in REF_CONSUMING:
            total += n
    return total


# ------------------------------------------------------- chromosome map ----
CHR_TO_ACC = {
    "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12",
    "4": "NC_000004.12", "5": "NC_000005.10", "6": "NC_000006.12",
    "7": "NC_000007.14", "8": "NC_000008.11", "9": "NC_000009.12",
    "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
    "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10",
    "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10",
    "19": "NC_000019.10", "20": "NC_000020.11", "21": "NC_000021.9",
    "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10",
    "M": "NC_012920.1", "MT": "NC_012920.1",
}
ACC_TO_CHR = {v: k for k, v in CHR_TO_ACC.items()}


def norm_chrom(c: str, mode: str) -> str:
    if mode == "exact":
        return c
    if mode == "strip-chr":
        return c[3:] if c.lower().startswith("chr") else c
    if mode == "accession":
        if c.lower().startswith("chr"):
            return CHR_TO_ACC.get(c[3:], c)
        return ACC_TO_CHR.get(c, c)
    raise ValueError(mode)


def pick_chrom_mode(junction_chroms, gene_chroms):
    """Pick the naming mode that best aligns junction chroms to gene chroms."""
    best, best_cov = "exact", -1.0
    for mode in ("exact", "strip-chr", "accession"):
        g = {norm_chrom(c, mode) for c in gene_chroms}
        cov = sum(1 for c in junction_chroms if norm_chrom(c, mode) in g) / max(1, len(junction_chroms))
        if cov > best_cov:
            best, best_cov = mode, cov
    return best, best_cov


# ------------------------------------------------------------ gene index ---
class Gene:
    __slots__ = ("chrom", "start", "end", "name", "strand")

    def __init__(self, chrom, start, end, name, strand):
        self.chrom, self.start, self.end = chrom, start, end
        self.name, self.strand = name, strand


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r", encoding="utf-8", errors="replace")


def parse_attributes(attr: str):
    out = {}
    for part in attr.strip().rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        m = re.match(r'(\S+)\s+"?([^";]+)"?', part)
        if m:
            out[m.group(1)] = m.group(2).strip().strip('"')
    return out


NCBI_DEFAULT_COLS = {"accession": 0, "begin": 1, "end": 2, "chromosome": 3,
                     "orientation": 4, "name": 5, "symbol": 6, "gene id": 7, "gene type": 8}


def split_row(line):
    """Split a table row: pipe-delimited (markdown paste) or tab/whitespace.
    Tolerates CRLF endings and quote-wrapped cells (Excel text export)."""
    s = line.rstrip("\r\n")
    if s.lstrip().startswith("|"):
        cells = [c.strip() for c in s.strip().strip("|").split("|")]
    else:
        f = s.split("\t")
        cells = f if len(f) > 1 else s.split()
    return [c.strip('"') if c.startswith('"') and c.endswith('"') and len(c) > 1 else c
            for c in cells]


def ncbi_gene(cells, idx):
    def g(key):
        i = idx.get(key)
        return cells[i] if i is not None and i < len(cells) else ""
    chrom = g("chromosome") or g("accession")
    start = int(g("begin").replace(",", "")) - 1   # NCBI Begin is 1-based -> BED 0-based
    end = int(g("end").replace(",", ""))
    name = g("symbol") or g("name") or f"{chrom}:{start}-{end}"
    strand = {"plus": "+", "minus": "-", "+": "+", "-": "-"}.get(g("orientation").lower(), ".")
    return Gene(chrom, start, end, name, strand)


def load_genes(path, feature):
    """Return (list of Gene, format). Auto-detects BED vs GTF/GFF vs NCBI gene table."""
    genes, fmt = [], None
    ncbi_idx = dict(NCBI_DEFAULT_COLS)
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = split_row(line)
            if f and f[0].lower() in ("track", "browser"):
                continue  # UCSC BED header lines
            if fmt is None:
                low = [c.lower() for c in f]
                is_gff = len(f) >= 9 and f[3].isdigit() and f[4].isdigit() and not f[2].isdigit()
                if is_gff:
                    fmt = "gff"
                elif "accession" in low and "begin" in low:
                    fmt = "ncbi"                       # header row: map columns by name
                    ncbi_idx = {c.strip().lower(): i for i, c in enumerate(f)}
                    continue
                elif f[0].startswith("NC_"):
                    fmt = "ncbi"                       # headerless NCBI table, positional
                else:
                    fmt = "bed"
            if fmt == "gff":
                if len(f) < 9 or f[2] not in feature:
                    continue
                attrs = parse_attributes(f[8])
                name = (attrs.get("gene_name") or attrs.get("Name")
                        or attrs.get("gene_id") or attrs.get("ID")
                        or attrs.get("transcript_name") or f[2])
                genes.append(Gene(f[0], int(f[3]) - 1, int(f[4]), name, f[6] if f[6] in "+-" else "."))
            elif fmt == "ncbi":
                if len(f) < 3 or not f[ncbi_idx["begin"]].replace(",", "").isdigit():
                    continue
                genes.append(ncbi_gene(f, ncbi_idx))
            else:
                chrom = f[0]
                start = int(f[1].replace(",", ""))   # tolerate thousands separators
                end = int(f[2].replace(",", ""))
                name = f[3] if len(f) > 3 else f"{chrom}:{start}-{end}"
                strand_raw = f[5].lower() if len(f) > 5 else "."
                strand = {"+": "+", "-": "-", "plus": "+", "minus": "-"}.get(strand_raw, ".")
                genes.append(Gene(chrom, start, end, name, strand))
    return genes, fmt


class GeneIndex:
    """Per-chromosome sorted interval index with prefix-max-end pruning."""

    def __init__(self, genes, chrom_mode):
        self.mode = chrom_mode
        by = defaultdict(list)
        for g in genes:
            by[norm_chrom(g.chrom, chrom_mode)].append(g)
        self.by_chrom = {}
        for chrom, gs in by.items():
            gs.sort(key=lambda g: (g.start, g.end))
            starts = [g.start for g in gs]
            pmax, m = [], -1
            for g in gs:
                m = max(m, g.end)
                pmax.append(m)
            self.by_chrom[chrom] = (starts, pmax, gs)

    def query(self, chrom, qstart, qend, strand=None):
        """All genes overlapping BED interval [qstart, qend); strand filter if given."""
        hit = self.by_chrom.get(norm_chrom(chrom, self.mode))
        if not hit:
            return []
        starts, pmax, gs = hit
        res = []
        hi = bisect_right(starts, qend - 1)  # gene.start < qend  <=>  start <= qend-1
        i = hi - 1
        while i >= 0 and pmax[i] > qstart:   # else no gene in [0..i] can overlap
            g = gs[i]
            if g.end > qstart and (strand is None or g.strand in (strand, ".")):
                res.append(g)
            i -= 1
        return res


# ------------------------------------------------------------ junctions ----
JUNCTION_SUFFIXES = (".junction", ".junction.gz")


def collect_junction_files(input_path):
    """Same convention as OGM-toolkit process_smap.collect_smap_files:
    a file -> [file]; a directory -> every *.junction / *.junction.gz found
    recursively, sorted; anything else -> error."""
    if os.path.isfile(input_path):
        return [input_path]
    if os.path.isdir(input_path):
        found = []
        for root, _, files in os.walk(input_path):
            for f in sorted(files):
                # skip header-added copies (*.with_header.junction) so a folder
                # containing both original and headed files is not double-counted
                if ".with_header." in f:
                    continue
                if f.lower().endswith(JUNCTION_SUFFIXES):
                    found.append(os.path.join(root, f))
        return sorted(found)
    raise ValueError(f"input path not found: {input_path}")


def sample_name(path):
    """Sample label from the junction FILE NAME: the part before
    '.Chimeric.out.junction' (e.g. 'SampleA.Chimeric.out.junction' -> 'SampleA').
    Bare 'Chimeric.out.junction' -> parent folder name; other '*.junction'
    names -> file stem."""
    base = os.path.basename(path)
    if base.lower().endswith(".gz"):
        base = base[:-3]
    low = base.lower()
    marker = ".chimeric.out.junction"
    if marker in low:
        prefix = base[: low.index(marker)]
        if prefix:
            return prefix
    if low == "chimeric.out.junction":
        # bare DRAGEN file name -> parent folder
        parent = os.path.basename(os.path.dirname(os.path.abspath(path)))
        return parent or base
    if low.endswith(".junction"):
        return base[: -len(".junction")]
    return base


def parse_junctions(path):
    reads = []
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line[0] in "#@":
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 14:
                f = line.split()
            if len(f) < 14:
                sys.stderr.write(f"WARNING: skipping short line ({len(f)} cols): {line[:80]}\n")
                continue
            reads.append({
                "d_chr": f[0], "d_pos": int(f[1]), "d_strand": f[2],
                "a_chr": f[3], "a_pos": int(f[4]), "a_strand": f[5],
                "read": f[9],
                "s1_start": int(f[10]), "s1_cigar": f[11],
                "s2_start": int(f[12]), "s2_cigar": f[13],
            })
    return reads


def names(genes):
    return ",".join(sorted({g.name for g in genes})) if genes else "."


def annotate_read(r, idx):
    """Return the 4 strand-aware + 4 strand-agnostic gene annotations."""
    s1_end = r["s1_start"] - 1 + ref_len(r["s1_cigar"])
    s2_end = r["s2_start"] - 1 + ref_len(r["s2_cigar"])
    q = {
        "donor":    (r["d_chr"], r["d_pos"] - 1, r["d_pos"], r["d_strand"]),
        "acceptor": (r["a_chr"], r["a_pos"] - 1, r["a_pos"], r["a_strand"]),
        "seg1":     (r["d_chr"], r["s1_start"] - 1, s1_end, r["d_strand"]),
        "seg2":     (r["a_chr"], r["s2_start"] - 1, s2_end, r["a_strand"]),
    }
    out = {}
    for k, (c, s, e, st) in q.items():
        out[k] = names(idx.query(c, s, e, st))
        out[k + "_any"] = names(idx.query(c, s, e, None))
    out["s1_end"], out["s2_end"] = s1_end, s2_end
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--junction", required=True,
                    help="a Chimeric.out.junction file (optionally .gz) OR a folder — for a "
                         "folder, all *.junction / *.junction.gz files are collected "
                         "recursively (OGM-toolkit collect_smap_files convention)")
    ap.add_argument("--list-files", action="store_true",
                    help="only list the junction file(s) found for --junction and exit")
    ap.add_argument("--genes", help="gene annotation: BED or GTF/GFF (optionally .gz); required unless --list-files")
    ap.add_argument("--outdir", help="output directory; required unless --list-files")
    ap.add_argument("--feature", default="gene",
                    help="GFF/GTF feature(s) to use, comma-separated (default: gene)")
    ap.add_argument("--chrom-mode", choices=["auto", "exact", "strip-chr", "accession"], default="auto",
                    help="chromosome naming reconciliation (default: auto)")
    ap.add_argument("--min-reads", type=int, default=1,
                    help="only keep junctions with >= this many supporting reads in junction_summary (default: 1)")
    ap.add_argument("--strand", choices=["any", "aware"], default="any",
                    help="primary annotation mode: 'any' = capture genes on BOTH strands (default); "
                         "'aware' = same-strand only. The other mode is always kept in *_aware/*_any columns")
    args = ap.parse_args()
    if not args.list_files and (not args.genes or not args.outdir):
        ap.error("--genes and --outdir are required unless --list-files")

    if args.list_files:
        try:
            jfiles = collect_junction_files(args.junction)
        except ValueError as e:
            sys.exit(f"ERROR: {e}")
        for p in jfiles:
            print(p)
        print(f"{len(jfiles)} junction file(s) found", file=sys.stderr)
        return

    if args.outdir:
        os.makedirs(args.outdir, exist_ok=True)
    feature = set(args.feature.split(","))

    print(f"[1/4] loading genes from {args.genes} ...", flush=True)
    genes, fmt = load_genes(args.genes, feature)
    print(f"      {len(genes):,} features loaded (format: {fmt})", flush=True)
    if not genes:
        sys.exit("ERROR: no gene records parsed — check --genes file format")

    print(f"[2/4] collecting junction files from {args.junction} ...", flush=True)
    try:
        jfiles = collect_junction_files(args.junction)
    except ValueError as e:
        sys.exit(f"ERROR: {e}")
    if not jfiles:
        sys.exit(f"ERROR: no *.junction / *.junction.gz files found under {args.junction}")
    reads = []
    for p in jfiles:
        rs = parse_junctions(p)
        for r in rs:
            r["source_file"] = p
            r["sample"] = sample_name(p)
        reads.extend(rs)
        print(f"        {p}: {len(rs):,} lines", flush=True)
    print(f"      {len(reads):,} junction lines loaded from {len(jfiles)} file(s)", flush=True)
    if not reads:
        sys.exit("ERROR: no junction lines parsed — check --junction input")

    jch = {r["d_chr"] for r in reads} | {r["a_chr"] for r in reads}
    gch = {g.chrom for g in genes}
    mode = args.chrom_mode if args.chrom_mode != "auto" else pick_chrom_mode(jch, gch)[0]
    cov = pick_chrom_mode(jch, gch)[1]
    print(f"[3/4] chromosome mode: {mode} (junction chroms found in gene file: {cov:.1%})", flush=True)
    idx = GeneIndex(genes, mode)

    print("[4/4] annotating (donor / acceptor / segment1 / segment2, strand-aware + any) ...", flush=True)
    read_path = os.path.join(args.outdir, "annotated_reads.tsv")
    hits = defaultdict(int)
    groups = defaultdict(list)
    sec_suffix = "_aware" if args.strand == "any" else "_any"

    def prim(a, k):
        """Primary annotation for coordinate k, per --strand mode."""
        return a[k + "_any"] if args.strand == "any" else a[k]

    def sec(a, k):
        """The other mode, always kept for comparison."""
        return a[k] if args.strand == "any" else a[k + "_any"]

    with open(read_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "source_file",
            "sample",
            "read_name",
            "donor_chr", "donor_pos", "donor_strand",
            "acceptor_chr", "acceptor_pos", "acceptor_strand",
            "seg1_chr", "seg1_start", "seg1_end",
            "seg2_chr", "seg2_start", "seg2_end",
            "gene_donor", "gene_acceptor", "gene_seg1", "gene_seg2",
            "gene_donor" + sec_suffix, "gene_acceptor" + sec_suffix,
            "gene_seg1" + sec_suffix, "gene_seg2" + sec_suffix,
            "same_chrom", "strand_match",
        ])
        for r in reads:
            a = annotate_read(r, idx)
            for k in ("donor", "acceptor", "seg1", "seg2"):
                if prim(a, k) != ".":
                    hits[k] += 1
            key = (r["sample"], r["d_chr"], r["d_pos"], r["d_strand"], r["a_chr"], r["a_pos"], r["a_strand"])
            groups[key].append((r, a))
            w.writerow([
                r["source_file"], r["sample"], r["read"], r["d_chr"], r["d_pos"], r["d_strand"],
                r["a_chr"], r["a_pos"], r["a_strand"],
                r["d_chr"], r["s1_start"], a["s1_end"],
                r["a_chr"], r["s2_start"], a["s2_end"],
                prim(a, "donor"), prim(a, "acceptor"), prim(a, "seg1"), prim(a, "seg2"),
                sec(a, "donor"), sec(a, "acceptor"), sec(a, "seg1"), sec(a, "seg2"),
                "yes" if r["d_chr"] == r["a_chr"] else "no",
                "yes" if r["d_strand"] == r["a_strand"] else "no",
            ])

    # ---- junction-level summary ----
    def event_class(r0, gd, ga):
        if r0["d_chr"] != r0["a_chr"]:
            return "interchromosomal"
        if gd == "." or ga == ".":
            return "unannotated"
        if gd == ga:
            return "intragenic_same_gene"
        return "same_chrom_diff_genes"

    summ_path = os.path.join(args.outdir, "junction_summary.tsv")
    dis_path = os.path.join(args.outdir, "disagreements.tsv")
    n_dis = 0
    with open(summ_path, "w", newline="") as fs, open(dis_path, "w", newline="") as fd:
        ws = csv.writer(fs, delimiter="\t")
        wd = csv.writer(fd, delimiter="\t")
        ws.writerow(["n_reads", "sample", "donor_chr", "donor_pos", "donor_strand",
                     "acceptor_chr", "acceptor_pos", "acceptor_strand",
                     "gene_donor", "gene_acceptor", "gene_seg1", "gene_seg2",
                     "strand_match", "event_class"])
        wd.writerow(["n_reads", "sample", "donor_chr", "donor_pos", "donor_strand",
                     "acceptor_chr", "acceptor_pos", "acceptor_strand",
                     "gene_donor", "gene_acceptor", "gene_seg1", "gene_seg2",
                     "gene_donor_any", "gene_acceptor_any", "gene_seg1_any", "gene_seg2_any",
                     "reasons"])
        rows = []
        for key, items in groups.items():
            if len(items) < args.min_reads:
                continue
            r0, a0 = items[0]
            gd = ",".join(sorted({prim(a, "donor") for _, a in items} - {"."})) or "."
            ga = ",".join(sorted({prim(a, "acceptor") for _, a in items} - {"."})) or "."
            g1 = ",".join(sorted({prim(a, "seg1") for _, a in items} - {"."})) or "."
            g2 = ",".join(sorted({prim(a, "seg2") for _, a in items} - {"."})) or "."
            strand_match = "yes" if r0["d_strand"] == r0["a_strand"] else "no"
            rows.append((len(items), key, gd, ga, g1, g2, strand_match, event_class(r0, gd, ga), a0))
        rows.sort(key=lambda x: -x[0])
        for n, key, gd, ga, g1, g2, sm, ec, a0 in rows:
            ws.writerow([n, *key, gd, ga, g1, g2, sm, ec])
            reasons = []
            if a0["donor"] != a0["donor_any"] or a0["acceptor"] != a0["acceptor_any"] \
               or a0["seg1"] != a0["seg1_any"] or a0["seg2"] != a0["seg2_any"]:
                reasons.append("strand_matters")
            if gd != ga:
                reasons.append("donor_vs_acceptor_gene_differ")
            if gd != "." and g1 != "." and gd not in g1:
                reasons.append("seg1_gene != donor_gene")
            if ga != "." and g2 != "." and ga not in g2:
                reasons.append("seg2_gene != acceptor_gene")
            if reasons:
                n_dis += 1
                wd.writerow([n, *key, gd, ga, g1, g2,
                             a0["donor_any"], a0["acceptor_any"], a0["seg1_any"], a0["seg2_any"],
                             ";".join(reasons)])

    print("\n===== done =====")
    print(f"junction lines annotated : {len(reads):,}")
    print(f"primary strand mode: {args.strand}")
    for k in ("donor", "acceptor", "seg1", "seg2"):
        print(f"  {k:>9}: {hits[k]:,} lines with >=1 gene hit ({hits[k]/len(reads):.1%})")
    print(f"unique breakpoints       : {len(groups):,}  (>= {args.min_reads} reads: {len(rows):,})")
    print(f"disagreement rows        : {n_dis:,}")
    unmatched = sorted({norm_chrom(c, mode) for c in jch} - set(idx.by_chrom))
    if unmatched:
        print(f"WARNING: junction chroms absent from gene file: {', '.join(unmatched)}")
    print(f"\noutputs in {args.outdir}/:")
    print("  annotated_reads.tsv   — one row per split read, all 4 annotations (strand-aware + any)")
    print("  junction_summary.tsv  — one row per breakpoint, sorted by supporting reads")
    print("  disagreements.tsv     — breakpoints where the 4 annotations / strand modes differ")


if __name__ == "__main__":
    main()
