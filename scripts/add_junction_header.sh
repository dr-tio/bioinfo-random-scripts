#!/usr/bin/env bash
# add_junction_header.sh — prepend a descriptive header to a headerless DRAGEN
# RNA `Chimeric.out.junction` file (the file itself has no header line).
#
# Column meanings (see Illumina DRAGEN docs, "Chimeric.out.junction File"):
#   1 donor_chr        chromosome of the donor
#   2 donor_pos        first base of the donor intron (1-based)
#   3 donor_strand     strand of the donor
#   4 acceptor_chr     chromosome of the acceptor
#   5 acceptor_pos     first base of the acceptor intron (1-based)
#   6 acceptor_strand  strand of the acceptor
#   7-9 ph1..ph3       placeholders (always 1, *, *)
#   10 read_name       read name
#   11 seg1_start      first base of the first segment (on the + strand)
#   12 seg1_cigar      CIGAR of the first segment (DRAGEN 'p' = chimera gap)
#   13 seg2_start      first base of the second segment
#   14 seg2_cigar      CIGAR of the second segment
#
# Usage:
#   add_junction_header.sh <in.junction> [out.junction]
# If no output path is given, writes <in>.with_header.junction next to the input.
# Lines starting with '#' are ignored by annotate_junctions.py, so both the
# original and the headed file work as input.
set -euo pipefail

in="${1:?usage: add_junction_header.sh <in.junction> [out.junction]}"
out="${2:-${in%.junction}.with_header.junction}"

{ printf '#donor_chr\tdonor_pos\tdonor_strand\tacceptor_chr\tacceptor_pos\tacceptor_strand\tph1\tph2\tph3\tread_name\tseg1_start\tseg1_cigar\tseg2_start\tseg2_cigar\n'
  cat "$in"
} > "$out"

echo "wrote $out ($(wc -l < "$out") lines incl. header)"
