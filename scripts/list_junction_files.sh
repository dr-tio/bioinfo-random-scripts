#!/usr/bin/env bash
# list_junction_files.sh — list DRAGEN Chimeric.out.junction files, given a
# single file OR a folder (searched recursively). Shell-only companion to
# annotate_junctions.py --list-files; same collect convention as
# OGM-toolkit's process_smap.collect_smap_files.
#
# Usage:
#   list_junction_files.sh <file-or-folder>
set -euo pipefail

input="${1:?usage: list_junction_files.sh <file-or-folder>}"

if [[ -f "$input" ]]; then
    echo "$input"
    exit 0
fi
if [[ ! -d "$input" ]]; then
    echo "error: not found: $input" >&2
    exit 1
fi

find "$input" -type f \( -name '*.junction' -o -name '*.junction.gz' \) | sort
