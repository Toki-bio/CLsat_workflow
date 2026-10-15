#!/usr/bin/env bash
set -euo pipefail

# Build BED from FASTA headers that contain "...|SCAFFOLD:START-END|..."
# Keeps full header as BED name.
#
# Usage:
#   make_bed_from_clsats_headers.sh in.fa > out.bed
#
# Output: BED6 (chrom, start0, end0, name, 0, .)

IN="${1:-}"
[[ -n "$IN" && -s "$IN" ]] || { echo "ERROR: input fasta missing/empty" >&2; exit 2; }

awk 'BEGIN{FS=OFS="\t"}
/^>/{
  hdr=substr($0,2)

  # Find first occurrence of  SCAFFOLD:START-END  anywhere in header
  # Example: arm|JAWWNF010000123.1:2760899-2800411|len=...
  if (match(hdr, /([A-Za-z0-9_.-]+):([0-9]+)-([0-9]+)/, m)) {
    scaf=m[1]; s=m[2]+0; e=m[3]+0;
    if (s<1 || e<1) next;
    # BED is 0-based start, half-open end
    print scaf, s-1, e, hdr, 0, "."
  } else {
    # no coords
    # print "WARN: no coords in header: " hdr > "/dev/stderr"
  }
}' "$IN"
