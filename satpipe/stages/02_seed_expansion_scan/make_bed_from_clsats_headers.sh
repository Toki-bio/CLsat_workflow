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

# The position is field 2 of the header (sp|SCAFFOLD:START-END|...), split at its LAST ':' so any contig
# name works. A header without it is an error: it used to be dropped silently, and the old pattern
# ([A-Za-z0-9_.-]+:start-end anywhere) misread contig names holding other characters (GLM audit 2026-10-24).
awk 'BEGIN{OFS="\t"}
/^>/{
  hdr=substr($0,2)
  n=split(hdr, f, "|")
  pos=(n>=2 ? f[2] : "")
  c=0; for(i=length(pos); i>0; i--) if(substr(pos,i,1)==":"){ c=i; break }
  scaf=(c>1 ? substr(pos,1,c-1) : ""); rng=(c ? substr(pos,c+1) : "")
  if (scaf=="" || rng !~ /^[0-9]+-[0-9]+$/) { print "FATAL: no SCAFFOLD:START-END in field 2 of header: " hdr > "/dev/stderr"; bad=1; exit }
  split(rng, r, "-"); s=r[1]+0; e=r[2]+0
  if (s<1 || e<s) { print "FATAL: bad coordinates in header: " hdr > "/dev/stderr"; bad=1; exit }
  # BED is 0-based start, half-open end
  print scaf, s-1, e, hdr, 0, "."
}
END{ exit bad }' "$IN"
