#!/usr/bin/env bash
# =============================================================================
# Retrospective screen of every singleton locus against the family consensuses
# (manuscript Methods 3.1).
#
# WHY THIS EXISTS
#   Step 6 clusters each assembly's rotated monomers with
#       vsearch --cluster_fast --id 0.90
#   which defaults to --strand plus, so a monomer and its reverse complement never
#   join the same cluster. Only the aligned members of MULTI-locus clusters are
#   pooled into the cross-assembly pass (step 8), so a cluster of one -- a
#   "singleton locus" -- can never become a seed, however strong its signal.
#
#   Note that "singleton" describes cluster membership, not array size: among the
#   979 Darevskia singleton loci monomer copy number runs from 10 to 418 (median 13.7),
#   and 39 of them carry 50 or more copies.
#
# WHAT IT PRODUCES
#   sing_best.tsv  - best family and bitscore per singleton locus
#   Hits above the MIN_BITS floor are then intersected against the catalogue
#   (if given) to separate: already inside a catalogued locus / absent from it.
#
# RESULT AS RUN (2026-07-27, Darevskia)
#   979 singleton loci; 14 above 50 bits, all above 77; 7 inside the 541-locus catalogue
#   (6 swallowed by neighbouring loci + L541 itself, added by hand); 7 outside -> C1-C7.
#
# Usage: singleton_screen.sh <clusters_by_file_dir> <rotated_references.fa> <out_dir> [catalogue.bed]
#   clusters_by_file_dir   Step-6 output: one subdirectory per species (<code>.fa.clean.fas/)
#   rotated_references.fa  the family consensuses (CLsat_cons.rot in the Darevskia run)
#   catalogue.bed          merged catalogue to intersect the hits with (optional)
# Environment: MIN_BITS [50]
# Tools: ssearch36, bedtools (if catalogue given) on PATH.
# =============================================================================
set -euo pipefail
[[ $# -ge 3 ]] || { echo "Usage: $0 <clusters_by_file_dir> <rotated_references.fa> <out_dir> [catalogue.bed]" >&2; exit 1; }
CB="$1"; REF="$2"; OUT="$3"; BED="${4:-}"
MIN_BITS="${MIN_BITS:-50}"
command -v ssearch36 >/dev/null || { echo "ssearch36 not in PATH" >&2; exit 1; }
[[ -s "$REF" ]] || { echo "reference FASTA missing: $REF" >&2; exit 1; }
mkdir -p "$OUT"; cd "$OUT"

# --- 1. collect every cluster holding exactly one locus ----------------------
: > singletons.fa
for d in "$CB"/*/; do
  [ -d "$d" ] || continue
  a=$(basename "$d"); a="${a%%.*}"
  n=0
  for f in "$d"/*cluster*.fa; do
    [ -f "$f" ] || continue
    if [ "$(grep -c '^>' "$f")" -le 1 ]; then cat "$f" >> singletons.fa; n=$((n+1)); fi
  done
  printf '%-4s %4d singleton loci\n' "$a" "$n"
done
echo "total singleton loci: $(grep -c '^>' singletons.fa)"

# --- 2. double the consensuses so the search is phase-insensitive ------------
awk '/^>/{if(s)print h"\n"s s; h=$0; s=""; next}{s=s$0}END{print h"\n"s s}' \
    "$REF" > ref_doubled.fa

# --- 3. search every singleton against all families --------------------------
ssearch36 -n -m 8C -E 100 -b 5000 singletons.fa ref_doubled.fa 2>/dev/null > sing_vs_ref.tsv

# best family per singleton locus
awk '!/^#/{if($12>b[$1]){b[$1]=$12; f[$1]=$2}}
     END{for(q in b) printf "%s\t%s\t%.1f\n", q, f[q], b[q]}' sing_vs_ref.tsv \
  | sort -k3,3gr > sing_best.tsv

echo
echo "singleton loci scoring >= $MIN_BITS bits against a family consensus:"
awk -v M="$MIN_BITS" '$3>=M' sing_best.tsv

# --- 4. which of those are already inside a catalogued locus? ----------------
if [ -n "$BED" ] && [ -f "$BED" ]; then
  command -v bedtools >/dev/null || { echo "bedtools not in PATH" >&2; exit 1; }
  awk -v M="$MIN_BITS" '$3>=M{print $1}' sing_best.tsv \
    | awk -F'|' '{split($2,a,":"); split(a[2],b,"-"); print a[1]"\t"b[1]-1"\t"b[2]"\t"$0}' \
    | sort -k1,1 -k2,2n > hits.bed
  echo
  echo "hits already covered by a catalogued locus:"
  bedtools intersect -a hits.bed -b <(sort -k1,1 -k2,2n "$BED") -u | cut -f4
  echo
  echo "hits ABSENT from the catalogue (candidates for inspection):"
  bedtools intersect -a hits.bed -b <(sort -k1,1 -k2,2n "$BED") -v | cut -f4 | tee candidates_absent.txt
fi
