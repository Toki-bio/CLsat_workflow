#!/usr/bin/env bash
# =============================================================================
# Retrospective screen of every singleton locus against the four CLsat consensuses
# (manuscript Methods 3.1).
#
# WHY THIS EXISTS
#   Step 6 clusters each assembly's rotated monomers with
#       vsearch --cluster_fast --id 0.90
#   which defaults to --strand plus, so a monomer and its reverse complement never
#   join the same cluster. Only the aligned members of MULTI-locus clusters are
#   pooled into the cross-assembly pass (step 8), so a cluster of one -- a
#   "singleton locus" -- can never become a seed, however strong its CLsat signal.
#
#   Note that "singleton" describes cluster membership, not array size: among the
#   979 singleton loci monomer copy number runs from 10 to 418 (median 13.7), and
#   39 of them carry 50 or more copies.
#
# WHAT IT PRODUCES
#   sing_best.tsv  - best CLsat family and bitscore per singleton locus
#   Hits above the 50-bit noise floor are then intersected against the catalogue
#   (loci_540.bed) to separate: already inside a catalogued locus / absent from it.
#
# RESULT AS RUN (2026-07-27)
#   979 singleton loci; 14 above 50 bits, all above 77; 7 inside catalogued loci
#   (including L541, added by hand); 7 outside -> the C1-C7 candidates.
# =============================================================================
set -euo pipefail

CLEAN=/staging/Sergei/Lacertidae/Darevskia/Satellites/motifs_output/clean
CB=$CLEAN/clusters_by_file
REF=/staging/Sergei/Lacertidae/Darevskia/Satellites/CLUSTERS/CLsat/consi_CLsat_original/CLsat_cons.rot
SS=/staging/conda/envs/bioinfo/bin/ssearch36
OUT=$CLEAN/singleton_screen
mkdir -p "$OUT"; cd "$OUT"

# --- 1. collect every cluster holding exactly one locus ----------------------
: > singletons.fa
for a in arm mix nai unm unp val; do
  n=0
  for f in "$CB/$a.fa.clean.fas"/*cluster*.fa; do
    [ -f "$f" ] || continue
    if [ "$(grep -c '^>' "$f")" -le 1 ]; then cat "$f" >> singletons.fa; n=$((n+1)); fi
  done
  printf '%-4s %4d singleton loci\n' "$a" "$n"
done
echo "total singleton loci: $(grep -c '^>' singletons.fa)"

# --- 2. double the CLsat consensuses so the search is phase-insensitive ------
awk '/^>/{if(s)print h"\n"s s; h=$0; s=""; next}{s=s$0}END{print h"\n"s s}' \
    "$REF" > clsat_doubled.fa

# --- 3. search every singleton against all four families --------------------
"$SS" -n -m 8C -E 100 -b 5000 singletons.fa clsat_doubled.fa 2>/dev/null > sing_vs_clsat.tsv

# best family per singleton locus
awk '!/^#/{if($12>b[$1]){b[$1]=$12; f[$1]=$2}}
     END{for(q in b) printf "%s\t%s\t%.1f\n", q, f[q], b[q]}' sing_vs_clsat.tsv \
  | sort -k3,3gr > sing_best.tsv

echo
echo "singleton loci scoring >= 50 bits against a CLsat consensus:"
awk '$3>=50' sing_best.tsv

# --- 4. which of those are already inside a catalogued locus? ----------------
# loci_540.bed is the 541-locus merged catalogue from step 06.
BED=/staging/Sergei/Lacertidae/Darevskia/Satellites/manuscript/loci_540.bed
if [ -f "$BED" ]; then
  awk '$3>=50{print $1}' sing_best.tsv \
    | awk -F'|' '{split($2,a,":"); split(a[2],b,"-"); print a[1]"\t"b[1]-1"\t"b[2]"\t"$0}' \
    | sort -k1,1 -k2,2n > hits.bed
  echo
  echo "hits already covered by a catalogued locus:"
  bedtools intersect -a hits.bed -b <(sort -k1,1 -k2,2n "$BED") -u | cut -f4
  echo
  echo "hits ABSENT from the catalogue (the C1-C7 candidates):"
  bedtools intersect -a hits.bed -b <(sort -k1,1 -k2,2n "$BED") -v | cut -f4
fi
