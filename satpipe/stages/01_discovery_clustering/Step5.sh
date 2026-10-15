#!/usr/bin/env bash
# Step 5 - canonical circular phase: rotate every motif to its lexicographically-minimal
# rotation so that identical monomers in different phases become identical strings before
# clustering. PLUS STRAND ONLY: the reverse complement is not considered, which matches
# vsearch's default --strand plus in Steps 6 and 8. (This is why the single L541 locus,
# whose TRF motif was reported in the opposite orientation, never clustered; see README.md.)
# Input motifs are uppercase (Step3.sh uppercases), so byte-order comparison is safe.
#
# Usage: Step5.sh [input_dir]   (default: current directory; reads *.clean.fas)
# Output: <input_dir>/clean/<file>.rot   (Steps 6-7 are then run from inside clean/)
set -euo pipefail
IN="${1:-.}"
mkdir -p "$IN/clean"
shopt -s nullglob
n=0
for fa in "$IN"/*.clean.fas; do
  awk '/^>/ {print; next}
       { s=$0; n=length(s); min=s; d=s s;
         for (i=2; i<=n; i++) { r=substr(d,i,n); if (r<min) min=r }
         print min }' "$fa" > "$IN/clean/$(basename "$fa").rot"
  n=$((n+1))
done
echo "Step 5: rotated $n file(s) into $IN/clean/"
