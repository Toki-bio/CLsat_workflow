#!/usr/bin/env bash
# Step 5 - canonical circular phase: rotate every motif to its lexicographically-minimal
# rotation so that identical monomers in different phases become identical strings before
# clustering.
#
# STRAND=plus (default, the published behaviour): the reverse complement is NOT considered,
# which matches vsearch's --strand plus in Steps 6 and 8. This is why the single L541 locus,
# whose TRF motif was reported in the opposite orientation, never clustered (see PIPELINE.md).
# STRAND=both: the minimum is taken over the rotations of the motif AND of its reverse
# complement, so opposite-strand reports of the same monomer become identical. Changes the
# published numbers; opt-in.
#
# Input motifs are uppercase (Step3.sh uppercases), so byte-order comparison is safe.
#
# Usage: Step5.sh [input_dir]   (default: current directory; reads *.clean.fas)
# Output: <input_dir>/clean/<file>.rot   (Steps 6-7 are then run from inside clean/)
set -euo pipefail
IN="${1:-.}"
STRAND="${STRAND:-plus}"
[[ "$STRAND" == plus || "$STRAND" == both ]] || { echo "STRAND must be plus or both, got '$STRAND'" >&2; exit 1; }
mkdir -p "$IN/clean"
shopt -s nullglob
n=0
for fa in "$IN"/*.clean.fas; do
  awk -v strand="$STRAND" '
       function minrot(s,    n, min, d, i, r) {
         n=length(s); min=s; d=s s
         for (i=2; i<=n; i++) { r=substr(d,i,n); if (r<min) min=r }
         return min
       }
       function revcomp(s,    n, i, c, out) {
         n=length(s); out=""
         for (i=n; i>=1; i--) {
           c=substr(s,i,1)
           if (c=="A") c="T"; else if (c=="T") c="A"; else if (c=="C") c="G"; else if (c=="G") c="C"
           out=out c
         }
         return out
       }
       /^>/ {print; next}
       { m=minrot($0)
         if (strand=="both") { r=minrot(revcomp($0)); if (r<m) m=r }
         print m }' "$fa" > "$IN/clean/$(basename "$fa").rot"
  n=$((n+1))
done
echo "Step 5: rotated $n file(s) into $IN/clean/ (strand=$STRAND)"
