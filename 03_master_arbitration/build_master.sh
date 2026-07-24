#!/usr/bin/env bash
set -euo pipefail
die(){ echo "ERROR: $*" >&2; exit 1; }
BASE="${PWD}"; mkdir -p "$BASE/out"

echo "=== step 1: concat 4 family BEDs + add family col ==="
for fam in 1 2 3 4; do
  f="out_CLsat${fam}_v4.expanded_merged.bed"
  [[ -s "$f" ]] || die "missing: $f"
  awk -v FAM="CLsat${fam}" 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,FAM}' "$f"
done | sort -k1,1 -k2,2n > "$BASE/out/all_fams.bed"
echo "  $(wc -l < "$BASE/out/all_fams.bed") rows pre-merge"

echo "=== step 2: split by species ==="
awk 'BEGIN{FS=OFS="\t"}{n=split($4,a,"|"); sp=a[1]
  print $1,$2,$3,$4,$5 > ("out/tmp." sp ".bed")}' "$BASE/out/all_fams.bed"
for f in "$BASE"/out/tmp.*.bed; do
  sp=$(basename "$f" .bed | sed 's/^tmp\.//')
  echo "  $sp: $(wc -l < "$f") rows"
done

echo "=== step 3: bedtools merge per species with provenance ==="
for f in "$BASE"/out/tmp.*.bed; do
  sp=$(basename "$f" .bed | sed 's/^tmp\.//')
  sort -k1,1 -k2,2n "$f" \
  | bedtools merge -i - -c 4,5 -o distinct,distinct -delim ',' \
  | awk -v sp="$sp" 'BEGIN{FS=OFS="\t"}{
      prov=($5?$5:"."); tag=sp "|" $1 ":" ($2+1) "-" $3 "|MASTER"
      print $1,$2,$3,tag,prov
    }' > "$BASE/out/MASTER.${sp}.merged.bed"
  echo "  $sp: $(wc -l < "$BASE/out/MASTER.${sp}.merged.bed") merged"
done

echo "=== step 4: recombine ==="
cat "$BASE"/out/MASTER.*.merged.bed \
  | sort -k1,1 -k2,2n > "$BASE/out/MASTER.merged.bed"

# cleanup temp
rm -f "$BASE"/out/all_fams.bed "$BASE"/out/tmp.*.bed

echo ""
echo "DONE — MASTER: $(wc -l < "$BASE/out/MASTER.merged.bed") loci"
echo "Provenance:"
awk -F'\t' '{print $5}' "$BASE/out/MASTER.merged.bed" | sort | uniq -c | sort -rn
