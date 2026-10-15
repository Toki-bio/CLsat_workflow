#!/usr/bin/env bash
set -euo pipefail
die(){ echo "ERROR: $*" >&2; exit 1; }
export LC_ALL=C   # byte-order sort everywhere, so locus numbering (L1..) does not depend on locale

# Step 12 - MASTER locus universe: union of every family's expanded intervals, merged per
# species with family provenance.
#
# Usage: build_master.sh [families_expansion.tsv] [out_dir]
#   families_expansion.tsv   family<TAB>expanded_merged.bed  (one row per family; '#' lines ignored),
#                            the PASS1+2 output of clisat_expand_then_scan_merged.sh
#                            <seeds> <genomes.tsv> <cons> <MONO> out_<family>.
#                            Without it, every out_*.expanded_merged.bed in the current
#                            directory is used and the family name is read from the file name.
#   out_dir                  [default: ./out]
# Output: <out_dir>/MASTER.merged.bed  (chrom start end code|SCF:START1-END|MASTER provenance)
# History: the published catalogue was built from v6full_CLsat1 (cli_v6, fixed --min-aln-frac)
# and the Jun-12 out_CLsat{2,3,4}_v4 files; see PIPELINE.md.
TAB="${1:-}"
OUTD="${2:-${PWD}/out}"
mkdir -p "$OUTD"

declare -a FAMS=() BEDS=()
if [[ -n "$TAB" ]]; then
  [[ -s "$TAB" ]] || die "families table missing/empty: $TAB"
  while IFS=$'\t' read -r fam bed _; do
    [[ -z "$fam" || "$fam" == \#* ]] && continue
    FAMS+=("$fam"); BEDS+=("$bed")
  done < "$TAB"
else
  shopt -s nullglob
  for f in out_*.expanded_merged.bed; do
    fam="${f#out_}"; fam="${fam%.expanded_merged.bed}"
    FAMS+=("$fam"); BEDS+=("$f")
  done
fi
(( ${#FAMS[@]} > 0 )) || die "no families: give families_expansion.tsv or run where out_<family>.expanded_merged.bed files are"

echo "=== step 1: concat ${#FAMS[@]} family BEDs + add family col ==="
for i in "${!FAMS[@]}"; do
  f="${BEDS[$i]}"
  [[ -s "$f" ]] || die "missing: $f (family ${FAMS[$i]})"
  awk -v FAM="${FAMS[$i]}" 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,FAM}' "$f"
done | sort -k1,1 -k2,2n > "$OUTD/all_fams.bed"
echo "  $(wc -l < "$OUTD/all_fams.bed") rows pre-merge"

echo "=== step 2: split by species ==="
rm -f "$OUTD"/tmp.*.bed
awk -v D="$OUTD" 'BEGIN{FS=OFS="\t"}{n=split($4,a,"|"); sp=a[1]
  print $1,$2,$3,$4,$5 > (D "/tmp." sp ".bed")}' "$OUTD/all_fams.bed"
for f in "$OUTD"/tmp.*.bed; do
  sp=$(basename "$f" .bed | sed 's/^tmp\.//')
  echo "  $sp: $(wc -l < "$f") rows"
done

echo "=== step 3: bedtools merge per species with provenance ==="
rm -f "$OUTD"/MASTER.*.merged.bed
for f in "$OUTD"/tmp.*.bed; do
  sp=$(basename "$f" .bed | sed 's/^tmp\.//')
  sort -k1,1 -k2,2n "$f" \
  | bedtools merge -i - -c 4,5 -o distinct,distinct -delim ',' \
  | awk -v sp="$sp" 'BEGIN{FS=OFS="\t"}{
      prov=($5?$5:"."); tag=sp "|" $1 ":" ($2+1) "-" $3 "|MASTER"
      print $1,$2,$3,tag,prov
    }' > "$OUTD/MASTER.${sp}.merged.bed"
  echo "  $sp: $(wc -l < "$OUTD/MASTER.${sp}.merged.bed") merged"
done

echo "=== step 4: recombine ==="
cat "$OUTD"/MASTER.*.merged.bed \
  | sort -k1,1 -k2,2n > "$OUTD/MASTER.merged.bed"

# cleanup temp
rm -f "$OUTD"/all_fams.bed "$OUTD"/tmp.*.bed

echo ""
echo "DONE — MASTER: $(wc -l < "$OUTD/MASTER.merged.bed") loci"
echo "Provenance:"
awk -F'\t' '{print $5}' "$OUTD/MASTER.merged.bed" | sort | uniq -c | sort -rn
