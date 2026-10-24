#!/bin/bash
set -euo pipefail
# Sub-subfamily assignment: the 2-monomer window at the array start (LB region; RB as fallback)
# is aligned against the doubled consensus of every sub-subfamily of the locus's family; the
# highest bitscore wins and the margin to the second best is recorded.
#
# Usage: assign_subfamily.sh <work_dir> <subfamily_consensus_dir> <members.tsv>
#   work_dir                 holds locus_master.tsv, manifest.tsv, all_regions.fa (from run_extract.sh);
#                            receives locus_subfamily.tsv
#   subfamily_consensus_dir  one FASTA per sub-subfamily, <name>.fa, already DOUBLED (the Darevskia
#                            run's subx/); a missing file skips that sub-subfamily
#   members.tsv              family<TAB>sub-subfamily, one row per sub-subfamily
# Tools: samtools, ssearch36 on PATH.
[[ $# -eq 3 ]] || { echo "Usage: $0 <work_dir> <subfamily_consensus_dir> <members.tsv>" >&2; exit 1; }
WORK="$1"; SUBX="$2"; MEMBERS="$3"
MASTER=$WORK/locus_master.tsv          # Lid scf start0 end species family MONO
FA=$WORK/all_regions.fa                 # has LB/RB regions per locus (sp|scf:start-end)
MANIFEST=$WORK/manifest.tsv             # locus kind sp scf rs re trunc explen
command -v samtools >/dev/null || { echo "samtools not in PATH" >&2; exit 1; }
command -v ssearch36 >/dev/null || { echo "ssearch36 not in PATH" >&2; exit 1; }
[[ -d "$SUBX" ]] || { echo "sub-subfamily consensus dir not found: $SUBX" >&2; exit 1; }

cd "$WORK"
samtools faidx "$FA"   # an unindexable regions FASTA is an error (was masked, and every locus became NA)

# Map: which sub-subfamilies belong to each family (space-separated list per family)
declare -A FAMSUBS
while IFS=$'\t' read -r fam sub _; do
  [[ -z "$fam" || "$fam" == \#* ]] && continue
  FAMSUBS[$fam]="${FAMSUBS[$fam]:-}${FAMSUBS[$fam]:+ }$sub"
done < "$MEMBERS"

OUT=$WORK/locus_subfamily.tsv
printf 'locus\tmajor\tsubfamily\tpct_id\taln_len\tbits\tsecond\tsecond_bits\tdelta\n' > "$OUT"

# For each locus, pull its LB region (array start, 2 monomers) and align vs the
# sub-consensuses of its major family. Pick best by bitscore.
tmpq=$(mktemp)
while IFS=$'\t' read -r lid scf s0 e sp major mono; do
  # find LB region header for this locus
  read -r _ kind rsp rscf rs re trunc explen < <(awk -F'\t' -v L="$lid" '$1==L && $2=="LB"{print; exit}' "$MANIFEST")
  if [ -z "${rs:-}" ] || [ "$rs" = "0" ]; then
    # fallback to RB
    read -r _ kind rsp rscf rs re trunc explen < <(awk -F'\t' -v L="$lid" '$1==L && $2=="RB"{print; exit}' "$MANIFEST")
  fi
  [ -z "${rs:-}" ] && { printf '%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' "$lid" "$major" >> "$OUT"; continue; }
  hdr="${rsp}|${rscf}:${rs}-${re}"
  # stage 16 writes every LB/RB region: a missing one is an error, not an NA row
  samtools faidx "$FA" "$hdr" > $tmpq.raw || { echo "FATAL: region $hdr of $lid not in $FA" >&2; exit 1; }
  awk 'NR>1{printf "%s",$0} END{print ""}' $tmpq.raw > $tmpq.seq
  qlen=$(tr -d '\n' < $tmpq.seq | wc -c)
  [ "$qlen" -lt 30 ] && { printf '%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' "$lid" "$major" >> "$OUT"; continue; }
  printf '>q\n%s\n' "$(cat $tmpq.seq)" > $tmpq.fa

  # align vs each candidate sub-consensus, collect best bits
  best=""; bestbits=-1; bestpid=0; bestaln=0; second=""; secondbits=-1
  for sub in ${FAMSUBS[$major]:-}; do
    consx=$SUBX/${sub}.fa
    # every listed sub-subfamily must have its consensus: a missing one would silently narrow the choice
    [ -s "$consx" ] || { echo "FATAL: sub-subfamily $sub of $major (members.tsv) has no consensus file $consx" >&2; exit 1; }
    # awk, not grep -v: no hit must give an empty result, not a pipefail exit
    res=$(ssearch36 -n -z 3 -m 8C -E 1000 $tmpq.fa "$consx" | awk '!/^#/' | sort -k12,12gr | awk 'NR==1')
    [ -z "$res" ] && continue
    bits=$(echo "$res" | cut -f12)
    pid=$(echo "$res" | cut -f3)
    aln=$(echo "$res" | cut -f4)
    # Compare bitscores as floats (the earlier integer truncation let a 0.4-bit-lower
    # candidate win on ties of the integer part; L299/L327/L396/L407 had negative delta).
    if awk -v a="$bits" -v b="$bestbits" 'BEGIN{exit !(a>b)}'; then
      second="$best"; secondbits="$bestbits"
      best="$sub"; bestbits="$bits"; bestpid="$pid"; bestaln="$aln"
    elif awk -v a="$bits" -v b="$secondbits" 'BEGIN{exit !(a>b)}'; then
      second="$sub"; secondbits="$bits"
    fi
  done
  if [ -z "$best" ]; then
    printf '%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' "$lid" "$major" >> "$OUT"
  else
    # delta = best - second (how confident)
    if [ "$secondbits" = "-1" ] || [ -z "$second" ]; then
      delta="NA"; second="NA"; secondbits="NA"
    else
      delta=$(awk -v a="$bestbits" -v b="$secondbits" 'BEGIN{printf "%.1f", a-b}')
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$lid" "$major" "$best" "$bestpid" "$bestaln" "$bestbits" "$second" "$secondbits" "$delta" >> "$OUT"
  fi
done < "$MASTER"
rm -f $tmpq $tmpq.fa $tmpq.seq $tmpq.raw

echo "rows: $(( $(wc -l < "$OUT") - 1 ))"
echo "=== sample ==="
head -8 "$OUT"
echo "=== subfamily distribution ==="
awk -F'\t' 'NR>1{print $3}' "$OUT" | sort | uniq -c | sort -rn
