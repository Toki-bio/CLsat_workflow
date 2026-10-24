#!/usr/bin/env bash
set -euo pipefail
# Sub-subfamily assignment on a corrected catalogue; same method as
# 05_subfamily_terminal_boundary/assign_subfamily.sh (2-monomer LB window vs doubled
# sub-subfamily consensuses, floats compared), with a progress counter and `bc` delta.
#
# Usage: run_assign.sh <work_dir> <subfamily_consensus_dir> <members.tsv>
#   work_dir                 holds locus_master.tsv, manifest.tsv, all_regions.fa (run_extract.sh);
#                            receives locus_subfamily.tsv
#   subfamily_consensus_dir  one DOUBLED FASTA per sub-subfamily, <name>.fa (subx/)
#   members.tsv              family<TAB>sub-subfamily
# Tools: samtools, ssearch36, bc on PATH.
[[ $# -eq 3 ]] || { echo "Usage: $0 <work_dir> <subfamily_consensus_dir> <members.tsv>" >&2; exit 1; }
WORK="$1"; SUBX="$2"; MEMBERS="$3"
MASTER=$WORK/locus_master.tsv
FA=$WORK/all_regions.fa
MANIFEST=$WORK/manifest.tsv
command -v samtools >/dev/null || { echo "samtools not in PATH" >&2; exit 1; }
command -v ssearch36 >/dev/null || { echo "ssearch36 not in PATH" >&2; exit 1; }
[[ -d "$SUBX" ]] || { echo "sub-subfamily consensus dir not found: $SUBX" >&2; exit 1; }

cd "$WORK"

declare -A FAMSUBS
while IFS=$'\t' read -r fam sub _; do
  [[ -z "$fam" || "$fam" == \#* ]] && continue
  FAMSUBS[$fam]="${FAMSUBS[$fam]:-}${FAMSUBS[$fam]:+ }$sub"
done < "$MEMBERS"

OUT=$WORK/locus_subfamily.tsv
printf 'locus\tmajor\tsubfamily\tpct_id\taln_len\tbits\tsecond\tsecond_bits\tdelta\n' > "$OUT"

tmpq=$(mktemp)
n=0
while IFS=$'\t' read -r lid scf s0 e sp major mono; do
  n=$((n+1))
  read -r _ kind rsp rscf rs re trunc explen < <(awk -F'\t' -v L="$lid" '$1==L && $2=="LB"{print; exit}' "$MANIFEST")
  if [ -z "${rs:-}" ] || [ "$rs" = "0" ]; then
    read -r _ kind rsp rscf rs re trunc explen < <(awk -F'\t' -v L="$lid" '$1==L && $2=="RB"{print; exit}' "$MANIFEST")
  fi
  if [ -z "${rs:-}" ]; then
    printf '%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' "$lid" "$major" >> "$OUT"
    continue
  fi
  hdr="${rsp}|${rscf}:${rs}-${re}"
  samtools faidx "$FA" "$hdr" 2>/dev/null | awk 'NR>1{printf "%s",$0} END{print ""}' > $tmpq.seq
  qlen=$(tr -d '\n' < $tmpq.seq | wc -c)
  if [ "$qlen" -lt 30 ]; then
    printf '%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' "$lid" "$major" >> "$OUT"
    continue
  fi
  printf '>q\n%s\n' "$(cat $tmpq.seq)" > $tmpq.fa

  best=""; bestbits=-1; bestpid=0; bestaln=0; second=""; secondbits=-1
  for sub in ${FAMSUBS[$major]:-}; do
    consx=$SUBX/${sub}.fa
    [ -s "$consx" ] || continue
    res=$(ssearch36 -n -z 3 -m 8C -E 1000 $tmpq.fa "$consx" 2>/dev/null | grep -v '^#' | sort -k12,12gr | head -1)
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
  delta="NA"
  if [ "$bestbits" != "-1" ] && [ -n "$secondbits" ] && [ "$secondbits" != "-1" ]; then
    delta=$(echo "$bestbits - $secondbits" | bc 2>/dev/null || echo NA)
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$lid" "$major" "${best:-NA}" "${bestpid:-NA}" "${bestaln:-NA}" "${bestbits:-NA}" \
    "${second:-NA}" "${secondbits:-NA}" "$delta" >> "$OUT"
  [ $((n % 100)) -eq 0 ] && echo "  processed $n loci..." >&2
done < "$MASTER"

touch "$WORK/subfamily.done"
echo ASSIGN_DONE
