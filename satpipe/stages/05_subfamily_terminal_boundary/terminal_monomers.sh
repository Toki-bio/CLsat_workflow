#!/bin/bash
set -euo pipefail
# Characterise the VERY FIRST and VERY LAST monomer of each array.
#
# Usage: terminal_monomers.sh <work_dir> <genomes.tsv> <families.tsv>
#   work_dir      holds locus_master.tsv (Lid scf start0 end species family MONO); receives terminal_monomers.tsv
#   genomes.tsv   species_code  genome.fa  (indexed: genome.fa.fai)
#   families.tsv  family  mono  consensus.fa   (single monomer, NOT doubled)
# Environment (defaults = published Darevskia run):
#   WINDOW_MONO   window into the array from each terminus, in monomers   [1.5]
#   EDGE_TOL      tolerated gap (bp) between alignment and array edge     [10]
# Tools: samtools, ssearch36 on PATH.
#
# Method (rewritten 2026-10-15): extract a window spanning WINDOW_MONO monomers INTO the array
# from each terminus and align it against the DOUBLED consensus with ssearch36 -m 8C, keeping the
# single highest-scoring alignment (both strands are searched). The doubled consensus lets one
# local alignment cover the whole window without a phase break, so the consensus coordinate at
# the array edge is read directly from the alignment and gives the phase of the terminal monomer.
# The strand of the hit is recorded and used to decide which consensus coordinates the terminal
# monomer spans.
#
# The previous implementation aligned the SINGLE consensus and chose among the (at most two,
# fwd/rc) rows by subject coordinate; ssearch36 reports one alignment per library sequence per
# strand, so that choice picked the best-scoring monomer inside the window, usually the
# complete internal one, not the terminal one, and discarded strand.
#
# Columns:
#   cons_start, cons_end : consensus coordinates (1-based, single monomer frame) spanned by the
#                          terminal monomer, in the consensus reading orientation
#   covered = cons_end-cons_start+1
#   partial = (covered < 0.9*MONO) ? yes : no
#   strand  : + if the array reads in consensus orientation along the scaffold, - otherwise
#   edge_gap: bases between the array edge and the start of the alignment (0 = alignment
#             reaches the edge). Rows with edge_gap > EDGE_TOL are reported as NA.
[[ $# -eq 3 ]] || { echo "Usage: $0 <work_dir> <genomes.tsv> <families.tsv>" >&2; exit 1; }
WORK="$1"; GENOMES="$2"; FAMTAB="$3"
MASTER="$WORK/locus_master.tsv"
WINDOW_MONO="${WINDOW_MONO:-1.5}"
EDGE_TOL="${EDGE_TOL:-10}"
command -v samtools >/dev/null || { echo "samtools not in PATH" >&2; exit 1; }
command -v ssearch36 >/dev/null || { echo "ssearch36 not in PATH" >&2; exit 1; }
[[ -s "$MASTER" ]] || { echo "missing $MASTER" >&2; exit 1; }

cd "$WORK"

# Single-monomer consensus per family (NOT doubled)
declare -A CONS
FIRSTFAM=""
while IFS=$'\t' read -r fam mono cons _; do
  [[ -z "$fam" || "$fam" == \#* ]] && continue
  CONS[$fam]="$cons"; [[ -z "$FIRSTFAM" ]] && FIRSTFAM="$fam"
done < "$FAMTAB"

declare -A G
while IFS=$'\t' read -r sp path _; do [[ -z "$sp" || "$sp" == \#* ]] && continue; G[$sp]="$path"; done < "$GENOMES"

get_scaflen(){ awk -v s="$2" -F'\t' '$1==s{print $2;exit}' "${1}.fai"; }

OUT=$WORK/terminal_monomers.tsv
printf 'locus\tsubtype\tmono\tend\tcons_start\tcons_end\tcovered\tfrac\tpartial\tpct_id\tstrand\tedge_gap\n' > "$OUT"

tmpq=$(mktemp)
# doubled consensuses, built once
declare -A CONS2
for fam in "${!CONS[@]}"; do
  CONS2[$fam]=$tmpq.$fam.x2.fa
  awk -v n="${fam}_x2" '/^>/{next}{s=s $0}END{print ">" n; print s s}' "${CONS[$fam]}" > "${CONS2[$fam]}"
done

while IFS=$'\t' read -r lid scf s0 e sp subt mono; do
  genome="${G[$sp]:-}"
  consx="${CONS2[$subt]:-}"
  [ -z "$genome" ] || [ -z "$consx" ] && { continue; }
  scaflen=$(get_scaflen "$genome" "$scf")
  [ -z "$scaflen" ] && continue
  span=$(awk -v m="$mono" -v w="$WINDOW_MONO" 'BEGIN{print int(m*w)}')   # WINDOW_MONO monomers

  for END in 5p 3p; do
    if [ "$END" = "5p" ]; then
      r_s=$(( s0 + 1 )); r_e=$(( s0 + span ))
    else
      r_s=$(( e - span + 1 )); r_e=$(( e ))
    fi
    (( r_s < 1 )) && r_s=1
    (( r_e > scaflen )) && r_e=$scaflen
    (( r_s > r_e )) && continue
    samtools faidx "$genome" "${scf}:${r_s}-${r_e}" 2>/dev/null | awk 'NR>1{printf "%s",$0}END{print ""}' > $tmpq.seq
    qlen=$(tr -d '\n' < $tmpq.seq | wc -c)
    [ "$qlen" -lt 20 ] && continue
    printf '>q\n%s\n' "$(cat $tmpq.seq)" > $tmpq.fa
    # best alignment of the window against the doubled consensus, either strand
    res=$(ssearch36 -n -m 8C -E 100 $tmpq.fa $consx 2>/dev/null | grep -v '^#' | sort -k12,12gr | head -1)
    if [ -z "$res" ]; then
      printf '%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' "$lid" "$subt" "$mono" "$END" >> "$OUT"
      continue
    fi
    echo "$res" | awk -F'\t' -v lid="$lid" -v subt="$subt" -v mono="$mono" -v which="$END" -v qlen="$qlen" -v tol="$EDGE_TOL" 'BEGIN{OFS="\t"}{
      pid=$3; q7=$7; q8=$8; s9=$9; s10=$10
      strand=(s9<=s10)?"+":"-"
      if (which=="5p") { gap=q7-1 } else { gap=qlen-q8 }
      if (gap>tol) { print lid,subt,mono,which,"NA","NA","NA","NA","NA",pid,strand,gap; next }
      # consensus position (doubled frame) at the array edge, extrapolated over the gap
      if (which=="5p") { pe=(strand=="+")?s9-gap:s9+gap } else { pe=(strand=="+")?s10+gap:s10-gap }
      p1=((pe-1)%mono+mono)%mono+1          # fold into the single-monomer frame, 1..mono
      # terminal monomer, in consensus reading orientation:
      #   array start, + strand : first monomer runs from p1 to mono
      #   array start, - strand : the array reads backwards, so this edge is the END of the
      #                           last monomer in reading orientation: it runs from 1 to p1
      #   array end,   + strand : last monomer runs from 1 to p1
      #   array end,   - strand : first monomer in reading orientation runs from p1 to mono
      if ((which=="5p" && strand=="+") || (which=="3p" && strand=="-")) { cs=p1; ce=mono } else { cs=1; ce=p1 }
      cov=ce-cs+1
      printf "%s\t%s\t%d\t%s\t%d\t%d\t%d\t%.2f\t%s\t%s\t%s\t%d\n", lid,subt,mono,which,cs,ce,cov,cov/mono,(cov<0.9*mono)?"yes":"no",pid,strand,gap
    }' >> "$OUT"
  done
done < "$MASTER"
rm -f $tmpq $tmpq.fa $tmpq.seq "${CONS2[@]}"

echo "rows: $(( $(wc -l < "$OUT") - 1 ))"
echo "=== sample ==="
head -10 "$OUT"
echo "=== 5p partial fraction distribution (how often first monomer is incomplete) ==="
awk -F'\t' '$4=="5p" && $9!="NA"{print $9}' "$OUT" | sort | uniq -c
echo "=== 5p covered-length histogram ($FIRSTFAM) ==="
awk -F'\t' -v F="$FIRSTFAM" '$4=="5p" && $2==F && $7!="NA"{b=int($7/15)*15; print b"-"b+14}' "$OUT" | sort -t- -k1n | uniq -c
