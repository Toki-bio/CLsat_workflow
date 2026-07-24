#!/bin/bash
set -euo pipefail
source /staging/conda/etc/profile.d/conda.sh 2>/dev/null || true
conda activate bioinfo 2>/dev/null || export PATH=/staging/conda/envs/bioinfo/bin:$PATH

DIR=/staging/Sergei/Lacertidae/Darevskia/Satellites/manuscript/CLsat_full_clusters_6genome
WORK=$DIR/flank_analysis
MASTER=$WORK/locus_master.tsv     # Lid scf start0 end species subtype MONO
GENOMES=$DIR/genomes.tsv

cd $WORK

# Single-monomer consensus per family (NOT doubled) - to find phase of terminal monomer
declare -A CONS
CONS[CLsat1]=$DIR/CLsat1.cons.fa
CONS[CLsat2]=$DIR/CLsat2.cons.fa
CONS[CLsat3]=$DIR/CLsat3.cons.fa
CONS[CLsat4]=$DIR/CLsat4.cons.fa

declare -A G
while IFS=$'\t' read -r sp path; do [ -z "$sp" ] && continue; G[$sp]="$path"; done < $GENOMES

get_scaflen(){ awk -v s="$2" -F'\t' '$1==s{print $2;exit}' "${1}.fai"; }

# Output: characterize the VERY FIRST and VERY LAST monomer of each array.
# Method: extract a window from the array boundary spanning ~1.5 monomers INTO
# the array from each terminus, align a single consensus monomer (m8C), and read:
#   - cons_start, cons_end : which consensus coords the terminal monomer spans
#   - covered = cons_end-cons_start+1 : how many consensus bp the terminal monomer covers
#   - partial = (covered < 0.9*MONO) ? yes : no  (incomplete terminal monomer)
# For the 5' terminus we look at the FIRST monomer; for 3' the LAST monomer.

OUT=$WORK/terminal_monomers.tsv
printf 'locus\tsubtype\tmono\tend\tcons_start\tcons_end\tcovered\tfrac\tpartial\tpct_id\n' > $OUT

tmpq=$(mktemp)
while IFS=$'\t' read -r lid scf s0 e sp subt mono; do
  genome="${G[$sp]}"
  consx="${CONS[$subt]}"
  [ -z "$genome" ] || [ -z "$consx" ] && { continue; }
  scaflen=$(get_scaflen "$genome" "$scf")
  [ -z "$scaflen" ] && continue
  span=$(( mono + mono/2 ))   # 1.5 monomers

  for END in 5p 3p; do
    if [ "$END" = "5p" ]; then
      # first 1.5 monomers of the array (genomic start)
      r_s=$(( s0 + 1 )); r_e=$(( s0 + span ))
    else
      # last 1.5 monomers of the array (genomic end)
      r_s=$(( e - span + 1 )); r_e=$(( e ))
    fi
    (( r_s < 1 )) && r_s=1
    (( r_e > scaflen )) && r_e=$scaflen
    (( r_s > r_e )) && continue
    samtools faidx "$genome" "${scf}:${r_s}-${r_e}" 2>/dev/null | awk 'NR>1{printf "%s",$0}END{print ""}' > $tmpq.seq
    qlen=$(tr -d '\n' < $tmpq.seq | wc -c)
    [ "$qlen" -lt 20 ] && continue
    printf '>q\n%s\n' "$(cat $tmpq.seq)" > $tmpq.fa
    # Align single consensus monomer; take the alignment whose subject (consensus)
    # coords reveal the phase. For 5p we want the FIRST aligned block (smallest qstart);
    # for 3p the LAST aligned block (largest qend).
    if [ "$END" = "5p" ]; then
      res=$(ssearch36 -n -m 8C -E 100 $tmpq.fa $consx 2>/dev/null | grep -v '^#' | sort -k7,7n | head -1)
    else
      res=$(ssearch36 -n -m 8C -E 100 $tmpq.fa $consx 2>/dev/null | grep -v '^#' | sort -k8,8nr | head -1)
    fi
    if [ -z "$res" ]; then
      printf '%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\n' "$lid" "$subt" "$mono" "$END" >> $OUT
      continue
    fi
    pid=$(echo "$res" | cut -f3)
    ss=$(echo "$res" | cut -f9)    # consensus start
    se=$(echo "$res" | cut -f10)   # consensus end
    # normalize orientation
    cs=$ss; ce=$se
    if [ "$se" -lt "$ss" ]; then cs=$se; ce=$ss; fi
    covered=$(( ce - cs + 1 ))
    frac=$(awk -v c="$covered" -v m="$mono" 'BEGIN{printf "%.2f", c/m}')
    partial=$(awk -v c="$covered" -v m="$mono" 'BEGIN{print (c < 0.9*m)?"yes":"no"}')
    printf '%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\n' "$lid" "$subt" "$mono" "$END" "$cs" "$ce" "$covered" "$frac" "$partial" "$pid" >> $OUT
  done
done < $MASTER
rm -f $tmpq $tmpq.fa $tmpq.seq

echo "rows: $(( $(wc -l < $OUT) - 1 ))"
echo "=== sample ==="
head -10 $OUT
echo "=== 5p partial fraction distribution (how often first monomer is incomplete) ==="
awk -F'\t' '$4=="5p" && $9!="NA"{print $9}' $OUT | sort | uniq -c
echo "=== 5p covered-length histogram (CLsat1) ==="
awk -F'\t' '$4=="5p" && $2=="CLsat1" && $7!="NA"{b=int($7/15)*15; print b"-"b+14}' $OUT | sort -t- -k1n | uniq -c
