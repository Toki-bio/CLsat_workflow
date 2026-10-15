#!/bin/bash
set -euo pipefail
DIR=/staging/Sergei/Lacertidae/Darevskia/Satellites/manuscript/CLsat_full_clusters_6genome
WORK=$DIR/flank_analysis
MASTER=$WORK/locus_master.tsv
FA=$WORK/all_regions.fa

# Consensus files (doubled, so a 2-monomer query can map to any phase)
declare -A CONS
CONS[CLsat1]=$DIR/CLsat1.cons.fa
CONS[CLsat2]=$DIR/CLsat2.cons.fa
CONS[CLsat3]=$DIR/CLsat3.cons.fa
CONS[CLsat4]=$DIR/CLsat4.cons.fa

mkdir -p $WORK/cons_x $WORK/bphase
# Build doubled consensus per subtype (so phase can be read 0..MONO-1)
for st in CLsat1 CLsat2 CLsat3 CLsat4; do
  awk '/^>/{next}{gsub(/[ \t\r\n]/,"");s=s $0} END{printf(">%s_x\n%s%s\n","'"$st"'",s,s)}' \
    ${CONS[$st]} > $WORK/cons_x/${st}.x.fa
done

# Index the regions FASTA so we can pull individual records
samtools faidx $FA

# Output: locus kind subtype border_pos_in_consensus aln_len pct_id strand raw_query_len
OUT=$WORK/locus_boundary_phase.tsv
printf 'locus\tkind\tsubtype\tcons_coord\taln_len\tpct_id\tstrand\tquery_len\n' > $OUT

# For each LB/RB manifest row, extract the region seq, align vs doubled consensus,
# report where in the consensus the array boundary monomer maps.
# cons_coord: for LB = consensus position where the array STARTS (5' phase)
#             for RB = consensus position where the array ENDS (3' phase)
tmpq=$(mktemp)
while IFS=$'\t' read -r lid kind sp scf rs re trunc explen; do
  [ "$kind" = "LB" ] || [ "$kind" = "RB" ] || continue
  [ "$rs" = "0" ] && continue
  # subtype + mono for this locus
  read -r _ _ _ _ _ subt mono < <(grep -P "^${lid}\t" $MASTER)
  hdr="${sp}|${scf}:${rs}-${re}"
  # extract region sequence
  samtools faidx $FA "$hdr" 2>/dev/null | awk 'NR>1{printf "%s",$0} END{print ""}' > $tmpq.seq || continue
  qlen=$(tr -d '\n' < $tmpq.seq | wc -c)
  [ "$qlen" -lt 20 ] && continue
  printf '>q\n%s\n' "$(cat $tmpq.seq)" > $tmpq.fa
  consx=$WORK/cons_x/${subt}.x.fa
  # ssearch36: query=region, library=doubled consensus
  # -m 8C tabular: qid sid pid alen mm go qs qe ss se ev bits
  res=$(ssearch36 -n -m 8C -E 1000 $tmpq.fa $consx 2>/dev/null | grep -v '^#' | sort -k12,12gr | head -1)
  [ -z "$res" ] && { printf '%s\t%s\t%s\tNA\tNA\tNA\tNA\t%d\n' "$lid" "$kind" "$subt" "$qlen" >> $OUT; continue; }
  pid=$(echo "$res" | cut -f3)
  alen=$(echo "$res" | cut -f4)
  ss=$(echo "$res" | cut -f9)    # subject (consensus) start
  se=$(echo "$res" | cut -f10)   # subject end
  # strand
  if [ "$se" -ge "$ss" ]; then strand="+"; else strand="-"; fi
  # cons_coord: the consensus position of the array boundary.
  # For LB (5' start of array): the array begins at the query 5' end -> subject coord = ss (mod MONO)
  # For RB (3' end of array):   the array ends at the query 3' end -> subject coord = se (mod MONO)
  if [ "$kind" = "LB" ]; then raw=$ss; else raw=$se; fi
  cc=$(( (raw - 1) % mono + 1 ))   # 1-based coordinate within a single monomer
  printf '%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\n' "$lid" "$kind" "$subt" "$cc" "$alen" "$pid" "$strand" "$qlen" >> $OUT
done < $WORK/manifest.tsv
rm -f $tmpq $tmpq.fa $tmpq.seq

echo "boundary phase rows: $(( $(wc -l < $OUT) - 1 ))"
echo "=== sample ==="
head -8 $OUT
echo "=== LB cons_coord distribution (CLsat1) ==="
awk -F'\t' '$2=="LB" && $3=="CLsat1" && $4!="NA"{print $4}' $OUT | sort -n | uniq -c | sort -rn | head -10
