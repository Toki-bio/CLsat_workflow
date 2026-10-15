#!/bin/bash
set -euo pipefail
source /staging/conda/etc/profile.d/conda.sh 2>/dev/null || true
conda activate bioinfo 2>/dev/null || export PATH=/staging/conda/envs/bioinfo/bin:$PATH

DIR=/staging/Sergei/Lacertidae/Darevskia/Satellites/manuscript/CLsat_full_clusters_6genome
WORK=$DIR/flank_analysis
MASTER=$WORK/locus_master.tsv          # Lid scf start0 end species subtype MONO
SUBCONS=$WORK/subfamily_consensus.fa
FA=$WORK/all_regions.fa                 # has LB/RB regions per locus (sp|scf:start-end)
MANIFEST=$WORK/manifest.tsv             # locus kind sp scf rs re trunc explen

cd $WORK
samtools faidx $FA 2>/dev/null || true

# Build doubled sub-consensus (so a 2-monomer query aligns cleanly at any phase)
awk '/^>/{name=substr($0,2); next}{print ">"name"\n"$0$0}' $SUBCONS > $WORK/subfamily_consensus.x.fa

# Map: which sub-subfamilies belong to each major family
# CLsat1 -> CLsat1.1..1.5 ; CLsat2 -> CLsat2 ; CLsat3 -> CLsat3.1,3.2 ; CLsat4 -> CLsat4
declare -A FAMSUBS
FAMSUBS[CLsat1]="CLsat1.1 CLsat1.2 CLsat1.3 CLsat1.4 CLsat1.5"
FAMSUBS[CLsat2]="CLsat2"
FAMSUBS[CLsat3]="CLsat3.1 CLsat3.2"
FAMSUBS[CLsat4]="CLsat4"

# Extract each sub-consensus (doubled) into its own file for targeted search

OUT=$WORK/locus_subfamily.tsv
printf 'locus\tmajor\tsubfamily\tpct_id\taln_len\tbits\tsecond\tsecond_bits\tdelta\n' > $OUT

# For each locus, pull its LB region (array start, 2 monomers) and align vs the
# sub-consensuses of its major family. Pick best by bitscore.
tmpq=$(mktemp)
while IFS=$'\t' read -r lid scf s0 e sp major mono; do
  # find LB region header for this locus
  read -r _ kind rsp rscf rs re trunc explen < <(awk -F'\t' -v L="$lid" '$1==L && $2=="LB"{print; exit}' $MANIFEST)
  if [ -z "${rs:-}" ] || [ "$rs" = "0" ]; then
    # fallback to RB
    read -r _ kind rsp rscf rs re trunc explen < <(awk -F'\t' -v L="$lid" '$1==L && $2=="RB"{print; exit}' $MANIFEST)
  fi
  [ -z "${rs:-}" ] && { printf '%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' "$lid" "$major" >> $OUT; continue; }
  hdr="${rsp}|${rscf}:${rs}-${re}"
  samtools faidx $FA "$hdr" 2>/dev/null | awk 'NR>1{printf "%s",$0} END{print ""}' > $tmpq.seq
  qlen=$(tr -d '\n' < $tmpq.seq | wc -c)
  [ "$qlen" -lt 30 ] && { printf '%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' "$lid" "$major" >> $OUT; continue; }
  printf '>q\n%s\n' "$(cat $tmpq.seq)" > $tmpq.fa

  # align vs each candidate sub-consensus, collect best bits
  best=""; bestbits=-1; bestpid=0; bestaln=0; second=""; secondbits=-1
  for sub in ${FAMSUBS[$major]}; do
    consx=$WORK/subx/${sub}.fa
    [ -s "$consx" ] || continue
    res=$(ssearch36 -n -m 8C -E 1000 $tmpq.fa $consx 2>/dev/null | grep -v '^#' | sort -k12,12gr | head -1)
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
    printf '%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' "$lid" "$major" >> $OUT
  else
    # delta = best - second (how confident)
    if [ "$secondbits" = "-1" ] || [ -z "$second" ]; then
      delta="NA"; second="NA"; secondbits="NA"
    else
      delta=$(awk -v a="$bestbits" -v b="$secondbits" 'BEGIN{printf "%.1f", a-b}')
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$lid" "$major" "$best" "$bestpid" "$bestaln" "$bestbits" "$second" "$secondbits" "$delta" >> $OUT
  fi
done < $MASTER
rm -f $tmpq $tmpq.fa $tmpq.seq

echo "rows: $(( $(wc -l < $OUT) - 1 ))"
echo "=== sample ==="
head -8 $OUT
echo "=== subfamily distribution ==="
awk -F'\t' 'NR>1{print $3}' $OUT | sort | uniq -c | sort -rn
