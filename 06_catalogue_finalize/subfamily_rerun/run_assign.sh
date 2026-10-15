#!/usr/bin/env bash
set -euo pipefail
source /staging/conda/etc/profile.d/conda.sh 2>/dev/null || true
conda activate bioinfo 2>/dev/null || export PATH=/staging/conda/envs/bioinfo/bin:$PATH

OLD=/staging/Sergei/Lacertidae/Darevskia/Satellites/manuscript/CLsat_full_clusters_6genome/flank_analysis
WORK=/staging/tmp/clsat_v6/subfamily
MASTER=$WORK/locus_master.tsv
FA=$WORK/all_regions.fa
MANIFEST=$WORK/manifest.tsv

cd $WORK
cp -r $OLD/subx $WORK/subx

declare -A FAMSUBS
FAMSUBS[CLsat1]="CLsat1.1 CLsat1.2 CLsat1.3 CLsat1.4 CLsat1.5"
FAMSUBS[CLsat2]="CLsat2"
FAMSUBS[CLsat3]="CLsat3.1 CLsat3.2"
FAMSUBS[CLsat4]="CLsat4"

OUT=$WORK/locus_subfamily.tsv
printf 'locus\tmajor\tsubfamily\tpct_id\taln_len\tbits\tsecond\tsecond_bits\tdelta\n' > $OUT

tmpq=$(mktemp)
n=0
while IFS=$'\t' read -r lid scf s0 e sp major mono; do
  n=$((n+1))
  read -r _ kind rsp rscf rs re trunc explen < <(awk -F'\t' -v L="$lid" '$1==L && $2=="LB"{print; exit}' $MANIFEST)
  if [ -z "${rs:-}" ] || [ "$rs" = "0" ]; then
    read -r _ kind rsp rscf rs re trunc explen < <(awk -F'\t' -v L="$lid" '$1==L && $2=="RB"{print; exit}' $MANIFEST)
  fi
  if [ -z "${rs:-}" ]; then
    printf '%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' "$lid" "$major" >> $OUT
    continue
  fi
  hdr="${rsp}|${rscf}:${rs}-${re}"
  samtools faidx $FA "$hdr" 2>/dev/null | awk 'NR>1{printf "%s",$0} END{print ""}' > $tmpq.seq
  qlen=$(tr -d '\n' < $tmpq.seq | wc -c)
  if [ "$qlen" -lt 30 ]; then
    printf '%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' "$lid" "$major" >> $OUT
    continue
  fi
  printf '>q\n%s\n' "$(cat $tmpq.seq)" > $tmpq.fa

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
  delta="NA"
  if [ "$bestbits" != "-1" ] && [ -n "$secondbits" ] && [ "$secondbits" != "-1" ]; then
    delta=$(echo "$bestbits - $secondbits" | bc 2>/dev/null || echo NA)
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$lid" "$major" "${best:-NA}" "${bestpid:-NA}" "${bestaln:-NA}" "${bestbits:-NA}" \
    "${second:-NA}" "${secondbits:-NA}" "$delta" >> $OUT
  [ $((n % 100)) -eq 0 ] && echo "  processed $n loci..." >&2
done < $MASTER

touch $WORK/subfamily.done
echo ASSIGN_DONE
