#!/bin/bash
set -euo pipefail
DIR=/staging/Sergei/Lacertidae/Darevskia/Satellites/manuscript/CLsat_full_clusters_6genome
WORK=$DIR/flank_analysis
GENOMES=$DIR/genomes.tsv
MASTER=$WORK/locus_master.tsv   # Lid scf start0 end species subtype MONO
FLANK=500

mkdir -p $WORK/regions

# species -> genome path
declare -A G
while IFS=$'\t' read -r sp path; do
  [ -z "$sp" ] && continue
  G[$sp]="$path"
done < $GENOMES

# Build a samtools region list per species, plus a manifest mapping region->(locus,kind)
# kinds: LF (left flank), RF (right flank), LB (left boundary 2*MONO), RB (right boundary 2*MONO)
: > $WORK/manifest.tsv   # locus  kind  species  scaffold  reg_start1  reg_end  truncated  expected_len

# per-species region files
declare -A RF_FILE
for sp in "${!G[@]}"; do
  : > $WORK/regions/${sp}.regions.txt
done

# scaffold length lookup per species (from .fai)
get_scaflen() {
  local genome="$1" scf="$2"
  awk -v s="$scf" -F'\t' '$1==s{print $2; exit}' "${genome}.fai"
}

while IFS=$'\t' read -r lid scf s0 e sp subt mono; do
  genome="${G[$sp]}"
  scaflen=$(get_scaflen "$genome" "$scf")
  [ -z "$scaflen" ] && { echo "WARN no scaflen $lid $scf" >&2; continue; }

  # ----- LEFT FLANK: [s0-FLANK, s0]  (1-based: s0-FLANK+1 .. s0) -----
  lf_s1=$(( s0 - FLANK + 1 )); lf_e=$(( s0 ))
  lf_trunc=0
  if (( lf_s1 < 1 )); then lf_s1=1; lf_trunc=1; fi
  if (( lf_e >= 1 && lf_e >= lf_s1 )); then
    echo "${scf}:${lf_s1}-${lf_e}" >> $WORK/regions/${sp}.regions.txt
    printf '%s\tLF\t%s\t%s\t%d\t%d\t%d\t%d\n' "$lid" "$sp" "$scf" "$lf_s1" "$lf_e" "$lf_trunc" "$FLANK" >> $WORK/manifest.tsv
  else
    printf '%s\tLF\t%s\t%s\t%d\t%d\t%d\t%d\n' "$lid" "$sp" "$scf" 0 0 1 "$FLANK" >> $WORK/manifest.tsv
  fi

  # ----- RIGHT FLANK: [e, e+FLANK] (1-based: e+1 .. e+FLANK) -----
  rf_s1=$(( e + 1 )); rf_e=$(( e + FLANK ))
  rf_trunc=0
  if (( rf_e > scaflen )); then rf_e=$scaflen; rf_trunc=1; fi
  if (( rf_s1 <= rf_e && rf_s1 <= scaflen )); then
    echo "${scf}:${rf_s1}-${rf_e}" >> $WORK/regions/${sp}.regions.txt
    printf '%s\tRF\t%s\t%s\t%d\t%d\t%d\t%d\n' "$lid" "$sp" "$scf" "$rf_s1" "$rf_e" "$rf_trunc" "$FLANK" >> $WORK/manifest.tsv
  else
    printf '%s\tRF\t%s\t%s\t%d\t%d\t%d\t%d\n' "$lid" "$sp" "$scf" 0 0 1 "$FLANK" >> $WORK/manifest.tsv
  fi

  # ----- LEFT BOUNDARY: first 2 monomers at the 5' end of the array, inside: [s0+1, s0+2*mono] -----
  span=$(( 2 * mono ))
  lb_s1=$(( s0 + 1 )); lb_e=$(( s0 + span ))
  if (( lb_e > e )); then lb_e=$e; fi          # don't exceed locus end (short locus)
  if (( lb_e > scaflen )); then lb_e=$scaflen; fi
  if (( lb_s1 <= lb_e )); then
    echo "${scf}:${lb_s1}-${lb_e}" >> $WORK/regions/${sp}.regions.txt
    printf '%s\tLB\t%s\t%s\t%d\t%d\t%d\t%d\n' "$lid" "$sp" "$scf" "$lb_s1" "$lb_e" 0 "$span" >> $WORK/manifest.tsv
  fi

  # ----- RIGHT BOUNDARY: last 2 monomers at the end of the array, inside: [e-2*mono, e] -----
  rb_s1=$(( e - span + 1 )); rb_e=$(( e ))
  if (( rb_s1 <= s0 )); then rb_s1=$(( s0 + 1 )); fi
  if (( rb_s1 < 1 )); then rb_s1=1; fi
  if (( rb_s1 <= rb_e )); then
    echo "${scf}:${rb_s1}-${rb_e}" >> $WORK/regions/${sp}.regions.txt
    printf '%s\tRB\t%s\t%s\t%d\t%d\t%d\t%d\n' "$lid" "$sp" "$scf" "$rb_s1" "$rb_e" 0 "$span" >> $WORK/manifest.tsv
  fi

done < $MASTER

echo "manifest rows: $(wc -l < $WORK/manifest.tsv)"
echo "per-species region counts:"
for sp in "${!G[@]}"; do
  echo "  $sp: $(wc -l < $WORK/regions/${sp}.regions.txt)"
done

# ----- Extract sequences per species in one faidx batch each -----
: > $WORK/all_regions.fa
for sp in "${!G[@]}"; do
  genome="${G[$sp]}"
  rf=$WORK/regions/${sp}.regions.txt
  [ -s "$rf" ] || continue
  # tag each header with species to keep regions unique across species
  samtools faidx "$genome" --region-file "$rf" 2>/dev/null \
    | awk -v sp="$sp" '/^>/{print ">"sp"|"substr($0,2); next}{print}' \
    >> $WORK/all_regions.fa
done

echo "extracted FASTA records: $(grep -c '^>' $WORK/all_regions.fa)"
echo "=== sample headers ==="
grep '^>' $WORK/all_regions.fa | head -4
