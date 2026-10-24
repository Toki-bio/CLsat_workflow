#!/usr/bin/env bash
set -euo pipefail
# Extract the regions the sub-subfamily / terminal / boundary-phase stages read:
#   LF/RF  REGION_FLANK bp outside each array edge
#   LB/RB  2 x MONO bp inside each array edge
# into <work_dir>/manifest.tsv (locus kind sp scf rs re trunc explen; 1-based inclusive) and
# <work_dir>/all_regions.fa (header sp|scf:start-end, indexed).
#
# Usage: run_extract.sh <work_dir> <genomes.tsv>
#   work_dir     holds locus_master.tsv (Lid scf start0 end species family MONO)
#   genomes.tsv  species_code  genome.fa (indexed)
# Environment: REGION_FLANK  outside-flank length, bp  [500]
# Tools: samtools on PATH.
[[ $# -eq 2 ]] || { echo "Usage: $0 <work_dir> <genomes.tsv>" >&2; exit 1; }
WORK="$1"; GENOMES="$2"
MASTER=$WORK/locus_master.tsv
FLANK="${REGION_FLANK:-500}"
command -v samtools >/dev/null || { echo "samtools not in PATH" >&2; exit 1; }
[[ -s "$MASTER" ]] || { echo "missing $MASTER" >&2; exit 1; }

mkdir -p "$WORK/regions"

declare -A G
while IFS=$'\t' read -r sp path _; do
  [[ -z "$sp" || "$sp" == \#* ]] && continue
  G[$sp]="$path"
done < "$GENOMES"

: > "$WORK/manifest.tsv"

get_scaflen() {
  local genome="$1" scf="$2"
  awk -v s="$scf" -F'\t' '$1==s{print $2; exit}' "${genome}.fai"
}

declare -A RSTR
for sp in "${!G[@]}"; do
  RSTR[$sp]=""
done

while IFS=$'\t' read -r lid scf s0 e sp subt mono; do
  genome="${G[$sp]:-}"
  # a locus whose genome, contig or coordinates do not fit its genome is a catalogue error (once: a locus of one
  # genome filed under another that shares the contig name); stop instead of silently skipping its regions
  [ -z "$genome" ] && { echo "FATAL: genome $sp of $lid not in $GENOMES" >&2; exit 1; }
  scaflen=$(get_scaflen "$genome" "$scf")
  [ -z "$scaflen" ] && { echo "FATAL: contig $scf of $lid not in genome $sp ($genome)" >&2; exit 1; }
  [ "$e" -gt "$scaflen" ] && { echo "FATAL: $lid ends at $e beyond contig $scf of genome $sp (length $scaflen)" >&2; exit 1; }

  lf_s1=$(( s0 - FLANK + 1 )); lf_e=$(( s0 ))
  lf_trunc=0
  [ "$lf_s1" -lt 1 ] && { lf_s1=1; lf_trunc=1; }
  if [ "$lf_e" -ge "$lf_s1" ]; then
    printf '%s\tLF\t%s\t%s\t%s\t%s\t%s\t%s\n' "$lid" "$sp" "$scf" "$lf_s1" "$lf_e" "$lf_trunc" "$FLANK" >> "$WORK/manifest.tsv"
    RSTR[$sp]="${RSTR[$sp]}${scf}:${lf_s1}-${lf_e}"$'\n'
  fi

  rf_s1=$(( e + 1 )); rf_e=$(( e + FLANK ))
  rf_trunc=0
  [ "$rf_e" -gt "$scaflen" ] && { rf_e=$scaflen; rf_trunc=1; }
  if [ "$rf_e" -ge "$rf_s1" ]; then
    printf '%s\tRF\t%s\t%s\t%s\t%s\t%s\t%s\n' "$lid" "$sp" "$scf" "$rf_s1" "$rf_e" "$rf_trunc" "$FLANK" >> "$WORK/manifest.tsv"
    RSTR[$sp]="${RSTR[$sp]}${scf}:${rf_s1}-${rf_e}"$'\n'
  fi

  bwin=$(( 2 * mono ))
  lb_s1=$(( s0 + 1 )); lb_e=$(( s0 + bwin ))
  lb_trunc=0
  [ "$lb_e" -gt "$e" ] && { lb_e=$e; lb_trunc=1; }
  if [ "$lb_e" -ge "$lb_s1" ]; then
    printf '%s\tLB\t%s\t%s\t%s\t%s\t%s\t%s\n' "$lid" "$sp" "$scf" "$lb_s1" "$lb_e" "$lb_trunc" "$bwin" >> "$WORK/manifest.tsv"
    RSTR[$sp]="${RSTR[$sp]}${scf}:${lb_s1}-${lb_e}"$'\n'
  fi

  rb_e=$e; rb_s1=$(( e - bwin + 1 ))
  rb_trunc=0
  [ "$rb_s1" -lt "$((s0+1))" ] && { rb_s1=$((s0+1)); rb_trunc=1; }
  if [ "$rb_e" -ge "$rb_s1" ]; then
    printf '%s\tRB\t%s\t%s\t%s\t%s\t%s\t%s\n' "$lid" "$sp" "$scf" "$rb_s1" "$rb_e" "$rb_trunc" "$bwin" >> "$WORK/manifest.tsv"
    RSTR[$sp]="${RSTR[$sp]}${scf}:${rb_s1}-${rb_e}"$'\n'
  fi
done < "$MASTER"

echo "manifest rows: $(wc -l < "$WORK/manifest.tsv")"

# Extract all regions per species into one FASTA (all_regions.fa), header sp|scf:start-end
: > "$WORK/all_regions.fa"
for sp in "${!G[@]}"; do
  [ -z "${RSTR[$sp]}" ] && continue
  printf '%s' "${RSTR[$sp]}" | sort -u > "$WORK/regions/${sp}.regions.txt"
  genome="${G[$sp]}"
  while read -r reg; do
    [ -z "$reg" ] && continue
    samtools faidx "$genome" "$reg" 2>/dev/null | \
      awk -v sp="$sp" 'NR==1{print ">"sp"|"substr($0,2); next}{print}'
  done < "$WORK/regions/${sp}.regions.txt" >> "$WORK/all_regions.fa"
done
samtools faidx "$WORK/all_regions.fa"
echo "all_regions.fa records: $(grep -c '^>' "$WORK/all_regions.fa")"

touch "$WORK/extract.done"
echo EXTRACT_DONE
