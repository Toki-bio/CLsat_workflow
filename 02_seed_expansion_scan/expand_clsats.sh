#!/usr/bin/env bash
set -euo pipefail

# expand_clsats.sh
# Expand tandem-repeat loci in BED to "true" ends by walking in 3-monomer steps
# and scoring each window against a repeated consensus with ssearch36.
#
# Inputs:
#   1) initial.bed   (BED3+; columns: chrom start0 end0 [name ...])
#   2) genome.fa     (FASTA indexed by samtools faidx; .fai must exist or will be created)
#   3) consensus.fa  (single-seq FASTA recommended)
#   4) monomer_len   (e.g. 146)
#   5) out_prefix    (e.g. out/dar_CLsat1)
#
# Outputs:
#   OUT.stop.tsv
#   OUT.expanded_per_locus.bed
#   OUT.expanded_per_locus.tsv

die(){ echo "ERROR: $*" >&2; exit 1; }
log(){ echo "[expand_clsats] $*" >&2; }

usage(){
  cat >&2 <<'EOF'
Usage:
  expand_clsats.sh initial.bed genome.fa consensus.fa MONO out_prefix [options]

Options:
  --step-monomers N     (default 3)
  --win-monomers  N     (default 3)
  --query-repeat  N     (default 2)   # consensus is repeated N times for search target
  --min-aln-frac  F     (default 1.0) # fraction of (query_repeat*MONO) for min alignment length
  --min-bits      X     (default 50)
  --max-steps     N     (default 500)
  --to-scaffold-ends    (default off) # if on: clamp to edges and keep walking until edge hit
  --progress-every N    (default 50)

Notes:
  - If a BED locus scaffold is missing from genome.fa.fai, it is recorded with reason=scaffold_missing
    and the locus is kept unmodified for that end (expanded_bp=0).
EOF
}

# -------------------------
# Parse args
# -------------------------
[[ $# -lt 5 ]] && usage && exit 1
BED=$1
GENOME=$2
CONS=$3
MONO=$4
OUT=$5
shift 5

STEP_MON=3
WIN_MON=3
QREP=2
MIN_ALN_FRAC=1.0
MIN_BITS=50
MAX_STEPS=500
TO_SCAFF_ENDS=0
PROG_EVERY=50

while [[ $# -gt 0 ]]; do
  case "$1" in
    --step-monomers) STEP_MON=$2; shift 2;;
    --win-monomers)  WIN_MON=$2; shift 2;;
    --query-repeat)  QREP=$2; shift 2;;
    --min-aln-frac)  MIN_ALN_FRAC=$2; shift 2;;
    --min-bits)      MIN_BITS=$2; shift 2;;
    --max-steps)     MAX_STEPS=$2; shift 2;;
    --to-scaffold-ends) TO_SCAFF_ENDS=1; shift;;
    --progress-every) PROG_EVERY=$2; shift 2;;
    -h|--help) usage; exit 0;;
    *) die "Unknown option: $1";;
  esac
done

# -------------------------
# Validate inputs
# -------------------------
[[ -s "$BED" ]]    || die "bed missing/empty: $BED"
[[ -s "$GENOME" ]] || die "genome.fa missing/empty: $GENOME"
[[ -s "$CONS" ]]   || die "consensus.fa missing/empty: $CONS"
[[ "$MONO" =~ ^[0-9]+$ ]] || die "MONO must be integer, got: $MONO"

command -v samtools >/dev/null 2>&1 || die "samtools not found"
command -v ssearch36 >/dev/null 2>&1 || die "ssearch36 not found"
command -v awk >/dev/null 2>&1 || die "awk not found"

OUTDIR=$(dirname "$OUT")
mkdir -p "$OUTDIR"

# Ensure genome index exists
if [[ ! -s "${GENOME}.fai" ]]; then
  log "Indexing genome with samtools faidx: $GENOME"
  samtools faidx "$GENOME"
fi
FAI="${GENOME}.fai"

STEP_BP=$(( STEP_MON * MONO ))
WIN_BP=$((  WIN_MON * MONO ))
QLEN=$(( QREP * MONO ))
# min alignment length (integer) = floor( MIN_ALN_FRAC * QLEN )
MIN_ALN=$(awk -v q="$QLEN" -v f="$MIN_ALN_FRAC" 'BEGIN{m=int(q*f); if(m<1)m=1; print m}')

log "Inputs:"
log "  BED=$BED"
log "  GENOME=$GENOME"
log "  CONS=$CONS"
log "  MONO=$MONO  step_bp=$STEP_BP  win_bp=$WIN_BP  min_aln=$MIN_ALN"
log "  min_bits=$MIN_BITS  to_scaffold_ends=$TO_SCAFF_ENDS  max_steps=$MAX_STEPS  progress_every=$PROG_EVERY"
log "  OUT=$OUT"

# Count loci
NLOCI=$(awk 'BEGIN{FS="[ \t]+"} NF>=3 && $1!~/^#/ && $1!="chrom" && $1!="CHROM"{n++} END{print n+0}' "$BED")
log "Parsed loci from BED: $NLOCI"

STOP_TSV="${OUT}.stop.tsv"
OUT_BED="${OUT}.expanded_per_locus.bed"
OUT_TSV="${OUT}.expanded_per_locus.tsv"

# Headers
printf "chrom\tseed_start0\tseed_end0\tname\tend\tDIR\tstop_step\twin_s1\twin_e1\treason\tclean_pos1\torig_endpos1\texpanded_bp\twalked_bp\tbest_bits\tbest_alen\tbest_pid\tbest_sstrand\n" > "$STOP_TSV"
printf "chrom\tstart0\tend0\tname\tseed_start0\tseed_end0\treason_L\treason_R\texpand_L_bp\texpand_R_bp\texpand_total_bp\twalk_L_bp\twalk_R_bp\twalk_total_bp\n" > "$OUT_TSV"
: > "$OUT_BED"

# Build repeated-consensus target once (CONS repeated QREP+1 times is fine; using QREP+1 gives a little slack)
TMPROOT=$(mktemp -d -p "${OUTDIR}" "tmp.expand_clsats.XXXXXXXX") || exit 1
trap 'rm -rf "$TMPROOT"' EXIT

CONS_X=$(awk -v R=$((QREP+1)) '
  BEGIN{seq=""}
  /^>/{next}
  {gsub(/[ \t\r\n]/,"",$0); seq=seq $0}
  END{
    printf(">cons_x\n");
    for(i=0;i<R;i++) printf("%s", seq);
    printf("\n");
  }' "$CONS" > "${TMPROOT}/cons.x.fa" && echo "${TMPROOT}/cons.x.fa")
CONS_X="${TMPROOT}/cons.x.fa"

# Helpers
scaf_len(){
  # prints length or empty if missing
  awk -v s="$1" 'BEGIN{FS="\t"} $1==s {print $2; exit}' "$FAI"
}

pass_bits_and_len(){
  # return 0 if pass
  local bits="$1" alen="$2"
  awk -v bs="$bits" -v mbs="$MIN_BITS" -v al="$alen" -v mal="$MIN_ALN" 'BEGIN{ if((bs+0)>=mbs && (al+0)>=mal) exit 0; exit 1 }'
}

best_hit(){
  # stdin: m8C; output: bits alen pid qs qe ss se OR "NA 0 NA 0 0 0 0"
  awk 'BEGIN{best=-1; al=0; pid="NA"; qs=0; qe=0; ss=0; se=0}
       $0~/^#/{next}
       {bs=$12+0; if(bs>best){best=bs; pid=$3; al=$4; qs=$7; qe=$8; ss=$9; se=$10}}
       END{
         if(best<0){print "NA 0 NA 0 0 0 0"}
         else{print best, al, pid, qs, qe, ss, se}
       }'
}

extract_window(){
  # extract_window chrom s1 e1 id outfa
  local chrom="$1" s1="$2" e1="$3" id="$4" outfa="$5"
  # If scaffold missing, samtools faidx returns non-zero (and may print nothing).
  if ! samtools faidx "$GENOME" "${chrom}:${s1}-${e1}" 2>/dev/null \
    | awk -v nm="$id" 'BEGIN{ORS=""} /^>/{print ">"nm"\n"; next} {gsub(/[ \t\r\n]/,"",$0); print $0} END{print "\n"}' \
    > "$outfa"
  then
    return 1
  fi
  # Also guard against empty output (rare but possible)
  if ! grep -q '^>' "$outfa"; then
    return 1
  fi
  return 0
}

# Main loop
idx=0
c_signal=0
c_flip=0
c_edge=0
c_max=0
c_scafmiss=0
c_faidx=0

# Read BED robustly (tabs/spaces). Keep seed name if present.
while IFS= read -r line; do
  [[ -z "$line" ]] && continue
  [[ "$line" =~ ^# ]] && continue

  chrom=$(awk 'BEGIN{FS="[ \t]+"} {print $1}' <<<"$line")
  [[ "$chrom" == "chrom" || "$chrom" == "CHROM" || "$chrom" == "Chr" ]] && continue

  seed_start0=$(awk 'BEGIN{FS="[ \t]+"} {print $2}' <<<"$line")
  seed_end0=$(awk   'BEGIN{FS="[ \t]+"} {print $3}' <<<"$line")
  name=$(awk 'BEGIN{FS="[ \t]+"} {print (NF>=4?$4:"")}' <<<"$line")
  [[ -z "$name" ]] && name="${chrom}:$((seed_start0+1))-${seed_end0}"

  # numeric safety
  seed_start0=$((seed_start0+0))
  seed_end0=$((seed_end0+0))

  idx=$((idx+1))
  if (( idx==1 || (PROG_EVERY>0 && idx%PROG_EVERY==0) )); then
    log "Processing locus ${idx}/${NLOCI} : ${name}"
  fi

  L_scaf_len=$(scaf_len "$chrom" || true)
  if [[ -z "$L_scaf_len" ]]; then
    # Scaffold not in genome index => mark both ends as scaffold_missing
    c_scafmiss=$((c_scafmiss+1))
    origL=$((seed_start0+1))
    origR=$((seed_end0))
    printf "%s\t%d\t%d\t%s\tL\t-\t0\tNA\tNA\tscaffold_missing\t%d\t%d\t0\t0\tNA\t0\tNA\tNA\n" \
      "$chrom" "$seed_start0" "$seed_end0" "$name" "$origL" "$origL" >> "$STOP_TSV"
    printf "%s\t%d\t%d\t%s\tR\t+\t0\tNA\tNA\tscaffold_missing\t%d\t%d\t0\t0\tNA\t0\tNA\tNA\n" \
      "$chrom" "$seed_start0" "$seed_end0" "$name" "$origR" "$origR" >> "$STOP_TSV"

    # keep locus unmodified
    printf "%s\t%d\t%d\t%s\n" "$chrom" "$seed_start0" "$seed_end0" "$name" >> "$OUT_BED"
    printf "%s\t%d\t%d\t%s\t%d\t%d\tscaffold_missing\tscaffold_missing\t0\t0\t0\t0\t0\t0\n" \
      "$chrom" "$seed_start0" "$seed_end0" "$name" "$seed_start0" "$seed_end0" >> "$OUT_TSV"
    continue
  fi

  # Process two ends
  cleanL=$((seed_start0+1))
  cleanR=$((seed_end0))
  reasonL="NA"
  reasonR="NA"
  walkL=0
  walkR=0
  expL=0
  expR=0

  for end in L R; do
    if [[ "$end" == "L" ]]; then
      DIR="-"
      anchor=$((seed_start0+1))   # 1-based
      max_i=$MAX_STEPS
    else
      DIR="+"
      anchor=$((seed_end0))       # 1-based
      max_i=$MAX_STEPS
    fi

    saw_pass=0
    last_pass_s1=0
    last_pass_e1=0
    last_bits="NA"
    last_alen=0
    last_pid="NA"
    last_sstrand="NA"

    stop_reason="NA"
    stop_step=0
    stop_win_s1="NA"
    stop_win_e1="NA"
    clean_pos1=$anchor

    # walk
    i=1
    while (( i <= max_i )); do
      if [[ "$DIR" == "-" ]]; then
        win_e1=$(( anchor - ( (i-1) * STEP_BP ) ))
        win_s1=$(( win_e1 - WIN_BP + 1 ))
      else
        win_s1=$(( anchor + ( (i-1) * STEP_BP ) ))
        win_e1=$(( win_s1 + WIN_BP - 1 ))
      fi

      # Edge conditions
      if (( win_s1 < 1 )); then
        if (( TO_SCAFF_ENDS==1 )); then
          win_s1=1
        else
          stop_reason="past_scaffold_start"
          stop_step=$i
          stop_win_s1=$win_s1
          stop_win_e1=$win_e1
          clean_pos1=1
          c_edge=$((c_edge+1))
          break
        fi
      fi
      if (( win_e1 > L_scaf_len )); then
        if (( TO_SCAFF_ENDS==1 )); then
          win_e1=$L_scaf_len
        else
          stop_reason="past_scaffold_end"
          stop_step=$i
          stop_win_s1=$win_s1
          stop_win_e1=$win_e1
          clean_pos1=$L_scaf_len
          c_edge=$((c_edge+1))
          break
        fi
      fi

      # If clamped window becomes invalid, stop
      if (( win_s1 > win_e1 )); then
        stop_reason="window_invalid"
        stop_step=$i
        stop_win_s1=$win_s1
        stop_win_e1=$win_e1
        clean_pos1=$anchor
        break
      fi

      # Extract + ssearch best hit
      wtmp=$(mktemp -d -p "$TMPROOT" "w.XXXXXXXX") || exit 1
      winfa="${wtmp}/win.fa"

      wid="${name}|${end}|step=${i}|${chrom}:${win_s1}-${win_e1}"
      if ! extract_window "$chrom" "$win_s1" "$win_e1" "$wid" "$winfa"; then
        rm -rf "$wtmp"
        stop_reason="faidx_failed"
        stop_step=$i
        stop_win_s1=$win_s1
        stop_win_e1=$win_e1
        clean_pos1=$anchor
        c_faidx=$((c_faidx+1))
        break
      fi

      read -r bits alen pid qs qe ss se < <(
        ssearch36 -m 8C -E 1000000 "$winfa" "$CONS_X" 2>/dev/null | best_hit
      )

      rm -rf "$wtmp"

      # Determine subject strand for flip detection
      sstrand=$(awk -v a="$ss" -v b="$se" 'BEGIN{print (b+0>=a+0?"+":"-")}')

      if [[ "$bits" != "NA" ]] && pass_bits_and_len "$bits" "$alen"; then
        if (( saw_pass==1 )) && [[ "$last_sstrand" != "NA" ]] && [[ "$sstrand" != "$last_sstrand" ]]; then
          stop_reason="flip"
          stop_step=$i
          stop_win_s1=$win_s1
          stop_win_e1=$win_e1
          # keep boundary at the last passing window
          if [[ "$end" == "L" ]]; then clean_pos1=$last_pass_s1; else clean_pos1=$last_pass_e1; fi
          c_flip=$((c_flip+1))
          break
        fi

        saw_pass=1
        last_pass_s1=$win_s1
        last_pass_e1=$win_e1
        last_bits=$bits
        last_alen=$alen
        last_pid=$pid
        last_sstrand=$sstrand

        i=$((i+1))
        continue
      else
        # signal lost
        stop_reason="signal_lost"
        stop_step=$i
        stop_win_s1=$win_s1
        stop_win_e1=$win_e1
        if (( saw_pass==0 )); then
          clean_pos1=$anchor
        else
          if [[ "$end" == "L" ]]; then
            clean_pos1=$last_pass_s1
          else
            clean_pos1=$last_pass_e1
          fi
        fi
        c_signal=$((c_signal+1))
        break
      fi
    done

    if [[ "$stop_reason" == "NA" ]]; then
      stop_reason="max_steps"
      stop_step=$max_i
      stop_win_s1="NA"
      stop_win_e1="NA"
      if (( saw_pass==0 )); then
        clean_pos1=$anchor
      else
        if [[ "$end" == "L" ]]; then clean_pos1=$last_pass_s1; else clean_pos1=$last_pass_e1; fi
      fi
      c_max=$((c_max+1))
    fi

    orig_endpos1=$anchor
    walked_bp=$(( (stop_step>0 ? (stop_step-1) : 0) * STEP_BP ))
    # expanded_bp relative to original end:
    if [[ "$end" == "L" ]]; then
      expanded_bp=$(( orig_endpos1 - clean_pos1 ))
      (( expanded_bp < 0 )) && expanded_bp=0
      cleanL=$clean_pos1
      reasonL=$stop_reason
      walkL=$((orig_endpos1 - clean_pos1)); (( walkL < 0 )) && walkL=0
      expL=$expanded_bp
    else
      expanded_bp=$(( clean_pos1 - orig_endpos1 ))
      (( expanded_bp < 0 )) && expanded_bp=0
      cleanR=$clean_pos1
      reasonR=$stop_reason
      walkR=$((clean_pos1 - orig_endpos1)); (( walkR < 0 )) && walkR=0
      expR=$expanded_bp
    fi

    printf "%s\t%d\t%d\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%s\n" \
      "$chrom" "$seed_start0" "$seed_end0" "$name" \
      "$end" "$DIR" "$stop_step" "$stop_win_s1" "$stop_win_e1" "$stop_reason" \
      "$clean_pos1" "$orig_endpos1" \
      "$expanded_bp" "$expanded_bp" \
      "$last_bits" "$last_alen" "$last_pid" "$last_sstrand" \
      >> "$STOP_TSV"
  done

  # Write per-locus expanded BED + TSV
  # Ensure L <= R
  if (( cleanL > cleanR )); then
    # fallback to seed if weird
    cleanL=$((seed_start0+1))
    cleanR=$((seed_end0))
    reasonL="bad_clean"
    reasonR="bad_clean"
    expL=0; expR=0; walkL=0; walkR=0
  fi

  new_start0=$((cleanL-1))
  new_end0=$((cleanR))

  printf "%s\t%d\t%d\t%s\n" "$chrom" "$new_start0" "$new_end0" "$name" >> "$OUT_BED"

  expTot=$((expL+expR))
  walkTot=$((walkL+walkR))

  printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n" \
    "$chrom" "$new_start0" "$new_end0" "$name" \
    "$seed_start0" "$seed_end0" \
    "$reasonL" "$reasonR" \
    "$expL" "$expR" "$expTot" \
    "$walkL" "$walkR" "$walkTot" \
    >> "$OUT_TSV"

done < "$BED"

log "Done. loci=$NLOCI ends=$((NLOCI*2)) stops: signal_lost=$c_signal flip=$c_flip scaffold_edge=$c_edge max_steps=$c_max scaffold_missing=$c_scafmiss faidx_failed=$c_faidx"
log "Outputs:"
log "  $STOP_TSV"
log "  $OUT_BED"
log "  $OUT_TSV"
