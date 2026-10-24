#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# clisat_expand_then_scan_merged.sh  (cli_v6: fixed --min-aln-frac base; referred to as
# cli_v6.sh / cli.sh elsewhere in this repository — 06_catalogue_finalize/cli_v6.sh links here)
# ============================================================================
#
# PURPOSE
#   Given TRF seed loci (BED4) for one CLsat family, expands each locus
#   outward until satellite signal is lost, merges overlapping expanded
#   intervals per species+scaffold, then scans every merged interval at
#   monomer-step resolution to produce a bits-score profile.
#
# TYPICAL USAGE
#   ./clisat_expand_then_scan_merged.sh CLsat1.bed genomes.tsv CLsat1.cons.fa 146 out_CLsat1 --jobs 16
#
# POSITIONAL ARGUMENTS
#   1  BED            BED3+ file; col4 (name) must start "species|..."
#   2  GENOMES_TSV    two-column TSV: species <TAB> /abs/path/genome.fa
#   3  CONS           consensus FASTA (single sequence, one monomer)
#   4  MONO           monomer length in bp (integer, e.g. 146)
#   5  OUT_PREFIX     prefix for all output files (directory created if needed)
#
# OPTIONS
#
#   Mode
#   --scan-only         Skip PASS1 (expand) and PASS2 (merge). Treat BED as
#                       already-merged intervals and run PASS3 (scan) only.
#                       Use for 2nd-scan runs on MASTER.merged.bed to avoid
#                       re-expansion beyond 1st-scan boundaries. Without this
#                       flag the script re-expands each interval outward by up
#                       to 3*MONO before scanning, which shifts boundaries
#                       beyond what the 1st scan established.
#
#   Expansion (PASS 1) -- ignored with --scan-only
#   --step-monomers N   step size between expansion windows, in monomers [3]
#   --win-monomers  N   window width for each expansion probe, in monomers [3]
#   --query-repeat  N   ssearch36 query = (N+1) tandem copies of consensus [2]
#                       Default N=2 -> 3-copy query of length 3*MONO.
#                       The extra copy avoids boundary edge effects. The doc
#                       value N is NOT the copy count; copy count = N+1.
#   --min-aln-frac  F   fraction of the WINDOW (WIN_BP = 3*MONO) that must align [0.67]
#   --min-bits      N   minimum ssearch36 bits score to accept a window [50]
#   --max-steps     N   hard cap on dense expansion steps per direction [500]
#                       At step=3*MONO=146 this gives ~217 kb max per side.
#                       Set to 0 for unlimited dense walk (no probe needed).
#   --to-scaffold-ends  Clamp windows to scaffold boundaries and keep walking
#                       (default: stop at scaffold boundary and record reason).
#
#   Three-phase expansion strategy (PASS 1)
#   ----------------------------------------
#   Phase 1: dense batched walk (up to --max-steps steps).
#     One ssearch36 call covers all windows in a direction. Stops when signal
#     lost, strand flips, scaffold boundary reached, or max-steps hit.
#
#   Phase 2: coarse binary probe (fires only when Phase 1 hits max-steps with
#     signal still present, i.e. array continues beyond the cap).
#     Advances from Phase 1 boundary in large jumps (--probe-monomers*MONO bp).
#     Signal present -> take the jump. Signal absent -> halve the step.
#     Stops when step < --probe-min-step*MONO. O(log N) calls.
#     Scaffold boundary during probe: clean_pos1 = last confirmed position
#     (probe_pos), NOT the scaffold edge -- avoids mis-reporting boundary.
#     Probe window is probe_step wide (>> WIN_BP); the MIN_BITS threshold is the
#     same number, but the TEST is weaker: a probe block passes if ANY hit in it
#     reaches MIN_BITS, with no MIN_ALN requirement, so a gap of up to PROBE_BP
#     (100*MONO = 14.6 kb) inside a long array is bridged, whereas the dense walk
#     (Phase 1) stops at the first failing 3-monomer window. Arrays that reach
#     Phase 2 (> MAX_STEPS*STEP_BP ~ 219 kb) therefore have boundaries that are
#     reproducible only to within a few kb. This is the likely cause of the
#     split-vs-merged differences on megabase arrays noted in PIPELINE.md.
#
#   Phase 3: dense pinpoint (fires after Phase 2).
#     Starts ONE STEP beyond probe_pos (probe_pos is already confirmed passing;
#     re-checking it would waste one ssearch36 call per direction).
#     Capped at 2*(PROBE_BP/STEP_BP)+4 steps for bounded cost.
#     Gives monomer-step boundary resolution regardless of array size.
#
#   --probe-monomers N  coarse probe initial step, in monomers [100]
#   --probe-min-step N  coarse probe minimum step / Phase 3 search width,
#                       in monomers [50]
#
#   Scan (PASS 3)
#   --scan-win-bp   N   scan window width in bp [2*MONO]
#   --scan-step-bp  N   scan step in bp [MONO]
#   --scan-min-bits N   minimum bits to emit a scan hit [0]
#
#   Parallelism
#   --jobs          N   GNU parallel workers for PASS1 and PASS3 [24]
#
# OUTPUTS
#   OUT.expanded_per_locus.bed   BED4: one expanded interval per seed locus
#                                (empty in --scan-only mode)
#   OUT.expanded_per_locus.tsv   Metadata per locus (header-only in --scan-only)
#   OUT.expanded_merged.bed      BED4: merged intervals per species+scaffold
#                                (copy of input BED in --scan-only mode)
#   OUT.bits.tsv                 Scan profile: scaffold pos bits name species
#                                pos is 1-based; use pos-1 as BED start
#
# STOP REASONS (reason_L / reason_R columns)
#   signal_lost          bits score fell below --min-bits
#   flip                 alignment strand reversed
#   past_scaffold_start  walk/probe reached scaffold start; clean_pos1 is
#                        last confirmed position (not forced to 1)
#   past_scaffold_end    walk/probe reached scaffold end; clean_pos1 is
#                        last confirmed position (not forced to L_scaf_len)
#   window_invalid       computed window coordinates were invalid
#   max_steps            dense walk hit cap; probe also exhausted
#   coarse_probe         boundary found by coarse probe + dense pinpoint
#
# WHY GIANT ARRAYS HIT max_steps
#   TRF reports short (~2-4 kb) windows even within megabase arrays, so seed
#   loci represent interior TRF fragments of large CLsat blocks. The default
#   --max-steps 500 caps expansion at ~217 kb for MONO=146. CLsat1 arrays in
#   val/unm/unp can exceed 1 Mb. The coarse+pinpoint probe automatically
#   extends to the true boundary at O(log N) extra cost.
#   Example: unm JAWDKO010000075.1 array ~1 Mb; dense walk stops at 248 kb;
#   coarse probe (~17 calls) reaches 15 kb; pinpoint (12 calls) finds exact
#   boundary at 10,439 bp -- 238 kb extra coverage, ~29 total extra calls.
#
# REQUIREMENTS
#   samtools >= 1.12 (for --region-file), ssearch36, bedtools, GNU parallel
# ============================================================================

die(){ echo "ERROR: $*" >&2; exit 1; }
log(){ echo "[expand_then_scan] $*" >&2; }

usage(){ awk 'NR>=4 && !/^#/{exit} NR>=4{sub(/^# ?/,""); print}' "$0"; }
[[ "${1:-}" == "-h" || "${1:-}" == "--help" ]] && { usage; exit 0; }
[[ $# -lt 5 ]] && die "Usage: $0 loci.bed genomes.tsv consensus.fa MONO out_prefix [options]  (--help for details)"

BED=$1; GENOMES_TSV=$2; CONS=$3; MONO=$4; OUT=$5
shift 5

STEP_MON=3; WIN_MON=3; QREP=2; MIN_ALN_FRAC=0.67; MIN_BITS=50
MAX_STEPS=500; TO_SCAFF_ENDS=0
PROBE_MON=100; PROBE_MIN_STEP=50
SCAN_WIN_BP=$((2*MONO)); SCAN_STEP_BP=$((1*MONO)); SCAN_MIN_BITS=0
JOBS=24
SCAN_ONLY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --step-monomers)    STEP_MON=$2;       shift 2;;
    --win-monomers)     WIN_MON=$2;        shift 2;;
    --query-repeat)     QREP=$2;           shift 2;;
    --min-aln-frac)     MIN_ALN_FRAC=$2;   shift 2;;
    --min-bits)         MIN_BITS=$2;       shift 2;;
    --max-steps)        MAX_STEPS=$2;      shift 2;;
    --to-scaffold-ends) TO_SCAFF_ENDS=1;   shift;;
    --probe-monomers)   PROBE_MON=$2;      shift 2;;
    --probe-min-step)   PROBE_MIN_STEP=$2; shift 2;;
    --scan-win-bp)      SCAN_WIN_BP=$2;    shift 2;;
    --scan-step-bp)     SCAN_STEP_BP=$2;   shift 2;;
    --scan-min-bits)    SCAN_MIN_BITS=$2;  shift 2;;
    --jobs)             JOBS=$2;           shift 2;;
    --scan-only)        SCAN_ONLY=1;       shift;;
    -h|--help) usage; exit 0;;
    *) die "Unknown option: $1";;
  esac
done

[[ -s "$BED" ]]         || die "bed missing/empty: $BED"
[[ -s "$GENOMES_TSV" ]] || die "genomes.tsv missing/empty: $GENOMES_TSV"
[[ -s "$CONS" ]]        || die "consensus missing/empty: $CONS"
[[ "$MONO" =~ ^[0-9]+$ ]] || die "MONO must be integer"

command -v samtools  >/dev/null 2>&1 || die "samtools not found"
command -v ssearch36 >/dev/null 2>&1 || die "ssearch36 not found"
command -v bedtools  >/dev/null 2>&1 || die "bedtools not found"
command -v parallel  >/dev/null 2>&1 || die "GNU parallel not found"

BED="$(readlink -f "$BED")"
GENOMES_TSV="$(readlink -f "$GENOMES_TSV")"
CONS="$(readlink -f "$CONS")"

OUTDIR=$(dirname "$OUT")
mkdir -p "$OUTDIR"

STEP_BP=$(( STEP_MON * MONO ))
WIN_BP=$((  WIN_MON  * MONO ))
# (the query is (QREP+1) copies of the consensus; its length is no longer used — see below)
# FIX: MIN_ALN is the minimum ALIGNED LENGTH a WINDOW must yield. It must be a fraction of
# WIN_BP (the window), not of QLEN (the query). With QLEN and MIN_ALN_FRAC=1.0 it demanded a
# 438-column alignment from a 438bp window -- unsatisfiable, so the expansion never ran.
MIN_ALN=$(awk -v w="$WIN_BP" -v f="$MIN_ALN_FRAC" 'BEGIN{m=int(w*f); if(m<1)m=1; print m}')
MIN_BITS_INT=${MIN_BITS%%.*}
SCAN_MIN_BITS_INT=${SCAN_MIN_BITS%%.*}
PROBE_BP=$(( PROBE_MON * MONO ))
PROBE_MIN_BP=$(( PROBE_MIN_STEP * MONO ))

EXP_BED="${OUT}.expanded_per_locus.bed"
EXP_TSV="${OUT}.expanded_per_locus.tsv"
MRG_BED="${OUT}.expanded_merged.bed"
BITS_TSV="${OUT}.bits.tsv"

printf "species\tchrom\tstart0\tend0\tname\tseed_start0\tseed_end0\treason_L\treason_R\texpand_L_bp\texpand_R_bp\texpand_total_bp\n" > "$EXP_TSV"
printf "scaffold\tpos\tbits\tname\tspecies\n" > "$BITS_TSV"

# genome map
declare -A GENOME_OF
while IFS=$'\t' read -r sp gp; do
  [[ -z "${sp:-}" || "$sp" =~ ^# ]] && continue
  [[ -n "${gp:-}" ]] || die "Bad genomes.tsv line: $sp"
  gp="$(readlink -f "$gp")" || die "Cannot resolve genome path: $sp"
  [[ -s "$gp" && -r "$gp" ]] || die "Genome missing/unreadable: $sp: $gp"
  GENOME_OF["$sp"]="$gp"
done < "$GENOMES_TSV"

for sp in "${!GENOME_OF[@]}"; do
  g="${GENOME_OF[$sp]}"
  [[ -s "${g}.fai" ]] || samtools faidx "$g"
done

TMPROOT=$(mktemp -d -p "${OUTDIR}" "tmp.expand_then_scan.XXXXXXXX")
trap 'rm -rf "$TMPROOT"' EXIT

CONS_X="${TMPROOT}/cons.x.fa"
awk -v R=$((QREP+1)) '
  BEGIN{seq=""}
  /^>/{next}
  {gsub(/[ \t\r\n]/,"",$0); seq=seq $0}
  END{
    printf(">cons_x\n");
    for(i=0;i<R;i++) printf("%s", seq);
    printf("\n");
  }' "$CONS" > "$CONS_X"

GENOME_MAP="${TMPROOT}/genome_map.tsv"
for sp in "${!GENOME_OF[@]}"; do
  printf "%s\t%s\n" "$sp" "${GENOME_OF[$sp]}"
done > "$GENOME_MAP"

# ss_retry <out> <ssearch36 args...>: ssearch36 under parallel load on a shared file system fails now and then
# ("cannot open library" with a garbled name). Those failures used to be hidden (2>/dev/null) and scored as no
# signal; now a call is retried 3 times and only then fails the job (seen 2026-10-24 on Monsoon).
ss_retry(){
  local out="$1"; shift; local k
  for k in 1 2 3; do
    if ssearch36 "$@" > "$out" 2> "$out.err"; then rm -f "$out.err"; return 0; fi
    sleep "$k"
  done
  echo "FATAL: ssearch36 $* failed 3 times: $(tail -2 "$out.err" | tr '\n' ' ')" >&2
  return 1
}
export -f ss_retry

# ── scan-only: skip PASS1+2 entirely ────────────────────────────────────────
if [[ "$SCAN_ONLY" == "1" ]]; then
  log "SCAN_ONLY mode: skipping PASS1+2, using input BED directly as merged regions"
  # every row must be tab-separated chrom start end sp|... with a genome of the genome table: a row that is not
  # would otherwise be dropped by scan_region without a trace (GLM audit 2026-10-24)
  awk -F'\t' 'NR==FNR{g[$1]=1; next} /^[ \t]*$|^#/{next}
    { i=index($4,"|"); sp=(i>1?substr($4,1,i-1):"")
      if(NF<4 || $2!~/^[0-9]+$/ || $3!~/^[0-9]+$/ || !(sp in g)){ print "FATAL: scan-only BED row " FNR " is not chrom<TAB>start<TAB>end<TAB>sp|... with a known genome: " $0 > "/dev/stderr"; bad=1; exit } }
    END{ exit bad }' "$GENOME_MAP" "$BED" || exit 1
  cp "$BED" "$MRG_BED"
  : > "$EXP_BED"
  # jump straight to PASS3
else

# ============================================================
# expand_locus: three-phase expansion worker
# ============================================================
expand_locus(){
  local TMPROOT="$1"  line="$2"  GENOME_MAP="$3"  CONS_X="$4"  CONS="$5"
  local STEP_BP="$6"  WIN_BP="$7"  MIN_BITS_INT="$8"  MIN_ALN="$9"
  local MAX_STEPS="${10}"  TO_SCAFF_ENDS="${11}"
  local PROBE_BP="${12}"  PROBE_MIN_BP="${13}"
  local SCAN_WIN_BP="${14}"  SCAN_STEP_BP="${15}"  SCAN_MIN_BITS_INT="${16}"

  [[ -z "$line" || "$line" =~ ^# ]] && return 0

  local chrom seed_start0 seed_end0 name sp genome fai
  chrom=$(      awk 'BEGIN{FS="[ \t]+"}{print $1}' <<<"$line")
  seed_start0=$(awk 'BEGIN{FS="[ \t]+"}{print $2}' <<<"$line")
  seed_end0=$(  awk 'BEGIN{FS="[ \t]+"}{print $3}' <<<"$line")
  name=$(       awk 'BEGIN{FS="[ \t]+"}{print (NF>=4?$4:"")}' <<<"$line")

  # Every locus that cannot be processed is reported on a DROP line (collected into
  # OUT.dropped.tsv) instead of vanishing silently.
  if [[ -z "$name" || "$name" != *"|"* ]]; then
    printf "DROP\t%s\t%s\t%s\t%s\tname_lacks_species_prefix\n" "$chrom" "$seed_start0" "$seed_end0" "$name"; return 0
  fi
  sp="${name%%|*}"
  genome=$(awk -v S="$sp" 'BEGIN{FS="\t"} $1==S{print $2; exit}' "$GENOME_MAP")
  if [[ -z "$genome" ]]; then
    printf "DROP\t%s\t%s\t%s\t%s\tspecies_not_in_genomes_tsv\n" "$chrom" "$seed_start0" "$seed_end0" "$name"; return 0
  fi
  fai="${genome}.fai"

  seed_start0=$((seed_start0+0)); seed_end0=$((seed_end0+0))
  local L_scaf_len
  L_scaf_len=$(awk -v s="$chrom" 'BEGIN{FS="\t"} $1==s{print $2; exit}' "$fai")
  if [[ -z "$L_scaf_len" ]]; then
    printf "DROP\t%s\t%s\t%s\t%s\tscaffold_not_in_fai\n" "$chrom" "$seed_start0" "$seed_end0" "$name"; return 0
  fi

  local cleanL=$((seed_start0+1))  cleanR=$((seed_end0))
  local reasonL="NA"  reasonR="NA"  expL=0  expR=0

  local WDIR="${TMPROOT}/locus_$$_${RANDOM}"
  mkdir -p "$WDIR"

  # Single-window query (used by phases 2 and 3)
  _qwin(){
    local s1=$1 e1=$2
    local wfa="${WDIR}/qwin.fa"
    # a failed fetch or search is an error (return 1), not "no signal"
    samtools faidx "$genome" "${chrom}:${s1}-${e1}" > "${WDIR}/qwin.raw" \
      || { echo "FATAL: samtools faidx $genome ${chrom}:${s1}-${e1} failed" >&2; return 1; }
    awk 'BEGIN{ORS=""} /^>/{print ">w\n";next} \
             {gsub(/[ \t\r\n]/,"",$0);print $0} END{print "\n"}' "${WDIR}/qwin.raw" > "$wfa"
    grep -q '^>' "$wfa" || { echo "FATAL: no sequence for ${chrom}:${s1}-${e1} in $genome" >&2; return 1; }
    ss_retry "${WDIR}/qwin.m8c" -n -m 8C -E 1000000 "$wfa" "$CONS_X" || return 1
    awk 'BEGIN{m=-1} /^[^#]/{b=$12+0; if(b>m)m=b} END{print m+0}' "${WDIR}/qwin.m8c"
  }

  for end in L R; do
    local DIR anchor isLeft
    if [[ "$end" == "L" ]]; then DIR="-"; anchor=$((seed_start0+1)); isLeft=1
    else                          DIR="+"; anchor=$((seed_end0));      isLeft=0; fi

    # ── Phase 1: dense batched walk ──────────────────────────────────────
    local -a WIN_S1 WIN_E1
    WIN_S1=(); WIN_E1=()
    local early_stop_reason="" early_stop_pos=0 i=1

    while (( MAX_STEPS == 0 || i <= MAX_STEPS )); do
      local win_s1 win_e1
      if [[ "$DIR" == "-" ]]; then
        win_e1=$(( anchor - (i-1)*STEP_BP )); win_s1=$(( win_e1 - WIN_BP + 1 ))
      else
        win_s1=$(( anchor + (i-1)*STEP_BP )); win_e1=$(( win_s1 + WIN_BP - 1 ))
      fi
      if (( win_s1 < 1 )); then
        if (( TO_SCAFF_ENDS==1 )); then
          win_s1=1
          # degenerate: entire window past scaffold start -> stop
          if (( win_e1 < 1 )); then
            early_stop_reason="past_scaffold_start"; early_stop_pos=1; break
          fi
        else
          early_stop_reason="past_scaffold_start"; early_stop_pos=1; break
        fi
      fi
      if (( win_e1 > L_scaf_len )); then
        if (( TO_SCAFF_ENDS==1 )); then
          win_e1=$L_scaf_len
          # degenerate: entire window past scaffold end -> stop
          if (( win_s1 > L_scaf_len )); then
            early_stop_reason="past_scaffold_end"; early_stop_pos=$L_scaf_len; break
          fi
        else
          early_stop_reason="past_scaffold_end"; early_stop_pos=$L_scaf_len; break
        fi
      fi
      (( win_s1 <= win_e1 )) || { early_stop_reason="window_invalid"; early_stop_pos=$anchor; break; }
      WIN_S1+=("$win_s1"); WIN_E1+=("$win_e1")
      i=$((i+1))
    done

    local n_wins=${#WIN_S1[@]}
    local regions_file="${WDIR}/regions_${end}.txt"
    local batch_fa="${WDIR}/batch_${end}.fa"
    : > "$regions_file"
    for (( j=0; j<n_wins; j++ )); do
      printf "%s:%d-%d\n" "$chrom" "${WIN_S1[$j]}" "${WIN_E1[$j]}" >> "$regions_file"
    done
    if [[ "$n_wins" -gt 0 ]]; then
      # a failed batch fetch is an error: it used to give an empty batch and 'signal_lost' at 0 bp (GLM audit #20)
      samtools faidx "$genome" --region-file "$regions_file" > "${WDIR}/batch_${end}.raw" \
        || { echo "FATAL: samtools faidx $genome --region-file failed for $name" >&2; rm -rf "$WDIR"; return 1; }
      awk '/^>/{idx++; print ">"idx; next}{print}' "${WDIR}/batch_${end}.raw" > "$batch_fa"
    else
      : > "$batch_fa"
    fi
    local search_out="${WDIR}/search_${end}.m8c"
    if [[ -s "$batch_fa" ]]; then
      ss_retry "$search_out" -n -m 8C -E 1000000 "$batch_fa" "$CONS_X" || { rm -rf "$WDIR"; return 1; }
    else
      : > "$search_out"
    fi

    local -a IDX_BITS IDX_ALEN IDX_SS IDX_SE
    IDX_BITS=(); IDX_ALEN=(); IDX_SS=(); IDX_SE=()
    for (( j=0; j<n_wins; j++ )); do IDX_BITS+=(-1); IDX_ALEN+=(0); IDX_SS+=(0); IDX_SE+=(0); done

    while IFS=$'\t' read -r qid _ _ alen _ _ _ _ ss se _ bits; do
      [[ "$qid" =~ ^[0-9]+$ ]] || continue
      local idx=$(( qid-1 ))
      (( idx < 0 || idx >= n_wins )) && continue
      local prev=${IDX_BITS[$idx]}
      if (( ${bits%%.*} > ${prev%%.*} )) || \
         { (( ${bits%%.*} == ${prev%%.*} )) && [[ "${bits##*.}" > "${prev##*.}" ]]; }; then
        IDX_BITS[$idx]=$bits; IDX_ALEN[$idx]=$alen; IDX_SS[$idx]=$ss; IDX_SE[$idx]=$se
      fi
    done < <(grep -v '^#' "$search_out" 2>/dev/null || true)

    local saw_pass=0 last_pass_s1=0 last_pass_e1=0 last_sstrand="NA"
    local stop_reason="NA" clean_pos1=$anchor

    for (( j=0; j<n_wins; j++ )); do
      local b=${IDX_BITS[$j]} al=${IDX_ALEN[$j]} ss_j=${IDX_SS[$j]} se_j=${IDX_SE[$j]}
      local ws1=${WIN_S1[$j]} we1=${WIN_E1[$j]}
      local pass=0
      if [[ "$b" != "-1" ]] && (( ${b%%.*} >= MIN_BITS_INT )) && (( al >= MIN_ALN )); then
        pass=1
      fi
      if (( pass==1 )); then
        local sstrand; (( se_j >= ss_j )) && sstrand="+" || sstrand="-"
        if (( saw_pass==1 )) && [[ "$last_sstrand" != "NA" && "$sstrand" != "$last_sstrand" ]]; then
          stop_reason="flip"; clean_pos1=$(( isLeft==1 ? last_pass_s1 : last_pass_e1 )); break
        fi
        saw_pass=1; last_pass_s1=$ws1; last_pass_e1=$we1; last_sstrand=$sstrand
      else
        stop_reason="signal_lost"
        clean_pos1=$(( saw_pass==0 ? anchor : (isLeft==1 ? last_pass_s1 : last_pass_e1) ))
        break
      fi
    done

    # -- resolve end of phase 1 --
    if [[ "$stop_reason" == "NA" ]]; then
      if [[ -n "$early_stop_reason" ]]; then
        stop_reason="$early_stop_reason"; clean_pos1=$early_stop_pos

      elif (( MAX_STEPS > 0 && n_wins >= MAX_STEPS )); then
        # phase 1 exhausted its budget with signal still on
        local phase1_end=$(( saw_pass==1 ? (isLeft==1 ? last_pass_s1 : last_pass_e1) : anchor ))

        # ── Phase 2: coarse binary probe ─────────────────────────────────
        local probe_pos=$phase1_end probe_step=$PROBE_BP probe_hit_boundary=0

        while (( probe_step >= PROBE_MIN_BP )); do
          local next_pos
          if [[ "$DIR" == "-" ]]; then next_pos=$(( probe_pos - probe_step ))
          else                          next_pos=$(( probe_pos + probe_step )); fi

          if (( next_pos < 1 )); then
            if (( TO_SCAFF_ENDS==1 )); then next_pos=1
            else
              # use probe_pos (last confirmed signal), not scaffold start
              stop_reason="past_scaffold_start"
              clean_pos1=$probe_pos
              probe_hit_boundary=1; break
            fi
          fi
          if (( next_pos > L_scaf_len )); then
            if (( TO_SCAFF_ENDS==1 )); then next_pos=$L_scaf_len
            else
              # use probe_pos (last confirmed signal), not scaffold end
              stop_reason="past_scaffold_end"
              clean_pos1=$probe_pos
              probe_hit_boundary=1; break
            fi
          fi

          # probe window: step-wide
          local pw_s1 pw_e1
          if [[ "$DIR" == "-" ]]; then
            pw_e1=$next_pos; pw_s1=$(( pw_e1 - probe_step + 1 ))
            (( pw_s1 < 1 )) && pw_s1=1
          else
            pw_s1=$next_pos; pw_e1=$(( pw_s1 + probe_step - 1 ))
            (( pw_e1 > L_scaf_len )) && pw_e1=$L_scaf_len
          fi

          local pb; pb=$(_qwin "$pw_s1" "$pw_e1") || { rm -rf "$WDIR"; return 1; }
          if [[ "$pb" != "-1" ]] && (( ${pb%%.*} >= MIN_BITS_INT )); then
            probe_pos=$next_pos   # signal present: advance
          else
            probe_step=$(( probe_step / 2 ))   # signal absent: halve step
          fi
        done

        if (( probe_hit_boundary == 0 )); then
          # probe_pos = last confirmed position, ±PROBE_MIN_BP from true boundary
          # ── Phase 3: dense pinpoint ──────────────────────────────────────
          # probe_pos is already confirmed passing; start ONE STEP beyond it
          # to avoid re-checking it (saves one ssearch36 call per direction)
          local pinpoint_cap=$(( 2 * PROBE_BP / STEP_BP + 4 ))
          local pp_last_pass=$probe_pos

          for (( k=1; k<=pinpoint_cap; k++ )); do
            local pp_s1 pp_e1
            if [[ "$DIR" == "-" ]]; then
              pp_e1=$(( probe_pos - k*STEP_BP )); pp_s1=$(( pp_e1 - WIN_BP + 1 ))
              (( pp_s1 < 1 )) && pp_s1=1
              if (( pp_e1 < 1 )); then
                pp_last_pass=1; stop_reason="past_scaffold_start"; break
              fi
            else
              pp_s1=$(( probe_pos + k*STEP_BP )); pp_e1=$(( pp_s1 + WIN_BP - 1 ))
              (( pp_e1 > L_scaf_len )) && pp_e1=$L_scaf_len
              if (( pp_s1 > L_scaf_len )); then
                pp_last_pass=$L_scaf_len; stop_reason="past_scaffold_end"; break
              fi
            fi

            local pb2; pb2=$(_qwin "$pp_s1" "$pp_e1") || { rm -rf "$WDIR"; return 1; }
            if [[ "$pb2" != "-1" ]] && (( ${pb2%%.*} >= MIN_BITS_INT )); then
              pp_last_pass=$(( isLeft==1 ? pp_s1 : pp_e1 ))
            else
              stop_reason="signal_lost"; break
            fi
          done

          stop_reason="coarse_probe"
          clean_pos1=$pp_last_pass
        fi
        # if probe_hit_boundary, stop_reason and clean_pos1 already set above

      else
        # MAX_STEPS==0 unlimited walk consumed everything
        stop_reason="max_steps"
        clean_pos1=$(( saw_pass==0 ? anchor : (isLeft==1 ? last_pass_s1 : last_pass_e1) ))
      fi
    fi

    if (( isLeft==1 )); then
      expL=$(( anchor - clean_pos1 )); ((expL<0))&&expL=0
      cleanL=$clean_pos1; reasonL=$stop_reason
    else
      expR=$(( clean_pos1 - anchor )); ((expR<0))&&expR=0
      cleanR=$clean_pos1; reasonR=$stop_reason
    fi
  done

  rm -rf "$WDIR"
  if (( cleanL > cleanR )); then
    printf "DROP\t%s\t%s\t%s\t%s\tempty_after_expansion\n" "$chrom" "$seed_start0" "$seed_end0" "$name"; return 0
  fi

  local new_start0=$((cleanL-1)) new_end0=$((cleanR)) expTot=$((expL+expR))
  printf "BED\t%s\t%d\t%d\t%s\n" "$chrom" "$new_start0" "$new_end0" "$name"
  printf "TSV\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\n" \
    "$sp" "$chrom" "$new_start0" "$new_end0" "$name" \
    "$seed_start0" "$seed_end0" "$reasonL" "$reasonR" "$expL" "$expR" "$expTot"
}
export -f expand_locus

export TMPROOT GENOME_MAP CONS_X CONS
export STEP_BP WIN_BP MIN_BITS_INT MIN_ALN MAX_STEPS TO_SCAFF_ENDS
export PROBE_BP PROBE_MIN_BP
export SCAN_WIN_BP SCAN_STEP_BP SCAN_MIN_BITS_INT

log "PASS1 expand (parallel --jobs $JOBS): step_bp=$STEP_BP win_bp=$WIN_BP min_bits=$MIN_BITS_INT min_aln=$MIN_ALN max_steps=$MAX_STEPS probe_bp=$PROBE_BP probe_min_bp=$PROBE_MIN_BP"

EXPAND_OUT="${TMPROOT}/expand_results.tsv"

# Skip blank lines, comments and a literal header row only. (An earlier version also
# dropped every row whose scaffold name begins with "Chr", i.e. all chromosome-named assemblies.)
grep -v -E '^\s*$|^#' "$BED" | awk 'BEGIN{FS="[ \t]+"} !($1=="chrom"||$1=="CHROM"||$1=="Chr"||$1=="chr")' \
  | parallel --jobs "$JOBS" --keep-order --halt never \
      expand_locus "$TMPROOT" {} "$GENOME_MAP" "$CONS_X" "$CONS" \
        "$STEP_BP" "$WIN_BP" "$MIN_BITS_INT" "$MIN_ALN" "$MAX_STEPS" "$TO_SCAFF_ENDS" \
        "$PROBE_BP" "$PROBE_MIN_BP" \
        "$SCAN_WIN_BP" "$SCAN_STEP_BP" "$SCAN_MIN_BITS_INT" \
  > "$EXPAND_OUT"

grep '^BED' "$EXPAND_OUT" | cut -f2- > "$EXP_BED" || true
grep '^TSV' "$EXPAND_OUT" | cut -f2- >> "$EXP_TSV" || true
DROP_TSV="${OUT}.dropped.tsv"
grep '^DROP' "$EXPAND_OUT" | cut -f2- > "$DROP_TSV" || true
n_in=$(grep -v -E '^\s*$|^#' "$BED" | awk 'BEGIN{FS="[ \t]+"} !($1=="chrom"||$1=="CHROM"||$1=="Chr"||$1=="chr")' | wc -l)
n_out=$(wc -l < "$EXP_BED"); n_drop=$(wc -l < "$DROP_TSV")
log "PASS1 reconciliation: $n_in seeds in, $n_out expanded, $n_drop dropped (see $DROP_TSV)"
(( n_in == n_out + n_drop )) || log "WARNING: seed count does not reconcile ($n_in != $n_out + $n_drop)"

# ── PASS 2 ───────────────────────────────────────────────────────────────────
log "PASS2 merge expanded intervals"

awk -F'\t' 'BEGIN{OFS="\t"}{
  n=$4; sp=n; sub(/\|.*/,"",sp); print sp,$1,$2,$3,n
}' "$EXP_BED" | sort -k2,2 -k3,3n > "${TMPROOT}/exp.sp.bed"

: > "$MRG_BED"
while read -r sp; do
  [[ -z "$sp" ]] && continue
  awk -F'\t' -v SP="$sp" 'BEGIN{OFS="\t"} $1==SP{print $2,$3,$4,$5}' "${TMPROOT}/exp.sp.bed" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - -c 4 -o count,collapse \
    | awk -v SP="$sp" 'BEGIN{OFS="\t"}{
        tag=SP "|" $1 ":" ($2+1) "-" $3 "|nSeeds=" $4;
        print $1,$2,$3,tag
      }' >> "$MRG_BED"
done < <(cut -f1 "${TMPROOT}/exp.sp.bed" | sort -u)

fi  # end non-scan-only block (PASS1+2)

# ============================================================================
# scan_region worker (PASS 3)
# Per-window ssearch36 matches original score behavior exactly.
# Batching shifts bits scores ~1-3 due to ssearch36 E-value calibration
# across multiple queries, so we keep per-window calls here.
# ============================================================================
scan_region(){
  local TMPROOT="$1"  line="$2"  GENOME_MAP="$3"  CONS="$4"
  local SCAN_WIN_BP="$5"  SCAN_STEP_BP="$6"  SCAN_MIN_BITS_INT="$7"

  local chrom start0 end0 tag sp genome
  IFS=$'\t' read -r chrom start0 end0 tag <<<"$line"
  sp="${tag%%|*}"
  genome=$(awk -v S="$sp" 'BEGIN{FS="\t"} $1==S{print $2; exit}' "$GENOME_MAP")
  # an unknown genome or a failed tool call is an error (parallel reports the job, pipefail stops the script);
  # it used to give "no windows" silently (GLM audit 2026-10-24)
  [[ -n "$genome" ]] || { echo "FATAL scan_region: genome '$sp' of '$tag' not in $GENOME_MAP" >&2; return 1; }

  local scan_start=$((start0+1))
  local scan_end=$(( end0 - SCAN_WIN_BP + 1 ))
  (( scan_end < scan_start )) && return 0

  local WDIR="${TMPROOT}/scan_$$_${RANDOM}"
  mkdir -p "$WDIR"
  local winfa="${WDIR}/win.fa"

  for ((pos=scan_start; pos<=scan_end; pos+=SCAN_STEP_BP)); do
    local s1=$pos e1=$((pos + SCAN_WIN_BP - 1))
    samtools faidx "$genome" "${chrom}:${s1}-${e1}" > "${WDIR}/raw.fa" \
      || { echo "FATAL scan_region: samtools faidx $genome ${chrom}:${s1}-${e1} failed" >&2; rm -rf "$WDIR"; return 1; }
    awk -v nm="${tag}|scan|${chrom}:${s1}-${e1}" \
          'BEGIN{ORS=""} /^>/{print ">"nm"\n"; next} \
           {gsub(/[ \t\r\n]/,"",$0); print $0} END{print "\n"}' "${WDIR}/raw.fa" > "$winfa"
    grep -q '^>' "$winfa" || { echo "FATAL scan_region: no sequence for ${chrom}:${s1}-${e1} in $genome" >&2; rm -rf "$WDIR"; return 1; }
    local bits
    ss_retry "${WDIR}/ss.out" -n -m 8C "$winfa" "$CONS" || { rm -rf "$WDIR"; return 1; }
    bits=$(awk 'BEGIN{m=0} /^[^#]/{if(($12+0)>m) m=$12+0} END{print m+0}' "${WDIR}/ss.out")
    if (( ${bits%%.*} >= SCAN_MIN_BITS_INT )); then
      printf "%s\t%d\t%s\t%s\t%s\n" "$chrom" "$pos" "$bits" "$tag" "$sp"
    fi
  done
  rm -rf "$WDIR"
}
export -f scan_region

# ── PASS 3 ───────────────────────────────────────────────────────────────────
log "PASS3 scan merged intervals (parallel --jobs $JOBS)"

grep -v -E '^\s*$|^#' "$MRG_BED" \
  | parallel --jobs "$JOBS" --keep-order --halt never \
      scan_region "$TMPROOT" {} "$GENOME_MAP" "$CONS" \
        "$SCAN_WIN_BP" "$SCAN_STEP_BP" "$SCAN_MIN_BITS_INT" \
  >> "$BITS_TSV"

log "DONE"
log "  $EXP_BED"
log "  $EXP_TSV"
log "  $MRG_BED"
log "  $BITS_TSV"