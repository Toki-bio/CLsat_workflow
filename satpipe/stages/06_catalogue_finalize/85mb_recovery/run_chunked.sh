#!/usr/bin/env bash
# TRF scan of one very long scaffold in overlapping chunks, so TRF does not hang on a giant
# sequence: extract the scaffold, split it, run TRF per chunk in parallel, translate the chunk
# coordinates back to scaffold coordinates and de-duplicate calls seen twice in an overlap.
# Overlap only needs to exceed a few monomers (arrays split across a chunk boundary reappear as
# two adjacent seeds and get reconnected by the same boundary-expansion pipeline used for
# everything else) -- no need for enormous overlap. The published run needed this for the
# 85.83 Mb unp scaffold.
#
# A chunk on which TRF does not finish within TRF_TIMEOUT is cut into NSPLIT overlapping pieces
# (overlap OVERLAP_BP, capped so the pieces shrink) and each piece is scanned again, recursively,
# until the pieces finish or are shorter than MIN_PIECE_BP; pieces that still time out are written
# to unscanned.tsv and left out, never imputed. This generalises the bisection cascade of the
# published run, whose one pathological chunk ended at a ~57 kb telomeric (TTAGGG)n tail that TRF
# never terminates on. Worst case per level = one TRF_TIMEOUT; 2 Mb -> 20 kb takes four levels.
#
# Usage: run_chunked.sh <genome.fa> <scaffold_name> <out_dir>
# Environment (defaults = published run):
#   CHUNK_BP      chunk length                                   [2000000]
#   OVERLAP_BP    overlap between consecutive chunks / halves    [50000]
#   TRF_PARAMS    TRF parameters                                 [2 5 7 80 10 50 2000]
#   TRF_TIMEOUT   seconds before a piece is bisected             [1800]
#   MIN_PIECE_BP  smallest piece still bisected                  [20000]
#   JOBS          parallel TRF jobs                              [16]
# Tools: samtools, trf, parallel, python3, timeout (coreutils) on PATH.
# Outputs in out_dir: scaf.fa, piece_coords.tsv (name start end, 1-based), chunks/<name>.{fa,dat,err},
#                     unscanned.tsv, scaf.trf.merged.dat, scaf.trf.dedup.dat (scaffold coordinates), trf_chunked.done
set -euo pipefail
[[ $# -eq 3 ]] || { echo "Usage: $0 <genome.fa> <scaffold_name> <out_dir>" >&2; exit 1; }
GENOME="$1"; SCFNAME="$2"; V="$3"
CHUNK="${CHUNK_BP:-2000000}"
OVERLAP="${OVERLAP_BP:-50000}"
TRF_PARAMS="${TRF_PARAMS:-2 5 7 80 10 50 2000}"
TRF_TIMEOUT="${TRF_TIMEOUT:-1800}"
MIN_PIECE="${MIN_PIECE_BP:-20000}"
NSPLIT="${NSPLIT:-4}"           # a timed-out piece is cut into this many overlapping pieces
JOBS="${JOBS:-16}"
for t in samtools trf parallel python3 timeout; do command -v $t >/dev/null || { echo "$t not in PATH" >&2; exit 1; }; done
[[ "$TRF_PARAMS" =~ ^[0-9]+(\ [0-9]+){6}$ ]] || { echo "TRF_PARAMS must be 7 integers" >&2; exit 1; }
# overlap >= chunk would advance the chunker by one base per piece (GLM audit #135)
[ "$OVERLAP" -lt "$CHUNK" ] || { echo "OVERLAP_BP ($OVERLAP) must be smaller than CHUNK_BP ($CHUNK)" >&2; exit 1; }

mkdir -p "$V"
cd "$V"
[[ -s "$GENOME.fai" ]] || samtools faidx "$GENOME"
samtools faidx "$GENOME" "$SCFNAME" > scaf.fa
samtools faidx scaf.fa
LEN=$(awk -F'\t' '{print $2}' scaf.fa.fai)
echo "scaffold $SCFNAME length: $LEN"

mkdir -p chunks
rm -f chunks/*.fa chunks/*.dat chunks/*.err chunks/*.timeout
: > piece_coords.tsv
: > unscanned.tsv

# extract_piece <name> <start1> <end1>  -> chunks/<name>.fa, records piece_coords.tsv
extract_piece(){
  samtools faidx scaf.fa "${SCFNAME}:$2-$3" | \
    awk -v nm="$1" 'BEGIN{ORS=""} /^>/{print ">"nm"\n"; next} {gsub(/[ \t\r\n]/,"",$0); print $0} END{print "\n"}' \
    > "chunks/$1.fa"
  printf "%s\t%d\t%d\n" "$1" "$2" "$3" >> piece_coords.tsv
}

start=1; n=0
while [ "$start" -le "$LEN" ]; do
  end=$(( start + CHUNK - 1 ))
  [ "$end" -gt "$LEN" ] && end=$LEN
  n=$((n+1))
  extract_piece "chunk_${n}" "$start" "$end"
  nextstart=$(( end - OVERLAP + 1 ))
  [ "$nextstart" -le "$start" ] && nextstart=$(( start + 1 ))
  start=$nextstart
  [ "$end" -eq "$LEN" ] && break
done
echo "chunks created: $n"

# run_one <chunks/name.fa>: TRF with a timeout; a timed-out piece leaves <name>.timeout
run_one() {
  cf="$1"
  base=$(basename "$cf" .fa)
  rm -f "chunks/${base}.dat" "chunks/${base}.timeout" "chunks/${base}.failed"
  timeout "$TRF_TIMEOUT" trf "$cf" $TRF_PARAMS -h -ngs > "chunks/${base}.dat" 2> "chunks/${base}.err"
  rc=$?
  if [ "$rc" -eq 124 ]; then rm -f "chunks/${base}.dat"; touch "chunks/${base}.timeout"
  elif [ "$rc" -ge 125 ]; then
    # killed / crashed (OOM 137, segfault 139, timeout's own errors 125-127): the .dat is partial, never accept it.
    # (Lower non-zero codes are not used: some TRF builds exit non-zero after a normal run.)
    rm -f "chunks/${base}.dat"; echo "$rc" > "chunks/${base}.failed"
  fi
  return 0
}
export -f run_one
export TRF_PARAMS TRF_TIMEOUT

# scan a list of pieces in parallel; bisect the ones that time out; repeat until none is left
scan_list(){ # scan_list <file with one chunks/name.fa per line>
  parallel --jobs "$JOBS" run_one {} < "$1"
  if compgen -G "chunks/*.failed" > /dev/null; then
    echo "FATAL: TRF crashed on: $(for f in chunks/*.failed; do printf '%s (exit %s) ' "$(basename "$f" .failed)" "$(cat "$f")"; done)" >&2
    exit 1
  fi
}
level=0
ls chunks/chunk_*.fa > pieces.level0.txt
echo "=== running TRF on each chunk (parallel, timeout ${TRF_TIMEOUT}s) ==="
scan_list pieces.level0.txt
while :; do
  next="pieces.level$((level+1)).txt"; : > "$next"
  for t in chunks/*.timeout; do
    [ -e "$t" ] || continue
    name=$(basename "$t" .timeout)
    read -r _ s e < <(awk -v n="$name" -F'\t' '$1==n{print; exit}' piece_coords.tsv)
    plen=$(( e - s + 1 ))
    rm -f "$t"
    if [ "$plen" -le "$MIN_PIECE" ]; then
      printf "%s\t%d\t%d\t%d\tTRF did not finish within %ss on a %d bp piece\n" "$name" "$s" "$e" "$plen" "$TRF_TIMEOUT" "$plen" >> unscanned.tsv
      echo "  UNSCANNED: $name ${SCFNAME}:$s-$e ($plen bp)"
      continue
    fi
    # split into NSPLIT pieces; the overlap is capped so that every piece is strictly shorter
    # than its parent (an uncapped 50 kb overlap on a 60 kb piece would never shrink)
    ov=$OVERLAP; [ "$ov" -gt $(( plen / (4 * NSPLIT) )) ] && ov=$(( plen / (4 * NSPLIT) ))
    step=$(( (plen + NSPLIT - 1) / NSPLIT ))
    for k in $(seq 0 $((NSPLIT - 1))); do
      ps_=$(( s + k * step - ov )); [ "$ps_" -lt "$s" ] && ps_=$s
      pe_=$(( s + (k + 1) * step - 1 + ov )); [ "$pe_" -gt "$e" ] && pe_=$e
      [ "$ps_" -le "$pe_" ] || continue
      sub="${name}_$((k + 1))"
      extract_piece "$sub" "$ps_" "$pe_"; echo "chunks/$sub.fa" >> "$next"
      echo "  split $name (${plen} bp) -> $sub ${ps_}-${pe_}"
    done
  done
  [ -s "$next" ] || break
  level=$((level+1))
  echo "=== level $level: $(wc -l < "$next") pieces ==="
  scan_list "$next"
  [ "$level" -lt 40 ] || { echo "bisection did not converge" >&2; exit 1; }
done

echo "=== translate piece-local coordinates back to scaffold coordinates ==="
python3 - << 'PYEOF'
import glob, os
coords = {}
for line in open('piece_coords.tsv'):
    n, s, e = line.split()
    coords[n] = int(s)
out = open('scaf.trf.merged.dat', 'w')
n_rows = 0
for f in sorted(glob.glob('chunks/*.dat'), key=lambda x: coords.get(os.path.basename(x)[:-4], 0)):
    n = os.path.basename(f)[:-4]
    offset = coords[n] - 1   # piece pos 1 == scaffold pos coords[n]
    for line in open(f):
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
            s, e = int(parts[0]), int(parts[1])
            parts[0] = str(s + offset)
            parts[1] = str(e + offset)
            out.write(' '.join(parts) + '\n'); n_rows += 1
out.close()
print('wrote scaf.trf.merged.dat: %d rows from %d pieces' % (n_rows, len(glob.glob('chunks/*.dat'))))
PYEOF

echo "=== dedupe calls seen twice in an overlap (identical start, end and period in scaffold coordinates) ==="
# Only exact duplicates are removed. The earlier rule (same period and copies, start within 100 bp of the previous
# row's end) also dropped distinct neighbouring arrays and let duplicates through when another call sorted between
# them (GLM audit 2026-10-24). A call cut at a chunk boundary stays as two seeds, which expansion joins.
sort -k1,1n -k2,2n -k3,3n scaf.trf.merged.dat | awk '!seen[$1" "$2" "$3]++' > scaf.trf.dedup.dat
wc -l scaf.trf.merged.dat scaf.trf.dedup.dat
if [ -s unscanned.tsv ]; then echo "unscanned pieces (see unscanned.tsv): $(wc -l < unscanned.tsv)"; fi
echo TRF_CHUNKED_DONE > trf_chunked.done
