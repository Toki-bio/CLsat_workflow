#!/usr/bin/env bash
# TRF scan of one very long scaffold in overlapping chunks, so TRF does not hang on a giant
# sequence: extract the scaffold, split it, run TRF per chunk in parallel, translate the chunk
# coordinates back to scaffold coordinates and de-duplicate calls seen twice in an overlap.
# Overlap only needs to exceed a few monomers (arrays split across a chunk boundary reappear as
# two adjacent seeds and get reconnected by the same boundary-expansion pipeline used for
# everything else) -- no need for enormous overlap. The published run needed this for the
# 85.83 Mb unp scaffold.
#
# Usage: run_chunked.sh <genome.fa> <scaffold_name> <out_dir>
# Environment (defaults = published run):
#   CHUNK_BP     chunk length                                  [2000000]
#   OVERLAP_BP   overlap between consecutive chunks            [50000]
#   TRF_PARAMS   TRF parameters                                [2 5 7 80 10 50 2000]
#   JOBS         parallel TRF jobs                             [16]
# Tools: samtools, trf, parallel, python3 on PATH.
# Outputs in out_dir: scaf.fa, chunk_coords.tsv, chunks/chunk_N.{fa,dat,err},
#                     scaf.trf.merged.dat, scaf.trf.dedup.dat (scaffold coordinates), trf_chunked.done
set -euo pipefail
[[ $# -eq 3 ]] || { echo "Usage: $0 <genome.fa> <scaffold_name> <out_dir>" >&2; exit 1; }
GENOME="$1"; SCFNAME="$2"; V="$3"
CHUNK="${CHUNK_BP:-2000000}"
OVERLAP="${OVERLAP_BP:-50000}"
TRF_PARAMS="${TRF_PARAMS:-2 5 7 80 10 50 2000}"
JOBS="${JOBS:-16}"
for t in samtools trf parallel python3; do command -v $t >/dev/null || { echo "$t not in PATH" >&2; exit 1; }; done
[[ "$TRF_PARAMS" =~ ^[0-9]+(\ [0-9]+){6}$ ]] || { echo "TRF_PARAMS must be 7 integers" >&2; exit 1; }

mkdir -p "$V"
cd "$V"
[[ -s "$GENOME.fai" ]] || samtools faidx "$GENOME"
samtools faidx "$GENOME" "$SCFNAME" > scaf.fa
samtools faidx scaf.fa
LEN=$(awk -F'\t' '{print $2}' scaf.fa.fai)
echo "scaffold $SCFNAME length: $LEN"

mkdir -p chunks
rm -f chunks/*.fa chunks/*.dat chunks/*.err
: > chunk_coords.tsv
start=1
n=0
while [ "$start" -le "$LEN" ]; do
  end=$(( start + CHUNK - 1 ))
  [ "$end" -gt "$LEN" ] && end=$LEN
  n=$((n+1))
  cf="chunks/chunk_${n}.fa"
  samtools faidx scaf.fa "${SCFNAME}:${start}-${end}" | \
    awk -v nm="chunk_${n}" 'BEGIN{ORS=""} /^>/{print ">"nm"\n"; next} {gsub(/[ \t\r\n]/,"",$0); print $0} END{print "\n"}' \
    > "$cf"
  printf "%d\t%d\t%d\n" "$n" "$start" "$end" >> chunk_coords.tsv
  nextstart=$(( end - OVERLAP + 1 ))
  [ "$nextstart" -le "$start" ] && nextstart=$(( start + 1 ))
  start=$nextstart
  [ "$end" -eq "$LEN" ] && break
done
echo "chunks created: $n"

echo "=== running TRF on each chunk (parallel) ==="
run_one() {
  cf="$1"
  base=$(basename "$cf" .fa)
  # shellcheck disable=SC2086
  trf "$cf" $TRF_PARAMS -h -ngs > "chunks/${base}.dat" 2> "chunks/${base}.err" || true
}
export -f run_one
export TRF_PARAMS
ls chunks/chunk_*.fa | parallel --jobs "$JOBS" run_one {}

echo "=== translate chunk-local coordinates back to scaffold coordinates ==="
python3 - << 'PYEOF'
import glob, re
coords = {}
for line in open('chunk_coords.tsv'):
    n, s, e = line.split()
    coords[int(n)] = int(s)

out = open('scaf.trf.merged.dat', 'w')
for f in sorted(glob.glob('chunks/chunk_*.dat'), key=lambda x: int(re.search(r'chunk_(\d+)', x).group(1))):
    n = int(re.search(r'chunk_(\d+)', f).group(1))
    offset = coords[n] - 1   # chunk pos 1 == scaffold pos coords[n]
    for line in open(f):
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
            s, e = int(parts[0]), int(parts[1])
            parts[0] = str(s + offset)
            parts[1] = str(e + offset)
            out.write(' '.join(parts) + '\n')
out.close()
print('wrote scaf.trf.merged.dat')
PYEOF

echo "=== dedupe overlapping calls (same period/copies within a few bp = same repeat seen twice) ==="
sort -k1,1n -k4,4n scaf.trf.merged.dat | \
  awk '{key=$4"_"$3; if(key!=prevkey || $1-prevend>100){print; prevkey=key} prevend=$2}' \
  > scaf.trf.dedup.dat
wc -l scaf.trf.merged.dat scaf.trf.dedup.dat
echo TRF_CHUNKED_DONE > trf_chunked.done
