#!/usr/bin/env bash
# 85 Mb TRF rescan, done properly this time: split into overlapping chunks so TRF doesn't
# hang on one giant sequence, run each chunk in parallel, then merge back to scaffold
# coordinates. Overlap only needs to exceed a few CLsat monomers (arrays split across a
# chunk boundary reappear as two adjacent seeds and get reconnected by the same
# boundary-expansion pipeline used for everything else) -- no need for enormous overlap.
set -euo pipefail

B=/staging/conda/envs/bioinfo/bin
V=/staging/tmp/clsat_v6/scan85
SCAF=$V/scaf.fa
CHUNK=2000000     # 2 Mb chunks
OVERLAP=50000     # 50 kb overlap (>> any single CLsat monomer, generous margin)

cd "$V"
LEN=$($B/samtools faidx "$SCAF" 2>/dev/null; awk '{print $2}' "$SCAF.fai")
echo "scaffold length: $LEN"

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
  $B/samtools faidx "$SCAF" "JAWDKN010000003.1:${start}-${end}" | \
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
  /staging/miniconda3/bin/trf "$cf" 2 5 7 80 10 50 2000 -h -ngs > "chunks/${base}.dat" 2> "chunks/${base}.err" || true
}
export -f run_one
ls chunks/chunk_*.fa | $B/parallel --jobs 16 run_one {}

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
echo TRF_CHUNKED_DONE > trf85.done
