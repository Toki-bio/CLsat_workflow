#!/usr/bin/env bash
set -euo pipefail
cd /staging/tmp/clsat_v6/scan85

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
    offset = coords[n] - 1
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

echo "=== dedupe overlapping calls ==="
sort -k1,1n -k4,4n scaf.trf.merged.dat | \
  awk '{key=$4"_"$3; if(key!=prevkey || $1-prevend>100){print; prevkey=key} prevend=$2}' \
  > scaf.trf.dedup.dat
wc -l scaf.trf.merged.dat scaf.trf.dedup.dat
echo TRF_CHUNKED_DONE > trf85.done
echo FINISHED
