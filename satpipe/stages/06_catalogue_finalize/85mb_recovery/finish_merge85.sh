#!/usr/bin/env bash
# Re-merge the per-piece TRF outputs of run_chunked.sh into scaffold coordinates after any
# piece .dat has been replaced by hand (the published run rebuilt chunk_44 from a bisection
# cascade; see legacy/published_run/85mb_recovery/). Pieces are listed in piece_coords.tsv.
#
# Usage: finish_merge85.sh <out_dir>     (the run_chunked.sh output directory)
set -euo pipefail
[[ $# -eq 1 ]] || { echo "Usage: $0 <out_dir>" >&2; exit 1; }
cd "$1"
[[ -s piece_coords.tsv ]] || { echo "piece_coords.tsv missing in $1" >&2; exit 1; }

echo "=== translate piece-local coordinates back to scaffold coordinates ==="
python3 - << 'PYEOF'
import glob, os
coords = {}
for line in open('piece_coords.tsv'):
    n, s, e = line.split()
    coords[n] = int(s)
out = open('scaf.trf.merged.dat', 'w')
for f in sorted(glob.glob('chunks/*.dat'), key=lambda x: coords.get(os.path.basename(x)[:-4], 0)):
    n = os.path.basename(f)[:-4]
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
echo TRF_CHUNKED_DONE > trf_chunked.done
echo FINISHED
