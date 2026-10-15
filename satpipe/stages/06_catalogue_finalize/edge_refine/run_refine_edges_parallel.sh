#!/usr/bin/env bash
set -euo pipefail
V=/staging/tmp/clsat_v6
IN=$V/refine_input_540.bed
OUT=$V/refined_540_20260716
NJOBS=16

WORK=$V/refine_chunks
rm -rf $WORK
mkdir -p $WORK

echo "=== splitting $(wc -l < $IN) loci into $NJOBS chunks ==="
split -n l/$NJOBS -d -a 2 $IN $WORK/chunk_

echo "=== running refine_edges.py on each chunk in parallel ==="
run_one() {
  cf="$1"
  base=$(basename "$cf")
  python3 $V/refine_edges.py "$cf" "$WORK/${base}.out" > "$WORK/${base}.log" 2>&1
}
export -f run_one
export V WORK
ls $WORK/chunk_* | /staging/conda/envs/bioinfo/bin/parallel --jobs $NJOBS run_one {}

echo "=== merging chunk outputs, renumbering loci to match original order ==="
python3 - <<'PYEOF'
import glob, re

WORK = '/staging/tmp/clsat_v6/refine_chunks'
OUT = '/staging/tmp/clsat_v6/refined_540_20260716'

# chunk files in original split order (chunk_00, chunk_01, ...)
chunk_bases = sorted(glob.glob(f'{WORK}/chunk_*.out.bed'))
chunk_bases = [c[:-4] for c in chunk_bases]  # strip .bed, keep the ".out" base

bed_out = open(f'{OUT}.bed', 'w')
rep_out = open(f'{OUT}.report.tsv', 'w')
header_written = False
global_i = 0
old_to_new = {}  # (chunk_base, old_Rn) -> new_Rn

for cb in chunk_bases:
    bed_lines = open(f'{cb}.bed').readlines()
    rep_lines = open(f'{cb}.report.tsv').readlines()
    rep_hdr, rep_lines = rep_lines[0], rep_lines[1:]
    if not header_written:
        rep_out.write(rep_hdr)
        header_written = True
    for bl, rl in zip(bed_lines, rep_lines):
        global_i += 1
        p = bl.rstrip('\n').split('\t')
        old_r = p[3]
        p[3] = f'R{global_i}'
        bed_out.write('\t'.join(p) + '\n')
        rp = rl.rstrip('\n').split('\t')
        rp[0] = f'R{global_i}'
        rep_out.write('\t'.join(rp) + '\n')

bed_out.close()
rep_out.close()
print(f'merged {global_i} loci -> {OUT}.bed, {OUT}.report.tsv')
PYEOF

touch $OUT.DONE
echo REFINE_PARALLEL_DONE
