#!/usr/bin/env bash
set -euo pipefail
# Run refine_edges.py on a locus BED in parallel chunks and merge the pieces back, renumbering
# loci R1..Rn in the original order.
#
# Usage: run_refine_edges_parallel.sh <loci.bed> <out_prefix> <genomes.tsv> <families.tsv>
#   loci.bed  chrom start0 end species family (5 columns)
# Environment: NJOBS [16]
# Output: <out_prefix>.bed, <out_prefix>.report.tsv, <out_prefix>.DONE; chunks in <out_prefix>.chunks/
[[ $# -eq 4 ]] || { echo "Usage: $0 <loci.bed> <out_prefix> <genomes.tsv> <families.tsv>" >&2; exit 1; }
IN="$1"; OUT="$2"; GENOMES="$3"; FAMILIES="$4"
NJOBS="${NJOBS:-16}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
command -v parallel >/dev/null || { echo "GNU parallel not in PATH" >&2; exit 1; }

WORK="$OUT.chunks"
rm -rf "$WORK"
mkdir -p "$WORK" "$(dirname "$OUT")"

echo "=== splitting $(wc -l < "$IN") loci into $NJOBS chunks ==="
split -n l/"$NJOBS" -d -a 2 "$IN" "$WORK/chunk_"

echo "=== running refine_edges.py on each chunk in parallel ==="
run_one() {
  cf="$1"
  base=$(basename "$cf")
  python3 "$HERE/refine_edges.py" "$cf" "$WORK/${base}.out" "$GENOMES" "$FAMILIES" > "$WORK/${base}.log" 2>&1
}
export -f run_one
export HERE WORK GENOMES FAMILIES
ls "$WORK"/chunk_?? | parallel --jobs "$NJOBS" run_one {}

echo "=== merging chunk outputs, renumbering loci to match original order ==="
python3 - "$WORK" "$OUT" <<'PYEOF'
import glob, sys

WORK, OUT = sys.argv[1], sys.argv[2]

# chunk files in original split order (chunk_00, chunk_01, ...)
chunk_bases = sorted(glob.glob(f'{WORK}/chunk_*.out.bed'))
chunk_bases = [c[:-4] for c in chunk_bases]  # strip .bed, keep the ".out" base

bed_out = open(f'{OUT}.bed', 'w')
rep_out = open(f'{OUT}.report.tsv', 'w')
header_written = False
global_i = 0

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
        p[3] = f'R{global_i}'
        bed_out.write('\t'.join(p) + '\n')
        rp = rl.rstrip('\n').split('\t')
        rp[0] = f'R{global_i}'
        rep_out.write('\t'.join(rp) + '\n')

bed_out.close()
rep_out.close()
print(f'merged {global_i} loci -> {OUT}.bed, {OUT}.report.tsv')
PYEOF

touch "$OUT.DONE"
echo REFINE_PARALLEL_DONE
