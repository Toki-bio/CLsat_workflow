#!/bin/bash
# End-to-end toy test of the stage-19 single-copy anchor rules with the real tools (TRF, ssearch36, minimap2, samtools).
# Usage: run_toy.sh <work_dir> <trf_binary>
# Runs rebuild_ortho_segmented.py twice on the toy (rules on / rules off), twice more for determinism, then checks.
set -euo pipefail
W="${1:?work dir}"; TRF="${2:?trf binary}"
HERE="$(cd "$(dirname "$0")" && pwd)"; S="$HERE/../../satpipe/stages/04_flank_orthology/rebuild_ortho_segmented.py"
mkdir -p "$W"; cd "$W" || exit 1
python3 "$HERE/make_toy.py" toy
for g in spa spb; do
  samtools faidx toy/genomes/$g.fa
  (cd toy/genomes && "$TRF" $g.fa 2 5 7 80 10 50 2000 -h -d > /dev/null 2>&1 || true)
  [ -s toy/genomes/$g.fa.2.5.7.80.10.50.2000.dat ] || { echo "FATAL: no TRF .dat for $g"; exit 1; }
done
TA="--trf-dat spa=toy/genomes/spa.fa.2.5.7.80.10.50.2000.dat --trf-dat spb=toy/genomes/spb.fa.2.5.7.80.10.50.2000.dat"
ON="$TA --single-copy --clear-margin 0.9 --reciprocal --threads 2"
COMMON="--loci toy/locus_master.tsv --genomes toy/genomes.tsv --families toy/families.tsv"
python3 "$S" $COMMON --out on1 $ON > on1.log
python3 "$S" $COMMON --out on2 $ON > on2.log
python3 "$S" $COMMON --out off1 > off1.log
python3 "$S" $COMMON --out off2 > off2.log
python3 "$HERE/check_toy.py" toy/expected.json on1 off1
for f in ortho_pairs.tsv flank_class.tsv flank_copies.tsv flank_orthologs.json flanks.json ortho_oneway.tsv; do
  cmp -s on1/$f on2/$f || { echo "FAIL determinism (rules on): $f differs"; exit 1; }
done
for f in ortho_pairs.tsv flank_class.tsv flank_orthologs.json flanks.json; do
  cmp -s off1/$f off2/$f || { echo "FAIL determinism (rules off): $f differs"; exit 1; }
done
echo "PASS determinism (two runs byte-identical, rules on and off)"
echo TOY_DONE
