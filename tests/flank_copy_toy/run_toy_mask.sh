#!/bin/bash
# End-to-end toy for --repeat-mask / --cut-at-locus with the real tools. Usage: run_toy_mask.sh <work_dir> <trf_binary>
set -euo pipefail
W="${1:?work dir}"; TRF="${2:?trf binary}"
HERE="$(cd "$(dirname "$0")" && pwd)"; S="$HERE/../../satpipe/stages/04_flank_orthology/rebuild_ortho_segmented.py"
mkdir -p "$W"; cd "$W" || exit 1
python3 "$HERE/make_toy_mask.py" toy
for g in spa spb; do
  samtools faidx toy/genomes/$g.fa
  (cd toy/genomes && "$TRF" $g.fa 2 5 7 80 10 50 2000 -h -d > /dev/null 2>&1 || true)
  [ -s toy/genomes/$g.fa.2.5.7.80.10.50.2000.dat ] || { echo "FATAL: no TRF .dat for $g"; exit 1; }
done
TA="--trf-dat spa=toy/genomes/spa.fa.2.5.7.80.10.50.2000.dat --trf-dat spb=toy/genomes/spb.fa.2.5.7.80.10.50.2000.dat"
ON="$TA --single-copy --clear-margin 0.9 --reciprocal --threads 2"
C="--loci toy/locus_master.tsv --genomes toy/genomes.tsv --families toy/families.tsv"
python3 "$S" $C --out rules $ON > rules.log
python3 "$S" $C --out mask600 $ON --repeat-mask > mask600.log
python3 "$S" $C --out long $ON --repeat-mask --cut-at-locus --flank 3000 > long.log
python3 "$S" $C --out long2 $ON --repeat-mask --cut-at-locus --flank 3000 > long2.log
python3 "$HERE/check_toy_mask.py" toy/expected.json rules mask600 long
for f in ortho_pairs.tsv flank_class.tsv flank_copies.tsv flanks.json; do
  cmp -s long/$f long2/$f || { echo "FAIL determinism (long): $f differs"; exit 1; }
done
echo "PASS determinism (long, two runs byte-identical)"
echo TOYMASK_DONE
