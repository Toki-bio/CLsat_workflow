#!/bin/bash
# Toy for satpipe.ortho_stack with the real tools (samtools, ssearch36, mafft on PATH). Usage: run_toy.sh <work_dir>
set -euo pipefail
W="${1:?work dir}"
HERE="$(cd "$(dirname "$0")" && pwd)"; REPO="$HERE/../.."
mkdir -p "$W"; cd "$W" || exit 1
python3 "$HERE/make_toy.py" toy
for c in a b c; do samtools faidx toy/genomes/$c.fa; done
for k in 1 2; do
  PYTHONPATH="$REPO" python3 - "$k" <<'PY'
import json, shutil, sys
from satpipe.ortho_stack import build_stacks
k = sys.argv[1]
loci = json.load(open('toy/flanks.json'))
pairs = [tuple(l.split('\t')[:2]) for l in open('toy/ortho_pairs.tsv').read().split('\n')[1:] if l]
cons = {'F': open('toy/cons.fa').read().split('\n')[1]}
build_stacks(loci, pairs, ['a', 'b', 'c'], cons, {'F': '#60a5fa'}, {c: 'toy/genomes/%s.fa' % c for c in 'abc'},
             'out%s' % k, {t: shutil.which(t) for t in ('samtools', 'ssearch36', 'mafft')})
PY
done
python3 "$HERE/check_toy.py" toy out1
for f in stacks.json stack/G001.fa stack/G001.bed stack/G001.mask.json stack/G002.fa; do
  cmp -s out1/$f out2/$f || { echo "FAIL determinism: $f differs"; exit 1; }
done
echo "PASS determinism (two runs byte-identical)"
echo TOYSTACK_DONE
