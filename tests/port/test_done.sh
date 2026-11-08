#!/usr/bin/env bash
# WP0: fingerprint / mark_done / is_done (lib/satpipe/done.sh) against the Python driver (satpipe/run.py) as oracle.
# Needs a POSIX host (symlinks), python3 that can import satpipe.run (PyYAML), and ~70 MB of free temp space.
# Usage: tests/port/test_done.sh [PYTHON]       (TMPDIR decides where the toy files go)
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd -- "$HERE/../.." && pwd)
source "$ROOT/lib/satpipe/common.sh"
source "$ROOT/lib/satpipe/flat.sh"
source "$ROOT/lib/satpipe/done.sh"
require_gawk
PY=${1:-python3}
TMP=$(mktemp -d "${TMPDIR:-/tmp}/satpipe_wp0.XXXXXX")
trap 'rm -rf -- "$TMP"' EXIT
fail=0
report() { printf '%-52s %s\n' "$1" "$2"; case $2 in PASS*) ;; *) fail=1 ;; esac; }
pyrun() { PYTHONPATH=$ROOT "$PY" -c "$@"; }

# ---- toy tree (spec 2.3): small files, a 70 MB file, a dir with symlinked file + dotfile + subdir, symlinked .dat
T=$TMP/toy
mkdir -p "$T/d/sub" "$T/empty_dir"
printf 'ACGT\n' > "$T/small.txt"; : > "$T/empty.txt"
head -c 73400320 /dev/zero > "$T/big.bin"                       # 70 MiB > 64 MiB -> size-mtime
printf 'a\n' > "$T/d/a.txt"; printf 'h\n' > "$T/d/.hidden"; printf 'B\n' > "$T/d/B"; printf '_\n' > "$T/d/_x"
printf 's\n' > "$T/d/sub/inner.txt"
ln -s a.txt "$T/d/link_to_a"; ln -s nowhere "$T/d/broken"; ln -s sub "$T/d/link_dir"; ln -s ../big.bin "$T/d/big_link"
ln -s small.txt "$T/sym.dat"; ln -s big.bin "$T/big.dat"; ln -s d "$T/dirlink"
cases=(small.txt empty.txt big.bin d d/sub empty_dir sym.dat big.dat dirlink d/link_to_a)
n=0; ok=0
for c in "${cases[@]}"; do
    n=$((n + 1))
    b=$(fingerprint "$T/$c")
    p=$(pyrun 'import sys; from satpipe.run import fingerprint; print(fingerprint(sys.argv[1]))' "$T/$c")
    if [ "$b" = "$p" ]; then ok=$((ok + 1)); else log "fingerprint $c: bash=$b python=$p"; fi
done
[ "$ok" = "$n" ] && report "fingerprint bash == Python" "PASS $ok/$n" || report "fingerprint bash == Python" "FAIL $ok/$n"
if fingerprint "$T/missing" 2>/dev/null; then report "fingerprint of a missing path fails" FAIL; else report "fingerprint of a missing path fails" PASS; fi

# ---- mark_done: record bytes equal Python's json.dump(rec, indent=1) with the same 'finished'
RUN=$TMP/run; mkdir -p "$RUN"
json_to_flat "$HERE/fixtures/done_14_arbitrate.json" params > "$TMP/params.flat"
printf '%s\n' "$T/small.txt" "$T/d" "$T/big.dat" "$T/small.txt" "$T/sym.dat" > "$TMP/inputs.lst"
printf '%s\n' "$T/empty.txt" "$T/d" "$T/not_made.tsv" "$T/big.bin" > "$TMP/outputs.lst"
mark_done "$RUN/14_arbitrate/DONE" 14_arbitrate "$TMP/inputs.lst" "$TMP/params.flat" "$TMP/outputs.lst"
pyrun '
import json, sys, os
from satpipe.run import fingerprint, Ctx
done, fx, ins, outs = sys.argv[1:5]
bash = json.load(open(done))
params = json.load(open(fx))["params"]
inputs = [l.rstrip("\n") for l in open(ins) if l.strip()]
outputs = [l.rstrip("\n") for l in open(outs) if l.strip()]
rec = {"stage": "14_arbitrate", "finished": bash["finished"], "params_hash": Ctx.params_hash(params), "params": params,
       "inputs": {f: fingerprint(f) for f in inputs},
       "outputs": {f: fingerprint(f) for f in outputs if os.path.isfile(f)}}
sys.stdout.write(json.dumps(rec, indent=1))' "$RUN/14_arbitrate/DONE" "$HERE/fixtures/done_14_arbitrate.json" "$TMP/inputs.lst" "$TMP/outputs.lst" > "$TMP/py_done.json"
cmp -s "$TMP/py_done.json" "$RUN/14_arbitrate/DONE" && report "mark_done record == Python json.dump(indent=1)" PASS || report "mark_done record" "FAIL (diff $TMP/py_done.json)"
grep -q '"params_hash": "a6e3f7b1e6fc11ea"' "$RUN/14_arbitrate/DONE" && report "mark_done params_hash = real run's a6e3f7b1e6fc11ea" PASS || report "mark_done params_hash" FAIL

# ---- is_done: bash vs Python Ctx.is_done on the same files, scenario by scenario
py_is_done() {   # RUN STAGE INPUTS PARAMS_JSON FORCE
    pyrun '
import json, sys
from satpipe.run import Ctx
run, stage, ins, pj, force = sys.argv[1:6]
inputs = [l.rstrip("\n") for l in open(ins) if l.strip()]
c = Ctx(cfg=None, run=run, force=(force == "1"))
print("1" if c.is_done(stage, inputs, json.load(open(pj))) else "0")' "$@"
}
flat_json "$TMP/params.flat" > "$TMP/params.json"
cmp_is_done() {   # label expected
    local b p
    if SATPIPE_FORCE=${FORCE:-0} is_done "$RUN/14_arbitrate/DONE" "$TMP/inputs.lst" "$TMP/params.flat"; then b=1; else b=0; fi
    p=$(py_is_done "$RUN" 14_arbitrate "$TMP/inputs.lst" "$TMP/params.json" "${FORCE:-0}")
    [ "$b" = "$p" ] && [ "$b" = "$2" ] && report "is_done: $1" "PASS (bash $b python $p)" || report "is_done: $1" "FAIL (bash $b python $p want $2)"
}
cmp_is_done "fresh record" 1
FORCE=1 cmp_is_done "--force" 0
cp "$TMP/inputs.lst" "$TMP/inputs.keep"; printf '%s\n' "$T/empty.txt" >> "$TMP/inputs.lst"
cmp_is_done "planned input set grew" 0
cp "$TMP/inputs.keep" "$TMP/inputs.lst"; head -n 2 "$TMP/inputs.keep" > "$TMP/inputs.lst"
cmp_is_done "planned input set shrank" 0
cp "$TMP/inputs.keep" "$TMP/inputs.lst"
cp "$TMP/params.flat" "$TMP/params.keep"; sed 's/^arbitration\/delta\tint\t20/arbitration\/delta\tfloat\t20.0/' "$TMP/params.keep" > "$TMP/params.flat"
flat_json "$TMP/params.flat" > "$TMP/params.json"
cmp_is_done "param 20 -> 20.0 (hash differs)" 0
cp "$TMP/params.keep" "$TMP/params.flat"; flat_json "$TMP/params.flat" > "$TMP/params.json"
cmp_is_done "params restored" 1
printf 'x' >> "$T/d/.hidden"
cmp_is_done "dotfile inside an input dir changed" 0
cp "$RUN/14_arbitrate/DONE" "$TMP/done.keep"
mark_done "$RUN/14_arbitrate/DONE" 14_arbitrate "$TMP/inputs.lst" "$TMP/params.flat" "$TMP/outputs.lst"
cmp_is_done "re-marked" 1
mv "$T/sym.dat" "$T/sym.dat.away"
cmp_is_done "recorded input missing" 0
mv "$T/sym.dat.away" "$T/sym.dat"
printf '{"stage": ' > "$RUN/14_arbitrate/DONE"
cmp_is_done "unparsable DONE" 0
rm -f -- "$RUN/14_arbitrate/DONE"
cmp_is_done "no DONE" 0

exit "$fail"
