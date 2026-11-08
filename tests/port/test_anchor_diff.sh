#!/usr/bin/env bash
# WP9 acceptance: differential test anchor.awk vs satpipe/anchor.py on N generated cases (default 1000, seed 9).
# Usage: tests/port/test_anchor_diff.sh [N] [SEED] [PYTHON]     (work dir: $TMPDIR; kept when a difference is found)
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd -- "$HERE/../.." && pwd)
source "$ROOT/lib/satpipe/common.sh"
require_gawk
N=${1:-1000}; SEED=${2:-9}; PY=${3:-python3}
W=$(mktemp -d "${TMPDIR:-/tmp}/satpipe_wp9.XXXXXX")
PYTHONIOENCODING=utf-8 "$PY" "$HERE/anchor_oracle.py" gen "$W" "$N" "$SEED"
t0=$(date +%s.%N)
PYTHONIOENCODING=utf-8 "$PY" "$HERE/anchor_oracle.py" run "$W/cases.tsv" | tr -d '\r' > "$W/py.out"
t1=$(date +%s.%N)
"$GAWK" -f "$PYJSON" -f "$SATPIPE_LIB/anchor.awk" -f "$HERE/anchor_diff.awk" "$W/cases.tsv" > "$W/gawk.out"
t2=$(date +%s.%N)
n=$(wc -l < "$W/py.out"); ng=$(wc -l < "$W/gawk.out")
nd=$(diff "$W/py.out" "$W/gawk.out" | grep -c '^<' || true)
"$GAWK" -F '\t' '{ k[$1]++; if ($0 ~ /\tERR\t/) e[$1]++ } END { for (x in k) printf "  %-3s %4d cases, %4d end in an exception\n", x, k[x], e[x] }' "$W/py.out" | LC_ALL=C sort
"$GAWK" '/\tOUT\t/ { o++ } /snapped/ { s++ } /after rotation the sequence starts/ { m++ } /ambiguous/ { a++ } /not within/ { w++ } /appears twice/ { d++ }
    /cannot place/ { c++ } /not found in the anchor/ { f++ }
    END { printf "  paths: P placed %d, snapped %d, motif-mismatch %d, frame ambiguous %d, motif not in anchor %d, cannot place %d, snap-not-within %d, duplicate name %d\n", o, s, m, a, f, c, w, d }' "$W/py.out"
printf 'cases %d (python lines %d, gawk lines %d), differing lines %d; time python %.1fs gawk %.1fs\n' \
    "$N" "$n" "$ng" "$nd" "$("$GAWK" -v a="$t0" -v b="$t1" "BEGIN{print b-a}")" "$("$GAWK" -v a="$t1" -v b="$t2" "BEGIN{print b-a}")"
if [ "$nd" = 0 ] && [ "$n" = "$ng" ] && [ "$n" = "$N" ]; then echo "PASS"; rm -rf -- "$W"; else echo "FAIL (kept $W)"; diff "$W/py.out" "$W/gawk.out" | head -20; exit 1; fi
