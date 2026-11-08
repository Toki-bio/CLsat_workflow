#!/usr/bin/env bash
# WP8: subx_dir (lib/satpipe/subx.sh = run.py:108-131), the three tests of tests/test_subx_dir.py, plus a differential
# check of the multi-FASTA case against Python's subx_dir on an awkward FASTA (CRLF, a lone CR, descriptions, wrapped and
# padded lines, a junk line before the first header, a repeated name).
# The config node subfamily.consensuses comes from a TYPED-flat fixture ({'subfamily': {'consensuses': cons}}).
# Usage: tests/sh/test_subx.sh      env: PYTHON (default python3; the differential part is skipped when it cannot import satpipe)
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd -- "$HERE/../.." && pwd)
source "$ROOT/lib/satpipe/common.sh"
source "$ROOT/lib/satpipe/flat.sh"
source "$ROOT/lib/satpipe/config.sh"
source "$ROOT/lib/satpipe/subx.sh"
require_gawk
PY=${PYTHON:-python3}
TMP=$(mktemp -d "${TMPDIR:-/tmp}/satpipe_subx.XXXXXX")
trap 'rm -rf -- "$TMP"' EXIT
N=0
check() { local label=$1; shift; N=$((N + 1)); if "$@"; then printf 'ok   %3d %s\n' "$N" "$label"; else printf 'FAIL %3d %s\n' "$N" "$label"; exit 1; fi; }
stage_die() { printf 'stage_die: %s\n' "$1" >&2; exit 1; }

# ctx_for(tmp, cons): RUN and a config whose subfamily.consensuses is CONS (or null when the argument is absent)
ctx_for() {
    RUN=$1/run
    CFG_FLAT=$1/cfg.flat
    if [ $# -ge 2 ]; then
        { flat_emit subfamily map 1 file; flat_emit subfamily/consensuses path "$2" file; } > "$CFG_FLAT"
    else
        { flat_emit subfamily map 1 file; flat_emit subfamily/consensuses null None default; } > "$CFG_FLAT"
    fi
}
same() { [ "$(cat "$1"; printf x)" = "$(printf "$2"; printf x)" ]; }

# test_fasta_becomes_doubled_files
T=$TMP/t1; mkdir -p "$T"
printf '>CLsat1.1 x\nAAAG\nCT\n>CLsat3.2\nGGTTA\n' > "$T/subs.fasta"
ctx_for "$T" "$T/subs.fasta"; subx_dir
check "fasta: returns run/tables/subx"            [ "$SUBX" = "$T/run/tables/subx" ]
check "fasta: CLsat1.1.fa doubled"                same "$T/run/tables/subx/CLsat1.1.fa" '>CLsat1.1\nAAAGCTAAAGCT\n'
check "fasta: CLsat3.2.fa doubled"                same "$T/run/tables/subx/CLsat3.2.fa" '>CLsat3.2\nGGTTAGGTTA\n'
check "fasta: exactly two files"                  [ "$(ls "$T/run/tables/subx" | wc -l)" -eq 2 ]

# test_directory_used_as_given
T=$TMP/t2; mkdir -p "$T/subx"
ctx_for "$T" "$T/subx"; subx_dir
check "directory: used as given"                  [ "$SUBX" = "$T/subx" ]
check "directory: nothing written"                [ ! -e "$T/run" ]

# test_missing_or_absent
T=$TMP/t3; mkdir -p "$T"
ctx_for "$T" "$T/nope.fa"; subx_dir
check "missing path: None"                        [ -z "$SUBX" ]
ctx_for "$T"; subx_dir
check "absent (null): None"                       [ -z "$SUBX" ]
check "missing/absent: nothing written"           [ ! -e "$T/run" ]

# differential against Python's subx_dir
T=$TMP/t4; mkdir -p "$T"
printf 'junk line\r\n>A.1 desc\twith tab\r\nAC GT\r\n  ac \r\n\r\n>B\rGG\rTT\n>A.1\nCCCC\n>C.3  \n\n>D\nA\x1fC\x0b\n' > "$T/odd.fa"
ctx_for "$T" "$T/odd.fa"; subx_dir
if PYTHONPATH=$ROOT "$PY" -c 'import satpipe.run' 2>/dev/null; then
    PYTHONPATH=$ROOT "$PY" - "$T/odd.fa" "$T/pyrun" <<'EOF'
import sys, types
from satpipe.run import subx_dir
print(subx_dir(types.SimpleNamespace(run=sys.argv[2], cfg=types.SimpleNamespace(data={'subfamily': {'consensuses': sys.argv[1]}}))))
EOF
    check "differential: same file set and bytes"  diff -r "$T/pyrun/tables/subx" "$T/run/tables/subx"
else
    echo "skip differential: $PY cannot import satpipe.run"
fi
echo "passed $N"
