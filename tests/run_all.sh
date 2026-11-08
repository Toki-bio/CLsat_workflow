#!/usr/bin/env bash
# tests/run_all.sh - the test suite of the bash driver, plain bash (no bats). Exit 0 only when every test passes.
#   1. pytest for the stage-script tests that stay Python (PORT_SPEC 5e): tests/test_monomerize.py, tests/test_flank_rules.py
#      (+ tests/test_flank_copy.py when present)
#   2. every tests/sh/test_*.sh (bash tests of the driver: schema, check, config, anchor, subx)
#   (the tests/port/ equivalence harness against the Python driver was removed at the switch-over, WP14; it lives at tag
#   python-driver-final)
# Usage: tests/run_all.sh [--only NAME_SUBSTRING]
# Env:   PYTHON (python3 for the stage scripts), PYTEST (default "$PYTHON -m pytest"; any python with pytest), TMPDIR (all work dirs go below it), KEEP_LOGS=dir (logs kept there)
# Needs a POSIX host with gawk and jq (servers; not Windows).
set -uo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
PY=${PYTHON:-python3}
export LC_ALL=C PYTHON=$PY
ONLY=""
while [ $# -gt 0 ]; do
    case $1 in
        --only) ONLY=$2; shift ;;
        -h|--help) sed -n 2,11p "$0"; exit 0 ;;
        *) echo "run_all.sh: unknown argument $1" >&2; exit 2 ;;
    esac
    shift
done
T=$(mktemp -d "${TMPDIR:-/tmp}/spt.XXXX") || exit 2
L=${KEEP_LOGS:-$T/logs}; mkdir -p "$L"
NP=0; NF=0; FAILED=()

run() {   # run NAME cmd...   (cwd = repo root; output to $L/NAME.log)
    local name=$1 rc=0 t0 log; shift
    if [ -n "$ONLY" ] && [[ $name != *"$ONLY"* ]]; then return 0; fi
    log=$L/${name//\//_}.log
    t0=$(date +%s)
    ( cd "$ROOT" && "$@" ) > "$log" 2>&1 || rc=$?
    if [ $rc = 0 ]; then NP=$((NP + 1)); printf 'PASS %-28s %4ss  %s\n' "$name" $(( $(date +%s) - t0 )) "$(grep -v '^[[:space:]]*$' "$log" | tail -n 1 | cut -c1-90)"
    else NF=$((NF + 1)); FAILED+=("$name"); printf 'FAIL %-28s %4ss  rc=%s  log %s\n' "$name" $(( $(date +%s) - t0 )) "$rc" "$log"
         tail -n 15 "$log" | sed 's/^/     | /'; fi
}

echo "satpipe tests  $(date +%FT%T)  $(hostname)  repo $(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo '?')"
echo "gawk: $(gawk --version 2>/dev/null | head -n 1)   jq: $(jq --version 2>/dev/null)   python: $("$PY" --version 2>&1)"

# 1. pytest (stage scripts)
PYT=(tests/test_monomerize.py tests/test_flank_rules.py)
[ -f "$ROOT/tests/test_flank_copy.py" ] && PYT+=(tests/test_flank_copy.py)
read -ra PYTEST_CMD <<< "${PYTEST:-$PY -m pytest}"
run pytest_stage_scripts "${PYTEST_CMD[@]}" -q -p no:cacheprovider --noconftest "${PYT[@]}"

# 2. bash tests of the driver
for f in "$ROOT"/tests/sh/test_*.sh; do
    n=${f##*/}; run "sh/${n%.sh}" bash "tests/sh/$n"
done


echo "passed $NP, failed $NF${FAILED[*]:+  (${FAILED[*]})}   logs: $L"
[ -z "${KEEP_LOGS:-}" ] && [ $NF = 0 ] && rm -rf -- "$T"
[ $NF = 0 ]
