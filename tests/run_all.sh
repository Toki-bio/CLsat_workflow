#!/usr/bin/env bash
# tests/run_all.sh - the test suite of the bash driver, plain bash (no bats). Exit 0 only when every test passes.
#   1. pytest for the stage-script tests that stay Python (PORT_SPEC 5e): tests/test_monomerize.py, tests/test_flank_rules.py
#      (+ tests/test_flank_copy.py when present)
#   2. every tests/sh/test_*.sh (bash tests of the driver: schema, check, config, anchor, subx)
#   3. the cheap tests/port/*.sh equivalence checks that need no bioinformatics tools (only bash, gawk, jq, python3 with
#      PyYAML as the Python-driver oracle): pyjson/hash fixtures, YAML differential, golden config, check/CLI diff, anchor
#      differential, fingerprint/DONE, driver on stub stages, stages 16-20 on fake scripts. These compare with the Python
#      driver and are dropped at switch-over (PORT_SPEC 8, WP14).
# Usage: tests/run_all.sh [--no-port] [--only NAME_SUBSTRING]
# Env:   PYTHON (python3 with PyYAML), PYTEST (default "$PYTHON -m pytest"; any python with pytest), TMPDIR (all work dirs go below it), KEEP_LOGS=dir (logs kept there)
# Needs a POSIX host with gawk and jq (servers; not Windows).
set -uo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
PY=${PYTHON:-python3}
export LC_ALL=C PYTHON=$PY
PORT=1; ONLY=""
while [ $# -gt 0 ]; do
    case $1 in
        --no-port) PORT=0 ;;
        --only) ONLY=$2; shift ;;
        -h|--help) sed -n 2,13p "$0"; exit 0 ;;
        *) echo "run_all.sh: unknown argument $1" >&2; exit 2 ;;
    esac
    shift
done
# short work path: `satpipe check` refuses an out_dir longer than 70 characters (ssearch36 path limit), and the stub runs
# of test_driver.sh / test_stages_wp8.sh put theirs at <work>/driver/proj/runs/t
T=$(mktemp -d "${TMPDIR:-/tmp}/spt.XXXX") || exit 2
[ ${#T} -le 40 ] || echo "WARNING: work dir $T is long; test_driver/test_stages_wp8 fail the 70-character out_dir check above ~45" >&2
export CONFIGS_PART=${CONFIGS_PART-}   # check_diff.sh reads it under set -u (unset -> 'unbound variable')
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
run pytest_stage_scripts "${PYTEST_CMD[@]}" -q -p no:cacheprovider --noconftest "${PYT[@]}"   # conftest.py imports the Python driver (PyYAML)

# 2. bash tests of the driver
for f in "$ROOT"/tests/sh/test_*.sh; do
    n=${f##*/}; run "sh/${n%.sh}" bash "tests/sh/$n"
done

# 3. cheap equivalence checks (no stage tools)
if [ $PORT = 1 ]; then
    run port/test_pyjson       bash tests/port/test_pyjson.sh
    run port/test_yaml_diff    bash tests/port/test_yaml_diff.sh
    run port/golden_config     bash tests/port/golden_config.sh
    run port/check_diff        bash tests/port/check_diff.sh
    run port/test_anchor_diff  bash tests/port/test_anchor_diff.sh 1000 9 "$PY"
    run port/test_done         bash tests/port/test_done.sh "$PY"
    run port/test_driver       bash tests/port/test_driver.sh "$T/driver"
    run port/test_stages_wp8   bash tests/port/test_stages_wp8.sh "$T/wp8"
fi

echo "passed $NP, failed $NF${FAILED[*]:+  (${FAILED[*]})}   logs: $L"
[ -z "${KEEP_LOGS:-}" ] && [ $NF = 0 ] && rm -rf -- "$T"
[ $NF = 0 ]
