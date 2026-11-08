#!/usr/bin/env bash
# WP11 molly: tests/molly/run_molly.sh unchanged, once with the Python driver and once with the bash driver as `satpipe`
# (a wrapper first on PATH), in the SAME work directory (moved aside after each), then the comparisons of e2e_toy.sh:
# rc/stdout/trace, every run-dir file (data cmp, DONE masked, manifest jq -S, logs own category), viewer/ (diff_viewer.sh),
# and cross reruns (each driver on a copy of the other's run dir runs nothing).
# Usage: tests/port/e2e_molly.sh WORKDIR MOLLY_TRF_DIR [JOBS]    env: PYTHON (python3 with PyYAML)
# Phases (env PHASES, default "py sh cmp"): run one side only with PHASES=py etc.
set -uo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
W=${1:?usage: e2e_molly.sh WORKDIR MOLLY_TRF_DIR [JOBS]}; M=${2:?molly_trf_dir}; JOBS=${3:-16}
PY=${PYTHON:-python3}; GAWK=${SATPIPE_GAWK:-gawk}
export LC_ALL=C
unset SATPIPE_TRACE SATPIPE_DRY SATPIPE_STAGE_FILES
mkdir -p "$W/out"; W=$(cd "$W" && pwd)
# shellcheck source=e2e_lib.sh
. "$ROOT/tests/port/e2e_lib.sh"
MW=$W/m            # run_molly.sh work dir (same path for both drivers)
mkdir -p "$W/bin_py" "$W/bin_sh"
printf '#!/bin/bash\nPYTHONPATH=%q exec %q -m satpipe.cli "$@"\n' "$ROOT" "$(command -v "$PY")" > "$W/bin_py/satpipe"
printf '#!/bin/bash\nexec %q "$@"\n' "$ROOT/bin/satpipe" > "$W/bin_sh/satpipe"
chmod +x "$W/bin_py/satpipe" "$W/bin_sh/satpipe"

for s in ${PHASES:-py sh cmp}; do
    case $s in
    py|sh)
        echo "### molly $s driver start $(date +%FT%T)"
        rm -rf -- "$MW" "$W/m.$s"
        : > "$W/out/$s.trace"; rc=0
        ( cd "$W" && PATH=$W/bin_$s:$PATH SATPIPE_TRACE="$W/out/$s.trace" bash "$ROOT/tests/molly/run_molly.sh" "$ROOT" "$M" "$MW" "$JOBS" ) \
            > "$W/out/$s.stdout" 2> "$W/out/$s.stderr" || rc=$?
        echo "$rc" > "$W/out/$s.rc"
        echo "### molly $s driver end $(date +%FT%T) rc=$rc"
        mv "$MW" "$W/m.$s" ;;
    cmp)
        echo "--- py vs sh"
        # the CMD records carry PATH, which holds this harness's own wrapper dir (bin_py / bin_sh): one name for both
        for x in py sh; do
            mkdir -p "$W/out/n"
            for k in rc stdout stderr; do cp "$W/out/$x.$k" "$W/out/n/$x.$k"; done
            sed "s#$W/bin_$x#$W/bin_X#g" "$W/out/$x.trace" > "$W/out/n/$x.trace"
        done
        compare_drive molly "$W/out/n/py" "$W/out/n/sh"
        compare_runs molly "$W/m.py/run" "$W/m.sh/run" "$W/out/cmp"
        for f in ref/cluster_family.tsv ref/families_note.txt molly.yaml; do
            cmp -s "$W/m.py/$f" "$W/m.sh/$f" && ok "molly $f identical" || fail "molly $f differs"
        done
        if bash "$ROOT/tests/port/diff_viewer.sh" "$W/m.py/run/viewer" "$W/m.sh/run/viewer" > "$W/out/viewer.diff" 2>&1; then ok "molly viewer dir: $(tail -n 1 "$W/out/viewer.diff")"
        else fail "molly viewer dir" "$W/out/viewer.diff"; fi
        echo "  manifest summary (sh): $(jq -c '{centroids_n, seed_counts, master_loci, segmented_loci, catalogue_loci, ortholog_pairs, candidates}' "$W/m.sh/run/manifest.json")"
        echo "--- cross reruns"
        for pair in "sh py" "py sh"; do
            set -- $pair
            rm -rf -- "$MW"; mkdir -p "$MW"
            for x in ref molly.yaml; do cp -a "$W/m.$2/$x" "$MW/"; done
            cross_rerun molly "$1" "$W/m.$2/run" "$MW/run" "$W/out/x_$1_on_$2" run "$MW/molly.yaml"
            rm -rf -- "$MW"
        done
        echo "E2E molly: passed $NPASS, failed $NFAIL, log-only differences $NLOG  ($(date +%FT%T))" ;;
    esac
done
[ "$NFAIL" = 0 ]
