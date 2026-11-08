#!/usr/bin/env bash
# WP11 (b): dry-run traces (SATPIPE_DRY=1, SATPIPE_TRACE) of both drivers on COPIES of a finished Darevskia run directory.
# The source run dir and its config are only read. The config is rewritten with out_dir = WORKDIR/run and every relative
# path made absolute against the ORIGINAL config's directory; each driver gets a fresh copy of the run dir at the same
# path (dry stage bodies read files their commands would have made, and some rewrite files, e.g. seeds/<fam>.bed), and runs
# `satpipe run CFG --only <stage>` for 02b_large .. 21_anchor (02_trf cannot dry-run, PORT_SPEC 5b).
# Compared per stage: rc, stdout, stderr, PLAN/CMD/OUT records.
# Usage: tests/port/e2e_darevskia_trace.sh WORKDIR SRC_RUN_DIR SRC_CONFIG     env: PYTHON (python3 with PyYAML)
set -uo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
W=${1:?usage}; SRC=${2:?src run dir}; CFG0=${3:?src config}
PY=${PYTHON:-python3}; GAWK=${SATPIPE_GAWK:-gawk}
export LC_ALL=C
unset SATPIPE_TRACE SATPIPE_DRY SATPIPE_STAGE_FILES
mkdir -p "$W/out"; W=$(cd "$W" && pwd)
# shellcheck source=e2e_lib.sh
. "$ROOT/tests/port/e2e_lib.sh"
RUN=$W/run; CD=$(cd "$(dirname "$CFG0")" && pwd)
# out_dir -> copy; `key: ../x` / `key: ./x` / bare relative file values of path keys -> absolute against the config dir
"$GAWK" -v RUN="$RUN" -v CD="$CD" '
    /^out_dir:/ { print "out_dir: " RUN; next }
    { while (match($0, /(: |\{consensus: )\.\.?\/[^ ,}#]*/)) {
          pre = substr($0, 1, RSTART - 1); m = substr($0, RSTART, RLENGTH); rest = substr($0, RSTART + RLENGTH)
          k = index(m, " "); $0 = pre substr(m, 1, k) CD "/" substr(m, k + 1) rest } print }' "$CFG0" > "$W/darevskia.yaml"
echo "config differences (source -> copy):"; diff "$CFG0" "$W/darevskia.yaml" | sed 's/^/  /'
# one stage at a time (--only), each on a fresh copy: a dry stage body cannot feed the next one (06_cluster's body removes
# clusters_by_file, which Step6 would refill; 07_mafft then lists it), so a single dry pass over 02b..21 stops at 07.
STAGES=${STAGES:-02b_large 03_motifs 04_filter 05_rotate 06_cluster 07_mafft 08_cross 09_assign 10_seeds 11_expand 12_master 13_rescan 14_arbitrate 15_catalogue 16_regions 17_terminal 18_subfamily 19_orthology 20_candidates 21_anchor}
for sid in $STAGES; do
    for s in py sh; do
        rm -rf -- "$RUN"; mkdir -p "$RUN"
        # everything but the viewer outputs (not read by any stage)
        for x in "$SRC"/* ; do case ${x##*/} in viewer|viewer_flankscan) ;; *) cp -a "$x" "$RUN/" ;; esac; done
        export SATPIPE_DRY=1
        drive "$s" "$W/out/$sid.$s" run "$W/darevskia.yaml" --only "$sid"
        unset SATPIPE_DRY
        rm -rf -- "$RUN"
    done
    o=$W/out/$sid
    echo "== $sid: rc py=$(cat "$o.py.rc") sh=$(cat "$o.sh.rc")"
    if [ "$(cat "$o.py.rc")" = "$(cat "$o.sh.rc")" ]; then ok "$sid rc=$(cat "$o.sh.rc") both"; else fail "$sid rc differs" "$o.sh.stderr"; fi
    cmp_text "$sid stdout" "$o.py.stdout" "$o.sh.stdout" "$o.diff.stdout"
    if diff <(norm < "$o.py.stderr") <(norm < "$o.sh.stderr") > "$o.diff.stderr"; then ok "$sid stderr equal"
    elif grep -q '^Traceback' "$o.py.stderr" && diff <(tail -n 1 "$o.py.stderr" | norm) <(sed -E 's/^satpipe run: //' "$o.sh.stderr" | norm) > "$o.diff.stderr"; then
        ok "$sid stderr: same exception text (Python traceback vs 'satpipe run: <msg>', PORT_LOG WP5)"
    else fail "$sid stderr" "$o.diff.stderr"; fi
    for k in PLAN CMD OUT; do
        if diff <(grep "^$k" "$o.py.trace" | norm | LC_ALL=C sort) <(grep "^$k" "$o.sh.trace" | norm | LC_ALL=C sort) > "$o.diff.$k"; then ok "$sid $k records equal ($(grep -c "^$k" "$o.sh.trace"))"
        else fail "$sid $k records differ" "$o.diff.$k"; fi
    done
    if diff <(norm < "$o.py.trace") <(norm < "$o.sh.trace") > /dev/null; then ok "$sid trace in order"; else echo "     ($sid trace equal as a set, order differs)"; fi
done
echo "DAREVSKIA TRACE: passed $NPASS, failed $NFAIL"
[ "$NFAIL" = 0 ]
