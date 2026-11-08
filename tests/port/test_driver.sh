#!/usr/bin/env bash
# WP5 acceptance: the bash driver (bin/satpipe + lib/satpipe/run.sh, stages.sh, manifest.sh) against the Python driver
# (satpipe/run.py with its test-only trace hook) on stub stages (tests/port/stub_stages.sh = tests/port/stub_run.py) and the
# real stage 21_anchor. Every step runs both drivers on the SAME run-directory path (the two states are swapped in and out)
# and compares: exit code, stdout, stderr, the trace (PLAN/CMD/OUT), every file of the run directory (DONE files and logs
# with timestamps masked) and manifest.json (byte layout with timestamps masked, and `jq -S`).
# Usage: tests/port/test_driver.sh WORKDIR      env: PYTHON (python3)   needs gawk, jq, bash >= 4 (servers, not Windows)
set -uo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
W=${1:?usage: test_driver.sh WORKDIR}
PY=${PYTHON:-python3}
export LC_ALL=C
rm -rf -- "$W"; mkdir -p "$W/proj" "$W/out"
W=$(cd "$W" && pwd)
P=$W/proj
unset SATPIPE_TRACE SATPIPE_DRY STUB_FAIL STUB_TAG

# ------------------------------------------------------------------ project
mkdir -p "$P/genomes" "$P/ref" "$P/configs" "$P/bin" "$P/seeds"
for c in spa spb; do printf '>%s_scf1\n%s\n' "$c" "$(printf 'ACGTTGCAAG%.0s' $(seq 1 30))" > "$P/genomes/$c.fna"; done
CONS=$ROOT/satpipe/stages/02_seed_expansion_scan/consensus
for f in CLsat1 CLsat2 CLsat3 CLsat4; do tr -d '\r' < "$CONS/$f.cons.fa" > "$P/ref/$f.cons.fa"; done
seqof() { grep -v '>' "$P/ref/$1.cons.fa" | tr -d '\n'; }
s1=$(seqof CLsat1); s2=$(seqof CLsat2)
printf '>CLsat1.1\n%s\n>CLsat2.1\n%s\n' "${s1:40}${s1:0:40}" "${s2:10}${s2:0:10}" > "$P/ref/subs.fa"
printf '# cluster\tfamily\n0\tCLsat1\n1\tCLsat2\n2\tCLsat3\n3\tCLsat4\n' > "$P/ref/cluster_family.tsv"
printf '#!/bin/sh\necho "fake ssearch36 $*"\nfor f in "$@"; do [ -f "$f" ] && cat "$f"; done\nexit 0\n' > "$P/bin/ssearch36"; chmod +x "$P/bin/ssearch36"
for f in CLsat1 CLsat2 CLsat3 CLsat4; do printf 'spa_scf1\t0\t100\t%s\n' "$f" > "$P/seeds/$f.bed"; done
cat > "$P/configs/run.yaml" <<'EOF'
run_name: wp5
out_dir: ../runs/t
genomes:
  spa: ../genomes/spa.fna
  spb: ../genomes/spb.fna
families:
  CLsat1: {consensus: ../ref/CLsat1.cons.fa, mono: 146}
  CLsat2: {consensus: ../ref/CLsat2.cons.fa, mono: 145, colour: "#ff0000"}
  CLsat3: {consensus: ../ref/CLsat3.cons.fa, mono: 146}
  CLsat4: {consensus: ../ref/CLsat4.cons.fa, mono: 116}
family_assignment: ../ref/cluster_family.tsv
expansion: {min_bits: 50.0}
discovery: {cluster_id: 0.90}
subfamily:
  consensuses: ../ref/subs.fa
  members: {CLsat1: [CLsat1.1], CLsat2: [CLsat2.1]}
anchor: {family: CLsat1, motif: AAAGAAG}
tools:
  bin_dir: ../bin
  ssearch36: ../bin/ssearch36
EOF
export SATPIPE_STAGE_FILES=$ROOT/tests/port/stub_stages.sh     # read by the bash driver only

# ------------------------------------------------------------------ one step, both drivers
NPASS=0; NFAIL=0
mask() { "$GAWK" '{ gsub(/[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9]:[0-9][0-9]:[0-9][0-9]/, "TS"); print }'; }
GAWK=${SATPIPE_GAWK:-gawk}
snapshot() {   # snapshot RUNDIR OUTDIR: every file, timestamps masked; manifest also as jq -S
    local r=$1 o=$2 f rel
    mkdir -p "$o"
    [ -d "$r" ] || { echo "(no run dir)" > "$o/_none"; return 0; }
    ( cd "$r" && find . -type f | LC_ALL=C sort ) > "$o/_files"
    while IFS= read -r rel; do
        f=$r/${rel#./}; mkdir -p "$o/$(dirname -- "$rel")"
        mask < "$f" > "$o/$rel"
    done < "$o/_files"
    if [ -f "$r/manifest.json" ]; then jq -S . "$r/manifest.json" | mask > "$o/_manifest.jqS"; fi
}
step() {   # step NAME [VAR=value ...] -- satpipe-args...
    local name=$1 side rc; shift
    local -a envs=()
    while [ "$1" != -- ]; do envs+=("$1"); shift; done; shift
    for side in py sh; do
        rm -rf -- "$P/runs"; mkdir -p "$P/runs"
        [ -d "$W/state_$side" ] && mv "$W/state_$side" "$P/runs/t"
        local tr=$W/out/$name.$side.trace; : > "$tr"
        rc=0
        if [ $side = py ]; then
            ( cd "$W" && env "${envs[@]}" SATPIPE_TRACE="$tr" "$PY" "$ROOT/tests/port/stub_run.py" "$@" ) > "$W/out/$name.$side.stdout" 2> "$W/out/$name.$side.stderr" || rc=$?
        else
            ( cd "$W" && env "${envs[@]}" SATPIPE_TRACE="$tr" "$ROOT/bin/satpipe" "$@" ) > "$W/out/$name.$side.stdout" 2> "$W/out/$name.$side.stderr" || rc=$?
        fi
        echo "$rc" > "$W/out/$name.$side.rc"
        rm -rf -- "$W/out/$name.$side.snap"; snapshot "$P/runs/t" "$W/out/$name.$side.snap"
        [ -d "$P/runs/t" ] && mv "$P/runs/t" "$W/state_$side"
    done
    local ok=1 k
    for k in rc stdout stderr trace; do
        if [ $k = stderr ] && [ "${TRACEBACK:-0}" = 1 ]; then
            # expected difference: Python raises a non-StageError (a traceback); bash prints `satpipe run: <its message>`
            # (list.index's ValueError text depends on the Python version: only the class is compared for it)
            diff <(tail -n 1 "$W/out/$name.py.stderr" | sed -E 's/^ValueError: .*/ValueError/; s/^[A-Za-z_.]*(Error|Exception): //') \
                 <(sed -E 's/^satpipe run: ValueError: .*/ValueError/; s/^satpipe run: //' "$W/out/$name.sh.stderr") > "$W/out/$name.diff.$k" || ok=0
            continue
        fi
        if ! diff <(mask < "$W/out/$name.py.$k") <(mask < "$W/out/$name.sh.$k") > "$W/out/$name.diff.$k"; then ok=0; fi
    done
    diff -r "$W/out/$name.py.snap" "$W/out/$name.sh.snap" > "$W/out/$name.diff.snap" || ok=0
    local nplan ncmd nout ndone
    nplan=$(grep -c '^PLAN' "$W/out/$name.sh.trace"); ncmd=$(grep -c '^CMD' "$W/out/$name.sh.trace"); nout=$(grep -c '^OUT' "$W/out/$name.sh.trace")
    ndone=$(grep -c '/DONE$' "$W/out/$name.sh.snap/_files" 2>/dev/null || true)
    if [ $ok = 1 ]; then NPASS=$((NPASS + 1)); printf 'ok   %-22s rc=%s PLAN=%s CMD=%s OUT=%s DONE=%s files=%s\n' "$name" "$(cat "$W/out/$name.sh.rc")" "$nplan" "$ncmd" "$nout" "$ndone" "$(grep -c '' "$W/out/$name.sh.snap/_files" 2>/dev/null || echo 0)"
    else NFAIL=$((NFAIL + 1)); printf 'FAIL %-22s (see %s/out/%s.diff.*)\n' "$name" "$W" "$name"; head -n 12 "$W"/out/"$name".diff.* | cut -c1-400; fi
}

CFG=proj/configs/run.yaml
step full        -- run "$CFG" --skip-tool-check
step rerun       -- run "$CFG" --skip-tool-check
# --plan (bash only): the PLAN records equal those Python traced in the up-to-date rerun, and nothing is written
plan_check() {
    mv "$W/state_sh" "$P/runs/t"; local before after
    before=$(cd "$P/runs/t" && find . -type f -newer "$P/configs/run.yaml" | LC_ALL=C sort | md5sum)
    ( cd "$W" && "$ROOT/bin/satpipe" run "$CFG" --skip-tool-check --plan ) > "$W/out/plan.sh.stdout" 2> "$W/out/plan.sh.stderr"
    after=$(cd "$P/runs/t" && find . -type f -newer "$P/configs/run.yaml" | LC_ALL=C sort | md5sum)
    mv "$P/runs/t" "$W/state_sh"
    grep '^PLAN' "$W/out/rerun.py.trace" | diff - "$W/out/plan.sh.stdout" > "$W/out/plan.diff" && [ ! -s "$W/out/plan.sh.stderr" ] && [ "$before" = "$after" ]
}
if plan_check; then NPASS=$((NPASS + 1)); echo "ok   plan_only              PLAN=$(grep -c '' "$W/out/plan.sh.stdout") (= Python rerun trace)"
else NFAIL=$((NFAIL + 1)); echo "FAIL plan_only"; head "$W/out/plan.diff" "$W/out/plan.sh.stderr"; fi
step force_only  -- run "$CFG" --skip-tool-check --force --only 05_rotate,12_master
step from_to     -- run "$CFG" --skip-tool-check --force --from 10_seeds --to 13_rescan
step only_anchor -- run "$CFG" --skip-tool-check --force --only 21_anchor,03_motifs
printf '>spb_scf2\nACGTACGTAC\n' >> "$P/genomes/spb.fna"
step input_changed -- run "$CFG" --skip-tool-check
sed -i 's/min_bits: 50.0/min_bits: 60/' "$P/configs/run.yaml"
step params_changed -- run "$CFG" --skip-tool-check
s3=$(seqof CLsat3); printf '>CLsat3\n%s\n' "${s3:5}${s3:0:5}" > "$P/ref/CLsat3.cons.fa"
step cons_changed -- run "$CFG" --skip-tool-check --to 09_assign
step stage_fails STUB_FAIL=13_rescan -- run "$CFG" --skip-tool-check --force --from 12_master
step after_fail  -- run "$CFG" --skip-tool-check
step anchor_force -- run "$CFG" --skip-tool-check --only 21_anchor --force
step unknown_from -- run "$CFG" --skip-tool-check --from 99_nothing
step unknown_to   -- run "$CFG" --skip-tool-check --to 5
TRACEBACK=1 step to_before_from -- run "$CFG" --skip-tool-check --from 10_seeds --to 05_rotate
step only_empty  -- run "$CFG" --skip-tool-check --only ''
step dry SATPIPE_DRY=1 -- run "$CFG" --skip-tool-check --from 09_assign --to 11_expand
cp "$P/configs/run.yaml" "$P/configs/seeds.yaml"
printf 'inputs:\n  seeds: {CLsat1: ../seeds/CLsat1.bed, CLsat2: ../seeds/CLsat2.bed, CLsat3: ../seeds/CLsat3.bed, CLsat4: ../seeds/CLsat4.bed}\n' >> "$P/configs/seeds.yaml"
step seeds_given -- run proj/configs/seeds.yaml --skip-tool-check --to 10_seeds
# stage 21 fails (motif absent): StageError text, manifest written
sed 's/motif: AAAGAAG/motif: CCCCCCCCCC/' "$P/configs/run.yaml" > "$P/configs/badmotif.yaml"
TRACEBACK=1 step anchor_error -- run proj/configs/badmotif.yaml --skip-tool-check --only 21_anchor
# a project that fails the check: nothing written
step no_config   -- run proj/configs/none.yaml
printf 'run_name: x\nout_dir: ../runs/t\n' > "$P/configs/invalid.yaml"
step invalid_cfg -- run proj/configs/invalid.yaml --skip-tool-check
sed 's#../genomes/spb.fna#../genomes/missing.fna#' "$P/configs/run.yaml" > "$P/configs/missing.yaml"
step missing_file -- run proj/configs/missing.yaml --skip-tool-check
step missing_file_tools -- run proj/configs/missing.yaml
rm -rf "$W/state_py" "$W/state_sh"
step empty_full  -- run "$CFG" --skip-tool-check --only 00_tables,21_anchor
echo "passed $NPASS, failed $NFAIL"
[ $NFAIL = 0 ]
