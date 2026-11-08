#!/usr/bin/env bash
# WP6 acceptance: stages 00_tables .. 08_cross (lib/satpipe/st_00_tables.sh .. st_08_cross.sh) run for real on the toy
# project (tests/toy/make_toy.py) by the Python driver (python -m satpipe.cli, with its test-only trace hook) and by the
# bash driver (bin/satpipe), each on the SAME run-directory path (states moved aside), and compared:
#   exit code; stdout and trace (PLAN/CMD/OUT; argv[0] of Python's sys.executable normalised to python3, see PORT_SPEC 5b;
#   trace also compared sorted because 07_mafft runs in parallel); every file of the run directory byte for byte (cmp)
#   except logs/* (timestamps masked, python path normalised, line order of parallel logs kept), DONE (finished masked)
#   and manifest.json (timestamps masked, jq -S); then each driver reruns on the OTHER driver's run directory and must
#   report every stage "up to date, skipped" (the DONE records are interchangeable).
# Variants: base (rotated_consensus given, monomerize, strand plus), large (large_scaffold_bp small: the 02b chunk path
# with several chunks), given (inputs.trf_dat for spa, stale trf/*.dat removed, refs.rot derived by minrot, no Step 3b,
# strand both), unscanned (both TRF tables given; tools.bin_dir puts a trf wrapper first on PATH that hangs on chunk_2,
# so 02b's timeout lists that chunk in unscanned.tsv and the manifest).
# Usage: tests/port/test_stages_wp6.sh WORKDIR [VARIANT...]   env: PYTHON (python3), SIZE (toy size, 1200000), JOBS (4)
# Needs trf, parallel, vsearch, mafft, samtools, gawk, jq, timeout on PATH (servers, not Windows).
set -uo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
W=${1:?usage: test_stages_wp6.sh WORKDIR [VARIANT...]}; shift
VARIANTS=("$@"); [ ${#VARIANTS[@]} -gt 0 ] || VARIANTS=(base large given unscanned)
PY=${PYTHON:-python3}
export LC_ALL=C
unset SATPIPE_TRACE SATPIPE_DRY
mkdir -p "$W"; W=$(cd "$W" && pwd)
GAWK=${SATPIPE_GAWK:-gawk}
TOY=$W/toy
CONS=$ROOT/satpipe/stages/02_seed_expansion_scan/consensus
if [ ! -s "$TOY/genomes/spb.fna" ]; then
    "$PY" "$ROOT/tests/toy/make_toy.py" "$TOY" "$CONS" --seed "${SEED:-7}" --size "${SIZE:-1200000}" > "$W/make_toy.log" 2>&1 || { echo "make_toy failed"; exit 1; }
fi
printf '0\tCLsat1\n1\tCLsat2\n2\tCLsat3\n3\tCLsat4\n' > "$W/cluster_family.tsv"
RUN=$W/run

mkcfg() {   # mkcfg VARIANT > yaml
    local v=$1
    cat <<EOF
run_name: wp6_$v
out_dir: $RUN
genomes:
  spa: $TOY/genomes/spa.fna
  spb: $TOY/genomes/spb.fna
families:
  CLsat1: {consensus: $CONS/CLsat1.cons.fa, mono: 146}
  CLsat2: {consensus: $CONS/CLsat2.cons.fa, mono: 145}
  CLsat3: {consensus: $CONS/CLsat3.cons.fa, mono: 146}
  CLsat4: {consensus: $CONS/CLsat4.cons.fa, mono: 116}
family_assignment: $W/cluster_family.tsv
subfamily:
  consensuses: $TOY/subx
  members: {CLsat1: [CLsat1.1, CLsat1.2], CLsat2: [CLsat2], CLsat3: [CLsat3.1, CLsat3.2], CLsat4: [CLsat4]}
expansion: {jobs: ${JOBS:-4}}
EOF
    case $v in
        base)  printf 'references: {rotated_consensus: %s}\ndiscovery: {monomerize: true, strand: plus}\n' "$TOY/refs.rot" ;;
        large) printf 'references: {rotated_consensus: %s}\ndiscovery: {monomerize: true, strand: plus}\n' "$TOY/refs.rot"
               printf 'finalisation: {large_scaffold_bp: 200000, chunk_bp: 400000, overlap_bp: 30000, min_piece_bp: 20000}\n' ;;
        given) printf 'discovery: {monomerize: false, strand: both, cluster_id: 0.90}\ninputs:\n  trf_dat: {spa: %s}\n' "$W/given_spa.dat" ;;
        unscanned)
               printf 'references: {rotated_consensus: %s}\ninputs:\n  trf_dat: {spa: %s, spb: %s}\n' "$TOY/refs.rot" "$W/given_spa.dat" "$W/given_spb.dat"
               printf 'finalisation: {large_scaffold_bp: 200000, chunk_bp: 400000, overlap_bp: 30000, min_piece_bp: 400000, trf_timeout: 60}\n'
               printf 'tools: {bin_dir: %s}\n' "$W/fakebin" ;;
    esac
}

prestate() {   # the state both drivers start from
    rm -rf -- "$RUN"; rm -f -- "$TOY"/genomes/*.fai
    if [ "$1" = given ]; then mkdir -p "$RUN/trf"; echo stale > "$RUN/trf/old.dat"; echo stale > "$RUN/trf/.dat"; fi
}

norm() {   # timestamps, Python's interpreter path (Step 3b argv[0]) -> python3
    "$GAWK" -v PYX="$PYX" '{ gsub(/[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9]:[0-9][0-9]:[0-9][0-9]/, "TS")
        if (PYX != "") { while ((i = index($0, PYX)) > 0) $0 = substr($0, 1, i - 1) "python3" substr($0, i + length(PYX)) }
        print }'
}
PYX=$("$PY" -c 'import sys; print(sys.executable)')

NPASS=0; NFAIL=0
ok()   { NPASS=$((NPASS + 1)); printf 'ok   %s\n' "$1"; }
fail() { NFAIL=$((NFAIL + 1)); printf 'FAIL %s\n' "$1"; [ -n "${2:-}" ] && head -n 20 "$2" | cut -c1-300; }

drive() {   # drive SIDE CFG OUTPREFIX args...
    local side=$1 cfg=$2 o=$3 rc=0; shift 3
    : > "$o.trace"
    if [ "$side" = py ]; then
        ( cd "$W" && PYTHONPATH=$ROOT SATPIPE_TRACE="$o.trace" "$PY" -m satpipe.cli run "$cfg" "$@" ) > "$o.stdout" 2> "$o.stderr" || rc=$?
    else
        ( cd "$W" && SATPIPE_TRACE="$o.trace" "$ROOT/bin/satpipe" run "$cfg" "$@" ) > "$o.stdout" 2> "$o.stderr" || rc=$?
    fi
    echo "$rc" > "$o.rc"
}

compare_runs() {   # compare_runs VARIANT  (run dirs $W/$v.py.run and $W/$v.sh.run)
    local v=$1 a=$W/$1.py.run b=$W/$1.sh.run o=$W/out/$1 f k bad=0 n=0 nd=0
    ( cd "$a" && find . \( -type f -o -type l \) | LC_ALL=C sort ) > "$o.files.py"
    ( cd "$b" && find . \( -type f -o -type l \) | LC_ALL=C sort ) > "$o.files.sh"
    if diff "$o.files.py" "$o.files.sh" > "$o.diff.files"; then ok "$v files: same list ($(grep -c '' "$o.files.py") files)"; else fail "$v files: lists differ" "$o.diff.files"; fi
    : > "$o.diff.content"
    while IFS= read -r f; do
        f=${f#./}
        [ -e "$b/$f" ] || continue
        case $f in
            logs/satpipe.log)   # the driver's own log mirrors stdout: 07_mafft runs in parallel, so line order may differ
                diff <(norm < "$a/$f") <(norm < "$b/$f") > /dev/null \
                || diff <(norm < "$a/$f" | LC_ALL=C sort) <(norm < "$b/$f" | LC_ALL=C sort) >> "$o.diff.content" \
                || { echo "  log differs: $f" >> "$o.diff.content"; bad=1; } ;;
            logs/*) diff <(norm < "$a/$f") <(norm < "$b/$f") >> "$o.diff.content" || { echo "  log differs: $f" >> "$o.diff.content"; bad=1; } ;;
            */DONE) nd=$((nd + 1)); diff <(norm < "$a/$f") <(norm < "$b/$f") >> "$o.diff.content" || { echo "  DONE differs: $f" >> "$o.diff.content"; bad=1; } ;;
            manifest.json) diff <(jq -S . "$a/$f" | norm) <(jq -S . "$b/$f" | norm) >> "$o.diff.content" || { echo "  manifest differs" >> "$o.diff.content"; bad=1; } ;;
            *) n=$((n + 1)); if [ -L "$a/$f" ]; then [ "$(readlink "$a/$f")" = "$(readlink "$b/$f")" ] || { echo "  link differs: $f" >> "$o.diff.content"; bad=1; }
               else cmp -s "$a/$f" "$b/$f" || { echo "  cmp differs: $f" >> "$o.diff.content"; bad=1; }; fi ;;
        esac
    done < "$o.files.py"
    if [ $bad = 0 ]; then ok "$v content: $n output files cmp-identical, $nd DONE equal (finished masked), logs equal (masked), manifest jq -S equal"
    else fail "$v content" "$o.diff.content"; fi
}

mkdir -p "$W/out"
for v in "${VARIANTS[@]}"; do
    mkcfg "$v" > "$W/$v.yaml"
    if [ "$v" = given ] || [ "$v" = unscanned ]; then
        for c in spa spb; do
            if [ -s "$W/base.py.run/trf/$c.dat" ]; then cp "$W/base.py.run/trf/$c.dat" "$W/given_$c.dat"
            else printf 'Tandem Repeats Finder Program\n' > "$W/given_$c.dat"; fi
        done
    fi
    if [ "$v" = unscanned ]; then
        mkdir -p "$W/fakebin"
        printf '#!/bin/sh\ncase "$1" in *chunk_2.fa) exec sleep 75 ;; esac\nexec %s "$@"\n' "$(command -v trf)" > "$W/fakebin/trf"
        chmod +x "$W/fakebin/trf"
    fi
    for side in py sh; do
        prestate "$v"
        drive "$side" "$W/$v.yaml" "$W/out/$v.$side" --to 08_cross --skip-tool-check
        rm -rf -- "$W/$v.$side.run"; [ -d "$RUN" ] && mv "$RUN" "$W/$v.$side.run"
    done
    o=$W/out/$v
    exp=0; [ "$v" = large ] && exp=1      # large: Step 4 fails on the empty motif set (scaffolds left are motif-free), in both drivers
    if [ "$(cat "$o.py.rc")" = "$(cat "$o.sh.rc")" ] && [ "$(cat "$o.sh.rc")" = $exp ]; then ok "$v rc=$(cat "$o.sh.rc") (expected $exp)"; else fail "$v rc py=$(cat "$o.py.rc") sh=$(cat "$o.sh.rc")" "$o.sh.stderr"; fi
    if diff <(norm < "$o.py.stdout") <(norm < "$o.sh.stdout") > "$o.diff.stdout"; then ok "$v stdout ($(grep -c '' "$o.sh.stdout") lines)"
    else
        # 07_mafft prints its command lines from parallel jobs: compare that stage's lines as a set
        if diff <(norm < "$o.py.stdout" | LC_ALL=C sort) <(norm < "$o.sh.stdout" | LC_ALL=C sort) > "$o.diff.stdout.sorted"; then ok "$v stdout (same lines, parallel order differs)"
        else fail "$v stdout" "$o.diff.stdout"; fi
    fi
    if diff <(norm < "$o.py.stderr") <(norm < "$o.sh.stderr") > "$o.diff.stderr"; then ok "$v stderr"; else fail "$v stderr" "$o.diff.stderr"; fi
    if diff <(norm < "$o.py.trace") <(norm < "$o.sh.trace") > "$o.diff.trace"; then
        ok "$v trace in order: PLAN=$(grep -c '^PLAN' "$o.sh.trace") CMD=$(grep -c '^CMD' "$o.sh.trace") OUT=$(grep -c '^OUT' "$o.sh.trace")"
    elif diff <(norm < "$o.py.trace" | LC_ALL=C sort) <(norm < "$o.sh.trace" | LC_ALL=C sort) > "$o.diff.trace.sorted"; then
        ok "$v trace sorted (parallel order differs): PLAN=$(grep -c '^PLAN' "$o.sh.trace") CMD=$(grep -c '^CMD' "$o.sh.trace") OUT=$(grep -c '^OUT' "$o.sh.trace")"
    else fail "$v trace" "$o.diff.trace.sorted"; fi
    [ -d "$W/$v.py.run" ] && [ -d "$W/$v.sh.run" ] && compare_runs "$v"
    # interchangeable DONE records: each driver reruns on the other driver's run directory and skips every stage
    for pair in "py sh" "sh py"; do
        set -- $pair
        [ "$(cat "$o.$1.rc")" = 0 ] || continue
        rm -rf -- "$RUN"; cp -a "$W/$v.$1.run" "$RUN"
        drive "$2" "$W/$v.yaml" "$o.rerun_$2_on_$1" --to 08_cross --skip-tool-check
        k=$(grep -c 'up to date, skipped' "$o.rerun_$2_on_$1.stdout")
        if [ "$(cat "$o.rerun_$2_on_$1.rc")" = 0 ] && [ "$k" = 9 ] && ! grep -q '^\[.*\] running' "$o.rerun_$2_on_$1.stdout"; then ok "$v $2 driver on $1 run dir: 9/9 stages up to date"
        else fail "$v $2 driver on $1 run dir: $k skipped" "$o.rerun_$2_on_$1.stdout"; fi
        rm -rf -- "$RUN"
    done
done
echo "passed $NPASS, failed $NFAIL"
[ $NFAIL = 0 ]
