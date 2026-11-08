#!/usr/bin/env bash
# WP8 acceptance: stages 16_regions .. 20_candidates of the bash driver (lib/satpipe/st_16..st_20, subx.sh) against the
# Python driver (satpipe/run.py with its test-only trace hook), both REAL drivers (`python -m satpipe.cli run` and
# bin/satpipe run). The stage scripts are replaced by a fake: tools.bin_dir holds `bash` and `python3` wrappers that, for
# the eleven stage scripts these stages call, write their output files from argv + the stage env variables (so a wrong
# argument or env value changes a fingerprint) and pass every other call to the real interpreter. Each step runs both
# drivers on the SAME run-directory path (states swapped in and out, as tests/port/test_driver.sh) and compares exit code,
# stdout, stderr, the trace (PLAN/CMD/OUT), every file of the run directory (DONE and logs with timestamps masked) and
# manifest.json. Variants: hapdup off/on, subfamily.consensuses multi-FASTA / directory / absent, seeds given for every
# family (20_candidates skip), discovery folder missing (error), zero candidates, missing hapdup.tsv (FileNotFoundError),
# dry mode, --plan.
# Usage: tests/port/test_stages_wp8.sh WORKDIR    env: PYTHON (python3)   needs gawk, jq, bash >= 4 (servers, not Windows)
set -uo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
W=${1:?usage: test_stages_wp8.sh WORKDIR}
PY=${PYTHON:-python3}
export LC_ALL=C
rm -rf -- "$W"; mkdir -p "$W/proj" "$W/out"
W=$(cd "$W" && pwd)
P=$W/proj
unset SATPIPE_TRACE SATPIPE_DRY SATPIPE_STAGE_FILES FAKE_ABSENT FAKE_NOHAPDUP SLURM_CPUS_PER_TASK
GAWK=${SATPIPE_GAWK:-gawk}
REAL_PY=$(command -v "$PY"); REAL_BASH=$(command -v bash)

# ------------------------------------------------------------------ project
mkdir -p "$P/genomes" "$P/ref/subx" "$P/configs" "$P/bin" "$P/seeds"
for c in spa spb; do printf '>%s_scf1\n%s\n' "$c" "$(printf 'ACGTTGCAAG%.0s' $(seq 1 30))" > "$P/genomes/$c.fna"; done
CONS=$ROOT/satpipe/stages/02_seed_expansion_scan/consensus
for f in CLsat1 CLsat2 CLsat3 CLsat4; do tr -d '\r' < "$CONS/$f.cons.fa" > "$P/ref/$f.cons.fa"; done
printf '# cluster\tfamily\n0\tCLsat1\n1\tCLsat2\n2\tCLsat3\n3\tCLsat4\n' > "$P/ref/cluster_family.tsv"
# multi-record FASTA with a description, wrapped lines, CRLF, a leading junk line and a repeated name (the later wins)
printf 'junk\r\n>CLsat1.1 desc x\r\nAAAG\r\nCT \r\n>CLsat2.1\nGGTTA\n>CLsat1.2\nTTTT\n>CLsat1.2\nCCCC\n' > "$P/ref/subs.fa"
printf '>CLsat1.1\nAAAGCTAAAGCT\n' > "$P/ref/subx/CLsat1.1.fa"; printf '>CLsat2.1\nGGTTAGGTTA\n' > "$P/ref/subx/CLsat2.1.fa"
for f in CLsat1 CLsat2 CLsat3 CLsat4; do printf 'spa_scf1\t0\t100\t%s\n' "$f" > "$P/seeds/$f.bed"; done

# fake interpreters: bash/python3 SCRIPT ARGS -> outputs of SCRIPT from argv + env; anything else -> the real one
cat > "$P/bin/fake" <<EOF
#!$REAL_BASH
self=\${0##*/}
case \$self in bash) real=$REAL_BASH ;; *) real=$REAL_PY ;; esac
[ \$# -gt 0 ] || exec "\$real"
name=\${1##*/}
envs="WINDOW_MONO=\${WINDOW_MONO-} EDGE_TOL=\${EDGE_TOL-} MIN_QCOV=\${MIN_QCOV-} MIN_IDENT=\${MIN_IDENT-} MIN_TOVL=\${MIN_TOVL-} THREADS=\${THREADS-} CONFIDENT_DELTA=\${CONFIDENT_DELTA-} MIN_BITS=\${MIN_BITS-} LC_ALL=\${LC_ALL-}"
w() { printf '%s %s | %s\n' "\$self" "\$*" "\$envs"; }   # one record line
shift
case \$name in
  build_locus_master.sh) { w "\$name" "\$@"; cat "\$1"; } > "\$3/locus_master.tsv" ;;
  hapdup_filter.sh)
    [ -n "\${FAKE_NOHAPDUP:-}" ] && exit 0
    { w "\$name" "\$@"; printf 'L1\tspa\tKEPT\nL2\tspb\tREMOVED\nL3 spb  REMOVED  \nL4\tREMOVED_NOT\r\nL5\tx\tREMOVED\rL6 REMOVED\n# REMOVED\n'; } > "\$3/hapdup.tsv"
    printf 'filtered\n' >> "\$1" ;;
  run_extract.sh) w "\$name" "\$@" > "\$1/manifest.tsv"; { w "\$name" "\$@"; cat "\$1/locus_master.tsv"; } > "\$1/all_regions.fa" ;;
  terminal_monomers.sh) w "\$name" "\$@" > "\$1/terminal_monomers.tsv" ;;
  boundary_phase.sh) w "\$name" "\$@" > "\$1/locus_boundary_phase.tsv" ;;
  build_json.py) w "\$name" "\$@" > "\$1/boundary_phase.json" ;;
  assign_subfamily.sh) { w "\$name" "\$@"; for f in "\$2"/*; do echo "== \${f##*/}"; cat "\$f"; done; } > "\$1/locus_subfamily.tsv" ;;
  build_subfam_term_json.py) { w "\$name" "\$@"; cat "\$1/locus_subfamily.tsv"; } > "\$1/subfamily.json"; w "\$name" "\$@" > "\$1/terminal_monomers.json" ;;
  rebuild_ortho_segmented.py)
    while [ \$# -gt 0 ]; do [ "\$1" = --out ] && o=\$2; shift; done
    w "\$name" > "\$o/flanks.json"; w "\$name" > "\$o/flank_orthologs.json"
    printf 'a\tb\nL1\tL2\nL3\tL4\r\nL5\tL6' > "\$o/ortho_pairs.tsv" ;;
  singleton_screen.sh)
    printf '# best\nspb|scf2:5-9\tCLsat1\t50.5\n\nspa|scf1:1-10\tCLsat2\t70\nspa|a:b:20-30\tCLsat1\t70.0\r\nspa|scf1:1-10\tCLsat2\t71\nshort\tx\n spb|z:1-2\tCLsat3\t99\n' > "\$3/sing_best.tsv"
    if [ -z "\${FAKE_ABSENT+x}" ]; then
      printf 'spb|scf2:5-9\n\n  spa|a:b:20-30  \nspa|scf1:1-10\nspb|scf9:3-4\nspa|scf0:3-4\r\nspb|z:1-2\n' > "\$3/candidates_absent.txt"
    elif [ -n "\$FAKE_ABSENT" ]; then printf '%s\n' "\$FAKE_ABSENT" > "\$3/candidates_absent.txt"; fi
    w "\$name" "\$@" > "\$3/screen.log" ;;
  candidate_scan.py) { w "\$name" "\$@"; cat "\$1"; } > "\$4/scan.txt" ;;
  candidate_flank5k.py) w "\$name" "\$@" > "\$1/candidates_summary.json" ;;
  candidate_subsubfamily.py) { w "\$name" "\$@"; for f in "\$4"/*; do echo "== \${f##*/}"; cat "\$f"; done; } >> "\$1/candidates_summary.json" ;;
  *) exec "\$real" "\$self_arg" "\$@" ;;
esac
EOF
# the passthrough must keep the original first argument
sed -i 's/exec "\$real" "\$self_arg" "\$@"/exec "$real" "$FAKE_FIRST" "$@"/; s/^shift$/FAKE_FIRST=$1; shift/' "$P/bin/fake"
chmod +x "$P/bin/fake"; ln -s fake "$P/bin/bash"; ln -s fake "$P/bin/python3"

cat > "$P/configs/run.yaml" <<'EOF'
run_name: wp8
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
terminal: {window_mono: 2}
orthology: {min_pid: 85, flank_bp: 700}
subfamily:
  consensuses: ../ref/subs.fa
  members: {CLsat1: [CLsat1.1, CLsat1.2], CLsat2: [CLsat2.1]}
  confident_delta: 12.5
tools:
  bin_dir: ../bin
EOF
C=$P/configs
sed 's#^tools:#hapdup:\n  pairs: [[spa, spb]]\n  min_qcov: 0.9\ntools:#' "$C/run.yaml" > "$C/hapdup.yaml"
sed 's#consensuses: ../ref/subs.fa#consensuses: ../ref/subx#' "$C/run.yaml" > "$C/subxdir.yaml"
"$GAWK" '/^subfamily:/ { skip = 1; next } skip && /^  / { next } { skip = 0; print }' "$C/run.yaml" > "$C/nosub.yaml"
sed 's#consensuses: ../ref/subs.fa#consensuses: ../ref/subs.fa\n  members: {}#; /members: {CLsat1/d' "$C/run.yaml" > "$C/nomembers.yaml"
cp "$C/run.yaml" "$C/seeds.yaml"
printf 'inputs:\n  seeds: {CLsat1: ../seeds/CLsat1.bed, CLsat2: ../seeds/CLsat2.bed, CLsat3: ../seeds/CLsat3.bed, CLsat4: ../seeds/CLsat4.bed}\n' >> "$C/seeds.yaml"
cp "$C/run.yaml" "$C/someseeds.yaml"
printf 'inputs:\n  seeds: {CLsat3: ../seeds/CLsat3.bed, CLsat1: ../seeds/CLsat1.bed}\n' >> "$C/someseeds.yaml"

# the state stages 00-15 leave behind (what 16-20 read)
seed_state() {
    local r=$1
    mkdir -p "$r/tables" "$r/catalogue" "$r/motifs/clean/clusters_by_file"
    printf 'spa\t%s\nspb\t%s\n' "$P/genomes/spa.fna" "$P/genomes/spb.fna" > "$r/tables/genomes.tsv"
    printf 'CLsat1\t146\nCLsat2\t145\n' > "$r/tables/families.tsv"
    printf 'CLsat1\tCLsat1.1\n' > "$r/tables/members.tsv"
    printf '>r\nACGT\n' > "$r/tables/refs.rot"
    printf 'spa_scf1\t0\t100\tL1\n' > "$r/catalogue/MASTER.merged.bed"
    printf 'spa_scf1\t0\t100\tL1\tspa\n' > "$r/catalogue/MASTER.loci.segmented.bed"
    printf 'L1\tCLsat1\t100\n' > "$r/catalogue/MASTER.final_assignment.tsv"
}
seed_state "$W/state_py"; seed_state "$W/state_sh"

# ------------------------------------------------------------------ one step, both drivers (as tests/port/test_driver.sh)
NPASS=0; NFAIL=0
mask() { "$GAWK" '{ gsub(/[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9]:[0-9][0-9]:[0-9][0-9]/, "TS"); print }'; }
snapshot() {
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
            ( cd "$W" && env "${envs[@]}" PYTHONPATH="$ROOT" SATPIPE_TRACE="$tr" "$PY" -m satpipe.cli "$@" ) > "$W/out/$name.$side.stdout" 2> "$W/out/$name.$side.stderr" || rc=$?
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
            # Python raises a non-StageError (traceback); bash prints `satpipe run: <its message>`
            # unpack messages: Python 3.14 adds ", got N" (3.13 does not); compare without it
            diff <(tail -n 1 "$W/out/$name.py.stderr" | sed -E 's/^[A-Za-z_.]*(Error|Exception): //; s/(\(expected [0-9]+), got [0-9]+\)/\1)/') \
                 <(sed -E 's/^satpipe run: //; s/(\(expected [0-9]+), got [0-9]+\)/\1)/' "$W/out/$name.sh.stderr") > "$W/out/$name.diff.$k" || ok=0
            continue
        fi
        if ! diff <(mask < "$W/out/$name.py.$k") <(mask < "$W/out/$name.sh.$k") > "$W/out/$name.diff.$k"; then ok=0; fi
    done
    diff -r "$W/out/$name.py.snap" "$W/out/$name.sh.snap" > "$W/out/$name.diff.snap" || ok=0
    local nplan ncmd nout ndone
    nplan=$(grep -c '^PLAN' "$W/out/$name.sh.trace"); ncmd=$(grep -c '^CMD' "$W/out/$name.sh.trace"); nout=$(grep -c '^OUT' "$W/out/$name.sh.trace")
    ndone=$(grep -c '/DONE$' "$W/out/$name.sh.snap/_files" 2>/dev/null || true)
    if [ $ok = 1 ]; then NPASS=$((NPASS + 1)); printf 'ok   %-18s rc=%s PLAN=%s CMD=%s OUT=%s DONE=%s files=%s\n' "$name" "$(cat "$W/out/$name.sh.rc")" "$nplan" "$ncmd" "$nout" "$ndone" "$(grep -c '' "$W/out/$name.sh.snap/_files" 2>/dev/null || echo 0)"
    else NFAIL=$((NFAIL + 1)); printf 'FAIL %-18s (see %s/out/%s.diff.*)\n' "$name" "$W" "$name"; head -n 12 "$W"/out/"$name".diff.* | cut -c1-400; fi
}
# a check of one file's content in the bash state (the expected values the fake data were built for)
expect() {   # expect LABEL FILE EXPECTED-TEXT
    if [ "$(cat "$W/state_sh/$2" 2>&1)" = "$3" ]; then NPASS=$((NPASS + 1)); echo "ok   $1"
    else NFAIL=$((NFAIL + 1)); echo "FAIL $1: $W/state_sh/$2 is:"; cat "$W/state_sh/$2"; fi
}

R="run --skip-tool-check --from 16_regions --to 20_candidates"
step base        -- $R proj/configs/run.yaml
expect "candidates order" candidates/candidates.tsv "$(printf 'C1\tspa\tscf1\t1\t10\nC2\tspa\ta:b\t20\t30\nC3\tspb\tscf2\t5\t9\nC4\tspa\tscf0\t3\t4\nC5\tspb\tscf9\t3\t4\nC6\tspb\tz\t1\t2')"
expect "subx files" tables/subx/CLsat1.2.fa "$(printf '>CLsat1.2\nCCCCCCCC')"
step rerun       -- $R proj/configs/run.yaml
step hapdup_on   -- $R proj/configs/hapdup.yaml
if [ "$(jq .hapdup_removed "$W/state_sh/manifest.json")" = 5 ]; then NPASS=$((NPASS + 1)); echo "ok   hapdup_removed = 5"
else NFAIL=$((NFAIL + 1)); echo "FAIL hapdup_removed: $(jq .hapdup_removed "$W/state_sh/manifest.json")"; fi
step hapdup_thr SLURM_CPUS_PER_TASK=3 -- $R proj/configs/hapdup.yaml --force --only 16_regions
step hapdup_off  -- $R proj/configs/run.yaml
step subx_dir    -- $R proj/configs/subxdir.yaml
step no_subx     -- $R proj/configs/nosub.yaml
step no_members  -- $R proj/configs/nomembers.yaml
step zero_cands FAKE_ABSENT=' ' -- $R proj/configs/run.yaml --force --only 20_candidates
step no_absent_file FAKE_ABSENT= -- $R proj/configs/run.yaml --force --only 20_candidates
TRACEBACK=1 step bad_cand FAKE_ABSENT='spa|nocolon' -- $R proj/configs/run.yaml --force --only 20_candidates
TRACEBACK=1 step bad_cand2 FAKE_ABSENT='spa|s:1-2-3' -- $R proj/configs/run.yaml --force --only 20_candidates
step cands_again -- $R proj/configs/run.yaml --force --only 20_candidates
step dry SATPIPE_DRY=1 -- $R proj/configs/hapdup.yaml
for s in py sh; do rm -f "$W/state_$s/regions/hapdup.tsv"; done
TRACEBACK=1 step no_hapdup_tsv FAKE_NOHAPDUP=1 -- $R proj/configs/hapdup.yaml --force --only 16_regions
# --plan (bash only) equals the PLAN records Python traced in an up-to-date rerun
step rerun_hd    -- $R proj/configs/hapdup.yaml
mv "$W/state_sh" "$P/runs/t"
( cd "$W" && "$ROOT/bin/satpipe" $R proj/configs/hapdup.yaml --plan ) > "$W/out/plan.sh.stdout" 2> "$W/out/plan.sh.stderr"
mv "$P/runs/t" "$W/state_sh"
if grep '^PLAN' "$W/out/rerun_hd.py.trace" | diff - "$W/out/plan.sh.stdout" > "$W/out/plan.diff" && [ ! -s "$W/out/plan.sh.stderr" ]; then
    NPASS=$((NPASS + 1)); echo "ok   plan_only          PLAN=$(grep -c '' "$W/out/plan.sh.stdout")"
else NFAIL=$((NFAIL + 1)); echo "FAIL plan_only"; head "$W/out/plan.diff" "$W/out/plan.sh.stderr"; fi
# discovery did not run here: seeds for every family -> skip message; seeds for some -> StageError
for s in py sh; do rm -rf "$W/state_$s/motifs"; done
step seeds_given -- $R proj/configs/seeds.yaml
step seeds_some  -- $R proj/configs/someseeds.yaml --force --only 20_candidates
step no_seeds    -- $R proj/configs/run.yaml --force --only 20_candidates
echo "passed $NPASS, failed $NFAIL"
[ $NFAIL = 0 ]
