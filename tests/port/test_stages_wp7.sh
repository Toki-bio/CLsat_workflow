#!/usr/bin/env bash
# WP7 acceptance: stages 09_assign .. 15_catalogue of the bash driver (lib/satpipe/st_09..15*.sh, expand_cmd.sh, assign.awk,
# manual.awk) against the Python driver (satpipe/run.py with its test-only trace hook), REAL runs on the synthetic toy
# (tests/toy/make_toy.py), modelled on tests/port/test_driver.sh: every step runs both drivers on the SAME run-directory path
# (the two states are swapped in and out) and compares exit code, stdout, stderr, the trace (PLAN/CMD/OUT), every file of the
# run directory except logs/ (timestamps masked: DONE `finished`, satpipe.log) and manifest.json (`jq -S`). logs/ is diffed
# too (timestamps masked) and reported separately: the stage scripts write tool output there.
# Stages 00-08 are run ONCE with the Python driver (WP6 ports them); every step starts from a copy of that run directory or
# of the state a previous step left, as stated per step.
# Usage: tests/port/test_stages_wp7.sh WORKDIR    env: PYTHON (python3), SEED (7), SIZE (1200000), JOBS (8)
# Needs the stage tools on PATH (sinederella env: ssearch36 vsearch mafft trf samtools bedtools parallel gawk) + jq; bash >= 4.
set -uo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
W=${1:?usage: test_stages_wp7.sh WORKDIR}
PY=${PYTHON:-python3}
export LC_ALL=C
GAWK=${SATPIPE_GAWK:-gawk}
rm -rf -- "$W"; mkdir -p "$W/proj" "$W/out" "$W/tmp"
W=$(cd "$W" && pwd)
export TMPDIR=$W/tmp
P=$W/proj
unset SATPIPE_TRACE SATPIPE_DRY SATPIPE_STAGE_FILES
CONS=$ROOT/satpipe/stages/02_seed_expansion_scan/consensus
FAMS="CLsat1 CLsat2 CLsat3 CLsat4"

# ------------------------------------------------------------------ toy project + config writer
"$PY" "$ROOT/tests/toy/make_toy.py" "$P/toy" "$CONS" --seed "${SEED:-7}" --size "${SIZE:-1200000}" > "$W/out/make_toy.log" 2>&1 || { echo "make_toy failed"; exit 1; }
mkdir -p "$P/ref" "$P/configs" "$P/seeds"
# cfg NAME [extra YAML lines...]   writes $P/configs/NAME.yaml: the toy config (run_driver_toy.sh:17-33) + extra top-level lines
cfgw() {
    local name=$1; shift
    {
        cat <<EOF
run_name: toy
out_dir: $P/runs/t
genomes:
  spa: $P/toy/genomes/spa.fna
  spb: $P/toy/genomes/spb.fna
families:
  CLsat1: {consensus: $CONS/CLsat1.cons.fa, mono: 146}
  CLsat2: {consensus: $CONS/CLsat2.cons.fa, mono: 145}
  CLsat3: {consensus: $CONS/CLsat3.cons.fa, mono: 146}
  CLsat4: {consensus: $CONS/CLsat4.cons.fa, mono: 116}
references: {rotated_consensus: $P/toy/refs.rot}
subfamily:
  consensuses: $P/toy/subx
  members: {CLsat1: [CLsat1.1, CLsat1.2], CLsat2: [CLsat2], CLsat3: [CLsat3.1, CLsat3.2], CLsat4: [CLsat4]}
discovery: {monomerize: true, strand: plus}
finalisation: {large_scaffold_bp: 50000000}
EOF
        local l has_fa=0 has_exp=0
        for l in "$@"; do
            case $l in family_assignment:*) has_fa=1 ;; expansion:*) has_exp=1 ;; esac
            printf '%s\n' "$l"
        done
        [ $has_fa = 1 ] || printf 'family_assignment: %s\n' "$P/ref/cluster_family.tsv"
        [ $has_exp = 1 ] || printf 'expansion: {jobs: %s}\n' "${JOBS:-8}"
    } > "$P/configs/$name.yaml"
}

# ------------------------------------------------------------------ stages 00-08 once (Python driver)
printf '0\tCLsat1\n1\tCLsat2\n2\tCLsat3\n3\tCLsat4\n' > "$P/ref/cluster_family.tsv"   # placeholder so `check` passes
cfgw base
( cd "$W" && PYTHONPATH=$ROOT "$PY" -m satpipe.cli run "$P/configs/base.yaml" --to 08_cross ) > "$W/out/base08.log" 2>&1 \
    || { echo "Python run to 08_cross failed: $W/out/base08.log"; tail -n 20 "$W/out/base08.log"; exit 1; }
# the curated table from the Step-8 .uc, as run_driver_toy.sh:38 (sorted by cluster number inside gawk)
"$GAWK" -F '\t' '($1 == "S" || $1 == "H") && $9 ~ /^CLsat[1-4]$/ { t[$2 + 0] = $2 "\t" $9 }
    END { PROCINFO["sorted_in"] = "@ind_num_asc"; for (k in t) print t[k] }' "$P/runs/t/cross/cross_species.uc" > "$P/ref/cluster_family.tsv"
echo "cluster_family.tsv: $(tr '\t\n' ': ' < "$P/ref/cluster_family.tsv")"
mv "$P/runs/t" "$W/base08"

# ------------------------------------------------------------------ one step, both drivers
NPASS=0; NFAIL=0
mask() { "$GAWK" '{ gsub(/[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9]:[0-9][0-9]:[0-9][0-9]/, "TS"); print }'; }
snapshot() {   # snapshot RUNDIR OUTDIR: every file, timestamps masked (logs/ under OUTDIR/logs); manifest also as jq -S
    local r=$1 o=$2 f rel
    mkdir -p "$o/files" "$o/logs"
    [ -d "$r" ] || { echo "(no run dir)" > "$o/files/_none"; return 0; }
    ( cd "$r" && find . -type f | LC_ALL=C sort ) > "$o/_files"
    while IFS= read -r rel; do
        f=$r/${rel#./}
        case $rel in ./logs/*) mkdir -p "$o/logs/$(dirname -- "$rel")"; mask < "$f" > "$o/logs/$rel" ;;
                     *) mkdir -p "$o/files/$(dirname -- "$rel")"; mask < "$f" > "$o/files/$rel" ;; esac
    done < "$o/_files"
    if [ -f "$r/manifest.json" ]; then jq -S . "$r/manifest.json" | mask > "$o/_manifest.jqS"; fi
}
init() { rm -rf -- "$W/state_py" "$W/state_sh"; cp -a "$1" "$W/state_py"; cp -a "$1" "$W/state_sh"; }
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
            ( cd "$W" && env "${envs[@]}" SATPIPE_TRACE="$tr" PYTHONPATH="$ROOT" "$PY" -m satpipe.cli "$@" ) > "$W/out/$name.$side.stdout" 2> "$W/out/$name.$side.stderr" || rc=$?
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
            # a Python exception other than StageError: Python prints a traceback, bash `satpipe run: <its message>`
            diff <(tail -n 1 "$W/out/$name.py.stderr" | sed -E 's/^[A-Za-z_.]*(Error|Exception): //') \
                 <(sed -E 's/^satpipe run: //' "$W/out/$name.sh.stderr") > "$W/out/$name.diff.$k" || ok=0
            continue
        fi
        if ! diff <(mask < "$W/out/$name.py.$k") <(mask < "$W/out/$name.sh.$k") > "$W/out/$name.diff.$k"; then ok=0; fi
    done
    diff -r "$W/out/$name.py.snap/files" "$W/out/$name.sh.snap/files" > "$W/out/$name.diff.snap" || ok=0
    diff "$W/out/$name.py.snap/_manifest.jqS" "$W/out/$name.sh.snap/_manifest.jqS" > "$W/out/$name.diff.jqS" 2>&1 || ok=0
    local nlog nplan ncmd nout ndone
    nlog=$(diff -rq "$W/out/$name.py.snap/logs" "$W/out/$name.sh.snap/logs" | grep -c '' || true)
    diff -r "$W/out/$name.py.snap/logs" "$W/out/$name.sh.snap/logs" > "$W/out/$name.diff.logs"
    nplan=$(grep -c '^PLAN' "$W/out/$name.sh.trace"); ncmd=$(grep -c '^CMD' "$W/out/$name.sh.trace"); nout=$(grep -c '^OUT' "$W/out/$name.sh.trace")
    ndone=$(grep -c '/DONE$' "$W/out/$name.sh.snap/_files" 2>/dev/null || true)
    if [ $ok = 1 ]; then NPASS=$((NPASS + 1)); printf 'ok   %-18s rc=%s PLAN=%s CMD=%s OUT=%s DONE=%s files=%s logs_differing=%s\n' "$name" "$(cat "$W/out/$name.sh.rc")" "$nplan" "$ncmd" "$nout" "$ndone" "$(grep -c '' "$W/out/$name.sh.snap/_files" 2>/dev/null || echo 0)" "$nlog"
    else NFAIL=$((NFAIL + 1)); printf 'FAIL %-18s (see %s/out/%s.diff.*)\n' "$name" "$W" "$name"; head -n 12 "$W"/out/"$name".diff.{rc,stdout,stderr,trace,snap,jqS} 2>/dev/null | cut -c1-400; fi
}
save() { rm -rf -- "$W/$1"; cp -a "$W/state_py" "$W/$1"; }   # a finished state (equal on both sides when the step passed)

C=$P/configs
# ------------------------------------------------------------------ 1. the toy through 09-15 (from the Python 00-08 run)
init "$W/base08"
step full         -- run "$C/base.yaml" --from 09_assign --to 15_catalogue
step rerun        -- run "$C/base.yaml" --from 09_assign --to 15_catalogue          # every stage up to date (DONE read back)
save base15
step dry   SATPIPE_DRY=1 -- run "$C/base.yaml" --from 09_assign --to 15_catalogue  # trace only, on the finished state

# ------------------------------------------------------------------ 2. config texts reaching argv/env (dry, from base15)
cfgw floats 'expansion: {min_bits: 50.0, min_aln_frac: 0.70, to_scaffold_ends: true, jobs: 4}' 'scan: {min_bits: 0.0, win_mono: 3}' \
            'arbitration: {delta: 20.0, delta_ambig: 15, bin_min_bits: 40, short_segments: merge}'
init "$W/base15"; step floats SATPIPE_DRY=1 -- run "$C/floats.yaml" --from 11_expand --to 15_catalogue
cfgw bmbnull 'expansion: {min_bits: 60.0, jobs: 4}'                # arbitration.bin_min_bits None -> BIN_MIN_BITS=60.0 (run.py:827)
init "$W/base15"; step bmb_null SATPIPE_DRY=1 -- run "$C/bmbnull.yaml" --only 14_arbitrate
cfgw bmbreal 'arbitration: {bin_min_bits: 40.0}'
init "$W/base15"; step bmb_real -- run "$C/bmbreal.yaml" --only 14_arbitrate,15_catalogue

# ------------------------------------------------------------------ 3. the curated table (real 09, from base15)
T=$P/ref/cluster_family.tsv
{ printf '# curated\r\n\r\n'; "$GAWK" -F '\t' '{ printf " %s \t%s\tnote\r\n", $1, $2 }' "$T"; } > "$P/ref/t_crlf.tsv"
cfgw tcrlf "family_assignment: $P/ref/t_crlf.tsv"
init "$W/base15"; step table_crlf -- run "$C/tcrlf.yaml" --only 09_assign --force
{ "$GAWK" -F '\t' '$2 != "CLsat4"' "$T"; printf '999\tCLsat1\n7\tCLsatX\n'; } > "$P/ref/t_bad.tsv"
cfgw tbad "family_assignment: $P/ref/t_bad.tsv"
init "$W/base15"; step table_bad -- run "$C/tbad.yaml" --only 09_assign
"$GAWK" -F '\t' -v OFS='\t' '{ if ($2 == "CLsat1") $2 = "CLsat2"; else if ($2 == "CLsat2") $2 = "CLsat1"; print }' "$T" > "$P/ref/t_swap.tsv"
cfgw tswap "family_assignment: $P/ref/t_swap.tsv"
init "$W/base15"; step table_swap -- run "$C/tswap.yaml" --only 09_assign
"$GAWK" -F '\t' -v OFS='\t' '{ if ($2 == "CLsat3") $2 = "CLsat4"; else if ($2 == "CLsat4") $2 = "CLsat3"; print }' "$T" > "$P/ref/t_swap34.tsv"
cfgw tswap34 "family_assignment: $P/ref/t_swap34.tsv"
init "$W/base15"; step table_swap34 -- run "$C/tswap34.yaml" --only 09_assign
# a cluster number that Step 8 did not make (check cannot see it; the table's other errors are refused by check first)
{ cat "$T"; printf '999\tCLsat1\n'; } > "$P/ref/t_999.tsv"
cfgw t999 "family_assignment: $P/ref/t_999.tsv"
init "$W/base15"; step table_999 -- run "$C/t999.yaml" --only 09_assign
# a fake ssearch36 (tools.ssearch36; --skip-tool-check) for the best-hit rule: comment lines, short lines, a bits field that is
# not a float, ties (first kept), a diverged mismatch (WARNING, run continues), a pid that is not a float (0.0), %.1f rounding
mkdir -p "$P/bin"
cat > "$P/bin/fake_ssearch36" <<'EOF'
#!/bin/sh
q=; for a in "$@"; do case $a in *centroids_query.fa) q=$a ;; esac; done
echo "# fake ssearch36"
gawk -F '\t' 'BEGIN { split("CLsat1 CLsat2 CLsat3 CLsat4", f, " ") }
    /^>/ { c = substr($1, 2); n++
           printf "%s\t%s\t%s\t1\t1\t1\t1\t1\t1\t1\t1\t%s\n", c, f[(n % 4) + 1], "80.00", "150.25"
           printf "%s\t%s\t%s\t1\t1\t1\t1\t1\t1\t1\t1\t%s\n", c, f[((n + 1) % 4) + 1], "99.0", "150.25"
           printf "%s\t%s\t%s\t1\t1\t1\t1\t1\t1\t1\t1\t%s\n", c, "CLsat3", "x", "1_0.05"
           printf "%s\tshort\n", c
           printf "%s\t%s\t%s\t1\t1\t1\t1\t1\t1\t1\t1\t%s\n", c, "CLsat2", "97", "abc" }' "$q"
EOF
chmod +x "$P/bin/fake_ssearch36"
cfgw fakess "tools: {ssearch36: $P/bin/fake_ssearch36}"
init "$W/base15"; step fake_ssearch -- run "$C/fakess.yaml" --only 09_assign --force --skip-tool-check

# ------------------------------------------------------------------ 4. manual loci (real 15, from base15)
SEG=$W/base15/arbitration/MASTER.loci.segmented.bed
# scaffold -> genome code from the .fai files (00_tables indexed the genomes)
{ "$GAWK" -F '\t' -v OFS='\t' '{ print $1, "spa" }' "$P/toy/genomes/spa.fna.fai"
  "$GAWK" -F '\t' -v OFS='\t' '{ print $1, "spb" }' "$P/toy/genomes/spb.fna.fai"; } > "$W/scf_code.tsv"
# loci of the pipeline catalogue re-entered by hand: inferred family (2), family in col 7 (1); track/browser/comment lines skipped
"$GAWK" -F '\t' -v OFS='\t' 'NR == FNR { sp[$1] = $2; next }
    FNR == 1 { print "track name=manual"; print "browser position x"; print "# hand-added" }
    FNR <= 2 { print $1, $2, $3, sp[$1] "|" $1 ":" $2 "-" $3; next }
    FNR == 3 { print $1, $2 + 100, $3 - 100, sp[$1] "|manual3", "0", "+", "CLsat2"; exit }' "$W/scf_code.tsv" "$SEG" > "$P/ref/manual_ok.bed"
cfgw manok "manual_loci: $P/ref/manual_ok.bed"
init "$W/base15"; step manual_ok -- run "$C/manok.yaml" --only 15_catalogue
save base15m
step manual_rerun -- run "$C/manok.yaml" --only 15_catalogue
# a 1-bp locus in background sequence (>= 20 kb from every MASTER locus; no scan window starts inside it): scores 0 with
# every family -> StageError
# after the earlier rows were appended
"$GAWK" -F '\t' -v OFS='\t' 'FILENAME ~ /MASTER/ { n[$1]++; s[$1, n[$1]] = $2; e[$1, n[$1]] = $3; next }
    { for (a = 10000; a + 2000 < $2; a += 10000) { ok = 1
        for (i = 1; i <= n[$1]; i++) if (a + 2000 > s[$1, i] - 20000 && a < e[$1, i] + 20000) ok = 0
        if (ok) { print $1, a + 3, a + 4, "spa|background"; exit } } }' "$W/base15/master/MASTER.merged.bed" "$P/toy/genomes/spa.fna.fai" > "$W/bg.bed"
echo "background locus: $(cat "$W/bg.bed")"
cat "$P/ref/manual_ok.bed" "$W/bg.bed" > "$P/ref/manual_zero.bed"
cfgw manzero "manual_loci: $P/ref/manual_zero.bed"
init "$W/base15"; step manual_zero -- run "$C/manzero.yaml" --only 15_catalogue
# an unknown family in col 7 after two good rows
{ sed -n 4,5p "$P/ref/manual_ok.bed"; "$GAWK" -F '\t' -v OFS='\t' 'NR == 6 { print $1, $2, $3, $4, 0, "+", "CLsat9" }' "$P/ref/manual_ok.bed"; } > "$P/ref/manual_badfam.bed"
cfgw manbadfam "manual_loci: $P/ref/manual_badfam.bed"
init "$W/base15"; step manual_badfam -- run "$C/manbadfam.yaml" --only 15_catalogue
# names without a configured genome (3 cols, no '|', unknown code x2): only the first three are listed
{ sed -n 4p "$P/ref/manual_ok.bed"; printf 'chrX\t1\t100\n'; printf "chrY\t1\t100\tnobar\nchrZ\t1\t100\tzz|q\nchrW\t1\t100\tww|q\n"; } > "$P/ref/manual_badgen.bed"
cfgw manbadgen "manual_loci: $P/ref/manual_badgen.bed"
init "$W/base15"; step manual_badgen -- run "$C/manbadgen.yaml" --only 15_catalogue
: > "$P/ref/manual_empty.bed"
cfgw manempty "manual_loci: $P/ref/manual_empty.bed"
init "$W/base15"; step manual_empty -- run "$C/manempty.yaml" --only 15_catalogue

# ------------------------------------------------------------------ 5. seeds given for every family (03-09 skipped; real 10-15 from base08)
for f in $FAMS; do cp "$W/base15/seeds/$f.bed" "$P/seeds/$f.bed"; done
cfgw seeds "inputs: {seeds: {CLsat1: $P/seeds/CLsat1.bed, CLsat2: $P/seeds/CLsat2.bed, CLsat3: $P/seeds/CLsat3.bed, CLsat4: $P/seeds/CLsat4.bed}}"
init "$W/base08"; step seeds_given -- run "$C/seeds.yaml" --from 10_seeds --to 15_catalogue
step seeds_rerun -- run "$C/seeds.yaml" --from 10_seeds --to 15_catalogue

echo "WP7 stages: $NPASS passed, $NFAIL failed  (work dir $W)"
[ $NFAIL -eq 0 ]
