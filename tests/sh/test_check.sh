#!/usr/bin/env bash
# tests/test_check.py (12 tests) for the bash driver: bin/satpipe check / run --list-stages / stages / config, tools.sh.
# The `project` fixture is tests/sh/fixture.sh (conftest.py:24-48); conftest's MAX_BASE_PATH lift is SATPIPE_MAX_BASE_PATH.
#   bash tests/sh/test_check.sh        (needs gawk on PATH or SATPIPE_GAWK; exit 1 on the first failure)
set -euo pipefail
export LC_ALL=C
ROOT=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../.." && pwd -P)
SP=$ROOT/bin/satpipe
# shellcheck source=fixture.sh
. "$ROOT/tests/sh/fixture.sh"
TMP=$(mktemp -d "${TMPDIR:-/tmp}/satpipe_tchk.XXXXXX"); TMP=$(cd "$TMP" && pwd -P)
trap 'rm -rf -- "$TMP"' EXIT
N=0
check() { local label=$1; shift; N=$((N + 1)); if "$@"; then printf 'ok   %2d %s\n' "$N" "$label"; else printf 'FAIL %2d %s\n' "$N" "$label"; exit 1; fi; }

# fresh project per test (pytest's tmp_path)
new_project() { D=$TMP/t$N; rm -rf -- "$D"; make_project "$D"; CFG=$D/configs/run.yaml; }
# report CFG ARGS...  -> $OUT (text report), $RC
report() { RC=0; OUT=$("$SP" check "$@") || RC=$?; }
# has LEVEL SECTION SUBSTRING    an item of that level (any level when LEVEL is -) and section whose message holds SUBSTRING
has() {
    local l
    while IFS= read -r l; do
        [[ $l == "["* ]] || continue
        local lv=${l:1:4} rest=${l:7}; lv=${lv%% *}
        local sec=${rest%%: *} msg=${rest#*: }
        if { [ "$1" = - ] || [ "$lv" = "$1" ]; } && { [ "$2" = - ] || [ "$sec" = "$2" ]; } && [[ $msg == *"$3"* ]]; then return 0; fi
    done <<< "$OUT"
    return 1
}
failed() { [ "$RC" = 1 ] && [[ $OUT == *"check FAILED:"* ]]; }
passed() { [ "$RC" = 0 ] && [[ $OUT == *"check passed:"* ]]; }

t_fixture_passes_without_tools() {
    new_project; report "$CFG" --no-tools
    passed && has - - "148 bp = mono" && has - - "2 cluster->family rows" && has - - "manual_loci: 1 intervals" \
        && has - - "stage script references present"
}
check test_fixture_passes_without_tools t_fixture_passes_without_tools

t_mono_mismatch_is_warning() {
    new_project; edit "$CFG" "mono: 148" "mono: 146"; report "$CFG" --no-tools
    passed && has WARN families "consensus length 148 != mono 146"
}
check test_mono_mismatch_is_warning t_mono_mismatch_is_warning

t_unknown_family_in_table_fails() {
    new_project; printf '0\tSatX\n7\tSatQ\n' > "$D/ref/cluster_family.tsv"; report "$CFG" --no-tools
    failed && has - family_assignment "family 'SatQ' is not in"
}
check test_unknown_family_in_table_fails t_unknown_family_in_table_fails

t_family_without_cluster_fails() {
    new_project; edit "$CFG" $'families:\n' $'families:\n  SatY: {consensus: ../ref/SatX.cons.fa, mono: 148}\n'; report "$CFG" --no-tools
    has FAIL family_assignment "family SatY has no cluster"
}
check test_family_without_cluster_fails t_family_without_cluster_fails

t_missing_genome_fails() {
    new_project; rm -f -- "$D/genomes/spb.fna"; report "$CFG" --no-tools
    failed || return 1
    local l; while IFS= read -r l; do [[ $l == "[FAIL] genomes: "*spb*"not found"* ]] && return 0; done <<< "$OUT"; return 1
}
check test_missing_genome_fails t_missing_genome_fails

t_bad_bed_is_reported() {
    new_project; printf 'scf\t10\t5\tspa|x\n' > "$D/ref/manual.bed"; report "$CFG" --no-tools
    has FAIL manual_loci "0 <= start < end"
}
check test_bad_bed_is_reported t_bad_bed_is_reported

t_tools_missing_are_reported() {
    new_project; mkdir -p "$TMP/empty"
    local bash_bin gawk_bin; bash_bin=$(command -v bash); gawk_bin=$(command -v "${SATPIPE_GAWK:-gawk}")
    RC=0; OUT=$(PATH=$TMP/empty SATPIPE_GAWK=$gawk_bin "$bash_bin" "$SP" check "$CFG" --no-files) || RC=$?
    failed || return 1
    local names=" " fails=" " l name t
    while IFS= read -r l; do
        [[ $l == "["*"] tools: "* ]] || continue
        name=${l#*] tools: }; name=${name%%:*}; names+="$name "
        [[ $l == "[FAIL]"* ]] && fails+="$name "
    done <<< "$OUT"
    for t in trf vsearch ssearch36 gawk samtools bedtools parallel mafft; do [[ $names == *" $t "* ]] || return 1; done
    while IFS=$'\037' read -r t _ _ _ req _; do                 # every required tool of tools.tsv failed
        case $t in ''|'#'*) continue ;; esac
        [ "$req" = True ] || continue
        [[ $fails == *" $t "* ]] || { echo "required $t did not fail"; return 1; }
    done < <(gawk '{ gsub(/\t/, "\037"); print }' "$ROOT/lib/satpipe/tools.tsv")
}
check test_tools_missing_are_reported t_tools_missing_are_reported

t_tool_override_and_bin_dir() {
    new_project; local fake=$TMP/bin$N; mkdir -p "$fake"
    printf '#!/bin/sh\necho '"'"'vsearch v2.30.0_linux_x86_64'"'"' >&2\n' > "$fake/vsearch"; chmod 755 "$fake/vsearch"
    append "$CFG" $'\ntools:\n  bin_dir: '"$fake"$'\n'
    report "$CFG" --no-files
    local l bin_ok=0 vs=""
    while IFS= read -r l; do
        [[ $l == "[OK  ] tools: bin_dir:"* ]] && bin_ok=1
        [[ $l == "["*"] tools: vsearch"[\ :]* ]] && vs=$l
    done <<< "$OUT"
    [ $bin_ok = 1 ] && [ -n "$vs" ] && [[ $vs != *"not found"* ]] || return 1
    local after=${vs#*: }; after=${after#*: }                   # msgs["vsearch"][1].split(": ")[1]
    [[ $after == "$fake"* ]] || [[ $vs == *v2.30.0* ]]
}
check test_tool_override_and_bin_dir t_tool_override_and_bin_dir

t_cli_exit_codes() {
    new_project; local o
    "$SP" check "$CFG" --no-tools > /dev/null || return 1
    o=$("$SP" check "$CFG" --no-tools --json) || return 1
    [[ $o == *'"failed": false'* ]] || return 1
    printf 'x\tSatX\n' > "$D/ref/cluster_family.tsv"
    local rc=0; "$SP" check "$CFG" --no-tools > /dev/null || rc=$?; [ $rc = 1 ] || return 1
    o=$("$SP" run --list-stages "$CFG") || return 1; [[ $o == *02_trf* ]] || return 1
    "$SP" stages > /dev/null || return 1
    "$SP" config --schema-doc > /dev/null || return 1
    "$SP" config "$CFG" > /dev/null || return 1
    rc=0; "$SP" config > /dev/null 2>&1 || rc=$?; [ $rc = 2 ]
}
check test_cli_exit_codes t_cli_exit_codes

t_probe_all_runs_without_crashing() {
    local n_tsv n_probe
    n_tsv=$(gawk -F '\t' '!/^#/ && $1 != "" { n++ } END { print n + 0 }' "$ROOT/lib/satpipe/tools.tsv")
    n_probe=$(bash -c '. "$1/lib/satpipe/common.sh"; . "$1/lib/satpipe/tools.sh"; tools_probe_all' _ "$ROOT" | wc -l)
    [ "$n_tsv" = 14 ] && [ "$n_probe" = "$n_tsv" ]
}
check test_probe_all_runs_without_crashing t_probe_all_runs_without_crashing

t_long_out_dir_fails() {
    new_project; local rc=0
    SATPIPE_MAX_BASE_PATH=20 "$SP" check "$CFG" --no-tools > /dev/null || rc=$?
    [ $rc = 1 ]
}
check test_long_out_dir_fails t_long_out_dir_fails

t_unknown_anchor_family_and_tool_key_are_config_errors() {
    new_project; local base e rc
    base=$(cat "$CFG"; printf x); base=${base%x}
    printf '%s\nanchor: {family: NoSuchFamily}\n' "$base" > "$CFG"
    rc=0; e=$("$SP" config "$CFG" 2>&1 >/dev/null) || rc=$?
    [ $rc = 1 ] && [[ $e == *anchor.family* ]] || return 1
    printf '%s\ntools: {ssearch_36: /opt/x/ssearch36}\n' "$base" > "$CFG"
    rc=0; e=$("$SP" config "$CFG" 2>&1 >/dev/null) || rc=$?
    [ $rc = 1 ] && [[ $e == *tools.ssearch_36* ]]
}
check test_unknown_anchor_family_and_tool_key_are_config_errors t_unknown_anchor_family_and_tool_key_are_config_errors

# --- stage scope of the check in `satpipe run` (check.sh _ck_scope, run.sh _run_scope; 2026-11-26)
# runplan ARGS...  -> $OUT (stderr of `satpipe run --plan`: the check report or its scoped WARNs), $RC; nothing is written
runplan() { RC=0; OUT=$("$SP" run "$@" --plan 2>&1 >/dev/null) || RC=$?; }
# warned SECTION SUBSTRING   a FAIL turned WARN by the scope
warned() { local l; while IFS= read -r l; do [[ $l == "satpipe run: [WARN] $1: "*"$2"*"(not a failure for this run: "* ]] && return 0; done <<< "$OUT"; return 1; }
failline() { local l; while IFS= read -r l; do [[ $l == "[FAIL] $1: "*"$2"* ]] && return 0; done <<< "$OUT"; return 1; }

t_early_range_family_table_is_warning() {   # (a) de novo discovery before the curated table exists
    new_project; : > "$D/ref/cluster_family.tsv"
    runplan "$CFG" --to 04_filter --skip-tool-check
    [ "$RC" = 0 ] && warned family_assignment "family SatX has no cluster in the table" && [[ $OUT == *"only stage 09_assign needs it"* ]] || return 1
    rm -f -- "$D/ref/cluster_family.tsv"                        # table not written yet at all
    runplan "$CFG" --to 08_cross --skip-tool-check
    [ "$RC" = 0 ] && warned family_assignment "table: file not found" || return 1
    runplan "$CFG" --only 02_trf,03_motifs --skip-tool-check
    [ "$RC" = 0 ] && warned family_assignment "table: file not found"
}
check test_early_range_family_table_is_warning t_early_range_family_table_is_warning

t_range_with_stage_09_still_fails() {        # (b) the same table fails whenever stage 09 runs, and in `satpipe check`
    new_project; : > "$D/ref/cluster_family.tsv"
    local a
    for a in "" "--to 09_assign" "--from 08_cross --to 10_seeds" "--only 04_filter,09_assign" "--from 00_tables"; do
        # shellcheck disable=SC2086
        runplan "$CFG" $a --skip-tool-check
        [ "$RC" = 1 ] && failline family_assignment "family SatX has no cluster in the table" \
            && [[ $OUT == *"refusing to start"* ]] && [[ $OUT != *"not a failure for this run"* ]] || { echo "range '$a'"; return 1; }
    done
    report "$CFG" --no-tools
    failed && has FAIL family_assignment "family SatX has no cluster in the table" && ! has - stages "check scoped"
}
check test_range_with_stage_09_still_fails t_range_with_stage_09_still_fails

t_later_stage_tool_missing_is_warning() {    # (c) a tool only later stages need: WARN for an early range, FAIL once a stage needs it
    new_project; local fb=$TMP/fakebin$N t
    mkdir -p "$fb"
    ln -s "$(command -v bash)" "$fb/bash"; ln -s "$(command -v "${SATPIPE_GAWK:-gawk}")" "$fb/gawk"
    _fake() { printf '#!/bin/sh\necho "%s"\n' "$2" > "$fb/$1"; chmod 755 "$fb/$1"; }
    _fake trf 'Tandem Repeats Finder, Version 4.10.0'; _fake parallel 'GNU parallel 20240122'
    _fake samtools 'samtools 1.22.1'; _fake python3 'Python 3.12.11'; _fake ssearch36 ' version: 36.3.8h'
    # missing on purpose: vsearch (06, 08), mafft (07), bedtools (11+)
    RC=0; OUT=$(PATH=$fb SATPIPE_GAWK=$fb/gawk "$fb/bash" "$SP" run "$CFG" --to 04_filter --plan 2>&1 >/dev/null) || RC=$?
    [ "$RC" = 0 ] || { printf '%s\n' "$OUT" | tail -n 20; return 1; }
    for t in vsearch mafft bedtools; do warned tools "$t: not found" || { echo "no WARN for $t"; return 1; }; done
    RC=0; OUT=$(PATH=$fb SATPIPE_GAWK=$fb/gawk "$fb/bash" "$SP" run "$CFG" --to 08_cross --plan 2>&1 >/dev/null) || RC=$?
    [ "$RC" = 1 ] && failline tools "vsearch: not found" && failline tools "mafft: not found" || return 1
    [[ $OUT == *"[WARN] tools: bedtools: not found"*"(not a failure for this run: only stages "*"12_master"* ]] || return 1
    rm -f -- "$fb/trf"                                          # a tool the early range itself needs still fails
    RC=0; OUT=$(PATH=$fb SATPIPE_GAWK=$fb/gawk "$fb/bash" "$SP" run "$CFG" --to 04_filter --plan 2>&1 >/dev/null) || RC=$?
    [ "$RC" = 1 ] && failline tools "trf: not found"
}
check test_later_stage_tool_missing_is_warning t_later_stage_tool_missing_is_warning

echo "test_check.sh: $N/15 passed"
