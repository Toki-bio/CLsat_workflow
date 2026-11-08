#!/usr/bin/env bash
# tests/test_schema.py (12 tests) for the bash driver: values via `satpipe config --flat`, errors via `satpipe config`
# (stderr `error: ...`), warnings via `satpipe check --no-files --no-tools` ([WARN] config lines).
#   bash tests/sh/test_schema.sh       (needs gawk on PATH or SATPIPE_GAWK; exit 1 on the first failure)
set -euo pipefail
export LC_ALL=C
ROOT=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../.." && pwd -P)
SP=$ROOT/bin/satpipe
# shellcheck source=fixture.sh
. "$ROOT/tests/sh/fixture.sh"
TMP=$(mktemp -d "${TMPDIR:-/tmp}/satpipe_tsch.XXXXXX"); TMP=$(cd "$TMP" && pwd -P)
trap 'rm -rf -- "$TMP"' EXIT
N=0
check() { local label=$1; shift; N=$((N + 1)); if "$@"; then printf 'ok   %2d %s\n' "$N" "$label"; else printf 'FAIL %2d %s\n' "$N" "$label"; exit 1; fi; }

new_project() { D=$TMP/t$N; rm -rf -- "$D"; make_project "$D"; CFG=$D/configs/run.yaml; }
# load CFG -> $FLAT (TYPED flat) or $ERR (error lines), $RC
load() { RC=0; FLAT=$("$SP" config --flat "$1" 2> "$TMP/err") || RC=$?; ERR=$(cat "$TMP/err"); }
# row PATH -> "type<TAB>text" of the TYPED flat row (status 1 when absent)
row() { local p=$1 l; while IFS= read -r l; do [ "${l%%$'\t'*}" = "$p" ] && { l=${l#*$'\t'}; printf '%s' "${l%$'\t'*}"; return 0; }; done <<< "$FLAT"; return 1; }
val() { local r; r=$(row "$1") || return 1; printf '%s' "${r#*$'\t'}"; }
typ() { local r; r=$(row "$1") || return 1; printf '%s' "${r%%$'\t'*}"; }
# numeq PATH NUMBER   Python ==: an int or float leaf numerically equal
numeq() { case $(typ "$1") in int|float) ;; *) return 1 ;; esac; gawk -v a="$(val "$1")" -v b="$2" 'BEGIN { exit !(a + 0 == b + 0) }'; }
streq() { [ "$(typ "$1")" = str ] && [ "$(val "$1")" = "$2" ]; }
isnull() { [ "$(typ "$1")" = null ]; }
# keys PATH  child keys, space-joined, document order
keys() { local p=$1/ l k out=""; while IFS= read -r l; do k=${l%%$'\t'*}; [[ $k == "$p"* ]] || continue; k=${k#"$p"}; [[ $k == */* ]] || out+=${out:+ }$k; done <<< "$FLAT"; printf '%s' "$out"; }
haserr() { [[ $ERR == *"$1"* ]]; }

t_defaults_are_published_values() {
    new_project; load "$CFG"; [ $RC = 0 ] || return 1
    streq discovery/trf "2 5 7 80 10 50 2000" || return 1
    [ "$(typ discovery/period)" = list ] && [ "$(val discovery/period)" = 2 ] && [ "$(typ discovery/period/0)" = int ] \
        && [ "$(val discovery/period/0)" = 80 ] && [ "$(val discovery/period/1)" = 300 ] || return 1
    streq discovery/strand plus && numeq discovery/cluster_id 0.90 && numeq discovery/cross_cluster_id 0.90 || return 1
    numeq expansion/min_bits 50 && numeq expansion/min_aln_frac 0.67 && numeq expansion/max_steps 500 && numeq expansion/probe_mono 100 || return 1
    [ "$(keys arbitration)" = "delta delta_ambig seg_bin seg_min_run seg_min_locus short_segments bin_min_bits max_gap_monomers bin_share" ] || return 1
    numeq arbitration/delta 20.0 && numeq arbitration/delta_ambig 15.0 && numeq arbitration/seg_bin 1000 && numeq arbitration/seg_min_run 10 \
        && numeq arbitration/seg_min_locus 2000 && streq arbitration/short_segments drop && isnull arbitration/bin_min_bits \
        && numeq arbitration/max_gap_monomers 100 && numeq arbitration/bin_share 0.75 || return 1
    [ "$(keys orthology)" = "flank_bp guard_bp min_aln min_pid" ] || return 1
    numeq orthology/flank_bp 600 && numeq orthology/guard_bp 0 && numeq orthology/min_aln 200 && numeq orthology/min_pid 80.0 || return 1
    numeq subfamily/confident_delta 10.0 || return 1
    [ "$(keys scan)" = "win_mono step_mono min_bits" ] && numeq scan/win_mono 2 && numeq scan/step_mono 1 && numeq scan/min_bits 0.0
}
check test_defaults_are_published_values t_defaults_are_published_values

t_paths_resolve_relative_to_config() {
    new_project; load "$CFG"
    [ "$(val genomes/spa)" = "$D/genomes/spa.fna" ] && [ "$(val out_dir)" = "$D/runs/toy" ] && [[ $(val families/SatX/consensus) == *ref/SatX.cons.fa ]]
}
check test_paths_resolve_relative_to_config t_paths_resolve_relative_to_config

t_unknown_key_is_error() {
    new_project; append "$CFG" $'\nexpansion:\n  min_bitz: 3\n'; load "$CFG"
    [ $RC = 1 ] && haserr "expansion.min_bitz: unknown key"
}
check test_unknown_key_is_error t_unknown_key_is_error

t_missing_required_and_bad_types_all_reported() {
    # validate() on this dict (base_dir None): the same content as a YAML file
    mkdir -p "$TMP/v"
    printf 'run_name: x y\ngenomes: {}\nfamilies:\n  F: {mono: "146"}\ndiscovery: {strand: minus, cluster_id: 0}\n' > "$TMP/v/c.yaml"
    load "$TMP/v/c.yaml"; [ $RC = 1 ] || return 1
    local s
    for s in "out_dir: required" "run_name: 'x y' does not match" "genomes: expected at least 1 entries" "families.F.consensus: required" \
             "families.F.mono: expected an integer" "discovery.strand: must be one of" "discovery.cluster_id: must be > 0.0" "family_assignment: required"; do
        haserr "$s" || { echo "missing: $s"; return 1; }
    done
}
check test_missing_required_and_bad_types_all_reported t_missing_required_and_bad_types_all_reported

t_cross_checks() {
    new_project; append "$CFG" $'\nsubfamily:\n  members:\n    SatY: [a]\n  consensuses: ../ref\nviewer:\n  species_names:\n    zzz: Z\n'; load "$CFG"
    [ $RC = 1 ] && haserr "subfamily.members.SatY: not a family" && haserr "viewer.species_names.zzz: not a code"
}
check test_cross_checks t_cross_checks

t_trf_dat_code_must_be_a_genome() {
    new_project; edit "$CFG" "spb: ../trf/spb.dat" "spc: ../trf/spb.dat"; load "$CFG"
    [ $RC = 1 ] && haserr "inputs.trf_dat.spc: not a code" || return 1
    edit "$CFG" $'    spc: ../trf/spb.dat\n' ""
    load "$CFG"; [ $RC = 0 ] || return 1                     # partial coverage is a warning, not an error
    local o; o=$("$SP" check "$CFG" --no-files --no-tools) || return 1
    [[ $o == *"[WARN] config: "*"covers only part"* ]]
}
check test_trf_dat_code_must_be_a_genome t_trf_dat_code_must_be_a_genome

t_schema_doc_lists_every_leaf() {
    local md k; md=$("$SP" config --schema-doc) || return 1
    for k in '`run_name`' '`genomes`' '`families.<key>.mono`' '`expansion.min_aln_frac`' '`arbitration.short_segments`' \
             '`references.rotated_consensus`' '`tools`'; do
        [[ $md == *"$k"* ]] || { echo "missing $k"; return 1; }
    done
}
check test_schema_doc_lists_every_leaf t_schema_doc_lists_every_leaf

t_darevskia_config_is_schema_valid() {
    local c=$ROOT/configs/darevskia_clsat.yaml f o
    load "$c"; [ $RC = 0 ] || return 1
    [ "$(keys genomes | tr ' ' '\n' | LC_ALL=C sort | tr '\n' ' ')" = "arm mix nai unm unp val " ] || return 1
    [ "$(val families/CLsat1/mono) $(val families/CLsat2/mono) $(val families/CLsat3/mono) $(val families/CLsat4/mono)" = "146 145 146 116" ] || return 1
    streq arbitration/short_segments drop || return 1
    [ "$(val tools/bin_dir)" = /staging/conda/envs/bioinfo/bin ] || return 1
    [[ $(val family_assignment) == */satpipe/stages/02_seed_expansion_scan/cluster_family.tsv ]] && [ -f "$(val family_assignment)" ] || return 1
    for f in $(keys families); do [ -f "$(val "families/$f/consensus")" ] || return 1; done
    [ -f "$(val manual_loci)" ] || return 1
    o=$("$SP" check "$c" --no-files --no-tools) || return 1
    [[ $o != *"[WARN] config: "* ]]                          # not c.warnings
}
check test_darevskia_config_is_schema_valid t_darevskia_config_is_schema_valid

t_example_minimal_is_schema_valid() {
    load "$ROOT/configs/example_minimal.yaml"; [ $RC = 0 ] && numeq expansion/jobs 16
}
check test_example_minimal_is_schema_valid t_example_minimal_is_schema_valid

t_hapdup_off_by_default() {
    new_project; load "$CFG"; [ $RC = 0 ] || return 1
    local t; t=$(typ hapdup/pairs || true)                    # not (x.get("hapdup") or {}).get("pairs")
    [ -z "$t" ] || [ "$t" = null ] || { [ "$t" = list ] && [ "$(val hapdup/pairs)" = 0 ]; }
}
check test_hapdup_off_by_default t_hapdup_off_by_default

t_hapdup_pairs_validate() {
    new_project; append "$CFG" $'\nhapdup:\n  pairs: [[spa, spa]]\n'; load "$CFG"; [ $RC = 0 ] || return 1
    [ "$(val hapdup/pairs)" = 1 ] && [ "$(val hapdup/pairs/0)" = 2 ] && streq hapdup/pairs/0/0 spa && streq hapdup/pairs/0/1 spa || return 1
    numeq hapdup/min_ident 0.99 && numeq hapdup/min_qcov 0.95 && numeq hapdup/edge_tol 100 && numeq hapdup/min_tovl 0.5
}
check test_hapdup_pairs_validate t_hapdup_pairs_validate

t_hapdup_pair_needs_two_codes() {
    new_project; append "$CFG" $'\nhapdup:\n  pairs: [[spa]]\n'; load "$CFG"
    [ $RC = 1 ] && haserr "hapdup.pairs[0]"
}
check test_hapdup_pair_needs_two_codes t_hapdup_pair_needs_two_codes

echo "test_schema.sh: $N/12 passed"
