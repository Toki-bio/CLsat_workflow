# stub_stages.sh - test-only stub stages for the bash driver (WP5 acceptance), the twin of tests/port/stub_run.py.
# Loaded through SATPIPE_STAGE_FILES; defines plan_<id>/st_<id> for every stage except 21_anchor (the real one).
# Written as WP6-8 stage code should be written: plan helpers, sh_run/sh_capture/sh_to, manifest_set, `if` not `&&`.

STUB_SCRIPT='o=$1; shift; cat "$@" > "$o"; printf "%s %s\n" "$0" "${STUB_TAG:-none}" >> "$o"; i=0; while [ $i -lt 20 ]; do i=$((i + 1)); echo "log line $i of $0"; done; [ "$0" != "${STUB_FAIL:-}" ]'

_stub_deps() {
    case $1 in
        02_trf|03_motifs) echo 00_tables ;; 02b_large) echo 02_trf ;; 04_filter) echo 03_motifs ;; 05_rotate) echo 04_filter ;;
        06_cluster) echo 05_rotate ;; 07_mafft) echo 06_cluster ;; 08_cross) echo 07_mafft ;; 09_assign) echo 08_cross ;;
        10_seeds) echo 09_assign ;; 11_expand) echo 10_seeds 02b_large ;; 12_master) echo 11_expand ;; 13_rescan) echo 12_master ;;
        14_arbitrate) echo 13_rescan ;; 15_catalogue) echo 14_arbitrate ;; 16_regions) echo 15_catalogue ;; 17_terminal) echo 16_regions ;;
        18_subfamily) echo 17_terminal ;; 19_orthology) echo 16_regions ;; 20_candidates) echo 15_catalogue 19_orthology ;;
    esac
}
# STUB_IN=() existing outputs of the stages SID depends on
_stub_in() { local d; STUB_IN=(); for d in $(_stub_deps "$1"); do if [ -f "$RUN/$d/stub.out" ]; then STUB_IN+=("$RUN/$d/stub.out"); fi; done; }
_stub_cons()    { local k; while IFS= read -r k; do STUB_X+=("$(cfg "families/$k/consensus")"); done < <(cfg_keys families); }
_stub_genomes() { local k; while IFS= read -r k; do STUB_X+=("$(cfg "genomes/$k")"); done < <(cfg_keys genomes); }

# stub_gen SID CWD(0|1) [sh_run options...]   uses STUB_X=() extra arguments
stub_gen() {
    local sid=$1 usecwd=$2 d o
    shift 2
    d=$(ctx_d "$sid"); o=$d/stub.out
    _stub_in "$sid"
    if [ "$usecwd" = 1 ]; then sh_run "$sid" "$d" "$@" -- bash -c "$STUB_SCRIPT" "$sid" stub.out "${STUB_IN[@]}" "${STUB_X[@]}"
    else sh_run "$sid" "" "$@" -- bash -c "$STUB_SCRIPT" "$sid" "$o" "${STUB_IN[@]}" "${STUB_X[@]}"; fi
    emit_outputs "$o" "$d"
}

_stub_plan_deps() { _stub_in "$1"; emit_inputs "${STUB_IN[@]}"; }

_stub_generic() { STUB_X=(); stub_gen "$1" 0; }

for _sid in "${STAGE_IDS[@]}"; do
    [ "$_sid" = 21_anchor ] && continue
    eval "plan_$_sid() { _stub_plan_deps $_sid; }"
    eval "st_$_sid() { _stub_generic $_sid; }"
done
unset _sid

plan_00_tables() {
    cons_files
    case $(cfg_type references/rotated_consensus 2>/dev/null) in path|str) emit_inputs "$(cfg references/rotated_consensus)" ;; esac
    cfg_subtree genomes genomes; cfg_subtree families families; param_cfg_or_null subfamily/members members
}
st_00_tables() { STUB_X=(); _stub_cons; stub_gen 00_tables 0; }

plan_02_trf() { _stub_plan_deps 02_trf; genome_files; cfg_subtree discovery/trf trf; cfg_subtree finalisation/large_scaffold_bp large; cfg_subtree genomes genomes; }
st_02_trf() { STUB_X=(); _stub_genomes; stub_gen 02_trf 0 log_name=02_trf_split; }

plan_03_motifs() { _stub_plan_deps 03_motifs; cfg_subtree discovery discovery; param_set step3b str runs-v3; }
st_03_motifs() {
    STUB_X=()
    stub_gen 03_motifs 1 PERIOD_MIN="$(cfg discovery/period/0)" PERIOD_MAX="$(cfg discovery/period/1)" MIN_COPIES="$(cfg discovery/min_copies)" \
        TRF_PARAMS="$(cfg discovery/trf)" SIMPLE_REPEAT_COV="$(cfg discovery/simple_repeat_cov)" STRAND="$(cfg discovery/strand)" \
        CLUSTER_ID="$(cfg discovery/cluster_id)" CROSS_CLUSTER_ID="$(cfg discovery/cross_cluster_id)"
}

st_07_mafft() { STUB_X=(); stub_gen 07_mafft 1 log_name=07_mafft_a; : > "$STAGE_OUTPUTS"; stub_gen 07_mafft 1 log_name=07_mafft_b STUB_TAG=b; }

plan_09_assign() { emit_inputs "$(cfg family_assignment)"; _stub_plan_deps 09_assign; cons_files; }
st_09_assign() {
    local d o out n m
    d=$(ctx_d 09_assign); o=$d/stub.out
    _stub_in 09_assign
    sh_capture out 09_assign -- "$(tool ssearch36)" -n -z 3 "${STUB_IN[@]}"
    printf '%s' "$out" > "$o"
    n=$("$GAWK" 'END { print NR }' "$o")
    m=$(md5sum < "$o"); m=${m%% *}
    manifest_set centroids_n "$n"
    manifest_set centroids_md5 "\"$m\""
    say "  centroids: $n"
    emit_outputs "$o"
}

st_10_seeds() {
    local d o f i=0 j="{" sep=""
    d=$(ctx_d 10_seeds); o=$d/stub.out
    _stub_in 10_seeds
    sh_to "$o" 10_seeds -- bash -c 'echo seeds; cat "$@"' mk "${STUB_IN[@]}"
    while IFS= read -r f; do i=$((i + 1)); j+="$sep\"$f\": $i"; sep=", "; done < <(cfg_keys families)
    manifest_set seed_counts "$j}"
    emit_outputs "$o"
}

plan_11_expand() {
    _stub_plan_deps 11_expand; cons_files; genome_files
    cfg_subtree expansion expansion; cfg_subtree scan scan
    local f; param_set mono map "$(cfg_keys families | "$GAWK" 'END { print NR }')"
    while IFS= read -r f; do cfg_subtree "families/$f/mono" "mono/$f"; done < <(cfg_keys families)
}
st_11_expand() {
    local f d
    d=$(ctx_d 11_expand)
    while IFS= read -r f; do
        local -a x=(--min-bits "$(cfg expansion/min_bits)" --scan-win-bp "$(( $(cfg scan/win_mono) * $(cfg "families/$f/mono") ))"
                    --jobs "$(cfg expansion/jobs)" --cluster-id "$(cfg discovery/cluster_id)")
        if [ "$(cfg expansion/to_scaffold_ends)" = True ]; then x+=(--to-scaffold-ends); fi
        sh_run 11_expand "" log_name="11_expand_$f" -- bash -c 'printf "%s\n" "$@" >> "$0"' "$d/argv.txt" "${x[@]}"
    done < <(cfg_keys families)
    STUB_X=("$d/argv.txt"); stub_gen 11_expand 0
}

st_12_master() { manifest_set master_loci 627; _stub_generic 12_master; }

plan_14_arbitrate() { _stub_plan_deps 14_arbitrate; cfg_subtree arbitration arbitration; cfg_subtree expansion/min_bits min_bits; param_set rule str bin-majority-v1; }
st_14_arbitrate() {
    local b
    if [ "$(cfg_type arbitration/bin_min_bits)" = null ]; then b=$(cfg expansion/min_bits); else b=$(cfg arbitration/bin_min_bits); fi
    STUB_X=(); stub_gen 14_arbitrate 0 BIN_MIN_BITS="$b" DELTA="$(cfg arbitration/delta)" STUB_TAG=arb
}

st_15_catalogue() {
    manifest_set segmented_loci 650
    manifest_set catalogue_loci 651
    manifest_set manual_loci '[{"locus": "L1", "bed": ["spa_scf1", "0", "10", "x"], "family": "CLsat1", "bits": {"CLsat1": 12.5, "CLsat2": 0.0}}]'
    STUB_X=(); stub_gen 15_catalogue 1
}

plan_16_regions() {
    _stub_plan_deps 16_regions; cfg_subtree finalisation fin
    if [ "$(cfg_type hapdup/pairs 2>/dev/null)" = list ] && [ "$(cfg hapdup/pairs)" != "[]" ]; then cfg_subtree hapdup hapdup; fi
}

plan_18_subfamily() { _stub_plan_deps 18_subfamily; subx_inputs; param_cfg_or_null subfamily subfamily; param_cfg_or_null anchor anchor; }

plan_19_orthology() { _stub_plan_deps 19_orthology; cfg_subtree orthology orthology; }
st_19_orthology() { manifest_set ortholog_pairs 53; _stub_generic 19_orthology; }

plan_20_candidates() { _stub_plan_deps 20_candidates; cfg_subtree expansion/min_bits min_bits; cfg_subtree orthology/flank_bp flank; }
