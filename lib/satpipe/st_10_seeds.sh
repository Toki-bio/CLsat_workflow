# st_10_seeds.sh - stage 10_seeds (run.py:720-748): seed BED per family (given, or from the Step-9 member headers)
# plus the large-scaffold recovered seeds. Stage contract: lib/satpipe/stages.sh header.

. "$SATPIPE_LIB/pyio.sh"

# _seeds_given   SEEDS_GIVEN[fam]=path for (inputs or {}).get('seeds') or {} (run.py:723); SEEDS_GIVEN_KEYS in config order
_seeds_given() {
    local k
    SEEDS_GIVEN=(); SEEDS_GIVEN_KEYS=()
    if [ "$(cfg_type inputs/seeds 2>/dev/null)" = map ]; then
        while IFS= read -r k; do SEEDS_GIVEN[$k]=$(cfg "inputs/seeds/$k"); SEEDS_GIVEN_KEYS+=("$k"); done < <(cfg_keys inputs/seeds)
    fi
    return 0
}

# run.py:745-748: members.fa of the families without a given seed BED, then the given BEDs (map order), then the
# recovered BEDs that exist; params {}
plan_10_seeds() {
    local -A SEEDS_GIVEN=(); local -a SEEDS_GIVEN_KEYS=(); local f
    _seeds_given
    while IFS= read -r f; do
        if [ -z "${SEEDS_GIVEN[$f]+x}" ]; then emit_inputs "$RUN/seeds/$f.members.fa"; fi
    done < <(cfg_keys families)
    for f in "${SEEDS_GIVEN_KEYS[@]}"; do emit_inputs "${SEEDS_GIVEN[$f]}"; done
    while IFS= read -r f; do
        if [ -f "$RUN/large/$f.recovered.bed" ]; then emit_inputs "$RUN/large/$f.recovered.bed"; fi
    done < <(cfg_keys families)
}

st_10_seeds() {
    local -A SEEDS_GIVEN=(); local -a SEEDS_GIVEN_KEYS=(); local -a fams=()
    local s f bed rec n j="{" msg="" sep=""
    _seeds_given
    s=$(ctx_d seeds)
    while IFS= read -r f; do fams+=("$f"); done < <(cfg_keys families)
    for f in "${fams[@]}"; do
        bed=$s/$f.bed
        if [ -n "${SEEDS_GIVEN[$f]+x}" ]; then
            py_copyfile "${SEEDS_GIVEN[$f]}" "$bed"
        else
            sh_to "$bed" 10_seeds -- bash "$SCRIPTS/02_seed_expansion_scan/make_bed_from_clsats_headers.sh" "$s/$f.members.fa"
        fi
        rec=$RUN/large/$f.recovered.bed
        if [ -f "$rec" ] && [ -s "$rec" ]; then cat -- "$rec" >> "$bed"; fi
        n=$(py_nlines "$bed")
        j+="$sep\"$f\": $n"; msg+="$sep$f $n"; sep=", "
        emit_outputs "$bed"
    done
    manifest_set seed_counts "$j}"
    say "  seeds: $msg"
}
