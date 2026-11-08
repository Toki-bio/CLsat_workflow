# st_19_orthology.sh - stage 19_orthology (run.py st_orthology, branch flank-single-copy). Contract: lib/satpipe/stages.sh header (WP5).
# The argv is built in an array in Python's order; option groups follow the four base options.

# shellcheck source=subx.sh
. "$SATPIPE_LIB/subx.sh"

# _orth_trf   (run.py _orth_trf) genome-wide TRF calls for the tandem mask, one "code<TAB>path" line each:
# trf/<code>.dat if it is a file, then sorted(glob(large/<code>/*/scaf.trf.dedup.dat)) (C-locale byte order)
_orth_trf() {
    local c q
    local -a codes=()
    mapfile -t codes < <(cfg_keys genomes)
    for c in "${codes[@]}"; do
        [ -f "$RUN/trf/$c.dat" ] && printf '%s\t%s\n' "$c" "$RUN/trf/$c.dat"
        for q in "$RUN/large/$c"/*/scaf.trf.dedup.dat; do
            [ -e "$q" ] && printf '%s\t%s\n' "$c" "$q"
        done
    done
    return 0
}

_orth_float_gt0() { V=$1 "$GAWK" 'BEGIN { exit !(ENVIRON["V"] + 0 > 0) }'; }

plan_19_orthology() {
    local -a trf=()
    emit_inputs "$RUN/regions/locus_master.tsv" "$RUN/tables/genomes.tsv" "$RUN/tables/families.tsv"
    cons_files
    genome_files
    if [ "$(cfg orthology/tandem_mask)" = True ]; then
        mapfile -t trf < <(_orth_trf | cut -f2-)
        [ ${#trf[@]} -eq 0 ] || emit_inputs "${trf[@]}"
    fi
    cfg_subtree orthology orthology                                    # {'orthology': cfg['orthology']}
}

st_19_orthology() {
    local out pairs c p lost="" steps
    local -a codes=() trf=()
    out=$(ctx_d ortho)
    local -a argv=(python3 "$SCRIPTS/04_flank_orthology/rebuild_ortho_segmented.py"
        --loci "$RUN/regions/locus_master.tsv" --genomes "$RUN/tables/genomes.tsv"
        --families "$RUN/tables/families.tsv" --out "$out"
        --flank "$(cfg orthology/flank_bp)" --min-aln "$(cfg orthology/min_aln)"
        --min-pid "$(cfg orthology/min_pid)" --guard "$(cfg orthology/guard_bp)")
    if [ "$(cfg orthology/tandem_mask)" = True ]; then
        mapfile -t trf < <(_orth_trf)
        mapfile -t codes < <(cfg_keys genomes)
        for c in "${codes[@]}"; do
            printf '%s\n' "${trf[@]}" | "$GAWK" -F '\t' -v c="$c" '$1 == c { f = 1 } END { exit !f }' \
                || lost=${lost:+$lost, }$c
        done
        [ -z "$lost" ] || stage_die "19_orthology: orthology.tandem_mask needs the TRF calls of every genome; none for $lost"
        for p in "${trf[@]}"; do argv+=(--trf-dat "${p%%$'\t'*}=${p#*$'\t'}"); done
    fi
    if [ "$(cfg orthology/single_copy)" = True ]; then
        argv+=(--single-copy --copy-min-pid "$(cfg orthology/copy_min_pid)" --copy-min-aln "$(cfg orthology/copy_min_aln)"
               --threads "$(cfg expansion/jobs)")
    fi
    if _orth_float_gt0 "$(cfg orthology/clear_margin)"; then
        argv+=(--clear-margin "$(cfg orthology/clear_margin)")
        [ "$(cfg orthology/reciprocal)" = True ] && argv+=(--reciprocal)
    fi
    steps=$(cfg orthology/extend_steps); steps=${steps#[}; steps=${steps%]}; steps=${steps//, /,}
    if [ "$(cfg orthology/repeat_mask)" = True ]; then
        if [ "$(cfg orthology/single_copy)" != True ] || [ "$(cfg orthology/tandem_mask)" != True ]; then
            stage_die "19_orthology: orthology.repeat_mask needs single_copy and tandem_mask"
        fi
        argv+=(--repeat-mask --mask-min-aln "$(cfg orthology/mask_min_aln)")
        [ -z "$steps" ] || argv+=(--extend-steps "$steps")
    elif [ -n "$steps" ]; then
        stage_die "19_orthology: orthology.extend_steps needs repeat_mask"
    fi
    sh_run 19_orthology "" -- "${argv[@]}"
    # lines of ortho_pairs.tsv minus the header: sum(1 for _ in open(..)) - 1, text-mode line count
    _py_need_file "$out/ortho_pairs.tsv"
    pairs=$("$GAWK" 'BEGIN { RS = "\r\n|\r|\n" } END { print NR - 1 }' "$out/ortho_pairs.tsv")
    manifest_set ortholog_pairs "$pairs"
    say "  ortholog locus pairs: $pairs"
    emit_outputs "$out/"{flanks.json,flank_orthologs.json,ortho_pairs.tsv}
}
