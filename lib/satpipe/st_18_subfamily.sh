# st_18_subfamily.sh - stage 18_subfamily (run.py:969-988), subx_dir in subx.sh. Contract: lib/satpipe/stages.sh header (WP5).

# shellcheck source=subx.sh
. "$SATPIPE_LIB/subx.sh"

plan_18_subfamily() {                                                  # run.py:986-988
    emit_inputs "$RUN/regions/"{locus_master.tsv,terminal_monomers.tsv,manifest.tsv,all_regions.fa} "$RUN/tables/members.tsv"
    subx_inputs
    param_cfg_or_null subfamily subfamily                              # {'subfamily': cfg.get('subfamily')}
}

st_18_subfamily() {
    local reg=$RUN/regions                                             # run.py:971 (no mkdir)
    subx_dir                                                           # run.py:973 (writes tables/subx/ for a FASTA)
    if [ -n "$SUBX" ]; then
        sh_run 18_subfamily "" -- bash "$SCRIPTS/05_subfamily_terminal_boundary/assign_subfamily.sh" "$reg" "$SUBX" \
            "$RUN/tables/members.tsv"                                                           # run.py:974-976
    else
        say '  no sub-subfamily consensus directory: writing an empty locus_subfamily.tsv'    # run.py:978
        [ -d "$reg" ] || stage_die "[Errno 2] No such file or directory: '$reg/locus_subfamily.tsv'"
        printf 'locus\tmajor\tsubfamily\tpct_id\taln_len\tbits\tsecond\tsecond_bits\tdelta\n' > "$reg/locus_subfamily.tsv" \
            || stage_die "cannot write $reg/locus_subfamily.tsv"                               # run.py:979-980
    fi
    # CONFIDENT_DELTA = str(sub.get('confident_delta', 10)); the schema default fills it (run.py:981-982)
    sh_run 18_subfamily "" CONFIDENT_DELTA="$(_p8_confident_delta)" -- \
        python3 "$SCRIPTS/05_subfamily_terminal_boundary/build_subfam_term_json.py" "$reg" "$RUN/tables/families.tsv"
    emit_outputs "$reg/"{locus_subfamily.tsv,subfamily.json,terminal_monomers.json}          # run.py:983
}
