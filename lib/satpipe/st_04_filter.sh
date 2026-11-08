# st_04_filter.sh - stage 04_filter (run.py:510-519): Step 4 per genome, cwd motifs/
. "$SATPIPE_LIB/st_common_wp6.sh"

plan_04_filter() {                                                             # run.py:517-519
    local c
    wp6_codes
    for c in "${WP6_CODES[@]}"; do emit_inputs "$RUN/motifs/$c.fa"; done
    cfg_subtree discovery/trf trf
    cfg_subtree discovery/simple_repeat_cov cov
}

st_04_filter() {
    local m c
    m=$(ctx_d motifs)
    wp6_codes; wp6_disc_env
    for c in "${WP6_CODES[@]}"; do
        sh_run 04_filter "$m" "${WP6_DISC_ENV[@]}" -- bash "$SCRIPTS/01_discovery_clustering/Step4.sh" "$c.fa"
    done
    for c in "${WP6_CODES[@]}"; do emit_outputs "$m/$c.fa.clean.fas"; done
}
