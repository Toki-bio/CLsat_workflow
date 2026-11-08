# st_05_rotate.sh - stage 05_rotate (run.py:522-530): Step 5 (minimal rotation, strand)
. "$SATPIPE_LIB/st_common_wp6.sh"

plan_05_rotate() {                                                             # run.py:529-530
    local c
    wp6_codes
    for c in "${WP6_CODES[@]}"; do emit_inputs "$RUN/motifs/$c.fa.clean.fas"; done
    cfg_subtree discovery/strand strand
}

st_05_rotate() {
    local m c
    m=$(ctx_d motifs)
    wp6_disc_env
    sh_run 05_rotate "" "${WP6_DISC_ENV[@]}" -- bash "$SCRIPTS/01_discovery_clustering/Step5.sh" "$m"
    wp6_codes
    for c in "${WP6_CODES[@]}"; do emit_outputs "$m/clean/$c.fa.clean.fas.rot"; done
}
