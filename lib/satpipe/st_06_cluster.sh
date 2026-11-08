# st_06_cluster.sh - stage 06_cluster (run.py:533-545): Step 6 (vsearch per genome) into a fresh clusters_by_file/
. "$SATPIPE_LIB/st_common_wp6.sh"

plan_06_cluster() {                                                            # run.py:544-545
    local c
    wp6_codes
    for c in "${WP6_CODES[@]}"; do emit_inputs "$RUN/motifs/clean/$c.fa.clean.fas.rot"; done
    cfg_subtree discovery/cluster_id id
    cfg_subtree discovery/strand strand
}

st_06_cluster() {
    local clean=$RUN/motifs/clean                                             # no ctx.d here (run.py:535)
    rm -rf -- "$clean/clusters_by_file" || true
    wp6_disc_env
    sh_run 06_cluster "" "${WP6_DISC_ENV[@]}" -- bash "$SCRIPTS/01_discovery_clustering/Step6.sh" "$clean" "$clean/clusters_by_file"
    return 0
}
