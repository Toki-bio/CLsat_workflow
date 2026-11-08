# st_08_cross.sh - stage 08_cross (run.py:571-585): Step 8 cross-genome clustering; centroid md5 and count
. "$SATPIPE_LIB/st_common_wp6.sh"

plan_08_cross() {                                                              # run.py:584-585
    local c
    emit_inputs "$RUN/tables/refs.rot"
    wp6_codes
    for c in "${WP6_CODES[@]}"; do emit_inputs "$RUN/motifs/clean/$c.fa.clean.fas.rot"; done
    cfg_subtree discovery/cross_cluster_id id
    cfg_subtree discovery/strand strand
}

st_08_cross() {
    local x cent md5 n
    x=$(ctx_d cross)
    wp6_disc_env
    sh_run 08_cross "" "${WP6_DISC_ENV[@]}" -- bash "$SCRIPTS/02_seed_expansion_scan/Step8_cross_cluster.sh" \
        "$RUN/motifs/clean/clusters_by_file" "$RUN/tables/refs.rot" "$x"
    cent=$x/cross_species_centroids.fa
    [ -f "$cent" ] || stage_die "FileNotFoundError: [Errno 2] No such file or directory: '$cent'"
    md5=$(md5sum < "$cent") || stage_die "cannot read $cent"
    md5=${md5%% *}
    # read_fasta_names (run.py:76-80): one name per '>' line; a '>' line without a name is an IndexError
    n=$("$GAWK" '/^>/ { h = substr($0, 2); sub(/^[ \t\n\r\f\v\034-\037]+/, "", h); if (h == "") { bad = 1; exit } n++ }
                 END { if (bad) exit 3; print n + 0 }' "$cent") || stage_die "IndexError: list index out of range"
    manifest_set centroids_md5 "\"$md5\""
    manifest_set centroids_n "$n"
    say "  centroids: $n, md5 $md5"
    emit_outputs "$cent"
}
