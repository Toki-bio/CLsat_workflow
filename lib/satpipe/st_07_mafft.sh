# st_07_mafft.sh - stage 07_mafft (run.py:548-566): Step 7 in every cluster folder, min(n, expansion.jobs) at a time
. "$SATPIPE_LIB/st_common_wp6.sh"

plan_07_mafft() {                                                              # run.py:563-566
    local c cb=$RUN/motifs/clean/clusters_by_file
    wp6_codes
    for c in "${WP6_CODES[@]}"; do emit_inputs "$RUN/motifs/clean/$c.fa.clean.fas.rot"; done
    for c in "${WP6_CODES[@]}"; do if [ -d "$cb/$c.fa.clean.fas" ]; then emit_inputs "$cb/$c.fa.clean.fas"; fi; done
    cfg_subtree discovery/cluster_id id
    cfg_subtree discovery/strand strand
}

st_07_mafft() {
    local cb=$RUN/motifs/clean/clusters_by_file script=$SCRIPTS/01_discovery_clustering/Step7.sh
    local jobs n i running=0 e old
    local -a all=() dirs=()
    [ -d "$cb" ] || stage_die "FileNotFoundError: [Errno 2] No such file or directory: '$cb'"
    old=$(shopt -p nullglob dotglob || true)
    shopt -s nullglob dotglob
    all=("$cb"/*)
    eval "$old"
    for e in "${all[@]}"; do if [ -d "$e" ]; then dirs+=("$e"); fi; done   # the glob is sorted bytewise (LC_ALL=C)
    n=${#dirs[@]}
    jobs=$(cfg expansion/jobs)
    if [ "$n" -lt "$jobs" ]; then jobs=$n; fi
    if [ "$jobs" -lt 1 ]; then jobs=1; fi
    # ThreadPoolExecutor.map: every job runs; the first failure in submission order is raised (run.py:557-560)
    for ((i = 0; i < n; i++)); do
        if [ "$running" -ge "$jobs" ]; then wait -n || true; running=$((running - 1)); fi
        rm -f -- "$RUN_TMP/wp6_mafft_err.$i"
        ( STAGE_ERRFILE=$RUN_TMP/wp6_mafft_err.$i
          sh_run 07_mafft "${dirs[$i]}" log_name="07_mafft_${dirs[$i]##*/}" -- bash "$script" ) &
        running=$((running + 1))
    done
    wait || true
    for ((i = 0; i < n; i++)); do
        if [ -f "$RUN_TMP/wp6_mafft_err.$i" ]; then stage_die_file "$RUN_TMP/wp6_mafft_err.$i"; fi
    done
    return 0
}
