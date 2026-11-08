# st_09_assign.sh - stage 09_assign (run.py:618-717): curated cluster -> family table check, cluster report, Step 9.
# Stage contract: lib/satpipe/stages.sh header. The table/report logic is lib/satpipe/assign.awk.

. "$SATPIPE_LIB/pyio.sh"

# run.py:717: inputs = [family_assignment, cross/cross_species_centroids.fa] + _cons_files; params {}
plan_09_assign() {
    emit_inputs "$(cfg family_assignment)" "$RUN/cross/cross_species_centroids.fa"
    cons_files
}

st_09_assign() {
    local x=$RUN/cross cc=$RUN/cross/cross_clusters table job res ssf out f s k
    table=$(cfg family_assignment)
    # os.listdir(cc) (run.py:626): errors as Python's
    if [ ! -e "$cc" ]; then stage_die "[Errno 2] No such file or directory: $(py_repr_q "$cc")"; fi
    if [ ! -d "$cc" ]; then stage_die "[Errno 20] Not a directory: $(py_repr_q "$cc")"; fi
    job=$RUN_TMP/assign.job; res=$RUN_TMP/assign.out; ssf=$RUN_TMP/assign.ssearch
    {
        job_row set cc "$cc"; job_row set table "$table"; job_row set refs "$RUN/tables/refs.rot"
        job_row set cent "$x/cross_species_centroids.fa"; job_row set uc "$x/cross_species.uc"
        job_row set q "$x/centroids_query.fa"; job_row set lib "$x/consensus_x2.fa"; job_row set rep "$x/cluster_report.tsv"
        job_row set ss "$ssf"
        local -a ents=()
        local old; old=$(shopt -p nullglob dotglob || true)
        shopt -s nullglob dotglob; ents=("$cc"/clust_*); eval "$old"
        for f in "${ents[@]}"; do s=${f##*/}; job_row clu "${s#clust_}" "$cc/$s"; done
        while IFS= read -r k; do job_row fam "$k" "$(cfg "families/$k/consensus")"; done < <(cfg_keys families)
    } > "$job"
    # query + library (run.py:625-667)
    { job_row set mode prep; cat -- "$job"; } > "$job.prep"
    ASSIGN_JOB=$job.prep "$GAWK" -f "$PYJSON" -f "$SATPIPE_LIB/pytext.awk" -f "$SATPIPE_LIB/assign.awk" < /dev/null > "$res" \
        || stage_die "stage 09_assign: assign.awk failed"
    py_gawk_result "$res"
    # one ssearch36 run of all centroids vs the doubled consensuses (run.py:668-669), stdout captured
    sh_capture out 09_assign -- "$(tool ssearch36)" -n -z 3 -m 8C -E 1000 -b 1 -d 1 "$x/centroids_query.fa" "$x/consensus_x2.fa"
    printf '%s' "$out" > "$ssf"
    { job_row set mode report; cat -- "$job"; } > "$job.report"
    ASSIGN_JOB=$job.report "$GAWK" -f "$PYJSON" -f "$SATPIPE_LIB/pytext.awk" -f "$SATPIPE_LIB/assign.awk" < /dev/null > "$res" \
        || stage_die "stage 09_assign: assign.awk failed"
    py_gawk_result "$res"
    s=$(ctx_d seeds)
    sh_run 09_assign "" -- bash "$SCRIPTS/02_seed_expansion_scan/Step9_family_fasta.sh" "$cc" "$table" "$s" "$RUN/tables/refs.rot"
    while IFS= read -r k; do emit_outputs "$s/$k.members.fa"; done < <(cfg_keys families)
    emit_outputs "$x/cluster_report.tsv"
}
