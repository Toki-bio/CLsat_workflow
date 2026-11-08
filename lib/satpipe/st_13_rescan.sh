# st_13_rescan.sh - stage 13_rescan (run.py:802-817): scan-only rescan of the MASTER loci with every family consensus.
# Stage contract: lib/satpipe/stages.sh header.

. "$SATPIPE_LIB/expand_cmd.sh"

# run.py:816-817: master/MASTER.merged.bed + _cons_files + _genome_files; params {'scan': block, 'jobs': expansion.jobs}
plan_13_rescan() {
    emit_inputs "$RUN/master/MASTER.merged.bed"
    cons_files; genome_files
    cfg_subtree scan scan
    cfg_subtree expansion/jobs jobs
}

st_13_rescan() {
    local r pre f master=$RUN/master/MASTER.merged.bed
    local -a fams=()
    r=$(ctx_d rescan)
    while IFS= read -r f; do fams+=("$f"); done < <(cfg_keys families)
    : > "$r/arbitration_inputs.tsv"               # opened for writing before the loop (run.py:807)
    for f in "${fams[@]}"; do
        pre=$r/MASTER_$f
        expand_cmd "$master" "$f" "$pre" 1
        sh_run 13_rescan "" log_name="13_rescan_$f" -- "${EXPAND_CMD[@]}"
        printf '%s\t%s\t%s\t%s\n' "$f" "$(cfg "families/$f/mono")" "$pre.bits.tsv" "$RUN/seeds/$f.bed" >> "$r/arbitration_inputs.tsv"
    done
    for f in "${fams[@]}"; do emit_outputs "$r/MASTER_$f.bits.tsv"; done
    emit_outputs "$r/arbitration_inputs.tsv"
}
