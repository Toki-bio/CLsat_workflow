# st_12_master.sh - stage 12_master (run.py:788-799): merge the expansions into MASTER.merged.bed.
# Stage contract: lib/satpipe/stages.sh header.

. "$SATPIPE_LIB/pyio.sh"

# run.py:798-799: expand/out_<fam>.expanded_merged.bed per family + expand/expansions.tsv when it is a file; params {}
plan_12_master() {
    local f
    while IFS= read -r f; do emit_inputs "$RUN/expand/out_$f.expanded_merged.bed"; done < <(cfg_keys families)
    if [ -f "$RUN/expand/expansions.tsv" ]; then emit_inputs "$RUN/expand/expansions.tsv"; fi
}

st_12_master() {
    local m n
    m=$(ctx_d master)
    sh_run 12_master "" -- bash "$SCRIPTS/03_master_arbitration/build_master.sh" "$RUN/expand/expansions.tsv" "$m"
    n=$(py_nlines "$m/MASTER.merged.bed")
    manifest_set master_loci "$n"
    say "  MASTER loci: $n"
    emit_outputs "$m/MASTER.merged.bed"
}
