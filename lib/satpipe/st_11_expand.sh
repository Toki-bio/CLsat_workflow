# st_11_expand.sh - stage 11_expand (run.py:771-785): seed expansion per family. Stage contract: lib/satpipe/stages.sh header.

. "$SATPIPE_LIB/expand_cmd.sh"

# run.py:784-785: seeds/<fam>.bed per family + _cons_files + _genome_files;
# params {'expansion': block, 'scan': block, 'mono': {fam: mono}}
plan_11_expand() {
    local f n=0
    while IFS= read -r f; do emit_inputs "$RUN/seeds/$f.bed"; n=$((n + 1)); done < <(cfg_keys families)
    cons_files; genome_files
    cfg_subtree expansion expansion
    cfg_subtree scan scan
    param_set mono map "$n"
    while IFS= read -r f; do cfg_subtree "families/$f/mono" "mono/$f"; done < <(cfg_keys families)
}

st_11_expand() {
    local x pre f
    local -a fams=()
    x=$(ctx_d expand)
    while IFS= read -r f; do fams+=("$f"); done < <(cfg_keys families)
    : > "$x/expansions.tsv"                       # opened for writing before the loop (run.py:775)
    for f in "${fams[@]}"; do
        pre=$x/out_$f
        expand_cmd "$RUN/seeds/$f.bed" "$f" "$pre" 0
        sh_run 11_expand "" log_name="11_expand_$f" -- "${EXPAND_CMD[@]}"
        printf '%s\t%s\n' "$f" "$pre.expanded_merged.bed" >> "$x/expansions.tsv"
    done
    for f in "${fams[@]}"; do emit_outputs "$x/out_$f.expanded_merged.bed"; done
    emit_outputs "$x/expansions.tsv"
}
