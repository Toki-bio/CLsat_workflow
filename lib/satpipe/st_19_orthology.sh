# st_19_orthology.sh - stage 19_orthology (run.py:991-1007). Contract: lib/satpipe/stages.sh header (WP5).
# The argv is built in an array in Python's order so that option groups (branch flank-single-copy, WP13) can be
# appended after the four base options without touching the rest.

# shellcheck source=subx.sh
. "$SATPIPE_LIB/subx.sh"

plan_19_orthology() {                                                  # run.py:1005-1007
    emit_inputs "$RUN/regions/locus_master.tsv" "$RUN/tables/genomes.tsv" "$RUN/tables/families.tsv"
    cons_files
    genome_files
    cfg_subtree orthology orthology                                    # {'orthology': cfg['orthology']}
}

st_19_orthology() {
    local out pairs
    out=$(ctx_d ortho)                                                 # run.py:994 (after reading the block)
    local -a argv=(python3 "$SCRIPTS/04_flank_orthology/rebuild_ortho_segmented.py"
        --loci "$RUN/regions/locus_master.tsv" --genomes "$RUN/tables/genomes.tsv"
        --families "$RUN/tables/families.tsv" --out "$out"
        --flank "$(cfg orthology/flank_bp)" --min-aln "$(cfg orthology/min_aln)"
        --min-pid "$(cfg orthology/min_pid)" --guard "$(cfg orthology/guard_bp)")              # run.py:995-998
    # (WP13: flank-single-copy option groups go here, appended to argv in the branch's order)
    sh_run 19_orthology "" -- "${argv[@]}"
    # lines of ortho_pairs.tsv minus the header: sum(1 for _ in open(..)) - 1 (run.py:999), text-mode line count
    _py_need_file "$out/ortho_pairs.tsv"
    pairs=$("$GAWK" 'BEGIN { RS = "\r\n|\r|\n" } END { print NR - 1 }' "$out/ortho_pairs.tsv")
    manifest_set ortholog_pairs "$pairs"                               # run.py:1000
    say "  ortholog locus pairs: $pairs"                               # run.py:1001
    emit_outputs "$out/"{flanks.json,flank_orthologs.json,ortho_pairs.tsv}                     # run.py:1002
}
