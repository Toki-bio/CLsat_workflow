# st_16_regions.sh - stage 16_regions (run.py:921-948), incl. the haplotype-duplicate rule (commit 8687b90, run.py:928-939).
# Contract: lib/satpipe/stages.sh header (WP5).

# shellcheck source=subx.sh
. "$SATPIPE_LIB/subx.sh"

plan_16_regions() {                                                    # run.py:945-948
    emit_inputs "$RUN/catalogue/"{MASTER.merged.bed,MASTER.loci.segmented.bed,MASTER.final_assignment.tsv}
    genome_files
    # hapdup enters the params only when hapdup.pairs is a non-empty list (run.py:947-948), as the whole block
    if [ -n "$(cfg_keys hapdup/pairs 2>/dev/null)" ]; then cfg_subtree hapdup hapdup; fi
}

st_16_regions() {
    local reg cat=$RUN/catalogue pairs="" k removed
    reg=$(ctx_d regions)                                               # run.py:923
    sh_run 16_regions "" -- bash "$SCRIPTS/05_subfamily_terminal_boundary/build_locus_master.sh" \
        "$cat/MASTER.merged.bed" "$cat" "$reg" "$RUN/tables/families.tsv"                       # run.py:925-926
    if [ -n "$(cfg_keys hapdup/pairs 2>/dev/null)" ]; then            # run.py:927-928 (hd.get('pairs') truthy)
        while IFS= read -r k; do                                       # ','.join('%s:%s' % (a, b)) (run.py:934)
            pairs+=${pairs:+,}$(cfg "hapdup/pairs/$k/0"):$(cfg "hapdup/pairs/$k/1")
        done < <(cfg_keys hapdup/pairs)
        # env: str() of the block values (schema defaults fill them, so hd.get(.., default) never falls back);
        # THREADS = os.environ.get('SLURM_CPUS_PER_TASK', '8') of the driver (run.py:935-936)
        sh_run 16_regions "" EDGE_TOL="$(cfg hapdup/edge_tol)" MIN_QCOV="$(cfg hapdup/min_qcov)" \
            MIN_IDENT="$(cfg hapdup/min_ident)" MIN_TOVL="$(cfg hapdup/min_tovl)" THREADS="${SLURM_CPUS_PER_TASK-8}" -- \
            bash "$SCRIPTS/06_catalogue_finalize/hapdup/hapdup_filter.sh" \
            "$reg/locus_master.tsv" "$RUN/tables/genomes.tsv" "$reg" "$pairs"                  # run.py:931-934
        # lines whose last whitespace-separated word is REMOVED (l.split()[-1:] == ["REMOVED"], run.py:937), text mode
        _py_need_file "$reg/hapdup.tsv"
        removed=$("$GAWK" 'BEGIN { RS = "\r\n|\r|\n"; n = 0 }
            { sub(/[ \t\v\f\034-\037]+$/, "") } $0 ~ /(^|[ \t\v\f\034-\037])REMOVED$/ { n++ } END { print n }' "$reg/hapdup.tsv")
        manifest_set hapdup_removed "$removed"                         # run.py:938
        say "  hapdup: $removed loci removed as haplotype duplicates (regions/hapdup.tsv)"     # run.py:939
    fi
    sh_run 16_regions "" -- bash "$SCRIPTS/06_catalogue_finalize/subfamily_rerun/run_extract.sh" "$reg" \
        "$RUN/tables/genomes.tsv"                                                               # run.py:940-941
    emit_outputs "$reg/"{locus_master.tsv,manifest.tsv,all_regions.fa}                         # run.py:942
}
