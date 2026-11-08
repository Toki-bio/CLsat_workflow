# st_03_motifs.sh - stage 03_motifs (run.py:484-507): Step 3 (motifs from the TRF tables) and Step 3b (monomerize)
. "$SATPIPE_LIB/st_common_wp6.sh"

plan_03_motifs() {                                                             # run.py:500-504
    local c
    wp6_codes
    for c in "${WP6_CODES[@]}"; do emit_inputs "$RUN/trf/$c.dat"; done        # _dats (run.py:480-481)
    if [ -e "$RUN/tables/genomes.tsv" ]; then emit_inputs "$RUN/tables/genomes.tsv"; fi
    cfg_subtree discovery/period period
    cfg_subtree discovery/min_copies min_copies
    cfg_subtree discovery/monomerize monomerize
    cfg_subtree discovery/monomer_min_id monomer_min_id
    param_set step3b str runs-v3
}

st_03_motifs() {
    local m c cons="" fam
    rm -rf -- "$RUN/motifs" || true                                           # shutil.rmtree(ignore_errors=True)
    m=$(ctx_d motifs)
    wp6_disc_env
    sh_run 03_motifs "" "${WP6_DISC_ENV[@]}" -- bash "$SCRIPTS/01_discovery_clustering/Step3.sh" "$RUN/trf" "$m"
    wp6_codes
    if [ "$(cfg discovery/monomerize)" = True ]; then
        wp6_families
        for fam in "${WP6_FAMS[@]}"; do cons+=${cons:+:}$(cfg "families/$fam/consensus"); done
        # Python runs sys.executable; the bash driver runs python3 (PORT_SPEC 5b: argv[0] normalised)
        for c in "${WP6_CODES[@]}"; do
            sh_run 03_motifs "" PERIOD_MIN="$(cfg discovery/period/0)" MONOMER_MIN_ID="$(cfg discovery/monomer_min_id)" \
                MONOMER_GENOMES="$RUN/tables/genomes.tsv" MONOMER_CONS="$cons" -- \
                python3 "$SCRIPTS/01_discovery_clustering/Step3b_monomerize.py" "$m/$c.fa"
        done
    fi
    for c in "${WP6_CODES[@]}"; do emit_outputs "$m/$c.fa"; done
}
