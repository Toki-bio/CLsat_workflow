# st_17_terminal.sh - stage 17_terminal (run.py:951-966). Contract: lib/satpipe/stages.sh header (WP5).

plan_17_terminal() {                                                   # run.py:964-966
    emit_inputs "$RUN/regions/"{locus_master.tsv,manifest.tsv,all_regions.fa} "$RUN/tables/families.tsv"
    cons_files
    genome_files
    cfg_subtree terminal terminal                                      # {'terminal': cfg['terminal']}
}

st_17_terminal() {
    local reg=$RUN/regions                                             # run.py:953 (no mkdir)
    # env = ctx.env(WINDOW_MONO=t['window_mono'], EDGE_TOL=t['edge_tol']) for the first command only (run.py:955-957)
    sh_run 17_terminal "" WINDOW_MONO="$(cfg terminal/window_mono)" EDGE_TOL="$(cfg terminal/edge_tol)" -- \
        bash "$SCRIPTS/05_subfamily_terminal_boundary/terminal_monomers.sh" "$reg" "$RUN/tables/genomes.tsv" "$RUN/tables/families.tsv"
    sh_run 17_terminal "" -- bash "$SCRIPTS/05_subfamily_terminal_boundary/boundary_phase.sh" "$reg" \
        "$RUN/tables/families.tsv"                                                              # run.py:958-959
    sh_run 17_terminal "" -- python3 "$SCRIPTS/05_subfamily_terminal_boundary/build_json.py" "$reg"   # run.py:960
    emit_outputs "$reg/"{terminal_monomers.tsv,locus_boundary_phase.tsv,boundary_phase.json}  # run.py:961
}
