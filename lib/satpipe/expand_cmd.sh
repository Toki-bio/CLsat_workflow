# expand_cmd.sh - _expand_cmd (run.py:753-768): the argv of clisat_expand_then_scan_merged.sh for one family.
# Used by stages 11_expand, 13_rescan, 15_catalogue and the viewer's flank scans (viewer.py:140).
#   expand_cmd BED FAM OUT_PREFIX SCAN_ONLY      sets the array EXPAND_CMD; SCAN_ONLY 1 appends --scan-only
# Needs cfg / cfg_type (config.sh), $RUN and $SCRIPTS. Values are the Python str() texts of the config leaves (cfg); the two
# window products are int x int (scan.win_mono / scan.step_mono x families.<fam>.mono, all ints in the schema), so bash
# arithmetic gives str(int) exactly.

[ -n "${_SATPIPE_EXPAND_CMD_SH:-}" ] && return 0
_SATPIPE_EXPAND_CMD_SH=1

expand_cmd() {
    local bed=$1 fam=$2 pre=$3 scan_only=$4 mono
    mono=$(cfg "families/$fam/mono")
    EXPAND_CMD=(bash "$SCRIPTS/02_seed_expansion_scan/clisat_expand_then_scan_merged.sh"
        "$bed" "$RUN/tables/genomes.tsv" "$(cfg "families/$fam/consensus")" "$mono" "$pre"
        --step-monomers "$(cfg expansion/step_mono)" --win-monomers "$(cfg expansion/win_mono)"
        --query-repeat "$(cfg expansion/query_repeat)"
        --min-aln-frac "$(cfg expansion/min_aln_frac)" --min-bits "$(cfg expansion/min_bits)" --max-steps "$(cfg expansion/max_steps)"
        --probe-monomers "$(cfg expansion/probe_mono)" --probe-min-step "$(cfg expansion/probe_min_step)"
        --scan-win-bp "$(( $(cfg scan/win_mono) * mono ))" --scan-step-bp "$(( $(cfg scan/step_mono) * mono ))"
        --scan-min-bits "$(cfg scan/min_bits)" --jobs "$(cfg expansion/jobs)")
    # e.get('to_scaffold_ends') is truthy only as True (a bool in the schema)
    if [ "$(cfg_type expansion/to_scaffold_ends 2>/dev/null)" = bool ] && [ "$(cfg expansion/to_scaffold_ends)" = True ]; then
        EXPAND_CMD+=(--to-scaffold-ends)
    fi
    if [ "$scan_only" = 1 ]; then EXPAND_CMD+=(--scan-only); fi
    return 0
}
