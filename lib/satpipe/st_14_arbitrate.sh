# st_14_arbitrate.sh - stage 14_arbitrate (run.py:820-841): segmentation and family arbitration.
# Stage contract: lib/satpipe/stages.sh header.

. "$SATPIPE_LIB/pyio.sh"

# run.py:838-841: rescan/arbitration_inputs.tsv, master/MASTER.merged.bed, rescan/MASTER_<fam>.bits.tsv per family;
# params {'arbitration': block, 'scan_win_mono': scan.win_mono, 'bin_min_bits': expansion.min_bits, 'rule': 'bin-majority-v1'}
# (the params' bin_min_bits is ALWAYS expansion.min_bits, whatever arbitration.bin_min_bits says; the block carries that one)
plan_14_arbitrate() {
    local f
    emit_inputs "$RUN/rescan/arbitration_inputs.tsv" "$RUN/master/MASTER.merged.bed"
    while IFS= read -r f; do emit_inputs "$RUN/rescan/MASTER_$f.bits.tsv"; done < <(cfg_keys families)
    cfg_subtree arbitration arbitration
    cfg_subtree scan/win_mono scan_win_mono
    cfg_subtree expansion/min_bits bin_min_bits
    param_set rule str bin-majority-v1
}

st_14_arbitrate() {
    local out n bmb
    # BIN_MIN_BITS = arbitration.bin_min_bits, or expansion.min_bits when it is None (run.py:827): str() of that value,
    # e.g. 50 (schema default, an int) or 50.0 (written in the YAML)
    if [ "$(cfg_type arbitration/bin_min_bits 2>/dev/null)" = null ] || ! cfg_type arbitration/bin_min_bits > /dev/null 2>&1; then
        bmb=$(cfg expansion/min_bits)
    else
        bmb=$(cfg arbitration/bin_min_bits)
    fi
    out=$(ctx_d arbitration)
    sh_run 14_arbitrate "" DELTA="$(cfg arbitration/delta)" DELTA_AMBIG="$(cfg arbitration/delta_ambig)" \
        SEG_BIN="$(cfg arbitration/seg_bin)" SEG_MIN_RUN="$(cfg arbitration/seg_min_run)" \
        SEG_MIN_LOCUS="$(cfg arbitration/seg_min_locus)" SEG_SHORT="$(cfg arbitration/short_segments)" \
        SCAN_WIN_MONO="$(cfg scan/win_mono)" BIN_MIN_BITS="$bmb" \
        BIN_SHARE="$(cfg arbitration/bin_share)" MAX_GAP_MONO="$(cfg arbitration/max_gap_monomers)" -- \
        bash "$SCRIPTS/03_master_arbitration/arbitration_patched.sh" "$RUN/rescan/arbitration_inputs.tsv" \
        "$RUN/master/MASTER.merged.bed" "$out"
    n=$(py_nlines "$out/MASTER.loci.segmented.bed")
    manifest_set segmented_loci "$n"
    say "  segmented loci: $n"
    emit_outputs "$out/MASTER.loci.segmented.bed" "$out/MASTER.final_assignment.tsv" "$out/MASTER.locus_bits.tsv"
}
