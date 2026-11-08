# st_15_catalogue.sh - stage 15_catalogue (run.py:844-916): final catalogue = arbitration output + manual loci (tagged
# MANUAL, numbered after the pipeline loci, scanned with every consensus for their bits). Stage contract: stages.sh header.
# The manual-loci logic is lib/satpipe/manual.awk.

. "$SATPIPE_LIB/pyio.sh"
. "$SATPIPE_LIB/expand_cmd.sh"

# ctx.cfg.data.get('manual_loci') when truthy (a non-empty path), else ""
_manual_loci_path() {
    case $(cfg_type manual_loci 2>/dev/null) in path|str) cfg manual_loci ;; esac
    return 0
}

# run.py:915-916: arbitration/MASTER.final_assignment.tsv + [manual_loci] when set + _cons_files + _genome_files; params {}
plan_15_catalogue() {
    local m
    emit_inputs "$RUN/arbitration/MASTER.final_assignment.tsv"
    m=$(_manual_loci_path)
    if [ -n "$m" ]; then emit_inputs "$m"; fi
    cons_files; genome_files
}

_manual_gawk() { MANUAL_JOB=$1 "$GAWK" -f "$PYJSON" -f "$SATPIPE_LIB/pytext.awk" -f "$SATPIPE_LIB/manual.awk" < /dev/null; }

st_15_catalogue() {
    local cat arb=$RUN/arbitration manual n_master job res f k pre json="[]" nseg nman
    local -a fams=()
    cat=$(ctx_d catalogue)
    for f in MASTER.loci.segmented.bed MASTER.final_assignment.tsv MASTER.locus_bits.tsv MASTER.loci.bed; do
        py_copyfile "$arb/$f" "$cat/$f"
    done
    py_copyfile "$RUN/master/MASTER.merged.bed" "$cat/MASTER.merged.bed"
    while IFS= read -r f; do fams+=("$f"); done < <(cfg_keys families)
    for f in "${fams[@]}"; do py_copyfile "$RUN/rescan/MASTER_$f.bits.tsv" "$cat/MASTER_$f.bits.tsv"; done
    manual=$(_manual_loci_path)
    if [ -n "$manual" ]; then
        # os.path.getsize(manual) (run.py:857)
        if [ ! -e "$manual" ]; then stage_die "[Errno 2] No such file or directory: $(py_repr_q "$manual")"; fi
    fi
    if [ -n "$manual" ] && [ -s "$manual" ]; then
        n_master=$(py_nlines "$cat/MASTER.loci.bed")
        job=$RUN_TMP/manual.job; res=$RUN_TMP/manual.out
        {
            job_row set manual "$manual"
            while IFS= read -r k; do job_row gen "$k"; done < <(cfg_keys genomes)
            for f in "${fams[@]}"; do job_row fam "$f"; done
        } > "$job"
        : > "$RUN_TMP/manual.state"
        { job_row set mode prep; job_row set mbed "$cat/manual.bed"; cat -- "$job"; } > "$job.prep"
        _manual_gawk "$job.prep" > "$res" || stage_die "stage 15_catalogue: manual.awk failed"
        py_gawk_result "$res"
        for f in "${fams[@]}"; do
            pre=$cat/manual_$f
            expand_cmd "$cat/manual.bed" "$f" "$pre" 1
            sh_run 15_catalogue "" log_name="15_catalogue_$f" -- "${EXPAND_CMD[@]}"
            { job_row set mode fam; job_row set fam "$f"; job_row set bits "$pre.bits.tsv"
              job_row set out "$cat/MASTER_$f.bits.tsv"; job_row set state "$RUN_TMP/manual.state"; cat -- "$job"; } > "$job.fam"
            _manual_gawk "$job.fam" > "$res" || stage_die "stage 15_catalogue: manual.awk failed"
            py_gawk_result "$res"
        done
        { job_row set mode final; job_row set state "$RUN_TMP/manual.state"; job_row set nmaster "$n_master"
          job_row set mb "$cat/MASTER.merged.bed"; job_row set lb "$cat/MASTER.loci.bed"
          job_row set sb "$cat/MASTER.loci.segmented.bed"; job_row set bb "$cat/MASTER.locus_bits.tsv"
          job_row set fa "$cat/MASTER.final_assignment.tsv"; cat -- "$job"; } > "$job.final"
        _manual_gawk "$job.final" > "$res" || stage_die "stage 15_catalogue: manual.awk failed"
        py_gawk_result "$res"
        json=$("$GAWK" -F '\t' '$1 == "JSON" { print substr($0, 6) }' "$res")
    fi
    manifest_set manual_loci "$json"
    nseg=$(py_nlines "$cat/MASTER.loci.segmented.bed")
    manifest_set catalogue_loci "$nseg"
    nman=$("$GAWK" -f "$PYJSON" -e '{ if (!json_parse($0, "J")) exit 1; print JN["J"] }' <<< "$json")
    say "  catalogue: $nseg segmented loci ($nman manual)"
    # a header-less copy of the assignment for the viewer (its name in the viewer repo)
    py_copyfile "$cat/MASTER.final_assignment.tsv" "$cat/MASTER.final_assignment.clean.tsv"
    emit_outputs "$cat/MASTER.merged.bed" "$cat/MASTER.loci.segmented.bed" "$cat/MASTER.final_assignment.tsv" "$cat/MASTER.locus_bits.tsv"
}
