# st_02b_large.sh - stage 02b_large (run.py:413-466): chunked TRF of the scaffolds 02_trf left out, classified against
# the family consensuses; kept candidates become large/<family>.recovered.bed seeds.
. "$SATPIPE_LIB/st_common_wp6.sh"

plan_02b_large() {                                                             # run.py:461-463
    emit_inputs "$RUN/trf/large_scaffolds.tsv" "$RUN/tables/families.tsv"
    cons_files
    genome_files
    cfg_subtree finalisation fin
    cfg_subtree discovery/trf trf
    cfg_subtree discovery/min_copies min_copies
}

st_02b_large() {
    local code scf ln d cand uns log rec row fam j n
    local -a rows=() outs=() envs=()
    rec=$RUN_TMP/wp6_recovered
    : > "$rec"
    wp6_families
    # rows = read_tsv(large_scaffolds.tsv, 3); `for code, scf, ln in rows` needs exactly 3 fields (run.py:416, 421)
    mapfile -t rows < <(wp6_read_tsv "$RUN/trf/large_scaffolds.tsv" 3)
    envs=(CHUNK_BP="$(cfg finalisation/chunk_bp)" OVERLAP_BP="$(cfg finalisation/overlap_bp)" TRF_PARAMS="$(cfg discovery/trf)"
          TRF_TIMEOUT="$(cfg finalisation/trf_timeout)" MIN_PIECE_BP="$(cfg finalisation/min_piece_bp)" JOBS="$(cfg expansion/jobs)")
    for row in "${rows[@]}"; do
        IFS=$'\t' read -r code scf ln j <<< "$row" || true
        if [ -n "$j" ]; then stage_die "ValueError: too many values to unpack (expected 3)"; fi
        d=$(ctx_d large "$code" "$scf")
        log=02b_large_${code}_$scf
        sh_run 02b_large "" log_name="$log" "${envs[@]}" -- \
            bash "$SCRIPTS/06_catalogue_finalize/85mb_recovery/run_chunked.sh" "$(cfg "genomes/$code")" "$scf" "$d"   # run.py:424-428
        cand=$d/candidates.tsv
        sh_run 02b_large "" log_name="$log" -- \
            python3 "$SCRIPTS/06_catalogue_finalize/85mb_recovery/classify_candidates.py" \
            "$d/scaf.trf.dedup.dat" "$RUN/tables/families.tsv" "$cand" --min-copies "$(cfg discovery/min_copies)"   # run.py:429-432
        # RECOVERED seed rows (run.py:433-438): family<TAB>bed line; an unknown family is a KeyError at once
        [ -f "$cand" ] || stage_die "[Errno 2] No such file or directory: '$cand'"
        wp6_read_tsv "$cand" 6 | C=$code S=$scf "$GAWK" -F '\t' '$1 != "start" {
                printf "%s\t%s\t%d\t%d\t%s|%s:%d-%d|RECOVERED\t0\t.\n", $5, ENVIRON["S"], $1 - 1, $2, ENVIRON["C"], ENVIRON["S"], $1, $2 }' > "$rec.new"
        fam=$(F="${WP6_FAMS[*]}" "$GAWK" -F '\t' 'BEGIN { n = split(ENVIRON["F"], a, " "); for (i = 1; i <= n; i++) K[a[i]] = 1 }
                !($1 in K) { print $1; exit }' "$rec.new")
        if [ -n "$fam" ]; then stage_die "KeyError: '$fam'"; fi
        cat -- "$rec.new" >> "$rec"
        outs+=("$cand")
        uns=$d/unscanned.tsv
        if [ -f "$uns" ] && [ -s "$uns" ]; then                                # run.py:440-446
            n=$(wp6_read_tsv "$uns" 3 | "$GAWK" 'END { print NR }')
            wp6_manifest_extend unscanned "$(wp6_read_tsv "$uns" 3 | C=$code S=$scf "$GAWK" -f "$PYJSON" -F '\t' -e '
                { o = o (NR > 1 ? ", " : "") sprintf("{\"code\": %s, \"scaffold\": %s, \"piece\": %s, \"start\": %d, \"end\": %d, \"note\": %s}",
                      json_str(ENVIRON["C"]), json_str(ENVIRON["S"]), json_str($1), $2, $3, json_str(NF > 4 ? $5 : "")) }
                END { printf "[%s]", o }')"
            say "  $code $scf: $n piece(s) TRF could not finish (listed in $uns)"
        fi
    done
    d=$(ctx_d large)
    for fam in "${WP6_FAMS[@]}"; do                                            # run.py:447-452
        FAM=$fam "$GAWK" -F '\t' '$1 == ENVIRON["FAM"] { print substr($0, length($1) + 2) }' "$rec" > "$d/$fam.recovered.bed"
        outs+=("$d/$fam.recovered.bed")
    done
    manifest_set large_scaffolds "$(for row in "${rows[@]}"; do printf '%s\n' "$row"; done | "$GAWK" -f "$PYJSON" -F '\t' -e '
        NF { o = o (n++ ? ", " : "") sprintf("{\"code\": %s, \"scaffold\": %s, \"length\": %d}", json_str($1), json_str($2), $3) }
        END { printf "[%s]", o }')"
    emit_outputs "${outs[@]}"
}
