# st_00_tables.sh - stage 00_tables (run.py:300-341): tables/genomes.tsv families.tsv members.tsv refs.rot, genome .fai
. "$SATPIPE_LIB/st_common_wp6.sh"

# rotated_consensus given (truthy path), run.py:323, :335
_st00_rot() { case $(cfg_type references/rotated_consensus 2>/dev/null) in path|str) cfg references/rotated_consensus ;; esac; }

plan_00_tables() {                                                             # run.py:333-338
    local rot
    cons_files
    rot=$(_st00_rot)
    if [ -n "$rot" ]; then emit_inputs "$rot"; fi
    cfg_subtree genomes genomes
    cfg_subtree families families
    param_cfg_or_null subfamily/members members
}

st_00_tables() {
    local t rot refs code g fam
    t=$(ctx_d tables)
    wp6_codes; wp6_families
    # genomes.tsv / families.tsv (mono is an int: %d = its text) / members.tsv (run.py:302-312)
    for code in "${WP6_CODES[@]}"; do printf '%s\t%s\n' "$code" "$(cfg "genomes/$code")"; done > "$t/genomes.tsv"
    for fam in "${WP6_FAMS[@]}"; do printf '%s\t%s\t%s\n' "$fam" "$(cfg "families/$fam/mono")" "$(cfg "families/$fam/consensus")"; done > "$t/families.tsv"
    : > "$t/members.tsv"
    if [ "$(cfg_type subfamily/members 2>/dev/null)" = map ]; then
        local -a mf=() mm=()
        local i
        mapfile -t mf < <(cfg_keys subfamily/members)
        for fam in "${mf[@]}"; do
            mapfile -t mm < <(cfg_keys "subfamily/members/$fam")
            for i in "${mm[@]}"; do printf '%s\t%s\n' "$fam" "$(cfg "subfamily/members/$fam/$i")" >> "$t/members.tsv"; done
        done
    fi
    # Step-8 references (run.py:314-326)
    rot=$(_st00_rot)
    refs=$t/refs.rot
    if [ -n "$rot" ]; then
        cat -- "$rot" > "$refs" || stage_die "[Errno 2] No such file or directory: '$rot'"     # shutil.copyfile
        manifest_set refs_rot "$(S=$rot pyjson_gawk 'BEGIN { printf "{\"source\": %s, \"derived\": false}", json_str(ENVIRON["S"]) }')"
    else
        : > "$refs"
        for fam in "${WP6_FAMS[@]}"; do
            # ''.join(l.strip() for l in open(cons) if not l.startswith('>')).upper(), then minrot (run.py:111-113)
            F=$fam "$GAWK" 'BEGIN { RS = "\n" }
                { sub(/\r$/, "") }
                /^>/ { next }
                { gsub(/^[ \t\n\r\f\v\034-\037]+|[ \t\n\r\f\v\034-\037]+$/, ""); s = s $0 }
                END {
                    s = toupper(s); n = length(s)
                    if (n == 0) { print "ValueError: min() iterable argument is empty" > "/dev/stderr"; exit 3 }
                    d = s s; best = s
                    for (i = 2; i <= n; i++) { r = substr(d, i, n); if (r < best) best = r }
                    printf ">%s\n%s\n", ENVIRON["F"], best
                }' "$(cfg "families/$fam/consensus")" >> "$refs" || stage_die "ValueError: min() iterable argument is empty"
        done
        manifest_set refs_rot '{"source": null, "derived": true, "note": "Step-5 (lexicographically minimal) rotation of each family consensus; the published Darevskia run used the file CLsat_cons.rot"}'
    fi
    # genome indexes (run.py:327-329)
    for code in "${WP6_CODES[@]}"; do
        g=$(cfg "genomes/$code")
        if [ ! -f "$g.fai" ]; then sh_run 00_tables "" -- "$(tool samtools)" faidx "$g"; fi
    done
    emit_outputs "$t/genomes.tsv" "$t/families.tsv" "$t/members.tsv" "$t/refs.rot"
}
