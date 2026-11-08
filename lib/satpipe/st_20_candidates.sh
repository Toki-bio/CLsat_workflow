# st_20_candidates.sh - stage 20_candidates (run.py:1010-1066). Contract: lib/satpipe/stages.sh header (WP5).

# shellcheck source=subx.sh
. "$SATPIPE_LIB/subx.sh"

plan_20_candidates() {                                                 # run.py:1064-1066
    emit_inputs "$RUN/catalogue/MASTER.merged.bed" "$RUN/tables/refs.rot"
    cons_files
    genome_files
    subx_inputs
    cfg_subtree expansion/min_bits min_bits                            # {'min_bits': .., 'flank': ..} in this order
    cfg_subtree orthology/flank_bp flank
}

st_20_candidates() {
    local cb=$RUN/motifs/clean/clusters_by_file out tables=$RUN/tables f n missing
    out=$(ctx_d candidates)                                            # run.py:1016-1017 (made before the check)
    if [ ! -d "$cb" ]; then
        # seeds given for every family: nothing to screen (run.py:1018-1022); else an error naming the families
        # discovery ran for, sorted (run.py:1023-1024)
        local -A given=()
        if [ "$(cfg_type inputs/seeds 2>/dev/null)" = map ]; then
            while IFS= read -r f; do given[$f]=1; done < <(cfg_keys inputs/seeds)
        fi
        missing=$(while IFS= read -r f; do [ -n "${given[$f]:-}" ] || printf '%s\n' "$f"; done < <(cfg_keys families) \
            | "$GAWK" '{ a[$0] = 1 } END { PROCINFO["sorted_in"] = "@ind_str_asc"; for (k in a) printf "%s%s", (n++ ? ", " : ""), k }')
        if [ -z "$missing" ]; then
            say '  no per-species clusters in this run (seeds were given for every family): no singleton screen'
            return 0
        fi
        stage_die "20_candidates: $cb is missing although discovery ran for $missing"
    fi
    sh_run 20_candidates "" MIN_BITS="$(cfg expansion/min_bits)" -- bash "$SCRIPTS/07_post_catalogue/singleton_screen.sh" \
        "$cb" "$tables/refs.rot" "$out" "$RUN/catalogue/MASTER.merged.bed"                      # run.py:1025-1028
    # candidates.tsv: absent loci ordered by (-best bits, name) (run.py:1029-1042)
    _py_need_file "$out/sing_best.tsv"
    n=$(_p8_candidates_table "$out/sing_best.tsv" "$out/candidates_absent.txt" "$out/candidates.tsv") \
        || stage_die "${n:-20_candidates: cannot write $out/candidates.tsv}"
    manifest_set candidates "$n"                                       # run.py:1043
    say "  singleton loci above the floor and outside the catalogue: $n"                       # run.py:1044
    if [ "$n" = 0 ]; then emit_outputs "$out/candidates.tsv"; return 0; fi                    # run.py:1045-1046
    sh_run 20_candidates "" -- python3 "$SCRIPTS/07_post_catalogue/candidate_scan.py" "$out/candidates.tsv" \
        "$tables/genomes.tsv" "$tables/families.tsv" "$out" --flank "$(cfg orthology/flank_bp)"  # run.py:1048-1050
    sh_run 20_candidates "" -- python3 "$SCRIPTS/07_post_catalogue/candidate_flank5k.py" "$out" \
        "$tables/genomes.tsv" "$tables/families.tsv" --sat-bits "$(cfg expansion/min_bits)"     # run.py:1051-1053
    subx_dir                                                           # run.py:1055 (always called: writes tables/subx/)
    if [ -n "$SUBX" ] && [ -n "$(cfg_keys subfamily/members 2>/dev/null)" ]; then              # run.py:1056
        sh_run 20_candidates "" -- python3 "$SCRIPTS/07_post_catalogue/candidate_subsubfamily.py" "$out" \
            "$tables/genomes.tsv" "$tables/families.tsv" "$SUBX" "$tables/members.tsv" --already-doubled \
            --confident-delta "$(_p8_confident_delta)"                                         # run.py:1057-1061
    fi
    emit_outputs "$out/candidates.tsv" "$out/candidates_summary.json"                          # run.py:1062
}

# _p8_candidates_table BEST ABSENT CTAB   writes CTAB, prints the row count; on a malformed name prints the Python
# exception text and fails (rows before it are written, as Python's loop does).
#   best: read_tsv(BEST, 3) (run.py:144-154: text mode, skip blank-after-strip and '#' lines, >= 3 tab fields),
#         best[col1] = float(col3), a later row replaces an earlier one (run.py:1031-1032)
#   absent: stripped non-blank lines of ABSENT, [] when it is not a file (run.py:1033-1034)
#   order: key (-best.get(h, 0.0), h) (run.py:1035): bits descending, then the name bytewise ascending
#   row i: C<i> sp scf start end from sp|scf:start-end ; h.split('|')[0], [1]; reg.rsplit(':', 1); coords.split('-')
_p8_candidates_table() {
    local absent=""
    if [ -f "$2" ]; then absent=$2; fi
    BEST=$1 ABSENT=$absent CTAB=$3 "$GAWK" '
        function strip(s) { gsub(/^[ \t\n\v\f\r\034-\037]+|[ \t\n\v\f\r\034-\037]+$/, "", s); return s }
        function cmp(i1, v1, i2, v2,   b1, b2) {
            b1 = (v1 in best) ? best[v1] : 0; b2 = (v2 in best) ? best[v2] : 0
            if (b1 > b2) return -1; if (b1 < b2) return 1
            if ((v1 "") < (v2 "")) return -1; if ((v1 "") > (v2 "")) return 1
            return 0
        }
        BEGIN {
            RS = "\r\n|\r|\n"
            f = ENVIRON["BEST"]
            while ((getline line < f) > 0) {
                if (strip(line) == "" || substr(line, 1, 1) == "#") continue
                if (split(line, p, "\t") >= 3) best[p[1]] = p[3] + 0
            }
            close(f)
            f = ENVIRON["ABSENT"]; n = 0
            if (f != "") {                              # os.path.isfile, tested by the caller
                while ((getline line < f) > 0) { s = strip(line); if (s != "") h[++n] = s }
                close(f)
            }
            o = ENVIRON["CTAB"]; printf "" > o
            PROCINFO["sorted_in"] = "cmp"; i = 0
            for (k in h) {
                x = h[k]; i++
                np = split(x, parts, "|")
                if (np < 2) { close(o); printf "list index out of range"; exit 1 }
                reg = parts[2]
                c = match(reg, /:[^:]*$/)
                if (!c) { close(o); printf "not enough values to unpack (expected 2, got 1)"; exit 1 }
                scf = substr(reg, 1, RSTART - 1); coords = substr(reg, RSTART + 1)
                nc = split(coords, se, "-")
                if (nc != 2) { close(o); printf (nc < 2 ? "not enough values to unpack (expected 2, got %d)" : "too many values to unpack (expected 2, got %d)"), nc; exit 1 }   # Python >= 3.14 wording (3.13 omits ", got N")
                printf "C%d\t%s\t%s\t%s\t%s\n", i, parts[1], scf, se[1], se[2] > o
            }
            close(o)
            printf "%d", n
        }'
}
