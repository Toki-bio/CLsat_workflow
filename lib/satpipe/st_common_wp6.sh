# st_common_wp6.sh - helpers shared by the stage files of WP6 (st_00_tables .. st_08_cross); sourced by them, not by
# stages_load. May be used by other stage files (source it with the same guard).
#
#   wp6_codes                 WP6_CODES=() genome codes in config order (ctx.codes, run.py:150-152)
#   wp6_families              WP6_FAMS=() family names in config order (ctx.families)
#   wp6_disc_env              WP6_DISC_ENV=(VAR=value ...) of _disc_env (run.py:473-477), ready for sh_run
#   wp6_py_list_repr STR...   Python repr() of a list of str, e.g. ['a', 'b'] (StageError texts, run.py:389, 401)
#   wp6_manifest_extend KEY JSONLIST   ctx.manifest.setdefault(KEY, []).extend(LIST) (run.py:450-452), via manifest_set
#   wp6_read_tsv FILE MINCOLS read_tsv (run.py:116-126): non-blank, non-'#' lines with >= MINCOLS tab fields, verbatim
# Child commands never read a loop's stdin: every list is read into an array before commands run.

[ -n "${_SATPIPE_ST_COMMON_WP6_SH:-}" ] && return 0
_SATPIPE_ST_COMMON_WP6_SH=1

wp6_codes()    { mapfile -t WP6_CODES < <(cfg_keys genomes); }
wp6_families() { mapfile -t WP6_FAMS < <(cfg_keys families); }

wp6_disc_env() {                                               # run.py:473-477, same keys as Python's ctx.env(**extra)
    WP6_DISC_ENV=(PERIOD_MIN="$(cfg discovery/period/0)" PERIOD_MAX="$(cfg discovery/period/1)"
                  MIN_COPIES="$(cfg discovery/min_copies)" TRF_PARAMS="$(cfg discovery/trf)"
                  SIMPLE_REPEAT_COV="$(cfg discovery/simple_repeat_cov)" STRAND="$(cfg discovery/strand)"
                  CLUSTER_ID="$(cfg discovery/cluster_id)" CROSS_CLUSTER_ID="$(cfg discovery/cross_cluster_id)")
}

wp6_py_list_repr() {
    if [ $# -eq 0 ]; then printf '[]'; return 0; fi
    printf '%s\0' "$@" | "$GAWK" -f "$PYJSON" -e 'BEGIN { RS = "\0" } { o = o (NR > 1 ? ", " : "") py_repr_str($0) }
        END { printf "[%s]", o }'
}

wp6_manifest_extend() {
    local cur
    cur=$(MAN_JSON=$MAN_JSON K=$1 V=$2 "$GAWK" -f "$PYJSON" -e '
        function cp(src, par, key,   c, i) {
            c = jt_set(par, key, JT[src], JV[src])
            for (i = 1; i <= JN[src]; i++) cp(src SUBSEP JK[src, i], c, JK[src, i])
            return c
        }
        BEGIN {
            f = ENVIRON["MAN_JSON"]; RS = "^$"; s = ""; if ((getline s < f) <= 0) s = ""; close(f)
            K = ENVIRON["K"]
            if (!json_parse(s, "M") || JT["M"] != "map") { if ("M" in JT) jt_del("M"); jt_init("M", "map") }
            if (!json_parse(ENVIRON["V"], "N") || JT["N"] != "list") { print "wp6_manifest_extend: bad list" > "/dev/stderr"; exit 1 }
            jt_init("L", "list")
            if (("M" SUBSEP K) in JT) {
                if (JT["M" SUBSEP K] != "list") { print "AttributeError: manifest " K " is not a list" > "/dev/stderr"; exit 1 }
                for (i = 1; i <= JN["M" SUBSEP K]; i++) cp("M" SUBSEP K SUBSEP JK["M" SUBSEP K, i], "L", JN["L"])
            }
            for (i = 1; i <= JN["N"]; i++) cp("N" SUBSEP JK["N", i], "L", JN["L"])
            printf "%s", json_dumps_ordered("L")
        }' < /dev/null) || stage_die "manifest: cannot extend $1"
    manifest_set "$1" "$cur"
}

wp6_read_tsv() {
    MC=$2 "$GAWK" -F '\t' 'BEGIN { mc = ENVIRON["MC"] + 0 }
        { l = $0 } /^[ \t\n\r\f\v]*$/ || /^#/ { next } NF >= mc { print l }' "$1"
}
