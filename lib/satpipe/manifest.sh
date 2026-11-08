# manifest.sh - <out_dir>/manifest.json (PORT_SPEC 2.4; run.py:1117-1132, 1156). Sourced after common.sh.
# The manifest is held in a working file ($MAN_JSON, compact JSON in insertion order) and written as json.dump(indent=1)
# by manifest_write. Keys keep their position when overwritten (Python dict semantics). Not a TYPED-flat file: the DONE
# records under `stages` have file paths ('/') as keys.
#
#   manifest_init FILE              load FILE when it parses as a JSON object, else {} (run.py:1118-1122)
#   manifest_set KEY JSON           manifest[KEY] = JSON
#   manifest_set_flat KEY FLAT      manifest[KEY] = the TYPED-flat tree in FLAT (config: types preserved)
#   manifest_set_file KEY SUB FILE  manifest.setdefault(KEY, {})[SUB] = JSON document in FILE (DONE record, run.py:218)
#   manifest_append KEY STRING      manifest.setdefault(KEY, []).append(STRING) (runs)
#   manifest_write FILE             json.dump(manifest, f, indent=1)
#   manifest_summary KEY...         `k=str(v)` for every KEY present, joined by "; " (cli.py:86-87)

[ -n "${_SATPIPE_MANIFEST_SH:-}" ] && return 0
_SATPIPE_MANIFEST_SH=1

# _man OP   one gawk pass over $MAN_JSON; arguments in MAN_* environment variables; rewrites $MAN_JSON unless OP prints
_man() {
    local op=$1 tmp=$MAN_JSON.tmp
    MAN_OP=$op "$GAWK" -f "$PYJSON" -e '
        function cp(src, par, key,   c, i) {
            c = jt_set(par, key, JT[src], JV[src])
            for (i = 1; i <= JN[src]; i++) cp(src SUBSEP JK[src, i], c, JK[src, i])
            return c
        }
        function slurp(f,   s, l, rs) { rs = RS; RS = "^$"; s = ""; if ((getline s < f) <= 0) s = ""; close(f); RS = rs; return s }
        BEGIN {
            op = ENVIRON["MAN_OP"]; K = ENVIRON["MAN_KEY"]
            if (!json_parse(slurp(ENVIRON["MAN_JSON"]), "M") || JT["M"] != "map") { if ("M" in JT) jt_del("M"); jt_init("M", "map") }
            if (op == "init") { }
            else if (op == "set") {
                if (!json_parse(ENVIRON["MAN_VAL"], "N")) { print "manifest_set " K ": " JERR > "/dev/stderr"; exit 1 }
                cp("N", "M", K)
            } else if (op == "set_flat") {
                f = ENVIRON["MAN_VAL"]; jt_init("N", "map")
                while ((getline line < f) > 0) if (!flat_line("N", line)) { print JERR > "/dev/stderr"; exit 1 }
                close(f); cp("N", "M", K)
            } else if (op == "set_file") {
                if (!json_parse(slurp(ENVIRON["MAN_VAL"]), "N")) { print "manifest_set_file " K ": " JERR > "/dev/stderr"; exit 1 }
                if (!(("M" SUBSEP K) in JT) || JT["M" SUBSEP K] != "map") jt_set("M", K, "map", "")
                cp("N", "M" SUBSEP K, ENVIRON["MAN_SUB"])
            } else if (op == "append") {
                if (!(("M" SUBSEP K) in JT) || JT["M" SUBSEP K] != "list") jt_set("M", K, "list", "")
                L = "M" SUBSEP K; jt_set(L, JN[L], "str", ENVIRON["MAN_VAL"])
            } else if (op == "write") { printf "%s", json_indent1("M", 0); exit 0 }
            else if (op == "summary") {
                n = split(ENVIRON["MAN_VAL"], ks, " "); out = ""; c = 0
                for (i = 1; i <= n; i++) if (("M" SUBSEP ks[i]) in JT) out = out (c++ ? "; " : "") ks[i] "=" py_str("M" SUBSEP ks[i])
                printf "%s", out; exit 0
            }
            printf "%s", json_dumps_ordered("M")
        }' < /dev/null
}

manifest_init() {
    MAN_JSON=${MAN_JSON:-$(mktemp "${TMPDIR:-/tmp}/satpipe_manifest.XXXXXX")}
    if [ -f "$1" ]; then cp -- "$1" "$MAN_JSON"; else : > "$MAN_JSON"; fi
    _man_apply init
}
_man_apply() { local t=$MAN_JSON.tmp; MAN_JSON=$MAN_JSON _man "$1" > "$t" || { rm -f -- "$t"; return 1; }; mv -f -- "$t" "$MAN_JSON"; }
manifest_set()      { MAN_KEY=$1 MAN_VAL=$2 _man_apply set; }
manifest_set_flat() { MAN_KEY=$1 MAN_VAL=$2 _man_apply set_flat; }
manifest_set_file() { MAN_KEY=$1 MAN_SUB=$2 MAN_VAL=$3 _man_apply set_file; }
manifest_append()   { MAN_KEY=$1 MAN_VAL=$2 _man_apply append; }
manifest_write() {
    local t=$1.tmp.$$
    MAN_JSON=$MAN_JSON _man write > "$t" || { rm -f -- "$t"; return 1; }
    mv -f -- "$t" "$1"
}
manifest_summary() { local IFS=' '; MAN_JSON=$MAN_JSON MAN_VAL="$*" _man summary; }
