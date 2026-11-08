# flat.sh - TYPED-flat reader/writer helpers (PORT_SPEC 2.1). Sourced after common.sh.
# Line format: path<TAB>type<TAB>text<TAB>origin; text escaped as in pyjson.awk (flat_esc): \\ \t \n \r.
# Root is an implicit map (no row). Container rows carry the item count.
# Arguments reach gawk through ENVIRON, never -v (which would interpret backslashes).

# flat_emit PATH TYPE TEXT [ORIGIN]      print one line (escapes TEXT); for a container TEXT is the item count
flat_emit() {
    local t=$3
    t=${t//\\/\\\\}; t=${t//$'\t'/\\t}; t=${t//$'\n'/\\n}; t=${t//$'\r'/\\r}
    printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$t" "${4:--}"
}

# flat_get FILE PATH        unescaped text of a scalar leaf (Python str() of the value); exit 1 when absent
flat_get() {
    P=$2 pyjson_gawk 'BEGIN { FS = "\t"; P = ENVIRON["P"] }
        $1 == P { sub(/\r$/, ""); printf "%s", flat_unesc($3); f = 1; exit } END { exit !f }' "$1"
}

# flat_type FILE PATH       type column of PATH; exit 1 when absent
flat_type() {
    P=$2 "$GAWK" -F '\t' 'BEGIN { P = ENVIRON["P"] } $1 == P { print $2; f = 1; exit } END { exit !f }' "$1"
}

# flat_keys FILE PATH       keys of the direct children of PATH ("" = root), in document order, one per line
flat_keys() {
    P=$2 "$GAWK" -F '\t' 'BEGIN { P = ENVIRON["P"]; pre = (P == "") ? "" : P "/"; n = length(pre) }
        length($1) > n && substr($1, 1, n) == pre && index(substr($1, n + 1), "/") == 0 { print substr($1, n + 1) }' "$1"
}

# flat_subtree FILE SRC DST   the rows of SRC and everything under it, re-rooted at DST (cfg_subtree uses this)
flat_subtree() {
    S=$2 D=$3 "$GAWK" -F '\t' -v OFS='\t' 'BEGIN { S = ENVIRON["S"]; D = ENVIRON["D"] }
        $1 == S { $1 = D; print; next }
        substr($1, 1, length(S) + 1) == S "/" { $1 = D substr($1, length(S) + 1); print }' "$1"
}

# flat_json FILE            json.dumps(tree, sort_keys=True) of the whole file, no trailing newline
flat_json() {
    pyjson_gawk '{ if (!flat_line("R", $0)) { print JERR > "/dev/stderr"; bad = 1; exit 1 } }
        END { if (bad) exit 1; if (!("R" in JT)) jt_init("R", "map"); printf "%s", json_dumps("R") }' "$1"
}

# flat_py_repr FILE PATH    Python repr() of the node at PATH ("" = root)
flat_py_repr() {
    P=$2 pyjson_gawk '{ if (!flat_line("R", $0)) { print JERR > "/dev/stderr"; bad = 1; exit 1 } }
        END { if (bad) exit 1; P = ENVIRON["P"]; if (!("R" in JT)) jt_init("R", "map")
              id = "R"; if (P != "") { n = split(P, a, "/"); for (i = 1; i <= n; i++) id = id SUBSEP a[i] }
              if (!(id in JT)) exit 1; printf "%s", py_repr(id) }' "$1"
}

# json_to_flat FILE [KEY]   TYPED flat of a JSON document (or of its top-level KEY, e.g. DONE "params").
#                           ints keep their digits, floats become Python repr, so json_dumps(flat) == json.dumps(json.loads)
json_to_flat() {
    K=${2:-} pyjson_gawk 'BEGIN { RS = "^$"; K = ENVIRON["K"] }
        { got = 1; if (!json_parse($0, "J")) { print "json_to_flat: " JERR > "/dev/stderr"; bad = 1; exit 1 } }
        END { if (bad) exit 1; if (!got) { print "json_to_flat: empty input" > "/dev/stderr"; exit 1 }
              id = "J"; if (K != "") { id = "J" SUBSEP K; if (!(id in JT)) { print "json_to_flat: no key " K > "/dev/stderr"; exit 1 } }
              if (JT[id] != "map") { print "json_to_flat: not an object" > "/dev/stderr"; exit 1 }
              flat_write(id, "", "-") }' "$1"
}
