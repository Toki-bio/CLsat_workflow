# yaml_dump.awk - `satpipe config FILE`: the normalised config as block-style YAML (config.py:111 dump_config).
# Not byte-identical to yaml.safe_dump (PORT_SPEC 7); the contract is semantic: yaml.safe_load of this output equals
# load_config(FILE).data, types included (int 50 stays 50, float 50.0 stays 50.0, None, True/False, '' ...).
# Usage: LC_ALL=C gawk -f pyjson.awk -f yaml_dump.awk FLAT      (TYPED flat; `!warning` lines are ignored)
# Layout: mappings in flat order, `key: value`; lists one `- item` per line; a non-empty container inside a list
# starts on the line after its `-`; empty containers {} and []; strings plain only when they cannot be read back as
# anything else, else single-quoted ('' for '), or double-quoted with escapes when they hold control characters.

BEGIN { FS = "\t" }
{ sub(/\r$/, "") }
/^!/ { next }
{ if (!flat_load("Y", $0)) { print "yaml_dump.awk: " JERR > "/dev/stderr"; exit 2 } }
END {
    if (!("Y" in JT)) { print "{}"; exit 0 }
    if (JN["Y"] == 0) { print "{}"; exit 0 }
    emit_map("Y", "")
}

function yscalar(id,   t, v) {
    t = JT[id]; v = JV[id]
    if (t == "null") return "null"
    if (t == "bool") return (v == "True") ? "true" : "false"
    if (t == "int") return v
    if (t == "float") {
        if (v == "inf") return ".inf"; if (v == "-inf") return "-.inf"; if (v == "nan") return ".nan"
        if (index(v, ".") == 0 && index(v, "e")) sub(/e/, ".0e", v)   # PyYAML reads 1e-05 as a string, 1.0e-05 as a float
        return v
    }
    return ystr(v)
}

function ystr(s,   out, i, n, c, b) {
    if (s ~ /^[A-Za-z_\/][A-Za-z0-9_.\/+-]*( [A-Za-z0-9_.\/+-]+)*$/ &&
        s !~ /^(yes|Yes|YES|no|No|NO|true|True|TRUE|false|False|FALSE|on|On|ON|off|Off|OFF|null|Null|NULL)$/)
        return s
    if (s !~ /[\001-\037\177]/ && s !~ /\302[\200-\237]|\342\200[\250\251]|\357\273\277/) { gsub(/'/, "''", s); return "'" s "'" }
    out = "\""; n = length(s)
    for (i = 1; i <= n; i++) {
        c = substr(s, i, 1); b = _PJ_ORD[c]
        if (c == "\"") out = out "\\\""
        else if (c == "\\") out = out "\\\\"
        else if (c == "\n") out = out "\\n"
        else if (c == "\t") out = out "\\t"
        else if (c == "\r") out = out "\\r"
        else if (b < 32 || b == 127) out = out sprintf("\\x%02x", b)
        else out = out c
    }
    return out "\""
}

function emit_child(c, ind, lead,   t) {
    t = JT[c]
    if (t == "map" || t == "list") {
        if (JN[c] == 0) { print lead " " (t == "map" ? "{}" : "[]"); return }
        print lead
        if (t == "map") emit_map(c, ind "  "); else emit_list(c, ind "  ")
        return
    }
    print lead " " yscalar(c)
}

function emit_map(id, ind,   i, k) {
    for (i = 1; i <= JN[id]; i++) { k = JK[id, i]; emit_child(id SUBSEP k, ind, ind ystr(k) ":") }
}

function emit_list(id, ind,   i) {
    for (i = 1; i <= JN[id]; i++) emit_child(id SUBSEP JK[id, i], ind, ind "-")
}

# one TYPED-flat line into the tree under root; unlike pyjson flat_line, path components are flat-unescaped too
function flat_load(root, line,   f, n, a, i, id) {
    n = split(line, f, "\t")
    if (n < 3 || f[2] !~ /^(map|list|int|float|str|path|bool|null)$/) { JERR = "malformed flat line: " line; return 0 }
    if (!(root in JT)) jt_init(root, "map")
    n = split(f[1], a, "/"); id = root
    for (i = 1; i < n; i++) { id = id SUBSEP flat_unesc(a[i]); if (!(id in JT)) { JERR = "flat line before its parent: " line; return 0 } }
    jt_set(id, flat_unesc(a[n]), f[2], (f[2] == "map" || f[2] == "list") ? "" : flat_unesc(f[3]))
    return 1
}
