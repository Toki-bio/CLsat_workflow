# validate.awk - RAW flat (yaml.awk) + schema.tsv -> TYPED flat with defaults filled and paths resolved
# (PORT_SPEC 2.1, section 3; port of schema.py validate(), schema.py:248-340, and config.py:96-102).
#
# Usage:  LC_ALL=C BASE_DIR=<dir of the absolute config path> \
#         gawk -f pyjson.awk -f normpath.awk -f validate.awk schema.tsv RAWFILE      (RAWFILE may be - = stdin)
#
# Output (stdout):
#   valid:   TYPED flat lines path<TAB>type<TAB>text<TAB>origin in validate()'s order (block fields in schema order,
#            maps and lists in config order); exit 0
#            type   int float str path bool null map list (enum -> str, block -> map; a float field absent from the YAML
#                   keeps its default OBJECT, so `expansion/min_bits` is int 50 then, float 50.0 when written)
#            text   Python str() of the value, flat-escaped (\\ \t \n \r); containers: item count
#            origin file     = the key is in the YAML with a non-null value (containers too)
#                   default  = absent or null, value from the schema default (or None, or {} for a block)
#                   derived  = absent or null, filled at run time (schema.tsv absent_origin; anchor/family)
#            path components are flat-escaped too (only pathological quoted keys can need it)
#   errors:  one line per error, `!error<TAB><message>` (message flat-escaped), in validate()'s order; exit 1.
#            Messages are the exact strings of schema.py:252-335 (dotted paths, [i] for list items, Python repr()).
#   a RAW `!error` line (YAML subset error, top level not a mapping) is passed through unchanged; exit 1.
# Known approximations (pathological input only): mapping keys that are not strings in the YAML (1: x, true: x) are
# quoted as strings inside a repr()'d mapping; Python's \d also matches non-ASCII digits, [0-9] does not.

BEGIN { FS = "\t"; NROW = 0; RAWSEEN = 0; NERR = 0 }

# ---------------------------------------------------------------- schema.tsv (first file)
FNR == NR {
    sub(/\r$/, "")
    if (/^#/) { sub(/^#/, ""); for (i = 1; i <= NF; i++) SC[$i] = i; next }
    sp = $SC["path"]
    S_type[sp] = $SC["type"]; S_req[sp] = ($SC["required"] == "True")
    S_dtype[sp] = $SC["default_type"]; S_djson[sp] = $SC["default_json"]; S_dstr[sp] = $SC["default_str"]
    S_min[sp] = $SC["min"]; S_max[sp] = $SC["max"]; S_excl[sp] = ($SC["exclusive_min"] == "True")
    S_pat[sp] = $SC["pattern"]; S_pere[sp] = $SC["pattern_ere"]
    S_kpat[sp] = $SC["key_pattern"]; S_kere[sp] = $SC["key_pattern_ere"]
    S_choices[sp] = $SC["choices"]; S_len[sp] = $SC["length"]; S_minit[sp] = $SC["min_items"]
    S_aorig[sp] = $SC["absent_origin"]
    if (sp != ".") {                                   # register as a field of its parent block
        par = sp; k = sp
        if (match(sp, /\/[^\/]*$/)) { par = substr(sp, 1, RSTART - 1); k = substr(sp, RSTART + 1) } else par = "."
        if (k != "<key>" && k != "<item>") { S_nf[par]++; S_field[par, S_nf[par]] = k; S_isfield[par, k] = 1 }
    }
    next
}

# ---------------------------------------------------------------- RAW flat (second file)
{
    sub(/\r$/, "")
    if (substr($0, 1, 7) == "!error\t") { print; PASSERR = 1; exit 1 }
    if (!RAWSEEN) { jt_init("R", "map"); RAWSEEN = 1 }
    n = split($1, comp, "/")
    id = "R"
    for (i = 1; i < n; i++) id = id SUBSEP raw_key(comp[i])
    jt_set(id, raw_key(comp[n]), $2, ($2 == "map" || $2 == "list") ? "" : flat_unesc($3))
}

END {
    if (PASSERR) exit 1
    if (!RAWSEEN) jt_init("R", "map")                 # empty document: raw = {} (config.py:96-97)
    if (!(("." in S_type))) { print "validate.awk: schema.tsv not read" > "/dev/stderr"; exit 2 }
    jt_init("O", "map")
    val_block_children("R", ".", "", "O")
    if (NERR) { for (i = 1; i <= NERR; i++) print "!error\t" flat_esc(ERR[i]); exit 1 }
    write_flat("O", "")
    exit 0
}

# RAW key component: flat escapes plus '/' as \x2f (yaml.awk header)
function raw_key(c,   out, i, n, ch) {
    if (index(c, "\\") == 0) return c
    out = ""; n = length(c)
    for (i = 1; i <= n; i++) {
        ch = substr(c, i, 1)
        if (ch == "\\" && substr(c, i + 1, 3) == "x2f") { out = out "/"; i += 3 }
        else if (ch == "\\" && i < n) { ch = substr(c, ++i, 1); out = out (ch == "t" ? "\t" : ch == "n" ? "\n" : ch == "r" ? "\r" : ch) }
        else out = out ch
    }
    return out
}

function err(msg) { ERR[++NERR] = msg }

function pytype(id,   t) {
    t = JT[id]
    return (t == "map") ? "dict" : (t == "null") ? "NoneType" : (t == "path") ? "str" : t
}

# re.match(pattern, s) for the schema's ^...$ patterns: Python's $ also matches before one trailing newline
function py_match(s, ere) {
    if (s ~ ere) return 1
    return (substr(s, length(s)) == "\n" && substr(s, 1, length(s) - 1) ~ ere)
}

# str.strip() is empty (ASCII and the common Unicode white space; controls are refused by the YAML reader)
function blank(s) {
    gsub(/[ \t\n\r\013\014\034-\037]/, "", s)
    gsub(/\302\240|\341\232\200|\342\200[\200-\212\250\251\257]|\342\201\237|\343\200\200/, "", s)
    return s == ""
}

function choice_ok(sp, v,   c, n, a, i) {
    c = S_choices[sp]; c = substr(c, 2, length(c) - 2)
    n = split(c, a, ", ")
    for (i = 1; i <= n; i++) if (substr(a[i], 2, length(a[i]) - 2) == v) return 1
    return 0
}

# copy a raw subtree under a new output node (values returned unchanged by validate on a type error)
function copy_node(rid, op, ok, origin,   oid, i, k) {
    oid = jt_set(op, ok, JT[rid], JV[rid]); ORIG[oid] = origin
    for (i = 1; i <= JN[rid]; i++) { k = JK[rid, i]; copy_node(rid SUBSEP k, oid, k, origin) }
    return oid
}

function check_number(text, sp, ep,   lo, hi, v) {
    lo = S_min[sp]; hi = S_max[sp]; v = text + 0
    if (lo != "") {
        if (S_excl[sp] && !(v > lo + 0)) err(ep ": must be > " lo ", got " text)
        else if (!S_excl[sp] && v < lo + 0) err(ep ": must be >= " lo ", got " text)
    }
    if (hi != "" && v > hi + 0) err(ep ": must be <= " hi ", got " text)
}

# default object of a field (schema.py:270): typed by the default's own Python type
function put_default(sp, op, ok,   oid, origin) {
    origin = S_aorig[sp]
    if (S_dtype[sp] == "list" || S_dtype[sp] == "map") {
        if (!json_parse(S_djson[sp], "D")) { print "validate.awk: bad default for " sp ": " JERR > "/dev/stderr"; exit 2 }
        oid = copy_node("D", op, ok, origin); jt_del("D")
        return oid
    }
    oid = jt_set(op, ok, S_dtype[sp], S_dstr[sp]); ORIG[oid] = origin
    return oid
}

# validate(node, spec, path) -> output child `ok` of `op`.  rid = raw node id or "" when the key is absent.
function val(rid, sp, ep, op, ok,   t, rt, oid, i, k, n, kid, v) {
    t = S_type[sp]
    if (rid == "" || JT[rid] == "null") {
        if (S_req[sp]) { err(ep ": required"); oid = jt_set(op, ok, "null", "None"); ORIG[oid] = S_aorig[sp]; return }
        if (t == "block") {
            oid = jt_set(op, ok, "map", ""); ORIG[oid] = S_aorig[sp]
            jt_init("E" SUBSEP oid, "map")
            val_block_children("E" SUBSEP oid, sp, ep, oid)
            return
        }
        put_default(sp, op, ok)
        return
    }
    rt = JT[rid]
    if (t == "str") {
        if (rt != "str") { err(ep ": expected a string, got " pytype(rid)); copy_node(rid, op, ok, "file"); return }
        if (S_pat[sp] != "" && !py_match(JV[rid], S_pere[sp])) err(ep ": '" JV[rid] "' does not match " S_pat[sp])
        oid = jt_set(op, ok, "str", JV[rid]); ORIG[oid] = "file"; return
    }
    if (t == "enum") {
        if (!(rt == "str" && choice_ok(sp, JV[rid]))) err(ep ": must be one of " S_choices[sp] ", got " py_repr(rid))
        copy_node(rid, op, ok, "file"); return
    }
    if (t == "bool") {
        if (rt != "bool") err(ep ": expected true/false, got " py_repr(rid))
        copy_node(rid, op, ok, "file"); return
    }
    if (t == "int") {
        if (rt != "int") { err(ep ": expected an integer, got " py_repr(rid)); copy_node(rid, op, ok, "file"); return }
        check_number(JV[rid], sp, ep)
        oid = jt_set(op, ok, "int", JV[rid]); ORIG[oid] = "file"; return
    }
    if (t == "float") {
        if (rt != "int" && rt != "float") { err(ep ": expected a number, got " py_repr(rid)); copy_node(rid, op, ok, "file"); return }
        check_number(JV[rid], sp, ep)
        oid = jt_set(op, ok, "float", pyrepr_float(JV[rid] + 0)); ORIG[oid] = "file"; return
    }
    if (t == "path") {
        if (rt != "str" || blank(JV[rid])) { err(ep ": expected a file path, got " py_repr(rid)); copy_node(rid, op, ok, "file"); return }
        oid = jt_set(op, ok, "path", satpipe_norm_path(JV[rid], ENVIRON["BASE_DIR"])); ORIG[oid] = "file"; return
    }
    if (t == "list") {
        if (rt != "list") { err(ep ": expected a list, got " pytype(rid)); copy_node(rid, op, ok, "file"); return }
        n = JN[rid]
        if (S_len[sp] != "" && n != S_len[sp] + 0) err(ep ": expected exactly " S_len[sp] " items, got " n)
        if (S_minit[sp] != "" && n < S_minit[sp] + 0) err(ep ": expected at least " S_minit[sp] " items, got " n)
        oid = jt_set(op, ok, "list", ""); ORIG[oid] = "file"
        for (i = 1; i <= n; i++) { k = JK[rid, i]; val(rid SUBSEP k, sp "/<item>", ep "[" (i - 1) "]", oid, i - 1) }
        return
    }
    if (t == "map") {
        if (rt != "map") { err(ep ": expected a mapping, got " pytype(rid)); copy_node(rid, op, ok, "file"); return }
        n = JN[rid]
        if (S_minit[sp] != "" && n < S_minit[sp] + 0) err(ep ": expected at least " S_minit[sp] " entries, got " n)
        oid = jt_set(op, ok, "map", ""); ORIG[oid] = "file"
        for (i = 1; i <= n; i++) {
            k = JK[rid, i]
            if (S_kere[sp] != "" && !py_match(k, S_kere[sp])) err(ep "." k ": key does not match " S_kpat[sp])
            val(rid SUBSEP k, sp "/<key>", ep "." k, oid, k)
        }
        return
    }
    if (t == "block") {
        if (rt != "map") { err(ep ": expected a mapping, got " pytype(rid)); copy_node(rid, op, ok, "file"); return }
        oid = jt_set(op, ok, "map", ""); ORIG[oid] = "file"
        val_block_children(rid, sp, ep, oid)
        return
    }
    print "validate.awk: bad spec type " t " at " sp > "/dev/stderr"; exit 2
}

# the block part of validate() for raw mapping rid (schema.py:332-339): unknown keys first, then fields in order
function val_block_children(rid, sp, ep, oid,   i, k, allowed, pre) {
    pre = (ep == "") ? "" : ep "."
    allowed = ""
    for (i = 1; i <= S_nf[sp]; i++) allowed = allowed (i > 1 ? ", " : "") S_field[sp, i]
    for (i = 1; i <= JN[rid]; i++) {
        k = JK[rid, i]
        if (!((sp, k) in S_isfield)) err(pre k ": unknown key (allowed: " allowed ")")
    }
    for (i = 1; i <= S_nf[sp]; i++) {
        k = S_field[sp, i]
        val(((rid SUBSEP k) in JT) ? rid SUBSEP k : "", (sp == "." ? "" : sp "/") k, pre k, oid, k)
    }
}

function write_flat(id, prefix,   i, k, c, p) {
    for (i = 1; i <= JN[id]; i++) {
        k = JK[id, i]; c = id SUBSEP k
        p = (prefix == "") ? flat_esc(k) : prefix "/" flat_esc(k)
        if (JT[c] == "map" || JT[c] == "list") { print p "\t" JT[c] "\t" JN[c] "\t" ORIG[c]; write_flat(c, p) }
        else print p "\t" JT[c] "\t" flat_esc(JV[c]) "\t" ORIG[c]
    }
}
