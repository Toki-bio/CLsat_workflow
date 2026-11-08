# pyjson.awk - typed-tree utilities that reproduce Python's text forms byte for byte (PORT_SPEC 2.1, 2.2).
# Library: load with `gawk -f pyjson.awk -e '...'`. Needs gawk and LC_ALL=C (bytes, bytewise ordering).
#
# Tree model (several trees may coexist under different root ids):
#   JT[id]       type: map list int float str path bool null
#   JV[id]       scalar text = Python str() of the value (float text is already Python repr)
#   JN[id]       number of children of a container
#   JK[id, i]    key of child i (1-based, insertion order); list children have keys 0..n-1
#   child id     = id SUBSEP key
# TYPED flat line (config.flat.tsv, PLAN_PARAMS): path<TAB>type<TAB>text<TAB>origin; path components joined by '/',
#   list items by index; container rows carry the item count in text; root is an implicit map (no row).
#   text escaping in the flat file: '\' -> '\\', TAB -> '\t', LF -> '\n', CR -> '\r' (flat_esc / flat_unesc).
#
# API
#   jt_init(id, type)                    make an empty node
#   jt_set(parent, key, type, text)      add/overwrite a child (overwrite keeps the position, like a Python dict); returns child id
#   jt_del(id)                           delete a subtree (node itself included)
#   flat_line(root, line)                add one TYPED-flat line under root; returns 0 on a malformed line (JERR set)
#   flat_write(id, prefix, origin)       print the subtree of id as TYPED-flat lines (prefix = path of id, "" for root)
#   json_parse(text, root)               JSON text -> tree; 1 ok, 0 error (JERR set). ints keep their digits, floats -> repr
#   json_dumps(id)                       = json.dumps(x, sort_keys=True)            (ensure_ascii)
#   json_dumps_ordered(id)               = json.dumps(x)  (insertion order)
#   json_indent1(id)                     = json.dumps(x, indent=1)  (insertion order, as json.dump(rec, f, indent=1))
#   json_str(s)                          JSON string literal with ensure_ascii escapes
#   py_str(id), py_repr(id), py_repr_str(s), pyrepr_float(x)
#   flat_esc(s), flat_unesc(s)

BEGIN {
    for (_pj_i = 0; _pj_i < 256; _pj_i++) { _PJ_CHR[_pj_i] = sprintf("%c", _pj_i); _PJ_ORD[_PJ_CHR[_pj_i]] = _pj_i }
    if (length(_PJ_CHR[200]) != 1) { print "pyjson.awk: needs LC_ALL=C (byte characters)" > "/dev/stderr"; exit 2 }
    _PJ_JESC["\""] = "\\\""; _PJ_JESC["\\"] = "\\\\"; _PJ_JESC["\n"] = "\\n"; _PJ_JESC["\r"] = "\\r"
    _PJ_JESC["\t"] = "\\t"; _PJ_JESC["\b"] = "\\b"; _PJ_JESC["\f"] = "\\f"
}

# ---------------------------------------------------------------- tree

function jt_init(id, type) { JT[id] = type; JV[id] = ""; JN[id] = 0 }

function jt_set(parent, key, type, text,   c) {
    c = parent SUBSEP key
    if (c in JT) { _pj_delkids(c) }
    else { JK[parent, ++JN[parent]] = key }
    JT[c] = type; JV[c] = text; JN[c] = 0
    return c
}

function _pj_delkids(id,   i) {
    for (i = 1; i <= JN[id]; i++) { jt_del(id SUBSEP JK[id, i]); delete JK[id, i] }
    JN[id] = 0
}

function jt_del(id) { _pj_delkids(id); delete JT[id]; delete JV[id]; delete JN[id] }

# ---------------------------------------------------------------- TYPED flat

function flat_esc(s) {
    if (s !~ /[\\\t\n\r]/) return s
    gsub(/\\/, "\\\\", s); gsub(/\t/, "\\t", s); gsub(/\n/, "\\n", s); gsub(/\r/, "\\r", s)
    return s
}

function flat_unesc(s,   out, i, c, n) {
    if (index(s, "\\") == 0) return s
    out = ""; n = length(s)
    for (i = 1; i <= n; i++) {
        c = substr(s, i, 1)
        if (c == "\\" && i < n) {
            c = substr(s, ++i, 1)
            out = out (c == "t" ? "\t" : c == "n" ? "\n" : c == "r" ? "\r" : c)
        } else out = out c
    }
    return out
}

function _pj_ensure(root, path,   parts, n, i, id, k) {
    id = root
    if (path == "") return id
    n = split(path, parts, "/")
    for (i = 1; i <= n; i++) {
        k = parts[i]
        if (!((id SUBSEP k) in JT)) jt_set(id, k, "map", "")
        id = id SUBSEP k
    }
    return id
}

function flat_line(root, line,   f, n, path, par, key, p, id) {
    sub(/\r$/, "", line)
    n = split(line, f, "\t")
    if (n < 3 || f[1] == "" || f[2] !~ /^(map|list|int|float|str|path|bool|null)$/) { JERR = "malformed flat line: " line; return 0 }
    path = f[1]
    p = match(path, /\/[^\/]*$/)
    if (p) { par = substr(path, 1, p - 1); key = substr(path, p + 1) } else { par = ""; key = path }
    if (!(root in JT)) jt_init(root, "map")
    id = _pj_ensure(root, par)
    if (JT[id] != "map" && JT[id] != "list") { JERR = "flat line under a scalar: " line; return 0 }
    if (((id SUBSEP key) in JT) && f[2] ~ /^(map|list)$/ && JT[id SUBSEP key] == f[2]) {
        return 1                                   # container row after an implicitly created parent: keep children
    }
    jt_set(id, key, f[2], (f[2] ~ /^(map|list)$/) ? "" : flat_unesc(f[3]))
    return 1
}

function flat_write(id, prefix, origin,   i, k, c, p) {
    if (origin == "") origin = "-"
    for (i = 1; i <= JN[id]; i++) {
        k = JK[id, i]; c = id SUBSEP k
        if (index(k, "/") || k ~ /[\t\n\r]/ || k == "") { print "pyjson.awk: key not representable in a flat path: " k > "/dev/stderr"; exit 2 }
        p = (prefix == "") ? k : prefix "/" k
        if (JT[c] == "map" || JT[c] == "list") { print p "\t" JT[c] "\t" JN[c] "\t" origin; flat_write(c, p, origin) }
        else print p "\t" JT[c] "\t" flat_esc(JV[c]) "\t" origin
    }
}

# ---------------------------------------------------------------- JSON parse

function json_parse(text, root) {
    _PJ_S = text; _PJ_P = 1; _PJ_L = length(text); JERR = ""
    if (root in JT) jt_del(root)
    jt_init(root, "null")
    if (!_pj_value(root)) return 0
    _pj_ws()
    if (_PJ_P <= _PJ_L) { JERR = "trailing data at " _PJ_P; return 0 }
    return 1
}

function _pj_ws() { while (_PJ_P <= _PJ_L && index(" \t\n\r", substr(_PJ_S, _PJ_P, 1))) _PJ_P++ }

function _pj_value(id,   c, rest, s, k, n, child) {
    _pj_ws()
    if (_PJ_P > _PJ_L) { JERR = "unexpected end"; return 0 }
    c = substr(_PJ_S, _PJ_P, 1)
    if (c == "{") {
        JT[id] = "map"; JV[id] = ""; JN[id] = 0; _PJ_P++; _pj_ws()
        if (substr(_PJ_S, _PJ_P, 1) == "}") { _PJ_P++; return 1 }
        while (1) {
            _pj_ws()
            if (substr(_PJ_S, _PJ_P, 1) != "\"") { JERR = "expected key at " _PJ_P; return 0 }
            if (!_pj_string()) return 0
            k = _PJ_STR; _pj_ws()
            if (substr(_PJ_S, _PJ_P, 1) != ":") { JERR = "expected ':' at " _PJ_P; return 0 }
            _PJ_P++
            child = jt_set(id, k, "null", "")
            if (!_pj_value(child)) return 0
            _pj_ws(); c = substr(_PJ_S, _PJ_P++, 1)
            if (c == "}") return 1
            if (c != ",") { JERR = "expected ',' or '}' at " (_PJ_P - 1); return 0 }
        }
    }
    if (c == "[") {
        JT[id] = "list"; JV[id] = ""; JN[id] = 0; _PJ_P++; _pj_ws()
        if (substr(_PJ_S, _PJ_P, 1) == "]") { _PJ_P++; return 1 }
        n = 0
        while (1) {
            child = jt_set(id, n++, "null", "")
            if (!_pj_value(child)) return 0
            _pj_ws(); c = substr(_PJ_S, _PJ_P++, 1)
            if (c == "]") return 1
            if (c != ",") { JERR = "expected ',' or ']' at " (_PJ_P - 1); return 0 }
        }
    }
    if (c == "\"") { if (!_pj_string()) return 0; JT[id] = "str"; JV[id] = _PJ_STR; return 1 }
    rest = substr(_PJ_S, _PJ_P, 400)
    if (match(rest, /^-?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][-+]?[0-9]+)?/)) {
        s = substr(rest, 1, RLENGTH); _PJ_P += RLENGTH
        if (s ~ /^-?[0-9]+$/) { if (s ~ /^-0+$/) s = "0"; JT[id] = "int"; JV[id] = s }
        else { JT[id] = "float"; JV[id] = pyrepr_float(s + 0) }
        return 1
    }
    if (substr(rest, 1, 4) == "true")  { JT[id] = "bool"; JV[id] = "True";  _PJ_P += 4; return 1 }
    if (substr(rest, 1, 5) == "false") { JT[id] = "bool"; JV[id] = "False"; _PJ_P += 5; return 1 }
    if (substr(rest, 1, 4) == "null")  { JT[id] = "null"; JV[id] = "None";  _PJ_P += 4; return 1 }
    JERR = "unexpected character at " _PJ_P; return 0
}

function _pj_string(   out, c, h, cp, lo, seg) {
    _PJ_P++; out = ""
    while (1) {
        if (match(substr(_PJ_S, _PJ_P), /^[^"\\]+/)) { out = out substr(_PJ_S, _PJ_P, RLENGTH); _PJ_P += RLENGTH }
        if (_PJ_P > _PJ_L) { JERR = "unterminated string"; return 0 }
        c = substr(_PJ_S, _PJ_P, 1)
        if (c == "\"") { _PJ_P++; _PJ_STR = out; return 1 }
        c = substr(_PJ_S, _PJ_P + 1, 1); _PJ_P += 2
        if (c == "n") out = out "\n"; else if (c == "t") out = out "\t"; else if (c == "r") out = out "\r"
        else if (c == "b") out = out "\b"; else if (c == "f") out = out "\f"
        else if (c == "\"" || c == "\\" || c == "/") out = out c
        else if (c == "u") {
            h = substr(_PJ_S, _PJ_P, 4)
            if (h !~ /^[0-9A-Fa-f][0-9A-Fa-f][0-9A-Fa-f][0-9A-Fa-f]$/) { JERR = "bad \\u escape"; return 0 }
            cp = strtonum("0x" h); _PJ_P += 4
            if (cp >= 55296 && cp <= 56319 && substr(_PJ_S, _PJ_P, 2) == "\\u") {
                seg = substr(_PJ_S, _PJ_P + 2, 4)
                if (seg ~ /^[0-9A-Fa-f][0-9A-Fa-f][0-9A-Fa-f][0-9A-Fa-f]$/) {
                    lo = strtonum("0x" seg)
                    if (lo >= 56320 && lo <= 57343) { cp = 65536 + (cp - 55296) * 1024 + (lo - 56320); _PJ_P += 6 }
                }
            }
            out = out utf8_chr(cp)
        } else { JERR = "bad escape \\" c; return 0 }
    }
}

function utf8_chr(cp) {
    if (cp < 128) return _PJ_CHR[cp]
    if (cp < 2048) return _PJ_CHR[192 + int(cp / 64)] _PJ_CHR[128 + cp % 64]
    if (cp < 65536) return _PJ_CHR[224 + int(cp / 4096)] _PJ_CHR[128 + int(cp / 64) % 64] _PJ_CHR[128 + cp % 64]
    return _PJ_CHR[240 + int(cp / 262144)] _PJ_CHR[128 + int(cp / 4096) % 64] _PJ_CHR[128 + int(cp / 64) % 64] _PJ_CHR[128 + cp % 64]
}

# decode the UTF-8 sequence starting at byte i of s: sets _PJ_CP (code point) and returns its byte length
# (invalid bytes are returned as single code points 0x80-0xFF, i.e. Latin-1, never seen in valid str values)
function _pj_utf8_at(s, i,   b, n, k, cp, c) {
    b = _PJ_ORD[substr(s, i, 1)]
    if (b < 128) { _PJ_CP = b; return 1 }
    if (b >= 240 && b < 248) { n = 4; cp = b - 240 } else if (b >= 224) { n = 3; cp = b - 224 } else if (b >= 192) { n = 2; cp = b - 192 } else { _PJ_CP = b; return 1 }
    for (k = 1; k < n; k++) {
        c = _PJ_ORD[substr(s, i + k, 1)]
        if (c < 128 || c >= 192) { _PJ_CP = b; return 1 }
        cp = cp * 64 + (c - 128)
    }
    _PJ_CP = cp; return n
}

# ---------------------------------------------------------------- JSON output

function json_str(s,   out, i, n, c, b, len, cp, hi) {
    if (s !~ /[^ !#-\[\]-~]/) return "\"" s "\""
    out = "\""; n = length(s); i = 1
    while (i <= n) {
        c = substr(s, i, 1); b = _PJ_ORD[c]
        if (b >= 32 && b < 127 && c != "\"" && c != "\\") { out = out c; i++; continue }
        if (c in _PJ_JESC) { out = out _PJ_JESC[c]; i++; continue }
        if (b < 128) { out = out sprintf("\\u%04x", b); i++; continue }
        len = _pj_utf8_at(s, i); cp = _PJ_CP; i += len
        if (cp >= 65536) {
            cp -= 65536; hi = 55296 + int(cp / 1024)
            out = out sprintf("\\u%04x\\u%04x", hi, 56320 + cp % 1024)
        } else out = out sprintf("\\u%04x", cp)
    }
    return out "\""
}

function _pj_scalar_json(id,   t) {
    t = JT[id]
    if (t == "int" || t == "float") return JV[id]
    if (t == "bool") return (JV[id] == "True") ? "true" : "false"
    if (t == "null") return "null"
    return json_str(JV[id])
}

function json_dumps(id,   t, i, n, out, S, keys) {
    t = JT[id]
    if (t == "map") {
        if (JN[id] == 0) return "{}"
        for (i = 1; i <= JN[id]; i++) S[JK[id, i]] = 1
        n = asorti(S, keys, "@ind_str_asc")
        out = "{"
        for (i = 1; i <= n; i++) out = out (i > 1 ? ", " : "") json_str(keys[i]) ": " json_dumps(id SUBSEP keys[i])
        return out "}"
    }
    if (t == "list") {
        out = "["
        for (i = 1; i <= JN[id]; i++) out = out (i > 1 ? ", " : "") json_dumps(id SUBSEP JK[id, i])
        return out "]"
    }
    return _pj_scalar_json(id)
}

function json_dumps_ordered(id,   t, i, out, k) {
    t = JT[id]
    if (t == "map") {
        out = "{"
        for (i = 1; i <= JN[id]; i++) { k = JK[id, i]; out = out (i > 1 ? ", " : "") json_str(k) ": " json_dumps_ordered(id SUBSEP k) }
        return out "}"
    }
    if (t == "list") {
        out = "["
        for (i = 1; i <= JN[id]; i++) out = out (i > 1 ? ", " : "") json_dumps_ordered(id SUBSEP JK[id, i])
        return out "]"
    }
    return _pj_scalar_json(id)
}

function json_indent1(id, level,   t, i, out, k, pad, pad1) {
    t = JT[id]; pad = sprintf("%" (level + 0) "s", ""); pad1 = pad " "
    if (t == "map") {
        if (JN[id] == 0) return "{}"
        out = "{"
        for (i = 1; i <= JN[id]; i++) { k = JK[id, i]; out = out (i > 1 ? "," : "") "\n" pad1 json_str(k) ": " json_indent1(id SUBSEP k, level + 1) }
        return out "\n" pad "}"
    }
    if (t == "list") {
        if (JN[id] == 0) return "[]"
        out = "["
        for (i = 1; i <= JN[id]; i++) out = out (i > 1 ? "," : "") "\n" pad1 json_indent1(id SUBSEP JK[id, i], level + 1)
        return out "\n" pad "]"
    }
    return _pj_scalar_json(id)
}

# ---------------------------------------------------------------- Python str()/repr()

# Python float repr: shortest %.{p}e (p = 0..16) that round-trips, laid out with Python's rules
function pyrepr_float(x,   p, s, sign, mant, ex, digits, out, i) {
    if (x == 0) return (sprintf("%g", x) ~ /^-/) ? "-0.0" : "0.0"
    for (p = 0; p <= 16; p++) { s = sprintf("%." p "e", x); if (s + 0 == x) break }
    sign = ""; if (substr(s, 1, 1) == "-") { sign = "-"; s = substr(s, 2) }
    mant = s; sub(/e.*/, "", mant); ex = s; sub(/.*e/, "", ex); ex = ex + 0
    digits = mant; gsub(/\./, "", digits); sub(/0+$/, "", digits); if (digits == "") digits = "0"
    if (ex < -4 || ex >= 16) {
        out = substr(digits, 1, 1) (length(digits) > 1 ? "." substr(digits, 2) : "") "e" (ex < 0 ? "-" : "+") sprintf("%02d", ex < 0 ? -ex : ex)
    } else if (ex < 0) {
        out = "0."; for (i = 1; i < -ex; i++) out = out "0"; out = out digits
    } else if (length(digits) <= ex + 1) {
        out = digits; for (i = length(digits); i <= ex; i++) out = out "0"; out = out ".0"
    } else {
        out = substr(digits, 1, ex + 1) "." substr(digits, ex + 2)
    }
    return sign out
}

# code points Python's str.isprintable() rejects, approximated: exact for ASCII and Latin-1; beyond that the
# separator/format/surrogate/private ranges below (unassigned code points are not known here)
function _pj_unprintable(cp) {
    return (cp < 32 || (cp >= 127 && cp <= 160) || cp == 173 || cp == 5760 || (cp >= 8192 && cp <= 8207) ||
            (cp >= 8232 && cp <= 8239) || (cp >= 8287 && cp <= 8303) || cp == 12288 || cp == 65279 ||
            (cp >= 55296 && cp <= 63743) || (cp >= 65529 && cp <= 65531) || cp >= 983040)
}

function py_repr_str(s,   q, out, i, n, c, b, len, cp) {
    q = (index(s, "'") && !index(s, "\"")) ? "\"" : "'"
    out = q; n = length(s); i = 1
    while (i <= n) {
        c = substr(s, i, 1); b = _PJ_ORD[c]
        if (c == "\\") { out = out "\\\\"; i++; continue }
        if (c == q) { out = out "\\" c; i++; continue }
        if (c == "\t") { out = out "\\t"; i++; continue }
        if (c == "\n") { out = out "\\n"; i++; continue }
        if (c == "\r") { out = out "\\r"; i++; continue }
        if (b >= 32 && b < 127) { out = out c; i++; continue }
        if (b < 128) { out = out sprintf("\\x%02x", b); i++; continue }
        len = _pj_utf8_at(s, i); cp = _PJ_CP
        if (!_pj_unprintable(cp)) out = out substr(s, i, len)
        else if (cp < 256) out = out sprintf("\\x%02x", cp)
        else if (cp < 65536) out = out sprintf("\\u%04x", cp)
        else out = out sprintf("\\U%08x", cp)
        i += len
    }
    return out q
}

function py_repr(id,   t, i, out, k) {
    t = JT[id]
    if (t == "map") {
        out = "{"
        for (i = 1; i <= JN[id]; i++) { k = JK[id, i]; out = out (i > 1 ? ", " : "") py_repr_str(k) ": " py_repr(id SUBSEP k) }
        return out "}"
    }
    if (t == "list") {
        out = "["
        for (i = 1; i <= JN[id]; i++) out = out (i > 1 ? ", " : "") py_repr(id SUBSEP JK[id, i])
        return out "]"
    }
    if (t == "str" || t == "path") return py_repr_str(JV[id])
    return JV[id]
}

function py_str(id,   t) {
    t = JT[id]
    if (t == "map" || t == "list") return py_repr(id)
    return JV[id]
}
