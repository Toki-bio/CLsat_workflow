# yaml.awk -- YAML-subset reader for satpipe configs (PORT_SPEC.md section 3, work package 2).
#
# Usage:   LC_ALL=C gawk -f lib/satpipe/yaml.awk CONFIG.yaml
#
# Output (stdout, LF line ends), the RAW flat form, one line per node below the top-level mapping, in document order
# (mapping keys in first-occurrence order, list items in order):
#     path<TAB>type<TAB>text
#   path  components joined by "/"; a mapping key is Python str(key) of the key PyYAML would build (True, None, 1.5 ...);
#         a list item is its 0-based index (hapdup/pairs/0/1).
#   type  map | list | str | int | float | bool | null
#   text  map/list: number of entries; int: Python str(int) (sign "+" dropped, -0 -> 0, no size limit);
#         float: Python repr(float) (50.0, 0.9, 1e-05, -0.0); bool: True/False; null: None; str: the decoded string.
#   Escapes, so that every node stays on one line: in text and in key components  \ -> \\  TAB -> \t  LF -> \n  CR -> \r;
#   in key components also  / -> \x2f.  Nothing else is escaped (non-ASCII is passed through as UTF-8 bytes).
# An empty document, a document of comments only, "---" alone or a top-level null gives no lines (config.py: None -> {}).
#
# Errors: exactly one line on stdout and exit status 1, no node lines:
#     !error<TAB>YAML subset: <reason> (line N)
#     !error<TAB>top level must be a mapping                       (text of config.py:99)
#
# Accepted: block mappings to any depth (spaces only, each mapping at one consistent column); flow mappings {k: v, ...}
# and flow sequences [a, b], nested ([[a, b]], {k: [..]}), trailing comma allowed, with line breaks and comments allowed
# BETWEEN their tokens (configs/darevskia_clsat.yaml on branch flank-single-copy splits `orthology: {...}` over two
# lines; PyYAML ignores the indentation of such continuation lines, and so does this reader); double-quoted scalars with
# \" \\ \n \t \uXXXX; single-quoted with ''; plain scalars (spaces, commas, parentheses, '#' not preceded by a space,
# ':' not followed by a space); comments after a space or after a quoted/flow value; blank lines; `key:` alone -> null;
# "---" as the first content line; CRLF line ends (the CR is dropped, as Python's universal newlines do); a leading BOM.
# Typing of plain scalars = PyYAML 6 (YAML 1.1 resolver) for the kept forms; duplicate keys: the last value wins at the
# position of the first key (PyYAML/dict behaviour, keys 1, 1.0 and true being the same dict key).
# Refused (PyYAML would accept them, often with a silently changed type): block sequences, multi-line scalars (| > and
# continuation lines also inside [..]/{..}, quotes left open), anchors/aliases/tags, complex keys, directives, several
# documents, a flow key and its ':' on different lines, single-pair mappings inside [..], tab characters outside
# quotes/comments, empty keys, octal 010, 0x/0b, 1_000 and 1_0.5, sexagesimal 1:30, .inf/.nan, timestamps (2026-11-08), '=' and '<<', escapes other
# than the five above, float overflow, carriage returns inside a line, control characters, invalid UTF-8, C1 controls and
# U+2028/U+2029 (YAML line breaks).
#
# Must run with gawk >= 5 and LC_ALL=C (bytes). Contains no gawk-5.3.2-only features.

{ RAWL[++NLINES] = $0 }

END { parse_doc(); if (ROOT) emit(ROOT, "", 0); exit 0 }

function fail(msg) { printf "!error\tYAML subset: %s (line %d)\n", msg, LN; exit 1 }
function top_fail() { printf "!error\ttop level must be a mapping\n"; exit 1 }

# ---------- line preparation ----------
function prep(i,   s) {
  LN = i; s = RAWL[i]
  sub(/\r$/, "", s)
  if (index(s, "\r")) fail("carriage return inside a line")
  if (i == 1 && substr(s, 1, 3) == "\357\273\277") s = substr(s, 4)
  if (s ~ /[\001-\010\013\014\016-\037\177]/) fail("control character in the file")
  if (s !~ /^([\001-\177]|[\302-\337][\200-\277]|\340[\240-\277][\200-\277]|[\341-\354\356\357][\200-\277][\200-\277]|\355[\200-\237][\200-\277]|\360[\220-\277][\200-\277][\200-\277]|[\361-\363][\200-\277][\200-\277][\200-\277]|\364[\200-\217][\200-\277][\200-\277])*$/)
    fail("invalid UTF-8")
  if (s ~ /\302[\200-\237]/) fail("C1 control character (U+0080-U+009F, incl. NEL) in the file")
  if (s ~ /\342\200[\250\251]/) fail("U+2028/U+2029 line separator in the file")
  if (s ~ /\357\277[\276\277]/) fail("non-character U+FFFE/U+FFFF in the file")
  S = s; LEN = length(s); P = 1
}

# ---------- scanning helpers (globals S, LEN, P) ----------
function skip_sp(   c) {           # block context: spaces; a tab is refused
  while (P <= LEN) {
    c = substr(S, P, 1)
    if (c == " ") P++
    else if (c == "\t") fail("tab character outside a quoted scalar or comment")
    else break
  }
}
function skip_sp_flow() {         # inside [..] / {..}: spaces, comments and line breaks between tokens
  skip_sp()
  while (P > LEN || substr(S, P, 1) == "#") {
    if (++CUR > NLINES) fail("flow collection not closed before the end of the file")
    prep(CUR)
    if (S ~ /^(---|\.\.\.)([ 	]|$)/) fail("document marker inside a flow collection")
    while (P <= LEN && substr(S, P, 1) == " ") P++
    if (substr(S, P, 1) == "	") fail("tab indentation")
    CROSSED = 1
  }
}
function is_ws(c) { return c == "" || c == " " || c == "\t" }

function refuse_indicator(c, nx, flow) {   # c cannot start a plain scalar: say why
  if (c == "-" && is_ws(nx)) fail(flow ? "block sequence entry inside a flow collection" : "block sequences ('- item') are not supported")
  if (c == "?" && (is_ws(nx) || flow)) fail("complex keys ('? ') are not supported")
  if (c == ":") fail(is_ws(nx) ? "empty key (':' without a key)" : "a value cannot start with ':' inside a flow collection")
  if (c == "&" || c == "*") fail("anchors and aliases are not supported")
  if (c == "!") fail("tags are not supported")
  if (c == "|" || c == ">") fail("multi-line (block) scalars are not supported")
  if (c == "%") fail("directives are not supported")
  if (c == "@" || c == "`") fail("reserved character '" c "' cannot start a value")
  if (c == "#") fail("comment where a value is expected")
  if (c == "," || c == "]" || c == "}" || c == "[" || c == "{") fail("unexpected '" c "'")
}
function plain_ok(c, nx, flow) {
  if (index("-?:,[]{}#&*!|>'\"%@`", c) == 0) return 1
  return !is_ws(nx) && (c == "-" || (!flow && (c == "?" || c == ":")))
}

function scan_plain(flow,   st, q, e, ch, nx, r) {
  st = P; e = P - 1; q = P
  while (q <= LEN) {
    ch = substr(S, q, 1)
    if (ch == " ") {
      r = q; while (substr(S, r, 1) == " ") r++
      if (r > LEN || substr(S, r, 1) == "#") break
      q = r; continue
    }
    if (ch == "\t") fail("tab character outside a quoted scalar or comment")
    nx = substr(S, q + 1, 1)
    if (ch == ":" && (is_ws(nx) || (flow && index(",[]{}", nx)))) break
    if (flow && index(",?[]{}", ch)) break
    e = q; q++
  }
  P = e + 1
  return substr(S, st, e - st + 1)
}

function hexval(h,   i, v) { v = 0; h = tolower(h); for (i = 1; i <= length(h); i++) v = v * 16 + index("0123456789abcdef", substr(h, i, 1)) - 1; return v }
function utf8(cp) {
  if (cp < 128) return sprintf("%c", cp)
  if (cp < 2048) return sprintf("%c%c", 192 + int(cp / 64), 128 + cp % 64)
  return sprintf("%c%c%c", 224 + int(cp / 4096), 128 + int(cp / 64) % 64, 128 + cp % 64)
}

function parse_quoted(   q, out, ch, e, h, cp) {
  q = substr(S, P, 1); P++; out = ""
  while (1) {
    if (P > LEN) fail("quoted scalar not closed on its line (multi-line scalars are not supported)")
    ch = substr(S, P, 1)
    if (q == "'") {
      if (ch == "'") { if (substr(S, P + 1, 1) == "'") { out = out "'"; P += 2; continue } P++; return out }
      out = out ch; P++; continue
    }
    if (ch == "\"") { P++; return out }
    if (ch == "\\") {
      e = substr(S, P + 1, 1)
      if (e == "\"") out = out "\""
      else if (e == "\\") out = out "\\"
      else if (e == "n") out = out "\n"
      else if (e == "t") out = out "\t"
      else if (e == "u") {
        h = substr(S, P + 2, 4)
        if (h !~ /^[0-9A-Fa-f][0-9A-Fa-f][0-9A-Fa-f][0-9A-Fa-f]$/) fail("\\u escape needs four hex digits")
        cp = hexval(h)
        if (cp >= 55296 && cp <= 57343) fail("\\u escape of a UTF-16 surrogate")
        out = out utf8(cp); P += 6; continue
      }
      else if (e == "") fail("backslash at the end of a line (multi-line scalars are not supported)")
      else fail("unsupported escape '\\" e "' (allowed: \\\" \\\\ \\n \\t \\uXXXX)")
      P += 2; continue
    }
    out = out ch; P++
  }
}

# ---------- PyYAML 1.1 implicit typing of plain scalars; sets RT (type) and RV (text) ----------
function resolve(tok,   x, kind) {
  if (tok ~ /^(~|null|Null|NULL)$/) { RT = "null"; RV = "None"; return }
  if (tok ~ /^(true|True|TRUE|yes|Yes|YES|on|On|ON)$/) { RT = "bool"; RV = "True"; return }
  if (tok ~ /^(false|False|FALSE|no|No|NO|off|Off|OFF)$/) { RT = "bool"; RV = "False"; return }
  if (tok ~ /^[-+]?(0|[1-9][0-9]*)$/) {
    x = tok; sub(/^\+/, "", x); if (x == "-0") x = "0"
    RT = "int"; RV = x; return
  }
  if (tok ~ /^[-+]?[0-9]+\.[0-9]*([eE][-+][0-9]+)?$/ || tok ~ /^\.[0-9]+([eE][-+][0-9]+)?$/) {
    x = tok + 0
    if (x > 1.7976931348623157e308 || x < -1.7976931348623157e308) fail("float '" tok "' overflows (PyYAML would read inf)")
    RT = "float"; RV = pyrepr(x, substr(tok, 1, 1) == "-"); return
  }
  kind = ""
  if (tok ~ /^[-+]?0b[01_]+$/) kind = "a binary integer"
  else if (tok ~ /^[-+]?0x[0-9a-fA-F_]+$/) kind = "a hexadecimal integer"
  else if (tok ~ /^[-+]?0[0-7_]+$/) kind = "an octal integer"
  else if (tok ~ /^[-+]?(0|[1-9][0-9_]*)$/) kind = "an integer with '_' separators"
  else if (tok ~ /^[-+]?[1-9][0-9_]*(:[0-5]?[0-9])+$/ || tok ~ /^[-+]?[0-9][0-9_]*(:[0-5]?[0-9])+\.[0-9_]*$/) kind = "a sexagesimal (base 60) number"
  else if (tok ~ /^[-+]?[0-9][0-9_]*\.[0-9_]*([eE][-+][0-9]+)?$/ || tok ~ /^\.[0-9][0-9_]*([eE][-+][0-9]+)?$/) kind = "a float with '_' separators"
  else if (tok ~ /^[-+]?\.(inf|Inf|INF)$/ || tok ~ /^\.(nan|NaN|NAN)$/) kind = "infinity/NaN"
  else if (tok ~ /^[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]$/ || tok ~ /^[0-9][0-9][0-9][0-9]-[0-9][0-9]?-[0-9][0-9]?([Tt]|[ \t]+)[0-9][0-9]?:[0-9][0-9]:[0-9][0-9](\.[0-9]*)?([ \t]*(Z|[-+][0-9][0-9]?(:[0-9][0-9])?))?$/) kind = "a date/timestamp"
  else if (tok == "=") kind = "the YAML 'value' tag"
  else if (tok == "<<") kind = "a merge key"
  if (kind != "") fail("'" tok "' would be read by PyYAML as " kind "; quote it if a string is meant")
  RT = "str"; RV = tok
}

# Python float repr: shortest %.{p}e (p = 0..16) that round-trips, laid out with Python's rules (PORT_SPEC 2.1)
function pyrepr(x, neg,   p, s, sign, mant, ex, digits, out, i) {
  if (x == 0) return neg ? "-0.0" : "0.0"
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

# ---------- tree ----------
function new_node(t, v) { NN++; T[NN] = t; V[NN] = v; CNT[NN] = 0; return NN }
function key_ident(t, v,   x) {                     # Python dict key equality: 1 == 1.0 == True
  if (t == "str") return "s" SUBSEP v
  if (t == "null") return "z"
  if (t == "bool") return "n" SUBSEP (v == "True" ? "1" : "0")
  if (t == "int") return "n" SUBSEP v
  if (v ~ /^-?[0-9]+\.0$/) { x = v; sub(/\.0$/, "", x); if (x == "-0") x = "0"; return "n" SUBSEP x }
  return "f" SUBSEP v
}
function map_put(m, kt, kv, child,   id, j) {
  id = key_ident(kt, kv)
  if ((m, id) in KIDX) { CH[m, KIDX[m, id]] = child; return }
  j = ++CNT[m]; KIDX[m, id] = j; KSTR[m, j] = kv; CH[m, j] = child
}
function list_add(l, child) { CH[l, ++CNT[l]] = child }

# ---------- flow collections (one line) ----------
function parse_flow_node(   c) {
  skip_sp_flow(); c = substr(S, P, 1)
  if (c == "{") return parse_flow_map()
  if (c == "[") return parse_flow_seq()
  if (c == "\"" || c == "'") return new_node("str", parse_quoted())
  if (!plain_ok(c, substr(S, P + 1, 1), 1)) refuse_indicator(c, substr(S, P + 1, 1), 1)
  resolve(scan_plain(1)); return new_node(RT, RV)
}
function parse_flow_seq(   l, c) {
  P++; l = new_node("list", 0)
  skip_sp_flow(); if (substr(S, P, 1) == "]") { P++; return l }
  while (1) {
    list_add(l, parse_flow_node())
    CROSSED = 0; skip_sp_flow(); c = substr(S, P, 1)
    if (c == ":") fail("single-pair mappings inside a flow sequence are not supported")
    if (CROSSED && c != "," && c != "]") fail("plain scalar continued on the next line (multi-line scalars are not supported)")
    if (c == "]") { P++; return l }
    if (c != ",") fail("expected ',' or ']' in a flow sequence")
    P++; skip_sp_flow(); c = substr(S, P, 1)
    if (c == "]") { P++; return l }
    if (c == ",") fail("empty entry in a flow sequence")
  }
}
function parse_flow_map(   m, c, nx, kt, kv, val) {
  P++; m = new_node("map", 0)
  skip_sp_flow(); if (substr(S, P, 1) == "}") { P++; return m }
  while (1) {
    c = substr(S, P, 1); nx = substr(S, P + 1, 1)
    if (c == "{" || c == "[") fail("flow collections as mapping keys are not supported")
    if (c == "\"" || c == "'") { kt = "str"; kv = parse_quoted() }
    else {
      if (!plain_ok(c, nx, 1)) refuse_indicator(c, nx, 1)
      resolve(scan_plain(1)); kt = RT; kv = RV
    }
    CROSSED = 0; skip_sp_flow(); c = substr(S, P, 1)
    if (CROSSED && c != "," && c != "}") fail("a flow mapping key must be followed by ':' on its own line (multi-line keys are not supported)")
    if (c == ":") {
      P++; skip_sp_flow(); c = substr(S, P, 1)
      if (c == "," || c == "}") val = new_node("null", "None")
      else {
        val = parse_flow_node(); CROSSED = 0; skip_sp_flow(); c = substr(S, P, 1)
        if (CROSSED && c != "," && c != "}") fail("plain scalar continued on the next line (multi-line scalars are not supported)")
      }
    } else val = new_node("null", "None")
    map_put(m, kt, kv, val)
    if (c == "}") { P++; return m }
    if (c != ",") fail("expected ',' or '}' in a flow mapping")
    P++; skip_sp_flow(); c = substr(S, P, 1)
    if (c == "}") { P++; return m }
    if (c == ",") fail("empty entry in a flow mapping")
  }
}
function expect_end(   c) {        # after a complete value: spaces, then a comment or the end of the line
  skip_sp(); if (P > LEN) return
  c = substr(S, P, 1)
  if (c == "#") return
  if (c == ":") fail("mapping values are not allowed here (a nested mapping must start on its own line)")
  fail("unexpected text after the value")
}

# ---------- block structure ----------
# parse_key: 1 = key read (KT, KV set, P after ':'); 0 = no 'key:' on this line; -2 = line starts a flow collection
function parse_key(   c, nx) {
  c = substr(S, P, 1); nx = substr(S, P + 1, 1)
  if (c == "[" || c == "{") return -2
  if (c == "\"" || c == "'") { KT = "str"; KV = parse_quoted() }
  else {
    if (!plain_ok(c, nx, 0)) refuse_indicator(c, nx, 0)
    resolve(scan_plain(0)); KT = RT; KV = RV
  }
  skip_sp()
  if (substr(S, P, 1) == ":" && is_ws(substr(S, P + 1, 1))) { P++; return 1 }
  return 0
}
function parse_value(m,   kt, kv, c, nx, node) {   # after 'key:' of mapping m
  kt = KT; kv = KV
  skip_sp()
  if (P > LEN || substr(S, P, 1) == "#") { PEND = 1; PM = m; PKT = kt; PKV = kv; return }
  c = substr(S, P, 1); nx = substr(S, P + 1, 1)
  if (c == "{") node = parse_flow_map()
  else if (c == "[") node = parse_flow_seq()
  else if (c == "\"" || c == "'") node = new_node("str", parse_quoted())
  else {
    if (!plain_ok(c, nx, 0)) refuse_indicator(c, nx, 0)
    resolve(scan_plain(0)); node = new_node(RT, RV)
  }
  expect_end()
  map_put(m, kt, kv, node)
}

function parse_doc(   ind, c, r, seen, marker, nulldoc, m) {
  ROOT = 0; D = 0; PEND = 0; seen = 0; marker = 0; nulldoc = 0
  for (CUR = 1; CUR <= NLINES; CUR++) {          # a flow collection spanning lines advances CUR itself
    prep(CUR)
    while (P <= LEN && substr(S, P, 1) == " ") P++
    if (P > LEN) continue
    c = substr(S, P, 1)
    if (c == "\t") fail("tab indentation")
    if (c == "#") continue
    ind = P - 1
    if (ind == 0 && (S ~ /^---([ \t]|$)/ || S ~ /^\.\.\.([ \t]|$)/)) {
      if (!seen && !marker && substr(S, 1, 3) == "---") { P = 4; skip_sp(); if (P > LEN || substr(S, P, 1) == "#") { marker = 1; continue } }
      fail("several documents or document markers are not supported")
    }
    if (nulldoc || ROOTFLOW) fail("content after the end of the top-level value")
    if (!ROOT) {
      seen = 1
      if (c == "-" && is_ws(substr(S, P + 1, 1))) top_fail()     # a top-level block sequence is a list (config.py:99)
      r = parse_key()
      if (r == -2) {
        if (c == "[") { parse_flow_seq(); expect_end(); top_fail() }
        ROOT = parse_flow_map(); expect_end(); ROOTFLOW = 1; continue
      }
      if (r == 0) {
        expect_end()
        if (KT == "null") { nulldoc = 1; continue }
        top_fail()
      }
      ROOT = new_node("map", 0); D = 1; SI[1] = ind; SM[1] = ROOT
      parse_value(ROOT); continue
    }
    if (PEND) {
      PEND = 0
      if (ind > SI[D]) {
        m = new_node("map", 0); map_put(PM, PKT, PKV, m)
        D++; SI[D] = ind; SM[D] = m
        r = parse_key()
        if (r != 1) fail("expected 'key: value' (block sequences, multi-line scalars and flow values on their own line are not supported)")
        parse_value(m); continue
      }
      map_put(PM, PKT, PKV, new_node("null", "None"))
    }
    while (D > 0 && SI[D] > ind) D--
    if (D == 0) fail("line indented less than the top-level mapping")
    if (SI[D] < ind) fail("unexpected indentation (multi-line scalars are not supported)")
    r = parse_key()
    if (r != 1) fail("expected 'key: value' (block sequences, multi-line scalars and flow values on their own line are not supported)")
    parse_value(SM[D])
  }
  if (PEND) map_put(PM, PKT, PKV, new_node("null", "None"))
}

# ---------- output ----------
function esc_text(s) {
  gsub(/\\/, "\\\\\\\\", s); gsub(/\t/, "\\t", s); gsub(/\n/, "\\n", s); gsub(/\r/, "\\r", s)
  return s
}
function esc_key(s) { s = esc_text(s); gsub(/\//, "\\x2f", s); return s }
function emit(n, prefix, depth,   j, c, p) {
  for (j = 1; j <= CNT[n]; j++) {
    c = CH[n, j]
    if (T[n] == "map") p = (depth ? prefix "/" : "") esc_key(KSTR[n, j])
    else p = prefix "/" (j - 1)
    if (T[c] == "map" || T[c] == "list") { print p "\t" T[c] "\t" CNT[c]; emit(c, p, depth + 1) }
    else print p "\t" T[c] "\t" esc_text(V[c])
  }
}
