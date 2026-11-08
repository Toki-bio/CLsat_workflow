# check.awk - the file readers and the report renderer of `satpipe check` (check.py), used by check.sh and tools.sh.
# Library: gawk -f pyjson.awk -f check.awk -e 'BEGIN { <entry>() }'. Needs LC_ALL=C (bytes).
# Arguments come through ENVIRON (never -v, which would interpret backslashes).
#
# Python reads these files in text mode, so the readers copy that: universal newlines (\r\n, \r and \n all end a line,
# RS below), UTF-8 with errors="replace" (ck_utf8_fix: one U+FFFD per maximal invalid subpart), str.strip()/split() with
# Python's whitespace set (ck_strip, ck_ws_index), lengths in code points (ck_clen).
#
# Entries (stdout = report items `LEVEL<TAB>SECTION<TAB>MESSAGE`, MESSAGE flat-escaped, unless noted)
#   ck_fasta_main()     P, MAX            read_fasta (check.py:53-70): `name<TAB>length` rows (name flat-escaped)
#   ck_famtab_main()    P, FAMS           check_family_assignment after its _exists (check.py:82-115); FAMS = names joined by \037, config order
#   ck_bed_main()       P, LABEL, SECTION, CODES   check_bed after its _exists (check.py:118-140); CODES = genome codes joined by \037
#   ck_version_probe()  ARGV[1]=stdout file, ARGV[2]=stderr file, ERE, MINV, NAME   tools.py:80-93: `version<TAB>error` (escaped)
#   ck_render_main()    P = items file, FORMAT = text | json      Report.to_text / json.dumps(to_json(), indent=2) (check.py:38-47)

BEGIN { CK_UNIRS = "\r\n|\r|\n" }

# ------------------------------------------------------------------ Python text helpers

# bytes -> valid UTF-8 as Python's decode(errors="replace") would give it
function ck_utf8_fix(s,   out, i, n, b, need, lo, hi, k, j, c) {
    if (s !~ /[\200-\377]/) return s
    out = ""; n = length(s); i = 1
    while (i <= n) {
        c = substr(s, i, 1); b = _PJ_ORD[c]
        if (b < 128) { out = out c; i++; continue }
        lo = 128; hi = 191
        if (b >= 194 && b <= 223) need = 1
        else if (b == 224) { need = 2; lo = 160 }
        else if (b >= 225 && b <= 236 || b == 238 || b == 239) need = 2
        else if (b == 237) { need = 2; hi = 159 }
        else if (b == 240) { need = 3; lo = 144 }
        else if (b >= 241 && b <= 243) need = 3
        else if (b == 244) { need = 3; hi = 143 }
        else { out = out "\357\277\275"; i++; continue }
        for (k = 1; k <= need; k++) {
            j = _PJ_ORD[substr(s, i + k, 1)]
            if (i + k > n || j < lo || j > hi) break
            lo = 128; hi = 191
        }
        if (k > need) { out = out substr(s, i, need + 1); i += need + 1 }
        else { out = out "\357\277\275"; i += k }
    }
    return out
}

# Python str.isspace() code points, as UTF-8 regex alternatives
function _ck_ws() {
    return "[\t\n\v\f\r\034-\037 ]|\302[\205\240]|\341\232\200|\342\200[\200-\212\250\251\257]|\342\201\237|\343\200\200"
}
function ck_strip(s,   w) {
    w = _ck_ws()
    sub("^(" w ")+", "", s); sub("(" w ")+$", "", s)
    return s
}
# position of the first whitespace character in s (0 = none)
function ck_ws_index(s) { return match(s, _ck_ws()) ? RSTART : 0 }
# number of code points of valid UTF-8 text
function ck_clen(s,   t) { t = s; return length(s) - gsub(/[\200-\277]/, "", t) }
# first n code points
function ck_head(s, n,   i, k, b) {
    k = 0
    for (i = 1; i <= length(s); i++) {
        b = _PJ_ORD[substr(s, i, 1)]
        if (b < 128 || b >= 192) { if (++k > n) return substr(s, 1, i - 1) }
    }
    return s
}

# Python int(x) for a BED column: "" when int() raises ValueError, else the canonical digits with sign
function ck_int(x,   t, neg) {
    t = ck_strip(x)
    if (t !~ /^[-+]?[0-9](_?[0-9])*$/) return ""
    gsub(/_/, "", t)
    neg = (substr(t, 1, 1) == "-"); sub(/^[-+]/, "", t)
    if (length(t) > 4300) return ""                                  # sys.int_info.default_max_str_digits
    sub(/^0+/, "", t); if (t == "") return "0"
    return (neg ? "-" : "") t
}
# compare canonical integers: -1, 0, 1 (no precision limit)
function ck_icmp(a, b,   na, nb, r) {
    na = (substr(a, 1, 1) == "-"); nb = (substr(b, 1, 1) == "-")
    if (na != nb) return na ? -1 : 1
    if (na) { a = substr(a, 2); b = substr(b, 2) }
    if (length(a) != length(b)) r = (length(a) < length(b)) ? -1 : 1
    else r = (a < b) ? -1 : (a > b) ? 1 : 0
    return na ? -r : r
}

function ck_item(level, section, msg) { printf "%s\t%s\t%s\n", level, section, flat_esc(msg) }

# sorted(xs) of \037-joined strings, joined by sep
function ck_sorted_join(list, sep,   a, s, n, i, out) {
    n = split(list, a, "\037"); n = asort(a, s, "@val_str_asc")
    out = ""; for (i = 1; i <= n; i++) out = out (i > 1 ? sep : "") s[i]
    return out
}

# ------------------------------------------------------------------ read_fasta (check.py:53-70)
# fills CKF_NAME[1..k], CKF_LEN[1..k]; returns k
function ck_read_fasta(path, maxrec,   line, k, name, have, n, hdr, w, ors) {
    delete CKF_NAME; delete CKF_LEN
    ors = RS; RS = CK_UNIRS
    k = 0; have = 0; n = 0
    while ((getline line < path) > 0) {
        line = ck_utf8_fix(line)
        if (substr(line, 1, 1) == ">") {
            if (have) {
                CKF_NAME[++k] = name; CKF_LEN[k] = n
                if (k >= maxrec) { close(path); RS = ors; return k }
            }
            hdr = ck_strip(substr(line, 2))
            w = ck_ws_index(hdr)
            name = w ? substr(hdr, 1, w - 1) : hdr
            n = 0; have = 1
        } else if (have) n += ck_clen(ck_strip(line))
    }
    close(path); RS = ors
    if (have) { CKF_NAME[++k] = name; CKF_LEN[k] = n }
    return k
}

function ck_fasta_main(   k, i) {
    k = ck_read_fasta(ENVIRON["P"], ENVIRON["MAX"] + 0)
    for (i = 1; i <= k; i++) printf "%s\t%d\n", flat_esc(CKF_NAME[i]), CKF_LEN[i]
}

# ------------------------------------------------------------------ check_family_assignment (check.py:82-115)
function ck_famtab_main(   path, sec, nf, F, i, FAM, SEEN, CLU, line, ln, p, np, clu, fam, nrows, sorted, out) {
    path = ENVIRON["P"]; sec = "family_assignment"
    nf = split(ENVIRON["FAMS"], F, "\037")
    for (i = 1; i <= nf; i++) { FAM[F[i]] = 1; SEEN[F[i]] = 0 }
    sorted = ck_sorted_join(ENVIRON["FAMS"], ", ")
    RS = CK_UNIRS; ln = 0; nrows = 0
    while ((getline line < path) > 0) {
        ln++
        line = ck_utf8_fix(line)
        if (ck_strip(line) == "" || substr(line, 1, 1) == "#") continue
        np = split(line, p, "\t")
        if (np < 2) { ck_item("FAIL", sec, "line " ln ": need `cluster<TAB>family`, got " py_repr_str(ck_strip(line))); continue }
        clu = ck_strip(p[1]); fam = ck_strip(p[2])
        if (clu !~ /^[0-9]+$/) { ck_item("FAIL", sec, "line " ln ": cluster must be a non-negative integer, got " py_repr_str(clu)); continue }
        if (!(fam in FAM)) { ck_item("FAIL", sec, "line " ln ": family " py_repr_str(fam) " is not in `families` (" sorted ")"); continue }
        if (clu in CLU) { ck_item("FAIL", sec, "line " ln ": cluster " clu " listed twice (" CLU[clu] " and " fam ")"); continue }
        CLU[clu] = fam; SEEN[fam]++; nrows++
    }
    close(path)
    # Python iterates a set here (hash order, varies with PYTHONHASHSEED); config order is used
    for (i = 1; i <= nf; i++) if (!SEEN[F[i]]) ck_item("FAIL", sec, "family " F[i] " has no cluster in the table")
    if (nrows) {
        split(ck_sorted_join(ENVIRON["FAMS"], "\037"), F, "\037"); out = ""
        for (i = 1; i <= nf; i++) out = out (i > 1 ? ", " : "") F[i] ": " SEEN[F[i]]
        ck_item("OK", sec, nrows " cluster->family rows; " out)
    }
    ck_item("OK", sec, "cluster existence is checked at stage 09 against the Step-8 output (needs a run)")
}

# ------------------------------------------------------------------ check_bed (check.py:118-140)
function ck_bed_main(   path, label, sec, nc, C, CODE, i, line, ln, p, np, s, e, n, nm) {
    path = ENVIRON["P"]; label = ENVIRON["LABEL"]; sec = ENVIRON["SECTION"]
    nc = split(ENVIRON["CODES"], C, "\037"); for (i = 1; i <= nc; i++) CODE[C[i]] = 1
    RS = CK_UNIRS; ln = 0; n = 0
    while ((getline line < path) > 0) {
        ln++
        line = ck_utf8_fix(line)
        if (ck_strip(line) == "" || line ~ /^(#|track|browser)/) continue
        np = split(line, p, "\t")
        if (np < 3) { ck_item("FAIL", sec, label " line " ln ": fewer than 3 columns"); continue }
        s = ck_int(p[2]); e = ck_int(p[3])
        if (s == "" || e == "") { ck_item("FAIL", sec, label " line " ln ": start/end not integers"); continue }
        if (!(ck_icmp("0", s) <= 0 && ck_icmp(s, e) < 0))
            ck_item("FAIL", sec, label " line " ln ": need 0 <= start < end (BED is 0-based half-open)")
        if (np >= 4) {
            nm = p[4]; i = index(nm, "|")
            if (!((i ? substr(nm, 1, i - 1) : nm) in CODE))
                ck_item("WARN", sec, label " line " ln ": name " py_repr_str(p[4]) " does not start with a known species code")
        }
        n++
    }
    close(path)
    ck_item("OK", sec, label ": " n " intervals")
}

# ------------------------------------------------------------------ tool version (tools.py:80-93)
function _ck_slurp(path,   ors, t, x) {
    ors = RS; RS = "^$"; t = ""
    if ((getline x < path) > 0) t = x
    close(path); RS = ors
    gsub(/\r\n/, "\n", t); gsub(/\r/, "\n", t)                       # text=True: universal newlines
    return ck_utf8_fix(t)                                             # errors="replace"
}
# tuple(int(x) for x in re.findall(r"\d+", v)[:3]) < min_version
function _ck_version_lt(v, minv,   a, b, na, nb, i, t) {
    na = 0; t = v
    while (na < 3 && match(t, /[0-9]+/)) { a[++na] = substr(t, RSTART, RLENGTH); t = substr(t, RSTART + RLENGTH) }
    nb = split(minv, b, ".")
    for (i = 1; i <= na && i <= nb; i++) {
        if (ck_icmp(ck_int(a[i]), b[i]) != 0) return ck_icmp(ck_int(a[i]), b[i]) < 0
    }
    return na < nb
}
function ck_version_probe(   out, m, found, version, err, first, w, sep) {
    out = _ck_slurp(ARGV[1]) "\n" _ck_slurp(ARGV[2])
    version = ""; err = ""; found = match(out, ENVIRON["ERE"], m)
    if (found) version = m[1]
    if (!found) {
        first = ck_strip(out)
        sep = "[\n\v\f\034\035\036]|\302\205|\342\200[\250\251]"         # str.splitlines() boundaries (\r is gone already)
        if (first == "") first = "(no output)"
        else if (match(first, sep)) first = substr(first, 1, RSTART - 1)
        err = "version not recognised: " ck_head(first, 80)
    } else if (ENVIRON["MINV"] != "" && _ck_version_lt(version, ENVIRON["MINV"]))
        err = "version " version " < required " ENVIRON["MINV"]
    if (ENVIRON["NAME"] == "gawk" && found && version != "" && !index(out, "GNU Awk")) err = "not GNU awk"
    printf "%s\t%s", flat_esc(version), flat_esc(err)
}

# ------------------------------------------------------------------ report (check.py:38-47)
function ck_render_main(   path, fmt, line, f, n, L, S, M, i, c, failed) {
    path = ENVIRON["P"]; fmt = ENVIRON["FORMAT"]
    RS = "\n"; n = 0; c["OK"] = c["WARN"] = c["FAIL"] = 0
    while ((getline line < path) > 0) {
        split(line, f, "\t")
        L[++n] = f[1]; S[n] = f[2]; M[n] = flat_unesc(f[3]); c[f[1]]++
    }
    close(path)
    failed = c["FAIL"] > 0
    if (fmt == "json") {
        printf "{\n  \"failed\": %s,\n  \"counts\": {\n    \"OK\": %d,\n    \"WARN\": %d,\n    \"FAIL\": %d\n  },\n", \
            failed ? "true" : "false", c["OK"], c["WARN"], c["FAIL"]
        if (n == 0) { printf "  \"items\": []\n}\n"; return }
        printf "  \"items\": [\n"
        for (i = 1; i <= n; i++)
            printf "    {\n      \"level\": %s,\n      \"section\": %s,\n      \"message\": %s\n    }%s\n", \
                json_str(L[i]), json_str(S[i]), json_str(M[i]), (i < n ? "," : "")
        printf "  ]\n}\n"
        return
    }
    for (i = 1; i <= n; i++) printf "[%-4s] %s: %s\n", L[i], S[i], M[i]
    printf "check %s: %d ok, %d warnings, %d failures\n", failed ? "FAILED" : "passed", c["OK"], c["WARN"], c["FAIL"]
}
