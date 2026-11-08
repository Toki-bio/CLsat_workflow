# viewer.awk - `satpipe viewer` data steps (viewer.py:48-271) in gawk. Load as
#   gawk -f pyjson.awk -f viewer_lists.awk -f viewer.awk      (BEGIN-only program; phase in ENVIRON["VW_PHASE"])
# viewer.sh runs the phases in viewer.py's order and does the file copies, gzip and flank scans itself.
#
# JSON here is read with a token-array parser (vj_parse; pyjson.awk's json_parse re-scans the rest of the text for every
# string, which is quadratic on the multi-MB flanks.json of real runs) into the same typed tree as pyjson.awk (JT/JV/JN/JK),
# so py_str/py_repr/json_str apply. Values follow Python json.load: ints keep their digits, floats become Python repr,
# NaN/Infinity/-Infinity are floats nan/inf/-inf; duplicate keys keep the first position and the last value (dict).
# vj_dump writes json.dump(x, f, separators=(',', ':')) (insertion order, ensure_ascii).
# Text files are read with Python's universal newlines (RS = \r\n | \r | \n), as open() in text mode does.
#
# Environment: VW_PHASE, VW_CFG (TYPED flat config, made by phase `manifest`), VW_RUN, VW_OUT, VW_TMP, VW_CAND (path of
# candidates_summary.json when it is a file, else empty), VW_VERSION, VW_SUBDIR_OK (1 = subfamily.consensuses is a dir),
# phase-specific ones below. Errors: message on stderr, exit 1 (Python: SystemExit texts verbatim; FileNotFoundError and
# other exceptions as `satpipe viewer: <msg>`, the bin/satpipe convention for tracebacks).

function vfail(msg) { printf "%s\n", msg > "/dev/stderr"; close("/dev/stderr"); exit 1 }
function vfail_nofile(p) { vfail("satpipe viewer: [Errno 2] No such file or directory: " py_repr_str(p)) }
function vkeyerr(k) { vfail("satpipe viewer: KeyError: " py_repr_str(k)) }

# ---------------------------------------------------------------- files

function vw_slurp(p,   s, rs, r) {
    rs = RS; RS = "^$"; s = ""
    r = (getline s < p)
    RS = rs
    if (r < 0) vfail_nofile(p)
    close(p)
    return s
}

# vw_lines(p, arr) -> n   lines of a text file with universal newlines (arr[1..n], no terminators); missing file -> error.
# VW_LT[i] = 1 when line i was terminated by a newline.
function vw_lines(p, arr,   rs, n, r, line) {
    rs = RS; RS = "\r\n|\r|\n"; n = 0
    delete arr; delete VW_LT
    while ((r = (getline line < p)) > 0) { arr[++n] = line; VW_LT[n] = (RT != "") }
    RS = rs
    if (r < 0 && n == 0) vfail_nofile(p)
    close(p)
    return n
}

function vw_isfile(p,   r, line) { r = (getline line < p); if (r >= 0) close(p); return r >= 0 }

# read_tsv (run.py:144-153): skip blank (strip() == '') and '#' lines; rows with >= min_cols tab fields
function vw_read_tsv(p, arr, min_cols,   L, n, i, k, m, f) {
    n = vw_lines(p, L); m = 0; delete arr
    for (i = 1; i <= n; i++) {
        if (L[i] ~ /^[ \t\n\r\v\f\034-\037]*$/ || substr(L[i], 1, 1) == "#") continue
        k = split(L[i], f, "\t")
        if (k < min_cols) continue
        arr[++m] = L[i]
    }
    return m
}

# ---------------------------------------------------------------- JSON (token parser, compact dumper)

function vj_parse(text, root,   n, i) {
    n = patsplit(text, _VJT, /"([^"\\]|\\.)*"|-?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][-+]?[0-9]+)?|true|false|null|NaN|-?Infinity|[][{}:,]/, _VJS)
    JERR = ""
    for (i = 0; i <= n; i++) if (_VJS[i] !~ /^[ \t\n\r]*$/) { JERR = "invalid JSON text near token " i; return 0 }
    _VJN = n; _VJI = 1
    if (root in JT) jt_del(root)
    jt_init(root, "null")
    if (!_vj_val(root)) return 0
    if (_VJI <= n) { JERR = "Extra data"; return 0 }
    return 1
}

function _vj_val(id,   t, c, k, child, n) {
    if (_VJI > _VJN) { JERR = "Expecting value"; return 0 }
    t = _VJT[_VJI++]
    if (t == "{") {
        JT[id] = "map"; JV[id] = ""; JN[id] = 0
        if (_VJT[_VJI] == "}") { _VJI++; return 1 }
        while (1) {
            t = _VJT[_VJI++]
            if (substr(t, 1, 1) != "\"") { JERR = "Expecting property name enclosed in double quotes"; return 0 }
            k = _vj_str(t)
            if (_VJT[_VJI++] != ":") { JERR = "Expecting ':' delimiter"; return 0 }
            child = jt_set(id, k, "null", "")
            if (!_vj_val(child)) return 0
            t = _VJT[_VJI++]
            if (t == "}") return 1
            if (t != ",") { JERR = "Expecting ',' delimiter"; return 0 }
        }
    }
    if (t == "[") {
        JT[id] = "list"; JV[id] = ""; JN[id] = 0
        if (_VJT[_VJI] == "]") { _VJI++; return 1 }
        n = 0
        while (1) {
            child = jt_set(id, n++, "null", "")
            if (!_vj_val(child)) return 0
            t = _VJT[_VJI++]
            if (t == "]") return 1
            if (t != ",") { JERR = "Expecting ',' delimiter"; return 0 }
        }
    }
    c = substr(t, 1, 1)
    if (c == "\"") { JT[id] = "str"; JV[id] = _vj_str(t); return 1 }
    if (t == "true")  { JT[id] = "bool"; JV[id] = "True"; return 1 }
    if (t == "false") { JT[id] = "bool"; JV[id] = "False"; return 1 }
    if (t == "null")  { JT[id] = "null"; JV[id] = "None"; return 1 }
    if (t == "NaN")   { JT[id] = "float"; JV[id] = "nan"; return 1 }
    if (t == "Infinity")  { JT[id] = "float"; JV[id] = "inf"; return 1 }
    if (t == "-Infinity") { JT[id] = "float"; JV[id] = "-inf"; return 1 }
    if (t ~ /^-?[0-9]/) {
        if (t ~ /^-?[0-9]+$/) { JT[id] = "int"; JV[id] = (t == "-0") ? "0" : t }
        else { JT[id] = "float"; JV[id] = vj_float_text(t + 0) }
        return 1
    }
    JERR = "Expecting value"; return 0
}

function vj_float_text(x,   s) {
    s = sprintf("%f", x)
    if (s ~ /nan/) return "nan"
    if (s ~ /^inf/) return "inf"
    if (s ~ /^-inf/) return "-inf"
    return pyrepr_float(x)
}

function _vj_str(t,   s, out, i, n, c, h, cp, lo, seg) {
    s = substr(t, 2, length(t) - 2)
    if (index(s, "\\") == 0) return s
    out = ""; n = length(s); i = 1
    while (i <= n) {
        c = substr(s, i, 1)
        if (c != "\\") { out = out c; i++; continue }
        c = substr(s, i + 1, 1); i += 2
        if (c == "n") out = out "\n"; else if (c == "t") out = out "\t"; else if (c == "r") out = out "\r"
        else if (c == "b") out = out "\b"; else if (c == "f") out = out "\f"
        else if (c == "u") {
            h = substr(s, i, 4); cp = strtonum("0x" h); i += 4
            if (cp >= 55296 && cp <= 56319 && substr(s, i, 2) == "\\u") {
                seg = substr(s, i + 2, 4); lo = strtonum("0x" seg)
                if (lo >= 56320 && lo <= 57343) { cp = 65536 + (cp - 55296) * 1024 + (lo - 56320); i += 6 }
            }
            out = out utf8_chr(cp)
        } else out = out c
    }
    return out
}

function vj_scalar(id,   t, v) {
    t = JT[id]; v = JV[id]
    if (t == "int") return v
    if (t == "float") return (v == "nan") ? "NaN" : (v == "inf") ? "Infinity" : (v == "-inf") ? "-Infinity" : v
    if (t == "bool") return (v == "True") ? "true" : "false"
    if (t == "null") return "null"
    return json_str(v)
}

function vj_dump(id, f,   t, i, k) {
    t = JT[id]
    if (t == "map") {
        printf "{" > f
        for (i = 1; i <= JN[id]; i++) { k = JK[id, i]; printf "%s%s:", (i > 1 ? "," : ""), json_str(k) > f; vj_dump(id SUBSEP k, f) }
        printf "}" > f; return
    }
    if (t == "list") {
        printf "[" > f
        for (i = 1; i <= JN[id]; i++) { if (i > 1) printf "," > f; vj_dump(id SUBSEP JK[id, i], f) }
        printf "]" > f; return
    }
    printf "%s", vj_scalar(id) > f
}

function vj_load(p, root,   s) {
    s = vw_slurp(p)
    if (!vj_parse(s, root)) vfail("satpipe viewer: " p ": " JERR)
}

function vj_write(root, p) { printf "" > p; vj_dump(root, p); close(p) }

# ---------------------------------------------------------------- typed-tree helpers

function vj_get(id, key,   c) { c = id SUBSEP key; return (id != "" && (JT[id] == "map") && (c in JT)) ? c : "" }
function vj_req(id, key,   c) { c = vj_get(id, key); if (c == "") vkeyerr(key); return c }

function vj_cp(src, par, key,   c, i) {
    c = jt_set(par, key, JT[src], JV[src])
    for (i = 1; i <= JN[src]; i++) vj_cp(src SUBSEP JK[src, i], c, JK[src, i])
    return c
}

function vj_truthy(id,   t) {
    if (id == "" || !(id in JT)) return 0
    t = JT[id]
    if (t == "map" || t == "list") return JN[id] > 0
    if (t == "null") return 0
    if (t == "bool") return JV[id] == "True"
    if (t == "int") return JV[id] != "0"
    if (t == "float") return JV[id] != "0.0" && JV[id] != "-0.0"
    return JV[id] != ""
}

function vj_isnum(id) { return JT[id] == "int" || JT[id] == "float" }

# Python int(x) of a node or a text (str() / int literal rules for the forms that occur); returns the digits
function py_int_text(s,   neg) {
    gsub(/^[ \t\n\r\v\f]+|[ \t\n\r\v\f]+$/, "", s)
    if (s !~ /^[-+]?[0-9]+(_[0-9]+)*$/) return ""
    gsub(/_/, "", s); neg = (substr(s, 1, 1) == "-"); sub(/^[-+]/, "", s); sub(/^0+/, "", s)
    if (s == "") return "0"
    return (neg ? "-" : "") s
}
function vw_int_s(s,   r) {
    r = py_int_text(s)
    if (r == "") vfail("satpipe viewer: ValueError: invalid literal for int() with base 10: " py_repr_str(s))
    return r
}
function vw_int(id,   t, v) {
    t = JT[id]; v = JV[id]
    if (t == "int") return v
    if (t == "bool") return (v == "True") ? "1" : "0"
    if (t == "float") { if (v ~ /nan|inf/) vfail("satpipe viewer: cannot convert float " v " to integer"); return sprintf("%d", int(v + 0)) }
    if (t == "str" || t == "path") return vw_int_s(v)
    vfail("satpipe viewer: TypeError: int() argument must be a string, a bytes-like object or a real number, not " py_repr_str(t == "null" ? "NoneType" : t == "map" ? "dict" : t))
}

# Python float(x) as a gawk number (ordering only)
function vw_float(id,   t, v) {
    t = JT[id]; v = JV[id]
    if (t == "int" || t == "float") return v + 0
    if (t == "bool") return (v == "True") ? 1 : 0
    if (t == "str" || t == "path") {
        gsub(/^[ \t\n\r\v\f]+|[ \t\n\r\v\f]+$/, "", v)
        if (v ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$/) return v + 0
        vfail("satpipe viewer: ValueError: could not convert string to float: " py_repr_str(JV[id]))
    }
    vfail("satpipe viewer: TypeError: float() argument must be a string or a real number, not " py_repr_str(t == "null" ? "NoneType" : t == "map" ? "dict" : t))
}

# '%g' % x and '%d' % x of config numbers
function vw_fmt_g(id,   v) { v = JV[id]; if (JT[id] == "float" && v ~ /nan|inf/) return v; return sprintf("%g", vw_float(id)) }
function vw_fmt_d(id) { return vw_int(id) }

# ---------------------------------------------------------------- config (TYPED flat from the manifest)

function vw_load_cfg(   f, line, fc, gc, i) {
    f = ENVIRON["VW_CFG"]; jt_init("C", "map")
    while ((getline line < f) > 0) if (!flat_line("C", line)) vfail("satpipe viewer: " JERR)
    close(f)
    fc = vj_req("C", "families"); gc = vj_req("C", "genomes")
    NFAM = JN[fc]; for (i = 1; i <= NFAM; i++) FAM[i] = JK[fc, i]
    NCODE = JN[gc]; for (i = 1; i <= NCODE; i++) { CODE[i] = JK[gc, i]; GPATH[i] = JV[gc SUBSEP CODE[i]] }
}

# rows of the viewer BED (written by phase rows): R_SCF R_S R_E R_SP R_LID, 1..NR_ROWS
function vw_load_rows(   f, line, p) {
    f = ENVIRON["VW_TMP"] "/rows.tsv"; NR_ROWS = 0
    while ((getline line < f) > 0) {
        split(line, p, "\t"); NR_ROWS++
        R_SCF[NR_ROWS] = flat_unesc(p[1]); R_S[NR_ROWS] = p[2]; R_E[NR_ROWS] = p[3]; R_SP[NR_ROWS] = flat_unesc(p[4]); R_LID[NR_ROWS] = flat_unesc(p[5])
    }
    close(f)
}

# candidates_summary.json -> tree "K"; CAND = 1 when loaded; CS = summary node; CSORT[1..NCS] = ids in viewer.py:76 order
function vw_load_cand(   p, i, j, k, key, n, tk, tp) {
    CAND = 0; p = ENVIRON["VW_CAND"]
    if (p == "") return
    vj_load(p, "K"); CAND = 1
    if (JT["K"] != "map") vfail("satpipe viewer: TypeError: list indices must be integers or slices, not str")
    CS = vj_req("K", "summary")
    n = JN[CS]; NCS = 0
    for (i = 1; i <= n; i++) {                       # stable insertion sort by int(id[1:]) (non-digit -> 0)
        k = JK[CS, i]; key = substr(k, 2)
        key = (key ~ /^[0-9]+$/) ? key + 0 : 0
        for (j = NCS; j >= 1 && CKEY[j] > key; j--) { CSORT[j + 1] = CSORT[j]; CKEY[j + 1] = CKEY[j] }
        CSORT[j + 1] = k; CKEY[j + 1] = key; NCS++
    }
}

# ---------------------------------------------------------------- phases

# manifest: TYPED flat of manifest['config'] on stdout (viewer.py:43-45; the legacy-manifest converter)
function phase_manifest(   p, c) {
    p = ENVIRON["VW_RUN"] "/manifest.json"
    vj_load(p, "M")
    if (JT["M"] != "map") vfail("satpipe viewer: " p ": manifest is not a JSON object")
    c = vj_req("M", "config")
    if (JT[c] != "map") vfail("satpipe viewer: manifest config is not an object")
    flat_write(c, "", "-")
}

# rows: steps 1-2 (viewer.py:61-89): MASTER.merged.bed; rows.tsv; candidate rows for the assignment and bits copies
function phase_rows(   line, seg, n, i, k, f, out, bed, s, rows, cid, sn, tmp, fa, fb, j, fam, b, bits) {
    vw_load_cfg()
    seg = ENVIRON["VW_RUN"] "/catalogue/MASTER.loci.segmented.bed"
    n = vw_read_tsv(seg, rows, 1)
    for (i = 1; i <= n; i++) if (split(rows[i], f, "\t") < 5)
        vfail("viewer: MASTER.loci.segmented.bed has no genome column (col 5); it was written by an older satpipe. Rerun from 14_arbitrate.")
    NR_ROWS = 0
    for (i = 1; i <= n; i++) {
        split(rows[i], f, "\t"); NR_ROWS++
        R_SCF[NR_ROWS] = f[1]; R_S[NR_ROWS] = vw_int_s(f[2]); R_E[NR_ROWS] = vw_int_s(f[3]); R_SP[NR_ROWS] = f[5]; R_LID[NR_ROWS] = f[4]
    }
    vw_load_cand()
    if (CAND) for (j = 1; j <= NCS; j++) {
        cid = CSORT[j]; s = CS SUBSEP cid; NR_ROWS++
        R_SCF[NR_ROWS] = py_str(vj_req(s, "scf")); R_S[NR_ROWS] = vw_int(vj_req(s, "start")); R_E[NR_ROWS] = vw_int(vj_req(s, "end"))
        R_SP[NR_ROWS] = py_str(vj_req(s, "sp")); R_LID[NR_ROWS] = cid
    }
    out = ENVIRON["VW_OUT"]; bed = out "/MASTER.merged.bed"; tmp = ENVIRON["VW_TMP"]
    printf "" > bed
    for (i = 1; i <= NR_ROWS; i++)
        printf "%s\t%d\t%d\t%s|%s:%d-%d|%s\n", R_SCF[i], R_S[i], R_E[i], R_SP[i], R_SCF[i], R_S[i] + 1, R_E[i], R_LID[i] > bed
    close(bed)
    printf "" > (tmp "/rows.tsv")
    for (i = 1; i <= NR_ROWS; i++)
        printf "%s\t%s\t%s\t%s\t%s\n", flat_esc(R_SCF[i]), R_S[i], R_E[i], flat_esc(R_SP[i]), flat_esc(R_LID[i]) > (tmp "/rows.tsv")
    close(tmp "/rows.tsv")
    if (!CAND) return
    fa = tmp "/cand_assign.tsv"; fb = tmp "/cand_bits.tsv"
    printf "" > fa; printf "" > fb
    for (j = 1; j <= NCS; j++) {
        cid = CSORT[j]; s = CS SUBSEP cid
        printf "%s\t%s\t%s\t%s\t%s\t%s\tCANDIDATE_UNCATALOGUED\t%s\n", cid, py_str(vj_req(s, "call")), py_str(vj_req(s, "call_bits")),
            (vj_truthy(vj_req(s, "comp")) ? py_str(s SUBSEP "comp") : ""), py_str(vj_req(s, "comp_bits")), py_str(vj_req(s, "delta")),
            py_str(vj_req(s, "call")) > fa
        bits = vj_req(s, "bits"); line = cid
        for (k = 1; k <= NFAM; k++) { b = vj_get(bits, FAM[k]); line = line "\t" (b != "" ? py_str(b) : "0") }
        print line > fb
    }
    close(fa); close(fb)
    # cand.get('summary') truthy (viewer.py:148)
    printf "%d\n", (JN[CS] > 0) > (tmp "/cand_summary_nonempty")
    close(tmp "/cand_summary_nonempty")
}

# fai lengths of every genome, in config order; (code, contig) -> int(col 2). FAI_N[i], FAI_NAME[i, j], FAI_LEN2[i, j] (raw col 2)
function vw_read_fais(   i, n, j, L, f, k) {
    for (i = 1; i <= NCODE; i++) {
        n = vw_lines(GPATH[i] ".fai", L); FAI_N[i] = n
        for (j = 1; j <= n; j++) {
            k = split(L[j], f, "\t")
            FAI_NAME[i, j] = f[1]
            FAI_K[i, j] = k
            # line.split('\t') keeps the newline in the last field: a two-column .fai gives 'len\n' as p[1]
            FAI_LEN2[i, j] = (k == 2 && VW_LT[j]) ? f[2] "\n" : f[2]
        }
    }
}

# flankbed: step 5 windows (viewer.py:104-121) -> $VW_FDIR/flanks.bed; prints repr(sorted(families)) for the signature (:125)
function phase_flankbed(   i, j, bed, a, b, L, key, srt, srt2, n, r) {
    vw_load_cfg(); vw_load_rows(); vw_read_fais()
    for (i = 1; i <= NCODE; i++) for (j = 1; j <= FAI_N[i]; j++) {
        if (FAI_K[i, j] < 2) vfail("satpipe viewer: IndexError: list index out of range")
        LEN[CODE[i] SUBSEP FAI_NAME[i, j]] = vw_int_s(FAI_LEN2[i, j])
    }
    bed = ENVIRON["VW_FDIR"] "/flanks.bed"; printf "" > bed
    for (i = 1; i <= NR_ROWS; i++) {
        key = R_SP[i] SUBSEP R_SCF[i]
        L = (key in LEN) ? LEN[key] + 0 : R_E[i] + 0
        a = R_S[i] - 5000; if (a < 0) a = 0; b = R_S[i] + 0
        if (b > a) printf "%s\t%d\t%d\t%s|%s:%d-%d|FLANK\n", R_SCF[i], a, b, R_SP[i], R_SCF[i], a + 1, b > bed
        a = R_E[i] + 0; b = R_E[i] + 5000; if (L < b) b = L
        if (b > a) printf "%s\t%d\t%d\t%s|%s:%d-%d|FLANK\n", R_SCF[i], a, b, R_SP[i], R_SCF[i], a + 1, b > bed
    }
    close(bed)
    for (i = 1; i <= NFAM; i++) srt[i] = FAM[i]
    n = asort(srt, srt2, "@val_str_asc"); r = "["
    for (i = 1; i <= n; i++) r = r (i > 1 ? ", " : "") py_repr_str(srt2[i])
    printf "%s", r "]"
}

# flankfam: the rows of one family's flank scan (viewer.py:141-146), appended to stdout
function phase_flankfam(   p, n, L, i, f, k, sp, x, fam) {
    p = ENVIRON["VW_PRE"] ".bits.tsv"; fam = ENVIRON["VW_FAM"]
    n = vw_lines(p, L)
    for (i = 2; i <= n; i++) {
        k = split(L[i], f, "\t")
        if (k >= 4) { sp = f[4]; x = index(sp, "|"); if (x) sp = substr(sp, 1, x - 1); printf "%s\t%s\t%s\t%s\t%s\n", f[1], f[2], f[3], fam, sp }
    }
}

# flankcand: candidate flank rows (viewer.py:151-160)
function phase_flankcand(   cf, n, rows, j, s, key, call, i, f, k, x, sp) {
    vw_load_cand()
    for (j = 1; j <= JN[CS]; j++) {
        s = CS SUBSEP JK[CS, j]
        key = py_str(vj_req(s, "sp")) "|" py_str(vj_req(s, "scf")) ":" (vw_int(vj_req(s, "start")) + 1) "-" vw_int(vj_req(s, "end"))
        CALL_OF[key] = py_str(vj_req(s, "call"))
    }
    cf = ENVIRON["VW_RUN"] "/candidates/CAND_FLANK.bits.tsv"
    n = vw_read_tsv(cf, rows, 4)
    for (i = 1; i <= n; i++) {
        split(rows[i], f, "\t")
        key = f[4]; if (match(key, /\|[^|]*$/)) key = substr(key, 1, RSTART - 1)
        if (!(key in CALL_OF)) vfail("viewer: candidate flank row " f[4] " matches no candidate")
        sp = f[4]; x = index(sp, "|"); if (x) sp = substr(sp, 1, x - 1)
        printf "%s\t%d\t%s\t%s\t%s\n", f[1], vw_int_s(f[2]) + 1, f[3], CALL_OF[key], sp
    }
}

# flanksjson: step 6 first half (viewer.py:165-174)
function phase_flanksjson(   j, cid, s, fl, fam0, mono, e0, s0, c, lf, rf, d, q, fl0, r) {
    vw_load_cfg(); vw_load_cand()
    vj_load(ENVIRON["VW_RUN"] "/ortho/flanks.json", "F")
    if (CAND) {
        fl = vj_req("K", "flanks")
        for (j = 1; j <= JN[fl]; j++) {
            cid = JK[fl, j]; s = vj_req(CS, cid)
            fam0 = py_str(vj_req(s, "call"))
            mono = vj_req(vj_req(vj_req("C", "families"), fam0), "mono")
            c = jt_set("F", cid, "map", "")
            vj_cp(vj_req(s, "sp"), c, "sp"); vj_cp(vj_req(s, "scf"), c, "scf")
            s0 = vj_req(s, "start"); e0 = vj_req(s, "end")
            vj_cp(s0, c, "start"); vj_cp(e0, c, "end")
            jt_set(c, "subtype", JT[vj_req(s, "call")], JV[s SUBSEP "call"])
            vj_cp(mono, c, "mono")
            if (JT[s0] == "int" && JT[e0] == "int") { d = JV[e0] - JV[s0]; jt_set(c, "length", "int", sprintf("%d", d)) }
            else { d = vw_float(e0) - vw_float(s0); jt_set(c, "length", "float", vj_float_text(d)) }
            # round((end - start) / mono): round half to even
            q = d / vw_float(mono); fl0 = int(q); if (fl0 > q) fl0--
            r = q - fl0
            if (r > 0.5 || (r == 0.5 && fl0 % 2 != 0)) fl0++
            jt_set(c, "n_monomers", "int", sprintf("%d", fl0))
            lf = vj_cp(vj_req(vj_req(fl, cid), "LF"), c, "LF"); jt_set(lf, "class", "str", "candidate")
            rf = vj_cp(vj_req(fl SUBSEP cid, "RF"), c, "RF"); jt_set(rf, "class", "str", "candidate")
        }
    }
    vj_write("F", ENVIRON["VW_OUT"] "/flanks.json")
}

# subfam: viewer.py:181-188 (called only when <out>/subfamily.json exists and candidates were loaded)
function phase_subfam(   p, j, cid, s, c, k, n, keys) {
    vw_load_cand()
    p = ENVIRON["VW_OUT"] "/subfamily.json"
    vj_load(p, "S")
    if (JT["S"] != "map") vfail("satpipe viewer: subfamily.json is not an object")
    split("pid bits second delta confident", keys, " ")
    for (j = 1; j <= JN[CS]; j++) {
        cid = JK[CS, j]; s = CS SUBSEP cid
        if (vj_get(s, "subfamily") == "") continue
        c = jt_set("S", cid, "map", "")
        vj_cp(vj_req(s, "call"), c, "major"); vj_cp(s SUBSEP "subfamily", c, "subfamily")
        for (k = 1; k <= 5; k++) { n = vj_get(s, "ss_" keys[k]); if (n != "") vj_cp(n, c, keys[k]); else jt_set(c, keys[k], "null", "None") }
    }
    vj_write("S", p)
}

# scaflen: step 7 scaffold_lengths.tsv (viewer.py:191-198)
function phase_scaflen(   f, i, j) {
    vw_load_cfg(); vw_load_rows()
    for (i = 1; i <= NR_ROWS; i++) WANT[R_SP[i] SUBSEP R_SCF[i]] = 1
    f = ENVIRON["VW_OUT"] "/scaffold_lengths.tsv"
    printf "scaffold\tlength\tspecies\n" > f
    vw_read_fais()
    for (i = 1; i <= NCODE; i++) for (j = 1; j <= FAI_N[i]; j++) {
        if (!((CODE[i] SUBSEP FAI_NAME[i, j]) in WANT)) continue
        if (FAI_K[i, j] < 2) vfail("satpipe viewer: IndexError: list index out of range")
        printf "%s\t%s\t%s\n", FAI_NAME[i, j], FAI_LEN2[i, j], CODE[i] > f
    }
    close(f)
}

# _read_fasta(path) values joined (viewer.py:33-40): records by first word of the header, a repeated name restarts its record
function vw_fasta_seq(p,   L, n, i, name, k, w, order, seq, out, s, nn) {
    n = vw_lines(p, L); name = ""; nn = 0
    for (i = 1; i <= n; i++) {
        if (substr(L[i], 1, 1) == ">") {
            s = substr(L[i], 2); sub(/^[ \t\n\r\v\f\034-\037]+/, "", s)
            if (s == "") vfail("satpipe viewer: IndexError: list index out of range")
            match(s, /^[^ \t\n\r\v\f\034-\037]+/); name = substr(s, 1, RLENGTH)
            if (!(name in seq)) order[++nn] = name
            seq[name] = ""
        } else if (name != "") {
            s = L[i]; gsub(/^[ \t\n\r\v\f\034-\037]+|[ \t\n\r\v\f\034-\037]+$/, "", s)
            seq[name] = seq[name] toupper(s)
        }
    }
    out = ""; for (k = 1; k <= nn; k++) out = out seq[order[k]]
    return out
}

# config: step 8 (viewer.py:205-261) -> viewer_config.json; prints "n_cat n_cand n_ortho n_flag" for the closing line
function phase_config(   out, fams, i, fam, fd, c, f, view, labels, t, sbk, cdir, mem, j, m, mlist, p, seq, half, h, o, n_cat, n_cand,
                          n_flag, n_contra, n_with, lists, notes, cnt, LO, LF, LU, best, bestv, x, sl, cand_txt, delta, dambig, mp, ma) {
    vw_load_cfg()
    out = ENVIRON["VW_OUT"]
    read_master(out "/MASTER.merged.bed")
    read_assign(out "/MASTER.final_assignment.clean.tsv")
    vj_load(out "/flanks.json", "F")
    vj_load(out "/flank_orthologs.json", "O")
    if (vw_isfile(out "/subfamily.json")) vj_load(out "/subfamily.json", "S"); else jt_init("S", "map")
    jt_init("V", "map")
    jt_set("V", "generator", "str", "satpipe " ENVIRON["VW_VERSION"])
    view = vj_get("C", "viewer"); if (!vj_truthy(view)) view = ""
    t = (view != "") ? vj_get(view, "title") : ""
    if (vj_truthy(t)) vj_cp(t, "V", "title"); else jt_set("V", "title", "str", py_str(vj_req("C", "run_name")) " locus viewer")
    vj_cp(vj_req("C", "run_name"), "V", "run_name")
    fams = jt_set("V", "families", "list", "")
    fd = vj_req("C", "families")
    best = ""
    for (i = 1; i <= NFAM; i++) {
        fam = FAM[i]; f = fd SUBSEP fam
        c = jt_set(fams, i - 1, "map", "")
        jt_set(c, "name", "str", fam)
        vj_cp(vj_req(f, "mono"), c, "mono")
        x = vj_get(f, "colour")
        if (vj_truthy(x)) vj_cp(x, c, "colour"); else jt_set(c, "colour", "str", PALETTE[(i - 1) % 10])
        jt_set(c, "consensus", "str", vw_fasta_seq(py_str(vj_req(f, "consensus"))))
        jt_set(c, "bits_file", "str", "MASTER_" fam ".bits.tsv")
        if (best == "" || vw_float(f SUBSEP "mono") > bestv) { best = f SUBSEP "mono"; bestv = vw_float(best) }
    }
    sl = jt_set("V", "species", "list", "")
    labels = (view != "") ? vj_get(view, "species_names") : ""; if (!vj_truthy(labels)) labels = ""
    for (i = 1; i <= NCODE; i++) {
        c = jt_set(sl, i - 1, "map", ""); jt_set(c, "code", "str", CODE[i])
        x = (labels != "") ? vj_get(labels, CODE[i]) : ""
        if (x != "") vj_cp(x, c, "label"); else jt_set(c, "label", "str", CODE[i])
    }
    c = jt_set("V", "subfamilies", "map", "")
    sbk = vj_get("C", "subfamily"); if (!vj_truthy(sbk)) sbk = ""
    cdir = (sbk != "") ? vj_get(sbk, "consensuses") : ""
    if (vj_truthy(cdir) && ENVIRON["VW_SUBDIR_OK"] == "1") {
        mem = vj_get(sbk, "members")
        if (vj_truthy(mem)) for (i = 1; i <= JN[mem]; i++) {
            fam = JK[mem, i]; mlist = mem SUBSEP fam
            for (j = 1; j <= JN[mlist]; j++) {
                m = py_str(mlist SUBSEP JK[mlist, j])
                p = (substr(m, 1, 1) == "/") ? m ".fa" : JV[cdir] (JV[cdir] ~ /\/$/ ? "" : "/") m ".fa"
                if (!vw_isfile(p)) continue
                seq = vw_fasta_seq(p); h = int(length(seq) / 2); half = substr(seq, 1, h)
                if (length(seq) % 2 == 0 && half == substr(seq, h + 1)) seq = half
                x = jt_set(c, m, "map", ""); jt_set(x, "family", "str", fam); jt_set(x, "seq", "str", seq)
            }
        }
    }
    vj_cp(best, "V", "window_step")
    vj_cp(vj_req(vj_req("C", "expansion"), "min_bits"), "V", "sat_bits")
    vj_cp(vj_req(vj_req("C", "finalisation"), "edge_tol"), "V", "edge_tol")
    jt_set("V", "flank_view", "int", "5000")
    o = vj_req("C", "orthology")
    c = jt_set("V", "orthology", "map", "")
    vj_cp(vj_req(o, "min_pid"), c, "min_pid"); vj_cp(vj_req(o, "min_aln"), c, "min_aln"); vj_cp(vj_req(o, "flank_bp"), c, "flank_bp")
    # counts and lists keep their dict positions; the lists are built in place, the counts filled from them afterwards
    cnt = jt_set("V", "counts", "map", "")
    lists = jt_set("V", "lists", "map", "")
    LO = jt_set(lists, "ortho", "list", ""); LF = jt_set(lists, "flagged", "list", ""); LU = jt_set(lists, "full", "list", "")
    build_lists(LO, LF, LU)
    n_cat = 0; n_cand = 0
    for (i = 1; i <= M_N; i++) { if (substr(M_ID[i], 1, 1) == "C") n_cand++; else n_cat++ }
    n_flag = 0; n_contra = 0
    for (i = 0; i < JN[LF]; i++) {
        x = LF SUBSEP i
        if (substr(JV[x SUBSEP "id"], 1, 1) != "C") n_flag++
        if (JV[x SUBSEP "cat"] == "UNIQUE_SEED_CONTRADICTION") n_contra++
    }
    n_with = 0; for (i = 0; i < JN[LO]; i++) if (JV[LO SUBSEP i SUBSEP "n"] + 0 > 0) n_with++
    jt_set(cnt, "catalogue_loci", "int", n_cat); jt_set(cnt, "candidates", "int", n_cand); jt_set(cnt, "flagged", "int", n_flag)
    jt_set(cnt, "contradictions", "int", n_contra); jt_set(cnt, "anchorable", "int", JN[LO]); jt_set(cnt, "with_ortholog", "int", n_with)
    mp = vw_fmt_g(vj_req(o, "min_pid")); ma = vw_fmt_d(vj_req(o, "min_aln"))
    x = vj_req("C", "arbitration"); delta = vw_fmt_g(vj_req(x, "delta")); dambig = vw_fmt_g(vj_req(x, "delta_ambig"))
    cand_txt = n_cand ? ", plus the " n_cand " uncatalogued candidates (not catalogue members)" : ""
    notes = jt_set("V", "notes", "map", "")
    jt_set(notes, "ortho", "str", "Positional orthology from unique flanking sequence (>=" mp "% identity, >=" ma " bp), single best hit per target " \
        "assembly. \"Orthologous to\" lists distinct partner assemblies then the total partner-locus count.")
    jt_set(notes, "flag", "str", n_flag " catalogue loci annotated by arbitration" cand_txt ". Only " n_contra " " (n_contra == 1 ? "is" : "are") \
        " a seed/rescan contradiction (the rescan scored another family >=" delta " bits above the seed; the seed family is retained). " \
        "The rest are close-margin calls (segment or seed-conflict resolved with a margin below " dambig " bits) or loci without a seed.")
    jt_set(notes, "full", "str", "Every locus (" n_cat ", all " NCODE " assemblies) with its sub-subfamily call and both flanks' cross-assembly orthology " \
        "status (ssearch36, >=" mp "% identity, >=" ma " bp, single best hit per target assembly). L/R badges show each flank's class.")
    vj_write("V", out "/viewer_config.json")
    printf "%d %d %d %d\n", n_cat, n_cand, JN[LO], n_flag
}

BEGIN {
    split("#60a5fa #34d399 #f472b6 #fbbf24 #a78bfa #f87171 #2dd4bf #fb923c #c084fc #4ade80", _vp, " ")
    for (_vi = 1; _vi <= 10; _vi++) PALETTE[_vi - 1] = _vp[_vi]
    _ph = ENVIRON["VW_PHASE"]
    if (_ph == "manifest") phase_manifest()
    else if (_ph == "rows") phase_rows()
    else if (_ph == "flankbed") phase_flankbed()
    else if (_ph == "flankfam") phase_flankfam()
    else if (_ph == "flankcand") phase_flankcand()
    else if (_ph == "flanksjson") phase_flanksjson()
    else if (_ph == "subfam") phase_subfam()
    else if (_ph == "scaflen") phase_scaflen()
    else if (_ph == "config") phase_config()
    else vfail("viewer.awk: unknown phase " _ph)
    exit 0
}
