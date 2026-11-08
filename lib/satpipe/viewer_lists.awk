# viewer_lists.awk - satpipe/viewer_lists.py in gawk: the viewer's three inspection lists (ORTHO_LOCI, FLAGGED_LOCI,
# FULL_TABLE) from the viewer data files. Library for viewer.awk (needs pyjson.awk's typed tree and viewer.awk's
# vj_* helpers, which are loaded in the same gawk run).
#
#   read_master(path)    viewer_lists.py:16-28   M_N; M_ID[i] M_SP[i] (catalogue order); SP_OF[lid] (last row wins)
#   read_assign(path)    viewer_lists.py:31-39   A_FAM[lid] A_BITS A_OTHER A_OB A_DELTA A_FLAG (last row wins)
#   build_lists(LO, LF, LU)   viewer_lists.py:54-142 into three existing empty list nodes; reads the trees
#                        "F" (flanks.json), "O" (flank_orthologs.json), "S" (subfamily.json or {}); VL_PAIRS = path of the
#                        decided ortho_pairs.tsv (pairs=read_pairs(..)) or "" (pairs=None, the published viewer)
# The threshold key is '80' (build_lists default; viewer.py:208 passes none).

function read_master(path,   L, n, i, k, p, parts, np, cnt) {
    n = vw_lines(path, L); M_N = 0; cnt = 0
    for (i = 1; i <= n; i++) {
        k = split(L[i], p, "\t")
        if (k < 4) continue
        cnt++
        np = split(p[4], parts, "|")
        M_N++
        M_ID[M_N] = (np > 2) ? parts[3] : "L" cnt
        M_SP[M_N] = parts[1]
        vw_int_s(p[2]); vw_int_s(p[3])                 # int(p[1]), int(p[2]) must parse
        SP_OF[M_ID[M_N]] = parts[1]
    }
}

function read_assign(path,   L, n, i, k, p, id) {
    n = vw_lines(path, L)
    for (i = 1; i <= n; i++) {
        k = split(L[i], p, "\t")
        if (k < 7 || tolower(p[1]) == "locus_id" || tolower(p[1]) == "locus") continue
        id = p[1]
        A_FAM[id] = p[2]; A_BITS[id] = p[3]; A_OTHER[id] = p[4]; A_OB[id] = p[5]; A_DELTA[id] = p[6]; A_FLAG[id] = p[7]
    }
}

# _fmt (viewer_lists.py:42-46): '%.1f' % float(x), or str(x) when float() fails
function vl_fmt(id,   t, v, s) {
    t = JT[id]; v = JV[id]
    if (t == "float" && v ~ /^-?(nan|inf)$/) return v
    if (t == "int" || t == "float") return sprintf("%.1f", v + 0)
    if (t == "bool") return (v == "True") ? "1.0" : "0.0"
    if (t == "str" || t == "path") {
        s = v; gsub(/^[ \t\n\r\v\f]+|[ \t\n\r\v\f]+$/, "", s)
        if (s ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$/) return sprintf("%.1f", s + 0)
        if (tolower(s) ~ /^[-+]?(inf|infinity)$/) return (substr(s, 1, 1) == "-") ? "-inf" : "inf"
        if (tolower(s) ~ /^[-+]?nan$/) return "nan"
        return v
    }
    return py_str(id)
}

# index of the best hit in a hit list: highest bits, then higher pid, then smallest locus id; first on full ties
function vl_best(hl,   i, h, b, p, l, bb, bp, bl, best) {
    best = ""
    for (i = 0; i < JN[hl]; i++) {
        h = hl SUBSEP i
        b = (vj_get(h, "bits") != "") ? vw_float(h SUBSEP "bits") : 0
        p = (vj_get(h, "pid") != "") ? vw_float(h SUBSEP "pid") : 0
        l = py_str(vj_req(h, "locus")) ""
        if (best == "" || b > bb || (b == bb && (p > bp || (p == bp && l < bl)))) { best = h; bb = b; bp = p; bl = l }
    }
    return best
}

function _vl_pair(a, b) { return ((a "") < (b "")) ? a SUBSEP b : b SUBSEP a }

function _vl_add_partner(a, b) { if (!((a, b) in PSET)) { PSET[a, b] = 1; PL[a, ++PN[a]] = b } }

function build_lists(LO, LF, LU,   FLAGGED, key, mask, tab, q, i, hl, j, h, t, lid, sp, fam, r, fl, x, lfc, rfc, lfcid, rfcid,
                     sf, sbn, side, sid, b, n, judged, nj, all, k, both, nuf, bysp, spk, sps, nsp, lst, lst2, s, pl, cnt, fg, nfg, prow, fln,
                     PLN, npl, pf) {
    split("UNIQUE_SEED_CONTRADICTION SEED_CONFLICT_RESOLVED_AMBIGUOUS SEGMENT_RESOLVED_AMBIGUOUS NO_SEED_RESOLVED NO_SEED_AMBIGUOUS CANDIDATE_UNCATALOGUED", fln, " ")
    for (i in fln) FLAGGED[fln[i]] = 1
    # LF / RF hit tables at the threshold
    for (side = 1; side <= 2; side++) {
        key = (side == 1) ? "LF" : "RF"
        t = vj_get("O", key); tab = ""
        if (vj_truthy(t)) { tab = vj_get(t, "80") }
        TAB[side] = tab
        if (tab == "" || JT[tab] != "map") continue
        mask = side
        for (i = 1; i <= JN[tab]; i++) {
            q = JK[tab, i]; hl = tab SUBSEP q
            for (j = 0; j < JN[hl]; j++) {
                h = hl SUBSEP j; t = py_str(vj_req(h, "locus"))
                _vl_add_partner(q, t); _vl_add_partner(t, q)
                k = _vl_pair(q, t); SIDES[k] = or(SIDES[k] + 0, mask)
            }
        }
    }
    # decided pairs (viewer_lists.py read_pairs; VL_PAIRS = path of ortho_pairs.tsv when stage 19 decided them): a locus's
    # partners are exactly these pairs; the hit tables above stay as evidence (lfp/rfp columns, both-sided sides)
    if (VL_PAIRS != "") {
        delete PSET; delete PN; delete PL
        npl = vw_lines(VL_PAIRS, PLN)
        for (i = 1; i <= npl; i++) {
            if (split(PLN[i], pf, "	") < 2 || pf[1] == "locusA") continue
            _vl_add_partner(pf[1], pf[2]); _vl_add_partner(pf[2], pf[1])
        }
    }
    nfg = 0
    for (i = 1; i <= M_N; i++) {
        lid = M_ID[i]; sp = M_SP[i]
        if (!(lid in A_FAM)) continue
        fam = A_FAM[lid]
        if (A_FLAG[lid] in FLAGGED) {
            nfg++; FG_ID[nfg] = lid; FG_SP[nfg] = sp; FG_FAM[nfg] = fam
        }
        if (substr(lid, 1, 1) == "C") continue
        fl = vj_get("F", lid); if (!vj_truthy(fl)) fl = ""
        lfcid = ""; rfcid = ""
        if (fl != "") {
            x = vj_get(fl, "LF"); if (vj_truthy(x)) { x = vj_get(x, "class"); if (vj_truthy(x)) lfcid = x }
            x = vj_get(fl, "RF"); if (vj_truthy(x)) { x = vj_get(x, "class"); if (vj_truthy(x)) rfcid = x }
        }
        r = jt_set(LU, JN[LU], "map", "")
        jt_set(r, "id", "str", lid); jt_set(r, "sp", "str", sp); jt_set(r, "fam", "str", fam)
        sf = vj_get("S", lid); if (!vj_truthy(sf)) sf = ""
        sbn = (sf != "") ? vj_get(sf, "subfamily") : ""
        if (sbn != "") vj_cp(sbn, r, "sub"); else jt_set(r, "sub", "null", "None")
        if (vj_truthy(sbn) && vj_truthy(vj_get(sf, "confident"))) jt_set(r, "subc", "bool", "True")
        if (lfcid != "") vj_cp(lfcid, r, "lfc"); else jt_set(r, "lfc", "str", "none")
        if (rfcid != "") vj_cp(rfcid, r, "rfc"); else jt_set(r, "rfc", "str", "none")
        for (side = 1; side <= 2; side++) {
            sid = (side == 1) ? "lf" : "rf"
            hl = (TAB[side] != "" && JT[TAB[side]] == "map") ? vj_get(TAB[side], lid) : ""
            if (!vj_truthy(hl)) continue
            b = vl_best(hl)
            vj_cp(b SUBSEP "locus", r, sid "p")
            vj_cp(vj_req(b, "sp"), r, sid "sp")
            x = py_str(b SUBSEP "locus")
            jt_set(r, sid "fam", "str", (x in A_FAM) ? A_FAM[x] : "")
            jt_set(r, sid "pid", "str", vl_fmt(vj_req(b, "pid")))
            jt_set(r, sid "aln", "str", py_str(vj_req(b, "aln")))
        }
        n = PN[lid] + 0
        if (n > 0) {
            jt_set(r, "ortho", "bool", "True")
            nj = 0; delete judged
            if (VL_PAIRS == "") for (k = 1; k <= 2; k++) { x = vj_get(r, (k == 1) ? "lfp" : "rfp"); if (x != "") judged[++nj] = py_str(x) }
            if (nj == 0) for (k = 1; k <= n; k++) judged[++nj] = PL[lid, k]
            all = 1
            for (k = 1; k <= nj; k++) if (!((judged[k]) in A_FAM) || A_FAM[judged[k]] != fam) { all = 0; break }
            jt_set(r, "pfam", "str", all ? "SAME" : "DIFFERENT")
        }
        both = 0
        for (k = 1; k <= n; k++) if (SIDES[_vl_pair(lid, PL[lid, k])] == 3) { both = 1; break }
        if (both) jt_set(r, "both", "bool", "True")
        nuf = (lfcid != "" && JT[lfcid] == "str" && JV[lfcid] == "unique") + (rfcid != "" && JT[rfcid] == "str" && JV[rfcid] == "unique")
        if (nuf >= 1 || n > 0) {
            delete bysp; delete sps; nsp = 0
            for (k = 1; k <= n; k++) {
                x = PL[lid, k]
                if (!(x in SP_OF))
                    vfail("satpipe viewer: ValueError: ortholog partner " x " of " lid " is not in the catalogue (flank_orthologs.json and " \
                          "MASTER.merged.bed are from different runs?)")
                s = SP_OF[x]
                if (!(s in bysp)) { sps[++nsp] = s; bysp[s] = 0 }
                BY[s, ++bysp[s]] = x
            }
            prow = jt_set(LO, JN[LO], "map", "")
            jt_set(prow, "id", "str", lid); jt_set(prow, "sp", "str", sp); jt_set(prow, "fam", "str", fam)
            jt_set(prow, "nuf", "int", nuf); jt_set(prow, "n", "int", n)
            pl = jt_set(prow, "partners", "map", "")
            if (nsp > 0) {
                asort(sps, spk, "@val_str_asc")
                for (k = 1; k <= nsp; k++) {
                    s = spk[k]; delete lst
                    for (j = 1; j <= bysp[s]; j++) lst[j] = BY[s, j]
                    cnt = asort(lst, lst2, "@val_str_asc")
                    x = jt_set(pl, s, "list", "")
                    for (j = 1; j <= cnt; j++) jt_set(x, j - 1, "str", lst2[j])
                }
            }
        }
    }
    # flagged: uncatalogued candidates first, then catalogue loci, each in catalogue order (stable sort)
    for (k = 0; k <= 1; k++) for (i = 1; i <= nfg; i++) {
        if ((substr(FG_ID[i], 1, 1) == "C") != (k == 0)) continue
        lid = FG_ID[i]; r = jt_set(LF, JN[LF], "map", "")
        jt_set(r, "id", "str", lid); jt_set(r, "sp", "str", FG_SP[i]); jt_set(r, "cat", "str", A_FLAG[lid])
        jt_set(r, "f1", "str", FG_FAM[i]); jt_set(r, "s1", "str", A_BITS[lid]); jt_set(r, "f2", "str", A_OTHER[lid])
        jt_set(r, "s2", "str", A_OB[lid]); jt_set(r, "d", "str", A_DELTA[lid]); jt_set(r, "fin", "str", FG_FAM[i])
    }
}
