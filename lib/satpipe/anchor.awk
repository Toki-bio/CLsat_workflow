# anchor.awk - anchor phasing (satpipe/anchor.py) and the bodies of stage 21_anchor (run.py:1039-1096) and
# `satpipe anchor` (cli.py:91-125), PORT_SPEC row 21_anchor. Load AFTER pyjson.awk (flat_esc/flat_unesc, json_str):
#   gawk -f pyjson.awk -f anchor.awk -e 'BEGIN { anchor_main() }'          (anchor.sh does this)
# Needs LC_ALL=C: sequences are handled as bytes; toupper() is ASCII-only, which equals str.upper() for ASCII input.
# Scope: ASCII FASTA. Python's str.upper()/split()/strip() on non-ASCII text (e.g. U+00A0 as whitespace) is not reproduced.
#
# Error model: a function that would raise in Python sets ANC_ERRT (exception class) and ANC_ERR (str(e)) and returns "";
# ANC_ERRT == "" means success. Only AnchorError is caught where Python catches it (anchor.py:167, :190).
#
# Library API (each mirrors the anchor.py function of the same name; file:line of every tie/edge rule is cited inline)
#   anc_clean(s)                                       anchor.py:53-54
#   anc_local_align(q, t)                              anchor.py:57-97   -> ANC_AL[1..5] = score, q_start, t_start, cols, matches
#   anc_hamming(seq, motif)                            anchor.py:100-101
#   anc_frame(anchor_seq, motif, max_mm)               anchor.py:104-119
#   anc_rotate(seq, anchor, mi, min_aligned, motif, max_mm, max_snap)   anchor.py:122-151 -> rot; ANC_SHIFT ANC_IDENT ANC_COLS ANC_SNAP
#   anc_check_frame(seq, motif, max_mm)                anchor.py:154-157
#   anc_records(R, A, anchor, motif, max_mm, mi)       anchor.py:160-179 (R: input records, A: output arrays, see below)
#   anc_verify(R, P, anchor, motif, max_mm, mi)        anchor.py:182-195 -> P[1..P["n"]] problem strings
#   anc_read_fasta(path, R)                            anchor.py:23-39
# Record arrays: R["n"], R["name",k], R["seq",k], R["idx",name] = k (dict with insertion order).

function anc_raise(t, m) { ANC_ERRT = t; ANC_ERR = m; return "" }

# Python's % for ints: result has the sign of the divisor (anchor.py:138, :144). n == 0 -> ZeroDivisionError.
function anc_pymod(a, n,   r) {
    if (n == 0) { anc_raise("ZeroDivisionError", "integer division or modulo by zero"); return 0 }
    r = a % n
    return (r < 0) ? r + n : r
}

# anchor.py:53-54: remove '-' and ' ', upper-case, U -> T (in this order; nothing else is removed, N and '.' stay)
function anc_clean(s) { gsub(/[- ]/, "", s); s = toupper(s); gsub(/U/, "T", s); return s }

# anchor.py:57-97: Smith-Waterman, MATCH 2, MISMATCH -3, GAP -5 (anchor.py:16), linear gaps.
# Cell score s = max(0, d, u, l) (:73); trace 0 when s == 0, else diag if s == d, else up if s == u, else left (:75-82);
# best cell only replaced on a strictly larger score (:83), rows outer, columns inner -> first maximum in row-major order;
# traceback until a 0 trace or an edge (:87); returned starts are the stop cell (:97).
function anc_local_align(q, t,    n, m, Q, T, H, tr, i, j, W, best, bi, bj, d, u, l, s, qi, k, cols, mat, b0, b1, base) {
    n = length(q); m = length(t)
    split(q, Q, ""); split(t, T, "")
    W = m + 1
    for (j = 0; j <= m; j++) { H[j] = 0; H[W + j] = 0 }
    best = 0; bi = 0; bj = 0
    for (i = 1; i <= n; i++) {
        b1 = (i % 2) * W; b0 = W - b1          # current / previous row in the two-row buffer
        H[b1] = 0; qi = Q[i]; base = i * W
        for (j = 1; j <= m; j++) {
            d = H[b0 + j - 1] + (qi == T[j] ? 2 : -3)
            u = H[b0 + j] - 5
            l = H[b1 + j - 1] - 5
            s = d; if (u > s) s = u; if (l > s) s = l; if (s < 0) s = 0
            H[b1 + j] = s
            if (s != 0) tr[base + j] = (s == d) ? 1 : (s == u) ? 2 : 3     # absent = 0 (stop)
            if (s > best) { best = s; bi = i; bj = j }
        }
    }
    i = bi; j = bj; cols = 0; mat = 0
    while (i > 0 && j > 0 && ((i * W + j) in tr)) {
        k = tr[i * W + j]; cols++
        if (k == 1) { if (Q[i] == T[j]) mat++; i--; j-- }
        else if (k == 2) i--
        else j--
    }
    ANC_AL[1] = best; ANC_AL[2] = i; ANC_AL[3] = j; ANC_AL[4] = cols; ANC_AL[5] = mat
    return best
}

# anchor.py:100-101: mismatches over the common prefix + the part of the motif beyond the sequence
function anc_hamming(seq, motif,   ls, lm, n, i, c) {
    ls = length(seq); lm = length(motif); n = (ls < lm) ? ls : lm; c = 0
    for (i = 1; i <= n; i++) if (substr(seq, i, 1) != substr(motif, i, 1)) c++
    if (lm > ls) c += lm - ls
    return c
}

# anchor.py:104-119. cands = sorted((mm, i)) (:112): cands[0] = lowest mm, lowest i among those; the ambiguity test
# (:116) fires when a second position has the same lowest mm, and reports the next-lowest such i (cands[1][1]).
# Empty anchor with a motif: cands is empty -> cands[0] raises IndexError (:113).
function anc_frame(anchor_seq, motif, max_mm,   a, d, la, lm, i, mm, best, bpos, spos) {
    ANC_ERRT = ""
    a = anc_clean(anchor_seq)
    if (motif == "") return a                                   # :108-109
    motif = toupper(motif)                                      # :110
    d = a a; la = length(a); lm = length(motif)
    if (la == 0) return anc_raise("IndexError", "list index out of range")
    best = -1; spos = -1
    for (i = 0; i < la; i++) {
        mm = anc_hamming(substr(d, i + 1, lm), motif)
        if (best < 0 || mm < best) { best = mm; bpos = i; spos = -1 }
        else if (mm == best && spos < 0) spos = i
    }
    if (best > max_mm + 0)                                      # :114
        return anc_raise("AnchorError", sprintf("anchor motif %s not found in the anchor consensus (best %d mismatches)", motif, best))
    if (spos >= 0)                                              # :116
        return anc_raise("AnchorError", sprintf("anchor motif %s occurs more than once in the anchor consensus at equal quality " \
                                                "(positions %d and %d): the frame is ambiguous", motif, bpos, spos))
    return substr(d, bpos + 1, la)                              # :119
}

# anchor.py:122-151. The anchor is used as given (not cleaned, :132); threshold min(min_aligned, len(anchor)) (:134);
# identity mat/cols, 0.0 without columns (:133); shift (ts - qs) mod n (:138, Python modulo); snap only when the
# alignment start does not carry the motif (:139); snap candidates in p order, kept on a strictly smaller (mm, dist) (:145)
# -> lowest mm, then nearest, then lowest p; dist = circular distance (:144); dist <= max_snap (:145).
# Note the motif in the not-within message is the motif as passed, not upper-cased (:148-149).
function anc_rotate(seq, anchor, mi, min_aligned, motif, max_mm, max_snap,   s, d, cols, mat, ident, thr, n, shift, snapped, M, lm, p, mm, dist, a1, a2, found, bmm, bdist, bp) {
    ANC_ERRT = ""
    s = anc_clean(seq); d = s s
    anc_local_align(anchor, d)
    cols = ANC_AL[4]; mat = ANC_AL[5]
    ident = cols ? mat / cols : 0.0
    thr = (min_aligned + 0 < length(anchor)) ? min_aligned + 0 : length(anchor)
    if (cols < thr || ident < mi + 0)
        return anc_raise("AnchorError", sprintf("cannot place the sequence on the anchor (aligned %d columns, %.0f%% identity)", cols, 100 * ident))
    n = length(s)
    shift = anc_pymod(ANC_AL[3] - ANC_AL[2], n); if (ANC_ERRT != "") return ""
    snapped = 0
    if (motif != "") {
        M = toupper(motif); lm = length(M)
        if (anc_hamming(substr(d, shift + 1, lm), M) > max_mm + 0) {
            found = 0
            for (p = 0; p < n; p++) {
                mm = anc_hamming(substr(d, p + 1, lm), M)
                if (mm <= max_mm + 0) {
                    a1 = anc_pymod(p - shift, n); a2 = anc_pymod(shift - p, n); dist = (a1 < a2) ? a1 : a2
                    if (dist <= max_snap + 0 && (!found || mm < bmm || (mm == bmm && dist < bdist))) { found = 1; bmm = mm; bdist = dist; bp = p }
                }
            }
            if (!found)
                return anc_raise("AnchorError", sprintf("the motif %s is not within %d bases of the alignment start", motif, max_snap))
            shift = bp; snapped = bdist
        }
    }
    ANC_SHIFT = shift; ANC_IDENT = ident; ANC_COLS = cols; ANC_SNAP = snapped
    return substr(d, shift + 1, n)
}

function anc_check_frame(seq, motif, max_mm) {                  # anchor.py:154-157
    if (motif == "") return 1
    return anc_hamming(anc_clean(seq), toupper(motif)) <= max_mm + 0
}

function anc_maxlen(motif) { return (length(motif) > 7) ? length(motif) : 7 }     # max(7, len(motif)), anchor.py:175, :194

# dict assignment with Python semantics: a new key is appended, an existing key keeps its position (sub_out.update, run.py:1076)
function anc_put(A, name, seq) {
    if (!(("idx", name) in A)) { A["n"]++; A["idx", name] = A["n"]; A["name", A["n"]] = name }
    A["seq", A["idx", name]] = seq
}
function anc_clear(A) { delete A; A["n"] = 0 }

# anchor.py:160-179. A gets the rotated records (also motif-mismatch ones, :174) and A["rows"], A["row",k] = the report
# tuple as tab-joined str() values (:169, :175-176). On failure ANC_ERRT = AnchorError with every bad record (:177-178);
# a non-AnchorError from rotate propagates at once (only AnchorError is caught, :167).
function anc_records(R, A, anchor, motif, max_mm, mi,   k, name, seq, rot, ok, nbad, bad, st) {
    anc_clear(A); A["rows"] = 0; nbad = 0
    for (k = 1; k <= R["n"]; k++) {
        name = R["name", k]; seq = R["seq", k]
        rot = anc_rotate(seq, anchor, mi, 40, motif, max_mm, 40)
        if (ANC_ERRT == "AnchorError") {
            bad[++nbad] = name ": " ANC_ERR
            A["row", ++A["rows"]] = name "\t" length(anc_clean(seq)) "\t\t\t\tunplaced"
            continue
        }
        if (ANC_ERRT != "") return ""
        ok = anc_check_frame(rot, motif, max_mm)
        if (!ok) bad[++nbad] = name ": after rotation the sequence starts " substr(rot, 1, length(motif)) ", not the anchor motif " motif
        anc_put(A, name, rot)
        st = (ok ? "ok" : "motif-mismatch") (ANC_SNAP ? "; snapped " ANC_SNAP " to motif" : "")
        A["row", ++A["rows"]] = name "\t" length(rot) "\t" ANC_SHIFT "\t" sprintf("%.1f", 100 * ANC_IDENT) "\t" substr(rot, 1, anc_maxlen(motif)) "\t" st
    }
    ANC_ERRT = ""
    if (nbad) {
        ANC_ERR = "anchor phasing failed:\n  " bad[1]
        for (k = 2; k <= nbad; k++) ANC_ERR = ANC_ERR "\n  " bad[k]
        ANC_ERRT = "AnchorError"; return ""
    }
    return 1
}

# anchor.py:182-195
function anc_verify(R, P, anchor, motif, max_mm, mi,   k, name, s) {
    delete P; P["n"] = 0
    for (k = 1; k <= R["n"]; k++) {
        name = R["name", k]; s = anc_clean(R["seq", k])
        anc_rotate(s, anchor, mi, 40, motif, max_mm, 40)
        if (ANC_ERRT == "AnchorError") { P[++P["n"]] = name ": " ANC_ERR; continue }
        if (ANC_ERRT != "") return ""
        if (ANC_SHIFT != 0 || !anc_check_frame(s, motif, max_mm))
            P[++P["n"]] = name ": out of the anchor frame (rotation needed: " ANC_SHIFT " bases; starts " substr(s, 1, anc_maxlen(motif)) ")"
    }
    ANC_ERRT = ""
    return 1
}

# anchor.py:23-39. Text mode = universal newlines (\r\n, \r, \n end a line); header name = first whitespace-separated
# token of line[1:] (:32; no token -> IndexError); sequence lines stripped (:36); lines before the first header ignored
# (:35); a name seen twice -> ValueError '<path>: record name <n> appears twice' (:33-34). Unreadable -> OSError.
function anc_read_fasta(path, R,   ors, line, name, buf, have, rc, t) {
    delete R; R["n"] = 0; ANC_ERRT = ""
    ors = RS; RS = "\r\n|\r|\n"; have = 0
    while ((rc = (getline line < path)) > 0) {
        if (substr(line, 1, 1) == ">") {
            if (have) anc_put(R, name, buf)
            t = substr(line, 2); sub(/^[ \t\n\v\f\r\034\035\036\037]+/, "", t)
            if (t == "") { close(path); RS = ors; return anc_raise("IndexError", "list index out of range") }
            match(t, /^[^ \t\n\v\f\r\034\035\036\037]+/); name = substr(t, 1, RLENGTH); buf = ""; have = 1
            if (("idx", name) in R) { close(path); RS = ors; return anc_raise("ValueError", path ": record name " name " appears twice") }
        } else if (have) {
            gsub(/^[ \t\n\v\f\r\034\035\036\037]+|[ \t\n\v\f\r\034\035\036\037]+$/, "", line); buf = buf line
        }
    }
    close(path); RS = ors
    if (rc < 0) return anc_raise("OSError", "cannot read " path)
    if (have) anc_put(R, name, buf)
    return 1
}

function anc_join_values(R,   k, s) { s = ""; for (k = 1; k <= R["n"]; k++) s = s R["seq", k]; return s }   # ''.join(values())

function anc_write_fasta(path, A,   k) {                     # anchor.py:42-50, width 0
    printf "" > path
    for (k = 1; k <= A["n"]; k++) printf ">%s\n%s\n", A["name", k], A["seq", A["idx", A["name", k]]] > path
    close(path)
}

# ---------------------------------------------------------------- job driver (anchor.sh writes the job file)
# Job file (ENVIRON["ANC_JOB"]): one row per line, fields TAB-separated and flat_esc-escaped:
#   set KEY VALUE      KEY in: mode (stage|cli), afam, motif, mm, mi, indep (True/False), outdir, check (0/1)
#   fam NAME PATH      families in config order
#   cons PATH EXISTS   stage: subfamily.consensuses as configured, EXISTS 0/1 (os.path.exists, run.py:1058)
#   subfile PATH       stage: fasta_files(cons) in order (anchor.py:198-204, expanded by anchor.sh)
#   member FAM NAME    stage: subfamily.members in config order
#   input FP OUT       cli: each FASTA file in order and its output path (cli.py:119)
# Result (stdout): stage -> "MANIFEST\t<json>" and "SAY\t<text>"; errors -> "ERROR\t<class>\t<escaped message>" and exit 3.
# cli -> its own stdout/stderr lines exactly as cli.py prints them; exit 0/1.

function _anc_fail(t, m) { printf "ERROR\t%s\t%s\n", t, flat_esc(m); exit 3 }

function anchor_main(   job, line, f, nf, S, nfam, FN, FP, cons, cons_ex, nsub, SF, nmem, MF, MN, k, afam, frame, R, FR, SUB,
                        FO, SO, rows, ROW, r, i, part, PA, OUTR, j, out, n, nin, IN_FP, IN_OUT, P, bad, gapped, man, rotated, indep) {
    job = ENVIRON["ANC_JOB"]
    while ((getline line < job) > 0) {
        nf = split(line, f, "\t"); for (i = 1; i <= nf; i++) f[i] = flat_unesc(f[i])
        if (f[1] == "set") S[f[2]] = f[3]
        else if (f[1] == "fam") { FN[++nfam] = f[2]; FP[nfam] = f[3] }
        else if (f[1] == "cons") { cons = f[2]; cons_ex = f[3] }
        else if (f[1] == "subfile") SF[++nsub] = f[2]
        else if (f[1] == "member") { MF[++nmem] = f[2]; MN[nmem] = f[3] }
        else if (f[1] == "input") { IN_FP[++nin] = f[2]; IN_OUT[nin] = f[3] }
    }
    close(job)
    afam = (S["afam"] != "") ? S["afam"] : FN[1]                                  # run.py:1046 / cli.py:100
    if (S["mode"] == "cli") exit anchor_cli(FN, FP, nfam, afam, S, IN_FP, IN_OUT, nin)

    # ---- stage 21_anchor (run.py:1040-1096)
    anc_clear(FR)                                                               # fam_recs (:1048)
    for (k = 1; k <= nfam; k++) {
        if (!anc_read_fasta(FP[k], R)) _anc_fail(ANC_ERRT, ANC_ERR)
        anc_put(FR, FN[k], anc_join_values(R))
    }
    if (!(("idx", afam) in FR)) _anc_fail("KeyError", "'" afam "'")
    frame = anc_frame(FR["seq", FR["idx", afam]], S["motif"], S["mm"])          # :1049 (also in independent mode)
    if (ANC_ERRT != "") _anc_fail(ANC_ERRT, ANC_ERR)
    anc_clear(SUB)
    if (cons != "") {                                                           # :1055-1067
        if (!cons_ex) _anc_fail("StageError", "21_anchor: subfamily.consensuses not found: " cons)
        for (k = 1; k <= nsub; k++) {
            if (!anc_read_fasta(SF[k], R)) _anc_fail(ANC_ERRT, ANC_ERR)
            for (j = 1; j <= R["n"]; j++) {
                if (("idx", R["name", j]) in SUB) _anc_fail("StageError", "21_anchor: sub-subfamily record " R["name", j] " appears twice (" SF[k] ")")
                anc_put(SUB, R["name", j], R["seq", j])
            }
        }
        part = ""
        for (k = 1; k <= nmem; k++) if (!(("idx", MN[k]) in SUB)) part = part (part == "" ? "" : ", ") MN[k]
        if (part != "") _anc_fail("StageError", "21_anchor: sub-subfamilies listed in subfamily.members without a consensus record: " part)
    } else if (nmem && S["indep"] == "True") _anc_fail("KeyError", "'" MN[1] "'")   # sub[n] on an empty dict (:1073); crosscheck forbids it
    rows = 0; anc_clear(FO); anc_clear(SO); indep = (S["indep"] == "True")
    if (indep) {                                                                # :1069-1077
        for (k = 1; k <= nfam; k++) {
            anc_put(FO, FN[k], anc_clean(FR["seq", k]))
            ROW[++rows] = "families\t" FN[k] "\t" length(FO["seq", k]) "\t0\t\t" substr(FO["seq", k], 1, 7) "\town frame"
        }
        for (k = 1; k <= nmem; k = j) {                                         # members grouped by family, config order
            anc_clear(PA)
            for (j = k; j <= nmem && MF[j] == MF[k]; j++) anc_put(PA, MN[j], SUB["seq", SUB["idx", MN[j]]])
            if (!(("idx", MF[k]) in FO)) _anc_fail("KeyError", "'" MF[k] "'")
            if (!anc_records(PA, OUTR, FO["seq", FO["idx", MF[k]]], "", S["mm"], S["mi"])) _anc_fail(ANC_ERRT, ANC_ERR)
            for (i = 1; i <= OUTR["n"]; i++) anc_put(SO, OUTR["name", i], OUTR["seq", i])
            for (i = 1; i <= OUTR["rows"]; i++) ROW[++rows] = "subfamilies\t" OUTR["row", i]
        }
    } else {                                                                    # :1078-1083
        if (!anc_records(FR, FO, frame, S["motif"], S["mm"], S["mi"])) _anc_fail(ANC_ERRT, ANC_ERR)
        for (i = 1; i <= FO["rows"]; i++) ROW[++rows] = "families\t" FO["row", i]
        if (SUB["n"]) {
            if (!anc_records(SUB, SO, frame, S["motif"], S["mm"], S["mi"])) _anc_fail(ANC_ERRT, ANC_ERR)
            for (i = 1; i <= SO["rows"]; i++) ROW[++rows] = "subfamilies\t" SO["row", i]
        }
    }
    out = S["outdir"]                                                           # :1084-1091
    anc_write_fasta(out "/families.anchored.fa", FO)
    if (SO["n"]) anc_write_fasta(out "/subfamilies.anchored.fa", SO)
    printf "set\tname\tlength\tshift\tidentity_pct\tstart\tstatus\n" > (out "/anchor_report.tsv")
    rotated = 0
    for (i = 1; i <= rows; i++) {
        printf "%s\n", ROW[i] > (out "/anchor_report.tsv")
        split(ROW[i], f, "\t"); if (f[4] != "0") rotated++                     # x[3] not in (0, '0') (:1093)
    }
    close(out "/anchor_report.tsv")
    anc_clear(PA); anc_put(PA, afam "_anchor", frame); anc_write_fasta(out "/anchor_frame.fa", PA)
    man = "{\"family\": " json_str(afam) ", \"motif\": " json_str(S["motif"]) ", \"frame_start\": " json_str(substr(frame, 1, 20)) \
          ", \"records\": " rows ", \"rotated\": " rotated "}"
    printf "MANIFEST\t%s\n", man
    printf "SAY\t%s\n", flat_esc(sprintf("  anchor %s (%s...): %d consensus records in frame, %d needed rotation", afam, substr(frame, 1, 12), rows, rotated))
    exit 0
}

# cli.py:91-125 after run_check: frame, then each file in order; --check verifies; otherwise gapped files stop the
# command (return 1, :115-117), records are rotated and written. AnchorError -> 'satpipe anchor: <e>' and 1 (:122-124).
# Other Python exceptions (duplicate record name ValueError, IndexError, unreadable file) are tracebacks there; here they
# print 'satpipe anchor: <message>' and also exit 1.
function anchor_cli(FN, FP, nfam, afam, S, IN_FP, IN_OUT, nin,   k, R, frame, idx, bad, P, i, A, fp, gapped, err) {
    err = "/dev/stderr"
    for (k = 1; k <= nfam; k++) if (FN[k] == afam) idx = k
    if (!idx) { printf "satpipe anchor: KeyError '%s'\n", afam > err; return 1 }
    if (!anc_read_fasta(FP[idx], R)) { printf "satpipe anchor: %s\n", ANC_ERR > err; return 1 }
    frame = anc_frame(anc_join_values(R), S["motif"], S["mm"])
    if (ANC_ERRT != "") { printf "satpipe anchor: %s\n", ANC_ERR > err; return 1 }
    bad = 0
    for (k = 1; k <= nin; k++) {
        fp = IN_FP[k]
        if (!anc_read_fasta(fp, R)) { printf "satpipe anchor: %s\n", ANC_ERR > err; return 1 }
        if (S["check"] == "1") {
            if (!anc_verify(R, P, frame, S["motif"], S["mm"], S["mi"])) { printf "satpipe anchor: %s\n", ANC_ERR > err; return 1 }
            for (i = 1; i <= P["n"]; i++) printf "%s: %s\n", fp, P[i] > err
            bad += P["n"]
            printf "%s: %d/%d records in the anchor frame\n", fp, R["n"] - P["n"], R["n"]
        } else {
            gapped = 0
            for (i = 1; i <= R["n"]; i++) if (index(R["seq", i], "-")) gapped = 1
            if (gapped) { printf "%s: gapped sequences (an alignment): rotate the unaligned sequences and re-align\n", fp > err; return 1 }
            if (!anc_records(R, A, frame, S["motif"], S["mm"], S["mi"])) { printf "satpipe anchor: %s\n", ANC_ERR > err; return 1 }
            anc_write_fasta(IN_OUT[k], A)
            printf "%s -> %s\n", fp, IN_OUT[k]
        }
    }
    return bad ? 1 : 0
}
