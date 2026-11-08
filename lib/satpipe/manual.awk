# manual.awk - the manual-loci part of stage 15_catalogue (run.py:855-907).
# gawk -f pyjson.awk -f pytext.awk -f manual.awk   with ENVIRON["MANUAL_JOB"] = job file (rows written by job_row):
#   set<TAB>mode<TAB>prep|fam|final
#   set<TAB>manual<TAB>PATH          the manual_loci BED
#   gen<TAB>CODE                     the config's genome codes
#   fam<TAB>NAME                     the config's families in config order
#   prep:  set mbed PATH             write catalogue/manual.bed; ERR when a locus has no `<genome>|` name (run.py:859-870)
#   fam:   set fam NAME, set bits PRE.bits.tsv, set out catalogue/MASTER_<fam>.bits.tsv, set state FILE
#          append the scan rows (header skipped) to OUT and the per-locus maximum bits of this family to STATE as
#          `i<TAB>fam<TAB>%.17g` (run.py:874-881); the keys of maxbits[i] are inserted in family order, as in Python
#   final: set state FILE, set nmaster N, set mb lb sb bb fa PATHS (the five catalogue files, appended)
#          the rows of the five files (run.py:882-906); stdout `JSON<TAB>[...]` = the manifest's manual_loci list
# stdout: SAY / ERR / DIE protocol lines (pytext.awk) and the JSON line.

function job_read(f,   line, n, a, i, rs) {
    rs = RS; RS = "\n"
    while ((getline line < f) > 0) {
        n = split(line, a, "\t")
        for (i = 1; i <= n; i++) a[i] = flat_unesc(a[i])
        if (a[1] == "set") S[a[2]] = a[3]
        else if (a[1] == "gen") IS_GEN[a[2]] = 1
        else if (a[1] == "fam") { NFAM++; FAM[NFAM] = a[2]; IS_FAM[a[2]] = 1 }
    }
    close(f); RS = rs
}

# rows = [p for p in read_tsv(manual, 3) if not p[0].startswith(('track', 'browser'))]; R[i, j] (1-based), RN[i] = len(p);
# MSP[i] = p[3].split('|', 1)[0] when len(p) >= 4 and '|' in p[3], else ''
function read_rows(path,   n, i, line, np, p, j) {
    n = txt_read(path, TL)
    NROW = 0
    for (i = 1; i <= n; i++) {
        line = TL[i]
        if (py_strip(line) == "" || substr(line, 1, 1) == "#") continue
        np = split(line, p, "\t")
        if (np < 3) continue
        if (substr(p[1], 1, 5) == "track" || substr(p[1], 1, 7) == "browser") continue
        NROW++; RN[NROW] = np
        for (j = 1; j <= np; j++) R[NROW, j] = p[j]
        MSP[NROW] = (np >= 4 && index(p[4], "|")) ? substr(p[4], 1, index(p[4], "|") - 1) : ""
    }
}

# repr(p[:4]) of row i
function row4_repr(i,   k, A) { for (k = 1; k <= 4 && k <= RN[i]; k++) A[k] = R[i, k]; return py_repr_strlist(A, k - 1) }

function need_int(s) { if (!py_int(s)) py_die("invalid literal for int() with base 10: " py_repr_str(s)); return PYI }

BEGIN {
    job_read(ENVIRON["MANUAL_JOB"])
    mode = S["mode"]
    read_rows(S["manual"])
    if (mode == "prep") {
        printf "" > S["mbed"]
        for (i = 1; i <= NROW; i++) {
            if (RN[i] >= 4) printf "%s\t%s\t%s\t%s\n", R[i, 1], R[i, 2], R[i, 3], R[i, 4] > S["mbed"]
            else printf "%s\t%s\t%s\tMANUAL\n", R[i, 1], R[i, 2], R[i, 3] > S["mbed"]
        }
        close(S["mbed"])
        nb = 0; bad = ""
        for (i = 1; i <= NROW; i++) if (!(MSP[i] in IS_GEN)) { if (++nb <= 3) bad = bad (nb > 1 ? ", " : "") row4_repr(i) }
        if (nb) py_err("15_catalogue: manual loci need a name `<genome>|...` with a genome of the config: [" bad "]")
        exit 0
    }
    if (mode == "fam") {
        fam = S["fam"]; f = S["bits"]; out = S["out"]
        rs = RS; RS = "\r\n|\r|\n"; ERRNO = ""
        r = (getline line < f)                                     # fh.readline(): the header
        if (r < 0) py_die(_pt_oserror(f, ERRNO))
        printf "" >> out
        if (r > 0) while ((getline line < f) > 0) {
            printf "%s%s", line, (RT == "" ? "" : "\n") >> out
            np = split(line, p, "\t")
            for (i = 1; i <= NROW; i++) {
                if (p[1] != R[i, 1]) continue
                if (np < 4) py_die("list index out of range")
                g = index(p[4], "|") ? substr(p[4], 1, index(p[4], "|") - 1) : p[4]
                if (g != MSP[i]) continue
                a = need_int(R[i, 2])
                if (np < 2) py_die("list index out of range")
                b = need_int(p[2])
                if (!(a < b)) continue
                if (!(b <= need_int(R[i, 3]))) continue
                prev = (i in MAXB) ? MAXB[i] : 0
                if (np < 3) py_die("list index out of range")
                if (!py_float(p[3])) py_die("could not convert string to float: " py_repr_str(p[3]))
                MAXB[i] = (PYF > prev) ? PYF : prev                 # max(prev, x): the first argument on ties
            }
        }
        close(f); RS = rs
        for (i = 1; i <= NROW; i++) if (i in MAXB) printf "%d\t%s\t%.17g\n", i, fam, MAXB[i] >> S["state"]
        close(S["state"]); close(out)
        exit 0
    }
    # final
    st = S["state"]; rs = RS; RS = "\n"
    while ((getline line < st) > 0) {
        split(line, a, "\t"); i = a[1] + 0
        BK[i, ++BN[i]] = a[2]; BV[i, a[2]] = a[3] + 0
    }
    close(st); RS = rs
    nm = S["nmaster"] + 0
    json = "["
    for (i = 1; i <= NROW; i++) {
        lid = "L" (nm + i)
        name = (RN[i] >= 4) ? R[i, 4] : "MANUAL"
        if (RN[i] >= 7 && R[i, 7] != "") {
            if (!(R[i, 7] in IS_FAM))
                py_err(sprintf("15_catalogue: manual locus %s has family %s, not a configured family", row4_repr(i), py_repr_str(R[i, 7])))
            fam = R[i, 7]
        } else {
            any = 0
            for (k = 1; k <= BN[i]; k++) if (BV[i, BK[i, k]] > 0) any = 1
            if (!any) py_err(sprintf("15_catalogue: manual locus %s scores 0 against every family; give its family in col 7", row4_repr(i)))
            # max(sorted(keys), key=get): the highest bits; ties -> the alphabetically first family
            fam = ""
            for (k = 1; k <= BN[i]; k++) {
                f = BK[i, k]
                if (fam == "" || BV[i, f] > BV[i, fam] || (BV[i, f] == BV[i, fam] && f < fam)) fam = f
            }
        }
        fb = ((i, fam) in BV) ? BV[i, fam] : 0
        # runner-up: sorted(((b, f) ...), reverse=True)[0]: the highest bits; ties -> the alphabetically last family
        bo = ""; bbits = 0; have = 0
        for (k = 1; k <= BN[i]; k++) {
            f = BK[i, k]
            if (f == fam) continue
            if (!have || BV[i, f] > bbits || (BV[i, f] == bbits && f > bo)) { bo = f; bbits = BV[i, f]; have = 1 }
        }
        printf "%s\t%s\t%s\t%s\tMANUAL\n", R[i, 1], R[i, 2], R[i, 3], name >> S["mb"]
        printf "%s\t%s\t%s\t%s\t%s\n", R[i, 1], R[i, 2], R[i, 3], lid, MSP[i] >> S["lb"]
        printf "%s\t%s\t%s\t%s\t%s\n", R[i, 1], R[i, 2], R[i, 3], lid, MSP[i] >> S["sb"]
        line = lid
        for (k = 1; k <= NFAM; k++) line = line "\t" sprintf("%g", ((i, FAM[k]) in BV) ? BV[i, FAM[k]] : 0)
        printf "%s\n", line >> S["bb"]
        printf "%s\t%s\t%g\t%s\t%g\t%g\tMANUAL\tMANUAL\n", lid, fam, fb, bo, bbits, fb - bbits >> S["fa"]
        # manifest entry {'locus', 'bed': p[:4], 'family', 'bits'} (floats as Python repr)
        e = "{\"locus\": " json_str(lid) ", \"bed\": ["
        for (k = 1; k <= 4 && k <= RN[i]; k++) e = e (k > 1 ? ", " : "") json_str(R[i, k])
        e = e "], \"family\": " json_str(fam) ", \"bits\": {"
        for (k = 1; k <= BN[i]; k++) e = e (k > 1 ? ", " : "") json_str(BK[i, k]) ": " pyrepr_float(BV[i, BK[i, k]])
        json = json (i > 1 ? ", " : "") e "}}"
    }
    printf "JSON\t%s]\n", json
}
