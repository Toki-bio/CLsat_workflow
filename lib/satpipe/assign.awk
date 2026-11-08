# assign.awk - the logic of stage 09_assign (run.py:618-710): validate the curated cluster -> family table against the
# Step-8 clusters, write the centroid query and doubled-consensus library for ssearch36, then the cluster report.
# gawk -f pyjson.awk -f pytext.awk -f assign.awk   with ENVIRON["ASSIGN_JOB"] = job file (rows written by job_row):
#   set<TAB>mode<TAB>prep|report      prep: read everything, write q and lib (run.py:625-667);
#                                     report: read everything again + the ssearch36 output, write the report (668-710)
#   set<TAB>KEY<TAB>VALUE             cc table refs cent uc q lib rep ss (paths)
#   clu<TAB>NAME<TAB>PATH             one per entry of cross_clusters/ whose name starts with clust_ (NAME = the rest)
#   fam<TAB>NAME<TAB>CONSENSUS        the config's families in config order
# stdout: SAY / ERR / DIE protocol lines (pytext.awk); nothing else.

function job_read(f,   line, n, a, i, rs) {
    rs = RS; RS = "\n"
    while ((getline line < f) > 0) {
        n = split(line, a, "\t")
        for (i = 1; i <= n; i++) a[i] = flat_unesc(a[i])
        if (a[1] == "set") S[a[2]] = a[3]
        else if (a[1] == "clu") { NCLU++; CLU_NAME[NCLU] = a[2]; CLU_PATH[a[2]] = a[3]; IS_CLU[a[2]] = 1 }
        else if (a[1] == "fam") { NFAM++; FAM[NFAM] = a[2]; FAM_CONS[a[2]] = a[3]; IS_FAM[a[2]] = 1 }
    }
    close(f); RS = rs
}

# sorted(clusters, key=int): by int value; ties (never seen: "1" vs "01") by name, where Python keeps os.listdir order
function cmp_int(i1, v1, i2, v2) {
    if (CLU_INT[v1] < CLU_INT[v2]) return -1
    if (CLU_INT[v1] > CLU_INT[v2]) return 1
    return (v1 < v2) ? -1 : (v1 > v2)
}

BEGIN {
    job_read(ENVIRON["ASSIGN_JOB"])
    mode = S["mode"]
    # --- the curated table (run.py:629-637); dict `assigned` keeps first-insertion order, last value
    n = txt_read(S["table"], T)
    for (i = 1; i <= n; i++) {
        line = T[i]
        if (py_strip(line) == "" || substr(line, 1, 1) == "#") continue
        np = split(line, p, "\t")
        if (np < 2) continue
        clu = py_strip(p[1]); fam = py_strip(p[2])
        if (!(clu in IS_CLU)) ERRS[++NERR] = sprintf("cluster %s (family %s) does not exist in %s", clu, fam, S["cc"])
        if (!(fam in IS_FAM)) ERRS[++NERR] = sprintf("family %s is not in the config", fam)
        if (!(clu in ASSIGNED)) ASS_KEY[++NASS] = clu
        ASSIGNED[clu] = fam
    }
    for (k = 1; k <= NASS; k++) HAS_FAM[ASSIGNED[ASS_KEY[k]]] = 1
    for (i = 1; i <= NFAM; i++) if (!(FAM[i] in HAS_FAM)) ERRS[++NERR] = sprintf("family %s has no cluster in %s", FAM[i], S["table"])
    # --- reference names (642), centroid sequences (643-649), centroid of each cluster from the .uc S records (651-655)
    n = txt_read(S["refs"], L)
    for (i = 1; i <= n; i++) if (substr(L[i], 1, 1) == ">") REFS[py_word1(substr(L[i], 2))] = 1
    n = txt_read(S["cent"], L); name = ""; have = 0
    for (i = 1; i <= n; i++) {
        if (substr(L[i], 1, 1) == ">") { name = py_word1(substr(L[i], 2)); have = 1; CENT_SEQ[name] = "" }
        else { if (!have) py_die("None"); CENT_SEQ[name] = CENT_SEQ[name] py_strip(L[i]) }
    }
    n = txt_read(S["uc"], L)
    for (i = 1; i <= n; i++) {
        np = split(L[i], p, "\t")
        if (p[1] == "S") {
            if (np < 9) py_die("list index out of range")
            if (!(p[2] in CENT_OF)) CO_KEY[++NCO] = p[2]
            CENT_OF[p[2]] = p[9]
        }
    }
    # --- cluster sizes (656)
    for (i = 1; i <= NCLU; i++) {
        c = CLU_NAME[i]; m = txt_read(CLU_PATH[c], L); s = 0
        for (j = 1; j <= m; j++) if (substr(L[j], 1, 1) == ">") { py_word1(substr(L[j], 2)); s++ }
        SIZE[c] = s
    }
    if (mode == "prep") {
        # query: centroids in .uc order (658-663); library: doubled upper-case consensus per family (664-667)
        printf "" > S["q"]
        for (k = 1; k <= NCO; k++) { c = CO_KEY[k]; if (CENT_OF[c] in CENT_SEQ) printf ">%s\n%s\n", c, CENT_SEQ[CENT_OF[c]] > S["q"] }
        close(S["q"])
        printf "" > S["lib"]
        for (i = 1; i <= NFAM; i++) {
            m = txt_read(FAM_CONS[FAM[i]], L); s = ""
            for (j = 1; j <= m; j++) if (substr(L[j], 1, 1) != ">") s = s py_strip(L[j])
            s = toupper(s)
            printf ">%s\n%s%s\n", FAM[i], s, s > S["lib"]
        }
        close(S["lib"])
        exit 0
    }
    # --- best hit per centroid (670-681): r.stdout.splitlines(); first maximum kept (strict >)
    rs = RS; RS = "\r\n|[\n\r\v\f\034\035\036]"
    while ((getline line < S["ss"]) > 0) {
        if (substr(line, 1, 1) == "#") continue
        np = split(line, p, "\t")
        if (np < 12) continue
        if (!py_float(p[12])) continue
        if (!(p[1] in BEST_BITS) || PYF > BEST_BITS[p[1]]) { BEST_FAM[p[1]] = p[2]; BEST_BITS[p[1]] = PYF; BEST_PID[p[1]] = p[3] }
    }
    close(S["ss"]); RS = rs
    # --- report sorted by int(cluster) (682-689)
    for (i = 1; i <= NCLU; i++) {
        c = CLU_NAME[i]
        if (!py_int(c)) py_die("invalid literal for int() with base 10: " py_repr_str(c))
        CLU_INT[c] = PYI; ORD[i] = c
    }
    asort(ORD, SORTED, "cmp_int")
    rep = S["rep"]
    printf "cluster\tsize\tassigned_family\tcentroid\tis_reference\tbest_family\tbest_bits\tbest_pid\n" > rep
    for (i = 1; i <= NCLU; i++) {
        c = SORTED[i]; ce = (c in CENT_OF) ? CENT_OF[c] : ""
        if (c in BEST_BITS) { bf = BEST_FAM[c]; bb = BEST_BITS[c]; bp = BEST_PID[c] } else { bf = ""; bb = 0; bp = "" }
        printf "%s\t%d\t%s\t%s\t%s\t%s\t%.1f\t%s\n", c, SIZE[c], (c in ASSIGNED) ? ASSIGNED[c] : "", ce,
               (ce in REFS) ? "yes" : "no", bf, bb, bp > rep
    }
    close(rep)
    py_say(sprintf("  cluster report: %s (%d clusters, %d assigned)", rep, NCLU, NASS))
    # --- every assigned cluster's centroid must score best against its family (691-707)
    for (k = 1; k <= NASS; k++) {
        c = ASS_KEY[k]; fam = ASSIGNED[c]
        if ((c in IS_CLU) && (c in BEST_BITS) && BEST_FAM[c] != fam) {
            pid = py_float(BEST_PID[c]) ? PYF : 0
            msg = sprintf("cluster %s is assigned to %s but its centroid scores best against %s (%.1f bits, %s%% id)",
                          c, fam, BEST_FAM[c], BEST_BITS[c], BEST_PID[c])
            if (pid >= 95.0 && BEST_BITS[c] >= 100.0) ERRS[++NERR] = msg ": the table does not describe this run's clustering"
            else py_say("  WARNING: " msg " (diverged centroid; keeping the curated call)")
        }
    }
    if (NERR) {
        msg = "09_assign: the curated table does not match the clusters:\n  " ERRS[1]
        for (i = 2; i <= NERR; i++) msg = msg "\n  " ERRS[i]
        py_err(msg "\n  see " rep)
    }
}
