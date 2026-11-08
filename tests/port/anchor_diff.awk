# WP9 differential driver: evaluates tests/port/anchor_oracle.py's cases.tsv with lib/satpipe/anchor.awk and prints
# the same result lines as `anchor_oracle.py run`. Usage:
#   LC_ALL=C gawk -f lib/satpipe/pyjson.awk -f lib/satpipe/anchor.awk -f tests/port/anchor_diff.awk cases.tsv
function d_err() {
    if (ANC_ERRT == "AnchorError" || ANC_ERRT == "ValueError") return "ERR\t" ANC_ERRT "\t" flat_esc(ANC_ERR)
    return "ERR\t" ANC_ERRT "\t"
}
function d_recs(A,   k, s) {
    s = ""
    for (k = 1; k <= A["n"]; k++) s = s (k > 1 ? "," : "") A["name", k] "=" A["seq", A["idx", A["name", k]]]
    return flat_esc(s)
}
function d_parse(f, R,   n, x, i, p) {
    delete R; R["n"] = 0
    if (f == "") return
    n = split(f, x, ",")
    for (i = 1; i <= n; i++) { p = index(x[i], "="); anc_put(R, substr(x[i], 1, p - 1), substr(x[i], p + 1)) }
}
function d_probs(P,   i, s) { s = ""; for (i = 1; i <= P["n"]; i++) s = s (i > 1 ? "|" : "") P[i]; return flat_esc(s) }
function d_rows(A,   i, s) { s = ""; for (i = 1; i <= A["rows"]; i++) s = s (i > 1 ? "|" : "") A["row", i]; return flat_esc(s) }

BEGIN { FS = "\t" }
{
    sub(/\r$/, "")
    k = $1; o = k "\t" $2
    if (k == "A") { anc_local_align($3, $4); o = o "\t" ANC_AL[1] " " ANC_AL[2] " " ANC_AL[3] " " ANC_AL[4] " " ANC_AL[5] }
    else if (k == "F") { r = anc_frame($3, $4, $5); o = o "\t" (ANC_ERRT == "" ? "OK\t" r : d_err()) }
    else if (k == "R") {
        r = anc_rotate($3, $4, $5, $6, $7, $8, $9)
        o = o "\t" (ANC_ERRT == "" ? sprintf("OK\t%s\t%d\t%.17g\t%d\t%d", r, ANC_SHIFT, ANC_IDENT, ANC_COLS, ANC_SNAP) : d_err())
    }
    else if (k == "C") o = o "\t" (anc_check_frame($3, $4, $5) ? "1" : "0")
    else if (k == "RF") o = o "\t" (anc_read_fasta($3, R) ? "OK\t" d_recs(R) : d_err())
    else if (k == "P") {
        d_parse($7, R)
        frame = anc_frame($3, $4, $5)
        if (ANC_ERRT != "") o = o "\t" d_err()
        else {
            o = o "\tFRAME\t" frame
            if (!anc_verify(R, P, frame, $4, $5, $6)) o = o "\t" d_err()
            else {
                o = o "\tVERIFY\t" d_probs(P)
                if (!anc_records(R, OUT, frame, $4, $5, $6)) o = o "\t" d_err()
                else {
                    o = o "\tOUT\t" d_recs(OUT) "\tROWS\t" d_rows(OUT)
                    if (!anc_verify(OUT, P, frame, $4, $5, $6)) o = o "\t" d_err()
                    else o = o "\tREVERIFY\t" d_probs(P)
                }
            }
        }
    }
    print o
}
