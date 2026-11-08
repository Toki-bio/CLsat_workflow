# crosscheck.awk - the checks that combine several keys (config.py:42-82, _cross_checks), run on the TYPED flat
# that validate.awk wrote (config.py:101-102: only when the schema produced no error).
#
# Usage:  LC_ALL=C KNOWN_TOOLS="bash gawk ..." gawk -f pyjson.awk -f crosscheck.awk FLAT     (FLAT may be - = stdin)
#   KNOWN_TOOLS = the names of satpipe/tools.py TOOLS (config.sh supplies them); bin_dir is added here.
# Output: valid -> the TYPED flat lines unchanged, then one `!warning<TAB>msg` line per warning; exit 0.
#         errors -> `!error<TAB>msg` lines only, in config.py order; exit 1.  Messages flat-escaped (\\ \t \n \r).
#         `!error` input lines (from yaml.awk or validate.awk) are passed through; exit 1.

BEGIN { FS = "\t"; N = 0; NERR = 0; NWARN = 0 }

{
    sub(/\r$/, "")
    if (substr($0, 1, 1) == "!") { print; PASSERR = 1; next }
    LINE[++N] = $0
    if (!flat_load("C", $0)) { print "crosscheck.awk: " JERR > "/dev/stderr"; exit 2 }
}

END {
    if (PASSERR) exit 1
    if (!("C" in JT)) jt_init("C", "map")
    cross()
    if (NERR) { for (i = 1; i <= NERR; i++) print "!error\t" flat_esc(ERR[i]); exit 1 }
    for (i = 1; i <= N; i++) print LINE[i]
    for (i = 1; i <= NWARN; i++) print "!warning\t" flat_esc(WARN[i])
    exit 0
}

function err(m) { ERR[++NERR] = m }

# node id of a slash path, "" when absent
function node(p,   a, n, i, id) {
    id = "C"; n = split(p, a, "/")
    for (i = 1; i <= n; i++) { id = id SUBSEP a[i]; if (!(id in JT)) return "" }
    return id
}
# id of a map node, "" when absent / not a map / empty (Python `x or {}` then iteration)
function mapnode(p,   id) { id = node(p); return (id != "" && JT[id] == "map") ? id : "" }

function cross(   g, f, codes, fams, trf, i, k, n, sub_id, mem, m, j, lst, seen, dup, cons, af, known, kn, a, sorted, s, per, p0, p1, sn, tl, rest) {
    if ((g = mapnode("genomes")) != "") for (i = 1; i <= JN[g]; i++) codes[JK[g, i]] = 1
    if ((f = mapnode("families")) != "") for (i = 1; i <= JN[f]; i++) fams[JK[f, i]] = 1
    trf = mapnode("inputs/trf_dat")
    n = 0
    if (trf != "") {
        for (i = 1; i <= JN[trf]; i++) { k = JK[trf, i]; if (!(k in codes)) err("inputs.trf_dat." k ": not a code in `genomes`") }
        n = JN[trf]
    }
    if (n > 0) {                                     # trf and set(trf) <= codes and set(trf) != codes
        sub_ok = 1
        for (i = 1; i <= JN[trf]; i++) { k = JK[trf, i]; intrf[k] = 1; if (!(k in codes)) sub_ok = 0 }
        rest = ""
        if (sub_ok) {
            PROCINFO["sorted_in"] = "@ind_str_asc"
            for (k in codes) if (!(k in intrf)) rest = rest (rest == "" ? "" : ", ") k
            delete PROCINFO["sorted_in"]
            if (rest != "")
                WARN[++NWARN] = "inputs.trf_dat covers only part of `genomes`; stage 02 (TRF) will run for the rest: " rest
        }
    }
    if ((s = mapnode("inputs/seeds")) != "")
        for (i = 1; i <= JN[s]; i++) { k = JK[s, i]; if (!(k in fams)) err("inputs.seeds." k ": not a family in `families`") }
    mem = mapnode("subfamily/members")
    if (mem != "") {
        for (i = 1; i <= JN[mem]; i++) {
            k = JK[mem, i]; lst = mem SUBSEP k
            if (!(k in fams)) err("subfamily.members." k ": not a family in `families`")
            if (JT[lst] == "list" && JN[lst] > 0) {
                delete seen; dup = 0
                for (j = 1; j <= JN[lst]; j++) { m = JV[lst SUBSEP JK[lst, j]]; if (m in seen) dup = 1; seen[m] = 1 }
                if (dup) err("subfamily.members." k ": duplicate names")
                for (j = 1; j <= JN[lst]; j++) {
                    m = JV[lst SUBSEP JK[lst, j]]
                    if ((m in fams) && m != k) err("subfamily.members." k ": '" m "' is also a family name")
                }
            }
        }
        cons = node("subfamily/consensuses")
        if (JN[mem] > 0 && (cons == "" || JT[cons] == "null" || JV[cons] == ""))
            err("subfamily.members given but subfamily.consensuses is missing")
    }
    if ((sn = mapnode("viewer/species_names")) != "")
        for (i = 1; i <= JN[sn]; i++) { k = JK[sn, i]; if (!(k in codes)) err("viewer.species_names." k ": not a code in `genomes`") }
    af = node("anchor/family")
    if (af != "" && JT[af] == "str" && JV[af] != "" && !(JV[af] in fams))
        err("anchor.family: '" JV[af] "' is not a family in `families`")
    n = split(ENVIRON["KNOWN_TOOLS"], a, " ")
    for (i = 1; i <= n; i++) if (a[i] != "") known[a[i]] = 1
    known["bin_dir"] = 1
    kn = ""
    PROCINFO["sorted_in"] = "@ind_str_asc"
    for (k in known) kn = kn (kn == "" ? "" : ", ") k
    delete PROCINFO["sorted_in"]
    if ((tl = mapnode("tools")) != "")
        for (i = 1; i <= JN[tl]; i++) { k = JK[tl, i]; if (!(k in known)) err("tools." k ": not a known tool (known: " kn ")") }
    per = node("discovery/period")
    if (per != "" && JT[per] == "list" && JN[per] == 2) {
        p0 = per SUBSEP JK[per, 1]; p1 = per SUBSEP JK[per, 2]
        if (JT[p0] == "int" && JT[p1] == "int" && JV[p0] + 0 > JV[p1] + 0)
            err("discovery.period: min " JV[p0] " > max " JV[p1])
    }
}

# one TYPED-flat line into the tree under root; unlike pyjson flat_line, path components are flat-unescaped too
function flat_load(root, line,   f, n, a, i, id) {
    n = split(line, f, "\t")
    if (n < 3 || f[2] !~ /^(map|list|int|float|str|path|bool|null)$/) { JERR = "malformed flat line: " line; return 0 }
    if (!(root in JT)) jt_init(root, "map")
    n = split(f[1], a, "/"); id = root
    for (i = 1; i < n; i++) { id = id SUBSEP flat_unesc(a[i]); if (!(id in JT)) { JERR = "flat line before its parent: " line; return 0 } }
    jt_set(id, flat_unesc(a[n]), f[2], (f[2] == "map" || f[2] == "list") ? "" : flat_unesc(f[3]))
    return 1
}
