# pytext.awk - Python text-file and literal semantics for the stage logic in gawk (WP7: assign.awk, manual.awk).
# Load after pyjson.awk:  gawk -f pyjson.awk -f pytext.awk -f <program>.awk   (LC_ALL=C)
#
#   txt_read(path, L)      list(open(path)) without line ends: universal newlines (\r\n, \r, \n), last line kept when it has
#                          no newline; L[1..n], returns n. An open failure calls py_die with the OSError text.
#   py_strip(s)            str.strip() (ASCII whitespace incl. \x1c-\x1f, as str.isspace)
#   py_word1(s)            s.split()[0]; py_die("list index out of range") when s has no word
#   py_int(s)              int(s): 1 and PYI set, or 0 (caller raises); accepts surrounding whitespace, sign, '_' separators
#   py_float(s)            float(s) for finite decimal literals: 1 and PYF set, or 0. inf/nan are NOT accepted (never written
#                          by ssearch36 or the scan script; Python would accept them: documented divergence)
#   py_repr_strlist(A, n)  repr() of a list of n str: ['a', 'b']
#   py_die(msg)            print DIE<TAB>msg (flat_esc'd) and exit (the bash side raises it as `satpipe run: msg`)
#   py_err(msg)            print ERR<TAB>msg (a StageError text) and exit
#   py_say(msg)            print SAY<TAB>msg (Ctx.say)

function py_die(msg) { printf "DIE\t%s\n", flat_esc(msg); _PT_DONE = 1; exit 0 }
function py_err(msg) { printf "ERR\t%s\n", flat_esc(msg); _PT_DONE = 1; exit 0 }
function py_say(msg) { printf "SAY\t%s\n", flat_esc(msg) }

function _pt_oserror(path, e,   n) {
    n = 5
    if (e ~ /No such file/) n = 2
    else if (e ~ /Is a directory/) n = 21
    else if (e ~ /Permission denied/) n = 13
    else if (e ~ /Not a directory/) n = 20
    return sprintf("[Errno %d] %s: %s", n, e, py_repr_str(path))
}

function txt_read(path, L,   n, line, rs, r) {
    delete L
    rs = RS; RS = "\r\n|\r|\n"; n = 0
    ERRNO = ""
    while ((r = (getline line < path)) > 0) L[++n] = line
    RS = rs
    if (r < 0) py_die(_pt_oserror(path, ERRNO))
    close(path)
    return n
}

function py_strip(s) { sub(/^[ \t\n\r\v\f\034\035\036\037]+/, "", s); sub(/[ \t\n\r\v\f\034\035\036\037]+$/, "", s); return s }

function py_word1(s) {
    if (match(s, /[^ \t\n\r\v\f\034\035\036\037]+/)) return substr(s, RSTART, RLENGTH)
    py_die("list index out of range")
}

function py_int(s,   t) {
    t = py_strip(s)
    if (t !~ /^[-+]?[0-9](_?[0-9])*$/) return 0
    gsub(/_/, "", t); PYI = t + 0
    return 1
}

function py_float(s,   t) {
    t = py_strip(s)
    if (t !~ /^[-+]?([0-9](_?[0-9])*(\.([0-9](_?[0-9])*)?)?|\.[0-9](_?[0-9])*)([eE][-+]?[0-9](_?[0-9])*)?$/) return 0
    gsub(/_/, "", t); PYF = t + 0
    return 1
}

function py_repr_strlist(A, n,   i, out) {
    out = "["
    for (i = 1; i <= n; i++) out = out (i > 1 ? ", " : "") py_repr_str(A[i])
    return out "]"
}
