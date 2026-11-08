# pyio.sh - Python file-operation equivalents used by the stage bodies of stages 09-15 (WP7) and the viewer.
# Sourced by st_09_assign.sh .. st_15_catalogue.sh and expand_cmd.sh, after common.sh (needs $GAWK, $PYJSON, stage_die).
#
#   py_repr_q STR               Python repr() of a str ('x' / "it's"), as embedded in exception texts
#   py_oserror PATH [w]         "[Errno N] <strerror>: 'PATH'" of open(PATH) for reading (or writing with w)
#   py_copyfile SRC DST         shutil.copyfile: content only, written through an existing DST (run.py:870-874, :751)
#   py_nlines FILE              sum(1 for _ in open(FILE)): text-mode line count (universal newlines: \r\n, \r, \n)
#   job_row FIELD...            one TAB-joined row of a gawk job file, every field flat_esc'd (\ TAB LF CR)
# A Python exception that is not a StageError (FileNotFoundError, ValueError, IndexError) is reported like WP5 does for
# tracebacks: `satpipe run: <exception text>` (stage_die), exit 1.

[ -n "${_SATPIPE_PYIO_SH:-}" ] && return 0
_SATPIPE_PYIO_SH=1

py_repr_q() { S=$1 "$GAWK" -f "$PYJSON" -e 'BEGIN { printf "%s", py_repr_str(ENVIRON["S"]) }' < /dev/null; }

py_oserror() {
    local p=$1 n t d
    if [ "${2:-}" = w ]; then
        d=$(dirname -- "$p")
        if [ ! -d "$d" ]; then n=2; t='No such file or directory'
        elif [ -d "$p" ]; then n=21; t='Is a directory'
        else n=13; t='Permission denied'; fi
    elif [ ! -e "$p" ]; then n=2; t='No such file or directory'
    elif [ -d "$p" ]; then n=21; t='Is a directory'
    else n=13; t='Permission denied'; fi
    printf '[Errno %d] %s: %s' "$n" "$t" "$(py_repr_q "$p")"
}

_py_readable() { [ -e "$1" ] && [ ! -d "$1" ] && [ -r "$1" ]; }

py_copyfile() {
    _py_readable "$1" || stage_die "$(py_oserror "$1")"
    { cat -- "$1" > "$2"; } 2> /dev/null || stage_die "$(py_oserror "$2" w)"
}

py_nlines() {
    _py_readable "$1" || stage_die "$(py_oserror "$1")"
    "$GAWK" 'BEGIN { RS = "\r\n|\r|\n" } END { print NR }' "$1"
}

job_row() {
    local out="" f sep=""
    for f in "$@"; do
        f=${f//\\/\\\\}; f=${f//$'\t'/\\t}; f=${f//$'\n'/\\n}; f=${f//$'\r'/\\r}
        out+=$sep$f; sep=$'\t'
    done
    printf '%s\n' "$out"
}

# py_gawk_result FILE   act on the protocol lines a stage gawk program printed: SAY<TAB>msg -> say; DIE<TAB>msg / ERR<TAB>msg
# -> stage_die (DIE = a Python exception text, ERR = a StageError text; both end the run the same way). Other lines are
# left to the caller. Messages are flat_esc'd.
py_gawk_result() {
    local k m
    while IFS=$'\t' read -r k m || [ -n "$k" ]; do
        case $k in
            SAY) say "$(S=$m "$GAWK" -f "$PYJSON" -e 'BEGIN { printf "%s", flat_unesc(ENVIRON["S"]) }' < /dev/null)" ;;
            DIE|ERR) stage_die "$(S=$m "$GAWK" -f "$PYJSON" -e 'BEGIN { printf "%s", flat_unesc(ENVIRON["S"]) }' < /dev/null)" ;;
        esac
    done < "$1"
    return 0
}
