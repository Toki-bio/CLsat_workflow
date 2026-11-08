# tools.sh - external tool registry (tools.tsv = tools.py:17-37) and version probe (tools.py:60-108).
# Sourced after common.sh (PORT_SPEC sections 1, 4 "tools probe").
#
# API
#   tools_probe_all [FLAT]   one line per tools.tsv row, in table order (probe_all, tools.py:105-108):
#                              name<TAB>required(True|False)<TAB>stages<TAB>note<TAB>found(1|0)<TAB>path<TAB>version<TAB>error
#                            fields flat-escaped (\\ \t \n \r). FLAT = the checked config as TYPED flat: its `tools` block gives
#                            the per-tool overrides and `bin_dir` (cfg.data.get("tools") or {}); without FLAT no overrides.
#   tools_locate NAME [OVERRIDE] [BIN_DIR]   locate() (tools.py:65-75); prints the path, status 1 when not found
#
# Environment the probe looks at (set by bin/satpipe before common.sh changes anything):
#   SATPIPE_USER_PATH / SATPIPE_USER_PATH_SET    the caller's PATH (shutil.which searches it; the probed tool runs with it)
#   SATPIPE_USER_LC_ALL / SATPIPE_USER_LC_ALL_SET the caller's LC_ALL (version output is read in the caller's locale, as Python does)
# Without them (library used on its own) the current PATH / LC_ALL are taken.

[ -n "${_SATPIPE_TOOLS_SH:-}" ] && return 0
_SATPIPE_TOOLS_SH=1

SATPIPE_TOOLS_TSV=${SATPIPE_TOOLS_TSV:-$SATPIPE_LIB/tools.tsv}
SATPIPE_PROBE_TIMEOUT=${SATPIPE_PROBE_TIMEOUT:-30}          # tools.py:78 timeout=30

# _tools_exe P      os.path.isfile(p) and os.access(p, X_OK)
_tools_exe() { [ -f "$1" ] && [ -x "$1" ]; }

# _tools_join A B   os.path.join for an absolute or relative A and a plain name B
_tools_join() { case $1 in ''|*/) printf '%s%s' "$1" "$2" ;; *) printf '%s/%s' "$1" "$2" ;; esac; }

# _tools_which NAME   shutil.which (POSIX): PATH unset -> CS_PATH /bin:/usr/bin; PATH "" -> None; empty entry -> NAME itself;
#                     a hit = exists, X_OK, not a directory
_tools_which() {
    local name=$1 p d c
    if [ -n "${SATPIPE_USER_PATH_SET+x}" ]; then
        if [ "$SATPIPE_USER_PATH_SET" = 1 ]; then p=$SATPIPE_USER_PATH; else p=/bin:/usr/bin; fi
    else
        if [ -n "${PATH+x}" ]; then p=$PATH; else p=/bin:/usr/bin; fi
    fi
    [ -n "$p" ] || return 1
    local -A seen=()
    local -a dirs=()
    local IFS=:
    read -r -d '' -a dirs < <(printf '%s:' "$p") || true     # p + ':': read drops one trailing ':', str.split keeps it
    for d in "${dirs[@]}"; do
        [ -n "${seen[x$d]+x}" ] && continue
        seen[x$d]=1
        c=$(_tools_join "$d" "$name")
        if [ -e "$c" ] && [ -x "$c" ] && [ ! -d "$c" ]; then printf '%s' "$c"; return 0; fi
    done
    return 1
}

tools_locate() {
    local name=$1 ov=${2:-} bin=${3:-} c
    if [ -n "$ov" ]; then
        if _tools_exe "$ov"; then printf '%s' "${ov//\\//}"; return 0; fi
        return 1
    fi
    if [ -n "$bin" ]; then
        c=$(_tools_join "$bin" "$name")
        if _tools_exe "$c"; then printf '%s' "${c//\\//}"; return 0; fi
    fi
    c=$(_tools_which "$name") || { [ "$name" = python3 ] && c=$(_tools_which python); } || return 1
    printf '%s' "${c//\\//}"
}

_tools_esc() { local t=$1; t=${t//\\/\\\\}; t=${t//$'\t'/\\t}; t=${t//$'\n'/\\n}; t=${t//$'\r'/\\r}; printf '%s' "$t"; }

# _tools_run PATH ARGV...  run ARGV (ARGV[0] = PATH) as the tool as subprocess.run(..., stdin=DEVNULL, timeout=30) does, in the caller's PATH and
#                          LC_ALL; stdout -> $_TP_OUT, stderr -> $_TP_ERR; sets _TP_OSERR to the OSError/TimeoutExpired text
_tools_run() {
    local path=$1 rc=0 t0 magic=""
    shift
    _TP_OSERR=""
    : > "$_TP_OUT"; : > "$_TP_ERR"
    # execve() of a file that is neither '#!' nor ELF fails with ENOEXEC in Python; timeout(1)'s execvp would hand it to sh
    IFS= read -r -n 4 -d '' magic < "$path" 2>/dev/null || true
    case $magic in
        '#!'*|$'\x7fELF') ;;
        *) _TP_OSERR="[Errno 8] Exec format error: $(_tools_pyrepr "$path")"; return 0 ;;
    esac
    t0=$SECONDS
    (
        if [ -n "${SATPIPE_USER_LC_ALL_SET+x}" ]; then
            if [ "$SATPIPE_USER_LC_ALL_SET" = 1 ]; then export LC_ALL=$SATPIPE_USER_LC_ALL; else unset LC_ALL; fi
        fi
        if [ -n "${SATPIPE_USER_PATH_SET+x}" ]; then
            if [ "$SATPIPE_USER_PATH_SET" = 1 ]; then export PATH=$SATPIPE_USER_PATH; else unset PATH; fi
        fi
        exec "$_TP_TIMEOUT" "$SATPIPE_PROBE_TIMEOUT" "$@"
    ) < /dev/null > "$_TP_OUT" 2> "$_TP_ERR" || rc=$?
    case $rc in
        124) if [ $((SECONDS - t0)) -ge $((SATPIPE_PROBE_TIMEOUT - 1)) ]; then
                 _TP_OSERR="Command $(_tools_pyrepr "$(_tools_pylist "$@")") timed out after $SATPIPE_PROBE_TIMEOUT seconds"
             fi ;;
        126|127) # timeout(1) could not exec it (its own message is the last stderr line): EACCES / ENOENT in Python
             if "$GAWK" 'END { exit !($0 ~ /^([^ ]*\/)?timeout: /) }' "$_TP_ERR"; then     # GNU prints argv[0]: /usr/bin/timeout:
                 if [ "$rc" = 126 ]; then _TP_OSERR="[Errno 13] Permission denied: $(_tools_pyrepr "$path")"
                 else _TP_OSERR="[Errno 2] No such file or directory: $(_tools_pyrepr "$path")"; fi
             fi ;;
    esac
    return 0
}

_tools_pyrepr() { S=$1 "$GAWK" -f "$PYJSON" -e 'BEGIN { printf "%s", py_repr_str(ENVIRON["S"]) }' < /dev/null; }
# _tools_pylist ARGS...   str() of a Python list of str (used inside TimeoutExpired's message, which quotes it again)
_tools_pylist() {
    local a out="[" sep=""
    for a in "$@"; do out+=$sep$(_tools_pyrepr "$a"); sep=", "; done
    printf '%s]' "$out"
}

# _tools_flat_overrides FLAT   rows `name<TAB>path` of the tools block with a non-empty str/path value
_tools_flat_overrides() {
    "$GAWK" -f "$PYJSON" -e 'BEGIN { FS = "\t" }
        { sub(/\r$/, "") }
        substr($1, 1, 6) == "tools/" && index(substr($1, 7), "/") == 0 && ($2 == "path" || $2 == "str") {
            v = flat_unesc($3); if (v != "") printf "%s\t%s\n", substr($1, 7), flat_esc(v) }' "$1"
}

# tools_tsv_us FILE   the TSV with TAB -> \037 (IFS=TAB would merge empty fields: TAB is IFS whitespace), CR dropped
tools_tsv_us() { "$GAWK" '{ sub(/\r$/, ""); gsub(/\t/, "\037"); print }' "$1"; }

tools_probe_all() {
    local flat=${1:-} name args re ere req stages note minv ov bin="" path ver err line
    local -A OV=()
    if [ -n "$flat" ]; then
        while IFS=$'\037' read -r name path; do
            OV[$name]=$(S=$path "$GAWK" -f "$PYJSON" -e 'BEGIN { printf "%s", flat_unesc(ENVIRON["S"]) }' < /dev/null)
        done < <(_tools_flat_overrides "$flat" | tools_tsv_us /dev/stdin)
        bin=${OV[bin_dir]:-}
    fi
    _TP_TIMEOUT=$(command -v timeout) || die "satpipe: timeout(1) not found (coreutils)"
    _TP_OUT=$(mktemp "${TMPDIR:-/tmp}/satpipe_probe_out.XXXXXX")
    _TP_ERR=$(mktemp "${TMPDIR:-/tmp}/satpipe_probe_err.XXXXXX")
    while IFS=$'\037' read -r name args re ere req stages note minv; do
        case $name in ''|'#'*) continue ;; esac
        ov=${OV[$name]:-}
        if ! path=$(tools_locate "$name" "$ov" "$bin"); then
            printf '%s\t%s\t%s\t%s\t0\t\t\t\n' "$name" "$req" "$(_tools_esc "$stages")" "$(_tools_esc "$note")"
            continue
        fi
        local -a av=()
        read -r -a av <<< "$args"
        _tools_run "$path" "$path" "${av[@]}"
        if [ -n "$_TP_OSERR" ]; then
            ver=""; err="could not run: $_TP_OSERR"
        else
            line=$(ERE=$ere MINV=$minv NAME=$name "$GAWK" -f "$PYJSON" -f "$SATPIPE_LIB/check.awk" -e 'BEGIN { ck_version_probe() }' "$_TP_OUT" "$_TP_ERR" < /dev/null)
            ver=${line%%$'\t'*}; err=${line#*$'\t'}
            ver=$(S=$ver "$GAWK" -f "$PYJSON" -e 'BEGIN { printf "%s", flat_unesc(ENVIRON["S"]) }' < /dev/null)
            err=$(S=$err "$GAWK" -f "$PYJSON" -e 'BEGIN { printf "%s", flat_unesc(ENVIRON["S"]) }' < /dev/null)
        fi
        printf '%s\t%s\t%s\t%s\t1\t%s\t%s\t%s\n' "$name" "$req" "$(_tools_esc "$stages")" "$(_tools_esc "$note")" \
            "$(_tools_esc "$path")" "$(_tools_esc "$ver")" "$(_tools_esc "$err")"
    done < <(tools_tsv_us "$SATPIPE_TOOLS_TSV")
    rm -f -- "$_TP_OUT" "$_TP_ERR"
}
