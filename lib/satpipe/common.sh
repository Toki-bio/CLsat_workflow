# common.sh - sourced first by bin/satpipe and every lib/satpipe/*.sh (PORT_SPEC section 1).
# Sets LC_ALL=C once (run.py:160), locates the checkout, provides die/say/log and the gawk check.

[ -n "${_SATPIPE_COMMON_SH:-}" ] && return 0
_SATPIPE_COMMON_SH=1

export LC_ALL=C

SATPIPE_LIB=$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]}")")
SATPIPE_ROOT=$(dirname -- "$(dirname -- "$SATPIPE_LIB")")
SCRIPTS=$SATPIPE_ROOT/satpipe/stages          # the stage scripts, called as today
PYJSON=$SATPIPE_LIB/pyjson.awk
GAWK=${SATPIPE_GAWK:-gawk}                     # risk 2: never plain `awk` (mawk on therioserver)
export SATPIPE_LIB SATPIPE_ROOT SCRIPTS PYJSON GAWK

# die MSG [CODE]   message to stderr, exit CODE (default 1)
die() { printf '%s\n' "$1" >&2; exit "${2:-1}"; }

# say MSG          Ctx.say: stdout, and the run log when SATPIPE_LOG is set
say() {
    printf '%s\n' "$1"
    if [ -n "${SATPIPE_LOG:-}" ]; then printf '%s\n' "$1" >> "$SATPIPE_LOG"; fi
}

# log MSG          stderr only (diagnostics that are not part of the Python driver's output)
log() { printf '%s\n' "$1" >&2; }

# require_gawk     GNU awk must be callable as $GAWK before anything else is read
require_gawk() {
    local v
    v=$("$GAWK" --version 2>/dev/null | head -n 1) || true
    case $v in
        "GNU Awk"*) return 0 ;;
        *) die "satpipe: GNU awk not found as '$GAWK' (set SATPIPE_GAWK or put gawk on PATH)" ;;
    esac
}

# pyjson_gawk PROGRAM [ARGS...]   run a gawk program with the pyjson.awk library loaded
pyjson_gawk() { local prog=$1; shift; "$GAWK" -f "$PYJSON" -e "$prog" "$@"; }
