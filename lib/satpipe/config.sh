# config.sh - load_config (config.py:85-108) as yaml.awk | validate.awk | crosscheck.awk, plus the cfg helpers and
# `satpipe config` (cli.py cmd_config). Sourced after common.sh and flat.sh (PORT_SPEC sections 1, 2.1, 3).
# Can also be run directly:  bash lib/satpipe/config.sh [--flat|--schema-doc] [FILE]   (= `satpipe config ...`)
#
# API
#   config_stream FILE           the whole pipeline's output on stdout: TYPED flat + `!warning<TAB>msg` lines (exit 0)
#                                or `!error<TAB>msg` lines only (exit 1); msg flat-escaped. Sets nothing.
#   config_load FILE [FLATOUT]   runs config_stream; sets CFG_PATH (abs config path, as Config.path), CFG_BASE_DIR,
#                                CFG_FLAT (= FLATOUT, default a mktemp file), CFG_ERRORS[], CFG_WARNINGS[] (unescaped).
#                                Returns 1 when CFG_ERRORS is not empty (ConfigError), else 0.
#   cfg PATH                     Python str() of the node at PATH (slash path) in CFG_FLAT (leaf text, repr for containers)
#   cfg_type PATH                type column of PATH
#   cfg_keys PATH                child keys of PATH in config order
#   cfg_subtree SRC DST          the rows of SRC re-rooted at DST (stdout; appended to $PLAN_PARAMS when that is set)
#   config_dump                  `satpipe config FILE` YAML of CFG_FLAT on stdout
#   config_main ARGS...          cmd_config: [--schema-doc] [--flat] [FILE]; exit codes 0 / 1 (errors) / 2 (no FILE)

[ -n "${_SATPIPE_CONFIG_SH:-}" ] && return 0
_SATPIPE_CONFIG_SH=1
# shellcheck source=common.sh
. "$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]}")")/common.sh"
# shellcheck source=flat.sh
. "$SATPIPE_LIB/flat.sh"

SATPIPE_SCHEMA_TSV=${SATPIPE_SCHEMA_TSV:-$SATPIPE_LIB/schema.tsv}   # test hook: another branch's schema.tsv

# Names of satpipe/tools.py TOOLS (the cross-check `tools.<k>: not a known tool`), column 1 of tools.tsv (WP4);
# tests/port/golden_config.sh checks it against tools.py.
config_known_tools() {
    "$GAWK" -F '\t' '!/^#/ && $1 != "" { printf "%s%s", (n++ ? " " : ""), $1 }' "$SATPIPE_LIB/tools.tsv"
}

# _config_abspath FILE   os.path.abspath: normpath(join(getcwd(), FILE)), no symlink resolution of FILE
_config_abspath() {
    local cwd; cwd=$(pwd -P)
    F=$1 CWD=$cwd "$GAWK" -f "$SATPIPE_LIB/normpath.awk" -e 'BEGIN { printf "%s", py_abspath(ENVIRON["F"], ENVIRON["CWD"]) }'
}

config_stream() {
    local file=$1 abs base
    abs=$(_config_abspath "$file")
    base=$(F=$abs "$GAWK" -f "$SATPIPE_LIB/normpath.awk" -e 'BEGIN { printf "%s", py_dirname(ENVIRON["F"]) }')
    if [ -d "$abs" ]; then printf '!error\t[Errno 21] Is a directory: %s\n' "$(py_repr_s "$abs")"; return 1; fi
    if [ ! -e "$abs" ]; then printf '!error\t[Errno 2] No such file or directory: %s\n' "$(py_repr_s "$abs")"; return 1; fi
    if [ ! -r "$abs" ]; then printf '!error\t[Errno 13] Permission denied: %s\n' "$(py_repr_s "$abs")"; return 1; fi
    (
        set -o pipefail
        "$GAWK" -f "$SATPIPE_LIB/yaml.awk" "$abs" \
        | BASE_DIR=$base "$GAWK" -f "$PYJSON" -f "$SATPIPE_LIB/normpath.awk" -f "$SATPIPE_LIB/validate.awk" "$SATPIPE_SCHEMA_TSV" - \
        | KNOWN_TOOLS=$(config_known_tools) "$GAWK" -f "$PYJSON" -f "$SATPIPE_LIB/crosscheck.awk" -
    )
}

py_repr_s() { S=$1 pyjson_gawk 'BEGIN { printf "%s", flat_esc(py_repr_str(ENVIRON["S"])) }'; }

config_load() {
    local file=$1 out=${2:-} stream line rc=0
    CFG_ERRORS=(); CFG_WARNINGS=()
    CFG_PATH=$(_config_abspath "$file")
    CFG_BASE_DIR=$(F=$CFG_PATH "$GAWK" -f "$SATPIPE_LIB/normpath.awk" -e 'BEGIN { printf "%s", py_dirname(ENVIRON["F"]) }')
    [ -n "$out" ] || out=$(mktemp "${TMPDIR:-/tmp}/satpipe_cfg.XXXXXX")
    CFG_FLAT=$out
    stream=$(mktemp "${TMPDIR:-/tmp}/satpipe_cfgs.XXXXXX")
    config_stream "$file" > "$stream" || rc=$?
    : > "$out"
    while IFS= read -r line; do
        case $line in
            $'!error\t'*)   CFG_ERRORS+=("$(_flat_unesc "${line#*$'\t'}")") ;;
            $'!warning\t'*) CFG_WARNINGS+=("$(_flat_unesc "${line#*$'\t'}")") ;;
            *) printf '%s\n' "$line" >> "$out" ;;
        esac
    done < "$stream"
    rm -f "$stream"
    if [ ${#CFG_ERRORS[@]} -gt 0 ]; then return 1; fi
    if [ "$rc" -ne 0 ]; then CFG_ERRORS+=("internal error: config pipeline exited $rc"); return 1; fi
    return 0
}

_flat_unesc() { S=$1 pyjson_gawk 'BEGIN { printf "%s", flat_unesc(ENVIRON["S"]) }'; }

cfg() {
    case $(flat_type "$CFG_FLAT" "$1") in
        map|list) flat_py_repr "$CFG_FLAT" "$1" ;;
        *) flat_get "$CFG_FLAT" "$1" ;;
    esac
}
cfg_type()    { flat_type "$CFG_FLAT" "$1"; }
cfg_keys()    { flat_keys "$CFG_FLAT" "$1"; }
cfg_subtree() {   # while a stage plan runs (PLAN_PARAMS set, stages.sh) the rows are appended to the params
    if [ -n "${PLAN_PARAMS:-}" ]; then flat_subtree "$CFG_FLAT" "$1" "$2" >> "$PLAN_PARAMS"; else flat_subtree "$CFG_FLAT" "$1" "$2"; fi
}
config_dump() { pyjson_gawk_file "$SATPIPE_LIB/yaml_dump.awk" "$CFG_FLAT"; }
pyjson_gawk_file() { local f=$1; shift; "$GAWK" -f "$PYJSON" -f "$f" "$@"; }

config_main() {
    local flat=0 doc=0 file="" a e
    for a in "$@"; do
        case $a in
            --flat) flat=1 ;;
            --schema-doc) doc=1 ;;
            *) file=$a ;;
        esac
    done
    if [ $doc -eq 1 ]; then "$GAWK" -f "$SATPIPE_LIB/schema_doc.awk" "$SATPIPE_SCHEMA_TSV"; return 0; fi
    if [ -z "$file" ]; then printf 'error: config path required unless --schema-doc\n' >&2; return 2; fi
    if ! config_load "$file"; then
        for e in "${CFG_ERRORS[@]}"; do printf 'error: %s\n' "$e" >&2; done
        rm -f "$CFG_FLAT"; return 1
    fi
    if [ $flat -eq 1 ]; then cat "$CFG_FLAT"; else config_dump; fi
    rm -f "$CFG_FLAT"
    return 0
}

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
    config_main "$@"; exit $?
fi
