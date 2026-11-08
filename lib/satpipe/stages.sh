# stages.sh - stage registry, run_stage and the helpers stage functions use (PORT_SPEC sections 1, 2.2, 5b;
# run.py:124-249). Sourced by run.sh after common.sh, flat.sh, done.sh, config.sh and manifest.sh.
#
# Registry: STAGE_IDS = column 1 of run_stages.tsv (run.py order). Stage code lives in lib/satpipe/st_<id>.sh (one file per
# stage, sourced when present) and anchor.sh (21_anchor); SATPIPE_STAGE_FILES (colon-separated, test hook) is sourced last
# and may override them. A stage with no st_<id> function fails with "stage <id> is not ported to the bash driver".
#
# Stage contract (= a @stage function + its .plan in run.py):
#   plan_<id>   optional (no plan = inputs [] and params {}, run.py:240). Runs under `set -e` with $RUN, $CFG_FLAT set.
#               emit_inputs PATH...         planned inputs, in Python's list order (duplicates allowed, as in the list)
#               cons_files | genome_files | subx_inputs    _cons_files/_genome_files/_subx_inputs (run.py:254-265)
#               cfg_subtree SRC DST         params[DST] = config node SRC (leaf or block, config order kept)
#               param_cfg_or_null SRC DST   same, or params[DST] = None when SRC is absent (dict.get)
#               param_set DST TYPE TEXT     a literal: param_set step3b str runs-v3; param_set members null None
#               Params rows must be emitted in the order of the Python dict literal (DONE files keep insertion order).
#   st_<id>     the body; runs under `set -e` in the driver shell (never call it in a subshell). Helpers:
#               cfg PATH / cfg_type PATH / cfg_keys PATH   config.sh (py_str text, type, child keys)
#               tool NAME                   Ctx.tool: tools.NAME from the config or NAME
#               ctx_d PART...               Ctx.d: mkdir -p $RUN/PART/...; prints the path
#               sh_run SID CWD [log_name=X] [VAR=value ...] -- ARGV...     Ctx.sh (CWD "" = none); fails the run on exit != 0
#               sh_capture OUTVAR SID [VAR=value ...] -- ARGV...           subprocess.run(capture_output) (run.py:638): stdout
#                                           into OUTVAR, stderr dropped, exit status ignored, no log, no say
#               sh_to FILE SID [VAR=value ...] -- ARGV...                  subprocess.run(stdout=f, check=True) (run.py:702)
#               say MSG                     Ctx.say
#               manifest_set KEY JSON       ctx.manifest[KEY] = JSON (top level; position kept when KEY exists)
#               stage_die MSG               raise StageError(MSG): `satpipe run: MSG` after the manifest is written, exit 1
#               emit_outputs PATH...        the returned output list (only regular files are fingerprinted, run.py:216)
#   VAR=value texts are Python str() texts (cfg gives them); the child environment is the caller's environment
#   (bin/satpipe's SATPIPE_USER_PATH / _LC_ALL / _TMPDIR restored, driver internals removed) + LC_ALL=C + PATH=bin_dir:PATH
#   + SATPIPE_TOOL_<NAME> + the VAR=values, as Ctx.env (run.py:158-169).
#   Parallel jobs (07_mafft): run `( STAGE_ERRFILE=<own file>; sh_run ... ) &` per job, wait for all, then
#   stage_die_file <first failed job's file> in submission order (ThreadPoolExecutor.map semantics, run.py:557-560).
#
# Trace (PORT_SPEC 5b): SATPIPE_TRACE=<file> appends PLAN<TAB>sid<TAB>json inputs<TAB>canonical params<TAB>hash,
# CMD<TAB>sid<TAB>log_name<TAB>cwd<TAB>json argv<TAB>json sorted [[k, v], ...] of env entries that differ from the caller's
# environment (log_name "" for sh_capture/sh_to), OUT<TAB>sid<TAB>json outputs. SATPIPE_DRY=1: commands are not run
# (success, empty output), is_done is false, no DONE is written.

[ -n "${_SATPIPE_STAGES_SH:-}" ] && return 0
_SATPIPE_STAGES_SH=1

STAGE_IDS=()
while IFS= read -r _l || [ -n "$_l" ]; do
    _l=${_l%$'\r'}; _l=${_l%%$'\t'*}
    case $_l in ''|'#'*) continue ;; esac
    STAGE_IDS+=("$_l")
done < "$SATPIPE_LIB/run_stages.tsv"
unset _l

stages_load() {
    local sid f
    for sid in "${STAGE_IDS[@]}"; do
        f=$SATPIPE_LIB/st_$sid.sh
        if [ -f "$f" ]; then . "$f"; fi
    done
    . "$SATPIPE_LIB/anchor.sh"
    if [ -n "${SATPIPE_STAGE_FILES:-}" ]; then
        local -a extra=()
        IFS=: read -r -a extra <<< "$SATPIPE_STAGE_FILES"
        for f in "${extra[@]}"; do [ -n "$f" ] && . "$f"; done
    fi
    return 0
}

# ------------------------------------------------------------------ plan helpers
emit_inputs()  { [ $# -eq 0 ] || printf '%s\n' "$@" >> "$PLAN_INPUTS"; }
emit_outputs() { [ $# -eq 0 ] || printf '%s\n' "$@" >> "$STAGE_OUTPUTS"; }
param_set()    { flat_emit "$1" "$2" "$3" >> "$PLAN_PARAMS"; }
param_cfg_or_null() { if cfg_type "$1" >/dev/null 2>&1; then cfg_subtree "$1" "$2"; else param_set "$2" null None; fi; }
genome_files() { local k; while IFS= read -r k; do emit_inputs "$(cfg "genomes/$k")"; done < <(cfg_keys genomes); }
cons_files()   { local k; while IFS= read -r k; do emit_inputs "$(cfg "families/$k/consensus")"; done < <(cfg_keys families); }
subx_inputs()  {
    local c
    case $(cfg_type subfamily/consensuses 2>/dev/null) in path|str) c=$(cfg subfamily/consensuses) ;; *) return 0 ;; esac
    if [ -n "$c" ] && [ -e "$c" ]; then emit_inputs "$c"; fi
    return 0
}

# ------------------------------------------------------------------ body helpers
ctx_d() { local p=$RUN a; for a in "$@"; do p=$p/$a; done; mkdir -p -- "$p" || stage_die "cannot create directory $p"; printf '%s\n' "$p"; }

tool() {
    case $(cfg_type "tools/$1" 2>/dev/null) in
        path|str) local v; v=$(cfg "tools/$1"); if [ -n "$v" ]; then printf '%s\n' "$v"; return 0; fi ;;
    esac
    printf '%s\n' "$1"
}

stage_die() {
    if [ -n "${STAGE_ERRFILE:-}" ]; then printf '%s' "$1" > "$STAGE_ERRFILE"; exit 1; fi
    die "satpipe run: $1"
}
stage_die_file() { local m; m=$(cat -- "$1"; printf x); stage_die "${m%x}"; }

# _env_setup   once per run: the tool part of Ctx.env (bin_dir, SATPIPE_TOOL_*) from the config
_env_setup() {
    local k v t
    ENV_BIN_DIR=""; ENV_TOOLS=()
    while IFS= read -r k; do
        t=$(cfg_type "tools/$k"); [ "$t" = path ] || [ "$t" = str ] || continue
        v=$(cfg "tools/$k")
        if [ "$k" = bin_dir ]; then ENV_BIN_DIR=$v; else ENV_TOOLS+=("SATPIPE_TOOL_${k^^}=$v"); fi
    done < <(cfg_keys tools 2>/dev/null)
    return 0
}

# _env_final VAR=value...   ENV_FINAL_K / ENV_FINAL_V: the variables Ctx.env sets, in order, later ones replacing earlier
_env_final() {
    local -A pos=()
    local a k i
    ENV_FINAL_K=(); ENV_FINAL_V=()
    _env_put() { if [ -n "${pos[$1]+x}" ]; then ENV_FINAL_V[${pos[$1]}]=$2; else pos[$1]=${#ENV_FINAL_K[@]}; ENV_FINAL_K+=("$1"); ENV_FINAL_V+=("$2"); fi; }
    _env_put LC_ALL C
    if [ -n "${ENV_BIN_DIR:-}" ]; then
        if _env_base PATH; then _env_put PATH "$ENV_BIN_DIR:$REPLY"; else _env_put PATH "$ENV_BIN_DIR:"; fi
    fi
    for a in "${ENV_TOOLS[@]}" "$@"; do k=${a%%=*}; _env_put "$k" "${a#*=}"; done
}

# _env_base NAME   REPLY = NAME's value in the caller's environment; status 1 when it was not set there
_SATPIPE_INTERNAL='SATPIPE_LIB SATPIPE_ROOT SCRIPTS PYJSON GAWK SATPIPE_USER_PATH_SET SATPIPE_USER_PATH SATPIPE_USER_LC_ALL_SET SATPIPE_USER_LC_ALL SATPIPE_USER_TMPDIR_SET SATPIPE_USER_TMPDIR'
_env_base() {
    REPLY=""
    case $1 in
        PATH|LC_ALL|TMPDIR)
            local s=SATPIPE_USER_$1_SET u=SATPIPE_USER_$1
            if [ -n "${!s+x}" ]; then [ "${!s}" = 1 ] || return 1; REPLY=${!u}; return 0; fi ;;
        SATPIPE_LIB|SATPIPE_ROOT|SCRIPTS|PYJSON|GAWK|SATPIPE_USER_*) return 1 ;;
    esac
    local v; v=$(printenv -- "$1" && printf x) || return 1
    v=${v%x}; REPLY=${v%$'\n'}
    return 0
}

# _env_restore   (in the child subshell) the caller's environment: driver internals out, PATH/LC_ALL/TMPDIR back
_env_restore() {
    local n
    if [ -n "${SATPIPE_USER_PATH_SET+x}" ]; then
        if [ "$SATPIPE_USER_PATH_SET" = 1 ]; then export PATH=$SATPIPE_USER_PATH; else unset PATH; fi
        if [ "$SATPIPE_USER_TMPDIR_SET" = 1 ]; then export TMPDIR=$SATPIPE_USER_TMPDIR; else unset TMPDIR; fi
    fi
    for n in $_SATPIPE_INTERNAL; do unset "$n"; done
}

# _sh_child CWD -- ARGV   run ARGV with the Ctx.env environment built by _env_final (in a subshell)
_sh_child() {
    local cwd=$1; shift; [ "${1:-}" = -- ] && shift
    (
        _env_restore
        local i
        for i in "${!ENV_FINAL_K[@]}"; do export "${ENV_FINAL_K[$i]}=${ENV_FINAL_V[$i]}"; done
        if [ -n "$cwd" ]; then cd -- "$cwd" || exit 1; fi
        exec "$@"
    )
}

# ------------------------------------------------------------------ trace
_json_list() {
    if [ $# -eq 0 ]; then printf '[]'; return 0; fi
    printf '%s\0' "$@" | "$GAWK" -f "$PYJSON" -e 'BEGIN { RS = "\0" } { a[++n] = $0 }
        END { printf "["; for (i = 1; i <= n; i++) printf "%s%s", (i > 1 ? ", " : ""), json_str(a[i]); printf "]" }'
}
_json_list_file() { local -a a=(); local l; while IFS= read -r l || [ -n "$l" ]; do a+=("${l%$'\r'}"); done < "$1"; _json_list "${a[@]}"; }
_trace() { [ -n "${SATPIPE_TRACE:-}" ] || return 0; local IFS=$'\t'; printf '%s\n' "$*" >> "$SATPIPE_TRACE"; }
_dry() { [ "${SATPIPE_DRY:-}" = 1 ]; }

# _trace_cmd SID LOG_NAME CWD -- ARGV   (uses ENV_FINAL_K/V)
_trace_cmd() {
    [ -n "${SATPIPE_TRACE:-}" ] || return 0
    local sid=$1 ln=$2 cwd=$3 i k argv envj; shift 3; [ "${1:-}" = -- ] && shift
    argv=$(_json_list "$@")
    local -a d=()
    for i in "${!ENV_FINAL_K[@]}"; do
        k=${ENV_FINAL_K[$i]}
        if _env_base "$k" && [ "$REPLY" = "${ENV_FINAL_V[$i]}" ]; then continue; fi
        d+=("$k" "${ENV_FINAL_V[$i]}")
    done
    if [ ${#d[@]} -eq 0 ]; then envj='[]'; else
        envj=$(printf '%s\0' "${d[@]}" | "$GAWK" -f "$PYJSON" -e 'BEGIN { RS = "\0" }
            NR % 2 { k = $0; next } { v[k] = $0 }
            END { PROCINFO["sorted_in"] = "@ind_str_asc"; printf "["
                  for (k in v) printf "%s[%s, %s]", (n++ ? ", " : ""), json_str(k), json_str(v[k]); printf "]" }')
    fi
    _trace CMD "$sid" "$ln" "$cwd" "$argv" "$envj"
}

# ------------------------------------------------------------------ commands
# _sh_opts ARGS...   splits [log_name=X] [VAR=value ...] -- ARGV; sets SH_LOG (or ""), SH_ENV=(), SH_ARGV=()
_sh_opts() {
    SH_LOG=""; SH_ENV=(); SH_ARGV=()
    while [ $# -gt 0 ] && [ "$1" != -- ]; do
        case $1 in
            log_name=*) SH_LOG=${1#log_name=} ;;
            [A-Za-z_]*=*) SH_ENV+=("$1") ;;
            *) stage_die "sh_run: bad argument before --: $1" ;;
        esac
        shift
    done
    [ $# -gt 0 ] && shift
    SH_ARGV=("$@")
    [ ${#SH_ARGV[@]} -gt 0 ] || stage_die "sh_run: no command"
}

sh_run() {
    local sid=$1 cwd=$2 logf line rc=0 tl
    shift 2
    _sh_opts "$@"
    logf=$(ctx_d logs)/${SH_LOG:-$sid}.log
    line=$(local IFS=' '; printf '%s' "${SH_ARGV[*]}")
    say "  \$ $line${cwd:+   (cwd $cwd)}"
    printf '\n### %s  %s\n$ %s\n' "$(date +%FT%T)" "$cwd" "$line" >> "$logf"
    _env_final "${SH_ENV[@]}"
    _trace_cmd "$sid" "${SH_LOG:-$sid}" "$cwd" -- "${SH_ARGV[@]}"
    _dry && return 0
    _sh_child "$cwd" -- "${SH_ARGV[@]}" >> "$logf" 2>&1 || rc=$?
    if [ $rc -ne 0 ]; then
        tl=$(tail -n 15 -- "$logf"; printf x); tl=${tl%x}
        stage_die "stage $sid: command failed (exit $rc): $line"$'\n'"--- last lines of $logf ---"$'\n'"$tl"
    fi
    return 0
}

sh_capture() {
    local __var=$1 sid=$2 __out
    shift 2
    _sh_opts "$@"
    _env_final "${SH_ENV[@]}"
    _trace_cmd "$sid" "" "" -- "${SH_ARGV[@]}"
    if _dry; then printf -v "$__var" '%s' ""; return 0; fi
    __out=$(_sh_child "" -- "${SH_ARGV[@]}" 2>/dev/null; printf x) || true
    printf -v "$__var" '%s' "${__out%x}"
    return 0
}

sh_to() {
    local file=$1 sid=$2 rc=0 line
    shift 2
    _sh_opts "$@"
    _env_final "${SH_ENV[@]}"
    _trace_cmd "$sid" "" "" -- "${SH_ARGV[@]}"
    if _dry; then : > "$file"; return 0; fi
    _sh_child "" -- "${SH_ARGV[@]}" > "$file" || rc=$?
    if [ $rc -ne 0 ]; then
        line=$(local IFS=' '; printf '%s' "${SH_ARGV[*]}")
        stage_die "stage $sid: command failed (exit $rc): $line"     # Python: CalledProcessError traceback
    fi
    return 0
}

# ------------------------------------------------------------------ run_stage (run.py:238-248)
# _stage_plan SID   fills $PLAN_INPUTS / $PLAN_PARAMS (files must exist and be empty); sets PLAN_HASH
_stage_plan() {
    local sid=$1
    declare -F "st_$sid" >/dev/null || stage_die "stage $sid is not ported to the bash driver"
    if declare -F "plan_$sid" >/dev/null; then set -e; "plan_$sid"; set +e; fi
    PLAN_HASH=$(params_hash "$PLAN_PARAMS") || stage_die "stage $sid: bad params"
}

_stage_plan_record() { printf 'PLAN\t%s\t%s\t%s\t%s\n' "$1" "$(_json_list_file "$PLAN_INPUTS")" "$(flat_json "$PLAN_PARAMS")" "$PLAN_HASH"; }

run_stage() {
    local sid=$1
    CUR_STAGE=$sid
    local PLAN_INPUTS=$RUN_TMP/inputs PLAN_PARAMS=$RUN_TMP/params STAGE_OUTPUTS=$RUN_TMP/outputs
    : > "$PLAN_INPUTS"; : > "$PLAN_PARAMS"; : > "$STAGE_OUTPUTS"
    _stage_plan "$sid"
    if [ -n "${SATPIPE_TRACE:-}" ]; then _stage_plan_record "$sid" >> "$SATPIPE_TRACE"; fi
    if ! _dry && is_done "$RUN/$sid/DONE" "$PLAN_INPUTS" "$PLAN_PARAMS"; then
        say "[$sid] up to date, skipped"
        return 0
    fi
    say "[$sid] running"
    local pp=$PLAN_PARAMS
    PLAN_PARAMS=""                       # cfg_subtree prints to stdout outside plan_*
    set -e; "st_$sid"; set +e
    PLAN_PARAMS=$pp
    _trace OUT "$sid" "$(_json_list_file "$STAGE_OUTPUTS")"
    if ! _dry; then
        mark_done "$RUN/$sid/DONE" "$sid" "$PLAN_INPUTS" "$PLAN_PARAMS" "$STAGE_OUTPUTS" || stage_die "stage $sid: cannot write DONE"
        manifest_set_file stages "$sid" "$RUN/$sid/DONE"
    fi
    say "[$sid] done"
    return 0
}
