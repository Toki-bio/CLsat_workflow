# run.sh - `satpipe run` (cli.py:67-88 cmd_run + run.py:1106-1158 run()). Sourced by bin/satpipe after check.sh.
#   run_main CONFIG   reads AP_from_stage, AP_to_stage, AP_only, AP_force, AP_skip_tool_check, AP_plan (bin/satpipe)
# --plan (bash only, PORT_SPEC 5b): after the check, print the PLAN record of every selected stage and exit; nothing is written.
# Errors: StageError texts as `satpipe run: MSG` on stderr, exit 1, after the manifest is written (the `finally` of run()).

[ -n "${_SATPIPE_RUN_SH:-}" ] && return 0
_SATPIPE_RUN_SH=1
. "$SATPIPE_LIB/flat.sh"
. "$SATPIPE_LIB/done.sh"
. "$SATPIPE_LIB/manifest.sh"
. "$SATPIPE_LIB/stages.sh"

_RUN_LIVE=0           # 1 inside run()'s try/finally: the manifest is written on any exit
RUN_TMP=""
STAGE_ERRFILE=""
CUR_STAGE=""

_run_cleanup() { [ -n "$RUN_TMP" ] && rm -rf -- "$RUN_TMP"; [ -n "${MAN_JSON:-}" ] && rm -f -- "$MAN_JSON" "$MAN_JSON.tmp"; check_cleanup; }

# EXIT trap while stages run: write the manifest (finally), then report the error as cmd_run does
_run_on_exit() {
    local rc=$? msg
    trap - EXIT INT TERM
    set +e
    if [ "$_RUN_LIVE" = 1 ]; then manifest_write "$RUN/manifest.json"; fi
    if [ -n "$STAGE_ERRFILE" ] && [ -f "$STAGE_ERRFILE" ]; then
        msg=$(cat -- "$STAGE_ERRFILE"; printf x); msg=${msg%x}
        printf 'satpipe run: %s\n' "$msg" >&2
        rc=1
    elif [ "$rc" -ne 0 ]; then
        printf 'satpipe run: stage %s: failed (exit %d)\n' "${CUR_STAGE:-?}" "$rc" >&2   # Python: a traceback
        rc=1
    fi
    _run_cleanup
    exit "$rc"
}

# _run_select   SEL=() selected stage ids (run.py:1133-1150); prints the error and returns 1 like StageError/ValueError
_run_select() {
    local -a sel=("${STAGE_IDS[@]}") only=()
    local i s found
    _idx() { local j; for ((j = 0; j < ${#sel[@]}; j++)); do [ "${sel[$j]}" = "$1" ] && { REPLY=$j; return 0; }; done; return 1; }
    _known() { local j; for j in "${STAGE_IDS[@]}"; do [ "$j" = "$1" ] && return 0; done; return 1; }
    _bad() { printf 'satpipe run: unknown stage %s; stages: %s\n' "$(S=$1 pyjson_gawk 'BEGIN { printf "%s", py_repr_str(ENVIRON["S"]) }')" "${STAGE_IDS[*]}" >&2; }
    if [ -n "${AP_from_stage:-}" ]; then
        _known "$AP_from_stage" || { _bad "$AP_from_stage"; return 1; }
        _idx "$AP_from_stage"; sel=("${sel[@]:$REPLY}")
    fi
    if [ -n "${AP_to_stage:-}" ]; then
        _known "$AP_to_stage" || { _bad "$AP_to_stage"; return 1; }
        if ! _idx "$AP_to_stage"; then   # --to before --from: list.index raises ValueError (a traceback in Python)
            printf 'satpipe run: ValueError: list.index(x): x not in list\n' >&2   # the traceback's last line
            return 1
        fi
        sel=("${sel[@]:0:$((REPLY + 1))}")
    fi
    if [ -n "${AP_only:-}" ]; then
        IFS=, read -r -a only <<< "$AP_only"
        local -a keep=()
        for s in "${sel[@]}"; do found=0; for i in "${only[@]}"; do [ "$i" = "$s" ] && found=1; done; [ $found = 1 ] && keep+=("$s"); done
        sel=("${keep[@]}")
    fi
    SEL=("${sel[@]}")
    # seeds given for every family: skip discovery (run.py:1146-1150)
    SEEDS_SKIP=0
    if [ "$(cfg_type inputs/seeds 2>/dev/null)" = map ] && [ "$(cfg inputs/seeds 2>/dev/null)" != "{}" ]; then
        local -A given=()
        local f
        while IFS= read -r f; do given[$f]=1; done < <(cfg_keys inputs/seeds)
        SEEDS_SKIP=1
        while IFS= read -r f; do [ -n "${given[$f]:-}" ] || SEEDS_SKIP=0; done < <(cfg_keys families)
    fi
    if [ $SEEDS_SKIP = 1 ]; then
        local -a keep=()
        for s in "${SEL[@]}"; do case $s in 03_motifs|04_filter|05_rotate|06_cluster|07_mafft|08_cross|09_assign) ;; *) keep+=("$s") ;; esac; done
        SEL=("${keep[@]}")
    fi
    return 0
}

_in_sel() { local s; for s in "${SEL[@]}"; do [ "$s" = "$1" ] && return 0; done; return 1; }

# _run_tools_json   manifest.tools: {name: {path, version, error}} (run.py:1124, tools.py:80-108)
_run_tools_json() {
    local -A ov=()
    local k t
    while IFS= read -r k; do
        t=$(cfg_type "tools/$k"); case $t in path|str) ov[$k]=$(cfg "tools/$k") ;; esac
    done < <(cfg_keys tools 2>/dev/null)
    local name req stages note found path ver err out="{" sep=""
    while IFS=$'\037' read -r name req stages note found path ver err; do
        out+=$sep$(N=$name F=$found P=$path V=$ver E=$err O=${ov[$name]:-} pyjson_gawk 'BEGIN {
            if (ENVIRON["F"] != 1) { p = "null"; v = "null"; o = ENVIRON["O"]; e = json_str("not found" (o != "" ? " (override " o ")" : "")) }
            else { p = json_str(flat_unesc(ENVIRON["P"])); v = flat_unesc(ENVIRON["V"]); v = (v == "") ? "null" : json_str(v)
                   e = flat_unesc(ENVIRON["E"]); e = (e == "") ? "null" : json_str(e) }
            printf "%s: {\"path\": %s, \"version\": %s, \"error\": %s}", json_str(ENVIRON["N"]), p, v, e }')
        sep=", "
    done < <(tools_probe_all "$CFG_FLAT" | tools_tsv_us /dev/stdin)
    printf '%s}' "$out"
}

run_main() {
    local cfgfile=$1 sid c
    SATPIPE_FORCE=${AP_force:-0}; export SATPIPE_FORCE
    check_run "$cfgfile" 1 "$((1 - ${AP_skip_tool_check:-0}))"
    if [ "$CHECK_FAILED" = 1 ]; then
        check_render text >&2
        printf 'satpipe run: refusing to start, fix the check failures above (or --skip-tool-check)\n' >&2
        check_cleanup; return 1
    fi
    stages_load
    RUN=$(cfg out_dir)
    RUN_TMP=$(mktemp -d "${TMPDIR:-/tmp}/satpipe_run.XXXXXX")
    STAGE_ERRFILE=$RUN_TMP/stage_error
    _env_setup

    if [ "${AP_plan:-0}" = 1 ]; then                  # plan only: no directory, log or manifest
        _run_select || { _run_cleanup; return 1; }
        trap _run_on_exit EXIT
        local PLAN_INPUTS=$RUN_TMP/inputs PLAN_PARAMS=$RUN_TMP/params
        for sid in "${STAGE_IDS[@]}"; do
            _in_sel "$sid" || continue
            CUR_STAGE=$sid
            : > "$PLAN_INPUTS"; : > "$PLAN_PARAMS"
            _stage_plan "$sid"
            _stage_plan_record "$sid"
        done
        trap - EXIT
        _run_cleanup; return 0
    fi

    mkdir -p -- "$RUN/logs" || { _run_cleanup; return 1; }
    SATPIPE_LOG=$RUN/logs/satpipe.log
    : >> "$SATPIPE_LOG"
    trap _run_on_exit EXIT
    trap 'exit 130' INT
    trap 'exit 143' TERM
    say "=== satpipe $SATPIPE_VERSION run $(cfg run_name)  $(date +%FT%T)"
    manifest_init "$RUN/manifest.json"
    manifest_set satpipe_version "$(S=$SATPIPE_VERSION pyjson_gawk 'BEGIN { printf "%s", json_str(ENVIRON["S"]) }')"
    manifest_set config_file "$(S=$CFG_PATH pyjson_gawk 'BEGIN { printf "%s", json_str(ENVIRON["S"]) }')"
    manifest_set_flat config "$CFG_FLAT"
    manifest_set tools "$(_run_tools_json)"
    c=$(git -C "$SATPIPE_ROOT" rev-parse HEAD 2>/dev/null) || c=""
    c=$(S=$c pyjson_gawk 'BEGIN { s = ENVIRON["S"]; gsub(/^[ \t\n\r\f\v]+|[ \t\n\r\f\v]+$/, "", s); printf "%s", (s == "" ? "null" : json_str(s)) }')
    manifest_set code_commit "$c"
    manifest_append runs "$(date +%FT%T)"
    if ! _run_select; then                            # raised before run()'s try: no manifest written
        trap - EXIT INT TERM
        _run_cleanup; return 1
    fi
    [ "$SEEDS_SKIP" = 1 ] && say "seeds given for every family: skipping 03_motifs 04_filter 05_rotate 06_cluster 07_mafft 08_cross 09_assign"
    _RUN_LIVE=1
    for sid in "${STAGE_IDS[@]}"; do
        _in_sel "$sid" || continue
        run_stage "$sid"
    done
    CUR_STAGE=""
    _RUN_LIVE=0
    manifest_write "$RUN/manifest.json"
    trap - EXIT INT TERM
    printf 'summary: %s\n' "$(manifest_summary centroids_n centroids_md5 seed_counts master_loci segmented_loci catalogue_loci ortholog_pairs)"
    _run_cleanup
    return 0
}
