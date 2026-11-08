# done.sh - DONE markers (PORT_SPEC 2.2, 2.3; run.py:51-66, 196-222). Sourced after common.sh and flat.sh.
#   fingerprint PATH                                       print the fingerprint (run.py:51-66); exit 1 when PATH is missing
#   params_hash PARAMS_FLAT                                sha256(json.dumps(params, sort_keys=True))[:16]
#   is_done DONE_FILE INPUTS_FILE PARAMS_FLAT              status 0 = up to date (run.py:196-209); SATPIPE_FORCE=1 forces 1
#   mark_done DONE_FILE STAGE INPUTS_FILE PARAMS_FLAT OUTPUTS_FILE   write the record, indent=1 (run.py:211-217)
# INPUTS_FILE / OUTPUTS_FILE: one path per line (as plan_* / st_* write them). PARAMS_FLAT: TYPED-flat params.

SATPIPE_SMALL=67108864      # run.py:40: files up to 64 MiB are sha256-hashed, larger ones size+mtime

_fp_file() {
    local sm s m h
    sm=$(stat -L -c '%s %Y' -- "$1") || return 1
    s=${sm% *}; m=${sm#* }
    if [ "$s" -le "$SATPIPE_SMALL" ]; then
        h=$(sha256sum < "$1") || return 1
        printf 'sha256:%s\n' "${h%% *}"
    else
        printf 'size-mtime:%s:%s\n' "$s" "$m"
    fi
}

fingerprint() {
    local p=$1
    if [ -d "$p" ]; then
        (
            shopt -s dotglob nullglob
            local e h
            # names sorted bytewise: the glob sorts with LC_ALL=C (common.sh)
            h=$(for e in "$p"/*; do
                    if [ -f "$e" ]; then printf '%s\t%s\n' "${e##*/}" "$(_fp_file "$e")" || exit 1; fi
                done | sha256sum) || exit 1
            printf 'dir-sha256:%s\n' "${h%% *}"
        )
        return
    fi
    if [ ! -e "$p" ]; then
        log "fingerprint: No such file or directory: '$p'"
        return 1
    fi
    _fp_file "$p"
}

params_hash() {
    local h
    h=$(flat_json "$1" | sha256sum) || return 1
    printf '%s\n' "${h:0:16}"
}

is_done() {
    local done=$1 inputs=$2 params=$3 rec want f fp
    [ "${SATPIPE_FORCE:-0}" = 1 ] && return 1
    [ -f "$done" ] || return 1
    # first line: params_hash (or a marker), then one "path<TAB>fingerprint" line per recorded input
    rec=$(pyjson_gawk 'BEGIN { RS = "^$" }
        { got = 1 }
        END {
            if (!got || !json_parse($0, "D") || JT["D"] != "map") { print "#UNPARSABLE"; exit }
            h = "D" SUBSEP "params_hash"
            print ((h in JT) && JT[h] == "str") ? JV[h] : "#NOHASH"
            i = "D" SUBSEP "inputs"
            if (!(i in JT)) exit
            if (JT[i] != "map") { print "#BADINPUTS"; exit }
            for (k = 1; k <= JN[i]; k++) { c = i SUBSEP JK[i, k]; print JK[i, k] "\t" (JT[c] == "str" ? JV[c] : "#NOTSTR") }
        }' "$done") || return 1
    [ "${rec%%$'\n'*}" = "#UNPARSABLE" ] && return 1
    want=$(params_hash "$params") || return 1
    [ "${rec%%$'\n'*}" = "$want" ] || return 1
    case $rec in *$'\n'*) rec=${rec#*$'\n'} ;; *) rec= ;; esac
    [ "$rec" = "#BADINPUTS" ] && return 1
    local -A recorded=() planned=()
    if [ -n "$rec" ]; then
        while IFS=$'\t' read -r f fp; do
            [ -e "$f" ] || return 1
            [ "$(fingerprint "$f")" = "$fp" ] || return 1
            recorded[$f]=1
        done <<< "$rec"
    fi
    while IFS= read -r f || [ -n "$f" ]; do
        f=${f%$'\r'}
        [ -n "$f" ] && planned[$f]=1
    done < "$inputs"
    [ "${#recorded[@]}" -eq "${#planned[@]}" ] || return 1
    for f in "${!planned[@]}"; do [ -n "${recorded[$f]:-}" ] || return 1; done
    return 0
}

# _fp_list FILE ONLY_FILES  "path<TAB>fingerprint" per line, first occurrence order (a dict keeps the first position)
_fp_list() {
    local f fp
    local -A seen=()
    while IFS= read -r f || [ -n "$f" ]; do
        f=${f%$'\r'}
        [ -n "$f" ] || continue
        [ -n "${seen[$f]:-}" ] && continue
        seen[$f]=1
        if [ "$2" = 1 ] && [ ! -f "$f" ]; then continue; fi
        case $f in *$'\t'*) die "mark_done: TAB in path: $f" ;; esac
        fp=$(fingerprint "$f") || return 1
        printf '%s\t%s\n' "$f" "$fp"
    done < "$1"
}

mark_done() {
    local done=$1 stage=$2 inputs=$3 params=$4 outputs=$5 hash fin tmp ins outs
    hash=$(params_hash "$params") || return 1
    fin=$(date +%FT%T)
    ins=$(_fp_list "$inputs" 0) || die "stage $stage: cannot fingerprint an input"
    outs=$(_fp_list "$outputs" 1) || die "stage $stage: cannot fingerprint an output"
    mkdir -p -- "$(dirname -- "$done")"
    tmp=$done.tmp.$$
    STAGE=$stage FIN=$fin HASH=$hash PARAMS=$params INS=$ins OUTS=$outs pyjson_gawk 'BEGIN {
            STAGE = ENVIRON["STAGE"]; FIN = ENVIRON["FIN"]; HASH = ENVIRON["HASH"]; PARAMS = ENVIRON["PARAMS"]
            INS = ENVIRON["INS"]; OUTS = ENVIRON["OUTS"]
            jt_init("D", "map")
            jt_set("D", "stage", "str", STAGE); jt_set("D", "finished", "str", FIN); jt_set("D", "params_hash", "str", HASH)
            jt_set("D", "params", "map", "")
            while ((getline line < PARAMS) > 0) if (!flat_line("D" SUBSEP "params", line)) { print JERR > "/dev/stderr"; exit 1 }
            jt_set("D", "inputs", "map", ""); _add("D" SUBSEP "inputs", INS)
            jt_set("D", "outputs", "map", ""); _add("D" SUBSEP "outputs", OUTS)
            printf "%s", json_indent1("D", 0)
        }
        function _add(id, s,   n, a, i, p) {
            n = split(s, a, "\n")
            for (i = 1; i <= n; i++) { if (a[i] == "") continue; p = index(a[i], "\t"); jt_set(id, substr(a[i], 1, p - 1), "str", substr(a[i], p + 1)) }
        }' > "$tmp" || { rm -f -- "$tmp"; return 1; }
    mv -f -- "$tmp" "$done"
}
