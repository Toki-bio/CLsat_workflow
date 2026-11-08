# check.sh - `satpipe check` (check.py; cli.py:16-22) on top of config.sh, tools.sh and check.awk.
# Sourced after common.sh (PORT_SPEC sections 1, 4 "check").
#
# API
#   check_run CONFIG DO_FILES DO_TOOLS   run_check (check.py:250-277). Writes the report items to $CHECK_ITEMS (a mktemp
#                                        file; lines LEVEL<TAB>SECTION<TAB>MESSAGE, message flat-escaped) and sets
#                                        CHECK_FAILED (0/1). When the config loaded, CFG_FLAT etc. stay set (config_load).
#   check_render text|json               Report.to_text() / json.dumps(to_json(), indent=2) of $CHECK_ITEMS on stdout
#   check_cleanup                        removes $CHECK_ITEMS and $CFG_FLAT
#   check_main CONFIG NO_FILES NO_TOOLS JSON   cmd_check: prints the report, returns 1 when it failed, else 0
#
# Test hook: SATPIPE_MAX_BASE_PATH (default 70) replaces check.py MAX_BASE_PATH (check.py:228), which the Python tests
# monkeypatch (conftest.py:12-17, test_check.py:112-117).

[ -n "${_SATPIPE_CHECK_SH:-}" ] && return 0
_SATPIPE_CHECK_SH=1
# shellcheck source=config.sh
. "$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]}")")/config.sh"
# shellcheck source=tools.sh
. "$SATPIPE_LIB/tools.sh"

SATPIPE_MAX_BASE_PATH=${SATPIPE_MAX_BASE_PATH:-70}

_ck() { "$GAWK" -f "$PYJSON" -f "$SATPIPE_LIB/check.awk" -e "BEGIN { $1() }" < /dev/null; }
_ck_unesc() { S=$1 "$GAWK" -f "$PYJSON" -e 'BEGIN { printf "%s", flat_unesc(ENVIRON["S"]) }' < /dev/null; }
# _ck_u S   REPLY = flat_unesc(S), trailing newlines kept (a $(...) would drop them)
_ck_u() { REPLY=$(S=$1 "$GAWK" -f "$PYJSON" -e 'BEGIN { printf "%s.", flat_unesc(ENVIRON["S"]) }' < /dev/null); REPLY=${REPLY%.}; }

# _ck_add LEVEL SECTION MESSAGE
_ck_add() {
    local m=$3
    m=${m//\\/\\\\}; m=${m//$'\t'/\\t}; m=${m//$'\n'/\\n}; m=${m//$'\r'/\\r}
    printf '%s\t%s\t%s\n' "$1" "$2" "$m" >> "$CHECK_ITEMS"
}

# _ck_exists SECTION LABEL PATH [dir] [nonfatal]   _exists (check.py:73-79)
_ck_exists() {
    local sec=$1 label=$2 p=$3 kind=${4:-file} fatal=${5:-1} ok=1
    if [ "$kind" = dir ]; then [ -d "$p" ] || ok=0; else [ -f "$p" ] || ok=0; fi
    if [ $ok = 1 ]; then _ck_add OK "$sec" "$label: $p"; return 0; fi
    _ck_add "$([ "$fatal" = 1 ] && echo FAIL || echo WARN)" "$sec" "$label: $kind not found: $p"
    return 1
}

# config accessors that tolerate absent/null nodes: _ck_val PATH prints the text of a non-empty str/path leaf, status 1 otherwise
_ck_val() {
    case $(cfg_type "$1" 2>/dev/null) in
        str|path) local v; v=$(cfg "$1"); [ -n "$v" ] || return 1; printf '%s' "$v" ;;
        *) return 1 ;;
    esac
}
# Map keys come out of the flat file escaped (flat_esc: a key may end in a newline, since re.match's `$` accepts one);
# lookups use the escaped key, messages and comparisons the real one (_ck_unesc).
_ck_keys() { [ "$(cfg_type "$1" 2>/dev/null)" = map ] && cfg_keys "$1"; }
# escaped keys of a map in Python sorted() order of the real keys (bytewise under LC_ALL=C)
_ck_sorted_keys() {
    _ck_keys "$1" | "$GAWK" -f "$PYJSON" -e '{ u[NR] = flat_unesc($0); e[flat_unesc($0)] = $0 }
        END { n = asort(u, s, "@val_str_asc"); for (i = 1; i <= n; i++) print e[s[i]] }'
}
# real keys of a map joined by \037 (the CODES / FAMS lists handed to check.awk)
_ck_ulist() { _ck_keys "$1" | "$GAWK" -f "$PYJSON" -e '{ printf "%s%s", (NR > 1 ? "\037" : ""), flat_unesc($0) }'; }

# _ck_fasta PATH MAX   read_fasta rows: name<TAB>length (name unescaped), one per line
_ck_fasta() { P=$1 MAX=$2 _ck ck_fasta_main; }

_ck_bed() {   # check_bed PATH LABEL CODES SECTION
    _ck_exists "$4" "$2" "$1" || return 0
    P=$1 LABEL=$2 CODES=$3 SECTION=$4 _ck ck_bed_main >> "$CHECK_ITEMS"
}

_ck_files() {
    local codes code g fam p n name ln mono rot names missing out parent sub what len T=$'\t'
    local k
    codes=$(_ck_ulist genomes)
    while IFS= read -r k; do
        [ -n "$k" ] || continue
        _ck_u "$k"; code=$REPLY; g=$(cfg "genomes/$k")
        if _ck_exists genomes "$code" "$g"; then
            if [ -f "$g.fai" ]; then _ck_add OK genomes "$code: index present ($g.fai)"
            else _ck_add WARN genomes "$code: no .fai index; stage 11 will run \`samtools faidx\` (needs write access next to the genome)"; fi
        fi
    done < <(_ck_sorted_keys genomes)
    while IFS= read -r k; do
        [ -n "$k" ] || continue
        _ck_u "$k"; _ck_exists inputs.trf_dat "$REPLY" "$(cfg "inputs/trf_dat/$k")" || true
    done < <(_ck_sorted_keys inputs/trf_dat)
    while IFS= read -r k; do
        [ -n "$k" ] || continue
        _ck_u "$k"; _ck_bed "$(cfg "inputs/seeds/$k")" "seeds $REPLY" "$codes" inputs.seeds
    done < <(_ck_sorted_keys inputs/seeds)
    while IFS= read -r k; do
        [ -n "$k" ] || continue
        _ck_u "$k"; fam=$REPLY; p=$(cfg "families/$k/consensus")
        if _ck_exists families "$fam consensus" "$p"; then
            local -a recs=()
            mapfile -t recs < <(_ck_fasta "$p" 5)
            if [ ${#recs[@]} -ne 1 ]; then
                _ck_add FAIL families "$fam: consensus must hold exactly one record, found ${#recs[@]}"
            else
                name=$(_ck_unesc "${recs[0]%$T*}"); ln=${recs[0]##*$T}; mono=$(cfg "families/$k/mono")
                if [ "$ln" != "$mono" ]; then
                    _ck_add WARN families "$fam: consensus length $ln != mono $mono (mono drives every window size; check which is intended)"
                else
                    _ck_add OK families "$fam: consensus $name, $ln bp = mono"
                fi
            fi
        fi
    done < <(_ck_sorted_keys families)
    if rot=$(_ck_val references/rotated_consensus); then
        if _ck_exists references rotated_consensus "$rot"; then
            local -A have=()
            while IFS= read -r n; do have[x$(_ck_unesc "${n%$T*}")]=1; done < <(_ck_fasta "$rot" 1000)
            missing=""
            while IFS= read -r k; do
                _ck_u "$k"; fam=$REPLY
                [ -n "${have[x$fam]+x}" ] || missing+=${missing:+, }$fam
            done < <(_ck_keys families)
            if [ -n "$missing" ]; then
                _ck_add WARN references "rotated_consensus lacks records named like families: $missing (Step 9 identifies references by header; the Darevskia file uses CLsat1..4)"
            fi
        fi
    else
        _ck_add WARN references "rotated_consensus not given: the driver will derive it by Step-5 rotation of each consensus and record the deviation (the published run used CLsat_cons.rot)"
    fi
    P=$(cfg family_assignment)
    if _ck_exists family_assignment table "$P"; then
        P=$P FAMS=$(_ck_ulist families) _ck ck_famtab_main >> "$CHECK_ITEMS"
    fi
    if p=$(_ck_val subfamily/consensuses); then
        local -A sn=()
        local nsn=0 f b m lost
        if [ -d "$p" ]; then
            local -a ent=()
            local old; old=$(shopt -p nullglob dotglob || true)
            shopt -s nullglob dotglob; ent=("$p"/*); eval "$old"
            n=0
            for f in "${ent[@]}"; do
                b=${f##*/}
                if [[ $b =~ \.(fa|fas|fasta|fna)$ ]]; then n=$((n + 1)); b=${b%.*}; [ -n "${sn[x$b]+x}" ] || { sn[x$b]=1; nsn=$((nsn + 1)); }; fi
            done
            if [ $n -gt 0 ]; then _ck_add OK subfamily "consensuses: directory with $n FASTA files"
            else _ck_add FAIL subfamily "consensuses: directory has no FASTA files: $p"; fi
        elif [ -f "$p" ]; then
            while IFS= read -r n; do
                b=$(_ck_unesc "${n%$T*}"); [ -n "${sn[x$b]+x}" ] || { sn[x$b]=1; nsn=$((nsn + 1)); }
            done < <(_ck_fasta "$p" 100000)
            _ck_add OK subfamily "consensuses: FASTA with $nsn records"
        else
            _ck_add WARN subfamily "consensuses: not found: $p (the sub-subfamily stage will write an empty table)"
        fi
        if [ $nsn -gt 0 ]; then
            while IFS= read -r k; do
                [ -n "$k" ] || continue
                [ "$(cfg_type "subfamily/members/$k")" = list ] || continue        # a null list makes Python raise TypeError
                lost=""
                while IFS= read -r m; do
                    m=$(cfg "subfamily/members/$k/$m")
                    [ -n "${sn[x$m]+x}" ] || lost+=${lost:+, }$m
                done < <(cfg_keys "subfamily/members/$k")
                _ck_u "$k"; [ -z "$lost" ] || _ck_add WARN subfamily "members of $REPLY without a consensus record/file: $lost"
            done < <(_ck_keys subfamily/members)
        fi
    else
        _ck_add OK subfamily "no sub-subfamily consensuses: stage skipped"
    fi
    if p=$(_ck_val manual_loci); then _ck_bed "$p" manual_loci "$codes" manual_loci; fi
    out=$(cfg out_dir)
    if [ -d "$out" ]; then
        local -a ent=()
        local old; old=$(shopt -p nullglob dotglob || true)
        shopt -s nullglob dotglob; ent=("$out"/*); eval "$old"
        if [ ${#ent[@]} -gt 0 ]; then _ck_add OK out_dir "exists: $out (non-empty; stages with DONE markers would be skipped)"
        else _ck_add OK out_dir "exists: $out"; fi
    else
        parent=$(F=${out%"${out##*[!/]}"} "$GAWK" -f "$SATPIPE_LIB/normpath.awk" -e 'BEGIN { printf "%s", py_dirname(ENVIRON["F"]) }' < /dev/null)
        [ -n "$parent" ] || parent=.
        if [ -d "$parent" ]; then _ck_add OK out_dir "will be created: $out"
        else _ck_add WARN out_dir "parent directory does not exist: $out"; fi
    fi
    # FASTA36 path limit (check.py:220-225)
    local cwd; cwd=$(pwd -P)
    for what in out_dir TMPDIR; do
        if [ $what = out_dir ]; then p=$out
        else   # os.environ.get("TMPDIR") or "/tmp", in the caller's environment (bin/satpipe keeps it in SATPIPE_USER_TMPDIR)
            if [ -n "${SATPIPE_USER_TMPDIR_SET+x}" ]; then p=$SATPIPE_USER_TMPDIR; else p=${TMPDIR:-}; fi
            [ -n "$p" ] || p=/tmp
        fi
        p=$(F=$p CWD=$cwd "$GAWK" -f "$SATPIPE_LIB/normpath.awk" -f "$PYJSON" -f "$SATPIPE_LIB/check.awk" \
              -e 'BEGIN { printf "%s", py_abspath(ENVIRON["F"], ENVIRON["CWD"]) }' < /dev/null)
        len=$(S=$p "$GAWK" -f "$PYJSON" -f "$SATPIPE_LIB/check.awk" -e 'BEGIN { printf "%d", ck_clen(ck_utf8_fix(ENVIRON["S"])) }' < /dev/null)
        if [ "$len" -gt "$SATPIPE_MAX_BASE_PATH" ]; then
            _ck_add FAIL paths "$what is $len characters ($p); keep it at most $SATPIPE_MAX_BASE_PATH, or ssearch36 cannot open the files below it (path limit 119)"
        fi
    done
}

# _ck_tools [FLAT]   check_tools (check.py:231-247)
_ck_tools() {
    local flat=${1:-} bin="" name req stages note found path ver err lv
    if [ -n "$flat" ] && bin=$(_ck_val tools/bin_dir); then
        if [ -d "$bin" ]; then _ck_add OK tools "bin_dir: $bin"; else _ck_add FAIL tools "bin_dir not found: $bin"; fi
    fi
    while IFS=$'\037' read -r name req stages note found path ver err; do
        if [ "$req" = True ]; then lv=FAIL; else lv=WARN; fi
        stages=$(_ck_unesc "$stages"); note=$(_ck_unesc "$note")
        if [ "$found" != 1 ]; then
            _ck_add "$lv" tools "$name: not found (stages $stages)${note:+; $note}"
            continue
        fi
        path=$(_ck_unesc "$path"); ver=$(_ck_unesc "$ver"); err=$(_ck_unesc "$err")
        if [ -n "$err" ]; then _ck_add "$lv" tools "$name: $path: $err"
        else _ck_add OK tools "$name $ver: $path"; fi
    done < <(tools_probe_all "$flat" | tools_tsv_us /dev/stdin)
}

# _ck_missing_scripts   stages/__init__.py:88-95 over stages_registry.tsv; prints `total<TAB>missing, comma-joined`
_ck_scripts() {
    local id name wired opt hand needs scripts s total=0 miss=""
    local -a ss=()
    while IFS=$'\037' read -r id name wired opt hand needs scripts; do
        case $id in ''|'#'*) continue ;; esac
        [ -n "$scripts" ] || continue
        IFS=, read -r -a ss <<< "$scripts"
        for s in "${ss[@]}"; do
            total=$((total + 1))
            [ -f "$SCRIPTS/$s" ] || miss+=${miss:+, }$s
        done
    done < <(tools_tsv_us "$SATPIPE_LIB/stages_registry.tsv")
    printf '%s\t%s' "$total" "$miss"
}

check_run() {
    local file=$1 do_files=$2 do_tools=$3 abs e w sc
    CHECK_ITEMS=$(mktemp "${TMPDIR:-/tmp}/satpipe_check.XXXXXX")
    CHECK_FAILED=0
    CFG_FLAT=""
    abs=$(_config_abspath "$file")
    if [ -d "$abs" ] || [ ! -e "$abs" ] || [ ! -r "$abs" ]; then
        # open() raises OSError: run_check reports it and returns without the tool check (check.py:260-262)
        config_stream "$file" > "$CHECK_ITEMS.s" || true
        while IFS= read -r e; do _ck_add FAIL config "$(_ck_unesc "${e#*$'\t'}")"; done < "$CHECK_ITEMS.s"
        rm -f -- "$CHECK_ITEMS.s"
    elif ! config_load "$file"; then
        rm -f -- "$CFG_FLAT"; CFG_FLAT=""
        for e in "${CFG_ERRORS[@]}"; do _ck_add FAIL config "$e"; done
        if [ "$do_tools" = 1 ]; then _ck_tools ""; fi
    else
        _ck_add OK config "${CFG_PATH//\\//}: schema valid; run_name=$(cfg run_name), $(flat_get "$CFG_FLAT" genomes) genomes, $(flat_get "$CFG_FLAT" families) families"
        for w in "${CFG_WARNINGS[@]}"; do _ck_add WARN config "$w"; done
        sc=$(_ck_scripts)
        if [ -n "${sc#*$'\t'}" ]; then _ck_add FAIL package "stage scripts missing from the installed package: ${sc#*$'\t'}"
        else _ck_add OK package "${sc%%$'\t'*} stage script references present"; fi
        if [ "$do_files" = 1 ]; then _ck_files; fi
        if [ "$do_tools" = 1 ]; then _ck_tools "$CFG_FLAT"; fi
    fi
    if "$GAWK" -F '\t' '$1 == "FAIL" { f = 1; exit } END { exit !f }' "$CHECK_ITEMS"; then CHECK_FAILED=1; fi
    return 0
}

check_render() { P=$CHECK_ITEMS FORMAT=$1 _ck ck_render_main; }

check_cleanup() { rm -f -- "${CHECK_ITEMS:-}" "${CFG_FLAT:-}"; }

check_main() {
    local cfgfile=$1 no_files=$2 no_tools=$3 json=$4 rc
    check_run "$cfgfile" "$((1 - no_files))" "$((1 - no_tools))"
    if [ "$json" = 1 ]; then check_render json; else check_render text; fi
    rc=$CHECK_FAILED
    check_cleanup
    return "$rc"
}
