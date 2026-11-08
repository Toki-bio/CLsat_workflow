# anchor.sh - stage 21_anchor (run.py:1039-1101) and `satpipe anchor` (cli.py:91-125) on top of anchor.awk.
# Sourced after common.sh and flat.sh (PORT_SPEC section 1). Also runnable standalone on a TYPED-flat config:
#   bash lib/satpipe/anchor.sh CONFIG_FLAT [--check] [-o OUT] FASTA...
#
# Functions
#   anchor_fasta_files PATH            fasta_files() (anchor.py:198-204): PATH, or the FASTA files of a directory, one per line
#   anchor_cli FLAT ARGS...            `satpipe anchor` after run_check: FLAT = the checked config as TYPED flat
#   plan_21_anchor / st_21_anchor      stage contract of PORT_SPEC section 1
#
# What the stage functions take from the driver (WP3/WP5; tests/sh/test_anchor.sh carries a minimal stub of exactly these):
#   cfg PATH (py_str text of a leaf), cfg_type PATH (type column; status 1 when absent), cfg_keys PATH (child keys in
#   document order), cfg_subtree SRC DST (appends the re-rooted rows to $PLAN_PARAMS), emit_inputs PATH... (appends to
#   $PLAN_INPUTS), emit_outputs PATH... (appends to $STAGE_OUTPUTS), say MSG, $RUN; optional: manifest_set KEY JSON
#   (manifest[KEY] = JSON; without it the JSON is kept in $ANCHOR_MANIFEST_JSON), stage_die MSG (default: die "satpipe run: MSG").

[ -n "${_SATPIPE_ANCHOR_SH:-}" ] && return 0
_SATPIPE_ANCHOR_SH=1

# one job-file row: fields TAB-joined, each escaped like flat_esc (\ TAB LF CR)
_anc_row() {
    local out="" f sep=""
    for f in "$@"; do
        f=${f//\\/\\\\}; f=${f//$'\t'/\\t}; f=${f//$'\n'/\\n}; f=${f//$'\r'/\\r}
        out+=$sep$f; sep=$'\t'
    done
    printf '%s\n' "$out"
}

anchor_fasta_files() {
    local p=$1 f n l
    if [ -d "$p" ]; then
        local -a all=()
        local old; old=$(shopt -p nullglob dotglob || true)
        shopt -s nullglob dotglob
        all=("$p"/*)
        eval "$old"
        for f in "${all[@]}"; do
            n=${f##*/}; l=${n,,}
            case $l in *.anchored.*) continue ;; esac
            case $l in *.fa|*.fas|*.fasta|*.fna) ;; *) continue ;; esac
            case $p in */) printf '%s\n' "$p$n" ;; *) printf '%s\n' "$p/$n" ;; esac   # os.path.join
        done
    else
        printf '%s\n' "$p"
    fi
}

# _anc_settings GET TYPE  - the anchor block as cli.py:98-101 / run.py:1044-1047 read it ('or' defaults included)
_anc_settings() {
    local get=$1 typ=$2 v
    if [ "$($typ anchor/family 2>/dev/null)" = str ]; then v=$($get anchor/family); else v=""; fi; _anc_row set afam "$v"
    if [ "$($typ anchor/motif 2>/dev/null)" = str ]; then v=$($get anchor/motif); else v=""; fi; _anc_row set motif "$v"
    case $($typ anchor/max_motif_mismatch 2>/dev/null) in int) v=$($get anchor/max_motif_mismatch) ;; *) v=1 ;; esac; _anc_row set mm "$v"
    case $($typ anchor/min_identity 2>/dev/null) in float|int) v=$($get anchor/min_identity) ;; *) v=0.6 ;; esac; _anc_row set mi "$v"
    if [ "$($typ anchor/independent 2>/dev/null)" = bool ]; then v=$($get anchor/independent); else v=False; fi; _anc_row set indep "$v"
}

# _anc_families GET KEYS
_anc_families() {
    local get=$1 keys=$2 f
    while IFS= read -r f; do _anc_row fam "$f" "$($get "families/$f/consensus")"; done < <($keys families)
}

_anc_unesc() { "$GAWK" -f "$PYJSON" -e 'BEGIN { printf "%s", flat_unesc(ENVIRON["S"]) }' </dev/null; }
_anc_gawk() { ANC_JOB=$1 "$GAWK" -f "$PYJSON" -f "$SATPIPE_LIB/anchor.awk" -e 'BEGIN { anchor_main() }'; }

# ------------------------------------------------------------------ satpipe anchor
_anc_flat_get()  { flat_get "$ANC_FLAT" "$1"; }
_anc_flat_type() { flat_type "$ANC_FLAT" "$1"; }
_anc_flat_keys() { flat_keys "$ANC_FLAT" "$1"; }

anchor_cli() {
    local ANC_FLAT=$1; shift
    local check=0 out="" have_out=0 a p fp job rc
    local -a fasta=()
    while [ $# -gt 0 ]; do
        case $1 in
            --check) check=1 ;;
            -o|--out) [ $# -ge 2 ] || { printf 'usage: satpipe anchor [-h] [--check] [-o OUT] config fasta [fasta ...]\nsatpipe anchor: error: argument -o/--out: expected one argument\n' >&2; return 2; }
                      out=$2; have_out=1; shift ;;
            --out=*) out=${1#--out=}; have_out=1 ;;
            -o?*) out=${1#-o}; have_out=1 ;;
            --) shift; fasta+=("$@"); break ;;
            -?*) printf 'usage: satpipe anchor [-h] [--check] [-o OUT] config fasta [fasta ...]\nsatpipe anchor: error: unrecognized arguments: %s\n' "$1" >&2; return 2 ;;
            *) fasta+=("$1") ;;
        esac
        shift
    done
    [ ${#fasta[@]} -gt 0 ] || { printf 'usage: satpipe anchor [-h] [--check] [-o OUT] config fasta [fasta ...]\nsatpipe anchor: error: the following arguments are required: fasta\n' >&2; return 2; }
    job=$(mktemp "${TMPDIR:-/tmp}/satpipe_anchor_job.XXXXXX")
    {
        _anc_row set mode cli; _anc_row set check "$check"
        _anc_settings _anc_flat_get _anc_flat_type
        _anc_families _anc_flat_get _anc_flat_keys
        for p in "${fasta[@]}"; do
            while IFS= read -r fp; do
                # cli.py:119: -o only for a single non-directory input, otherwise <file>.anchored.fa
                if [ -n "$out" ] && [ ${#fasta[@]} -eq 1 ] && [ ! -d "$p" ]; then _anc_row input "$fp" "$out"; else _anc_row input "$fp" "$fp.anchored.fa"; fi
            done < <(anchor_fasta_files "$p")
        done
    } > "$job"
    rc=0; _anc_gawk "$job" || rc=$?
    rm -f -- "$job"
    return "$rc"
}

# ------------------------------------------------------------------ stage 21_anchor
_anc_cfg_get()  { cfg "$1"; }
_anc_cfg_type() { cfg_type "$1"; }
_anc_cfg_keys() { cfg_keys "$1" 2>/dev/null || true; }

# the configured subfamily.consensuses, or nothing when null/absent
_anc_cons() { case $(cfg_type subfamily/consensuses 2>/dev/null) in path|str) cfg subfamily/consensuses ;; esac; }

plan_21_anchor() {                                                     # run.py:1099-1101
    local f c
    while IFS= read -r f; do emit_inputs "$(cfg "families/$f/consensus")"; done < <(_anc_cfg_keys families)
    c=$(_anc_cons)
    if [ -n "$c" ] && [ -e "$c" ]; then emit_inputs "$c"; fi
    if cfg_type anchor >/dev/null 2>&1; then cfg_subtree anchor anchor; else flat_emit anchor null None >> "$PLAN_PARAMS"; fi
    if cfg_type subfamily/members >/dev/null 2>&1; then cfg_subtree subfamily/members members; else flat_emit members null None >> "$PLAN_PARAMS"; fi
}

st_21_anchor() {
    local out=$RUN/anchored job c fp fam i res rc line cls msg man sayt T=$'\t'
    mkdir -p -- "$out"                                                 # ctx.d('anchored'), run.py:1050
    job=$(mktemp "${TMPDIR:-/tmp}/satpipe_anchor_job.XXXXXX")
    {
        _anc_row set mode stage; _anc_row set outdir "$out"
        _anc_settings _anc_cfg_get _anc_cfg_type
        _anc_families _anc_cfg_get _anc_cfg_keys
        c=$(_anc_cons)
        if [ -n "$c" ]; then
            if [ -e "$c" ]; then
                _anc_row cons "$c" 1
                while IFS= read -r fp; do _anc_row subfile "$fp"; done < <(anchor_fasta_files "$c")
            else _anc_row cons "$c" 0; fi
        fi
        while IFS= read -r fam; do
            while IFS= read -r i; do _anc_row member "$fam" "$(cfg "subfamily/members/$fam/$i")"; done < <(_anc_cfg_keys "subfamily/members/$fam")
        done < <(_anc_cfg_keys subfamily/members)
    } > "$job"
    rc=0; res=$(_anc_gawk "$job") || rc=$?
    rm -f -- "$job"
    man=""; sayt=""; cls=""; msg=""
    while IFS= read -r line; do
        case $line in
            MANIFEST$T*) man=${line#MANIFEST$T} ;;
            SAY$T*) sayt=${line#SAY$T} ;;
            ERROR$T*) [ -n "$cls" ] || { line=${line#ERROR$T}; cls=${line%%$T*}; msg=${line#*$T}; } ;;
        esac
    done <<< "$res"
    if [ "$rc" != 0 ]; then
        # StageError texts as run.py raises them; other classes (AnchorError, ValueError, ...) are a traceback in Python
        if [ -n "$cls" ]; then msg=$(S=$msg _anc_unesc); else msg="21_anchor: anchor.awk failed (exit $rc)"; fi
        if declare -F stage_die >/dev/null; then stage_die "$msg"; else die "satpipe run: $msg"; fi
        return 1
    fi
    sayt=$(S=$sayt _anc_unesc)
    if declare -F manifest_set >/dev/null; then manifest_set anchor "$man"; else ANCHOR_MANIFEST_JSON=$man; fi
    say "$sayt"
    emit_outputs "$out/families.anchored.fa" "$out/anchor_report.tsv" "$out/anchor_frame.fa"   # run.py:1096 (not subfamilies.anchored.fa)
}

# ------------------------------------------------------------------ standalone
if [ "${BASH_SOURCE[0]}" = "$0" ]; then
    set -euo pipefail
    source "$(dirname -- "$(readlink -f -- "$0")")/common.sh"
    source "$SATPIPE_LIB/flat.sh"
    require_gawk
    [ $# -ge 1 ] || die "usage: anchor.sh CONFIG_FLAT [--check] [-o OUT] FASTA..." 2
    anchor_cli "$@"
fi
