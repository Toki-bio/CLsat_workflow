# viewer.sh - `satpipe viewer RUN_DIR` (cli.py:59-67, viewer.py:48-271, viewer_lists.py). Sourced by bin/satpipe after
# common.sh; viewer_main reads the parsed options AP_out/AP_HAS_out, AP_viewer_html, AP_no_flank_scan.
#
# Steps in viewer.py's order: the config comes from manifest.json's `config` object (phase `manifest` of viewer.awk turns
# it into a TYPED-flat file: the legacy-manifest -> flat converter; every manifest, Python- or bash-written, carries the
# typed config, ints and floats as written, so no jq and no jq version matters); data steps run as viewer.awk phases;
# copies, gzip (-9 -n: member bytes differ from Python's GzipFile header, payload identical, PORT_SPEC 7) and the flank
# scans (expand_cmd + sh_run, logs/viewer_flankscan_<fam>.log) run here.
# Errors: Python's SystemExit texts verbatim (`viewer: ...`), StageError / FileNotFoundError as `satpipe viewer: <msg>`,
# other Python exceptions (tracebacks there) as `satpipe viewer: <msg>`; exit 1.

[ -n "${_SATPIPE_VIEWER_SH:-}" ] && return 0
_SATPIPE_VIEWER_SH=1

# shellcheck source=flat.sh
. "$SATPIPE_LIB/flat.sh"
# shellcheck source=config.sh
. "$SATPIPE_LIB/config.sh"
# shellcheck source=stages.sh
. "$SATPIPE_LIB/stages.sh"
# shellcheck source=expand_cmd.sh
. "$SATPIPE_LIB/expand_cmd.sh"

# the viewer's StageError text (cli.py:65), not the run's
stage_die() { die "satpipe viewer: $1"; }

_vw_awk() { VW_PHASE=$1 "$GAWK" -f "$PYJSON" -f "$SATPIPE_LIB/viewer_lists.awk" -f "$SATPIPE_LIB/viewer.awk"; }
_vw_repr() { S=$1 pyjson_gawk 'BEGIN { printf "%s", py_repr_str(ENVIRON["S"]) }'; }
_vw_nofile() { die "satpipe viewer: [Errno 2] No such file or directory: $(_vw_repr "$1")"; }
# shutil.copyfile SRC DST
_vw_copy() { [ -e "$1" ] || _vw_nofile "$1"; cp -- "$1" "$2" || die "satpipe viewer: cannot copy $1 to $2"; }
# _gzip SRC DST (viewer.py:28-30)
_vw_gzip() { gzip -9 -n -c < "$1" > "$2" || die "satpipe viewer: gzip failed on $1"; }

viewer_main() {
    local arg=$1 out html flank_scan=1 fam fams=() cand="" rc cdir t v sig old sigf fdir bed tmpf pre cf counts f k
    [ "${AP_no_flank_scan:-0}" = 1 ] && flank_scan=0
    RUN=$(F=$arg CWD=$(pwd -P) "$GAWK" -f "$SATPIPE_LIB/normpath.awk" -e 'BEGIN { printf "%s", py_abspath(ENVIRON["F"], ENVIRON["CWD"]) }')
    VW_TMP=$(mktemp -d "${TMPDIR:-/tmp}/satpipe_viewer.XXXXXX") || die "satpipe viewer: cannot create a temporary directory"
    trap 'rm -rf -- "$VW_TMP"' EXIT
    export VW_RUN=$RUN VW_TMP VW_CFG=$VW_TMP/config.flat.tsv VW_VERSION=$SATPIPE_VERSION

    # load_run_config (viewer.py:43-45)
    [ -e "$RUN/manifest.json" ] || _vw_nofile "$RUN/manifest.json"
    _vw_awk manifest > "$VW_CFG" || exit 1
    CFG_FLAT=$VW_CFG
    while IFS= read -r fam; do fams+=("$fam"); done < <(cfg_keys families)
    if [ "${AP_HAS_out:-0}" = 1 ] && [ -n "${AP_out:-}" ]; then out=$AP_out; else out=$RUN/viewer; fi
    mkdir -p -- "$out" || die "satpipe viewer: cannot create $out"
    export VW_OUT=$out
    _env_setup

    # 1. catalogue BED (+ candidates); 2. assignment + bits copies with candidate rows
    if [ -f "$RUN/candidates/candidates_summary.json" ]; then cand=$RUN/candidates/candidates_summary.json; fi
    export VW_CAND=$cand
    _vw_awk rows || exit 1
    _vw_copy "$RUN/catalogue/MASTER.final_assignment.tsv" "$out/MASTER.final_assignment.clean.tsv"
    _vw_copy "$RUN/catalogue/MASTER.locus_bits.tsv" "$out/MASTER.locus_bits.tsv"
    if [ -n "$cand" ]; then
        cat -- "$VW_TMP/cand_assign.tsv" >> "$out/MASTER.final_assignment.clean.tsv"
        cat -- "$VW_TMP/cand_bits.tsv" >> "$out/MASTER.locus_bits.tsv"
    fi

    # 4. per-family bits tracks
    for fam in "${fams[@]}"; do
        t=$out/MASTER_$fam.bits.tsv
        _vw_copy "$RUN/catalogue/MASTER_$fam.bits.tsv" "$t"
        if [ -n "$cand" ] && [ -f "$RUN/candidates/CAND_$fam.bits.tsv" ]; then cat -- "$RUN/candidates/CAND_$fam.bits.tsv" >> "$t"; fi
        _vw_gzip "$t" "$t.gz"
        rm -f -- "$t"
    done

    # 5. flank windows (background track)
    if [ $flank_scan = 1 ]; then
        fdir=$RUN/viewer_flankscan; bed=$fdir/flanks.bed; sigf=$fdir/flanks.bed.sha256
        mkdir -p -- "$fdir" || die "satpipe viewer: cannot create $fdir"
        export VW_FDIR=$fdir
        v=$(_vw_awk flankbed) || exit 1
        sig=$({ cat -- "$bed"; printf '%s' "$v"; } | sha256sum) || die "satpipe viewer: sha256sum failed"
        sig=${sig%% *}
        old=""
        if [ -f "$sigf" ]; then
            old=$("$GAWK" 'BEGIN { RS = "^$" } { gsub(/^[ \t\n\r\v\f\034-\037]+|[ \t\n\r\v\f\034-\037]+$/, ""); printf "%s", $0 }' "$sigf")
        fi
        if [ ! -f "$sigf" ] || [ "$old" != "$sig" ]; then
            # stale scans are dropped only for other flank windows (GLM audit #56, viewer.py:122-132)
            for fam in "${fams[@]}"; do
                for f in "$fdir"/* "$fdir"/.[!.]* "$fdir"/..?*; do
                    [ -e "$f" ] || [ -L "$f" ] || continue
                    k=${f##*/}
                    [[ $k == "flank_$fam."* ]] && rm -f -- "$f"
                done
            done
            printf '%s\n' "$sig" > "$sigf"
        fi
        tmpf=$out/MASTER_flank.bits.tsv
        printf 'scaffold\tpos\tbits\tclsat\tspecies\n' > "$tmpf"
        for fam in "${fams[@]}"; do
            pre=$fdir/flank_$fam
            if [ ! -f "$pre.bits.tsv" ]; then
                say "  flank scan $fam"
                for k in families/$fam/mono families/$fam/consensus expansion/step_mono expansion/win_mono expansion/query_repeat \
                         expansion/min_aln_frac expansion/min_bits expansion/max_steps expansion/probe_mono expansion/probe_min_step \
                         scan/win_mono scan/step_mono scan/min_bits expansion/jobs; do
                    cfg_type "$k" > /dev/null 2>&1 || die "satpipe viewer: KeyError: $(_vw_repr "${k##*/}")"
                done
                expand_cmd "$bed" "$fam" "$pre" 1
                sh_run viewer_flankscan "" "log_name=viewer_flankscan_$fam" -- "${EXPAND_CMD[@]}"
            fi
            VW_PRE=$pre VW_FAM=$fam _vw_awk flankfam >> "$tmpf" || exit 1
        done
        cf=$RUN/candidates/CAND_FLANK.bits.tsv
        if [ -n "$cand" ] && [ "$(cat -- "$VW_TMP/cand_summary_nonempty")" = 1 ] && [ ! -f "$cf" ]; then
            die "viewer: candidates exist but $cf is missing; rerun stage 20_candidates"
        fi
        if [ -n "$cand" ] && [ -f "$cf" ]; then _vw_awk flankcand >> "$tmpf" || exit 1; fi
        _vw_gzip "$tmpf" "$tmpf.gz"
        rm -f -- "$tmpf"
    fi

    # 6. JSON products
    _vw_awk flanksjson || exit 1
    _vw_copy "$RUN/ortho/flank_orthologs.json" "$out/flank_orthologs.json"
    # stage 19 decided the pairs (orthology.clear_margin > 0): the viewer lists exactly those (viewer.py, flank-single-copy)
    if CM=$(cfg orthology/clear_margin) "$GAWK" 'BEGIN { exit !(ENVIRON["CM"] + 0 > 0) }'; then
        _vw_copy "$RUN/ortho/ortho_pairs.tsv" "$out/ortho_pairs.tsv"
    fi
    for f in subfamily.json terminal_monomers.json; do
        if [ -f "$RUN/regions/$f" ]; then _vw_copy "$RUN/regions/$f" "$out/$f"; fi
    done
    if [ -f "$out/subfamily.json" ] && [ -n "$cand" ]; then _vw_awk subfam || exit 1; fi

    # 7. scaffold lengths, consensuses
    _vw_awk scaflen || exit 1
    for fam in "${fams[@]}"; do _vw_copy "$(cfg "families/$fam/consensus")" "$out/$fam.cons.fa"; done

    # 8. viewer_config.json (subfamily.consensuses must be a directory, os.path.isdir)
    export VW_SUBDIR_OK=0
    case $(cfg_type subfamily/consensuses 2>/dev/null) in
        str|path) cdir=$(cfg subfamily/consensuses); [ -n "$cdir" ] && [ -d "$cdir" ] && VW_SUBDIR_OK=1 ;;
    esac
    counts=$(_vw_awk config) || exit 1

    # 9. the page
    if [ "${AP_HAS_viewer_html:-0}" = 1 ] && [ -n "${AP_viewer_html:-}" ]; then html=$AP_viewer_html; else html=$SATPIPE_ROOT/satpipe/viewer_html/index.html; fi
    if [ -f "$html" ]; then _vw_copy "$html" "$out/index.html"
    else say "  WARNING: viewer page not found ($html); data files and viewer_config.json written without index.html"; fi
    : > "$out/.nojekyll"
    set -- $counts
    say "viewer: $1 catalogue loci, $2 candidates, $3 anchorable, $4 flagged -> $out"
    return 0
}
