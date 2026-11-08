# subx.sh - subx_dir (run.py:108-132) and _py_need_file, shared by the WP8 stage files st_16..st_20, which source it
# after common.sh, flat.sh and config.sh (cfg, cfg_type). Test: tests/sh/test_subx.sh.
#
#   subx_dir        SUBX = the directory of DOUBLED per-sub-subfamily FASTA files, or "" (Python None):
#                   subfamily.consensuses absent or empty -> "" (run.py:113-114); a directory (os.path.isdir follows
#                   symlinks) -> the path as given (:115-116); not a regular file -> "" (:117-118); else the multi-record
#                   FASTA is written out as $RUN/tables/subx/<record>.fa = ">name\nseqseq\n" (:119-130) and SUBX is that dir.
#                   Sets a variable instead of printing, so it is never run in a subshell (it writes files).
# FASTA reading as Python's text-mode `for line in open(cons)` (:122-128): universal newlines (\n, \r\n, \r end a line);
# a header's name = first whitespace-separated word after '>'; sequence lines are str.strip()ed and joined; lines before
# the first header are dropped; a record name seen twice: the later record's file wins (files are written after reading,
# in record order). ASCII whitespace set of str.split/strip: space \t \n \v \f \r \x1c-\x1f.

[ -n "${_SATPIPE_SUBX_SH:-}" ] && return 0
_SATPIPE_SUBX_SH=1

# _py_need_file FILE   a file a Python stage body open()s for reading: FileNotFoundError's text when it is absent
_py_need_file() { [ -f "$1" ] || stage_die "[Errno 2] No such file or directory: '$1'"; }

subx_dir() {
    local cons out
    SUBX=""
    case $(cfg_type subfamily/consensuses 2>/dev/null) in path|str) cons=$(cfg subfamily/consensuses) ;; *) return 0 ;; esac
    [ -n "$cons" ] || return 0
    if [ -d "$cons" ]; then SUBX=$cons; return 0; fi
    [ -f "$cons" ] || return 0
    out=$RUN/tables/subx
    mkdir -p -- "$out" || stage_die "[Errno 13] Permission denied: '$out'"
    OUT=$out "$GAWK" '
        BEGIN { RS = "\r\n|\r|\n"; WS = "[ \t\n\v\f\r\034-\037]" }
        /^>/ {
            if (name != "") { recs[++n] = name; seqs[n] = seq }
            h = substr($0, 2); sub("^" WS "+", "", h)
            if (h == "") { print "list index out of range" > "/dev/stderr"; bad = 1; exit 1 }   # line[1:].split()[0]
            match(h, WS); name = RSTART ? substr(h, 1, RSTART - 1) : h; seq = ""
            next
        }
        { s = $0; gsub("^" WS "+|" WS "+$", "", s); seq = seq s }
        END {
            if (bad) exit 1
            if (name != "") { recs[++n] = name; seqs[n] = seq }
            for (i = 1; i <= n; i++) {
                f = ENVIRON["OUT"] "/" recs[i] ".fa"
                printf ">%s\n%s\n", recs[i], seqs[i] seqs[i] > f
                close(f)
            }
        }' "$cons" || stage_die "subx_dir: cannot write $out from $cons"
    SUBX=$out
    return 0
}

# _p8_confident_delta   str(sub.get('confident_delta', 10)) (run.py:982, :1061)
_p8_confident_delta() {
    case $(cfg_type subfamily/confident_delta 2>/dev/null) in int|float) cfg subfamily/confident_delta ;; *) printf '10' ;; esac
}
