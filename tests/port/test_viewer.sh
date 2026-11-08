#!/usr/bin/env bash
# test_viewer.sh REPO SRC_RUN WORK [seed] [errors]
# `satpipe viewer` Python vs bash on a finished run directory (WP10 acceptance). SRC_RUN is only read: each driver gets
# WORK/{py,sh}/run with symlinks to SRC_RUN's catalogue candidates regions ortho tables manifest.json (the viewer writes
# viewer/, viewer_flankscan/ and logs/ into the run directory). `seed`: copy SRC_RUN/viewer_flankscan in first (real
# copy), so earlier scans are reused when their signature matches. `errors`: also the error cases (toy-sized runs only).
# Variants: A default (flank scan, <run>/viewer); R the same again (scan reuse); N --no-flank-scan --out --viewer-html FILE;
# W --viewer-html missing (WARNING line). Compared: viewer dirs (diff_viewer.sh), stdout, stderr, exit code,
# viewer_flankscan/ files, logs/viewer_flankscan_*.log (timestamps masked).
# `addgenome`: runs from before the genome column (col 5 of catalogue/MASTER.loci.segmented.bed, which both drivers require)
# get a copy of that file with col 5 = the genome whose .fai holds the contig (contig names must be unique across genomes);
# the other catalogue files stay symlinks. Same fixture for both drivers.
# Env: PY (default python3) runs `python -m satpipe.cli` with PYTHONPATH=REPO.
set -uo pipefail
export LC_ALL=C
REPO=$(cd "${1:?repo}" && pwd); SRC=${2:?src_run}; W=${3:?work}; shift 3
SEED=0; ERR=0; ADDG=0
for a in "$@"; do case $a in seed) SEED=1 ;; errors) ERR=1 ;; addgenome) ADDG=1 ;; esac; done
PY=${PY:-python3}
mkdir -p "$W"; W=$(cd "$W" && pwd)
fail=0
note() { printf '%s\n' "$*"; }

# mkrun DIR SRC [PART...]   a fresh run dir (DIR must not exist): symlinks to SRC's catalogue candidates regions ortho tables
# manifest.json, except each PART, which becomes an empty real directory (an error case copies files in; nothing is ever
# written through a symlink into SRC); PART `no<name>` leaves <name> out. With `addgenome` the catalogue is a real
# directory of per-file symlinks plus the rebuilt segmented BED (unless it is a PART).
mkrun() {
    local d=$1 s=$2 x p own; shift 2
    [ ! -e "$d" ] || { echo "mkrun: $d exists"; exit 1; }
    mkdir -p "$d"
    for x in catalogue candidates regions ortho tables manifest.json; do
        own=0
        for p in "$@"; do [ "$p" = "$x" ] && own=1; [ "$p" = "no$x" ] && own=2; done
        if [ $own = 0 ] && [ $x = catalogue ] && [ $ADDG = 1 ]; then own=3; fi
        case $own in
            0) [ -e "$s/$x" ] && ln -s "$s/$x" "$d/$x" ;;
            1) mkdir "$d/$x" ;;
            3) mkdir "$d/catalogue"
               for x in "$s"/catalogue/*; do [ "${x##*/}" = MASTER.loci.segmented.bed ] || ln -s "$x" "$d/catalogue/${x##*/}"; done
               "$PY" - "$s/manifest.json" "$s/catalogue/MASTER.loci.segmented.bed" > "$d/catalogue/MASTER.loci.segmented.bed" <<'PYX' || { echo "addgenome failed"; exit 1; }
import json, sys
m = json.load(open(sys.argv[1]))
own = {}
for code, g in m['config']['genomes'].items():
    for line in open(g + '.fai'):
        c = line.split(chr(9))[0]
        if c in own and own[c] != code:
            sys.exit('contig %s in two genomes' % c)
        own[c] = code
for line in open(sys.argv[2]):
    p = line.rstrip(chr(10)).split(chr(9))
    if len(p) == 4:
        p.append(own[p[0]])
    print(chr(9).join(p))
PYX
               ;;
        esac
    done
    if [ $SEED = 1 ] && [ -d "$s/viewer_flankscan" ]; then cp -r "$s/viewer_flankscan" "$d/viewer_flankscan"; fi
    return 0
}

# drive VARIANT DRIVER ARGS...   (@D@ in ARGS = the driver's work dir)
drive() {
    local v=$1 drv=$2 d=$W/$2; shift 2
    local -a args=() ; local a
    for a in "$@"; do args+=("${a//@D@/$d}"); done
    if [ "$drv" = py ]; then (cd "$d" && PYTHONPATH=$REPO "$PY" -m satpipe.cli viewer "${args[@]}") > "$W/$v.$drv.out" 2> "$W/$v.$drv.err"
    else (cd "$d" && "$REPO/bin/satpipe" viewer "${args[@]}") > "$W/$v.$drv.out" 2> "$W/$v.$drv.err"; fi
    echo $? > "$W/$v.$drv.rc"
}

norm() { sed -E "s#$W/(py|sh)#WORK#g; s/[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}/TS/g"; }

cmp_text() {   # cmp_text LABEL FILE_PY FILE_SH
    if diff <(norm < "$2") <(norm < "$3") > "$W/diff.tmp"; then return 0; fi
    note "  DIFF $1:"; sed 's/^/    /' "$W/diff.tmp" | head -20; fail=1
}

check() {      # check VARIANT VIEWER_SUBDIR(relative to the driver dir, or '')
    local v=$1 vd=$2 r f
    cmp_text "$v stdout" "$W/$v.py.out" "$W/$v.sh.out"
    cmp_text "$v stderr" "$W/$v.py.err" "$W/$v.sh.err"
    [ "$(cat "$W/$v.py.rc")" = "$(cat "$W/$v.sh.rc")" ] || { note "  DIFF $v exit: py $(cat "$W/$v.py.rc") sh $(cat "$W/$v.sh.rc")"; fail=1; }
    if [ -n "$vd" ]; then
        r=$(bash "$REPO/tests/port/diff_viewer.sh" "$W/py/$vd" "$W/sh/$vd") || { fail=1; }
        printf '%s\n' "$r" | sed "s/^/  $v viewer: /"
    fi
    if [ -d "$W/py/run/viewer_flankscan" ] || [ -d "$W/sh/run/viewer_flankscan" ]; then
        r=$(bash "$REPO/tests/port/diff_viewer.sh" "$W/py/run/viewer_flankscan" "$W/sh/run/viewer_flankscan") || { fail=1; }
        printf '%s\n' "$r" | sed "s/^/  $v viewer_flankscan: /"
    fi
    for f in "$W"/py/run/logs/viewer_flankscan_*.log; do
        [ -f "$f" ] || continue
        cmp_text "$v ${f##*/}" "$f" "$W/sh/run/logs/${f##*/}"
    done
    note "  $v: exit $(cat "$W/$v.sh.rc"); stdout $(wc -l < "$W/$v.sh.out") lines; last: $(tail -n 1 "$W/$v.sh.out" | norm)"
}

for d in py sh; do mkrun "$W/$d/run" "$SRC"; done
printf '<!doctype html><title>custom page</title>\n' > "$W/custom.html"

note "== $SRC"
for d in py sh; do drive A $d @D@/run; done; check A run/viewer
[ "$(cat "$W/A.py.rc")" = 0 ] || { note "  A FAILED in both drivers (no viewer built): $(tail -n 1 "$W/A.py.err")"; fail=1; }
for d in py sh; do drive R $d @D@/run; done; check R run/viewer
for d in py sh; do drive N $d @D@/run --no-flank-scan --out @D@/noflank --viewer-html "$W/custom.html"; done; check N noflank
for d in py sh; do drive W $d @D@/run --no-flank-scan --out @D@/nohtml --viewer-html "$W/missing.html"; done; check W nohtml

if [ $ERR = 1 ]; then
    note "== error cases"
    for d in py sh; do drive E1 $d @D@/no_such_run; done; check E1 ""
    for d in py sh; do                                  # an old segmented catalogue without the genome column
        mkrun "$W/$d/e2" "$SRC" catalogue
        cp "$SRC"/catalogue/* "$W/$d/e2/catalogue/"
        cut -f 1-4 "$SRC/catalogue/MASTER.loci.segmented.bed" > "$W/$d/e2/catalogue/MASTER.loci.segmented.bed"
        drive E2 $d @D@/e2 --no-flank-scan
    done; check E2 ""
    for d in py sh; do                                  # candidates without their flank scan
        mkrun "$W/$d/e3" "$SRC" candidates
        [ -d "$SRC/candidates" ] && cp "$SRC"/candidates/* "$W/$d/e3/candidates/"
        rm -f "$W/$d/e3/candidates/CAND_FLANK.bits.tsv"
        drive E3 $d @D@/e3
    done; check E3 ""
    for d in py sh; do                                  # no candidates at all, no subfamily.json
        mkrun "$W/$d/e4" "$SRC" regions nocandidates
        cp "$SRC"/regions/* "$W/$d/e4/regions/" 2>/dev/null; rm -f "$W/$d/e4/regions/subfamily.json"
        drive E4 $d @D@/e4
    done; check E4 e4/viewer
fi
note "RESULT $SRC: $([ $fail = 0 ] && echo IDENTICAL || echo DIFFERENT)"
exit $fail
