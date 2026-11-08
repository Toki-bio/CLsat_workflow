# e2e_lib.sh - shared functions of the WP11 end-to-end harnesses (tests/port/e2e_toy.sh, tests/port/e2e_molly.sh).
# Sourced, not run. Needs: ROOT (repo), W (work dir), PY (python with PyYAML), GAWK, jq, cmp, diff.
# Model: tests/port/test_stages_wp6.sh (same masking rules: ISO timestamps; Python's sys.executable -> python3; parallel
# 07_mafft line order compared as a set), extended to whole runs: the script stamps `[YYYY-MM-DD HH:MM:SS]` of
# arbitration logs are masked too (PORT_LOG WP7), and log files are reported in their own category.

NPASS=${NPASS:-0}; NFAIL=${NFAIL:-0}; NLOG=${NLOG:-0}
ok()   { NPASS=$((NPASS + 1)); printf 'ok   %s\n' "$1"; }
fail() { NFAIL=$((NFAIL + 1)); printf 'FAIL %s\n' "$1"; [ -n "${2:-}" ] && [ -f "$2" ] && head -n 30 "$2" | cut -c1-300; return 0; }
logd() { NLOG=$((NLOG + 1)); printf 'LOG  %s\n' "$1"; [ -n "${2:-}" ] && [ -f "$2" ] && head -n 12 "$2" | cut -c1-300; return 0; }

PYX=$("$PY" -c 'import sys; print(sys.executable)')
norm() {   # timestamps, Python's interpreter path (Step 3b argv[0]) -> python3
    "$GAWK" -v PYX="$PYX" '{
        gsub(/[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9][T ][0-9][0-9]:[0-9][0-9]:[0-9][0-9](\.[0-9]+)?/, "TS")
        if (PYX != "") { while ((i = index($0, PYX)) > 0) $0 = substr($0, 1, i - 1) "python3" substr($0, i + length(PYX)) }
        print }'
}

# drive SIDE OUTPREFIX satpipe-args...   (cwd $W; SATPIPE_TRACE = OUTPREFIX.trace)
drive() {
    local side=$1 o=$2 rc=0; shift 2
    : > "$o.trace"
    if [ "$side" = py ]; then
        ( cd "$W" && PYTHONPATH=$ROOT SATPIPE_TRACE="$o.trace" "$PY" -m satpipe.cli "$@" ) > "$o.stdout" 2> "$o.stderr" || rc=$?
    else
        ( cd "$W" && SATPIPE_TRACE="$o.trace" "$ROOT/bin/satpipe" "$@" ) > "$o.stdout" 2> "$o.stderr" || rc=$?
    fi
    echo "$rc" > "$o.rc"
}

# cmp_text LABEL A B OUT : in order after masking, else as a sorted set (parallel stages)
cmp_text() {
    local l=$1 a=$2 b=$3 o=$4
    if diff <(norm < "$a") <(norm < "$b") > "$o"; then ok "$l ($(grep -c '' "$b") lines, in order)"
    elif diff <(norm < "$a" | LC_ALL=C sort) <(norm < "$b" | LC_ALL=C sort) > "$o.sorted"; then ok "$l ($(grep -c '' "$b") lines, same set, parallel order differs)"
    else fail "$l" "$o.sorted"; fi
}

# compare_runs LABEL RUN_A RUN_B OUTPREFIX
compare_runs() {
    local v=$1 a=$2 b=$3 o=$4 f bad=0 n=0 nd=0 nl=0 nm=0 ng=0 nlogd=0
    ( cd "$a" && find . \( -type f -o -type l \) | LC_ALL=C sort ) > "$o.files.a"
    ( cd "$b" && find . \( -type f -o -type l \) | LC_ALL=C sort ) > "$o.files.b"
    if diff "$o.files.a" "$o.files.b" > "$o.diff.files"; then ok "$v files: same list ($(grep -c '' "$o.files.a") files)"; else fail "$v files: lists differ" "$o.diff.files"; fi
    : > "$o.diff.content"; : > "$o.diff.logs"
    while IFS= read -r f; do
        f=${f#./}
        [ -e "$b/$f" ] || [ -L "$b/$f" ] || continue
        if [ -L "$a/$f" ]; then
            n=$((n + 1)); [ "$(readlink "$a/$f")" = "$(readlink "$b/$f")" ] || { echo "  link differs: $f" >> "$o.diff.content"; bad=1; }
            continue
        fi
        case $f in
            logs/satpipe.log)
                nl=$((nl + 1))
                diff <(norm < "$a/$f") <(norm < "$b/$f") > /dev/null \
                || diff <(norm < "$a/$f" | LC_ALL=C sort) <(norm < "$b/$f" | LC_ALL=C sort) >> "$o.diff.content" \
                || { echo "  satpipe.log differs (as a set): $f" >> "$o.diff.content"; bad=1; } ;;
            logs/*|*.log)
                nl=$((nl + 1))
                cmp -s "$a/$f" "$b/$f" || diff <(norm < "$a/$f") <(norm < "$b/$f") > /dev/null \
                || { nlogd=$((nlogd + 1)); { echo "=== $f"; diff <(norm < "$a/$f") <(norm < "$b/$f") | head -n 8; } >> "$o.diff.logs"; } ;;
            */DONE) nd=$((nd + 1)); diff <(norm < "$a/$f") <(norm < "$b/$f") >> "$o.diff.content" || { echo "  DONE differs: $f" >> "$o.diff.content"; bad=1; } ;;
            manifest.json) nm=1; diff <(jq -S 'del(.runs)' "$a/$f" | norm) <(jq -S 'del(.runs)' "$b/$f" | norm) >> "$o.diff.content" || { echo "  manifest differs" >> "$o.diff.content"; bad=1; } ;;
            *.gz) ng=$((ng + 1)); cmp -s "$a/$f" "$b/$f" || cmp -s <(gzip -dc < "$a/$f") <(gzip -dc < "$b/$f") || { echo "  gz payload differs: $f" >> "$o.diff.content"; bad=1; } ;;
            *) n=$((n + 1)); cmp -s "$a/$f" "$b/$f" || { echo "  cmp differs: $f" >> "$o.diff.content"; bad=1; } ;;
        esac
    done < "$o.files.a"
    if [ $bad = 0 ]; then ok "$v content: $n data files cmp-identical, $ng gz payloads equal, $nd DONE equal (finished masked), manifest jq -S equal ($nm), satpipe.log equal"
    else fail "$v content" "$o.diff.content"; fi
    if [ $nlogd = 0 ]; then ok "$v logs: $nl log files equal after masking"
    else logd "$v logs: $nlogd of $nl log files differ after masking (see $o.diff.logs)" "$o.diff.logs"; fi
}

# compare_drive LABEL PREFIX_A PREFIX_B : rc, stdout, stderr, trace of two drive() calls
compare_drive() {
    local v=$1 a=$2 b=$3
    if [ "$(cat "$a.rc")" = "$(cat "$b.rc")" ]; then ok "$v rc=$(cat "$b.rc") both"; else fail "$v rc a=$(cat "$a.rc") b=$(cat "$b.rc")" "$b.stderr"; fi
    cmp_text "$v stdout" "$a.stdout" "$b.stdout" "$b.diff.stdout"
    cmp_text "$v stderr" "$a.stderr" "$b.stderr" "$b.diff.stderr"
    cmp_text "$v trace PLAN=$(grep -c '^PLAN' "$b.trace") CMD=$(grep -c '^CMD' "$b.trace") OUT=$(grep -c '^OUT' "$b.trace")" "$a.trace" "$b.trace" "$b.diff.trace"
}

# cross_rerun LABEL SIDE RUNCOPY_SRC RUN_PATH OUTPREFIX satpipe-run-args...
#   SIDE's driver reruns on a copy of the other driver's finished run dir: must run nothing (no CMD/OUT records), exit 0,
#   leave every DONE byte-identical, and print 'up to date, skipped' once per DONE.
cross_rerun() {
    local v=$1 side=$2 src=$3 run=$4 o=$5 k nd ch; shift 5
    rm -rf -- "$run"; cp -a "$src" "$run"
    ( cd "$run" && find . -name DONE | LC_ALL=C sort | xargs md5sum ) > "$o.done_before"
    drive "$side" "$o" "$@"
    ( cd "$run" && find . -name DONE | LC_ALL=C sort | xargs md5sum ) > "$o.done_after"
    k=$(grep -c 'up to date, skipped' "$o.stdout"); nd=$(grep -c '' "$o.done_before")
    ch=$(grep -c '^\(CMD\|OUT\)' "$o.trace")
    if [ "$(cat "$o.rc")" = 0 ] && [ "$k" = "$nd" ] && [ "$ch" = 0 ] && cmp -s "$o.done_before" "$o.done_after"; then
        ok "$v $side driver on the other run dir: $k/$nd stages up to date, 0 CMD/OUT, DONE unchanged"
    else fail "$v $side driver on the other run dir: rc=$(cat "$o.rc") skipped=$k DONE=$nd CMD/OUT=$ch" "$o.stdout"; fi
    rm -rf -- "$run"
}
