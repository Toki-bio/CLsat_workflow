#!/usr/bin/env bash
# WP0 server probe (PORT_SPEC 6 row 0, 8 risk 1-3): the tool forms the bash driver relies on, as they behave on this host.
# Prints a report; every check line is "OK|FAIL|INFO  <what>  <detail>". Exit 1 when any FAIL.
# Usage: tests/port/probe_servers.sh      (temp files under ${TMPDIR:-/tmp}; set TMPDIR on hosts where /tmp is not allowed)
set -uo pipefail
export LC_ALL=C
GAWK=${SATPIPE_GAWK:-gawk}
TMP=$(mktemp -d "${TMPDIR:-/tmp}/satpipe_probe.XXXXXX") || { echo "FAIL mktemp -d"; exit 1; }
trap 'rm -rf -- "$TMP"' EXIT
fails=0
ok()   { printf 'OK    %-34s %s\n' "$1" "${2:-}"; }
bad()  { printf 'FAIL  %-34s %s\n' "$1" "${2:-}"; fails=$((fails + 1)); }
info() { printf 'INFO  %-34s %s\n' "$1" "${2:-}"; }
chk()  { if [ "$2" = "$3" ]; then ok "$1" "$2"; else bad "$1" "got [$2] want [$3]"; fi; }
first() { "$@" 2>&1 | head -n 1; }

info host "$(hostname 2>/dev/null) $(uname -srm)"
info date "$(date)"
info bash "$BASH_VERSION"
info "gawk (\$SATPIPE_GAWK or PATH)" "$(command -v "$GAWK") :: $(first "$GAWK" --version)"
info "awk on PATH" "$(command -v awk) :: $(first awk --version || first awk -W version)"
info jq "$(command -v jq) :: $(first jq --version)"
info python3 "$(command -v python3) :: $(first python3 --version)"
for t in stat sha256sum md5sum date timeout xargs readlink mktemp tail cp ln rm sort head; do
    info "$t" "$(command -v $t) :: $(first $t --version)"
done

# gawk features used by pyjson.awk and the driver
case $(first "$GAWK" --version) in "GNU Awk"*) ok "gawk is GNU awk" ;; *) bad "gawk is GNU awk" "$(first "$GAWK" --version)" ;; esac
chk "gawk %c 200 is one byte (LC_ALL=C)" "$("$GAWK" 'BEGIN { printf "%d", length(sprintf("%c", 200)) }')" 1
chk "gawk -f lib -e prog" "$(printf 'function f() { return 7 }\n' > "$TMP/l.awk"; "$GAWK" -f "$TMP/l.awk" -e 'BEGIN { print f() }')" 7
chk "gawk asorti @ind_str_asc bytewise" "$("$GAWK" 'BEGIN { a["b"]; a["B"]; a["_"]; a["a"]; a["."]; n = asorti(a, k, "@ind_str_asc"); for (i = 1; i <= n; i++) printf "%s", k[i] }')" '.B_ab'
chk "gawk strtonum/and/RS=^\$" "$(printf 'x\ny' | "$GAWK" 'BEGIN { RS = "^$" } { printf "%d %d %d", strtonum("0xD83D"), and(255, 15), length($0) }')" '55357 15 3'
chk "gawk %.17e round trip of 0.1" "$("$GAWK" 'BEGIN { s = sprintf("%.16e", 0.1); print (s + 0 == 0.1) ? "yes" : "no" }')" yes

# coreutils forms (risk 1: uutils on therioserver)
printf 'ACGT\n' > "$TMP/f.txt"; ln -s f.txt "$TMP/l.txt"; touch -d '2026-01-02 03:04:05' "$TMP/f.txt" 2>/dev/null || true
py_st=$(python3 -c 'import os,sys; s=os.stat(sys.argv[1]); print(s.st_size, int(s.st_mtime))' "$TMP/l.txt" 2>/dev/null)
chk "stat -L -c '%s %Y' (symlink) == os.stat" "$(stat -L -c '%s %Y' "$TMP/l.txt")" "$py_st"
chk "stat -c '%s %Y' (file) == stat -L of link" "$(stat -c '%s %Y' "$TMP/f.txt")" "$py_st"
info "sha256sum < file (raw line)" "[$(sha256sum < "$TMP/f.txt")]"
x=$(sha256sum < "$TMP/f.txt"); chk "sha256sum < file (field 1)" "${x%% *}" "$(python3 -c 'import hashlib,sys; print(hashlib.sha256(open(sys.argv[1],"rb").read()).hexdigest())' "$TMP/f.txt")"
x=$(printf '' | sha256sum); chk "sha256sum of empty stdin" "${x%% *}" e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
x=$(md5sum < "$TMP/f.txt"); chk "md5sum < file (field 1)" "${x%% *}" "$(python3 -c 'import hashlib,sys; print(hashlib.md5(open(sys.argv[1],"rb").read()).hexdigest())' "$TMP/f.txt")"
d=$(date +%FT%T)
if [[ $d =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}$ ]]; then ok "date +%FT%T" "$d (python $(python3 -c 'import datetime as d; print(d.datetime.now().isoformat(timespec="seconds"))'))"; else bad "date +%FT%T" "$d"; fi
timeout 1 sleep 5; chk "timeout 1 sleep 5 -> exit 124" "$?" 124
timeout 30 true; chk "timeout 30 true -> exit 0" "$?" 0
s=$(date +%s)
printf '%s\n' 1 2 3 4 5 6 | xargs -P 3 -I{} sh -c 'sleep 1; echo {}' > "$TMP/x.out"
e=$(( $(date +%s) - s ))
chk "xargs -P 3 -I{} ran all 6" "$("$GAWK" '{ s += $1 } END { print NR, s }' "$TMP/x.out")" "6 21"
if [ "$e" -le 3 ]; then ok "xargs -P 3 parallel (6 x 1 s)" "${e} s"; else bad "xargs -P 3 parallel (6 x 1 s)" "${e} s"; fi
printf '%s\n' 1 2 3 | xargs -P 2 -I{} sh -c 'exit {}' ; chk "xargs -P exit status on failure" "$?" 123
chk "readlink -f symlink" "$(readlink -f "$TMP/l.txt")" "$(cd "$TMP" && pwd -P)/f.txt"
mkdir -p "$TMP/a/b"; ln -sfn a "$TMP/ln"; ln -sfn a/b "$TMP/ln"; chk "ln -sfn replaces dir link" "$(readlink "$TMP/ln")" a/b
seq 1 40 > "$TMP/s.txt"; chk "tail -n 15" "$(tail -n 15 "$TMP/s.txt" | head -n 1)" 26
chk "head -c 5" "$(printf 'abcdefgh' | head -c 5)" abcde
cp "$TMP/f.txt" "$TMP/g.txt"; chk "cp" "$(cat "$TMP/g.txt")" ACGT
rm -rf -- "$TMP/a"; [ ! -e "$TMP/a" ] && ok "rm -rf" || bad "rm -rf"
chk "mkdir -p nested" "$(mkdir -p "$TMP/x/y/z" && [ -d "$TMP/x/y/z" ] && echo y)" y
mkdir "$TMP/glob"; for n in b B _ a .h; do : > "$TMP/glob/$n"; done
chk "bash glob order (dotglob, LC_ALL=C)" "$(cd "$TMP/glob" && shopt -s dotglob && printf '%s,' *)" '.h,B,_,a,b,'
info "sort --version" "$(first sort --version)  (driver never calls sort)"

# jq number-literal preservation (spec 2.4: viewer legacy path needs jq >= 1.7)
if command -v jq >/dev/null 2>&1; then
    j=$(printf '{"b":50.0,"a":0.90,"c":100000000000000000001,"d":1.5e300,"e":"x/\\u00e9"}' | jq -c .)
    info "jq -c . on {b:50.0,a:0.90,c:1e20+1,d:1.5e300}" "$j"
    case $j in *'"b":50.0'*) ok "jq preserves 50.0" ;; *) bad "jq preserves 50.0" "$j" ;; esac
    case $j in *'"c":100000000000000000001'*) ok "jq preserves big int literal" ;; *) bad "jq preserves big int literal" "$j" ;; esac
    chk "jq -S key sort" "$(printf '{"b":1,"a":2}' | jq -S -c .)" '{"a":2,"b":1}'
    chk "jq -r to_entries" "$(printf '{"inputs":{"/x y":"sha256:1"}}' | jq -r '.inputs | to_entries[] | "\(.key)\t\(.value)"')" $'/x y\tsha256:1'
    chk "jq --argjson keeps 50.0" "$(jq -n -c --argjson v '{"x":50.0}' '{v:$v}')" '{"v":{"x":50.0}}'
else
    bad jq "not on PATH"
fi

echo "fails: $fails"
[ "$fails" = 0 ]
