#!/usr/bin/env bash
# WP0 unit tests of lib/satpipe/pyjson.awk + flat.sh against Python-made fixtures (no Python needed to run).
#   1. pyrepr_float on fixtures/floats.txt == fixtures/floats.python_repr.txt          (3027 values)
#   2. json_str / py_repr_str / flat escaping on fixtures/pyjson_strings.tsv            (gen_pyjson_cases.py)
#   3. the 88 PLAN rows of fixtures/probe_params.out: JSON -> TYPED flat -> params_hash == Python's hash
#   4. TYPED-flat helpers (flat_get/keys/subtree/py_repr) on a hand-made params file
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd -- "$HERE/../.." && pwd)
source "$ROOT/lib/satpipe/common.sh"
source "$ROOT/lib/satpipe/flat.sh"
source "$ROOT/lib/satpipe/done.sh"
require_gawk
FX=$HERE/fixtures
TMP=$(mktemp -d "${TMPDIR:-/tmp}/satpipe_wp0.XXXXXX")
trap 'rm -rf -- "$TMP"' EXIT
fail=0
report() { printf '%-44s %s\n' "$1" "$2"; case $2 in PASS*) ;; *) fail=1 ;; esac; }

# 1. floats
pyjson_gawk '{ print pyrepr_float($1 + 0) }' "$FX/floats.txt" > "$TMP/floats.out"
n=$(wc -l < "$FX/floats.txt" | tr -d ' ')
bad=$(pyjson_gawk 'NR == FNR { a[FNR] = $0; next } a[FNR] != $0 { b++ } END { print b + 0 }' "$TMP/floats.out" "$FX/floats.python_repr.txt")
[ "$bad" = 0 ] && report "pyrepr_float vs Python repr" "PASS $n/$n" || report "pyrepr_float vs Python repr" "FAIL $bad of $n differ"

# 2. strings
read -r n bj br nr < <(pyjson_gawk 'BEGIN { FS = "\t" }
    { s = flat_unesc($1); n++
      if (json_str(s) != $2) { bj++; if (bj <= 3) print "json_str mismatch line " NR ": " json_str(s) " vs " $2 > "/dev/stderr" }
      if (flat_esc(s) != $1) bj++
      if ($3 != "-") { nr++; if (flat_esc(py_repr_str(s)) != $3) { br++; if (br <= 3) print "py_repr mismatch line " NR ": " py_repr_str(s) " vs " $3 > "/dev/stderr" } } }
    END { print n + 0, bj + 0, br + 0, nr + 0 }' "$FX/pyjson_strings.tsv")
[ "$bj" = 0 ] && report "json_str (ensure_ascii) vs json.dumps" "PASS $n/$n" || report "json_str vs json.dumps" "FAIL $bj of $n"
[ "$br" = 0 ] && report "py_repr_str vs repr" "PASS $nr/$nr" || report "py_repr_str vs repr" "FAIL $br of $nr"

# 3. probe hashes
n=0; ok=0; okj=0
while IFS=$'\t' read -r tag sid hash js; do
    [ "$tag" = PLAN ] || continue
    n=$((n + 1))
    printf '%s' "$js" > "$TMP/p.json"
    json_to_flat "$TMP/p.json" > "$TMP/p.flat"
    [ "$(params_hash "$TMP/p.flat")" = "$hash" ] && ok=$((ok + 1)) || log "probe hash mismatch: $sid $hash"
    [ "$(flat_json "$TMP/p.flat")" = "$js" ] && okj=$((okj + 1))
done < "$FX/probe_params.out"
[ "$n" = 88 ] && [ "$ok" = 88 ] && report "probe params_hash (88 stage/config pairs)" "PASS $ok/$n" || report "probe params_hash" "FAIL $ok/$n"
[ "$okj" = "$n" ] && report "probe canonical JSON bytes" "PASS $okj/$n" || report "probe canonical JSON bytes" "FAIL $okj/$n"

# 4. flat helpers
{
    flat_emit terminal map 2; flat_emit terminal/window_mono int 3; flat_emit terminal/edge_tol float 0.5
    flat_emit families map 1; flat_emit families/CLsat1.a map 2; flat_emit families/CLsat1.a/colour null None
    flat_emit families/CLsat1.a/members list 2; flat_emit families/CLsat1.a/members/0 str "x'y"; flat_emit families/CLsat1.a/members/1 str $'tab\there\\'
    flat_emit period list 2; flat_emit period/0 int 80; flat_emit period/1 int 300; flat_emit empty map 0; flat_emit on bool True
} > "$TMP/h.flat"
got=$(flat_json "$TMP/h.flat")
want='{"empty": {}, "families": {"CLsat1.a": {"colour": null, "members": ["x'"'"'y", "tab\there\\"]}}, "on": true, "period": [80, 300], "terminal": {"edge_tol": 0.5, "window_mono": 3}}'
[ "$got" = "$want" ] && report "flat_json sort_keys/nesting/escapes" PASS || report "flat_json" "FAIL got $got"
[ "$(flat_get "$TMP/h.flat" families/CLsat1.a/members/1)" = $'tab\there\\' ] && report "flat_get unescape" PASS || report "flat_get" FAIL
[ "$(flat_keys "$TMP/h.flat" "" | tr '\n' ,)" = "terminal,families,period,empty,on," ] && report "flat_keys document order" PASS || report "flat_keys" FAIL
[ "$(flat_subtree "$TMP/h.flat" terminal t2 | cut -f1 | tr '\n' ,)" = "t2,t2/window_mono,t2/edge_tol," ] && report "flat_subtree re-root" PASS || report "flat_subtree" FAIL
[ "$(flat_py_repr "$TMP/h.flat" "")" = "{'terminal': {'window_mono': 3, 'edge_tol': 0.5}, 'families': {'CLsat1.a': {'colour': None, 'members': [\"x'y\", 'tab\\there\\\\']}}, 'period': [80, 300], 'empty': {}, 'on': True}" ] \
    && report "flat_py_repr dict/list/None/True" PASS || report "flat_py_repr" "FAIL $(flat_py_repr "$TMP/h.flat" "")"
: > "$TMP/e.flat"
[ "$(params_hash "$TMP/e.flat")" = 44136fa355b3678a ] && report "empty params hash 44136fa355b3678a" PASS || report "empty params hash" FAIL

exit "$fail"
