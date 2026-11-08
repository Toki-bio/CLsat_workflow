#!/usr/bin/env bash
# WP1 tests: schema.tsv, schema_doc.awk, stages_registry.tsv, run_stages.tsv against the Python driver.
#   bash tests/sh/test_schema_tsv.sh            (needs gawk, python3 or python, git; test code may use sort/diff)
# Env FLANK_REF (default origin/flank-single-copy): a second schema.py version to check; skipped if the ref is absent.
set -euo pipefail
export LC_ALL=C
ROOT=$(cd "$(dirname "$(readlink -f "$0")")/../.." && pwd)
PY=${PYTHON:-$(command -v python3 >/dev/null 2>&1 && echo python3 || echo python)}
LIB=$ROOT/lib/satpipe
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
cd "$ROOT"
pass=0
check() { local label=$1; shift
          if "$@"; then pass=$((pass+1)); echo "ok   $label"; else echo "FAIL $label"; exit 1; fi; }
nocr() { tr -d '\r' < "$1" > "$2"; }          # Windows python may write CRLF on stdout
same() { cmp -s "$1" "$2" || { diff "$1" "$2" | head -20; return 1; }; }

# 1. schema.tsv is the generator's output for this checkout's schema.py (never hand-edited)
"$PY" tools/gen_schema_tsv.py schema -o "$TMP/schema.tsv"
check "schema.tsv up to date with satpipe/schema.py" same "$TMP/schema.tsv" "$LIB/schema.tsv"
check "schema.tsv has no CR" bash -c "! grep -q \$'\r' '$LIB/schema.tsv'"

# 2. --schema-doc: gawk renderer == docs/config.md == Python to_markdown()
gawk -f "$LIB/schema_doc.awk" "$LIB/schema.tsv" > "$TMP/doc.md"
nocr docs/config.md "$TMP/config.md"
check "schema_doc.awk == docs/config.md ($(wc -c < "$TMP/doc.md") bytes)" same "$TMP/doc.md" "$TMP/config.md"
"$PY" tests/port/check_schema_tsv.py --markdown > "$TMP/py.md.raw"; nocr "$TMP/py.md.raw" "$TMP/py.md"
check "schema_doc.awk == Python to_markdown()" same "$TMP/doc.md" "$TMP/py.md"

# 3. field-by-field round trip, type strict (50 vs 50.0)
check "schema.tsv rebuilds SCHEMA field by field" "$PY" tests/port/check_schema_tsv.py "$LIB/schema.tsv"

# 4. ERE columns match exactly the strings the Python regexes match
"$PY" tests/port/check_schema_tsv.py --patterns "$LIB/schema.tsv" > "$TMP/pat.raw"; nocr "$TMP/pat.raw" "$TMP/pat.tsv"
gawk -F'\t' '{ g = ($3 ~ $2) ? 1 : 0; if (g != $4) { print "mismatch:", $1, $2, "[" $3 "]", "py=" $4, "gawk=" g; bad++ } n++ }
             END { print n " pattern/string cases, " bad+0 " mismatches" > "/dev/stderr"; exit bad > 0 }' "$TMP/pat.tsv" > "$TMP/pat.out" \
  || { cat "$TMP/pat.out"; false; }
check "pattern_ere/key_pattern_ere agree with Python re.match" true

# 5. a second schema version (flank-single-copy branch): generated table renders that version's to_markdown
REF=${FLANK_REF:-origin/flank-single-copy}
if git rev-parse -q --verify "$REF" >/dev/null; then
    git show "$REF:satpipe/schema.py" > "$TMP/schema_flank.py"
    "$PY" tools/gen_schema_tsv.py schema --schema "$TMP/schema_flank.py" -o "$TMP/schema_flank.tsv"
    gawk -f "$LIB/schema_doc.awk" "$TMP/schema_flank.tsv" > "$TMP/doc_flank.md"
    "$PY" tests/port/check_schema_tsv.py --markdown "$TMP/schema_flank.py" > "$TMP/py_flank.raw"; nocr "$TMP/py_flank.raw" "$TMP/py_flank.md"
    check "$REF: schema_doc.awk == its to_markdown() ($(wc -c < "$TMP/doc_flank.md") bytes)" same "$TMP/doc_flank.md" "$TMP/py_flank.md"
    git show "$REF:docs/config.md" | tr -d '\r' > "$TMP/config_flank.md"
    check "$REF: schema_doc.awk == its docs/config.md" same "$TMP/doc_flank.md" "$TMP/config_flank.md"
    check "$REF: schema.tsv rebuilds its SCHEMA" "$PY" tests/port/check_schema_tsv.py "$TMP/schema_flank.tsv" "$TMP/schema_flank.py"
else
    echo "skip $REF not available (git fetch origin flank-single-copy)"
fi

# 6. stage tables are the generator's output, and reproduce the Python driver
"$PY" tools/gen_schema_tsv.py stages -o "$TMP/stages_registry.tsv"
check "stages_registry.tsv up to date" same "$TMP/stages_registry.tsv" "$LIB/stages_registry.tsv"
"$PY" tools/gen_schema_tsv.py runorder -o "$TMP/run_stages.tsv"
check "run_stages.tsv up to date" same "$TMP/run_stages.tsv" "$LIB/run_stages.tsv"

# `satpipe stages` rendered from stages_registry.tsv (cli.py:42-56 layout) == Python output
gawk -F'\t' '
    function lst(s,   a, k, i, o) { k = split(s, a, ","); o = ""; for (i = 1; i <= k; i++) o = o (i > 1 ? ", " : "") a[i]; return o }
    NR == 1 { printf "%-12s %-5s %-46s %s\n", "stage", "wired", "name", "scripts"; next }
    { sc = lst($7); if (sc == "") sc = "(command in docs/files.md)"
      printf "%-12s %-5s %-46s %s\n", $1, ($3 == "True" ? "yes" : "no"), $2, sc
      nn = 0; delete notes
      if ($4 == "True") notes[++nn] = "optional"
      if ($5 != "") notes[++nn] = "hand-curated input: " $5
      if ($6 != "") notes[++nn] = "needs " lst($6)
      if (nn) { s = notes[1]; for (i = 2; i <= nn; i++) s = s "; " notes[i]; printf "%-12s %-5s %-46s [%s]\n", "", "", "", s } }
' "$LIB/stages_registry.tsv" > "$TMP/stages.sh.txt"
"$PY" -m satpipe.cli stages > "$TMP/stages.py.raw"; nocr "$TMP/stages.py.raw" "$TMP/stages.py.txt"
check "satpipe stages from stages_registry.tsv == Python ($(wc -l < "$TMP/stages.py.txt") lines)" same "$TMP/stages.sh.txt" "$TMP/stages.py.txt"

# run order + needs == run.py STAGES; ids == `satpipe run --list-stages`
"$PY" -c 'import sys; sys.path.insert(0, "."); from satpipe.run import STAGES
print("#id\tneeds"); [print(f.sid + "\t" + ",".join(f.needs)) for f in STAGES]' > "$TMP/run.py.raw"; nocr "$TMP/run.py.raw" "$TMP/run.py.tsv"
check "run_stages.tsv == run.py STAGES ids/order/needs ($(($(wc -l < "$TMP/run.py.tsv")-1)) stages)" same "$LIB/run_stages.tsv" "$TMP/run.py.tsv"
"$PY" -m satpipe.cli run nonexistent.yaml --list-stages > "$TMP/list.raw"; nocr "$TMP/list.raw" "$TMP/list.txt"
gawk -F'\t' 'NR > 1 { print $1 }' "$LIB/run_stages.tsv" > "$TMP/ids.txt"
check "run_stages.tsv ids == satpipe run --list-stages" same "$TMP/ids.txt" "$TMP/list.txt"

echo "all $pass checks passed"
