#!/usr/bin/env bash
# Schema and stage tables (lib/satpipe/schema.tsv, stages_registry.tsv, run_stages.tsv) - self-contained since the switch-over
# (WP14). schema.tsv is now the source of every parameter; it was generated from the Python driver's schema.py, last at
# tag python-driver-final, and the golden files in tests/sh/golden/ are that driver's own output at the tag:
#   stages.txt (`satpipe stages`), list_stages.txt (`satpipe run X --list-stages`), patterns.tsv (pattern/string/re.match cases).
#   bash tests/sh/test_schema_tsv.sh            (needs gawk)
set -euo pipefail
export LC_ALL=C
ROOT=$(cd "$(dirname "$(readlink -f "$0")")/../.." && pwd)
LIB=$ROOT/lib/satpipe; G=$ROOT/tests/sh/golden
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
cd "$ROOT"
pass=0
check() { local label=$1; shift
          if "$@"; then pass=$((pass+1)); echo "ok   $label"; else echo "FAIL $label"; exit 1; fi; }
same() { cmp -s "$1" "$2" || { diff "$1" "$2" | head -20; return 1; }; }

check "schema.tsv has no CR" bash -c "! grep -q \$'\r' '$LIB/schema.tsv'"
check "schema.tsv rows have 20 columns" gawk -F'\t' 'NF != 20 { print FILENAME ":" NR ": " NF " columns"; bad = 1 } END { exit bad }' "$LIB/schema.tsv"

# --schema-doc renders docs/config.md (after a schema.tsv change: bin/satpipe config --schema-doc > docs/config.md)
gawk -f "$LIB/schema_doc.awk" "$LIB/schema.tsv" > "$TMP/doc.md"
check "schema_doc.awk == docs/config.md ($(wc -c < "$TMP/doc.md") bytes)" same "$TMP/doc.md" docs/config.md

# ERE columns match exactly the strings Python's re.match matched (golden cases)
gawk -F'\t' '{ g = ($3 ~ $2) ? 1 : 0; if (g != $4) { print "mismatch:", $1, $2, "[" $3 "]", "py=" $4, "gawk=" g; bad++ } n++ }
             END { print n " pattern/string cases, " bad+0 " mismatches" > "/dev/stderr"; exit bad > 0 }' "$G/patterns.tsv" > "$TMP/pat.out" \
  || { cat "$TMP/pat.out"; false; }
check "pattern_ere/key_pattern_ere agree with the golden re.match cases" true
gawk -F'\t' 'NR > 1 && $12 != "" { print $12 "\t" $13 } NR > 1 && $14 != "" { print $14 "\t" $15 }' "$LIB/schema.tsv" | sort -u > "$TMP/pat_now"
gawk -F'\t' '{ print $1 "\t" $2 }' "$G/patterns.tsv" | sort -u > "$TMP/pat_gold"
check "every schema pattern has golden cases (new pattern: add cases to tests/sh/golden/patterns.tsv)" \
    bash -c "[ -z \"\$(comm -23 '$TMP/pat_now' '$TMP/pat_gold')\" ] || { comm -23 '$TMP/pat_now' '$TMP/pat_gold'; false; }"

# `satpipe stages` and `satpipe run --list-stages` (bash driver) == the tables / golden output
"$ROOT/bin/satpipe" stages > "$TMP/stages.txt"
check "satpipe stages == golden ($(wc -l < "$G/stages.txt") lines)" same "$TMP/stages.txt" "$G/stages.txt"
"$ROOT/bin/satpipe" run nonexistent.yaml --list-stages > "$TMP/list.txt"
check "satpipe run --list-stages == golden" same "$TMP/list.txt" "$G/list_stages.txt"
gawk -F'\t' 'NR > 1 { print $1 }' "$LIB/run_stages.tsv" > "$TMP/ids.txt"
check "run_stages.tsv ids == --list-stages" same "$TMP/ids.txt" "$G/list_stages.txt"

echo "all $pass checks passed"
