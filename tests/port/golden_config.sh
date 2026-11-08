#!/usr/bin/env bash
# WP3 golden config test (PORT_SPEC 5a): the bash config loader (lib/satpipe/config.sh = yaml.awk | validate.awk |
# crosscheck.awk) against the Python loader (tools/py_oracle.py flat = load_config + TYPED flat).
#   bash tests/port/golden_config.sh          exit 1 on any unexpected difference
# env: PYTHON (default python3), SATPIPE_GAWK (default gawk), KEEP=dir to keep the work dir,
#      FLANK_SCHEMA_PY = schema.py of branch flank-single-copy (default: `git show origin/flank-single-copy:satpipe/schema.py`)
# Checks, each config under both schemas (this checkout's, and flank-single-copy's via a temporary schema.tsv):
#   stream   config_stream == oracle flat, byte for byte (TYPED flat + !warning lines, or the !error list in order);
#            a YAML error compares the first word only (`YAML subset:` vs PyYAML's `YAML parse error:`, PORT_SPEC 3);
#            configs: configs/*.yaml, tests/port/yaml_fixtures/**, tests/port/bad/*.yaml (make_bad_configs.py)
#   cwd      the same with a relative config path from another directory and through a symlinked directory
#   flat     `config.sh --flat FILE` stdout == the oracle's TYPED flat (valid configs)
#   reload   `config.sh FILE` YAML, re-loaded by PyYAML, equals load_config(FILE).data with types (valid configs)
#   tools    the built-in tool-name list == satpipe/tools.py TOOLS
# Expected differences: tests/port/golden_config.expected (config base names, refused-by-decision forms).
set -euo pipefail
export LC_ALL=C
ROOT=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../.." && pwd -P)
PY=${PYTHON:-python3}
W=${KEEP:-$(mktemp -d)}; mkdir -p "$W"; W=$(cd "$W" && pwd -P)
[ -n "${KEEP:-}" ] || trap 'rm -rf "$W"' EXIT
# shellcheck source=../../lib/satpipe/config.sh
. "$ROOT/lib/satpipe/config.sh"
require_gawk

declare -A XFAIL=()
while IFS= read -r l; do l=${l%%#*}; l=${l//[[:space:]]/}; [ -n "$l" ] && XFAIL[$l]=1; done < "$ROOT/tests/port/golden_config.expected"

n_stream=0 n_first=0 n_xfail=0 n_valid=0 n_err=0 n_cwd=0 n_flat=0 n_reload=0 fail=0
bad() { fail=$((fail + 1)); echo "FAIL $*"; }

# schemas: name -> schema.tsv and schema.py
FLANK_PY=${FLANK_SCHEMA_PY:-}
if [ -z "$FLANK_PY" ]; then
    FLANK_PY=$W/flank_schema.py
    git -C "$ROOT" show origin/flank-single-copy:satpipe/schema.py > "$FLANK_PY" 2>/dev/null || FLANK_PY=""
fi
SCHEMAS=("main|$SATPIPE_LIB/schema.tsv|")
if [ -n "$FLANK_PY" ]; then
    "$PY" "$ROOT/tools/gen_schema_tsv.py" schema --schema "$FLANK_PY" -o "$W/flank_schema.tsv"
    tr -d '\r' < "$W/flank_schema.tsv" > "$W/flank_schema.tsv.lf" && mv "$W/flank_schema.tsv.lf" "$W/flank_schema.tsv"
    SCHEMAS+=("flank|$W/flank_schema.tsv|$FLANK_PY")
    echo "flank-single-copy schema: $(grep -c . "$W/flank_schema.tsv") rows (main: $(grep -c . "$SATPIPE_LIB/schema.tsv"))"
else
    bad "flank-single-copy schema.py not available (set FLANK_SCHEMA_PY)"
fi

# compare FILE with the oracle under schema S; CWD = directory to run in, ARG = config argument as given
compare() {   # $1 label, $2 schema tsv, $3 schema py, $4 cwd, $5 config arg, $6 base name
    local a=$W/a.out b=$W/b.out ra=0 rb=0
    ( cd "$4" && SATPIPE_SCHEMA_TSV=$2 config_stream "$5" ) > "$a" || ra=$?
    ( cd "$4" && SATPIPE_SCHEMA_PY=$3 "$PY" "$ROOT/tools/py_oracle.py" flat "$5" ) > "$b.raw" || rb=$?
    tr -d '\r' < "$b.raw" > "$b"
    if [ -n "${XFAIL[$6]:-}" ]; then
        n_xfail=$((n_xfail + 1))
        if cmp -s "$a" "$b"; then bad "[$1] $6: listed as expected difference but equal"; fi
        return 0
    fi
    if [ $ra -ne $rb ]; then bad "[$1] $5: rc bash=$ra python=$rb"; diff "$a" "$b" | head -8 || true; return 0; fi
    if cmp -s "$a" "$b"; then
        n_stream=$((n_stream + 1)); if [ $ra -eq 0 ]; then n_valid=$((n_valid + 1)); else n_err=$((n_err + 1)); fi
        return 0
    fi
    if [ "$(wc -l < "$a")" -eq 1 ] && [ "$(wc -l < "$b")" -eq 1 ] && grep -q $'^!error\tYAML subset: ' "$a" \
       && grep -q $'^!error\tYAML parse error: ' "$b"; then
        n_first=$((n_first + 1)); return 0
    fi
    bad "[$1] $5 (cwd $4)"; diff "$a" "$b" | head -12 || true
}

CONFIGS=()
for f in "$ROOT"/configs/*.yaml; do CONFIGS+=("$f"); done
while IFS= read -r f; do CONFIGS+=("$f"); done < <(find "$ROOT/tests/port/yaml_fixtures" "$ROOT/tests/port/bad" -name '*.yaml' | LC_ALL=C sort)
echo "configs: ${#CONFIGS[@]}"

for s in "${SCHEMAS[@]}"; do
    IFS='|' read -r sname stsv spy <<< "$s"
    for f in "${CONFIGS[@]}"; do compare "$sname/abs" "$stsv" "$spy" "$ROOT" "$f" "$(basename "$f" .yaml)"; done
done
echo "stream: $n_stream equal ($n_valid valid, $n_err error lists), $n_first YAML-error first-word equal, $n_xfail expected differences"

# relative paths from another directory and through a symlinked directory (main schema)
mkdir -p "$W/cwd/sub"; ln -sfn "$ROOT/tests/port" "$W/cwd/link"
for f in "$ROOT"/configs/*.yaml "$ROOT"/tests/port/bad/ok_*.yaml "$ROOT"/tests/port/yaml_fixtures/generated/*.yaml; do
    rel=${f#"$ROOT"/}
    b=$(basename "$f" .yaml)
    s0=$n_stream; f0=$fail
    compare main/rel "$SATPIPE_LIB/schema.tsv" "" "$ROOT/tests" "../$rel" "$b"
    compare main/rel2 "$SATPIPE_LIB/schema.tsv" "" "$W/cwd/sub" "$ROOT/tests/../$rel" "$b"
    case $rel in tests/port/*)
        compare main/link "$SATPIPE_LIB/schema.tsv" "" "$W/cwd" "link/${rel#tests/port/}" "$b"
        compare main/inlink "$SATPIPE_LIB/schema.tsv" "" "$W/cwd/link" "${rel#tests/port/}" "$b" ;;
    esac
    [ $fail -eq $f0 ] && n_cwd=$((n_cwd + n_stream - s0))
done
echo "cwd: $n_cwd relative/symlink comparisons equal"

# config --flat and the YAML dump re-load (valid configs, main schema)
for f in "${CONFIGS[@]}"; do
    b=$(basename "$f" .yaml)
    [ -n "${XFAIL[$b]:-}" ] && continue
    "$PY" "$ROOT/tools/py_oracle.py" flat "$f" > "$W/o.raw" 2>/dev/null || continue
    tr -d '\r' < "$W/o.raw" | grep -v '^!warning' > "$W/o.flat" || true
    bash "$ROOT/lib/satpipe/config.sh" --flat "$f" > "$W/c.flat" 2> "$W/c.err" || { bad "--flat $f: rc != 0"; continue; }
    if cmp -s "$W/c.flat" "$W/o.flat" && [ ! -s "$W/c.err" ]; then n_flat=$((n_flat + 1)); else bad "--flat $f"; diff "$W/c.flat" "$W/o.flat" | head -5 || true; fi
    bash "$ROOT/lib/satpipe/config.sh" "$f" > "$W/dump.yaml" || { bad "config $f: rc != 0"; continue; }
    if "$PY" "$ROOT/tools/py_oracle.py" reload "$f" "$W/dump.yaml" > "$W/r.out" 2>&1; then n_reload=$((n_reload + 1))
    else bad "reload $f"; head -c 600 "$W/r.out"; echo; fi
done
echo "flat: $n_flat equal; reload: $n_reload equal"

# error-path CLI: `config.sh FILE` on an invalid config prints `error: msg` lines on stderr, exit 1; no FILE -> exit 2
r=0; bash "$ROOT/lib/satpipe/config.sh" "$ROOT/tests/port/bad/cross_many.yaml" > "$W/e.out" 2> "$W/e.err" || r=$?
"$PY" -c "
import sys; sys.path.insert(0, sys.argv[1])
from satpipe.config import load_config, ConfigError
try: load_config(sys.argv[2])
except ConfigError as e: sys.stdout.write(''.join('error: ' + m + '\n' for m in e.errors))" "$ROOT" "$ROOT/tests/port/bad/cross_many.yaml" | tr -d '\r' > "$W/e.py"
if [ $r -eq 1 ] && [ ! -s "$W/e.out" ] && cmp -s "$W/e.err" "$W/e.py"; then echo "ok   config FILE error output ($(wc -l < "$W/e.err") lines)"; else bad "config FILE error output (rc $r)"; diff "$W/e.err" "$W/e.py" || true; fi
r=0; bash "$ROOT/lib/satpipe/config.sh" > /dev/null 2> "$W/e.err" || r=$?
if [ $r -eq 2 ] && [ "$(cat "$W/e.err")" = "error: config path required unless --schema-doc" ]; then echo "ok   config without a path: exit 2"; else bad "config without a path (rc $r)"; fi

t_py=$("$PY" "$ROOT/tools/py_oracle.py" tools | tr -d '\r')
if [ "$t_py" = "$(config_known_tools)" ]; then echo "ok   tool names = satpipe/tools.py TOOLS"; else bad "tool names: python '$t_py' bash '$(config_known_tools)'"; fi

echo "TOTAL: stream $n_stream equal + $n_first first-word + $n_xfail expected; cwd $n_cwd; flat $n_flat; reload $n_reload; failures $fail"
[ $fail -eq 0 ]
