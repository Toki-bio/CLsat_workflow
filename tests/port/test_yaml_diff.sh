#!/usr/bin/env bash
# Differential test of lib/satpipe/yaml.awk against PyYAML (tests/port/py_oracle_raw.py), PORT_SPEC.md section 6 WP2.
#   bash tests/port/test_yaml_diff.sh            # all checks below; exit 1 if any fails
# env: PYTHON (default python3), GAWK (default gawk), FUZZ_SEEDS (default "1"), FUZZ_N (default 200), KEEP=dir to keep work
# Checks:
#   accept: configs/*.yaml, tests/port/yaml_fixtures/**/*.yaml (flank-single-copy and Monsoon configs, test-generated
#           configs from make_yaml_fixtures.py) and the fuzz documents: yaml.awk output must equal the oracle byte for byte
#           and neither may report an error (an equal 'top level must be a mapping' line also counts as equal).
#   refuse: gen_yaml_fuzz.py --refused: yaml.awk must exit 1 with one '!error<TAB>YAML subset: ...' line.
set -euo pipefail
export LC_ALL=C
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
PY=${PYTHON:-python3}; GAWK=${GAWK:-gawk}
W=${KEEP:-$(mktemp -d)}; mkdir -p "$W"
[ -n "${KEEP:-}" ] || trap 'rm -rf "$W"' EXIT
nacc=0; nref=0; fail=0

accept() {   # $1 label, $2 file
  local a="$W/a.out" b="$W/b.out" ra=0 rb=0
  "$GAWK" -f "$ROOT/lib/satpipe/yaml.awk" "$2" > "$a" || ra=$?
  "$PY" "$ROOT/tests/port/py_oracle_raw.py" "$2" > "$b" || rb=$?
  nacc=$((nacc + 1))
  if [ $ra -ne $rb ] || ! cmp -s "$a" "$b" || { [ $ra -ne 0 ] && ! grep -q $'^!error\ttop level must be a mapping$' "$a"; }; then
    fail=$((fail + 1)); echo "DIFF [$1] $2 (gawk rc=$ra, python rc=$rb)"; diff "$a" "$b" | head -10 || true
  fi
}

for f in "$ROOT"/configs/*.yaml; do accept repo-config "$f"; done
while IFS= read -r f; do accept fixture "$f"; done < <(find "$ROOT/tests/port/yaml_fixtures" -name '*.yaml' | LC_ALL=C sort)
echo "fixtures: $nacc compared, $fail differ"
fx=$fail
for seed in ${FUZZ_SEEDS:-1}; do
  mkdir -p "$W/fuzz$seed"
  "$PY" "$ROOT/tests/port/gen_yaml_fuzz.py" --seed "$seed" --n "${FUZZ_N:-200}" "$W/fuzz$seed"
  for f in "$W/fuzz$seed"/*.yaml; do accept "fuzz$seed" "$f"; done
done
echo "fuzz: $((nacc - $(find "$ROOT/configs" "$ROOT/tests/port/yaml_fixtures" -name '*.yaml' | wc -l))) compared, $((fail - fx)) differ"

mkdir -p "$W/refused"
"$PY" "$ROOT/tests/port/gen_yaml_fuzz.py" --refused "$W/refused"
for f in "$W/refused"/*.yaml; do
  nref=$((nref + 1)); r=0
  out=$("$GAWK" -f "$ROOT/lib/satpipe/yaml.awk" "$f") || r=$?
  py=$("$PY" "$ROOT/tests/port/py_oracle_raw.py" "$f" | head -c 70 | tr '\t\n' '  ') || true
  if [ $r -eq 1 ] && [ "$(printf '%s\n' "$out" | wc -l)" -eq 1 ] && [[ $out == $'!error\tYAML subset: '* ]]; then
    printf 'refused %-26s %s   | PyYAML: %s\n' "$(basename "$f" .yaml)" "${out#*$'\t'}" "$py"
  else
    fail=$((fail + 1)); echo "NOT REFUSED $(basename "$f") rc=$r: $out"
  fi
done
echo "refused: $nref cases"
echo "TOTAL: $nacc accept-compared, $nref refuse-checked, $fail failures"
[ $fail -eq 0 ]
