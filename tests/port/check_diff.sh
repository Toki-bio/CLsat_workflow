#!/usr/bin/env bash
# WP4 equivalence: `bin/satpipe` vs `python -m satpipe.cli` on the same machine, same environment (PORT_SPEC 5a, 6 WP4).
#   bash tests/port/check_diff.sh [GROUP...]      groups: configs fixtures project tools cli  (default: all)
# env: PYTHON (default python3; needs PyYAML), SATPIPE_GAWK (default gawk), KEEP=dir to keep the work dir,
#      SYS_PATH (default /usr/bin:/bin) = a PATH without the bioinformatics tools (the "tools missing" runs)
# Compared per case: stdout byte for byte (CR stripped) and exit code; stderr only where the CLI prints there.
# Allowed differences (counted, not failures): a YAML error line compares by its first word (`YAML subset:` vs
# `YAML parse error:`, PORT_SPEC 3) and the refusals listed in tests/port/golden_config.expected.
# Groups:
#   configs   configs/*.yaml, tests/port/yaml_fixtures/**, tests/port/bad/*.yaml, the make_yaml_fixtures.py output,
#             each with: --no-files --no-tools | + --json | --no-tools (files) | --no-files (tools) | all, run from the
#             repo root with repo-relative config paths
#   fixtures  the conftest project and every edit test_check.py / test_schema.py make to it, files + tools, text + json
#   project   check with real files: FASTA/BED/table edge cases (CRLF, blank and # lines, bad ints, >5 records, subfamily dir)
#   tools     tool probe: missing tools (SYS_PATH and an empty PATH), fake tools in bin_dir/overrides (stderr version, old
#             version, unrecognised output, no shebang, missing interpreter, CRLF, no output, timeout when SLOW=1)
#   cli       stages, run --list-stages, config --schema-doc / FILE / --flat (vs tools/py_oracle.py flat), argparse errors
set -uo pipefail
export LC_ALL=C
ROOT=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../.." && pwd -P)
PY=${PYTHON:-python3}
SP=$ROOT/bin/satpipe
SYS_PATH=${SYS_PATH:-/usr/bin:/bin}
W=${KEEP:-$(mktemp -d "${TMPDIR:-/tmp}/wp4diff.XXXXXX")}; mkdir -p "$W"; W=$(cd "$W" && pwd -P)
[ -n "${KEEP:-}" ] || trap 'rm -rf "$W"' EXIT
GROUPS_=${*:-configs fixtures project tools cli}
export PYTHONPATH=$ROOT

declare -A XFAIL=()
while IFS= read -r l; do l=${l%%#*}; l=${l//[[:space:]]/}; [ -n "$l" ] && XFAIL[$l]=1; done < "$ROOT/tests/port/golden_config.expected"

n_eq=0 n_yaml=0 n_xfail=0 n_tb=0 fail=0
bad() { fail=$((fail + 1)); echo "FAIL $*"; }
norm() {   # CR strip + YAML error compared by its first word (PyYAML's text spans several lines in the text report)
    tr -d '\r' < "$1" | gawk '
        skip && !/^\[(OK  |WARN|FAIL)\] / && !/^check (passed|FAILED): / { next }
        { skip = 0 }
        /^\[FAIL\] config: YAML (subset|parse error): / { print "[FAIL] config: YAML <error>"; skip = 1; next }
        /^      "message": "YAML (subset|parse error): / { print "      \"message\": \"YAML <error>\""; next }
        { print }'
}
# pair LABEL XFAIL_KEY CWD ENVSPEC -- ARGS...   ENVSPEC: "-" or VAR=value[,VAR=value]  (applied to both drivers)
pair() {
    local label=$1 key=$2 cwd=$3 envs=$4; shift 5
    local -a ev=()
    [ "$envs" = - ] || IFS=, read -r -a ev <<< "$envs"
    local rp=0 rs=0
    ( cd "$cwd" && env "${ev[@]}" "$PY" -m satpipe.cli "$@" ) > "$W/p.out" 2> "$W/p.err" || rp=$?
    ( cd "$cwd" && env "${ev[@]}" bash "$SP" "$@" ) > "$W/s.out" 2> "$W/s.err" || rs=$?
    norm "$W/p.out" > "$W/p.n"; norm "$W/s.out" > "$W/s.n"
    norm "$W/p.err" > "$W/p.en"; norm "$W/s.err" > "$W/s.en"
    if [ -n "$key" ] && [ -n "${XFAIL[$key]:-}" ]; then n_xfail=$((n_xfail + 1)); return 0; fi
    if [ $rp = $rs ] && cmp -s "$W/p.n" "$W/s.n" && { [ "${CMP_ERR:-0}" = 0 ] || cmp -s "$W/p.en" "$W/s.en"; }; then
        if cmp -s "$W/p.out" "$W/s.out" || ! grep -q 'YAML <error>' "$W/p.n"; then n_eq=$((n_eq + 1)); else n_yaml=$((n_yaml + 1)); fi
        return 0
    fi
    if [ $rp = $rs ] && grep -q '^Traceback' "$W/p.err"; then         # Python crashed (not a report); bash exits the same way
        n_tb=$((n_tb + 1)); echo "PYTRACEBACK $label: $(tail -n 1 "$W/p.err")  | bash: $(tail -n 1 "$W/s.err")"; return 0
    fi
    bad "$label: rc py=$rp sh=$rs  ($*)"
    diff "$W/p.n" "$W/s.n" | head -10
    if [ "${CMP_ERR:-0}" = 1 ]; then diff "$W/p.en" "$W/s.en" | head -6; fi
    if grep -q Traceback "$W/p.err"; then echo "  (python traceback: $(tail -n 1 "$W/p.err"))"; fi
}

MODES=("--no-files --no-tools" "--no-files --no-tools --json" "--no-tools" "--no-files" "")
group_configs() {
    local f rel m; local -a cf=()
    local gen=tests/port/_gen_wp4_${CONFIGS_PART//\//_}
    "$PY" "$ROOT/tests/port/make_yaml_fixtures.py" "$ROOT/$gen" >/dev/null || bad "make_yaml_fixtures.py"
    for f in "$ROOT"/configs/*.yaml "$ROOT"/tests/port/yaml_fixtures/*/*.yaml "$ROOT"/tests/port/bad/*.yaml "$ROOT/$gen"/*.yaml; do cf+=("${f#"$ROOT/"}"); done
    if [ -n "${CONFIGS_PART:-}" ]; then          # i/n: every n-th config from the i-th (parallel runs)
        local i=${CONFIGS_PART%/*} n=${CONFIGS_PART#*/} j; local -a sel=()
        for ((j = i - 1; j < ${#cf[@]}; j += n)); do sel+=("${cf[$j]}"); done; cf=("${sel[@]}")
    fi
    echo "configs: ${#cf[@]} files x ${#MODES[@]} modes"
    for rel in "${cf[@]}"; do
        for m in "${MODES[@]}"; do
            # shellcheck disable=SC2086
            pair "$rel [$m]" "$(basename "$rel" .yaml)" "$ROOT" - -- check "$rel" $m
        done
    done
    rm -rf "${ROOT:?}/$gen"
}

# the conftest project (tests/sh/fixture.sh) and the edits of tests/test_check.py / test_schema.py
group_fixtures() {
    . "$ROOT/tests/sh/fixture.sh"; unset SATPIPE_MAX_BASE_PATH
    local d=$W/fx k m
    local -a edits=(base mono26 table34 family42 genome50 bed58 table97 bindir82 anchor125 tools129 nested36 cross59 trfcode70 trfpartial75 hapdup115 hapdup124)
    for k in "${edits[@]}"; do
        rm -rf "$d"; make_project "$d"; local c=$d/configs/run.yaml
        case $k in
            mono26) edit "$c" "mono: 148" "mono: 146" ;;
            table34) printf '0\tSatX\n7\tSatQ\n' > "$d/ref/cluster_family.tsv" ;;
            family42) edit "$c" $'families:\n' $'families:\n  SatY: {consensus: ../ref/SatX.cons.fa, mono: 148}\n' ;;
            genome50) rm -f "$d/genomes/spb.fna" ;;
            bed58) printf 'scf\t10\t5\tspa|x\n' > "$d/ref/manual.bed" ;;
            table97) printf 'x\tSatX\n' > "$d/ref/cluster_family.tsv" ;;
            bindir82) mkdir -p "$d/fakebin"; printf '#!/bin/sh\necho '"'"'vsearch v2.30.0_linux_x86_64'"'"' >&2\n' > "$d/fakebin/vsearch"
                      chmod 755 "$d/fakebin/vsearch"; append "$c" $'\ntools:\n  bin_dir: '"$d/fakebin"$'\n' ;;
            anchor125) append "$c" $'\nanchor: {family: NoSuchFamily}\n' ;;
            tools129) append "$c" $'\ntools: {ssearch_36: /opt/x/ssearch36}\n' ;;
            nested36) append "$c" $'\nexpansion:\n  min_bitz: 3\n' ;;
            cross59) append "$c" $'\nsubfamily:\n  members:\n    SatY: [a]\n  consensuses: ../ref\nviewer:\n  species_names:\n    zzz: Z\n' ;;
            trfcode70) edit "$c" "spb: ../trf/spb.dat" "spc: ../trf/spb.dat" ;;
            trfpartial75) edit "$c" $'    spb: ../trf/spb.dat\n' "" ;;
            hapdup115) append "$c" $'\nhapdup:\n  pairs: [[spa, spa]]\n' ;;
            hapdup124) append "$c" $'\nhapdup:\n  pairs: [[spa]]\n' ;;
        esac
        for m in "${MODES[@]}" "--json"; do
            # shellcheck disable=SC2086
            pair "fixture $k [$m]" "" "$d/configs" - -- check run.yaml $m
        done
        pair "fixture $k MAX_BASE_PATH" "" "$d" "TMPDIR=$d/a_rather_long_temporary_directory_name/that_goes_on_and_on_beyond_seventy" -- check configs/run.yaml --no-tools
    done
}

# real-file edge cases on a project copy
group_project() {
    . "$ROOT/tests/sh/fixture.sh"; unset SATPIPE_MAX_BASE_PATH
    local d=$W/pj k c
    local -a cases=(crlf blank_hash short_line bad_cluster dup_cluster two_fams_missing bed_edge bed_track bed_bigint bed_underscore
                    bed_ws cons_multi cons_none cons_crlf fai rot_ok rot_missing sub_dir sub_dir_empty sub_file sub_missing seeds
                    outdir_exists outdir_empty outdir_noparent nonutf8 trfdat_missing tab_name)
    for k in "${cases[@]}"; do
        rm -rf "$d"; make_project "$d"; c=$d/configs/run.yaml
        case $k in
            crlf) printf '0\tSatX\r\n5\tSatX\r\n' > "$d/ref/cluster_family.tsv"; printf 'spa_scf1\t0\t10\tspa|a\r\n' > "$d/ref/manual.bed" ;;
            blank_hash) printf '\n   \n#x\n0\tSatX\n\t\n' > "$d/ref/cluster_family.tsv" ;;
            short_line) printf '0\tSatX\n17\n  x y  \n' > "$d/ref/cluster_family.tsv" ;;
            bad_cluster) printf '0\tSatX\n-1\tSatX\n1.5\tSatX\n 7 \t SatX \n007\tSatX\n' > "$d/ref/cluster_family.tsv" ;;
            dup_cluster) printf '0\tSatX\n0\tSatX\n' > "$d/ref/cluster_family.tsv" ;;
            two_fams_missing) edit "$c" $'families:\n' $'families:\n  SatY: {consensus: ../ref/SatX.cons.fa, mono: 148}\n'
                              printf '0\tSatX\n3\tSatY\n4\tSatQ\n' > "$d/ref/cluster_family.tsv" ;;
            bed_edge) printf 'a\t0\t0\tspa|x\nb\t-1\t5\tzz|y\nc\t5\t6\nd\t1\n\ne\tx\t5\tspb\n' > "$d/ref/manual.bed" ;;
            bed_track) printf 'track name=x\nbrowser position\n#c\nspa_scf1\t1\t2\tspa\n' > "$d/ref/manual.bed" ;;
            bed_bigint) printf 'a\t99999999999999999998\t99999999999999999999\tspa|x\nb\t99999999999999999999\t99999999999999999998\tspa|x\n' > "$d/ref/manual.bed" ;;
            bed_underscore) printf 'a\t1_0\t2_0\tspa|x\nb\t1__0\t20\tspa|x\nc\t+1\t+2\tspa\n' > "$d/ref/manual.bed" ;;
            bed_ws) printf 'a\t 1 \t 2\t spa|x\n' > "$d/ref/manual.bed" ;;
            cons_multi) printf '>a\nAC\n>b\nAC\n>c\n>d\n>e\n>f\n>g\n' > "$d/ref/SatX.cons.fa" ;;
            cons_none) printf 'ACGT\n' > "$d/ref/SatX.cons.fa" ;;
            cons_crlf) printf '>SatX desc here\r\nACGT ACGT\r\n\r\n  ACG  \r\n' > "$d/ref/SatX.cons.fa" ;;
            fai) : > "$d/genomes/spa.fna.fai" ;;
            rot_ok) printf '>SatX x\nACGT\n' > "$d/ref/rot.fa"; append "$c" $'references:\n  rotated_consensus: ../ref/rot.fa\n' ;;
            rot_missing) printf '>Other\nACGT\n' > "$d/ref/rot.fa"
                         edit "$c" $'families:\n' $'families:\n  SatY: {consensus: ../ref/SatX.cons.fa, mono: 148}\n'
                         printf '0\tSatX\n1\tSatY\n' > "$d/ref/cluster_family.tsv"
                         append "$c" $'references:\n  rotated_consensus: ../ref/rot.fa\n' ;;
            sub_dir) mkdir -p "$d/subx"; for x in a.fa b.fas c.fasta d.fna e.txt .h.fa f.FA; do printf '>x\nA\n' > "$d/subx/$x"; done
                     append "$c" $'subfamily:\n  consensuses: ../subx\n  members:\n    SatX: [a, b, z, h, .h, f]\n' ;;
            sub_dir_empty) mkdir -p "$d/subx"; append "$c" $'subfamily:\n  consensuses: ../subx\n  members:\n    SatX: [a]\n' ;;
            sub_file) printf '>a x\nA\n>b\nC\n>a\nG\n' > "$d/subs.fa"; append "$c" $'subfamily:\n  consensuses: ../subs.fa\n  members:\n    SatX: [a, q]\n' ;;
            sub_missing) append "$c" $'subfamily:\n  consensuses: ../nope\n' ;;
            seeds) printf 'spa_scf1\t0\t10\tspa|x\nspa_scf1\t5\t3\tqq|x\n' > "$d/ref/seed.bed"
                   append "$c" $'\n' ; edit "$c" $'inputs:\n' $'inputs:\n  seeds:\n    SatX: ../ref/seed.bed\n    ZZ: ../ref/none.bed\n' ;;
            outdir_exists) mkdir -p "$d/runs/toy"; : > "$d/runs/toy/.hidden" ;;
            outdir_empty) mkdir -p "$d/runs/toy" ;;
            outdir_noparent) edit "$c" "out_dir: ../runs/toy" "out_dir: ../runs/deep/toy" ;;
            nonutf8) printf 'spa_scf1\t0\t10\tq\xff\xfex\xe2\x82|y\n' > "$d/ref/manual.bed"; printf '>Sat\xc3\xa9X\xff\nAC\xc3\xa9GT\n' > "$d/ref/SatX.cons.fa" ;;
            trfdat_missing) rm -f "$d/trf/spa.dat" ;;
            tab_name) printf 'spa_scf1\t0\t10\tsp\x07a\\x\t"q'"'"'\n' > "$d/ref/manual.bed" ;;
        esac
        pair "project $k" "" "$d/configs" - -- check run.yaml --no-tools
        pair "project $k json" "" "$d/configs" - -- check run.yaml --no-tools --json
    done
}

group_tools() {
    . "$ROOT/tests/sh/fixture.sh"; unset SATPIPE_MAX_BASE_PATH
    local d=$W/tl c gawk_abs bash_abs
    gawk_abs=$(command -v "${SATPIPE_GAWK:-gawk}"); bash_abs=$(command -v bash)
    rm -rf "$d"; make_project "$d"; c=$d/configs/run.yaml
    pair "tools: current PATH" "" "$d/configs" - -- check run.yaml --no-files
    pair "tools: current PATH json" "" "$d/configs" - -- check run.yaml --no-files --json
    pair "tools: SYS_PATH" "" "$d/configs" "PATH=$SYS_PATH,SATPIPE_GAWK=$gawk_abs" -- check run.yaml --no-files
    mkdir -p "$d/empty"
    # an empty PATH: python by absolute path, the bash driver by absolute bash (its shebang needs env + PATH); stdout only
    # (a site BASH_ENV, e.g. Monsoon's lmod init, complains on stderr when PATH is empty)
    local rp=0 rs=0 pyabs; pyabs=$(command -v "$PY")
    ( cd "$d/configs" && PATH=$d/empty "$pyabs" -m satpipe.cli check run.yaml --no-files ) > "$W/p.out" 2> "$W/p.err" || rp=$?
    ( cd "$d/configs" && PATH=$d/empty SATPIPE_GAWK=$gawk_abs "$bash_abs" "$SP" check run.yaml --no-files ) > "$W/s.out" 2> "$W/s.err" || rs=$?
    if [ $rp = $rs ] && cmp -s "$W/p.out" "$W/s.out"; then n_eq=$((n_eq + 1)); else bad "tools: empty PATH rc $rp/$rs"; diff "$W/p.out" "$W/s.out" | head; fi
    # fake tools
    local f=$d/fake; mkdir -p "$f"
    mk() { printf '%s\n' "$2" > "$f/$1"; chmod 755 "$f/$1"; }
    mk vsearch $'#!/bin/sh\necho "vsearch v2.30.0_linux_x86_64" >&2'
    mk samtools $'#!/bin/sh\necho "samtools 1.9"\necho "Using htslib 1.9"'
    mk mafft $'#!/bin/sh\nprintf "MAFFT\\r\\nv7.526 (2024/Apr/26)\\r\\n" >&2'
    mk bedtools $'#!/bin/sh\necho "   "; echo "  this tool prints a long line that is not a version string at all, it goes on and on past eighty characters"'
    mk parallel $'#!/bin/sh\nexit 3'
    mk trf 'not a script, no shebang'
    mk ssearch36 $'#!/nonexistent/interpreter\necho x'
    mk python3 $'#!/bin/sh\necho "Python 3.9.2"'
    mk gawk $'#!/bin/sh\necho "GNU Awk 3.1.8"'
    mk bash $'#!/bin/sh\necho "GNU bash, version 5.2.15(1)-release (x86_64-pc-linux-gnu)"; exit 1'
    mk cons $'#!/bin/sh\nprintf "a\\vb\\fc\\034d"'
    mk water $'#!/bin/sh\nprintf "EMBOSS:6.6.0.0\\n" >&2; exit 0'
    mk blastn $'#!/bin/sh\necho "blastn: 2.17.0+"; echo " Package: blast 2.17.0, build"'
    printf 'x\n' > "$f/makeblastdb"; chmod 644 "$f/makeblastdb"            # not executable: located on PATH instead
    append "$c" $'\ntools:\n  bin_dir: '"$f"$'\n  water: '"$f"$'/water\n  blastn: /nonexistent/blastn\n'
    pair "tools: fake bin_dir" "" "$d/configs" - -- check run.yaml --no-files
    pair "tools: fake bin_dir json" "" "$d/configs" - -- check run.yaml --no-files --json
    pair "tools: fake bin_dir SYS_PATH" "" "$d/configs" "PATH=$SYS_PATH,SATPIPE_GAWK=$gawk_abs" -- check run.yaml --no-files
    rm -rf "$d"; make_project "$d"
    append "$c" $'\ntools:\n  bin_dir: /nonexistent/bin\n  vsearch: '"$f"$'/vsearch\n'
    pair "tools: bin_dir missing + override" "" "$d/configs" - -- check run.yaml --no-files
    pair "tools: config error + tools" "" "$d/configs" - -- check "$ROOT/tests/port/bad/tools_unknown.yaml"
    if [ "${SLOW:-0}" = 1 ]; then
        mkdir -p "$f"; mk water $'#!/bin/sh\nexec sleep 40'
        rm -rf "$d"; make_project "$d"; append "$c" $'\ntools:\n  water: '"$f"$'/water\n'
        pair "tools: timeout" "" "$d/configs" - -- check run.yaml --no-files
    fi
}

group_cli() {
    . "$ROOT/tests/sh/fixture.sh"; unset SATPIPE_MAX_BASE_PATH
    local d=$W/cl c f
    rm -rf "$d"; make_project "$d"; c=$d/configs/run.yaml
    CMP_ERR=1
    pair "stages" "" "$ROOT" - -- stages
    pair "run --list-stages" "" "$ROOT" - -- run --list-stages "$c"
    pair "run CFG --list-stages --from x" "" "$ROOT" - -- run "$c" --from=x --list-stages
    pair "config --schema-doc" "" "$ROOT" - -- config --schema-doc
    pair "config --sch x" "" "$ROOT" - -- config --sch x
    pair "config (no file)" "" "$ROOT" - -- config
    pair "config bad" "" "$ROOT" - -- config tests/port/bad/unknown_several.yaml
    pair "--version" "" "$ROOT" - -- --version
    pair "-h" "" "$ROOT" - -- -h
    for f in check config stages run viewer anchor; do pair "$f -h" "" "$ROOT" - -- "$f" -h; done
    pair "no args" "" "$ROOT" - --
    pair "bogus cmd" "" "$ROOT" - -- bogus
    pair "--bogus" "" "$ROOT" - -- --bogus
    pair "--bogus check" "" "$ROOT" - -- --bogus check "$c" --no-tools
    pair "check" "" "$ROOT" - -- check
    pair "check a b" "" "$ROOT" - -- check a b
    pair "check --no" "" "$ROOT" - -- check a --no
    pair "check --js=1" "" "$ROOT" - -- check a --js=1
    pair "check -- a" "" "$ROOT" - -- check -- a
    pair "check abbrev" "" "$ROOT" - -- check "$c" --no-t --no-f --js
    pair "check missing file" "" "$ROOT" - -- check nonexistent.yaml
    pair "check dir" "" "$ROOT" - -- check configs
    pair "stages x" "" "$ROOT" - -- stages x
    pair "run" "" "$ROOT" - -- run
    pair "run --list-stages" "" "$ROOT" - -- run --list-stages
    pair "run --from" "" "$ROOT" - -- run cfg --from
    pair "run --f x" "" "$ROOT" - -- run cfg --f x
    pair "run --from -x" "" "$ROOT" - -- run cfg --from -x
    pair "run cfg -5" "" "$ROOT" - -- run cfg -5
    pair "run check fails" "" "$ROOT" - -- run tests/port/bad/unknown_top.yaml --skip-tool-check
    pair "anchor" "" "$ROOT" - -- anchor cfg
    pair "anchor extra" "" "$ROOT" - -- anchor cfg a --check b
    pair "anchor -o" "" "$ROOT" - -- anchor cfg -o
    pair "anchor bad cfg" "" "$ROOT" - -- anchor tests/port/bad/unknown_top.yaml x.fa
    pair "anchor missing files" "" "$ROOT" - -- anchor "$c" x.fa
    # anchor on the Darevskia consensuses (configs/darevskia_clsat.yaml has real consensus files; genomes are absent -> check fails)
    pair "anchor darevskia" "" "$ROOT" - -- anchor configs/darevskia_clsat.yaml --check satpipe/stages/02_seed_expansion_scan/consensus
    # a config whose files all exist: the anchor rotation itself
    local an=$W/an; rm -rf "$an"; make_project "$an"
    cp "$ROOT/satpipe/stages/02_seed_expansion_scan/consensus/CLsat1.cons.fa" "$an/ref/SatX.cons.fa"
    local mono; mono=$(gawk '!/^>/ { n += length($0) } END { print n }' "$an/ref/SatX.cons.fa")
    edit "$an/configs/run.yaml" "mono: 148" "mono: $mono"
    append "$an/configs/run.yaml" $'anchor: {motif: AAAGAAG}\n'
    cp "$ROOT/satpipe/stages/02_seed_expansion_scan/consensus/CLsat2.cons.fa" "$an/q.fa"
    pair "anchor --check real" "" "$an" - -- anchor configs/run.yaml --check q.fa ref/SatX.cons.fa
    CMP_ERR=0
    # config FILE (YAML dump: semantic, PORT_SPEC 4 'config dump') and config --flat (vs tools/py_oracle.py flat)
    local n_dump=0 n_flat=0
    for f in "$ROOT"/configs/*.yaml "$c"; do
        (cd "$ROOT" && bash "$SP" config "$f") > "$W/s.yaml" || { bad "config $f rc"; continue; }
        if "$PY" - "$f" "$W/s.yaml" <<'EOF'
import sys, yaml
from satpipe.config import load_config
d = load_config(sys.argv[1]).data
with open(sys.argv[2], encoding="utf-8") as fh:
    e = yaml.safe_load(fh)
def same(a, b):
    if type(a) is not type(b): return False
    if isinstance(a, dict): return list(a) == list(b) and all(same(a[k], b[k]) for k in a)
    if isinstance(a, list): return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    return a == b
sys.exit(0 if same(d, e) else 1)
EOF
        then n_dump=$((n_dump + 1)); else bad "config dump $f"; fi
        (cd "$ROOT" && bash "$SP" config --flat "$f") > "$W/s.flat" || { bad "config --flat $f rc"; continue; }
        (cd "$ROOT" && "$PY" tools/py_oracle.py flat "$f") | tr -d '\r' | grep -v '^!' > "$W/p.flat"
        if cmp -s "$W/s.flat" "$W/p.flat"; then n_flat=$((n_flat + 1)); else bad "config --flat $f"; diff "$W/p.flat" "$W/s.flat" | head -5; fi
    done
    echo "config dump re-load equal: $n_dump; config --flat equal: $n_flat"
}

for g in $GROUPS_; do echo "== $g"; "group_$g"; done
echo "check_diff: $n_eq equal, $n_yaml YAML-error first-word equal, $n_xfail expected differences, $n_tb Python tracebacks (same exit), $fail failures"
[ $fail = 0 ]
