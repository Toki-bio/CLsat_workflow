#!/usr/bin/env bash
# WP11 (c): whole-pipeline equivalence on the synthetic toy. For one VARIANT, the Python driver (python -m satpipe.cli,
# test-only trace hook) and the bash driver (bin/satpipe) each do ONE fresh `satpipe run CFG` (00_tables .. 21_anchor) on
# the SAME run-directory path (run moved aside after each), then:
#   rc / stdout / stderr / trace (PLAN, CMD, OUT) compared (masked; parallel 07_mafft order as a set);
#   every run-dir file: data cmp, DONE with `finished` masked, manifest jq -S (runs dropped), logs masked (own category);
#   cross rerun: each driver on a copy of the other's run dir must run nothing (every stage 'up to date, skipped');
#   `satpipe viewer RUN` by each driver on its own run dir: viewer/ via diff_viewer.sh, viewer_flankscan/ + stdout.
# variant `default` also runs the Python driver a second time (control: the same driver twice, shows what masking must
# cover; any difference there is nondeterminism of the stage scripts, not of the port).
# The curated cluster->family table comes from a separate prep run (Python driver, --to 08_cross, own out_dir, deleted
# afterwards), derived as tests/toy/run_driver_toy.sh does; `seeds` also takes seeds/*.bed of a prep run to 10_seeds.
# Variants: default | manual (manual_loci BED: two planted arrays re-entered, one with family in col 7) | seeds (inputs.seeds
#   for every family: stages 03-10 skipped) | hapdup (third genome spa2 = spa renamed + 3 array-only contigs cut from
#   planted spa arrays; hapdup.pairs [[spa, spa2]]) | samenames (make_toy.py --same-names).
# Usage: tests/port/e2e_toy.sh WORKDIR VARIANT     env: PYTHON (python3 with PyYAML), JOBS (8), SIZE (1200000)
# Needs the stage tools on PATH (trf parallel vsearch mafft ssearch36 samtools bedtools minimap2 gawk jq); servers only.
set -uo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
W=${1:?usage: e2e_toy.sh WORKDIR VARIANT}; V=${2:?variant}
PY=${PYTHON:-python3}; GAWK=${SATPIPE_GAWK:-gawk}; JOBS=${JOBS:-8}
export LC_ALL=C
unset SATPIPE_TRACE SATPIPE_DRY SATPIPE_STAGE_FILES
rm -rf -- "$W"; mkdir -p "$W/out"; W=$(cd "$W" && pwd)
# shellcheck source=e2e_lib.sh
. "$ROOT/tests/port/e2e_lib.sh"
TOY=$W/toy; RUN=$W/run; TAB=$W/cluster_family.tsv
CONS=$ROOT/satpipe/stages/02_seed_expansion_scan/consensus
echo "### e2e toy variant $V  ($(date +%FT%T), $(hostname), repo $(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null))"

"$PY" "$ROOT/tests/toy/make_toy.py" "$TOY" "$CONS" --seed "${SEED:-7}" --size "${SIZE:-1200000}" ${DIMER:+--dimer} \
    $([ "$V" = samenames ] && echo --same-names) > "$W/make_toy.log" 2>&1 || { echo "make_toy failed"; cat "$W/make_toy.log"; exit 1; }
GEN3=""
if [ "$V" = hapdup ]; then
    # spa2: every spa contig renamed h2_*, plus h2_arr1..3 = exact copies of three planted spa arrays (array-only contigs)
    "$PY" - "$TOY" <<'PYEOF'
import sys, os
toy = sys.argv[1]
seqs, name = {}, None
for line in open(os.path.join(toy, 'genomes', 'spa.fna')):
    line = line.strip()
    if line.startswith('>'): name = line[1:].split()[0]; seqs[name] = []
    elif name: seqs[name].append(line)
seqs = {k: ''.join(v) for k, v in seqs.items()}
arrays = []
for line in open(os.path.join(toy, 'truth.tsv')):
    p = line.rstrip('\n').split('\t')
    if p[0] == 'spa' and p[4].startswith('array') and int(p[6]) >= 20: arrays.append(p)
with open(os.path.join(toy, 'genomes', 'spa2.fna'), 'w') as f:
    for k, s in seqs.items():
        f.write('>h2_%s\n' % k)
        for i in range(0, len(s), 80): f.write(s[i:i + 80] + '\n')
    for j, p in enumerate(arrays[:3]):
        s = seqs[p[1]][int(p[2]):int(p[3])]
        f.write('>h2_arr%d\n' % (j + 1))
        for i in range(0, len(s), 80): f.write(s[i:i + 80] + '\n')
        print('array-only contig h2_arr%d = spa %s:%s-%s %s %s copies' % (j + 1, p[1], p[2], p[3], p[5], p[6]))
PYEOF
    GEN3="  spa2: $TOY/genomes/spa2.fna"
fi

mkcfg() {   # mkcfg OUT_DIR [extra yaml...]
    local out=$1; shift
    cat <<EOF
run_name: e2e_$V
out_dir: $out
genomes:
  spa: $TOY/genomes/spa.fna
  spb: $TOY/genomes/spb.fna
EOF
    [ -n "$GEN3" ] && printf '%s\n' "$GEN3"
    cat <<EOF
families:
  CLsat1: {consensus: $CONS/CLsat1.cons.fa, mono: 146}
  CLsat2: {consensus: $CONS/CLsat2.cons.fa, mono: 145}
  CLsat3: {consensus: $CONS/CLsat3.cons.fa, mono: 146}
  CLsat4: {consensus: $CONS/CLsat4.cons.fa, mono: 116}
references: {rotated_consensus: $TOY/refs.rot}
family_assignment: $TAB
subfamily:
  consensuses: $TOY/subx
  members: {CLsat1: [CLsat1.1, CLsat1.2], CLsat2: [CLsat2], CLsat3: [CLsat3.1, CLsat3.2], CLsat4: [CLsat4]}
expansion: {jobs: $JOBS}
discovery: {monomerize: true, strand: plus}
EOF
    [ $# -gt 0 ] && printf '%s\n' "$@"
    return 0
}

# ---- prep: the curated table (and seeds) from a separate Python run
printf '0\tCLsat1\n1\tCLsat2\n2\tCLsat3\n3\tCLsat4\n' > "$TAB"
mkcfg "$W/prep" > "$W/prep.yaml"
drive py "$W/out/prep" run "$W/prep.yaml" --to 08_cross
[ "$(cat "$W/out/prep.rc")" = 0 ] || { echo "prep run failed"; tail -n 20 "$W/out/prep.stdout" "$W/out/prep.stderr"; exit 1; }
"$GAWK" -F'\t' '($1=="S"||$1=="H") && $9 ~ /^CLsat[1-4]$/ {print $2"\t"$9}' "$W/prep/cross/cross_species.uc" | LC_ALL=C sort -k1,1n > "$TAB"
echo "cluster_family.tsv:"; sed 's/^/  /' "$TAB"
EXTRA=()
case $V in
    seeds)
        drive py "$W/out/prep10" run "$W/prep.yaml" --to 10_seeds
        [ "$(cat "$W/out/prep10.rc")" = 0 ] || { echo "prep to 10 failed"; tail -n 20 "$W/out/prep10.stdout"; exit 1; }
        mkdir -p "$W/seeds"; cp "$W"/prep/seeds/CLsat?.bed "$W/seeds/"
        EXTRA=("inputs:" "  seeds: {CLsat1: $W/seeds/CLsat1.bed, CLsat2: $W/seeds/CLsat2.bed, CLsat3: $W/seeds/CLsat3.bed, CLsat4: $W/seeds/CLsat4.bed}")
        wc -l "$W"/seeds/*.bed ;;
    manual)
        "$GAWK" -F'\t' -v OFS='\t' 'NR == 1 { print "track name=manual"; print "# hand-added" }
            $1 == "spa" && $5 ~ /^array/ && !a { a = 1; print $2, $3 + 50, $4 - 50, "spa|manual1"; next }
            $1 == "spb" && $5 ~ /^array/ && !b { b = 1; print $2, $3, $4, "spb|manual2", "0", "+", $6 }' "$TOY/truth.tsv" > "$W/manual.bed"
        cat "$W/manual.bed"
        EXTRA=("manual_loci: $W/manual.bed") ;;
    hapdup) EXTRA=("hapdup: {pairs: [[spa, spa2]]}") ;;
esac
rm -rf -- "$W/prep"
mkcfg "$RUN" "${EXTRA[@]}" > "$W/run.yaml"

# ---- one full run per driver (python twice for the control in `default`)
SIDES=(py sh); [ "$V" = default ] && SIDES=(py sh py2)
for s in "${SIDES[@]}"; do
    rm -rf -- "$RUN"; rm -f -- "$TOY"/genomes/*.fai
    t0=$(date +%s)
    drive "${s%2}" "$W/out/$s" run "$W/run.yaml"
    echo "  $s run: rc=$(cat "$W/out/$s.rc") $(( $(date +%s) - t0 )) s"
    rm -rf -- "$W/run.$s"; [ -d "$RUN" ] && mv "$RUN" "$W/run.$s"
done
echo "--- py vs sh"
compare_drive "$V" "$W/out/py" "$W/out/sh"
compare_runs "$V" "$W/run.py" "$W/run.sh" "$W/out/cmp"
echo "  manifest summary (sh): $(jq -c '{centroids_n, seed_counts, master_loci, segmented_loci, manual_loci, catalogue_loci, hapdup_removed, ortholog_pairs, candidates}' "$W/run.sh/manifest.json")"
if [ "$V" = default ]; then
    echo "--- control: py vs py (same driver twice)"
    compare_drive "$V-control" "$W/out/py" "$W/out/py2"
    compare_runs "$V-control" "$W/run.py" "$W/run.py2" "$W/out/cmp_control"
fi

# ---- cross reruns
echo "--- cross reruns"
cross_rerun "$V" sh "$W/run.py" "$RUN" "$W/out/x_sh_on_py" run "$W/run.yaml"
cross_rerun "$V" py "$W/run.sh" "$RUN" "$W/out/x_py_on_sh" run "$W/run.yaml"

# ---- viewer, each driver on its own run dir
echo "--- viewer"
for s in py sh; do
    rm -rf -- "$RUN"; mv "$W/run.$s" "$RUN"
    drive "$s" "$W/out/viewer_$s" viewer "$RUN"
    mv "$RUN" "$W/run.$s"
done
compare_drive "$V viewer" "$W/out/viewer_py" "$W/out/viewer_sh"
if [ -d "$W/run.py/viewer" ] && bash "$ROOT/tests/port/diff_viewer.sh" "$W/run.py/viewer" "$W/run.sh/viewer" > "$W/out/viewer.diff" 2>&1; then ok "$V viewer dir: $(tail -n 1 "$W/out/viewer.diff")"
else fail "$V viewer dir" "$W/out/viewer.diff"; fi
if [ -d "$W/run.py/viewer_flankscan" ]; then
    if diff -r "$W/run.py/viewer_flankscan" "$W/run.sh/viewer_flankscan" > "$W/out/viewer_flankscan.diff"; then ok "$V viewer_flankscan identical ($(find "$W/run.py/viewer_flankscan" -type f | wc -l) files)"
    else fail "$V viewer_flankscan" "$W/out/viewer_flankscan.diff"; fi
fi
[ "$V" = hapdup ] && echo "  hapdup_removed py=$(jq .hapdup_removed "$W/run.py/manifest.json") sh=$(jq .hapdup_removed "$W/run.sh/manifest.json")"
echo "E2E $V: passed $NPASS, failed $NFAIL, log-only differences $NLOG  ($(date +%FT%T))"
[ "$NFAIL" = 0 ]
