#!/usr/bin/env bash
# Run the whole driver (`satpipe run`) on the synthetic project and, when a run_toy_equivalence.sh
# work directory is given, compare the driver's outputs with that harness's NEW-chain outputs.
#
# Usage: run_driver_toy.sh <repo_dir> <work_dir> [equivalence_work_dir]
# Needs: satpipe installed (pip install -e <repo_dir>) and the tools on PATH.
set -euo pipefail
REPO="$(cd "${1:?repo_dir}" && pwd)"; W="${2:?work_dir}"; EQ="${3:-}"
mkdir -p "$W"; W="$(cd "$W" && pwd)"
export LC_ALL=C
TOY="$W/toy"
python3 "$REPO/tests/toy/make_toy.py" "$TOY" "$REPO/satpipe/stages/02_seed_expansion_scan/consensus" --seed "${SEED:-7}" --size "${SIZE:-1200000}"
# curated table: the harness derives it from the Step-8 .uc; here we take it from the equivalence run if present,
# else run the discovery stages first and derive it the same way
TAB="$W/cluster_family.tsv"
cat > "$W/toy.yaml" <<EOF
run_name: toy
out_dir: $W/run
genomes:
  spa: $TOY/genomes/spa.fna
  spb: $TOY/genomes/spb.fna
families:
  CLsat1: {consensus: $REPO/satpipe/stages/02_seed_expansion_scan/consensus/CLsat1.cons.fa, mono: 146}
  CLsat2: {consensus: $REPO/satpipe/stages/02_seed_expansion_scan/consensus/CLsat2.cons.fa, mono: 145}
  CLsat3: {consensus: $REPO/satpipe/stages/02_seed_expansion_scan/consensus/CLsat3.cons.fa, mono: 146}
  CLsat4: {consensus: $REPO/satpipe/stages/02_seed_expansion_scan/consensus/CLsat4.cons.fa, mono: 116}
references: {rotated_consensus: $TOY/refs.rot}
family_assignment: $TAB
subfamily:
  consensuses: $TOY/subx
  members: {CLsat1: [CLsat1.1, CLsat1.2], CLsat2: [CLsat2], CLsat3: [CLsat3.1, CLsat3.2], CLsat4: [CLsat4]}
expansion: {jobs: ${JOBS:-8}}
finalisation: {large_scaffold_bp: ${LARGE:-50000000}}
EOF
if [ -n "$EQ" ] && [ -s "$EQ/cluster_family.tsv" ]; then
  cp "$EQ/cluster_family.tsv" "$TAB"
else
  printf '0\tCLsat1\n1\tCLsat2\n2\tCLsat3\n3\tCLsat4\n' > "$TAB"   # placeholder (one row per family) so `check` passes; replaced after stage 08
  satpipe run "$W/toy.yaml" --to 08_cross
  awk -F'\t' '($1=="S"||$1=="H") && $9 ~ /^CLsat[1-4]$/ {print $2"\t"$9}' "$W/run/cross/cross_species.uc" | sort -k1,1n > "$TAB"
fi
echo "cluster_family.tsv:"; sed 's/^/  /' "$TAB"
satpipe run "$W/toy.yaml"
echo "--- second invocation must skip everything"
satpipe run "$W/toy.yaml" | grep -c 'skipped'
[ -z "$EQ" ] && { echo "DRIVER TOY: ran (no comparison requested)"; exit 0; }

echo "### compare driver run with the equivalence harness NEW chain in $EQ"
fail=0
cmpf(){ if cmp -s "$2" "$3"; then echo "  same   $1"; else echo "  DIFF   $1"; fail=1; fi; }
R="$W/run"; N="$EQ/newrun"
for sp in spa spb; do cmpf "motifs $sp" "$N/dat/motifs_output/$sp.fa" "$R/motifs/$sp.fa"; done
cmpf "centroids" "$N/dat/motifs_output/clean/cross_species_centroids.fa" "$R/cross/cross_species_centroids.fa"
for fam in CLsat1 CLsat2 CLsat3 CLsat4; do
  cmpf "seeds $fam" "$N/seeds/$fam.bed" "$R/seeds/$fam.bed"
  cmpf "expansion $fam" "$N/out_$fam.expanded_merged.bed" "$R/expand/out_$fam.expanded_merged.bed"
done
cmpf "MASTER" "$N/master/MASTER.merged.bed" "$R/master/MASTER.merged.bed"
cmpf "segmented" "$N/arb/MASTER.loci.segmented.bed" "$R/arbitration/MASTER.loci.segmented.bed"
cmpf "locus_master" "$N/fa/locus_master.tsv" "$R/regions/locus_master.tsv"
cmpf "manifest" "$N/fa/manifest.tsv" "$R/regions/manifest.tsv"
cmpf "terminal" "$N/fa/terminal_monomers.tsv" "$R/regions/terminal_monomers.tsv"
cmpf "boundary phase" "$N/fa/locus_boundary_phase.tsv" "$R/regions/locus_boundary_phase.tsv"
cmpf "flank_class" "$N/ortho/flank_class.tsv" "$R/ortho/flank_class.tsv"
python3 "$REPO/tests/toy/compare_noise.py" assignment "$N/arb/MASTER.final_assignment.tsv" "$R/arbitration/MASTER.final_assignment.tsv" && echo "  same~  final_assignment" || { echo "  DIFF   final_assignment"; fail=1; }
python3 "$REPO/tests/toy/compare_noise.py" subfamily "$N/fa/locus_subfamily.tsv" "$R/regions/locus_subfamily.tsv" && echo "  same~  locus_subfamily" || { echo "  DIFF   locus_subfamily"; fail=1; }
echo "manifest summary:"; python3 -c "import json;m=json.load(open('$R/manifest.json'));print({k:m.get(k) for k in ('centroids_n','seed_counts','master_loci','segmented_loci','catalogue_loci','ortholog_pairs','refs_rot')})"
[ "$fail" = 0 ] && echo "DRIVER TOY: PASS" || echo "DRIVER TOY: FAIL"
exit $fail
