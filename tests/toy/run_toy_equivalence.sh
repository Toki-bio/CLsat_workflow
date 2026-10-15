#!/usr/bin/env bash
# End-to-end equivalence test: run the ORIGINAL stage scripts (a given git commit, paths patched
# by sed the way the published re-runs did) and the PARAMETERISED scripts (working tree) on the
# same synthetic project, and compare every intermediate file.
#
# Usage: run_toy_equivalence.sh <repo_dir> <work_dir> [old_commit]
#   repo_dir    a clone of CLsat_workflow (working tree = the scripts under test)
#   work_dir    scratch directory (created); the toy project goes to <work_dir>/toy
#   old_commit  commit holding the pre-parameterisation scripts   [4d70b56]
# Needs on PATH: trf vsearch mafft ssearch36 gawk parallel samtools bedtools python3 bc
# Exit status 0 = every compared file identical.
set -euo pipefail
REPO="$(cd "${1:?repo_dir}" && pwd)"; W="${2:?work_dir}"; OLDC="${3:-4d70b56}"
mkdir -p "$W"; W="$(cd "$W" && pwd)"
export LC_ALL=C
NEW="$REPO/satpipe/stages"
TOY="$W/toy"
FAMS="CLsat1 CLsat2 CLsat3 CLsat4"
declare -A MONO=([CLsat1]=146 [CLsat2]=145 [CLsat3]=146 [CLsat4]=116)
JOBS="${JOBS:-8}"
fail=0
# the original scripts call `awk` and need gawk features (asorti); make sure awk IS gawk here
mkdir -p "$W/bin"; ln -sf "$(command -v gawk)" "$W/bin/awk"; export PATH="$W/bin:$PATH"
for t in trf vsearch mafft ssearch36 gawk parallel samtools bedtools python3; do command -v $t >/dev/null || { echo "$t not on PATH" >&2; exit 2; }; done
cmpf(){ # cmpf <label> <old> <new>
  if cmp -s "$2" "$3"; then echo "  same   $1"; else echo "  DIFF   $1 ($2 vs $3)"; fail=1; fi
}
cmpjson(){ python3 - "$2" "$3" <<'PY' && echo "  same   $1" || { echo "  DIFF   $1"; fail=1; }
import json,sys
a=json.load(open(sys.argv[1])); b=json.load(open(sys.argv[2]))
sys.exit(0 if a==b else 1)
PY
}

# ---------------------------------------------------------------- 0. toy project + TRF
echo "### toy project"
python3 "$REPO/tests/toy/make_toy.py" "$TOY" "$NEW/02_seed_expansion_scan/consensus" --seed "${SEED:-7}" --size "${SIZE:-1200000}"
for sp in spa spb; do samtools faidx "$TOY/genomes/$sp.fna"; done
mkdir -p "$TOY/dat"; ( cd "$TOY/dat" && for sp in spa spb; do
  [ -s "$sp.dat" ] || { trf "$TOY/genomes/$sp.fna" 2 5 7 80 10 50 2000 -h > /dev/null 2>&1 || true; mv "$sp.fna.2.5.7.80.10.50.2000.dat" "$sp.dat"; }
done )
# make the four consensus files and genomes.tsv visible under one dir "D" (the old scripts' layout)
D="$W/D"; mkdir -p "$D/out"
for fam in $FAMS; do cp "$NEW/02_seed_expansion_scan/consensus/$fam.cons.fa" "$D/"; done
cp "$TOY/genomes.tsv" "$D/genomes.tsv"
cp "$TOY/refs.rot" "$D/CLsat_cons.rot"

# ---------------------------------------------------------------- old scripts, path-patched
OLD="$W/old_scripts"; mkdir -p "$OLD"
getold(){ git -C "$REPO" show "$OLDC:$1" > "$OLD/$(basename "$1")"; }
for p in 01_discovery_clustering/Step3.sh 01_discovery_clustering/Step4.sh 01_discovery_clustering/Step5.sh \
         01_discovery_clustering/Step6.sh 01_discovery_clustering/Step7.sh 02_seed_expansion_scan/Step8_cross_cluster.sh \
         02_seed_expansion_scan/Step9_family_fasta.sh 02_seed_expansion_scan/make_bed_from_clsats_headers.sh \
         02_seed_expansion_scan/clisat_expand_then_scan_merged.sh 03_master_arbitration/build_master.sh \
         03_master_arbitration/arbitration_patched.sh 05_subfamily_terminal_boundary/build_locus_master.sh \
         05_subfamily_terminal_boundary/terminal_monomers.sh 05_subfamily_terminal_boundary/boundary_phase.sh \
         05_subfamily_terminal_boundary/assign_subfamily.sh 05_subfamily_terminal_boundary/build_json.py \
         05_subfamily_terminal_boundary/build_subfam_term_json.py 06_catalogue_finalize/subfamily_rerun/run_extract.sh \
         04_flank_orthology/rebuild_ortho_segmented.py; do getold "$p"; done
SD="$W/oldrun/seeds"
sed -i -e "s#^BASE=.*#BASE=\"$D\"#" -e "s#^OUT=\"\$BASE/out\"#OUT=\"$D/out\"#" \
       -e "s#\$OUT/MASTER_CLsat\([1-4]\).bits.tsv#$W/oldrun/scan/MASTER_CLsat\1.bits.tsv#" \
       -e "s#\$BASE/CLsat\([1-4]\).bed#$SD/CLsat\1.bed#" "$OLD/arbitration_patched.sh"
for f in build_locus_master.sh terminal_monomers.sh boundary_phase.sh assign_subfamily.sh; do
  sed -i -e "s#^DIR=.*#DIR=$D#" -e '/^source \/staging/d' -e '/^conda activate/d' "$OLD/$f"; done
sed -i -e "s#^D=.*#D=$D#" -e "s#^WORK=.*#WORK=$D/flank_analysis#" -e '/^source \/staging/d' -e '/^conda activate/d' "$OLD/run_extract.sh"
sed -i -e "s#^WORK = .*#WORK = \"$D/flank_analysis\"#" "$OLD/build_json.py"
sed -i -e "s#^W = .*#W = \"$D/flank_analysis\"#" "$OLD/build_subfam_term_json.py"
BIN="$(dirname "$(command -v ssearch36)")"
sed -i -e "s#^BIN = .*#BIN = '$BIN'#" -e "s#^SD = (.*#SD = ('$D'#" -e "s#^      'cross_clusters/selected/subclusters/seeds')#      '')#" \
       -e "s#^CAT = .*#CAT = '$D/flank_analysis/segmented_catalog.tsv'#" -e "s#^OUT = .*#OUT = '$W/oldrun/ortho'#" "$OLD/rebuild_ortho_segmented.py"
grep -n '/staging' "$OLD"/* && { echo "unpatched /staging path in old scripts"; exit 2; } || true

# ---------------------------------------------------------------- 1. OLD chain
echo "### OLD chain ($OLDC)"
O="$W/oldrun"; mkdir -p "$O/dat" "$O/scan" "$SD"
cp "$TOY/dat/"*.dat "$O/dat/"
( cd "$O/dat" && bash "$OLD/Step3.sh" > step3.log && cd motifs_output && for sp in spa spb; do bash "$OLD/Step4.sh" $sp.fa > /dev/null; done \
  && bash "$OLD/Step5.sh" . > /dev/null && cd clean && bash "$OLD/Step6.sh" > /dev/null \
  && for d in clusters_by_file/*/; do ( cd "$d" && bash "$OLD/Step7.sh" > step7.log 2>&1 ); done \
  && bash "$OLD/Step8_cross_cluster.sh" clusters_by_file "$D/CLsat_cons.rot" . > step8.log )
# curated table: the cluster whose centroid is the reference record
# (S = the reference is the centroid itself; H = it is a member of another centroid's cluster)
awk -F'\t' '($1=="S"||$1=="H") && $9 ~ /^CLsat[1-4]$/ {print $2"\t"$9}' "$O/dat/motifs_output/clean/cross_species.uc" | sort -k1,1n > "$W/cluster_family.tsv"
echo "  cluster_family.tsv:"; sed 's/^/    /' "$W/cluster_family.tsv"
( cd "$O/dat/motifs_output/clean" && bash "$OLD/Step9_family_fasta.sh" cross_clusters "$W/cluster_family.tsv" "$SD" > step9.log )
for fam in $FAMS; do bash "$OLD/make_bed_from_clsats_headers.sh" "$SD/$fam.members.fa" > "$SD/$fam.bed"; done
( cd "$O" && for fam in $FAMS; do
    bash "$OLD/clisat_expand_then_scan_merged.sh" "$SD/$fam.bed" "$D/genomes.tsv" "$D/$fam.cons.fa" "${MONO[$fam]}" "$O/out_$fam" --jobs "$JOBS" > "$O/expand_$fam.log" 2>&1
  done && bash "$OLD/build_master.sh" > build_master.log )   # writes $O/out/MASTER.merged.bed
cp "$O/out/MASTER.merged.bed" "$D/out/MASTER.merged.bed"
for fam in $FAMS; do
  bash "$OLD/clisat_expand_then_scan_merged.sh" "$D/out/MASTER.merged.bed" "$D/genomes.tsv" "$D/$fam.cons.fa" "${MONO[$fam]}" "$O/scan/MASTER_$fam" --scan-only --jobs "$JOBS" > "$O/scan_$fam.log" 2>&1
done
bash "$OLD/arbitration_patched.sh" > "$O/arbitration.log"           # writes $D/out/*
bash "$OLD/build_locus_master.sh" > "$O/locus_master.log"           # $D/flank_analysis/locus_master.tsv
bash "$OLD/run_extract.sh" > "$O/extract.log"                        # manifest.tsv all_regions.fa
mkdir -p "$D/flank_analysis/subx"; cp "$TOY/subx/"*.fa "$D/flank_analysis/subx/"; cp "$D/CLsat1.cons.fa" "$D/flank_analysis/subfamily_consensus.fa"
bash "$OLD/assign_subfamily.sh" > "$O/assign.log"
bash "$OLD/terminal_monomers.sh" > "$O/terminal.log"
bash "$OLD/boundary_phase.sh" > "$O/bphase.log"
python3 "$OLD/build_json.py" > "$O/build_json.log"
python3 "$OLD/build_subfam_term_json.py" > "$O/build_subfam.log"
cut -f1-6 "$D/flank_analysis/locus_master.tsv" > "$D/flank_analysis/segmented_catalog.tsv"
python3 "$OLD/rebuild_ortho_segmented.py" > "$O/ortho.log"
mv "$D/flank_analysis" "$O/flank_analysis"; mv "$D/out" "$O/arb"; mkdir -p "$D/out"

# ---------------------------------------------------------------- 2. NEW chain
echo "### NEW chain (working tree)"
N="$W/newrun"; mkdir -p "$N/dat" "$N/scan" "$N/seeds" "$N/arb" "$N/fa"
cp "$TOY/dat/"*.dat "$N/dat/"
( cd "$N/dat" && bash "$NEW/01_discovery_clustering/Step3.sh" > step3.log && cd motifs_output && for sp in spa spb; do bash "$NEW/01_discovery_clustering/Step4.sh" $sp.fa > /dev/null; done \
  && bash "$NEW/01_discovery_clustering/Step5.sh" . > /dev/null && cd clean && bash "$NEW/01_discovery_clustering/Step6.sh" > /dev/null \
  && for d in clusters_by_file/*/; do ( cd "$d" && bash "$NEW/01_discovery_clustering/Step7.sh" > step7.log 2>&1 ); done \
  && bash "$NEW/02_seed_expansion_scan/Step8_cross_cluster.sh" clusters_by_file "$TOY/refs.rot" . > step8.log \
  && bash "$NEW/02_seed_expansion_scan/Step9_family_fasta.sh" cross_clusters "$W/cluster_family.tsv" "$N/seeds" "$TOY/refs.rot" > step9.log )
for fam in $FAMS; do bash "$NEW/02_seed_expansion_scan/make_bed_from_clsats_headers.sh" "$N/seeds/$fam.members.fa" > "$N/seeds/$fam.bed"; done
: > "$N/expansions.tsv"; : > "$N/arbitration_inputs.tsv"
for fam in $FAMS; do
  bash "$NEW/02_seed_expansion_scan/clisat_expand_then_scan_merged.sh" "$N/seeds/$fam.bed" "$TOY/genomes.tsv" "$D/$fam.cons.fa" "${MONO[$fam]}" "$N/out_$fam" --jobs "$JOBS" > "$N/expand_$fam.log" 2>&1
  printf '%s\t%s\n' "$fam" "$N/out_$fam.expanded_merged.bed" >> "$N/expansions.tsv"
done
bash "$NEW/03_master_arbitration/build_master.sh" "$N/expansions.tsv" "$N/master" > "$N/build_master.log"
for fam in $FAMS; do
  bash "$NEW/02_seed_expansion_scan/clisat_expand_then_scan_merged.sh" "$N/master/MASTER.merged.bed" "$TOY/genomes.tsv" "$D/$fam.cons.fa" "${MONO[$fam]}" "$N/scan/MASTER_$fam" --scan-only --jobs "$JOBS" > "$N/scan_$fam.log" 2>&1
  printf '%s\t%s\t%s\t%s\n' "$fam" "${MONO[$fam]}" "$N/scan/MASTER_$fam.bits.tsv" "$N/seeds/$fam.bed" >> "$N/arbitration_inputs.tsv"
done
bash "$NEW/03_master_arbitration/arbitration_patched.sh" "$N/arbitration_inputs.tsv" "$N/master/MASTER.merged.bed" "$N/arb" > "$N/arbitration.log"
bash "$NEW/05_subfamily_terminal_boundary/build_locus_master.sh" "$N/master/MASTER.merged.bed" "$N/arb" "$N/fa" "$TOY/families.tsv" > "$N/locus_master.log"
bash "$NEW/06_catalogue_finalize/subfamily_rerun/run_extract.sh" "$N/fa" "$TOY/genomes.tsv" > "$N/extract.log"
bash "$NEW/05_subfamily_terminal_boundary/assign_subfamily.sh" "$N/fa" "$TOY/subx" "$TOY/members.tsv" > "$N/assign.log"
bash "$NEW/05_subfamily_terminal_boundary/terminal_monomers.sh" "$N/fa" "$TOY/genomes.tsv" "$TOY/families.tsv" > "$N/terminal.log"
bash "$NEW/05_subfamily_terminal_boundary/boundary_phase.sh" "$N/fa" "$TOY/families.tsv" > "$N/bphase.log"
python3 "$NEW/05_subfamily_terminal_boundary/build_json.py" "$N/fa" > "$N/build_json.log"
python3 "$NEW/05_subfamily_terminal_boundary/build_subfam_term_json.py" "$N/fa" "$TOY/families.tsv" > "$N/build_subfam.log"
python3 "$NEW/04_flank_orthology/rebuild_ortho_segmented.py" --loci "$N/fa/locus_master.tsv" --genomes "$TOY/genomes.tsv" \
        --families "$TOY/families.tsv" --out "$N/ortho" > "$N/ortho.log"

# ---------------------------------------------------------------- 3. compare
echo "### compare"
for sp in spa spb; do
  cmpf "step3 $sp"  "$O/dat/motifs_output/$sp.fa" "$N/dat/motifs_output/$sp.fa"
  cmpf "step4 $sp"  "$O/dat/motifs_output/$sp.fa.clean.fas" "$N/dat/motifs_output/$sp.fa.clean.fas"
  cmpf "step5 $sp"  "$O/dat/motifs_output/clean/$sp.fa.clean.fas.rot" "$N/dat/motifs_output/clean/$sp.fa.clean.fas.rot"
  cmpf "step6 $sp (cluster count)" <(ls "$O/dat/motifs_output/clean/clusters_by_file/$sp.fa.clean.fas" | sort) <(ls "$N/dat/motifs_output/clean/clusters_by_file/$sp.fa.clean.fas" | sort)
done
cmpf "step8 centroids" "$O/dat/motifs_output/clean/cross_species_centroids.fa" "$N/dat/motifs_output/clean/cross_species_centroids.fa"
for fam in $FAMS; do
  cmpf "step9 $fam members" "$SD/$fam.members.fa" "$N/seeds/$fam.members.fa"
  cmpf "step10 $fam seeds"  "$SD/$fam.bed" "$N/seeds/$fam.bed"
  cmpf "step11 $fam expansion" "$O/out_$fam.expanded_merged.bed" "$N/out_$fam.expanded_merged.bed"
  cmpf "step13 $fam bits" "$O/scan/MASTER_$fam.bits.tsv" "$N/scan/MASTER_$fam.bits.tsv"
done
cmpf "step12 MASTER" "$O/out/MASTER.merged.bed" "$N/master/MASTER.merged.bed"
for f in MASTER.loci.bed MASTER.loci.segmented.bed MASTER.locus_bits.tsv MASTER.final_assignment.tsv MASTER.locus_provenance.segmented.tsv; do
  cmpf "step14 $f" "$O/arb/$f" "$N/arb/$f"
done
for f in locus_master.tsv manifest.tsv all_regions.fa locus_subfamily.tsv terminal_monomers.tsv locus_boundary_phase.tsv; do
  cmpf "05 $f" "$O/flank_analysis/$f" "$N/fa/$f"
done
for f in subfamily.json terminal_monomers.json terminal_hist.json boundary_phase.json boundary_phase_hist.json; do
  cmpjson "05 $f" "$O/flank_analysis/$f" "$N/fa/$f"
done
cmpjson "ortho flanks.json" "$O/ortho/flanks_627.json" "$N/ortho/flanks.json"
cmpjson "ortho flank_orthologs.json" "$O/ortho/flank_orthologs_627.json" "$N/ortho/flank_orthologs.json"
cmpf "ortho pairs" "$O/ortho/ortho_pairs.tsv" "$N/ortho/ortho_pairs.tsv"
cmpf "ortho flank_class" "$O/ortho/flank_class.tsv" "$N/ortho/flank_class.tsv"
echo
echo "loci: $(wc -l < "$N/master/MASTER.merged.bed") MASTER, $(wc -l < "$N/arb/MASTER.loci.segmented.bed") segmented; ortholog pairs: $(( $(wc -l < "$N/ortho/ortho_pairs.tsv") - 1 ))"
if [ "$fail" = 0 ]; then echo "TOY EQUIVALENCE: PASS"; else echo "TOY EQUIVALENCE: FAIL"; fi
exit $fail
