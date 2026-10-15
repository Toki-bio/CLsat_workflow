#!/usr/bin/env bash
# Second test case (work-plan step 9): the Amazon molly, Poecilia formosa -- the two 2021
# parental-haplotype assemblies (lat = GCA_013036135.2, P. latipinna-derived; mex =
# GCA_013036085.2, P. mexicana-derived), with the TRF output and the Step 3-8 families that the
# CLsat-revision session produced on therioserver (HANDOFF.md, "Second test organism").
#
# Families for the test = three pooled Step-8 families of that session, consensus = EMBOSS cons of
# the family's Step-7 rotated alignment, mono = consensus length. The curated cluster->family table
# is derived automatically here (the Step-8 clusters that contain each reference record), which is
# what a curator would do for an unambiguous reference; for a new genome the curator inspects
# cross/cluster_report.tsv instead.
#
# Usage: run_molly.sh <repo_dir> <molly_trf_dir> <work_dir> [jobs]
#   molly_trf_dir  /home/toki/clsat_revision/molly_trf (genomes/, wf37/dat/{lat,mex}.dat, wf37/step8/s7/)
set -euo pipefail
REPO="$(cd "${1:?repo}" && pwd)"; M="${2:?molly_trf_dir}"; W="${3:?work_dir}"; JOBS="${4:-16}"
mkdir -p "$W/ref"; W="$(cd "$W" && pwd)"
for t in cons samtools satpipe; do command -v $t >/dev/null || { echo "$t not on PATH" >&2; exit 2; }; done

# genomes (indexed next to the originals' copies in the work dir so nothing is written into the source tree)
for pair in "lat:GCA_013036135.2" "mex:GCA_013036085.2"; do
  code="${pair%%:*}"; acc="${pair#*:}"
  [ -s "$W/ref/$code.fna" ] || ln -sf "$M/genomes/$acc.fna" "$W/ref/$code.fna"
  [ -s "$W/ref/$code.fna.fai" ] || samtools faidx "$W/ref/$code.fna"
done
# consensus per family from the Step-7 rotated alignment (degapped EMBOSS majority consensus)
declare -A FAM=([MolSat1]=molly_fam011_n29_lat14_mex15_unit185bp [MolSat2]=molly_fam004_n4_lat4_mex0_unit236bp [MolSat3]=molly_fam000_n3_lat0_mex3_unit253bp)
: > "$W/ref/families_note.txt"
for fam in MolSat1 MolSat2 MolSat3; do
  src="$M/wf37/step8/s7/${FAM[$fam]}_rotated/${FAM[$fam]}.rotated.fa"
  [ -s "$src" ] || { echo "missing alignment $src" >&2; exit 2; }
  cons -sequence "$src" -outseq "$W/ref/$fam.cons.raw.fa" -name "$fam" -auto 2>/dev/null
  awk -v n="$fam" '/^>/{next}{gsub(/[-nN]/,""); s=s toupper($0)} END{print ">"n; print s}' "$W/ref/$fam.cons.raw.fa" > "$W/ref/$fam.cons.fa"
  echo "$fam	${FAM[$fam]}	$(awk '!/^>/{print length($0)}' "$W/ref/$fam.cons.fa") bp" >> "$W/ref/families_note.txt"
done
cat "$W/ref/families_note.txt"
mono(){ awk '!/^>/{print length($0)}' "$W/ref/$1.cons.fa"; }

cat > "$W/molly.yaml" <<EOF
# Amazon molly test case (work-plan step 9); see tests/molly/run_molly.sh
run_name: molly_2021_halves
out_dir: $W/run
genomes:
  lat: $W/ref/lat.fna
  mex: $W/ref/mex.fna
inputs:
  trf_dat: {lat: $M/wf37/dat/lat.dat, mex: $M/wf37/dat/mex.dat}
families:
  MolSat1: {consensus: $W/ref/MolSat1.cons.fa, mono: $(mono MolSat1), colour: "#60a5fa"}
  MolSat2: {consensus: $W/ref/MolSat2.cons.fa, mono: $(mono MolSat2), colour: "#34d399"}
  MolSat3: {consensus: $W/ref/MolSat3.cons.fa, mono: $(mono MolSat3), colour: "#f472b6"}
family_assignment: $W/ref/cluster_family.tsv
expansion: {jobs: $JOBS}
finalisation: {large_scaffold_bp: 80000000}
tools: {bin_dir: $(dirname "$(command -v samtools)")}
viewer:
  title: "Satellite loci in the two Amazon molly haplotype assemblies (2021)"
  species_names: {lat: "P. formosa, latipinna-derived half (GCA_013036135.2)", mex: "P. formosa, mexicana-derived half (GCA_013036085.2)"}
EOF
printf '0\tMolSat1\n' > "$W/ref/cluster_family.tsv"     # placeholder until stage 08 has run
satpipe check "$W/molly.yaml" | tail -1
satpipe run "$W/molly.yaml" --to 08_cross
# curated table: every Step-8 cluster that holds a reference record
awk -F'\t' '($1=="S"||$1=="H") && $9 ~ /^MolSat[0-9]+$/ {print $2"\t"$9}' "$W/run/cross/cross_species.uc" | sort -k1,1n > "$W/ref/cluster_family.tsv"
echo "cluster_family.tsv:"; sed 's/^/  /' "$W/ref/cluster_family.tsv"
satpipe run "$W/molly.yaml"
satpipe viewer "$W/run"
echo "=== manifest"; python3 -c "import json;m=json.load(open('$W/run/manifest.json'));print({k:m.get(k) for k in ('centroids_n','seed_counts','master_loci','segmented_loci','catalogue_loci','ortholog_pairs')})"
echo "MOLLY DONE: $W/run/viewer"
