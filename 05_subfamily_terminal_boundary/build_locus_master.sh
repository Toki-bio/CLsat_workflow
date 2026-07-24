#!/bin/bash
set -euo pipefail
DIR=/staging/Sergei/Lacertidae/Darevskia/Satellites/manuscript/CLsat_full_clusters_6genome
OUT=$DIR/out
WORK=$DIR/flank_analysis
mkdir -p $WORK

# 1. scaffold -> species
awk -F'\t' '{split($4,a,"|"); print $1"\t"a[1]}' $OUT/MASTER.merged.bed | sort -u > $WORK/scf2sp.tsv

# 2. locus -> subtype
awk -F'\t' '{print $1"\t"$2}' $OUT/MASTER.final_assignment.tsv > $WORK/locus2subtype.tsv

# 3. MONO table
printf 'CLsat1\t146\nCLsat2\t145\nCLsat3\t146\nCLsat4\t116\n' > $WORK/mono.tsv

# 4a. segmented bed + species
awk -F'\t' 'FNR==NR{sp[$1]=$2;next}{print $4"\t"$1"\t"$2"\t"$3"\t"(($1 in sp)?sp[$1]:"NA")}' \
  $WORK/scf2sp.tsv $OUT/MASTER.loci.segmented.bed > $WORK/tmp1.tsv

# 4b. + subtype  (use subt[] - 'sub' is a reserved awk function)
awk -F'\t' 'FNR==NR{subt[$1]=$2;next}{print $0"\t"(($1 in subt)?subt[$1]:"NA")}' \
  $WORK/locus2subtype.tsv $WORK/tmp1.tsv > $WORK/tmp2.tsv

# 4c. + MONO
awk -F'\t' 'FNR==NR{mono[$1]=$2;next}{print $0"\t"(($6 in mono)?mono[$6]:0)}' \
  $WORK/mono.tsv $WORK/tmp2.tsv > $WORK/locus_master.tsv

echo "=== locus_master.tsv (first 5) ==="
head -5 $WORK/locus_master.tsv
echo "total: $(wc -l < $WORK/locus_master.tsv)"
echo "NA/zero rows: $(awk -F'\t' '$5=="NA"||$6=="NA"||$7==0' $WORK/locus_master.tsv | wc -l)"
echo "cols in row 1: $(head -1 $WORK/locus_master.tsv | awk -F'\t' '{print NF}')"
echo "subtype distribution:"
awk -F'\t' '{print $6}' $WORK/locus_master.tsv | sort | uniq -c
