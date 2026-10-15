#!/usr/bin/env bash
# Step 9 - build one FASTA per CLsat family from the curated cluster->family table.
# Usage: Step9_family_fasta.sh <cross_clusters_dir> [cluster_family.tsv] [out_dir]
# Output: <out_dir>/CLsat{1,2,3,4}.members.fa  -> Step 10: make_bed_from_clsats_headers.sh
set -euo pipefail
CC="${1:?cross_clusters dir}"; TAB="${2:-$(dirname "$0")/cluster_family.tsv}"; OUT="${3:-.}"
mkdir -p "$OUT"
for fam in CLsat1 CLsat2 CLsat3 CLsat4; do : > "$OUT/$fam.members.fa"; done
grep -v '^#' "$TAB" | awk -F'\t' 'NF>=2{print $1"\t"$2}' | while IFS=$'\t' read -r clu fam; do
  f="$CC/clust_$clu"
  [[ -s "$f" ]] || { echo "ERROR: $f missing" >&2; exit 1; }
  # drop the reference consensus itself (header starts with CLsat) so only genomic members become seeds
  awk '/^>/{keep=($0 !~ /^>CLsat/)} keep' "$f" >> "$OUT/$fam.members.fa"
done
for fam in CLsat1 CLsat2 CLsat3 CLsat4; do
  echo "$fam: $(grep -c '^>' "$OUT/$fam.members.fa") members"
done
echo "published seed counts: 750 / 288 / 101 / 44"
