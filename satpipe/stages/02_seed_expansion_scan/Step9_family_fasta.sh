#!/usr/bin/env bash
# Step 9 - build one FASTA per family from the curated cluster->family table.
#
# Usage: Step9_family_fasta.sh <cross_clusters_dir> [cluster_family.tsv] [out_dir] [rotated_references.fa]
#   cluster_family.tsv   curated: cluster<TAB>family[<TAB>note]; '#' lines ignored
#                        [default: cluster_family.tsv next to this script]
#   out_dir              [default: .]
#   rotated_references   the reference FASTA prepended in Step 8; its record names are
#                        removed from the member lists so that only genomic members become
#                        seeds. Alternatively set REF_NAMES="name1 name2 ...". With neither,
#                        nothing is removed and a warning is printed.
# Output: <out_dir>/<family>.members.fa for every family in the table
#         -> Step 10: make_bed_from_clsats_headers.sh
set -euo pipefail
CC="${1:?cross_clusters dir}"; TAB="${2:-$(dirname "$0")/cluster_family.tsv}"; OUT="${3:-.}"; REF="${4:-}"
[[ -s "$TAB" ]] || { echo "ERROR: cluster->family table missing/empty: $TAB" >&2; exit 1; }
mkdir -p "$OUT"

# reference record names to drop
REF_NAMES="${REF_NAMES:-}"
if [[ -n "$REF" ]]; then
  [[ -s "$REF" ]] || { echo "ERROR: reference FASTA missing/empty: $REF" >&2; exit 1; }
  REF_NAMES="$(awk '/^>/{sub(/^>/,""); print $1}' "$REF" | tr '\n' ' ')"
fi
[[ -n "$REF_NAMES" ]] || echo "WARNING: no reference names given; reference records (if any) stay in the member lists" >&2

# families, in order of first appearance in the table
mapfile -t FAMS < <(grep -v '^#' "$TAB" | awk -F'\t' 'NF>=2 && !seen[$2]++{print $2}')
(( ${#FAMS[@]} > 0 )) || { echo "ERROR: no cluster->family rows in $TAB" >&2; exit 1; }
for fam in "${FAMS[@]}"; do : > "$OUT/$fam.members.fa"; done

grep -v '^#' "$TAB" | awk -F'\t' 'NF>=2{print $1"\t"$2}' | while IFS=$'\t' read -r clu fam; do
  f="$CC/clust_$clu"
  [[ -s "$f" ]] || { echo "ERROR: $f missing" >&2; exit 1; }
  # drop the reference consensus records themselves so only genomic members become seeds
  awk -v names="$REF_NAMES" 'BEGIN{n=split(names,a," "); for(i=1;i<=n;i++) drop[a[i]]=1}
       /^>/{nm=substr($1,2); keep=!(nm in drop)} keep' "$f" >> "$OUT/$fam.members.fa"
done
for fam in "${FAMS[@]}"; do
  echo "$fam: $(grep -c '^>' "$OUT/$fam.members.fa" || true) members"
done
