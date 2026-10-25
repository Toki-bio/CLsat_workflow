#!/usr/bin/env bash
# Step 8 - cross-species clustering. Pool all Step-7 rotated per-species cluster sequences,
# prepend the family reference consensuses (rotated to the same phase convention; in the
# Darevskia run the file CLsat_cons.rot), degap, and cluster. Prepending the references
# makes them win centroid election for their families: --cluster_fast sorts by decreasing
# length and processes equal lengths in input order.
#
# Usage: Step8_cross_cluster.sh <clusters_by_file_dir> <rotated_references.fa> [out_dir]
# Environment (defaults = published Darevskia run):
#   CROSS_CLUSTER_ID  vsearch --id               [0.90]
#   STRAND            vsearch --strand           [plus]
# Output (in out_dir): all7_pooled.degap.fa, cross_species_centroids.fa, cross_species.uc,
#                      cross_clusters/clust_<N>  (N = vsearch cluster number, used in cluster_family.tsv)
set -euo pipefail
CB="${1:?clusters_by_file dir}"; REF="${2:?rotated reference consensuses FASTA}"; OUT="${3:-.}"
ID="${CROSS_CLUSTER_ID:-0.90}"
STRAND="${STRAND:-plus}"
command -v vsearch >/dev/null 2>&1 || { echo "Error: vsearch not found in PATH" >&2; exit 1; }
[[ -s "$REF" ]] || { echo "Error: reference FASTA missing/empty: $REF" >&2; exit 1; }
mkdir -p "$OUT/cross_clusters"
rm -f "$OUT"/cross_clusters/clust_*
# every per-genome cluster of 2+ loci must have its rotated alignment (Step 7); a missing one would drop its loci
# from the pool without a word (GLM audit #48)
miss=0
for f in "$CB"/*/*_cluster_*.fa; do
  [ -e "$f" ] || continue
  [ "$(grep -c '^>' "$f")" -gt 1 ] || continue
  b=$(basename "$f" .fa); r="$(dirname "$f")/${b}_rotated/${b}.rotated.fa"
  [ -s "$r" ] || { echo "Error: no rotated alignment for $f (expected $r)" >&2; miss=$((miss+1)); }
done
[ "$miss" -eq 0 ] || { echo "Error: $miss clusters without a rotated alignment; rerun stage 07" >&2; exit 1; }
cat "$CB"/*/*_rotated/*.rotated.fa > "$OUT/all7_pooled.rot.fa"
cat "$REF" "$OUT/all7_pooled.rot.fa" > "$OUT/all7_pooled.rot.fas"
awk '/^>/{print} !/^>/{gsub(/-/,""); if(length($0)>0) print}' "$OUT/all7_pooled.rot.fas" > "$OUT/all7_pooled.degap.fa"
vsearch --cluster_fast "$OUT/all7_pooled.degap.fa" --id "$ID" --strand "$STRAND" \
  --centroids "$OUT/cross_species_centroids.fa" --clusters "$OUT/cross_clusters/clust_" \
  --uc "$OUT/cross_species.uc" --quiet
echo "Step 8: $(grep -c '^>' "$OUT/cross_species_centroids.fa") centroids (id=$ID strand=$STRAND)"
