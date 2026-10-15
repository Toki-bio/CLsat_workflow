#!/usr/bin/env bash
# Step 8 - cross-species clustering. Pool all Step-7 rotated per-species cluster sequences,
# prepend the four CLsat reference consensuses (rotated to the same phase convention:
# CLsat_cons.rot), degap, and cluster at 90% identity. Prepending the references makes them
# win centroid election for their families: --cluster_fast sorts by decreasing length and
# processes equal lengths in input order.
#
# Usage: Step8_cross_cluster.sh <clusters_by_file_dir> <CLsat_cons.rot> [out_dir]
# Output (in out_dir): all7_pooled.degap.fa, cross_species_centroids.fa, cross_species.uc,
#                      cross_clusters/clust_<N>  (N = vsearch cluster number, used in cluster_family.tsv)
set -euo pipefail
CB="${1:?clusters_by_file dir}"; REF="${2:?CLsat_cons.rot}"; OUT="${3:-.}"
mkdir -p "$OUT/cross_clusters"
cat "$CB"/*/*_rotated/*.rotated.fa > "$OUT/all7_pooled.rot.fa"
cat "$REF" "$OUT/all7_pooled.rot.fa" > "$OUT/all7_pooled.rot.fas"
awk '/^>/{print} !/^>/{gsub(/-/,""); if(length($0)>0) print}' "$OUT/all7_pooled.rot.fas" > "$OUT/all7_pooled.degap.fa"
vsearch --cluster_fast "$OUT/all7_pooled.degap.fa" --id 0.90 \
  --centroids "$OUT/cross_species_centroids.fa" --clusters "$OUT/cross_clusters/clust_" \
  --uc "$OUT/cross_species.uc" --quiet
echo "Step 8: $(grep -c '^>' "$OUT/cross_species_centroids.fa") centroids (published run: 331)"
