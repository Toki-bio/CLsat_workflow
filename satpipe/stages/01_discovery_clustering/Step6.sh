#!/usr/bin/env bash
set -euo pipefail
# Step 6 - per-species clustering of the rotated motifs with vsearch.
#
# Usage: Step6.sh [input_dir] [output_dir]
#   input_dir   where the .rot files live                     [default: .]
#   output_dir  receives one subfolder per .rot input         [default: clusters_by_file]
# Environment (defaults = published Darevskia run):
#   CLUSTER_ID  vsearch --id                                  [0.90]
#   STRAND      vsearch --strand (plus | both)                [plus]

INPUT_DIR="${1:-.}"
OUTPUT_DIR="${2:-clusters_by_file}"
ID="${CLUSTER_ID:-0.90}"
STRAND="${STRAND:-plus}"

mkdir -p "$OUTPUT_DIR"
shopt -s nullglob

# require vsearch
command -v vsearch >/dev/null 2>&1 || { echo "Error: vsearch not found in PATH"; exit 1; }

for infa in "$INPUT_DIR"/*.rot; do
  # skip if no .rot files
  [ -e "$infa" ] || continue

  base=$(basename "$infa" .rot)
  base_dir="$OUTPUT_DIR/$base"
  mkdir -p "$base_dir"

  echo "=== Clustering '$infa' -> $base_dir (id=$ID strand=$STRAND) ==="

  # run vsearch and ask it to write cluster files with the given prefix inside base_dir.
  # vsearch numbers them cluster_${base}_0, cluster_${base}_1, ... (no zero padding) in
  # the order clusters were created. See PIPELINE.md Step 5/6 for the strand choice.
  rm -f "$base_dir"/cluster_"$base"_* "$base_dir/${base}"_cluster_*.fa
  # Step 7 writes one <base>_cluster_N_rotated/ folder per cluster and never removes them; folders left by an earlier
  # clustering (other motif set, other strand, other cluster numbers) would be pooled by Step 8 as if current
  rm -rf "$base_dir"/*_rotated
  vsearch --cluster_fast "$infa" --id "$ID" --strand "$STRAND" --clusters "$base_dir/cluster_${base}_" --quiet

  # rename to ${base}_cluster_N.fa keeping vsearch's own cluster number + 1, so the name
  # identifies the cluster regardless of lexical file order.
  idx=0
  for cf in "$base_dir"/cluster_"$base"_*; do
    [ -e "$cf" ] || continue
    k="${cf##*_}"
    mv -- "$cf" "$base_dir/${base}_cluster_$((k + 1)).fa"
    idx=$((idx + 1))
  done

  if (( idx == 0 )); then
    echo "WARNING: vsearch produced no cluster files for '$infa'."
    echo "  - Check vsearch version / options (--clusters must be supported by this build)."
  else
    echo "Created $idx cluster file(s) in: $base_dir"
  fi
done

echo "All done. Top-level output: $OUTPUT_DIR"
