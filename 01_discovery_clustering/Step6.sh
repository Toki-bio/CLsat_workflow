#!/usr/bin/env bash
set -euo pipefail

# INPUT_DIR: where your .rot files live
# OUTPUT_DIR: top-level folder that will receive one subfolder per .rot input
INPUT_DIR="."
OUTPUT_DIR="clusters_by_file"
# clustering identity (tune as needed)
ID=0.90

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

  echo "=== Clustering '$infa' -> $base_dir (id=$ID) ==="

  # run vsearch and ask it to write cluster files with the given prefix inside base_dir
  # created files will look like: cluster_${base}_00001.fa, cluster_${base}_00002.fa, ...
  vsearch --cluster_fast "$infa" --id "$ID" --clusters "$base_dir/cluster_${base}_" --quiet

  # rename cluster files to a neat, predictable pattern: ${base}_cluster_1.fa, ${base}_cluster_2.fa ...
  idx=0
  for cf in "$base_dir"/cluster_"$base"_*; do
    # nullglob ensures loop is skipped if no match
    [ -e "$cf" ] || continue
    idx=$((idx + 1))
    # use -n to avoid overwriting any existing file
    mv -n -- "$cf" "$base_dir/${base}_cluster_${idx}.fa"
  done

  if (( idx == 0 )); then
    echo "WARNING: vsearch produced no cluster files for '$infa'."
    echo "  - Check vsearch version / options. Some vsearch builds may not write cluster files."
    echo "  - If that happens I can provide a fallback that uses --uc to map sequences to clusters."
  else
    echo "Created $idx cluster file(s) in: $base_dir"
  fi
done

echo "All done. Top-level output: $OUTPUT_DIR"
