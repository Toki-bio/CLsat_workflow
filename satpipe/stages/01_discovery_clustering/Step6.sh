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
#   ORIENT      1 = write every member in its centroid's orientation (members that vsearch matched on the minus
#               strand are reverse-complemented, from the .uc table); 0 = as vsearch writes them  [1]
#               (with STRAND=plus no member is on the minus strand, so ORIENT changes nothing there)

INPUT_DIR="${1:-.}"
OUTPUT_DIR="${2:-clusters_by_file}"
ID="${CLUSTER_ID:-0.90}"
STRAND="${STRAND:-plus}"
ORIENT="${ORIENT:-1}"

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
  uc="$base_dir/${base}.uc"
  vsearch --cluster_fast "$infa" --id "$ID" --strand "$STRAND" --clusters "$base_dir/cluster_${base}_" --uc "$uc" --quiet

  # rename to ${base}_cluster_N.fa keeping vsearch's own cluster number + 1, so the name
  # identifies the cluster regardless of lexical file order.
  idx=0
  for cf in "$base_dir"/cluster_"$base"_*; do
    [ -e "$cf" ] || continue
    k="${cf##*_}"
    mv -- "$cf" "$base_dir/${base}_cluster_$((k + 1)).fa"
    idx=$((idx + 1))
  done

  # VSEARCH writes each member as it was in the input; a member matched on the minus strand would enter the Step 7
  # alignment reverse-complemented to the rest of its cluster. Turn those members round (uc: H rows, column 5 '-',
  # column 2 = cluster number, column 9 = member label up to the first blank).
  if [ "$ORIENT" = 1 ] && (( idx > 0 )); then
    nrev=$(awk -F'	' '$1=="H" && $5=="-"' "$uc" | wc -l)
    if (( nrev > 0 )); then
      awk -F'	' -v D="$base_dir" -v B="$base" '
        function rc(s,  i, o, c) { o = ""; for (i = length(s); i > 0; i--) { c = substr(s, i, 1)
            o = o (c=="A"?"T":c=="T"?"A":c=="C"?"G":c=="G"?"C":c=="a"?"t":c=="t"?"a":c=="c"?"g":c=="g"?"c":c) } return o }
        function flush() { if (h != "") { print h > out; print (h in R ? rc(sq) : sq) > out }; h = ""; sq = "" }
        FNR == NR { if ($1 == "H" && $5 == "-") { split($9, w, " "); R[">" w[1]] = 1; F[$2 + 1] = 1 }; next }
        END {
          for (k in F) {
            f = D "/" B "_cluster_" k ".fa"; out = f ".tmp"; h = ""; sq = ""
            while ((getline l < f) > 0) {
              if (l ~ /^>/) { flush(); split(l, w, " "); h = l; key = w[1]; if (key in R) R[h] = 1 }
              else sq = sq l
            }
            flush(); close(f); close(out)
            if (system("mv -- \"" out "\" \"" f "\"") != 0) { print "FATAL: cannot replace " f > "/dev/stderr"; exit 1 }
          }
        }' "$uc" /dev/null
      echo "Turned $nrev minus-strand member(s) into their centroid's orientation"
    fi
  fi

  if (( idx == 0 )); then
    echo "WARNING: vsearch produced no cluster files for '$infa'."
    echo "  - Check vsearch version / options (--clusters must be supported by this build)."
  else
    echo "Created $idx cluster file(s) in: $base_dir"
  fi
done

echo "All done. Top-level output: $OUTPUT_DIR"
