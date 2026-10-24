#!/bin/bash
set -euo pipefail
# Step 3 - extract candidate tandem-repeat motifs from TRF .dat files.
#
# Usage: Step3.sh [dat_dir] [out_dir]
#   dat_dir   directory searched (recursively) for *.dat files; the basename of each .dat file
#             is the species code written into every FASTA header            [default: .]
#   out_dir   output directory, one <code>.fa per .dat                          [default: motifs_output]
# Environment (defaults = published Darevskia run):
#   PERIOD_MIN  minimum TRF period accepted                                    [80]
#   PERIOD_MAX  maximum TRF period accepted                                    [300]
#   MIN_COPIES  minimum TRF copy number accepted                               [10]
#
# Header written per motif: >code|SCAFFOLD:START-END|len=PERIOD|copies=COPIES|locuslen=N
# (START-END are TRF's 1-based inclusive coordinates.)

IN_DIR="${1:-.}"
OUTPUT_DIR="${2:-motifs_output}"
MIN_LEN="${PERIOD_MIN:-80}"
MAX_LEN="${PERIOD_MAX:-300}"
MIN_COPIES="${MIN_COPIES:-10}"
mkdir -p "$OUTPUT_DIR"

# -P => never follow symlinks; include regular files and symlinks
find -P "$IN_DIR" \( -type f -o -type l \) -name '*.dat' -print0 |
while IFS= read -r -d '' dat_file; do
  # Always use symlink/file's own name for output file and FASTA headers
  output_base=$(basename -- "$dat_file" .dat)
  output_file="${OUTPUT_DIR}/${output_base}.fa"

  # Resolve symlink target for processing
  if [[ -L "$dat_file" ]]; then
    # Get absolute path of target
    target=$(readlink -f -- "$dat_file" 2>/dev/null)
    if [[ ! -f "$target" || ! -r "$target" ]]; then
      # a genome without readable TRF output must stop the run, not drop out of it silently
      echo "FATAL: Symlink $dat_file points to inaccessible or non-existent file: $target" >&2
      exit 1
    fi
  else
    if [[ ! -r "$dat_file" ]]; then
      echo "FATAL: File $dat_file is not readable" >&2
      exit 1
    fi
    target="$dat_file"
  fi

  echo "Processing $dat_file (target: $target)..."

  awk -v MIN_LEN="$MIN_LEN" -v MAX_LEN="$MAX_LEN" -v MIN_COPIES="$MIN_COPIES" -v base_name="$output_base" '
    BEGIN { found_data=0 }
    /^[ \t]*Sequence:/ { scaffold=$2; next }   # accession right after "Sequence:"
    /^[0-9]/ {
      # TRF numeric rows; be robust to column count differences
      if (NF < 6) next
      start   = $1+0
      endpos  = $2+0
      period  = $3+0
      copies  = $4+0
      motif   = toupper($(NF-1))   # consensus unit (penultimate field)
      locusln = length($NF)        # length of full repeated sequence (last field)

      if (period >= MIN_LEN && period <= MAX_LEN && copies >= MIN_COPIES) {
        printf(">%s|%s:%d-%d|len=%s|copies=%s|locuslen=%d\n",
               base_name, scaffold, start, endpos, period, copies, locusln)
        print motif
        found_data=1
      }
    }
    END {
      if (!found_data) {
        print "Warning: No valid motifs found in file" > "/dev/stderr"
      }
    }
  ' "$target" > "$output_file"

  motif_count=$(grep -c '^>' "$output_file" || true)   # grep -c prints 0 itself on no match
  echo "  Saved to $output_file ($motif_count motifs)"
done

echo "Processing complete. Results in $OUTPUT_DIR/"
