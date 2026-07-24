#!/bin/bash
set -euo pipefail

# ==========================
# Config
# ==========================
INPUT_DIR="."
OUTPUT_DIR="."

# ==========================
# Dependency Check
# ==========================
command -v mafft >/dev/null 2>&1 || { echo "Error: mafft not found. Please install MAFFT."; exit 1; }
command -v awk >/dev/null 2>&1 || { echo "Error: awk not found. Please install awk."; exit 1; }

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# ==========================
# Process each *.fa file with >1 sequence
# ==========================
for infa in "$INPUT_DIR"/*.fa; do
  [ -e "$infa" ] || { echo "No *.fa files found in $INPUT_DIR"; exit 1; }

  # Check if file has more than one sequence
  seq_count=$(grep -c '^>' "$infa" || echo 0)
  if [ "$seq_count" -le 1 ]; then
    echo "Skipping $infa: contains $seq_count sequence(s) (requires >1)"
    continue
  fi

  # Extract base filename without path and extension
  base=$(basename "$infa" .fa)

  # Create subfolder for rotated output
  ROTATED_DIR="$OUTPUT_DIR/${base}_rotated"
  mkdir -p "$ROTATED_DIR"

  # Step 0: MAFFT alignment
  echo "Aligning $infa with MAFFT..."
  aligned_fasta="$ROTATED_DIR/${base}.aligned.fasta"
  mafft --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --preservecase --quiet "$infa" > "$aligned_fasta" 2>/dev/null || {
    echo "Error: MAFFT alignment failed for $infa"
    continue
  }
  echo "Alignment done: $aligned_fasta"

  # Step 1: Find rightmost first non-gap position across all sequences
  cons_start=$(awk '
    BEGIN { max_start = 0 }
    /^>/ { 
      if (seq != "") { 
        for (i=1; i<=length(seq); i++) { 
          if (substr(seq, i, 1) != "-") { 
            if (i > max_start) max_start = i
            break 
          } 
        } 
        seq = "" 
      } 
    } 
    !/^>/ { seq = seq $0 } 
    END { 
      if (seq != "") { 
        for (i=1; i<=length(seq); i++) { 
          if (substr(seq, i, 1) != "-") { 
            if (i > max_start) max_start = i
            break 
          } 
        } 
      } 
      print max_start ? max_start : 1 
    }' "$aligned_fasta" 2>/dev/null || {
      echo "Error: Failed to find consensus start in $aligned_fasta"
      continue
    })

  # Check if cons_start is valid
  if [ "$cons_start" -eq 0 ]; then
    echo "Error: No non-gap characters found in $aligned_fasta. Skipping rotation."
    continue
  fi

  # Step 2: Rotate sequences and realign with MAFFT
  echo "Rotating sequences for $infa (consensus start: $cons_start) and realigning..."
  outfa="$ROTATED_DIR/${base}.rotated.fa"
  awk -v cs="$cons_start" '
    BEGIN { header = ""; seq = "" }
    /^>/ {
      if (header != "" && seq != "") {
        pre = substr(seq, 1, cs - 1)
        suf = substr(seq, cs)
        print header
        print suf pre
      }
      header = $0
      seq = ""
      next
    }
    { seq = seq $0 }
    END {
      if (header != "" && seq != "") {
        pre = substr(seq, 1, cs - 1)
        suf = substr(seq, cs)
        print header
        print suf pre
      }
    }' "$aligned_fasta" | mafft --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --preservecase --quiet - > "$outfa" 2>/dev/null || {
    echo "Error: MAFFT realignment failed for rotated sequences of $infa"
    continue
  }
  echo "Rotation and realignment done: $outfa"
done

echo "Pipeline completed successfully!"