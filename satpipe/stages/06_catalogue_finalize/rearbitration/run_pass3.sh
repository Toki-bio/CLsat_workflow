#!/usr/bin/env bash
set -euo pipefail
# Symmetric rescan (stage 13 / "PASS3"): scan a fixed locus BED with EVERY family consensus at
# monomer resolution (window 2 x MONO, step MONO), without re-expanding.
#
# Usage: run_pass3.sh <loci.bed> <genomes.tsv> <families.tsv> <out_dir> [--jobs N]
#   loci.bed      BED3+ with col4 "species|..." (e.g. MASTER.merged.bed or a corrected catalogue)
#   genomes.tsv   species_code  genome.fa
#   families.tsv  family  mono  consensus.fa
#   out_dir       receives scan_<family>.bits.tsv (+ the other --scan-only outputs) and pass3.done
# Environment: JOBS [16]  (overridden by --jobs)
[[ $# -ge 4 ]] || { echo "Usage: $0 <loci.bed> <genomes.tsv> <families.tsv> <out_dir> [--jobs N]" >&2; exit 1; }
BED="$1"; GENOMES="$2"; FAMTAB="$3"; OUT="$4"; shift 4
JOBS="${JOBS:-16}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --jobs) JOBS="$2"; shift 2;;
    *) echo "Unknown option: $1" >&2; exit 1;;
  esac
done
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLI="$HERE/../../02_seed_expansion_scan/clisat_expand_then_scan_merged.sh"
[[ -s "$CLI" ]] || { echo "expansion script not found: $CLI" >&2; exit 1; }
mkdir -p "$OUT"

while IFS=$'\t' read -r fam mono cons _; do
  [[ -z "$fam" || "$fam" == \#* ]] && continue
  echo "=== PASS3 scan: $fam (MONO=$mono) ==="
  bash "$CLI" "$BED" "$GENOMES" "$cons" "$mono" "$OUT/scan_$fam" --scan-only --jobs "$JOBS" 2>&1 \
    | grep -E '^\[expand_then_scan\] (PASS|DONE|SCAN_ONLY)' || true
  echo "  bits rows: $(wc -l < "$OUT/scan_${fam}.bits.tsv")"
done < "$FAMTAB"
echo ALL_PASS3_DONE > "$OUT/pass3.done"
