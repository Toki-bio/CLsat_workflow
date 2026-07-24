#!/usr/bin/env bash
set -euo pipefail
V=/staging/tmp/clsat_v6
SD=/staging/Sergei/Lacertidae/Darevskia/Satellites/motifs_output/clean/cross_clusters/selected/subclusters/seeds
BED=$V/pass3/corrected_loci.tagged.bed
declare -A MONO=([CLsat1]=146 [CLsat2]=145 [CLsat3]=146 [CLsat4]=116)
for fam in CLsat1 CLsat2 CLsat3 CLsat4; do
  echo "=== PASS3 scan: $fam (MONO=${MONO[$fam]}) ==="
  bash $SD/cli.sh $BED $SD/genomes.tsv $SD/$fam.cons.fa ${MONO[$fam]}        $V/pass3/scan_$fam --scan-only --jobs 32 2>&1 | grep -E '^\[expand_then_scan\] (PASS|DONE|SCAN_ONLY)'
  echo "  bits rows: $(wc -l < $V/pass3/scan_${fam}.bits.tsv)"
done
echo ALL_PASS3_DONE > $V/pass3/pass3.done
