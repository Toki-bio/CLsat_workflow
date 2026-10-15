#!/usr/bin/env bash
set -euo pipefail
SD=/staging/Sergei/Lacertidae/Darevskia/Satellites/motifs_output/clean/cross_clusters/selected/subclusters/seeds
V=/staging/tmp/clsat_v6
OUT=$V/clsat1_verify_rescan_20260716
mkdir -p $OUT

echo "=== Full CLsat1 boundary-expansion re-scan, from raw seeds, corrected pipeline ==="
echo "Input seeds : $SD/CLsat1.bed (750 rows, dated 2026-06-09, pre-expansion)"
echo "Consensus   : $SD/CLsat1.cons.fa"
echo "Script      : $V/cli_v6.sh (corrected MIN_ALN, MIN_ALN_FRAC=0.67)"
echo "Started     : $(date -Iseconds)"
echo

bash $V/cli_v6.sh $SD/CLsat1.bed $SD/genomes.tsv $SD/CLsat1.cons.fa 146 \
    $OUT/CLsat1_full --jobs 24 2>&1

echo
echo "Finished    : $(date -Iseconds)"
echo
echo "=== summary ==="
wc -l $OUT/CLsat1_full.expanded_per_locus.tsv $OUT/CLsat1_full.expanded_merged.bed
echo "zero-expansion loci (of 750 seeds):"
awk -F'\t' 'NR>1 && $NF==0' $OUT/CLsat1_full.expanded_per_locus.tsv | wc -l
echo "expansion distance distribution (bp, top 10):"
awk -F'\t' 'NR>1{print $NF}' $OUT/CLsat1_full.expanded_per_locus.tsv | sort -n | uniq -c | sort -rn | head -10

touch $OUT/DONE
echo CLSAT1_FULL_RESCAN_DONE
