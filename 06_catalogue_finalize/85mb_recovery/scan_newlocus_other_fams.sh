#!/usr/bin/env bash
set -euo pipefail
V=/staging/tmp/clsat_v6
SD=/staging/Sergei/Lacertidae/Darevskia/Satellites/motifs_output/clean/cross_clusters/selected/subclusters/seeds
echo -e "JAWDKN010000003.1\t85344126\t85768940\tunp|JAWDKN010000003.1:85344127-85768940|MASTER\tunp" \
    > $V/scan85/newlocus.tagged.bed

declare -A MONO=([CLsat2]=145 [CLsat3]=146 [CLsat4]=116)
for fam in CLsat2 CLsat3 CLsat4; do
  bash $V/cli_v6.sh $V/scan85/newlocus.tagged.bed $SD/genomes.tsv $SD/$fam.cons.fa ${MONO[$fam]} \
      $V/scan85/newlocus_$fam --scan-only --jobs 16
done
touch $V/scan85/newlocus_other_fams.done
echo DONE_OTHER_FAMS
