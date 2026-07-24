#!/usr/bin/env bash
set -euo pipefail
V=/staging/tmp/clsat_v6
SD=/staging/Sergei/Lacertidae/Darevskia/Satellites/motifs_output/clean/cross_clusters/selected/subclusters/seeds
bash $V/cli_v6.sh $V/scan85/new_candidates_seed.bed $SD/genomes.tsv $SD/CLsat1.cons.fa 146 \
    $V/scan85/new_candidates_out --jobs 16
touch $V/scan85/new_candidates.done
echo DONE_NEW_CANDIDATES
