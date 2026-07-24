#!/usr/bin/env bash
set -euo pipefail
cd /staging/tmp/clsat_v6/scan85/chunks/sub44/sub44_4_split
samtools faidx half1.fa
samtools faidx half1.fa "half1:1-16500" | sed "1s/.*/>q1/" > q1.fa
samtools faidx half1.fa "half1:16501-33000" | sed "1s/.*/>q2/" > q2.fa
ls q1.fa q2.fa | parallel --jobs 2 'trf {} 2 5 7 80 10 50 2000 -h -ngs > {/.}.log 2>&1'
touch quarters.done
echo DONE_QUARTERS
