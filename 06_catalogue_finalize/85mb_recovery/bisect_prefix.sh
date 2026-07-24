#!/usr/bin/env bash
set -euo pipefail
cd /staging/tmp/clsat_v6/scan85/chunks/sub44/sub44_4_split
samtools faidx sub44_4_4_prefix.fa "sub44_4_4_prefix:1-33000" | sed "1s/.*/>half1/" > half1.fa
samtools faidx sub44_4_4_prefix.fa "sub44_4_4_prefix:33001-66500" | sed "1s/.*/>half2/" > half2.fa
ls half1.fa half2.fa | parallel --jobs 2 'trf {} 2 5 7 80 10 50 2000 -h -ngs > {/.}.log 2>&1'
touch halves.done
echo DONE_HALVES
