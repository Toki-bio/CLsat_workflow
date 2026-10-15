#!/usr/bin/env bash
# STAGE 1 of the re-run: re-expand CLsat1 ONLY, all 6 genomes, with the fixed cli_v6.sh.
#
# Scope is CLsat1 alone because the Jun-14 QLEN edit broke only that family:
#   cli.sh.bak .. cli.sh_щдв  (Jun 9-12)  QLEN = QREP*MONO     = 292  -> MIN_ALN 292/438  OK
#   cli.sh                    (Jun 14)    QLEN = (QREP+1)*MONO = 438  -> MIN_ALN 438/438  IMPOSSIBLE
# and only out_CLsat1_v4.* was regenerated after that edit (Jun 14 05:30).
# CLsat2/3/4 still carry their Jun-12 results from the working script and are left alone.
#
# Writes ONLY to /staging/tmp/clsat_v6/full — nothing published is touched.
set -euo pipefail

SD=/staging/Sergei/Lacertidae/Darevskia/Satellites/motifs_output/clean/cross_clusters/selected/subclusters/seeds
V6=/staging/tmp/clsat_v6
OUT=$V6/full
BIN=/staging/conda/envs/bioinfo/bin
export PATH="$BIN:/staging/miniconda3/bin:$PATH"

mkdir -p "$OUT"; cd "$OUT"
[[ -s $V6/cli_v6.sh ]] || { echo "cli_v6.sh missing"; exit 1; }

echo "=== inputs ==="
printf '  CLsat1 seeds : %d (all species)\n' "$(wc -l < "$SD/CLsat1.bed")"
awk -F'\t' '{split($4,a,"|"); n[a[1]]++} END{for(k in n) printf "     %-4s %3d\n", k, n[k]}' "$SD/CLsat1.bed"
printf '  genomes      : %d\n' "$(wc -l < "$SD/genomes.tsv")"
echo "  MIN_ALN check: $(grep -n 'MIN_ALN=\$(awk' $V6/cli_v6.sh)"

echo
echo "=== PASS1+2: expand + merge (CLsat1, MONO=146) ==="
# --scan-only is NOT used: we want PASS1. PASS3 (the bit-score rescan) is a separate,
# much more expensive stage and is deliberately deferred.
bash "$V6/cli_v6.sh" "$SD/CLsat1.bed" "$SD/genomes.tsv" "$SD/CLsat1.cons.fa" 146 \
     "$OUT/v6full_CLsat1" --jobs 32 2>&1 | grep -E '^\[expand_then_scan\] (PASS|DONE)' || true

echo
echo "=== RESULT: v6 vs v4, CLsat1, per species ==="
printf '%-5s %9s %9s %9s %9s %9s\n' sp seeds v4_loci v6_loci v4_Mb v6_Mb
for s in arm mix nai unm unp val; do
  seed=$(awk -F'\t' -v S=$s '{split($4,a,"|"); if(a[1]==S) print}' "$SD/CLsat1.bed" | wc -l)
  [[ $seed -eq 0 ]] && continue
  v4=$(awk -F'\t' -v S=$s '{split($4,a,"|"); if(a[1]==S) print}' "$SD/out_CLsat1_v4.expanded_merged.bed" \
       | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - 2>/dev/null)
  v6=$(awk -F'\t' -v S=$s '{split($4,a,"|"); if(a[1]==S) print}' "$OUT/v6full_CLsat1.expanded_merged.bed" \
       | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - 2>/dev/null)
  printf '%-5s %9d %9d %9d %9.3f %9.3f\n' "$s" "$seed" \
    "$(echo "$v4" | grep -c . )" "$(echo "$v6" | grep -c . )" \
    "$(echo "$v4" | awk '{bp+=$3-$2} END{print bp/1e6}')" \
    "$(echo "$v6" | awk '{bp+=$3-$2} END{print bp/1e6}')"
done

echo
echo "=== stop reasons (v6) ==="
awk -F'\t' 'NR>1{R[$8]++; R[$9]++} END{for(k in R) printf "  %-20s %d\n", k, R[k]}' \
  "$OUT/v6full_CLsat1.expanded_per_locus.tsv" | sort -k2,2nr

echo
echo "=== sanity: did v6 DROP anything v4 had? (must be 0) ==="
cut -f1-3 "$SD/out_CLsat1_v4.expanded_merged.bed" | sort -k1,1 -k2,2n | bedtools merge -i - > /tmp/f_v4.bed
cut -f1-3 "$OUT/v6full_CLsat1.expanded_merged.bed" | sort -k1,1 -k2,2n | bedtools merge -i - > /tmp/f_v6.bed
bedtools subtract -a /tmp/f_v4.bed -b /tmp/f_v6.bed | awk '{bp+=$3-$2} END{printf "  v4 minus v6 : %d bp\n", bp+0}'
bedtools subtract -a /tmp/f_v6.bed -b /tmp/f_v4.bed | awk '{bp+=$3-$2} END{printf "  v6 minus v4 : %d bp (newly called)\n", bp+0}'
echo "DONE — $OUT"
