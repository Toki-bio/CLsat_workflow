#!/usr/bin/env bash
set -euo pipefail
V=/staging/tmp/clsat_v6
SD=/staging/Sergei/Lacertidae/Darevskia/Satellites/motifs_output/clean/cross_clusters/selected/subclusters/seeds
D=/staging/Sergei/Lacertidae/Darevskia/Satellites/manuscript/CLsat_full_clusters_6genome
OUT=$V/table1_540
mkdir -p $OUT

echo "=== build MASTER_BED (5-col) from the 540-locus catalog ==="
awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$3,$5"|"$1":"($2+1)"-"$3"|MASTER",$6}' \
    $V/scan85/loci_540.bed > $OUT/MASTER.merged.bed
wc -l $OUT/MASTER.merged.bed

echo
echo "=== patch arbitration_patched.sh ==="
sed \
  -e "s#^BASE=.*#BASE=\"$OUT\"#" \
  -e "s#^OUT=.*#OUT=\"$OUT\"#" \
  -e "s#\$OUT/MASTER_CLsat1.bits.tsv#$V/pass3/scan_CLsat1.bits.tsv#" \
  -e "s#\$OUT/MASTER_CLsat2.bits.tsv#$V/pass3/scan_CLsat2.bits.tsv#" \
  -e "s#\$OUT/MASTER_CLsat3.bits.tsv#$V/pass3/scan_CLsat3.bits.tsv#" \
  -e "s#\$OUT/MASTER_CLsat4.bits.tsv#$V/pass3/scan_CLsat4.bits.tsv#" \
  -e "s#\$BASE/CLsat1.bed#$SD/CLsat1.bed#" \
  -e "s#\$BASE/CLsat2.bed#$SD/CLsat2.bed#" \
  -e "s#\$BASE/CLsat3.bed#$SD/CLsat3.bed#" \
  -e "s#\$BASE/CLsat4.bed#$SD/CLsat4.bed#" \
  $D/arbitration_patched.sh > $OUT/arbitration_540.sh
chmod +x $OUT/arbitration_540.sh

echo
echo "=== run it ==="
bash $OUT/arbitration_540.sh

echo
echo "=== Table 1 major-family x species (540-locus catalog) ==="
python3 - "$OUT/MASTER.final_assignment.tsv" "$OUT/MASTER.merged.bed" << 'PYEOF'
import sys, collections
fa = sys.argv[1]; mb = sys.argv[2]
MONO = {'CLsat1': 146, 'CLsat2': 145, 'CLsat3': 146, 'CLsat4': 116}
sp_of = {}
with open(mb) as f:
    for i, line in enumerate(f, 1):
        p = line.rstrip('\n').split('\t')
        sp_of[f'L{i}'] = p[3].split('|')[0]
seg_len = {}
with open(f'{sys.argv[1].rsplit("/",1)[0]}/MASTER.loci.segmented.bed') as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        seg_len[p[3]] = int(p[2]) - int(p[1])
loci = collections.Counter(); copies = collections.Counter()
with open(fa) as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        lid, fam = p[0], p[1]
        base = lid.split('.')[0]
        sp = sp_of.get(base, '?')
        length = seg_len.get(lid, 0)
        loci[(fam, sp)] += 1
        copies[(fam, sp)] += round(length / MONO[fam])
species = ['mix', 'arm', 'val', 'unp', 'unm', 'nai']
fams = ['CLsat1', 'CLsat2', 'CLsat3', 'CLsat4']
print('%-8s' % 'Family', *['%14s' % s for s in species], '%14s' % 'Total')
grand_loci = 0; grand_bp = 0
for fam in fams:
    row_loci = [loci[(fam, s)] for s in species]
    row_cn = [copies[(fam, s)] for s in species]
    cells = ['%d / %s' % (l, format(c, ',')) for l, c in zip(row_loci, row_cn)]
    tot_loci, tot_cn = sum(row_loci), sum(row_cn)
    grand_loci += tot_loci
    print('%-8s' % fam, *['%14s' % c for c in cells],
          '%14s' % ('%d / %s' % (tot_loci, format(tot_cn, ','))))
print('GRAND TOTAL segmented loci:', grand_loci)
PYEOF
