#!/usr/bin/env python3
"""Rescan each candidate locus with 5 kb of context on both sides.

The 7 candidates were scanned over their raw TRF interval only. Catalogue loci are
seed-expanded first (Methods 4.1), walking outward while CLsat similarity holds, so a
TRF-interval-only scan understates any array that continues past its seed. This scans
locus +/- 5000 bp (FLANK_VIEW) at the same window/step, and reports the contiguous run
of windows above the 50-bit satellite floor that contains the locus -- i.e. the extent
the array would have had under boundary expansion.

Also emits flank-track windows so the viewer shows scanned-and-empty rather than
unexamined background around these loci.
"""
import json
import os
import subprocess

BASE = '/staging/Sergei/Lacertidae/Darevskia'
OUT = BASE + '/Satellites/motifs_output/clean/candidates'
SS = '/staging/conda/envs/bioinfo/bin/ssearch36'
SAT_BITS = 50.0
VIEW = 5000

GEN = {
    'unp': BASE + '/unp/GCA_032766555.1_D.uni.pHiFi.1.0_genomic.fna',
    'val': BASE + '/val/GCA_034642135.1_Dar_val_2.0_genomic.fna',
    'nai': BASE + '/nai/GCA_034642065.1_Dar_nair_1.0_genomic.fna',
    'unm': BASE + '/unm/GCA_032766585.1_D.uni.mHiFi.1.0_genomic.fna',
    'arm': BASE + '/arm/GCA_034641395.1_Dar_arm_1.0_genomic.fna',
}
MONO = {'CLsat1': 146, 'CLsat2': 145, 'CLsat3': 146, 'CLsat4': 116}

data = json.load(open('%s/candidates_summary.json' % OUT))
summary = data['summary']


def scaf_len(sp, scaf):
    for line in open(GEN[sp] + '.fai'):
        p = line.split('\t')
        if p[0] == scaf:
            return int(p[1])


def faidx(sp, scaf, a, b):
    r = subprocess.run(['samtools', 'faidx', GEN[sp], '%s:%d-%d' % (scaf, a, b)],
                       capture_output=True, text=True)
    return ''.join(l.strip() for l in r.stdout.splitlines()[1:]).upper()


def bits(seq, fam):
    with open('%s/_q5.fa' % OUT, 'w') as f:
        f.write('>q\n%s\n' % seq)
    o = subprocess.run([SS, '-n', '-m', '8C', '-E', '1000', '-b', '1',
                        '%s/_q5.fa' % OUT, '%s/%s.doubled.fa' % (OUT, fam)],
                       capture_output=True, text=True).stdout
    best = 0.0
    for line in o.splitlines():
        if not line.startswith('#'):
            p = line.split('\t')
            if len(p) >= 12:
                try:
                    best = max(best, float(p[11]))
                except ValueError:
                    pass
    return best


flank_rows = []
print('%-3s %-4s %-11s %-11s  %-9s %-11s %s'
      % ('id', 'asm', 'TRF start', 'TRF end', 'TRF len', 'signal span', 'change'))
report = {}
for cid in sorted(summary, key=lambda x: int(x[1:])):
    s = summary[cid]
    sp, scaf = s['sp'], s['scf']
    a, b = s['start'] + 1, s['end']
    slen = scaf_len(sp, scaf)
    ra, rb = max(1, a - VIEW), min(slen, b + VIEW)
    seq = faidx(sp, scaf, ra, rb)
    fam = s['call']
    m = MONO[fam]
    win, step = 2 * m, m
    wins = []
    pos = 0
    while pos + win <= len(seq):
        gpos = ra + pos                      # 1-based genomic start of window
        bt = bits(seq[pos:pos + win], fam)
        wins.append((gpos, bt))
        if gpos + win - 1 < a or gpos > b:   # window lies wholly in the flank
            flank_rows.append((scaf, gpos - 1, bt,
                               '%s|%s:%d-%d|CANDFLANK' % (sp, scaf, a, b), 'FLANK', sp))
        pos += step
    # contiguous run above the floor that overlaps the TRF interval
    above = [i for i, (g, bt) in enumerate(wins) if bt >= SAT_BITS]
    lo = hi = None
    for i in above:
        g = wins[i][0]
        if g + win - 1 >= a and g <= b:
            lo = i if lo is None else lo
            hi = i
    if lo is not None:
        while lo - 1 in above:
            lo -= 1
        while hi + 1 in above:
            hi += 1
        span_a, span_b = wins[lo][0], wins[hi][0] + win - 1
    else:
        span_a, span_b = a, b
    grew = (a - span_a) + (span_b - b)
    report[cid] = dict(trf=[a, b], span=[span_a, span_b],
                       trf_len=b - a + 1, span_len=span_b - span_a + 1,
                       left_ext=a - span_a, right_ext=span_b - b,
                       scaffold_len=slen, clipped=(ra == 1 or rb == slen))
    print('%-3s %-4s %11d %11d  %9d %11d %s'
          % (cid, sp, a, b, b - a + 1, span_b - span_a + 1,
             ('+%d bp (L+%d R+%d)' % (grew, a - span_a, span_b - b)) if grew else 'unchanged'))

with open('%s/CAND_FLANK.bits.tsv' % OUT, 'w') as f:
    for r in flank_rows:
        f.write('%s\t%d\t%.1f\t%s\t%s\t%s\n' % r)
data['extent'] = report
with open('%s/candidates_summary.json' % OUT, 'w') as f:
    json.dump(data, f, indent=1)
print('\nflank-track windows written: %d' % len(flank_rows))
