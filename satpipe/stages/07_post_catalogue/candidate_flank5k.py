#!/usr/bin/env python3
"""Rescan each candidate locus with VIEW bp of context on both sides.

The candidates were scanned over their raw TRF interval only. Catalogue loci are
seed-expanded first (Methods 4.1), walking outward while similarity holds, so a
TRF-interval-only scan understates any array that continues past its seed. This scans
locus +/- VIEW bp at the same window/step, and reports the contiguous run of windows above
the SAT_BITS satellite floor that contains the locus -- i.e. the extent the array would have
had under boundary expansion.

Also emits flank-track windows so the viewer shows scanned-and-empty rather than
unexamined background around these loci.

Usage: candidate_flank5k.py <out_dir> <genomes.tsv> <families.tsv> [--view 5000] [--sat-bits 50]
  out_dir  holds candidates_summary.json and <family>.doubled.fa from candidate_scan.py
Tools: ssearch36, samtools on PATH.
"""
import argparse
import json
import shutil
import subprocess
import sys

ap = argparse.ArgumentParser()
ap.add_argument('out')
ap.add_argument('genomes')
ap.add_argument('families')
ap.add_argument('--view', type=int, default=5000)
ap.add_argument('--sat-bits', type=float, default=50.0)
a = ap.parse_args()

SS = shutil.which('ssearch36') or sys.exit('ssearch36 not on PATH')
SAM = shutil.which('samtools') or sys.exit('samtools not on PATH')
OUT = a.out
SAT_BITS = a.sat_bits
VIEW = a.view


def read_tsv(path, min_cols):
    rows = []
    for line in open(path):
        if not line.strip() or line.startswith('#'):
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) < min_cols:
            sys.exit('%s: need %d columns, got %r' % (path, min_cols, line.strip()))
        rows.append(p)
    return rows


GEN = {p[0]: p[1] for p in read_tsv(a.genomes, 2)}
MONO = {p[0]: int(p[1]) for p in read_tsv(a.families, 3)}

data = json.load(open('%s/candidates_summary.json' % OUT))
summary = data['summary']


def scaf_len(sp, scaf):
    for line in open(GEN[sp] + '.fai'):
        p = line.split('\t')
        if p[0] == scaf:
            return int(p[1])


def faidx(sp, scaf, x, y):
    r = subprocess.run([SAM, 'faidx', GEN[sp], '%s:%d-%d' % (scaf, x, y)],
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
for cid in sorted(summary, key=lambda c: int(c[1:]) if c[1:].isdigit() else 0):
    s = summary[cid]
    sp, scaf = s['sp'], s['scf']
    x, y = s['start'] + 1, s['end']
    slen = scaf_len(sp, scaf)
    ra, rb = max(1, x - VIEW), min(slen, y + VIEW)
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
        if gpos + win - 1 < x or gpos > y:   # window lies wholly in the flank
            flank_rows.append((scaf, gpos - 1, bt,
                               '%s|%s:%d-%d|CANDFLANK' % (sp, scaf, x, y), 'FLANK', sp))
        pos += step
    # contiguous run above the floor that overlaps the TRF interval
    above = [i for i, (g, bt) in enumerate(wins) if bt >= SAT_BITS]
    lo = hi = None
    for i in above:
        g = wins[i][0]
        if g + win - 1 >= x and g <= y:
            lo = i if lo is None else lo
            hi = i
    if lo is not None:
        while lo - 1 in above:
            lo -= 1
        while hi + 1 in above:
            hi += 1
        span_a, span_b = wins[lo][0], wins[hi][0] + win - 1
    else:
        span_a, span_b = x, y
    grew = (x - span_a) + (span_b - y)
    report[cid] = dict(trf=[x, y], span=[span_a, span_b],
                       trf_len=y - x + 1, span_len=span_b - span_a + 1,
                       left_ext=x - span_a, right_ext=span_b - y,
                       scaffold_len=slen, clipped=(ra == 1 or rb == slen))
    print('%-3s %-4s %11d %11d  %9d %11d %s'
          % (cid, sp, x, y, y - x + 1, span_b - span_a + 1,
             ('+%d bp (L+%d R+%d)' % (grew, x - span_a, span_b - y)) if grew else 'unchanged'))

with open('%s/CAND_FLANK.bits.tsv' % OUT, 'w') as f:
    for r in flank_rows:
        f.write('%s\t%d\t%.1f\t%s\t%s\t%s\n' % r)
data['extent'] = report
with open('%s/candidates_summary.json' % OUT, 'w') as f:
    json.dump(data, f, indent=1)
print('\nflank-track windows written: %d' % len(flank_rows))
