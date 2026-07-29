#!/usr/bin/env python3
"""Sub-subfamily assignment for the 7 candidate loci, following Methods 3.6:
a two-monomer window from the array boundary is aligned (ssearch36 -m 8C) against the
doubled sub-subfamily consensuses of the assigned family only; highest bitscore wins,
and the margin to the second best is recorded."""
import json
import os
import subprocess

BASE = '/staging/Sergei/Lacertidae/Darevskia'
OUT = BASE + '/Satellites/motifs_output/clean/candidates'
SS = '/staging/conda/envs/bioinfo/bin/ssearch36'
SUB = '/tmp/subsubfamily_consensuses.fasta'

GEN = {
    'unp': BASE + '/unp/GCA_032766555.1_D.uni.pHiFi.1.0_genomic.fna',
    'val': BASE + '/val/GCA_034642135.1_Dar_val_2.0_genomic.fna',
    'nai': BASE + '/nai/GCA_034642065.1_Dar_nair_1.0_genomic.fna',
    'unm': BASE + '/unm/GCA_032766585.1_D.uni.mHiFi.1.0_genomic.fna',
    'arm': BASE + '/arm/GCA_034641395.1_Dar_arm_1.0_genomic.fna',
}


def read_fa(p):
    s, n = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith('>'):
            n = l[1:].split()[0]
            s[n] = []
        elif n:
            s[n].append(l)
    return {k: ''.join(v).upper() for k, v in s.items()}


ss = read_fa(SUB)
print('sub-subfamily consensuses: %s' % {k: len(v) for k, v in ss.items()})

# doubled subjects grouped by parent family
FAMOF = {}
for name in ss:
    FAMOF[name] = 'CLsat' + name.replace('CLsat', '').split('.')[0]
for name, seq in ss.items():
    with open('%s/ss_%s.fa' % (OUT, name), 'w') as f:
        f.write('>%s\n%s\n' % (name, seq + seq))

data = json.load(open('%s/candidates_summary.json' % OUT))
summary = data['summary']

for cid, s in sorted(summary.items()):
    sp, scaf, a, b = s['sp'], s['scf'], s['start'] + 1, s['end']
    fam = s['call']
    mono = 146 if fam in ('CLsat1', 'CLsat3') else (145 if fam == 'CLsat2' else 116)
    win = min(2 * mono, b - a + 1)
    r = subprocess.run(['samtools', 'faidx', GEN[sp], '%s:%d-%d' % (scaf, a, a + win - 1)],
                       capture_output=True, text=True)
    q = ''.join(l.strip() for l in r.stdout.splitlines()[1:]).upper()
    with open('%s/_ssq.fa' % OUT, 'w') as f:
        f.write('>q\n%s\n' % q)
    scores = {}
    for name in ss:
        if FAMOF[name] != fam:
            continue
        o = subprocess.run([SS, '-n', '-m', '8C', '-E', '1000', '-b', '1',
                            '%s/_ssq.fa' % OUT, '%s/ss_%s.fa' % (OUT, name)],
                           capture_output=True, text=True).stdout
        best = 0.0
        pid = 0.0
        for line in o.splitlines():
            if line.startswith('#'):
                continue
            p = line.split('\t')
            if len(p) >= 12:
                try:
                    if float(p[11]) > best:
                        best = float(p[11])
                        pid = float(p[2])
                except ValueError:
                    pass
        scores[name] = (best, pid)
    if not scores:
        continue
    order = sorted(scores.items(), key=lambda kv: -kv[1][0])
    top, (tb, tp) = order[0]
    second, sb = (order[1][0], order[1][1][0]) if len(order) > 1 else (None, 0.0)
    delta = round(tb - sb, 1)
    s['subfamily'] = top
    s['ss_bits'] = tb
    s['ss_pid'] = tp
    s['ss_second'] = second
    s['ss_delta'] = delta
    s['ss_confident'] = bool(delta >= 10 or second is None)
    print('  %-3s %-7s -> %-9s %6.1f bits  pid=%.1f  second=%-9s delta=%5.1f  confident=%s'
          % (cid, fam, top, tb, tp, second, delta, s['ss_confident']))

with open('%s/candidates_summary.json' % OUT, 'w') as f:
    json.dump(data, f, indent=1)
print('\nupdated candidates_summary.json')
