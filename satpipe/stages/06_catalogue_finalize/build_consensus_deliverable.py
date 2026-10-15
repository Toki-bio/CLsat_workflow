# -*- coding: utf-8 -*-
"""Reference-free per-assembly consensus comparison for one family, saved as deliverable files
(the manuscript's numbers came from this method: seed=5, 600 monomers per assembly, 300
bootstrap resamples; results reproduce exactly with the same seed).

Usage: build_consensus_deliverable.py <genomes.tsv> <seeds.bed> <consensus.fa> <out_dir>
                                      [--pairs a:b,c:d,...] [--n-monomers 600] [--n-boot 300] [--seed 5]
  genomes.tsv   species_code  genome.fa (indexed); the species order of the tables
  seeds.bed     the family's seed BED (col4 = species|...)
  consensus.fa  the family consensus (one record; its length is the monomer length)
  --pairs       species pairs for the bootstrap CI table; default: every pair. The summary
                compares the FIRST TWO pairs when at least two are given (the published run
                used unp:nai vs unp:val).
Tools: samtools on PATH.
Outputs: <family>_per_assembly_consensus.fa, <family>_consensus_pairwise_identity.tsv,
         <family>_consensus_bootstrap_CI.tsv, <family>_consensus_summary.txt
"""
import argparse
import collections
import itertools
import os
import random
import shutil
import statistics
import subprocess
import sys
import tempfile

ap = argparse.ArgumentParser()
ap.add_argument('genomes')
ap.add_argument('seeds')
ap.add_argument('consensus')
ap.add_argument('out')
ap.add_argument('--pairs', default='')
ap.add_argument('--n-monomers', type=int, default=600)
ap.add_argument('--n-boot', type=int, default=300)
ap.add_argument('--seed', type=int, default=5)
a = ap.parse_args()

SAM = shutil.which('samtools') or sys.exit('samtools not on PATH')
OUT = a.out
os.makedirs(OUT, exist_ok=True)
NMON = a.n_monomers
NBOOT = a.n_boot

GEN = {}
SPS = []
for line in open(a.genomes):
    if not line.strip() or line.startswith('#'):
        continue
    p = line.rstrip('\n').split('\t')
    if len(p) >= 2:
        GEN[p[0]] = p[1]
        SPS.append(p[0])
FAM = None
C1 = ''
for line in open(a.consensus):
    if line.startswith('>'):
        FAM = line[1:].split()[0]
    else:
        C1 += line.strip().upper()
FAM = FAM or 'family'
M = len(C1)
D1 = C1 + C1
TR = str.maketrans('ACGTN', 'TGCAN')
rc = lambda s: s.translate(TR)[::-1]


def ident(x, y):
    n = min(len(x), len(y))
    return sum(1 for i in range(n) if x[i] == y[i]) / n if n else 0.0


seeds = collections.defaultdict(list)
for line in open(a.seeds):
    p = line.rstrip('\n').split('\t')
    if len(p) >= 4:
        sp = p[3].split('|')[0]
        if sp in GEN:
            seeds[sp].append((p[0], int(p[1]), int(p[2])))


def monomers(sp, n_target, seed):
    SL = {}
    for line in open(GEN[sp] + '.fai'):
        q = line.split('\t')
        SL[q[0]] = int(q[1])
    R = [r for r in seeds[sp] if r[2] - r[1] > 5 * M]
    if not R:
        return []
    random.seed(seed)
    regs = []
    for _ in range(n_target * 4):
        scf, x, y = random.choice(R)
        if y > SL.get(scf, 0):
            continue
        p = random.randint(x, y - 3 * M)
        regs.append((scf, p + 1, p + 3 * M))
        if len(regs) >= n_target * 2:
            break
    with tempfile.NamedTemporaryFile('w', suffix='.txt', delete=False) as f:
        for scf, s, e in regs:
            f.write('%s:%d-%d\n' % (scf, s, e))
        rf = f.name
    out = subprocess.run([SAM, 'faidx', GEN[sp], '--region-file', rf],
                         capture_output=True, text=True).stdout
    os.unlink(rf)
    segs = []
    cur = None
    buf = []
    for line in out.splitlines():
        if line.startswith('>'):
            if cur and buf:
                segs.append(''.join(buf).upper())
            cur = line; buf = []
        elif cur:
            buf.append(line.strip())
    if cur and buf:
        segs.append(''.join(buf).upper())
    mons = []
    for sg in segs:
        if len(sg) < 3 * M or sg.count('N'):
            continue
        best = (None, -1, 0)
        for s2 in (sg, rc(sg)):
            for k in range(M):
                sc = ident(s2[:M], D1[k:k + M])
                if sc > best[1]:
                    best = (s2, sc, k)
        s2, sc, k = best
        if sc < 0.55:
            continue
        mon = s2[(M - k) % M:][:M]
        if len(mon) == M:
            mons.append(mon)
        if len(mons) >= n_target:
            break
    return mons


def consensus(ms):
    return ''.join((collections.Counter(m[i] for m in ms if m[i] in 'ACGT').most_common(1)
                    or [('N',)])[0][0] for i in range(M))


POOL = {sp: monomers(sp, NMON, a.seed) for sp in SPS}
SPS = [sp for sp in SPS if POOL[sp]]
if not SPS:
    sys.exit('no monomers sampled for any species (seed BED empty or arrays shorter than 5 monomers)')
CONS = {sp: consensus(POOL[sp]) for sp in SPS}

# --- deliverable 1: per-assembly consensus FASTA ---
with open(f'{OUT}/{FAM}_per_assembly_consensus.fa', 'w') as f:
    for sp in SPS:
        f.write(f'>{sp}_{FAM}_consensus n_monomers={len(POOL[sp])}\n{CONS[sp]}\n')

# --- deliverable 2: full pairwise identity matrix ---
with open(f'{OUT}/{FAM}_consensus_pairwise_identity.tsv', 'w') as f:
    f.write('assembly\t' + '\t'.join(SPS) + '\n')
    for x in SPS:
        row = [f'{100*ident(CONS[x], CONS[y]):.1f}' if x != y else '-' for y in SPS]
        f.write(x + '\t' + '\t'.join(row) + '\n')

# --- deliverable 3: bootstrap CI table ---
if a.pairs:
    PAIRS = [tuple(kp.split(':')) for kp in a.pairs.split(',')]
    for x, y in PAIRS:
        if x not in CONS or y not in CONS:
            sys.exit('pair %s:%s names an assembly without sampled monomers' % (x, y))
else:
    PAIRS = list(itertools.combinations(SPS, 2))
pairs = {pq: [] for pq in PAIRS}
random.seed(a.seed)
for _ in range(NBOOT):
    cons = {}
    for sp in SPS:
        ms = [random.choice(POOL[sp]) for _ in range(len(POOL[sp]))]
        cons[sp] = consensus(ms)
    for (x, y) in pairs:
        pairs[(x, y)].append(100 * ident(cons[x], cons[y]))

with open(f'{OUT}/{FAM}_consensus_bootstrap_CI.tsv', 'w') as f:
    f.write('pair\tmedian_pct_identity\tCI95_low\tCI95_high\n')
    for (x, y), v in pairs.items():
        v.sort()
        f.write(f'{x}-{y}\t{statistics.median(v):.1f}\t{v[int(0.025*NBOOT)]:.1f}\t{v[int(0.975*NBOOT)]:.1f}\n')

with open(f'{OUT}/{FAM}_consensus_summary.txt', 'w') as f:
    f.write(f'Reference-free per-assembly {FAM} consensus comparison\n')
    f.write(f'Method: {NMON} bootstrap-resampled monomers per assembly, seed={a.seed}, {NBOOT} resamples\n\n')
    if len(PAIRS) >= 2:
        (p1, p2) = PAIRS[0], PAIRS[1]
        d = [x - y for x, y in zip(pairs[p1], pairs[p2])]
        d.sort()
        frac = 100.0 * sum(1 for x in d if x > 0) / len(d)
        f.write(f'{p1[0]} consensus vs {p1[1]}: {ident(CONS[p1[0]], CONS[p1[1]])*100:.1f}%\n')
        f.write(f'{p2[0]} consensus vs {p2[1]}: {ident(CONS[p2[0]], CONS[p2[1]])*100:.1f}%\n')
        f.write(f'{p1[0]}-{p1[1]} minus {p2[0]}-{p2[1]}: median {statistics.median(d):+.1f} points, '
                f'95% CI {d[int(0.025*NBOOT)]:+.1f} to {d[int(0.975*NBOOT)]:+.1f}\n')
        f.write(f'fraction of resamples where {p1[0]}-{p1[1]} identity exceeds {p2[0]}-{p2[1]}: {frac:.1f}%\n')

print('wrote deliverables to', OUT)
print(open(f'{OUT}/{FAM}_consensus_summary.txt').read())
