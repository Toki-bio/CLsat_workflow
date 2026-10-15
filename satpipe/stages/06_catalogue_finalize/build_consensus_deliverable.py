# -*- coding: utf-8 -*-
"""Regenerate the reference-free per-assembly CLsat1 consensus comparison and SAVE it as
real deliverable files (this was previously only ever printed to a terminal and never
saved) -- same deterministic method (seed=5, 600 monomers/assembly) already used for the
manuscript's numbers, so results should reproduce exactly."""
import subprocess, random, collections, statistics

BIN = '/staging/conda/envs/bioinfo/bin'
SD = ('/staging/Sergei/Lacertidae/Darevskia/Satellites/motifs_output/clean/'
      'cross_clusters/selected/subclusters/seeds')
OUT = '/staging/tmp/clsat_v6/consensus_deliverable'
import os
os.makedirs(OUT, exist_ok=True)

M = 146
NMON = 600
NBOOT = 300

GEN = {}
for line in open(SD + '/genomes.tsv'):
    p = line.rstrip('\n').split('\t')
    if len(p) >= 2:
        GEN[p[0]] = p[1]
C1 = ''.join(l.strip() for l in open(SD + '/CLsat1.cons.fa') if not l.startswith('>')).upper()
D1 = C1 + C1
TR = str.maketrans('ACGTN', 'TGCAN')
rc = lambda s: s.translate(TR)[::-1]
def ident(a, b):
    n = min(len(a), len(b))
    return sum(1 for i in range(n) if a[i] == b[i]) / n if n else 0.0

seeds = collections.defaultdict(list)
for line in open(SD + '/CLsat1.bed'):
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
    random.seed(seed)
    regs = []
    for _ in range(n_target * 4):
        scf, a, b = random.choice(R)
        if b > SL.get(scf, 0):
            continue
        p = random.randint(a, b - 3 * M)
        regs.append((scf, p + 1, p + 3 * M))
        if len(regs) >= n_target * 2:
            break
    rf = '/tmp/bs.txt'
    with open(rf, 'w') as f:
        for scf, s, e in regs:
            f.write('%s:%d-%d\n' % (scf, s, e))
    out = subprocess.run([BIN + '/samtools', 'faidx', GEN[sp], '--region-file', rf],
                         capture_output=True, text=True).stdout
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

SPS = ['arm', 'val', 'mix', 'unm', 'nai', 'unp']
POOL = {sp: monomers(sp, NMON, 5) for sp in SPS}
CONS = {sp: consensus(POOL[sp]) for sp in SPS}

# --- deliverable 1: per-assembly consensus FASTA ---
with open(f'{OUT}/CLsat1_per_assembly_consensus.fa', 'w') as f:
    for sp in SPS:
        f.write(f'>{sp}_CLsat1_consensus n_monomers={len(POOL[sp])}\n{CONS[sp]}\n')

# --- deliverable 2: full pairwise identity matrix ---
with open(f'{OUT}/CLsat1_consensus_pairwise_identity.tsv', 'w') as f:
    f.write('assembly\t' + '\t'.join(SPS) + '\n')
    for a in SPS:
        row = [f'{100*ident(CONS[a], CONS[b]):.1f}' if a != b else '-' for b in SPS]
        f.write(a + '\t' + '\t'.join(row) + '\n')

# --- deliverable 3: bootstrap CI table (regenerated, same method/seed as manuscript numbers) ---
pairs = {('unp', 'nai'): [], ('unp', 'val'): [], ('unp', 'unm'): [], ('unp', 'arm'): [],
         ('unm', 'nai'): [], ('unm', 'val'): []}
random.seed(5)
for _ in range(NBOOT):
    cons = {}
    for sp in SPS:
        ms = [random.choice(POOL[sp]) for _ in range(len(POOL[sp]))]
        cons[sp] = consensus(ms)
    for (a, b) in pairs:
        pairs[(a, b)].append(100 * ident(cons[a], cons[b]))

with open(f'{OUT}/CLsat1_consensus_bootstrap_CI.tsv', 'w') as f:
    f.write('pair\tmedian_pct_identity\tCI95_low\tCI95_high\n')
    for (a, b), v in pairs.items():
        v.sort()
        f.write(f'{a}-{b}\t{statistics.median(v):.1f}\t{v[int(0.025*NBOOT)]:.1f}\t{v[int(0.975*NBOOT)]:.1f}\n')

d = [x - y for x, y in zip(pairs[('unp', 'nai')], pairs[('unp', 'val')])]
d.sort()
frac = 100.0 * sum(1 for x in d if x > 0) / len(d)
with open(f'{OUT}/CLsat1_consensus_summary.txt', 'w') as f:
    f.write('Reference-free per-assembly CLsat1 consensus comparison\n')
    f.write(f'Method: {NMON} bootstrap-resampled monomers per assembly, seed=5, {NBOOT} resamples\n\n')
    f.write(f'unp consensus vs nai: {ident(CONS["unp"], CONS["nai"])*100:.1f}%\n')
    f.write(f'unp consensus vs val: {ident(CONS["unp"], CONS["val"])*100:.1f}%\n')
    f.write(f'unp-nai minus unp-val: median {statistics.median(d):+.1f} points, '
            f'95% CI {d[int(0.025*NBOOT)]:+.1f} to {d[int(0.975*NBOOT)]:+.1f}\n')
    f.write(f'fraction of resamples where unp is closer to nai than to val: {frac:.1f}%\n')

print('wrote deliverables to', OUT)
print(open(f'{OUT}/CLsat1_consensus_summary.txt').read())
