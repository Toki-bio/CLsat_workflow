"""Sensitivity of the own-genome copy search (stage 19 single_copy rule) at the decision identity.
A 50 Mb random genome holds 100 random 600-letter flanks once each (the 'self' place) and, for 80 of them, one more
copy with substitutions and small indels: 20 each at 90, 85, 82, 80 % identity; 20 flanks get a 75 % copy (must
NOT count as a copy at the 80 % rule). Each search method must report, per flank, the number of OTHER places with
>= 80 % identity over >= 200 letters.
Usage: bench_copy_search.py <work_dir> <threads>"""
import os
import random
import subprocess
import sys
import time

W, T = sys.argv[1], sys.argv[2]
os.makedirs(W, exist_ok=True)
R = random.Random(7)
A = 'ACGT'


def rnd(n):
    return ''.join(R.choice(A) for _ in range(n))


def mut(s, ident):
    out = []
    for c in s:
        x = R.random()
        if x < (1 - ident) * 0.9:
            out.append(R.choice([b for b in A if b != c]))
        elif x < (1 - ident):                      # 10 % of differences are 1-letter indels
            if R.random() < 0.5:
                continue
            out.append(c + R.choice(A))
        else:
            out.append(c)
    return ''.join(out)


fl = [rnd(600) for _ in range(100)]
level = {}
for i in range(100):
    level[i] = [0.90, 0.85, 0.82, 0.80, 0.75][i // 20] if i < 100 else None
pieces, places = [], {}
pos = 0


def add(s):
    global pos
    pieces.append(s); pos += len(s)


add(rnd(200000))
for i, f in enumerate(fl):
    places[i] = (pos + 1, pos + 600); add(f); add(rnd(R.randint(100000, 300000)))
for i, f in enumerate(fl):
    c = mut(f, level[i])
    if R.random() < 0.5:
        c = c[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    add(c); add(rnd(R.randint(100000, 300000)))
g = ''.join(pieces)
gp = os.path.join(W, 'bench.fa')
with open(gp, 'w') as f:
    f.write('>chr1\n' + '\n'.join(g[i:i + 80] for i in range(0, len(g), 80)) + '\n')
qp = os.path.join(W, 'q.fa')
with open(qp, 'w') as f:
    for i, s in enumerate(fl):
        f.write('>f%d\n%s\n' % (i, s))
print('genome %.1f Mb' % (len(g) / 1e6), flush=True)


def count_paf(out):
    n = {}
    for l in out.splitlines():
        p = l.split('\t')
        q = int(p[0][1:]); ts, te, m, b = int(p[7]), int(p[8]), int(p[9]), int(p[10])
        s, e = places[q]
        if ts < e and te > s - 1:
            continue
        if b >= 200 and m / b >= 0.80:
            n[q] = n.get(q, 0) + 1
    return n


def count_blast(out):
    n = {}
    for l in out.splitlines():
        p = l.split('\t')
        q = int(p[0][1:]); pid, aln, ss, se = float(p[2]), int(p[3]), int(p[8]), int(p[9])
        lo, hi = min(ss, se), max(ss, se)
        s, e = places[q]
        if lo <= e and hi >= s:
            continue
        if aln >= 200 and pid >= 80:
            n[q] = n.get(q, 0) + 1
    return n


methods = {
    'mm2 asm20 -p0.2 (night T3 / first fix)': ['minimap2', '-x', 'asm20', '-c', '-N', '200', '-p', '0.2'],
    'mm2 asm20 -p0': ['minimap2', '-x', 'asm20', '-c', '-N', '200', '-p', '0'],
    'mm2 k15 w5 -p0': ['minimap2', '-k15', '-w5', '-s50', '-m30', '-A2', '-B4', '-O4,24', '-E2,1', '-c', '-N', '200', '-p', '0'],
    'mm2 k13 w4 -p0': ['minimap2', '-k13', '-w4', '-s40', '-m25', '-A2', '-B4', '-O4,24', '-E2,1', '-c', '-N', '200', '-p', '0'],
}
res = {}
for name, cmd in methods.items():
    t = time.time()
    o = subprocess.run(cmd + ['-t', T, gp, qp], capture_output=True, text=True).stdout
    res[name] = (count_paf(o), time.time() - t)
subprocess.run(['makeblastdb', '-in', gp, '-dbtype', 'nucl', '-out', os.path.join(W, 'bench')], capture_output=True)
for name, task in (('blastn word 11', ['-task', 'blastn', '-word_size', '11']), ('blastn dc-megablast', ['-task', 'dc-megablast'])):
    t = time.time()
    o = subprocess.run(['blastn', '-query', qp, '-db', os.path.join(W, 'bench'), '-outfmt', '6', '-evalue', '1e-10',
                        '-max_target_seqs', '1000', '-num_threads', T] + task, capture_output=True, text=True).stdout
    res[name] = (count_blast(o), time.time() - t)
print('method\tfound 90%\t85%\t82%\t80%\tfalse at 75%\tseconds')
for name, (n, sec) in res.items():
    f = lambda lv: sum(1 for i in range(100) if level[i] == lv and n.get(i, 0) >= 1)
    print('%s\t%d/20\t%d/20\t%d/20\t%d/20\t%d/20\t%.0f' % (name, f(0.90), f(0.85), f(0.82), f(0.80), f(0.75), sec), flush=True)
