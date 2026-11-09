"""Toy genomes with planted answers for satpipe.ortho_stack. Usage: make_toy.py <dir>

Genomes A (code a), B (code b), C (code c); one family F, 100-letter consensus.
Group 1 (loci A1, B1, C1): flanks 5,000 letters shared (1 % substitutions per genome); array = 10 monomer variants
  v1..v10 (each the consensus with 8 own substitutions), starting at consensus position 31 (first monomer 70 letters).
  B1 lacks v6 and v7. C1: the whole locus is reverse-complemented in genome C; its left flank also carries a
  200-letter insertion next to nothing in A/B.
Group 2 (A2, B2): homogenised arrays, A2 8 identical monomers, B2 6 -> no clear monomer match, all left-shifted.
A1's left-flank letters 201..300 from the array edge are marked masked (repeat) in flanks.json.
Writes genomes/*.fa, flanks.json (stage-19 format, start 0-based, flank runs 1-based), ortho_pairs.tsv,
cons.fa, expected.json."""
import json
import os
import random
import sys

R = random.Random(7)
D = sys.argv[1]
os.makedirs(D + '/genomes', exist_ok=True)
ACGT = 'ACGT'
comp = str.maketrans('ACGT', 'TGCA')
rc = lambda s: s.translate(comp)[::-1]
rnd = lambda n: ''.join(R.choice(ACGT) for _ in range(n))


def mut(s, rate):
    s = list(s)
    for i in range(len(s)):
        if R.random() < rate:
            s[i] = R.choice([x for x in ACGT if x != s[i]])
    return ''.join(s)


def subst(s, pos):
    s = list(s)
    for p in pos:
        s[p] = R.choice([x for x in ACGT if x != s[p]])
    return ''.join(s)


cons = rnd(100)
var = [subst(cons, R.sample(range(100), 8)) for _ in range(10)]
T = ''.join(var) + var[0]
arrA = T[30:1000]                       # 970 letters, monomer starts at 70, 170, ... 970 -> 10 segments
arrB = ''.join(v for i, v in enumerate(var) if i not in (5, 6)) + var[0]
arrB = arrB[30:800]                     # 8 variants
LF1, RF1 = rnd(5000), rnd(5000)
homo = cons
LF2, RF2 = rnd(5000), rnd(5000)
pad = lambda: rnd(3000)
loci, seqs = {}, {'a': '', 'b': '', 'c': ''}


def place(code, lid, lf, arr, rf, reverse=False):
    g = seqs[code] + pad()
    if reverse:
        block = rc(lf + arr + rf)
        s = len(g) + len(rf)            # in the genome, rc(rf) comes first
        e = s + len(arr)
        lflank, rflank = rc(rf), rc(lf)
    else:
        block = lf + arr + rf
        s = len(g) + len(lf); e = s + len(arr)
        lflank, rflank = lf, rf
    seqs[code] = g + block
    loci[lid] = {'sp': code, 'scf': 'ctg_' + code, 'start': s, 'end': e, 'subtype': 'F', 'mono': 100,
                 'length': e - s, 'n_monomers': round((e - s) / 100),
                 'LF': {'start': s - 600 + 1, 'end': s, 'len': 600, 'seq': lflank[-600:], 'truncated': False, 'class': 'unique'},
                 'RF': {'start': e + 1, 'end': e + 600, 'len': 600, 'seq': rflank[:600], 'truncated': False, 'class': 'unique'}}


place('a', 'L1', LF1, arrA, RF1)
place('b', 'L2', mut(LF1, 0.01), arrB, mut(RF1, 0.01))
lfc = mut(LF1, 0.01)
lfc = lfc[:2000] + rnd(200) + lfc[2000:]
place('c', 'L3', lfc, mut(arrA, 0.003), mut(RF1, 0.01), reverse=True)
place('a', 'L4', LF2, homo * 8, RF2)
# Group 3 (A L6, B L7): neighbouring monomers differ by 2 letters only (each variant = consensus + 1 own letter);
# A's 7th monomer is a second copy of its 3rd. B = A. B's 3rd and 7th have two exact partners -> not clear; the
# other 10 are clear (0 letters against >= 2).
pos = R.sample(range(100), 12)
near = [subst(cons, [pos[k]]) for k in range(12)]
near[6] = near[2]
LF3, RF3 = rnd(5000), rnd(5000)
place('a', 'L6', LF3, ''.join(near), RF3)
place('b', 'L5', mut(LF2, 0.01), homo * 6, mut(RF2, 0.01))
place('b', 'L7', mut(LF3, 0.01), ''.join(near), mut(RF3, 0.01))
for c in seqs:
    seqs[c] += pad()
    with open('%s/genomes/%s.fa' % (D, c), 'w') as f:
        f.write('>ctg_%s\n' % c)
        for i in range(0, len(seqs[c]), 80):
            f.write(seqs[c][i:i + 80] + '\n')
s = loci['L1']['start']
loci['L1']['LF']['mask'] = {'repeat': [[s - 299, s - 200]]}      # 1-based genomic: edge indices 200..299
json.dump(loci, open(D + '/flanks.json', 'w'))
with open(D + '/ortho_pairs.tsv', 'w') as f:
    f.write('locusA\tlocusB\tspA\tspB\tpid\taln\tbits\n')
    for a, b in (('L1', 'L2'), ('L1', 'L3'), ('L4', 'L5'), ('L6', 'L7')):
        f.write('%s\t%s\t%s\t%s\t99.0\t600\t1000.0\n' % (a, b, loci[a]['sp'], loci[b]['sp']))
open(D + '/cons.fa', 'w').write('>F\n%s\n' % cons)
json.dump({'G001': {'rows': ['L1', 'L2', 'L3'], 'ori': ['+', '+', '-'],
                    'ref_bounds': list(range(70, 970, 100)), 'n_ref': 10,
                    # ref segment index of each B monomer: B = v1(part) v2..v5 v8 v9 v10 (v1 tail)
                    'b_cols': [0, 1, 2, 3, 4, 7, 8, 9], 'c_cols': list(range(10)),
                    'mask_window': 2},
           'G002': {'rows': ['L4', 'L5'], 'anchored_b': 0, 'b_n': 6},
           'G003': {'rows': ['L6', 'L7'], 'b_anchored': [1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1]}}, open(D + '/expected.json', 'w'))
print('toy written')
