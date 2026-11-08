"""Toy genomes with known answers for the stage-19 refinements --repeat-mask and --cut-at-locus (author 2026-11-08).
Two genomes, one satellite family (150-letter unit x 20). Locus names = case + copy + genome (a / b).
  F  L flank = 300 unique letters + 300 letters of a dispersed repeat (5 other copies per genome);
     R flank = 600 letters of that repeat                                -> rules: no pair; mask (600): pair via L
  G  three loci per genome inside a carrier: 1,200 carrier letters on each side of the array, the same carrier at
     all three loci (copied together with the array); unique sequence beyond      -> 600 letters: no pair;
     3,000 letters + mask + cut: pairs G1, G2, G3 (anchored outside the carrier)
  H  two loci 1,000 letters apart; their outer flanks are 3,000 letters of a 3-copy repeat, the gap is unique
     -> pairs H1, H2 via the gap in every mode; with --cut-at-locus the gap is split at its middle (no letter used
     by two flanks, no flank reaching into the other locus)
  I  each flank starts with 3 worn satellite units (12 % changed: the array continues past the catalogue edge),
     then unique sequence                                                  -> 600 letters: no pair; long: pair
Writes <out>/genomes/*.fa, genomes.tsv, families.tsv, cons.fa, locus_master.tsv, expected.json."""
import json
import os
import random
import sys

out = sys.argv[1]
R = random.Random(20261109)
ACGT = 'ACGT'


def rnd(n):
    return ''.join(R.choice(ACGT) for _ in range(n))


def mut(s, rate):
    s = list(s)
    for i in range(len(s)):
        if R.random() < rate:
            s[i] = R.choice([c for c in ACGT if c != s[i]])
    return ''.join(s)


UNIT = rnd(150)
SAT = ''.join(mut(UNIT, 0.03) for _ in range(20))
REP = rnd(900)                     # dispersed repeat (case F)
CL, CR = rnd(1200), rnd(1200)      # carrier (case G)
R3 = rnd(3000)                     # repeat at H's outer sides
FU = rnd(300)                      # F's unique part
GX = [(rnd(1800), rnd(1800)) for _ in range(3)]
HG = rnd(1000)
IU = (rnd(2600), rnd(2600))
genome = {'spa': [], 'spb': []}
loci = []


def pos(g):
    return sum(len(x) for x in genome[g])


def add(g, s):
    genome[g].append(s)


def array(g, name):
    off = pos(g)
    add(g, SAT)
    loci.append((name, 'scf1', off, off + len(SAT), g, 'SatX', 150))


for g, r in (('spa', 0.0), ('spb', 0.05)):
    m = (lambda s, r=r: mut(s, r)) if r else (lambda s: s)
    t = g[-1]
    add(g, rnd(5000))
    # F
    add(g, m(FU) + mut(REP[:300], 0.02)); array(g, 'F' + t); add(g, mut(REP[300:], 0.02)); add(g, rnd(5000))
    # G1..G3
    for n in range(3):
        add(g, m(GX[n][0]) + mut(CL, 0.02)); array(g, 'G%d%s' % (n + 1, t)); add(g, mut(CR, 0.02) + m(GX[n][1]))
        add(g, rnd(5000))
    # H1, gap, H2
    add(g, mut(R3, 0.02)); array(g, 'H1' + t); add(g, m(HG)); array(g, 'H2' + t); add(g, mut(R3, 0.02))
    add(g, rnd(5000))
    # I
    add(g, m(IU[0]) + ''.join(mut(UNIT, 0.12) for _ in range(3))); array(g, 'I' + t)
    add(g, ''.join(mut(UNIT, 0.12) for _ in range(3)) + m(IU[1])); add(g, rnd(5000))
    # other copies of the repeats
    for _ in range(5):
        add(g, mut(REP, 0.03)); add(g, rnd(3000))
    add(g, mut(R3, 0.02)); add(g, rnd(3000))

os.makedirs(os.path.join(out, 'genomes'), exist_ok=True)
with open(os.path.join(out, 'genomes.tsv'), 'w') as gt:
    for g, parts in genome.items():
        p = os.path.join(out, 'genomes', g + '.fa')
        s = ''.join(parts)
        with open(p, 'w') as f:
            f.write('>scf1\n' + '\n'.join(s[i:i + 80] for i in range(0, len(s), 80)) + '\n')
        gt.write('%s\t%s\n' % (g, os.path.abspath(p)))
with open(os.path.join(out, 'cons.fa'), 'w') as f:
    f.write('>SatX\n%s\n' % UNIT)
with open(os.path.join(out, 'families.tsv'), 'w') as f:
    f.write('SatX\t150\t%s\n' % os.path.abspath(os.path.join(out, 'cons.fa')))
with open(os.path.join(out, 'locus_master.tsv'), 'w') as f:
    for x in loci:
        f.write('\t'.join(map(str, x)) + '\n')


def P(*c):
    return [[x + 'a', x + 'b'] for x in c]


exp = {'rules': P('H1', 'H2'),
       'mask600': P('F', 'H1', 'H2'),
       'long': P('F', 'G1', 'G2', 'G3', 'H1', 'H2', 'I'),
       'loci': {x[0]: [x[2], x[3]] for x in loci}}
json.dump(exp, open(os.path.join(out, 'expected.json'), 'w'), indent=1)
print('toy written:', out, len(loci), 'loci')
