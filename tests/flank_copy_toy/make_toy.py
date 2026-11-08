"""Toy genomes with known answers for the stage-19 single-copy anchor rules (author 2026-11-08).
Two genomes, one satellite family (150-letter unit x 20). Cases (locus names = case letter + genome):
  A  ordinary ortholog: unique flanks, copies in spb at 95 %                    -> pair Aa-Ab (rules on and off)
  B  left flank of Aa duplicated elsewhere in spa at 90 % (not at a locus)       -> Ba L 'repeated'; pair Ba-Bb via R
  C  spb holds two loci C1b, C2b whose flanks are independent 15 % mutants of
     Ca's flanks (C1b vs C2b ~73 %: each single copy) -> Ca hits both ~equally   -> no pair with rules on
  D  left flank carries a 450-letter tandem repeat (60-letter unit)              -> Da L 'tandem'; pair Da-Db via R
  E  right flank of Ea has a 75 % copy elsewhere in spa                         -> stays single; pair Ea-Eb
Writes <out>/genomes/*.fa, genomes.tsv, families.tsv, cons.fa, locus_master.tsv, expected.json."""
import json
import os
import random
import sys

out = sys.argv[1]
R = random.Random(20261108)
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
FL = {k: (rnd(600), rnd(600)) for k in 'ABCDE'}
tandem = rnd(60) * 8    # 480 letters; a 60-letter unit passes the low-complexity rule, so only the TRF mask catches it
FL['D'] = (FL['D'][0][:70] + tandem[:450] + FL['D'][0][520:600], FL['D'][1])
assert len(FL['D'][0]) == 600

genome = {'spa': [], 'spb': []}
loci = []


def spacer():
    return rnd(5000)


def put(g, name, left, right):
    seq = genome[g]
    seq.append(spacer())
    off = sum(len(x) for x in seq) + len(left)
    seq.append(left + SAT + right)
    loci.append((name, 'scf1', off, off + len(SAT), g, 'SatX', 150))


put('spa', 'Aa', *FL['A'])
put('spa', 'Ba', *FL['B'])
put('spa', 'Ca', *FL['C'])
put('spa', 'Da', *FL['D'])
put('spa', 'Ea', *FL['E'])
genome['spa'] += [spacer(), mut(FL['B'][0], 0.10), spacer(), mut(FL['E'][1], 0.25), spacer()]
put('spb', 'Ab', mut(FL['A'][0], 0.05), mut(FL['A'][1], 0.05))
put('spb', 'Bb', mut(FL['B'][0], 0.05), mut(FL['B'][1], 0.05))
put('spb', 'C1b', mut(FL['C'][0], 0.15), mut(FL['C'][1], 0.15))
put('spb', 'C2b', mut(FL['C'][0], 0.15), mut(FL['C'][1], 0.15))
put('spb', 'Db', mut(FL['D'][0][:70], 0.05) + rnd(450) + mut(FL['D'][0][520:], 0.05), mut(FL['D'][1], 0.05))
put('spb', 'Eb', mut(FL['E'][0], 0.05), mut(FL['E'][1], 0.05))
genome['spb'].append(spacer())

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
exp = {
    'on': {'pairs': [['Aa', 'Ab'], ['Ba', 'Bb'], ['Da', 'Db'], ['Ea', 'Eb']],
           'class': {'Ba|L': 'repeated', 'Da|L': 'tandem', 'Ea|R': 'unique', 'C1b|L': 'unique', 'C2b|L': 'unique',
                     'Aa|L': 'unique', 'Ab|R': 'unique'}},
    'off_must_include': [['Aa', 'Ab'], ['Ea', 'Eb']],
}
json.dump(exp, open(os.path.join(out, 'expected.json'), 'w'), indent=1)
print('toy written:', out, len(loci), 'loci')
