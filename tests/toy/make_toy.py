#!/usr/bin/env python3
"""Make a small synthetic satellite-DNA project for end-to-end tests.

Two fake genomes (codes spa, spb) with planted tandem arrays of the four real CLsat consensuses
(146/145/146/116 bp) in random background sequence, including:
  * ordinary arrays of 12-60 mutated copies per family and species,
  * one composite array (CLsat1 block followed by a CLsat3 block) to exercise segmentation,
  * one orthologous pair: the same unique 1.5 kb flanks around an array in spa and in spb.
Also writes genomes.tsv, families.tsv (family mono consensus), members.tsv, a subx/ directory
of DOUBLED sub-subfamily consensuses, and refs.rot (Step-5 rotated references for Step 8).

Usage: make_toy.py <out_dir> <consensus_dir> [--seed 7] [--size 1200000]
  consensus_dir  holds CLsat1.cons.fa .. CLsat4.cons.fa (satpipe/stages/02_seed_expansion_scan/consensus)
Deterministic for a given seed. Needs only the standard library.
"""
import argparse
import os
import random

ap = argparse.ArgumentParser()
ap.add_argument('out')
ap.add_argument('consensus_dir')
ap.add_argument('--seed', type=int, default=7)
ap.add_argument('--size', type=int, default=1200000, help='background length per genome, bp')
ap.add_argument('--dimer', action='store_true',
                help='also plant one 2-unit higher-order CLsat4 array per genome (TRF reports it at 232; Step 3b test). '
                     'Off by default so the equivalence toy is unchanged.')
ap.add_argument('--same-names', action='store_true',
                help='give both genomes the SAME contig names (chr1, scaf2) and plant a twin pair of arrays of different '
                     'families at the same coordinates in both, so any step that matches by contig name without the genome '
                     'mixes or mislabels them (the split-assembly case: hifiasm h1tg/h2tg names in two animals)')
a = ap.parse_args()
random.seed(a.seed)
os.makedirs(a.out, exist_ok=True)
os.makedirs(os.path.join(a.out, 'genomes'), exist_ok=True)
os.makedirs(os.path.join(a.out, 'subx'), exist_ok=True)

FAMS = ['CLsat1', 'CLsat2', 'CLsat3', 'CLsat4']
CONS = {}
for fam in FAMS:
    CONS[fam] = ''.join(l.strip() for l in open(os.path.join(a.consensus_dir, fam + '.cons.fa'))
                        if not l.startswith('>')).upper().replace('N', 'A')
MONO = {f: len(CONS[f]) for f in FAMS}
BASES = 'ACGT'


def rand_seq(n):
    return ''.join(random.choice(BASES) for _ in range(n))


def mutate(s, rate):
    out = []
    for c in s:
        r = random.random()
        if r < rate * 0.8:
            out.append(random.choice([b for b in BASES if b != c]))
        elif r < rate * 0.9:
            continue                      # deletion
        elif r < rate:
            out.append(c); out.append(random.choice(BASES))   # insertion
        else:
            out.append(c)
    return ''.join(out)


def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def array(fam, copies, rate, phase=0, minus=False):
    c = CONS[fam]
    c = c[phase:] + c[:phase]
    s = ''.join(mutate(c, rate) for _ in range(copies))
    return rc(s) if minus else s


def minrot(s):
    d = s + s
    return min(d[i:i + len(s)] for i in range(len(s)))


# ---- planted features --------------------------------------------------------------------
# each genome: one chromosome-sized scaffold (so expansion has room) + a short second scaffold
plan = {'spa': [], 'spb': []}
for sp in plan:
    for fam in FAMS:
        for k in range(3):
            plan[sp].append(('array', fam, random.randint(12, 60), random.uniform(0.02, 0.06),
                             random.randrange(MONO[fam]), random.random() < 0.3))
    plan[sp].append(('composite',))
    if a.dimer:
        plan[sp].append(('dimer',))
ortho_L = rand_seq(1500)
ortho_R = rand_seq(1500)

truth = open(os.path.join(a.out, 'truth.tsv'), 'w')
truth.write('species\tscaffold\tstart0\tend\tkind\tfamily\tcopies\n')
with open(os.path.join(a.out, 'genomes.tsv'), 'w') as gt:
    for sp in ('spa', 'spb'):
        feats = plan[sp][:]
        random.shuffle(feats)
        parts = []
        pos = 0
        gap = a.size // (len(feats) + 2)
        chrom = 'chr1' if a.same_names else f'{sp}_chr1'
        if a.same_names:
            # twin arrays at identical coordinates in both genomes, different families (CLsat1 in spa, CLsat4 in spb)
            random.seed(a.seed + 1000); bg = rand_seq(20000); random.seed(a.seed + (1 if sp == 'spa' else 2))
            tw = array('CLsat1' if sp == 'spa' else 'CLsat4', 40, 0.03)
            parts.append(bg); pos += len(bg)
            truth.write(f'{sp}\t{chrom}\t{pos}\t{pos+len(tw)}\ttwin\t{"CLsat1" if sp == "spa" else "CLsat4"}\t40\n')
            parts.append(tw); pos += len(tw)
        for feat in feats:
            bg = rand_seq(gap)
            parts.append(bg); pos += len(bg)
            if feat[0] == 'array':
                _, fam, copies, rate, phase, minus = feat
                s = array(fam, copies, rate, phase, minus)
                truth.write(f'{sp}\t{chrom}\t{pos}\t{pos+len(s)}\tarray{"-" if minus else "+"}\t{fam}\t{copies}\n')
                parts.append(s); pos += len(s)
            elif feat[0] == 'dimer':
                # a 2-unit higher-order array: CLsat4 unit A and its 4%-variant B alternate (A B A B ...), so the
                # 232-letter pair repeats better than the 116-letter unit and TRF reports period 232 (Step 3b test)
                A = CONS['CLsat4']; B = mutate(A, 0.04)
                s = ''.join(mutate(A + B, 0.005) for _ in range(30))
                truth.write(f'{sp}\t{chrom}\t{pos}\t{pos+len(s)}\tdimer\tCLsat4\t60\n')
                parts.append(s); pos += len(s)
            else:
                s1 = array('CLsat1', 40, 0.03); s3 = array('CLsat3', 40, 0.03)
                truth.write(f'{sp}\t{chrom}\t{pos}\t{pos+len(s1)}\tcomposite\tCLsat1\t40\n')
                truth.write(f'{sp}\t{chrom}\t{pos+len(s1)}\t{pos+len(s1)+len(s3)}\tcomposite\tCLsat3\t40\n')
                parts.append(s1 + s3); pos += len(s1) + len(s3)
        # orthologous pair: shared unique flanks around a CLsat2 array
        bg = rand_seq(gap); parts.append(bg); pos += len(bg)
        s = array('CLsat2', 25, 0.03)
        parts.append(ortho_L); pos += len(ortho_L)
        truth.write(f'{sp}\t{chrom}\t{pos}\t{pos+len(s)}\tortholog\tCLsat2\t25\n')
        parts.append(s); pos += len(s)
        parts.append(ortho_R); pos += len(ortho_R)
        parts.append(rand_seq(gap))
        seq = ''.join(parts)
        fa = os.path.join(a.out, 'genomes', f'{sp}.fna')
        with open(fa, 'w') as f:
            f.write(f'>{chrom}\n')
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + '\n')
            small = rand_seq(20000)
            f.write('>scaf2\n' if a.same_names else f'>{sp}_scaf2\n')
            for i in range(0, len(small), 80):
                f.write(small[i:i + 80] + '\n')
        gt.write(f'{sp}\t{os.path.abspath(fa)}\n')
truth.close()

# ---- reference files -------------------------------------------------------------------------
with open(os.path.join(a.out, 'families.tsv'), 'w') as f:
    for fam in FAMS:
        f.write(f'{fam}\t{MONO[fam]}\t{os.path.abspath(os.path.join(a.consensus_dir, fam + ".cons.fa"))}\n')
with open(os.path.join(a.out, 'refs.rot'), 'w') as f:        # Step-5 rotation of each consensus
    for fam in FAMS:
        f.write(f'>{fam}\n{minrot(CONS[fam])}\n')
members = {'CLsat1': ['CLsat1.1', 'CLsat1.2'], 'CLsat2': ['CLsat2'], 'CLsat3': ['CLsat3.1', 'CLsat3.2'], 'CLsat4': ['CLsat4']}
with open(os.path.join(a.out, 'members.tsv'), 'w') as f:
    for fam, subs in members.items():
        for sub in subs:
            f.write(f'{fam}\t{sub}\n')
            base = CONS[fam] if sub.endswith('.1') or '.' not in sub else mutate(CONS[fam], 0.10)
            with open(os.path.join(a.out, 'subx', sub + '.fa'), 'w') as g:   # DOUBLED, as subx/ is
                g.write(f'>{sub}\n{base}{base}\n')
print('toy project written to', a.out)
