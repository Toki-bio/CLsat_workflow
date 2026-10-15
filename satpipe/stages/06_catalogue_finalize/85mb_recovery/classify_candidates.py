# -*- coding: utf-8 -*-
"""Classify chunked-rescan TRF candidates (period near a family monomer length, copies >= N)
against the family consensus sequences by best rotation-aware alignment identity. Keep only
genuine matches; discard coincidental period matches to unrelated repeat families.

Usage: classify_candidates.py <scaf.trf.dedup.dat> <families.tsv> <out.tsv>
                              [--period-tol 2] [--min-copies 10] [--min-identity 0.70]
  families.tsv  family  mono  consensus.fa
A candidate's period must lie within +/- period-tol of some family's mono (the published
run used 114-118 and 143-148 for monos 116, 145, 146 = tol 2).
"""
import argparse
import sys

ap = argparse.ArgumentParser()
ap.add_argument('dat')
ap.add_argument('families')
ap.add_argument('out')
ap.add_argument('--period-tol', type=int, default=2)
ap.add_argument('--min-copies', type=float, default=10.0)
ap.add_argument('--min-identity', type=float, default=0.70)
a = ap.parse_args()

cons = {}
monos = []
for line in open(a.families):
    if not line.strip() or line.startswith('#'):
        continue
    p = line.rstrip('\n').split('\t')
    fam, mono, fa = p[0], int(p[1]), p[2]
    cons[fam] = ''.join(l.strip() for l in open(fa) if not l.startswith('>')).upper()
    monos.append(mono)
if not cons:
    sys.exit('no families in %s' % a.families)


def best_rotation_identity(motif, consensus):
    """Best identity of motif against any rotation of consensus (and its revcomp), same-length
    sliding compare -- handles TRF's arbitrary phase and strand for the reported motif."""
    def revcomp(s):
        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(comp.get(c, 'N') for c in reversed(s))

    best = 0.0
    for base in (consensus, revcomp(consensus)):
        doubled = base + base
        L = len(base)
        m = motif[:L] if len(motif) >= L else motif + motif[:L - len(motif)]
        for offset in range(L):
            window = doubled[offset:offset + L]
            if len(window) != len(m):
                continue
            matches = sum(1 for x, y in zip(window, m) if x == y)
            ident = matches / L
            if ident > best:
                best = ident
    return best


def period_ok(period):
    return any(abs(period - m) <= a.period_tol for m in monos)


candidates = []
for line in open(a.dat):
    p = line.split()
    if len(p) < 14:
        continue
    s, e, period, copies = int(p[0]), int(p[1]), int(p[2]), float(p[3])
    motif = p[13].upper()
    if not period_ok(period):
        continue
    if copies < a.min_copies:
        continue
    candidates.append((s, e, period, copies, motif))

print(f'{len(candidates)} candidates with period within {a.period_tol} of a family mono and copies>={a.min_copies:g}\n')
kept = []
for s, e, period, copies, motif in candidates:
    scores = {fam: best_rotation_identity(motif, c) for fam, c in cons.items()}
    best_fam = max(scores, key=scores.get)
    best_score = scores[best_fam]
    verdict = 'KEEP' if best_score >= a.min_identity else 'discard'
    print(f'{s}-{e}  period={period}  copies={copies:.1f}  best={best_fam} ({best_score:.2f})  {verdict}')
    if verdict == 'KEEP':
        kept.append((s, e, period, copies, best_fam, best_score))

print(f'\n{len(kept)} of {len(candidates)} candidates kept as genuine matches (>={a.min_identity:.0%} identity)')
with open(a.out, 'w') as out:
    out.write('start\tend\tperiod\tcopies\tfamily\tidentity\n')
    for s, e, period, copies, fam, score in kept:
        out.write(f'{s}\t{e}\t{period}\t{copies}\t{fam}\t{score:.3f}\n')
print('wrote', a.out)
