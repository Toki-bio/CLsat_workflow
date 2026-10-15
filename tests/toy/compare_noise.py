#!/usr/bin/env python3
"""Noise-aware comparison of bits-derived tables for the toy equivalence test.

ssearch36 estimates bit scores from shuffled-library statistics, so identical alignments get
different bits from run to run (up to ~15% on the toy). Quantities derived from bits inherit
that noise: delta = best - second, and any flag or boolean that thresholds delta. This
comparer requires every structural column to be identical, allows bits within a relative
tolerance, allows deltas within tol x (|a| + |b|) of the two bits they come from, and accepts a
flag/boolean difference only when the two deltas straddle the threshold within that same
noise band.

Usage: compare_noise.py assignment  <old.tsv> <new.tsv> [--tol 0.25] [--delta 20] [--delta-ambig 15]
       compare_noise.py subfamily   <old.tsv> <new.tsv> [--tol 0.25] [--confident 10]
       compare_noise.py subfamily-json <old.json> <new.json> [--tol 0.25] [--confident 10]
Exit 0 = equivalent within noise; 1 = a real difference (printed).
"""
import argparse
import json
import sys

ap = argparse.ArgumentParser()
ap.add_argument('mode', choices=['assignment', 'subfamily', 'subfamily-json'])
ap.add_argument('old')
ap.add_argument('new')
ap.add_argument('--tol', type=float, default=0.25)
ap.add_argument('--delta', type=float, default=20)
ap.add_argument('--delta-ambig', type=float, default=15)
ap.add_argument('--confident', type=float, default=10)
a = ap.parse_args()
TOL = a.tol
problems = []
noise_flips = 0


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def near(u, v, scale=None):
    if u is None or v is None:
        return u == v
    scale = max(abs(u), abs(v), 1.0) if scale is None else max(scale, 1.0)
    return abs(u - v) <= TOL * scale


def straddles(d1, d2, thr, band):
    """the two deltas lie on different sides of thr, both within band of it"""
    if d1 is None or d2 is None:
        return False
    return (d1 - thr) * (d2 - thr) <= 0 and abs(d1 - thr) <= band and abs(d2 - thr) <= band


def rows(path):
    return [l.rstrip('\n').split('\t') for l in open(path) if l.strip()]


if a.mode == 'assignment':
    # locus_id family assigned_bits best_other best_other_bits delta flag provenance
    A, B = rows(a.old), rows(a.new)
    if len(A) != len(B):
        problems.append('row count %d vs %d' % (len(A), len(B)))
    for i, (x, y) in enumerate(zip(A, B), 1):
        for c in (0, 1, 3, 7):
            if x[c] != y[c]:
                problems.append('line %d col %d: %r vs %r' % (i, c + 1, x[c], y[c]))
        b1, b2 = num(x[2]), num(y[2]); o1, o2 = num(x[4]), num(y[4]); d1, d2 = num(x[5]), num(y[5])
        if not near(b1, b2) or not near(o1, o2):
            problems.append('line %d bits: %s/%s vs %s/%s' % (i, x[2], x[4], y[2], y[4]))
        band = TOL * ((abs(b1) if b1 else 0) + (abs(o1) if o1 else 0))
        if not near(d1, d2, scale=(abs(b1) if b1 else 0) + (abs(o1) if o1 else 0)):
            problems.append('line %d delta: %s vs %s' % (i, x[5], y[5]))
        if x[6] != y[6]:
            if straddles(d1, d2, -a.delta, band) or straddles(d1, d2, a.delta_ambig, band) or straddles(d1, d2, 0, band):
                noise_flips += 1
                print('    noise flip line %d: %s (delta %s) vs %s (delta %s)' % (i, x[6], x[5], y[6], y[5]))
            else:
                problems.append('line %d flag: %r vs %r (delta %s vs %s)' % (i, x[6], y[6], x[5], y[5]))

elif a.mode == 'subfamily':
    # locus major subfamily pct_id aln_len bits second second_bits delta
    A, B = rows(a.old), rows(a.new)
    if len(A) != len(B):
        problems.append('row count %d vs %d' % (len(A), len(B)))
    for i, (x, y) in enumerate(zip(A, B), 1):
        if i == 1:
            continue
        for c in (0, 1, 2, 3, 4, 6):
            if x[c] != y[c]:
                problems.append('line %d col %d: %r vs %r' % (i, c + 1, x[c], y[c]))
        b1, b2 = num(x[5]), num(y[5]); s1, s2 = num(x[7]), num(y[7]); d1, d2 = num(x[8]), num(y[8])
        if not near(b1, b2) or not near(s1, s2):
            problems.append('line %d bits: %s/%s vs %s/%s' % (i, x[5], x[7], y[5], y[7]))
        if not near(d1, d2, scale=(abs(b1) if b1 else 0) + (abs(s1) if s1 else 0)):
            problems.append('line %d delta: %s vs %s' % (i, x[8], y[8]))

elif a.mode == 'subfamily-json':
    A, B = json.load(open(a.old)), json.load(open(a.new))
    if A.keys() != B.keys():
        problems.append('locus sets differ')
    for lid in A:
        x, y = A[lid], B.get(lid, {})
        for k in ('major', 'subfamily', 'second', 'pid'):
            if x.get(k) != y.get(k):
                problems.append('%s %s: %r vs %r' % (lid, k, x.get(k), y.get(k)))
        b1, b2 = x.get('bits'), y.get('bits'); d1, d2 = x.get('delta'), y.get('delta')
        if not near(b1, b2):
            problems.append('%s bits: %s vs %s' % (lid, b1, b2))
        scale = 2 * (abs(b1) if b1 else 0)
        if not near(d1, d2, scale=scale):
            problems.append('%s delta: %s vs %s' % (lid, d1, d2))
        if x.get('confident') != y.get('confident'):
            if straddles(d1, d2, a.confident, TOL * scale):
                noise_flips += 1
                print('    noise flip %s: confident %s (delta %s) vs %s (delta %s)' % (lid, x.get('confident'), d1, y.get('confident'), d2))
            else:
                problems.append('%s confident: %s vs %s (delta %s vs %s)' % (lid, x.get('confident'), y.get('confident'), d1, d2))

for p in problems[:10]:
    print('    ' + p)
if noise_flips:
    print('    (%d threshold flip(s) within ssearch36 bits noise, tolerated)' % noise_flips)
sys.exit(1 if problems else 0)
