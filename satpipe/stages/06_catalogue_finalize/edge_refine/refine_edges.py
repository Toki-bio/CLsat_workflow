# -*- coding: utf-8 -*-
"""Refine every locus boundary to NUCLEOTIDE resolution.

Needed for every family: with STEP_MON=3 the expansion walk discards a whole failing window
(3 x MONO), so edges are quantised and can sit up to ~3 monomers inside the true array end.
The expansion finds the right WINDOW; this finds the right BASE.

The rule (validated against a manual call on V1, which the user checked by hand):
  edge = the LAST BASE THAT MATCHES the locus's own end-local consensus, scanning outward,
         stopping at the first run of MAXMISS consecutive mismatches,
         BOUNDED ABOVE by the ssearch36 window coverage.

Why each piece is there, all of them learned the hard way:
  * end-local consensus, not the family consensus: terminal monomers can belong to a different
    subfamily (V2 is a CLsat1 locus whose 3' end is CLsat3).
  * strict last-match, not match density: density carries mismatches and over-ran V1 by 2 bp.
    The user's manual call is the strict one.
  * run-break (MAXMISS): past V1's true edge the genome reads TTTTC and the consensus TTTTC
    again by chance, so an unbounded "last match anywhere" walks off into the flank.
  * bounded by SW coverage: strict alone runs away in low-complexity flank -- at V517's 5' end
    the telomeric (CCCTAA)n never produces 4 consecutive mismatches and it ate 25 bp of
    telomere. SW coverage over-extends but never under-extends, so it is a safe upper bound.
  * the window test thresholds are on ALIGNED LENGTH and BITS, never on identity: flank windows
    reach 87.5% identity because SW aligns 8 bases perfectly. Identity carries no signal here.

Emits a refined BED plus a per-end report, and flags contig edges (which are assembly
truncations, not biological array ends, and must never be treated as boundaries).

Usage: refine_edges.py <loci.bed> <out_prefix> <genomes.tsv> <families.tsv>
  loci.bed      chrom  start0  end  species  family   (5 columns; family = the locus's call)
  genomes.tsv   species_code  genome.fa (indexed)
  families.tsv  family  mono  consensus.fa
Environment: EDGE_TOL [50]  MAXMISS [4]  GAP_TOL [150]  LEAD_TOL [30]  MIN_ALEN_FRAC [0.85]
             MIN_BITS_FRAC [0.55]  MIN_BLOCK [60]  PAD [1500]
Tools: ssearch36, samtools on PATH.
"""
import collections
import os
import shutil
import subprocess
import sys
import tempfile

if len(sys.argv) != 5:
    sys.exit('Usage: refine_edges.py <loci.bed> <out_prefix> <genomes.tsv> <families.tsv>')
BED_IN, OUT, GENOMES_TSV, FAMILIES_TSV = sys.argv[1:5]

SS = shutil.which('ssearch36') or sys.exit('ssearch36 not on PATH')
SAM = shutil.which('samtools') or sys.exit('samtools not on PATH')

MAXMISS = int(os.environ.get('MAXMISS', 4))
GAP_TOL = int(os.environ.get('GAP_TOL', 150))
LEAD_TOL = int(os.environ.get('LEAD_TOL', 30))
MIN_ALEN_FRAC = float(os.environ.get('MIN_ALEN_FRAC', 0.85))
MIN_BITS_FRAC = float(os.environ.get('MIN_BITS_FRAC', 0.55))
MIN_BLOCK = int(os.environ.get('MIN_BLOCK', 60))
EDGE_TOL = int(os.environ.get('EDGE_TOL', 50))
PAD = int(os.environ.get('PAD', 1500))          # how far outside the called edge we look


def read_tsv(path, min_cols):
    rows = []
    for line in open(path):
        if not line.strip() or line.startswith('#'):
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) < min_cols:
            sys.exit('%s: need %d columns, got %r' % (path, min_cols, line.strip()))
        rows.append(p)
    return rows


MONO = {}
FAMCONS = {}
for p in read_tsv(FAMILIES_TSV, 3):
    MONO[p[0]] = int(p[1])
    FAMCONS[p[0]] = ''.join(l.strip() for l in open(p[2]) if not l.startswith('>')).upper()
GEN = {p[0]: p[1] for p in read_tsv(GENOMES_TSV, 2)}
SL = {}
for sp, g in GEN.items():
    for line in open(g + '.fai'):
        p = line.split('\t')
        SL[(sp, p[0])] = int(p[1])

TR = str.maketrans('ACGTN', 'TGCAN')


def rc(s):
    return s.translate(TR)[::-1]


def ident(x, y):
    n = min(len(x), len(y))
    return sum(1 for i in range(n) if x[i] == y[i]) / n if n else 0.0


def faidx(sp, scf, s, e):
    s = max(1, s)
    e = min(SL[(sp, scf)], e)
    if e < s:
        return ''
    r = subprocess.run([SAM, 'faidx', GEN[sp], '%s:%d-%d' % (scf, s, e)],
                       capture_output=True, text=True)
    if r.returncode:
        return ''
    return ''.join(l.strip() for l in r.stdout.splitlines() if not l.startswith('>')).upper()


def sat_cov(seq, lib, wlen, step, cov=None):
    if cov is None:
        cov = bytearray(len(seq))
    if len(seq) < wlen:
        return cov
    with tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False) as f:
        for i in range(0, len(seq) - wlen + 1, step):
            f.write('>w%d\n%s\n' % (i, seq[i:i + wlen]))
        q = f.name
    o = subprocess.run([SS, '-n', '-z', '3', '-m', '8C', '-E', '100000', '-b', '1', '-d', '1', q, lib],
                       capture_output=True, text=True).stdout
    os.unlink(q)
    best = {}
    for line in o.splitlines():
        if line.startswith('#'):
            continue
        p = line.split('\t')
        if len(p) < 12:
            continue
        try:
            off = int(p[0][1:]); alen = int(p[3])
            qs, qe = int(p[6]), int(p[7]); bits = float(p[11])
        except ValueError:
            continue
        if off not in best or bits > best[off][3]:
            best[off] = (alen, min(qs, qe), max(qs, qe), bits)
    for off, (alen, qs, qe, bits) in best.items():
        if alen >= MIN_ALEN_FRAC * wlen and bits >= MIN_BITS_FRAC * wlen:
            for i in range(off + qs - 1, min(off + qe, len(seq))):
                cov[i] = 1
    return cov


def end_consensus(arr_out, M):
    """arr_out runs ARRAY->OUTWARD; build a consensus from the terminal monomers,
    each phased individually, and say which family they look like."""
    mons = []
    for k in range(30):
        x = k * M
        if x + 3 * M > len(arr_out):
            break
        seg = arr_out[x:x + 3 * M]
        b = (0, -1, None)
        for f2, c in FAMCONS.items():
            if MONO[f2] != M:
                continue
            d = c + c
            for r in range(M):
                sc = ident(seg[:M], d[r:r + M])
                if sc > b[1]:
                    b = (r, sc, f2)
        if b[1] < 0.60:
            continue
        off = (M - b[0]) % M
        m = seg[off:off + M]
        if len(m) == M:
            mons.append((m, b[2]))
    if not mons:
        return None, None, 0.0
    fam = collections.Counter(f for _, f in mons).most_common(1)[0][0]
    ms = [m for m, _ in mons]
    con = ''.join((collections.Counter(m[i] for m in ms if m[i] in 'ACGT').most_common(1)
                   or [('N',)])[0][0] for i in range(M))
    d = FAMCONS[fam] + FAMCONS[fam]
    kb = max(range(M), key=lambda k: ident(con, d[k:k + M]))
    con = con[(M - kb) % M:] + con[:(M - kb) % M] if kb else con
    mid = sorted(ident(m, con) for m in ms)[len(ms) // 2]
    return con, fam, mid


def refine_one_end(sp, scf, edge, direction, fam_hint):
    """edge = 1-based coord of the last base INSIDE the array on this side.
    direction = +1 (array is to the LEFT, we walk right) or -1 (array is to the RIGHT).
    returns (new_edge, shift, reason, fam, conf)"""
    M = MONO[fam_hint]
    inside = 32 * M
    L = SL[(sp, scf)]
    if direction > 0:
        g1, g2 = edge - inside + 1, edge + PAD
        seq = faidx(sp, scf, g1, g2)
        arr_out = seq                       # array first, then outward
        e0 = edge - g1 + 1                  # index one past the last inside base
        gedge = edge
    else:
        g1, g2 = edge - PAD, edge + inside - 1
        seq = faidx(sp, scf, g1, g2)
        arr_out = rc(seq)                   # now array-first, outward
        e0 = (g2 - edge + 1)
        gedge = edge
    if len(arr_out) < inside // 2:
        return edge, 0, 'too_short', fam_hint, 0.0
    if gedge <= EDGE_TOL or gedge >= L - EDGE_TOL:
        return edge, 0, 'CONTIG_EDGE', fam_hint, 0.0

    con, fam, conf = end_consensus(arr_out[max(0, e0 - 30 * M):e0], M)
    if not con:
        return edge, 0, 'no_consensus', fam_hint, 0.0
    # per-call temp library (the former fixed /tmp path collided under parallel runs)
    with tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False) as lf:
        lf.write('>c\n%s%s\n' % (con, con))
        lib = lf.name
    try:
        cov = sat_cov(arr_out, lib, 60, 20)
        d = 0
        while e0 + d < len(arr_out) and cov[e0 + d]:
            d += 1
        lo = max(0, e0 + d - 300)
        hi = min(len(arr_out), e0 + d + 300)
        sub = arr_out[lo:hi]
        fc = sat_cov(sub, lib, 40, 1)
    finally:
        os.unlink(lib)
    for i in range(len(sub)):
        if fc[i]:
            cov[lo + i] = 1

    # walk outward with the coverage (upper bound)
    if not any(cov[e0:e0 + LEAD_TOL]):
        resid_cov = 0
    else:
        last = 0
        j = 0
        n = len(arr_out)
        while e0 + j < n:
            if cov[e0 + j]:
                last = j + 1; j += 1; continue
            g = j
            while e0 + g < n and not cov[e0 + g]:
                g += 1
            blk = 0
            while e0 + g + blk < n and cov[e0 + g + blk]:
                blk += 1
            if e0 + g < n and (g - j) <= GAP_TOL and blk >= MIN_BLOCK:
                j = g; continue
            break
        resid_cov = last

    # strict last-match, bounded by resid_cov
    D = con * ((PAD + 2 * M) // M + 4)
    S = arr_out[max(0, e0 - 60):min(len(arr_out), e0 + resid_cov + 60)]
    base = max(0, e0 - 60)
    anchor = e0 - base
    if anchor < 10 or len(S) <= anchor:
        strict = 0
    else:
        kb = max(range(M), key=lambda k: ident(S[:anchor], D[k:k + anchor]))
        R = D[kb:kb + len(S)]
        last_match = anchor - 1
        miss = 0
        for i in range(anchor, len(S)):
            if S[i] == R[i]:
                last_match = i; miss = 0
            else:
                miss += 1
                if miss >= MAXMISS:
                    break
        strict = max(0, (base + last_match + 1) - e0)
    resid = min(strict, resid_cov)
    new_edge = edge + direction * resid
    return new_edge, resid, ('refined' if resid else 'unchanged'), fam, conf


loci = []
for line in open(BED_IN):
    p = line.rstrip('\n').split('\t')
    if not p or not p[0] or p[0].startswith('#'):
        continue
    if len(p) < 5:
        sys.exit('%s: need 5 columns (chrom start end species family), got %r' % (BED_IN, line.strip()))
    if p[4] not in MONO:
        sys.exit('%s: family %r is not in %s' % (BED_IN, p[4], FAMILIES_TSV))
    loci.append((p[0], int(p[1]), int(p[2]), p[3], p[4]))
print('loci to refine: %d' % len(loci), flush=True)

fb = open(OUT + '.bed', 'w')
fr = open(OUT + '.report.tsv', 'w')
fr.write('locus\tscaffold\tsp\told_start\told_end\tnew_start\tnew_end\t'
         'shift_L\tshift_R\treason_L\treason_R\tfam_L\tfam_R\tconf_L\tconf_R\n')
stat = collections.Counter()
shifts = []
for i, (scf, s, e, sp, fam) in enumerate(loci):
    if (sp, scf) not in SL:
        continue
    nL, shL, rL, fL, cL = refine_one_end(sp, scf, s + 1, -1, fam)
    nR, shR, rR, fR, cR = refine_one_end(sp, scf, e, +1, fam)
    ns, ne = min(nL - 1, s), max(nR, e)
    stat[rL] += 1; stat[rR] += 1
    if shL:
        shifts.append(shL)
    if shR:
        shifts.append(shR)
    fb.write('%s\t%d\t%d\tR%d\t%s\t%s\n' % (scf, ns, ne, i + 1, sp, fam))
    fr.write('R%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%.3f\t%.3f\n'
             % (i + 1, scf, sp, s, e, ns, ne, shL, shR, rL, rR, fL, fR, cL, cR))
    if (i + 1) % 50 == 0:
        print('  %d/%d' % (i + 1, len(loci)), flush=True)
fb.close(); fr.close()

n = sum(stat.values())
print('\n=== boundary outcomes (%d ends) ===' % n)
for k, v in stat.most_common():
    print('  %-14s %5d (%.1f%%)' % (k, v, 100.0 * v / n))
if shifts:
    shifts.sort()
    print('\nedges that MOVED: %d | median %d bp, mean %d, max %d'
          % (len(shifts), shifts[len(shifts) // 2], sum(shifts) // len(shifts), shifts[-1]))
print('\nwrote %s.bed and %s.report.tsv' % (OUT, OUT))
