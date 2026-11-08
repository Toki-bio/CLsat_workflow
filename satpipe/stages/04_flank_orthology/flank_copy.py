# -*- coding: utf-8 -*-
"""Single-copy anchors for flank orthology (stage 19). Author's rules, 2026-11-08:

1. Tandem mask: letters of a flank that lie inside any TRF array of its own genome (the run's genome-wide TRF
   calls, every period) are masked with N. A flank with fewer than MIN_ALN unmasked letters is class 'tandem'
   and is no anchor. So a satellite inside a flank can neither make the flank look repeated nor link two loci.
2. Single copy: the masked flank must occur once in its own whole genome. A second place with >= COPY_MIN_PID %
   identity over >= COPY_MIN_ALN letters makes the flank 'repeated' (no anchor). Two searches: the flanks of the
   other loci of the same genome (ssearch36, exact, the stage's own all-vs-all) and the whole genome (minimap2).
3. Clear best hit: in a target genome, the best anchor hit must beat the second-best hit at a DIFFERENT locus:
   second < CLEAR_MARGIN x best (bits). Ties never decide.
4. Reciprocal: a locus pair is kept only when a flank of each locus has the other locus as its clear best hit.

Pure functions only (no tools), so they can be tested on their own.
"""
import bisect
import collections


def parse_trf_dat(path, wanted):
    """TRF .dat (1-based inclusive) -> {scaffold: [(start, end), ...]} for scaffolds in `wanted`, merged and sorted.
    Streams the file; scaffolds not in `wanted` are skipped without storing."""
    out = collections.defaultdict(list)
    scf = None
    with open(path) as fh:
        for line in fh:
            if line.startswith('Sequence:'):
                scf = line.split(None, 1)[1].split()[0]
                continue
            if scf not in wanted or not line[:1].isdigit():
                continue
            p = line.split(None, 2)
            try:
                s, e = int(p[0]), int(p[1])
            except (ValueError, IndexError):
                continue
            out[scf].append((s, e))
    return {k: merge(v) for k, v in out.items()}


def merge(iv):
    iv = sorted(iv)
    m = []
    for s, e in iv:
        if m and s <= m[-1][1] + 1:
            if e > m[-1][1]:
                m[-1] = (m[-1][0], e)
        else:
            m.append((s, e))
    return m


def mask(seq, fs, iv):
    """Mask (N) the letters of `seq` (genomic, plus strand, first letter at 1-based position fs) that fall inside the
    merged intervals `iv` (1-based inclusive). Returns (masked_seq, masked_letters)."""
    if not iv:
        return seq, 0
    fe = fs + len(seq) - 1
    s = list(seq)
    n = 0
    k = bisect.bisect_right([a for a, _ in iv], fe)
    for a, b in iv[:k]:
        if b < fs:
            continue
        lo, hi = max(a, fs), min(b, fe)
        for x in range(lo - fs, hi - fs + 1):
            if s[x] != 'N':
                s[x] = 'N'
                n += 1
    return ''.join(s), n


def unmasked(seq):
    return sum(1 for c in seq if c not in 'Nn')


def paf_other_copies(paf_lines, self_place, min_pid, min_aln):
    """Count minimap2 PAF hits per query that are NOT the query's own place.
    self_place: {query: (scaffold, start1, end1)} (1-based inclusive genomic place of the flank).
    A hit counts when aligned block >= min_aln and identity (matches / block, in %) >= min_pid.
    A hit overlapping the query's own place on its own scaffold is the flank itself and is not counted."""
    n = collections.Counter()
    for line in paf_lines:
        p = line.rstrip('\n').split('\t')
        if len(p) < 12:
            continue
        q, t, ts, te, match, block = p[0], p[5], int(p[7]), int(p[8]), int(p[9]), int(p[10])
        if block < min_aln or 100.0 * match / block < min_pid:
            continue
        sp = self_place.get(q)
        if sp and t == sp[0] and ts < sp[2] and te > sp[1] - 1:      # PAF target coords 0-based half-open
            continue
        n[q] += 1
    return n


def clear_best(hits, margin):
    """hits: [(bits, target_locus, payload...)] of ONE query flank in ONE target genome.
    Best per target locus first (a locus can be hit by both of its flanks); returns the best entry when the
    second-best locus scores < margin x best, else None (tie or near tie: undecided)."""
    per = {}
    for h in hits:
        if h[1] not in per or h[0] > per[h[1]][0]:
            per[h[1]] = h
    if not per:
        return None
    r = sorted(per.values(), key=lambda h: (-h[0], h[1]))
    if len(r) > 1 and r[1][0] >= margin * r[0][0]:
        return None
    return r[0]


def reciprocal_pairs(directed):
    """directed: set of (locusA, locusB) meaning a flank of A has B as clear best hit.
    Returns (pairs kept: both directions present, one-way pairs)."""
    keep, oneway = set(), set()
    for a, b in directed:
        key = tuple(sorted((a, b)))
        if (b, a) in directed:
            keep.add(key)
        else:
            oneway.add((a, b))
    return keep, oneway


# ---- refinements, author 2026-11-08 (after the rules left molly with 84 anchors and Darevskia with 7) ----

def n_union(a, b):
    """Two versions of the same flank (same length): N wherever either has N."""
    if len(a) != len(b):
        raise ValueError('n_union: lengths differ (%d vs %d)' % (len(a), len(b)))
    return ''.join('N' if (x in 'Nn' or y in 'Nn') else x for x, y in zip(a, b))


def paf_other_pieces(paf_lines, self_place, min_pid, min_aln):
    """Like paf_other_copies, but returns the query pieces that have another copy:
    {query: [(start1, end1), ...]} (1-based inclusive, on the query as written), for hits with aligned block
    >= min_aln and identity >= min_pid that are not the query's own place."""
    out = collections.defaultdict(list)
    for line in paf_lines:
        p = line.rstrip('\n').split('\t')
        if len(p) < 12:
            continue
        q, qs, qe, t, ts, te, match, block = p[0], int(p[2]), int(p[3]), p[5], int(p[7]), int(p[8]), int(p[9]), int(p[10])
        if block < min_aln or 100.0 * match / block < min_pid:
            continue
        sp = self_place.get(q)
        if sp and t == sp[0] and ts < sp[2] and te > sp[1] - 1:
            continue
        out[q].append((qs + 1, qe))
    return dict(out)


def cut_flank(side, fs, fe, s, e, contig_len, loci_on_contig, self_i, guard):
    """Shorten a flank window instead of rejecting it (--cut-at-locus).
    Flank [fs, fe] is 1-based inclusive; loci are (s2, e2, j), 0-based half-open, so a locus occupies the 1-based
    letters s2+1..e2. When another locus lies on the flank's side, the letters between the two loci are split at
    their middle: the left locus' R flank ends at the middle letter, the right locus' L flank starts after it, so
    no letter serves as flank for two loci and a flank never reaches into another locus (guard letters kept free).
    Returns (fs, fe, why): why = 'locus' if another locus cut it, 'end' if the contig end did, '' otherwise.
    fe < fs means nothing is left."""
    why = ''
    if side == 'L':
        for (s2, e2, j) in loci_on_contig:
            if j == self_i or s2 >= s:
                continue
            mid = (e2 + 1 + s) // 2
            lo = max(mid + 1, e2 + 1 + guard)
            if lo > fs:
                fs, why = lo, 'locus'
        if fs < 1:
            fs, why = 1, why or 'end'
    else:
        for (s2, e2, j) in loci_on_contig:
            if j == self_i or e2 <= e:
                continue
            mid = (e + 1 + s2) // 2
            hi = min(mid, s2 - guard)
            if hi < fe:
                fe, why = hi, 'locus'
        if fe > contig_len:
            fe, why = contig_len, why or 'end'
    return fs, fe, why
