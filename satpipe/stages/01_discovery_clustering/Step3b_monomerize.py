#!/usr/bin/env python3
"""Step 3b - monomerize: replace a TRF motif that is k tandem copies of a shorter unit by one unit.

TRF sometimes reports a whole array at a multiple of the real unit (e.g. period 292 for a 146-letter CLsat1
array). Step 3 accepts it whenever the multiple still lies inside the period range, and Step 4 does not catch it
(TRF finds nothing inside a sequence that is exactly two units long). Such a motif then enters clustering as a
dimer.

Test, for k = 2 .. len // PERIOD_MIN (so a unit is never shorter than the accepted minimum period): align the
motif with itself rotated by round(len / k) letters, globally, with indels allowed inside a band of BAND letters
(unit edit distance). If the identity is >= MIN_ID, the motif is k units. The k with the highest identity wins;
the motif is replaced by its first round(len / k) letters. No consensus, family or external tool is used.

len= becomes the unit length, copies= is multiplied by k, and |monomerized=k:<old period>:<identity> is appended.
The TRF interval of such a motif is not trusted as an array outline (a multi-unit TRF row can bridge separate arrays
and the DNA between them): if it is longer than 2 x MONOMER_SEED_UNITS units, the motif is written twice, as a start
seed and an end seed of MONOMER_SEED_UNITS units each (|seed=start, |seed=end); seed expansion then finds the real
array ends. Every change is written to <out>.monomerize.tsv.

Usage: Step3b_monomerize.py <motifs.fa> [out.fa]   (default: rewrite in place, original kept as .before_monomerize)
Environment (defaults):  PERIOD_MIN [80]  MONOMER_MIN_ID [0.70]  MONOMER_BAND [12]  MONOMER_SEED_UNITS [10]
"""
import os
import re
import sys


def banded_identity(a, b, band):
    """1 - (unit edit distance of a global alignment restricted to |i - j| <= band) / max(len)."""
    n, m = len(a), len(b)
    if abs(n - m) > band:
        return 0.0
    inf = n + m + 1
    prev = [inf] * (m + 1)
    for j in range(0, min(m, band) + 1):
        prev[j] = j
    for i in range(1, n + 1):
        cur = [inf] * (m + 1)
        lo, hi = max(0, i - band), min(m, i + band)
        if lo == 0:
            cur[0] = i
        for j in range(max(1, lo), hi + 1):
            best = prev[j - 1] + (a[i - 1] != b[j - 1])
            if prev[j] + 1 < best:
                best = prev[j] + 1
            if cur[j - 1] + 1 < best:
                best = cur[j - 1] + 1
            cur[j] = best
        prev = cur
    return 1.0 - prev[m] / max(n, m)


def unit_count(seq, period_min, min_id, band):
    """(k, identity): k tandem units if the motif matches itself rotated by one unit; (1, 0.0) otherwise."""
    s = seq.upper()
    best = (1, 0.0)
    for k in range(2, len(s) // period_min + 1):
        sh = round(len(s) / k)
        idv = banded_identity(s, (s + s)[sh:sh + len(s)], band)
        if idv >= min_id and idv > best[1]:
            best = (k, idv)
    return best


def monomerize_header(h, k, unit_len, old_len, idv):
    h = re.sub(r'\|len=[^|]*', '|len=%d' % unit_len, h, count=1)
    def _copies(mo):
        try:
            return '|copies=%s' % ('%g' % (float(mo.group(1)) * k))
        except ValueError:
            return mo.group(0)
    h = re.sub(r'\|copies=([^|]*)', _copies, h, count=1)
    return h + '|monomerized=%d:%d:%.2f' % (k, old_len, idv)


COORD = re.compile(r'([A-Za-z0-9_.-]+):([0-9]+)-([0-9]+)')


def seed_headers(h, unit_len, seed_units):
    """A motif TRF reported at a multiple of its unit can span several arrays and the non-satellite DNA between them
    (seen in Darevskia: one period-292 row over 387 kb bridged a 4.75 kb gap between two CLsat1 arrays). Its TRF
    interval is therefore not used as a seed outline: it becomes a start seed and an end seed of seed_units units
    each, and seed expansion finds the real array ends from there. Short intervals keep one seed."""
    m = COORD.search(h)
    if not m:
        return [h]
    s, e = int(m.group(2)), int(m.group(3))
    w = seed_units * unit_len
    if e - s + 1 <= 2 * w:
        return [h]
    out = []
    for tag, a, b in (('start', s, s + w - 1), ('end', e - w + 1, e)):
        x = h[:m.start(2)] + '%d-%d' % (a, b) + h[m.end(3):]
        x = re.sub(r'\|locuslen=[^|]*', '|locuslen=%d' % (b - a + 1), x, count=1)
        out.append(x + '|seed=%s' % tag)
    return out


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    src = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else src
    period_min = int(os.environ.get('PERIOD_MIN', 80))
    min_id = float(os.environ.get('MONOMER_MIN_ID', 0.70))
    band = int(os.environ.get('MONOMER_BAND', 12))
    seed_units = int(os.environ.get('MONOMER_SEED_UNITS', 10))
    recs, h, buf = [], None, []
    for line in open(src):
        line = line.rstrip('\n')
        if line.startswith('>'):
            if h is not None:
                recs.append((h, ''.join(buf)))
            h, buf = line[1:], []
        elif line:
            buf.append(line.strip())
    if h is not None:
        recs.append((h, ''.join(buf)))
    if out == src:
        os.replace(src, src + '.before_monomerize')
    changed = []
    with open(out, 'w') as f:
        for h, s in recs:
            k, idv = unit_count(s, period_min, min_id, band)
            if k > 1:
                unit = s[:round(len(s) / k)]
                changed.append((h, len(s), k, len(unit), idv))
                for hh in seed_headers(monomerize_header(h, k, len(unit), len(s), idv), len(unit), seed_units):
                    f.write('>%s\n%s\n' % (hh, unit))
                continue
            f.write('>%s\n%s\n' % (h, s))
    with open(out + '.monomerize.tsv', 'w') as f:
        f.write('header\tmotif_len\tunits\tunit_len\tself_identity\n')
        for h, L, k, ul, idv in changed:
            f.write('%s\t%d\t%d\t%d\t%.3f\n' % (h, L, k, ul, idv))
    print('Step 3b: %s: %d motifs, %d monomerized (k=2: %d, k>=3: %d)'
          % (os.path.basename(out), len(recs), len(changed),
             sum(1 for c in changed if c[2] == 2), sum(1 for c in changed if c[2] >= 3)))


if __name__ == '__main__':
    main()
