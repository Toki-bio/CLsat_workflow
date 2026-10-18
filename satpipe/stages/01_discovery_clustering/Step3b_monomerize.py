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
and the DNA between them): if it is longer than 2 x MONOMER_SEED_UNITS units, the interval is scanned for CLsat (the
manuscript's definition: 150-letter windows with >= 80 letters aligned at >= 80% to a family consensus; windows
within ~1 kb join into a run; MONOMER_GENOMES = genomes.tsv, MONOMER_CONS = consensus files joined by ':') and every
run becomes its own seed (|seed=run<i>). Without genomes or tools, or when no run is found, the motif falls back to a
start seed and an end seed of MONOMER_SEED_UNITS units each (|seed=start, |seed=end) and seed expansion has to find the
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


def scan_runs(genome, scaf, s, e, consensus_files, win=150, min_aln=80, min_id=80.0, max_gap_windows=7, min_windows=3):
    """Runs of CLsat-like sequence inside the 1-based interval [s, e] of a genome, found the way the manuscript
    defines CLsat: 150-letter windows count when >= min_aln letters align to a family consensus at >= min_id % (ssearch36
    against the consensuses doubled); windows within max_gap_windows (7 = ~1 kb) of each other form a run; runs of fewer
    than min_windows windows are dropped. Returns [(start, end)] 1-based inclusive, or None if the tools are missing."""
    import shutil
    import subprocess
    import tempfile
    if not (shutil.which('samtools') and shutil.which('ssearch36')):
        return None
    out = subprocess.run(['samtools', 'faidx', genome, '%s:%d-%d' % (scaf, s, e)], capture_output=True, text=True).stdout
    seq = ''.join(out.split('\n')[1:]).upper()
    if not seq:
        return None
    d = tempfile.mkdtemp()
    qf, lf = os.path.join(d, 'q.fa'), os.path.join(d, 'lib.fa')
    with open(qf, 'w') as f:
        for i in range(0, len(seq) - win + 1, win):
            f.write('>w%d\n%s\n' % (i // win, seq[i:i + win]))
    with open(lf, 'w') as f:
        for p in consensus_files:
            c = ''.join(l.strip() for l in open(p) if not l.startswith('>')).upper()
            f.write('>%s\n%s\n' % (os.path.basename(p).split('.')[0], c + c))
    res = subprocess.run(['ssearch36', '-n', '-m', '8', '-E', '1000', '-b', '1', '-d', '1', qf, lf],
                         capture_output=True, text=True).stdout
    hit = set()
    for line in res.splitlines():
        p = line.split('\t')
        if len(p) >= 12 and not line.startswith('#'):
            try:
                if int(p[3]) >= min_aln and float(p[2]) >= min_id:
                    hit.add(int(p[0][1:]))
            except ValueError:
                pass
    runs, cur = [], None
    for w in sorted(hit):
        if cur and w - cur[1] <= max_gap_windows + 1:
            cur[1] = w; cur[2] += 1
        else:
            if cur:
                runs.append(cur)
            cur = [w, w, 1]
    if cur:
        runs.append(cur)
    shutil.rmtree(d, ignore_errors=True)
    return [(s + a * win, s + (b + 1) * win - 1) for a, b, n in runs if n >= min_windows]


def run_headers(h, runs):
    """One seed header per run (coordinates replaced, |seed=run<i> appended)."""
    m = COORD.search(h)
    out = []
    for i, (a, b) in enumerate(runs, 1):
        x = h[:m.start(2)] + '%d-%d' % (a, b) + h[m.end(3):]
        x = re.sub(r'\|locuslen=[^|]*', '|locuslen=%d' % (b - a + 1), x, count=1)
        out.append(x + '|seed=run%d' % i)
    return out


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
    genomes, cons = {}, []
    if os.environ.get('MONOMER_GENOMES') and os.environ.get('MONOMER_CONS'):
        for line in open(os.environ['MONOMER_GENOMES']):
            p = line.rstrip('\n').split('\t')
            if len(p) >= 2 and not line.startswith('#'):
                genomes[p[0]] = p[1]
        cons = [c for c in os.environ['MONOMER_CONS'].split(':') if c]
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
                mh = monomerize_header(h, k, len(unit), len(s), idv)
                heads = None
                mm = COORD.search(h)
                if genomes and cons and mm and (int(mm.group(3)) - int(mm.group(2)) + 1) > 2 * seed_units * len(unit):
                    code = h.split('|')[0]
                    if code in genomes:
                        runs = scan_runs(genomes[code], mm.group(1), int(mm.group(2)), int(mm.group(3)), cons)
                        if runs:
                            heads = run_headers(mh, runs)
                if heads is None:
                    heads = seed_headers(mh, len(unit), seed_units)
                for hh in heads:
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
