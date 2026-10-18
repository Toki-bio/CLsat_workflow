"""Score a toy run against its truth (tests/toy): every planted array (truth.tsv) found by a catalogue locus
(MASTER.loci.segmented.bed + final assignment), with the right family and covering >= 90% of it.
Usage: toy_score.py <work_dir>"""
import sys, os, collections
w = sys.argv[1]
truth = [l.rstrip('\n').split('\t') for l in open(os.path.join(w, 'toy', 'truth.tsv'))][1:]
cat = os.path.join(w, 'run', 'catalogue')
seg = [l.rstrip('\n').split('\t') for l in open(os.path.join(cat, 'MASTER.loci.segmented.bed')) if l.strip()]
fam = {}
for l in open(os.path.join(cat, 'MASTER.final_assignment.clean.tsv')):
    p = l.rstrip('\n').split('\t')
    if len(p) >= 2: fam[p[0]] = p[1]
res = collections.Counter(); miss = []
for sp, chrom, s, e, kind, f, copies in truth:
    s, e = int(s), int(e)
    best = None
    for r in seg:
        if r[0] != chrom: continue
        a, b = int(r[1]), int(r[2]); ov = max(0, min(b, e) - max(a, s))
        if ov and (best is None or ov > best[0]): best = (ov, fam.get(r[3].split('.')[0], fam.get(r[3], '?')), r[3])
    if best is None: res['missed'] += 1; miss.append((kind, f, chrom, s, e, 'no locus')); continue
    cov = best[0] / (e - s)
    ok = cov >= 0.9 and best[1] == f
    res['found, right family, >=90% covered' if ok else 'found but wrong family or <90%'] += 1
    if not ok: miss.append((kind, f, chrom, s, e, 'cov %.2f, family %s (%s)' % (cov, best[1], best[2])))
print(os.path.basename(w.rstrip('/')), 'planted', len(truth), dict(res))
for m in miss: print('  ', m)
