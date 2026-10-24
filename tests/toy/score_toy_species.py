"""Genome-aware score of a toy run (tests/toy, best with make_toy.py --same-names): every planted array must be
found by a catalogue locus of the SAME genome with the right family covering >= 90 % of it, no locus may be filed
under a genome where nothing was planted at its place, and every locus must have a sub-subfamily row.
Reads regions/locus_master.tsv (Lid scf start0 end species family MONO) and regions/locus_subfamily.tsv.
Usage: score_toy_species.py <work_dir>     exit 0 = pass"""
import collections
import os
import sys

w = sys.argv[1]
truth = [l.rstrip('\n').split('\t') for l in open(os.path.join(w, 'toy', 'truth.tsv'))][1:]
reg = os.path.join(w, 'run', 'regions')
loci = [l.rstrip('\n').split('\t') for l in open(os.path.join(reg, 'locus_master.tsv')) if l.strip()]
problems = []
res = collections.Counter()
for sp, chrom, s, e, kind, fam, copies in truth:
    s, e = int(s), int(e)
    best = None
    for lid, scf, a, b, lsp, lfam, mono in loci:
        if lsp != sp or scf != chrom:
            continue
        ov = max(0, min(int(b), e) - max(int(a), s))
        if ov and (best is None or ov > best[0]):
            best = (ov, lfam, lid)
    if best is None:
        res['missed'] += 1; problems.append('missed %s %s %s:%d-%d %s' % (kind, fam, sp, s, e, chrom)); continue
    ok = best[0] / (e - s) >= 0.9 and best[1] == fam
    res['ok' if ok else 'wrong'] += 1
    if not ok:
        problems.append('%s %s %s:%d-%d -> %s cov %.2f family %s' % (kind, fam, sp, s, e, best[2], best[0] / (e - s), best[1]))
for lid, scf, a, b, lsp, lfam, mono in loci:
    a, b = int(a), int(b)
    if not any(t[0] == lsp and t[1] == scf and min(b, int(t[3])) > max(a, int(t[2])) for t in truth):
        problems.append('locus %s filed under %s %s:%d-%d where nothing was planted' % (lid, lsp, scf, a, b))
sub = {l.split('\t')[0] for l in open(os.path.join(reg, 'locus_subfamily.tsv')) if not l.startswith('locus\t')}
lost = [l[0] for l in loci if l[0] not in sub]
if lost:
    problems.append('loci without a sub-subfamily row: %d (%s)' % (len(lost), ', '.join(lost[:5])))
print('planted', len(truth), dict(res), 'loci', len(loci))
for p in problems:
    print('  ', p)
print('SPECIES TOY: PASS' if not problems else 'SPECIES TOY: FAIL')
sys.exit(1 if problems else 0)
