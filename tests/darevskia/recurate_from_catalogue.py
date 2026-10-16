#!/usr/bin/env python3
"""Transfer the published curated cluster -> family assignment to a new Step-8 clustering.

The curated table names clusters by their vsearch number, which changes when the motif set changes
(e.g. after the Step-4 coverage fix). This script re-derives the table from the published catalogue:
every cluster member is a TRF motif locus (header sp|scf:start-end|...); if it lies inside a catalogue
locus, that locus's family votes for the cluster. Rules, in order: a cluster holding a reference record takes
that family; a cluster whose centroid is >= 95 % identical to a family consensus over >= 100 bits
(cross/cluster_report.tsv) takes that family; otherwise a cluster takes the majority family of its
members' catalogue loci when at least 60 % of its members lie in catalogue loci and 60 % of those
agree. Everything else stays unassigned (as the published table left it).

Usage: recurate_from_catalogue.py <run_dir> <published MASTER.merged.bed (viewer format)>
                                  <published MASTER.final_assignment.clean.tsv> <out cluster_family.tsv>
"""
import collections
import os
import sys

run, bed, assign, out = sys.argv[1:5]
fam_of = {}
for line in open(assign):
    p = line.rstrip('\n').split('\t')
    if len(p) >= 2:
        fam_of[p[0]] = p[1]
loci = collections.defaultdict(list)   # (sp, scf) -> [(s, e, fam)]
for line in open(bed):
    p = line.rstrip('\n').split('\t')
    if len(p) < 4:
        continue
    sp, _, lid = p[3].split('|')[:3]
    if lid.startswith('C'):
        continue
    loci[(sp, p[0])].append((int(p[1]), int(p[2]), fam_of.get(lid, '')))

report = {}
for line in open(os.path.join(run, 'cross', 'cluster_report.tsv')):
    p = line.rstrip('\n').split('\t')
    if p[0] == 'cluster':
        continue
    report[p[0]] = {'size': int(p[1]), 'best': p[5], 'bits': float(p[6]), 'pid': p[7], 'is_ref': p[4] == 'yes', 'centroid': p[3]}

cc = os.path.join(run, 'cross', 'cross_clusters')
rows = []
for fn in sorted(os.listdir(cc), key=lambda x: int(x.split('_')[1])):
    clu = fn.split('_')[1]
    votes = collections.Counter()
    n = 0
    for line in open(os.path.join(cc, fn)):
        if not line.startswith('>'):
            continue
        n += 1
        h = line[1:].strip()
        parts = h.split('|')
        if len(parts) < 2 or ':' not in parts[1]:
            continue
        sp, reg = parts[0], parts[1]
        scf, coords = reg.rsplit(':', 1)
        s, e = (int(x) for x in coords.split('-'))
        for ls, le, fam in loci.get((sp, scf), ()):
            if s - 1 < le and e > ls:
                votes[fam] += 1
                break
    r = report.get(clu, {})
    fam, nv = (votes.most_common(1)[0] if votes else ('', 0))
    total = sum(votes.values())
    decision = ''
    try:
        pid = float(r.get('pid') or 0)
    except ValueError:
        pid = 0.0
    if r.get('is_ref'):
        decision = r['best']          # the reference consensus sits in this cluster
    elif pid >= 95.0 and (r.get('bits') or 0) >= 100.0:
        decision = r['best']          # centroid practically identical to a family consensus (over a real alignment)
    elif fam and total >= 0.6 * n and nv >= 0.6 * total:
        decision = fam                # most members lie in catalogue loci of one family
    rows.append((clu, n, decision, dict(votes), r.get('best'), r.get('bits'), r.get('pid')))

with open(out, 'w') as f:
    f.write('# cluster -> family, re-derived from the published catalogue (majority family of the catalogue loci the\n')
    f.write('# cluster members lie in, confirmed by the centroid\'s best consensus hit); see tests/darevskia/recurate_from_catalogue.py\n')
    f.write('# cluster\tfamily\tnote\n')
    for clu, n, decision, votes, best, bits, pid in rows:
        if decision:
            f.write('%s\t%s\tsize=%d votes=%s centroid_best=%s %.1f bits %s%%\n' % (clu, decision, n, votes, best, bits or 0, pid))
assigned = [r for r in rows if r[2]]
print('clusters: %d, assigned: %d' % (len(rows), len(assigned)))
for fam in sorted({r[2] for r in assigned}):
    cl = [r for r in assigned if r[2] == fam]
    print('  %s: %d clusters, %d members: %s' % (fam, len(cl), sum(r[1] for r in cl), ' '.join('%s(%d)' % (r[0], r[1]) for r in cl)))
print('unassigned clusters with catalogue votes (member loci inside catalogue loci but no clear call):')
for clu, n, decision, votes, best, bits, pid in rows:
    if not decision and votes:
        print('  %s size=%d votes=%s centroid_best=%s %s%%' % (clu, n, votes, best, pid))
