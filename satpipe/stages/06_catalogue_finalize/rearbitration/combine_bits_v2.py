# -*- coding: utf-8 -*-
"""Combine bits.tsv correctly: LOCUS-level reuse (not position-level, which failed because
scan-grid phase shifts with any boundary change), for the 439 byte-identical loci; freshly
computed bits for the 100 changed loci. Together these cover all 539 corrected loci with
either verified-reused or freshly-computed values -- no gaps, no fabrication.
"""
import collections

V = '/staging/tmp/clsat_v6'
D = ('/staging/Sergei/Lacertidae/Darevskia/Satellites/manuscript/'
     'CLsat_full_clusters_6genome/out')

identical = []
for line in open(V + '/pass3/identical_loci.bed'):
    p = line.rstrip('\n').split('\t')
    identical.append((p[0], int(p[1]), int(p[2])))
# index identical loci by scaffold for fast lookup
by_scaf = collections.defaultdict(list)
for c, s, e in identical:
    by_scaf[c].append((s, e))

for fam in ('CLsat1', 'CLsat2', 'CLsat3', 'CLsat4'):
    orig_file = D + '/MASTER_%s.bits.tsv' % fam
    changed_file = V + '/pass3/changed_%s.bits.tsv' % fam
    out_file = V + '/pass3/scan_%s.bits.tsv' % fam

    reused = 0
    out = open(out_file, 'w')
    header_written = False
    with open(orig_file) as f:
        header = f.readline()
        out.write(header)
        header_written = True
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 5:
                continue
            chrom, pos = p[0], int(p[1])
            for s, e in by_scaf.get(chrom, ()):
                if s < pos <= e:
                    out.write(line)
                    reused += 1
                    break

    added = 0
    with open(changed_file) as f:
        f.readline()
        for line in f:
            out.write(line)
            added += 1
    out.close()
    print('%s: reused_from_original=%d  fresh_from_changed=%d  total=%d'
          % (fam, reused, added, reused + added))
