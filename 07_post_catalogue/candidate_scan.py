#!/usr/bin/env python3
"""Process the 7 CLsat-matching singleton loci that are absent from the catalogue,
using the same steps the pipeline applies to MASTER loci, so they can be inspected
in the viewer alongside catalogued loci.

For each candidate:
  * extract the locus plus 600 bp of flank on each side
  * four-way sliding-window ssearch36 scan against the CLsat1-4 consensuses
    (window = 2 monomers, step = 1 monomer, per-family monomer length)
  * per-locus maximum bitscore per family -> family call and competitor
  * 600 bp flank extraction and classification
  * sub-subfamily assignment against the doubled sub-subfamily consensuses

These are CANDIDATES for inspection, not catalogue members. They are given C-prefixed
identifiers and a distinct flag so they can never be confused with the 627/541 set.
"""
import os
import subprocess
import json
import collections

BASE = '/staging/Sergei/Lacertidae/Darevskia'
SAT = BASE + '/Satellites'
OUT = SAT + '/motifs_output/clean/candidates'
SS = '/staging/conda/envs/bioinfo/bin/ssearch36'
CONS = SAT + '/CLUSTERS/CLsat/consi_CLsat_original/CLsat_cons.rot'

os.makedirs(OUT, exist_ok=True)

CAND = [
    ('C1', 'unp', 'JAWDKN010000111.1', 1046086, 1047540),
    ('C2', 'val', 'JAWWNG010000434.1', 5035624, 5037236),
    ('C3', 'nai', 'JAWWNI010000232.1',  487484,  489487),
    ('C4', 'unm', 'JAWDKO010000033.1', 6419919, 6421922),
    ('C5', 'val', 'JAWWNG010000434.1', 5030863, 5032532),
    ('C6', 'val', 'JAWWNG010000953.1',  378254,  380111),
    ('C7', 'arm', 'JAWWNF010000348.1', 1559350, 1561206),
]

MONO = {'CLsat1': 146, 'CLsat2': 145, 'CLsat3': 146, 'CLsat4': 116}
FLANK = 600


def read_fa(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: ''.join(v).upper() for k, v in seqs.items()}


cons = read_fa(CONS)
print('consensuses: %s' % {k: len(v) for k, v in cons.items()})

# doubled subjects, one file per family
for fam, seq in cons.items():
    with open('%s/%s.doubled.fa' % (OUT, fam), 'w') as f:
        f.write('>%s\n%s\n' % (fam, seq + seq))

# locate genome fasta per assembly
GEN = {}
for sp in set(c[1] for c in CAND):
    hits = [p for p in
            subprocess.run(['bash', '-lc', 'ls %s/%s/*.fna 2>/dev/null' % (BASE, sp)],
                           capture_output=True, text=True).stdout.split()
            if p.endswith('.fna')]
    assert hits, 'no genome fasta for %s' % sp
    GEN[sp] = hits[0]
    print('%-4s %s' % (sp, hits[0]))


def faidx(sp, scaf, a, b):
    """1-based inclusive region"""
    r = subprocess.run(['samtools', 'faidx', GEN[sp], '%s:%d-%d' % (scaf, a, b)],
                       capture_output=True, text=True)
    return ''.join(l.strip() for l in r.stdout.splitlines()[1:]).upper()


def scaf_len(sp, scaf):
    for line in open(GEN[sp] + '.fai'):
        p = line.split('\t')
        if p[0] == scaf:
            return int(p[1])
    return None


def best_bits(query_seq, fam):
    """max bitscore of query vs doubled family consensus"""
    q = '%s/_q.fa' % OUT
    with open(q, 'w') as f:
        f.write('>q\n%s\n' % query_seq)
    r = subprocess.run([SS, '-n', '-m', '8C', '-E', '1000', '-b', '1',
                        q, '%s/%s.doubled.fa' % (OUT, fam)],
                       capture_output=True, text=True)
    best = 0.0
    for line in r.stdout.splitlines():
        if line.startswith('#'):
            continue
        p = line.split('\t')
        if len(p) >= 12:
            try:
                best = max(best, float(p[11]))
            except ValueError:
                pass
    return best


bits_rows = {f: [] for f in MONO}
summary = {}
flanks = {}

for cid, sp, scaf, a, b in CAND:
    slen = scaf_len(sp, scaf)
    seq = faidx(sp, scaf, a, b)
    print('\n%s  %s %s:%d-%d  len=%d  scaffold=%d' % (cid, sp, scaf, a, b, len(seq), slen))
    permax = {}
    for fam in ('CLsat1', 'CLsat2', 'CLsat3', 'CLsat4'):
        m = MONO[fam]
        win, step = 2 * m, m
        mx = 0.0
        n = 0
        pos = 0
        while pos + win <= len(seq):
            sub = seq[pos:pos + win]
            bt = best_bits(sub, fam)
            bits_rows[fam].append((scaf, a + pos, bt,
                                   '%s|%s:%d-%d|CAND' % (sp, scaf, a, b), fam, sp))
            mx = max(mx, bt)
            n += 1
            pos += step
        permax[fam] = mx
        print('   %-7s windows=%-3d max=%.1f' % (fam, n, mx))
    order = sorted(permax.items(), key=lambda kv: -kv[1])
    summary[cid] = dict(sp=sp, scf=scaf, start=a - 1, end=b, len=len(seq),
                        bits=permax, call=order[0][0], call_bits=order[0][1],
                        comp=order[1][0], comp_bits=order[1][1],
                        delta=round(order[0][1] - order[1][1], 1))
    # flanks
    lf_a, lf_b = max(1, a - FLANK), a - 1
    rf_a, rf_b = b + 1, min(slen, b + FLANK)
    lf = faidx(sp, scaf, lf_a, lf_b) if lf_b >= lf_a else ''
    rf = faidx(sp, scaf, rf_a, rf_b) if rf_b >= rf_a else ''
    flanks[cid] = dict(LF=dict(start=lf_a - 1, end=lf_b, len=len(lf), seq=lf,
                               truncated=(a - FLANK) < 1),
                       RF=dict(start=rf_a - 1, end=rf_b, len=len(rf), seq=rf,
                               truncated=(b + FLANK) > slen))

with open('%s/candidates_summary.json' % OUT, 'w') as f:
    json.dump(dict(summary=summary, flanks=flanks), f, indent=1)
for fam, rows in bits_rows.items():
    with open('%s/CAND_%s.bits.tsv' % (OUT, fam), 'w') as f:
        for r in rows:
            f.write('%s\t%d\t%.1f\t%s\t%s\t%s\n' % r)
print('\nwrote %s' % OUT)
for cid, s in summary.items():
    print('  %-3s %-4s %-8s call=%-7s %6.1f  comp=%-7s %6.1f  delta=%.1f'
          % (cid, s['sp'], s['len'], s['call'], s['call_bits'],
             s['comp'], s['comp_bits'], s['delta']))
