#!/usr/bin/env python3
"""Process candidate loci that are absent from the catalogue (singleton loci with a strong
family signal), using the same steps the pipeline applies to MASTER loci, so they can be
inspected in the viewer alongside catalogued loci.

For each candidate:
  * extract the locus plus FLANK bp of flank on each side
  * sliding-window ssearch36 scan against every family consensus
    (window = 2 monomers, step = 1 monomer, per-family monomer length)
  * per-locus maximum bitscore per family -> family call and competitor
  * flank extraction

These are CANDIDATES for inspection, not catalogue members. They are given C-prefixed
identifiers and a distinct flag so they can never be confused with the catalogue.

Usage: candidate_scan.py <candidates.tsv> <genomes.tsv> <families.tsv> <out_dir> [--flank 600]
  candidates.tsv  id  species  scaffold  start1  end1   (1-based inclusive; '#' lines ignored)
  genomes.tsv     species_code  genome.fa (indexed)
  families.tsv    family  mono  consensus.fa
Tools: ssearch36, samtools on PATH.
Outputs in out_dir: candidates_summary.json, CAND_<family>.bits.tsv, <family>.doubled.fa
Note: flank start/end in candidates_summary.json are 0-based (BED frame), unlike the L-locus
flanks.json from rebuild_ortho_segmented.py (1-based) -- see docs/files.md.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

ap = argparse.ArgumentParser()
ap.add_argument('candidates')
ap.add_argument('genomes')
ap.add_argument('families')
ap.add_argument('out')
ap.add_argument('--flank', type=int, default=600)
a = ap.parse_args()

SS = shutil.which('ssearch36') or sys.exit('ssearch36 not on PATH')
SAM = shutil.which('samtools') or sys.exit('samtools not on PATH')
OUT = a.out
FLANK = a.flank
os.makedirs(OUT, exist_ok=True)


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


CAND = [(p[0], p[1], p[2], int(p[3]), int(p[4])) for p in read_tsv(a.candidates, 5)]
GEN = {p[0]: p[1] for p in read_tsv(a.genomes, 2)}
MONO = {}
cons = {}
for p in read_tsv(a.families, 3):
    MONO[p[0]] = int(p[1])
    cons[p[0]] = ''.join(read_fa(p[2]).values())
FAMS = list(MONO)
print('consensuses: %s' % {k: len(v) for k, v in cons.items()})

# doubled subjects, one file per family
for fam, seq in cons.items():
    with open('%s/%s.doubled.fa' % (OUT, fam), 'w') as f:
        f.write('>%s\n%s\n' % (fam, seq + seq))

for cid, sp, scaf, x, y in CAND:
    if sp not in GEN:
        sys.exit('candidate %s: species %r not in %s' % (cid, sp, a.genomes))


def run(cmd):
    """Run a tool; a failure stops the stage with the tool's message (it used to give an empty result silently)."""
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit('%s failed (exit %d): %s' % (' '.join(str(c) for c in cmd), r.returncode, r.stderr.strip()[-500:]))
    return r


def fetch(sp, scaf, x, y):
    """1-based inclusive region; an empty result is an error (missing contig or coordinates past its end)"""
    r = run([SAM, 'faidx', GEN[sp], '%s:%d-%d' % (scaf, x, y)])
    seq = ''.join(l.strip() for l in r.stdout.splitlines()[1:]).upper()
    if not seq:
        sys.exit('samtools faidx %s %s:%d-%d returned no sequence: %s' % (GEN[sp], scaf, x, y, r.stderr.strip()[-300:]))
    return seq


def faidx(sp, scaf, x, y):
    """1-based inclusive region"""
    return fetch(sp, scaf, x, y)


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
    r = run([SS, '-n', '-m', '8C', '-E', '1000', '-b', '1', q, '%s/%s.doubled.fa' % (OUT, fam)])
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

for cid, sp, scaf, x, y in CAND:
    slen = scaf_len(sp, scaf)
    if slen is None:
        sys.exit('candidate %s: scaffold %s not in %s.fai' % (cid, scaf, GEN[sp]))
    seq = faidx(sp, scaf, x, y)
    print('\n%s  %s %s:%d-%d  len=%d  scaffold=%d' % (cid, sp, scaf, x, y, len(seq), slen))
    permax = {}
    for fam in FAMS:
        m = MONO[fam]
        win, step = 2 * m, m
        mx = 0.0
        n = 0
        pos = 0
        while pos + win <= len(seq):
            sub = seq[pos:pos + win]
            bt = best_bits(sub, fam)
            bits_rows[fam].append((scaf, x + pos, bt,
                                   '%s|%s:%d-%d|CAND' % (sp, scaf, x, y), fam, sp))
            mx = max(mx, bt)
            n += 1
            pos += step
        permax[fam] = mx
        print('   %-7s windows=%-3d max=%.1f' % (fam, n, mx))
    order = sorted(permax.items(), key=lambda kv: -kv[1])
    comp = order[1] if len(order) > 1 else (None, 0.0)
    summary[cid] = dict(sp=sp, scf=scaf, start=x - 1, end=y, len=len(seq),
                        bits=permax, call=order[0][0], call_bits=order[0][1],
                        comp=comp[0], comp_bits=comp[1],
                        delta=round(order[0][1] - comp[1], 1))
    # flanks (0-based start in the JSON)
    lf_a, lf_b = max(1, x - FLANK), x - 1
    rf_a, rf_b = y + 1, min(slen, y + FLANK)
    lf = faidx(sp, scaf, lf_a, lf_b) if lf_b >= lf_a else ''
    rf = faidx(sp, scaf, rf_a, rf_b) if rf_b >= rf_a else ''
    flanks[cid] = dict(LF=dict(start=lf_a - 1, end=lf_b, len=len(lf), seq=lf,
                               truncated=(x - FLANK) < 1),
                       RF=dict(start=rf_a - 1, end=rf_b, len=len(rf), seq=rf,
                               truncated=(y + FLANK) > slen))

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
