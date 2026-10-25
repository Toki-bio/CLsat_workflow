#!/usr/bin/env python3
"""Sub-subfamily assignment for the candidate loci, following Methods 3.6:
a two-monomer window from the array boundary is aligned (ssearch36 -m 8C) against the
doubled sub-subfamily consensuses of the assigned family only; highest bitscore wins,
and the margin to the second best is recorded.

Usage: candidate_subsubfamily.py <out_dir> <genomes.tsv> <families.tsv> <subfamily_consensuses> <members.tsv>
                                 [--confident-delta 10]
  out_dir                holds candidates_summary.json from candidate_scan.py (updated in place)
  genomes.tsv            species_code  genome.fa
  families.tsv           family  mono  consensus.fa
  subfamily_consensuses  one multi-record FASTA, or a directory of <name>.fa files (subx/);
                         sequences are doubled here, so give SINGLE monomers in the FASTA case
                         and note that the Darevskia subx/ files are already doubled (then pass
                         --already-doubled)
  members.tsv            family<TAB>sub-subfamily
Tools: ssearch36, samtools on PATH.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

ap = argparse.ArgumentParser()
ap.add_argument('out')
ap.add_argument('genomes')
ap.add_argument('families')
ap.add_argument('subcons')
ap.add_argument('members')
ap.add_argument('--confident-delta', type=float, default=10.0)
ap.add_argument('--already-doubled', action='store_true')
a = ap.parse_args()

SS = shutil.which('ssearch36') or sys.exit('ssearch36 not on PATH')
SAM = shutil.which('samtools') or sys.exit('samtools not on PATH')
OUT = a.out



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


def read_fa(p):
    s, n = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith('>'):
            n = l[1:].split()[0]
            s[n] = []
        elif n:
            s[n].append(l)
    return {k: ''.join(v).upper() for k, v in s.items()}


GEN = {p[0]: p[1] for p in read_tsv(a.genomes, 2)}
MONO = {p[0]: int(p[1]) for p in read_tsv(a.families, 3)}
FAMOF = {}
for p in read_tsv(a.members, 2):
    if p[1] in FAMOF and FAMOF[p[1]] != p[0]:   # the same sub-subfamily name under two families: last one won (#97)
        sys.exit('%s: sub-subfamily %s listed under %s and %s' % (a.members, p[1], FAMOF[p[1]], p[0]))
    FAMOF[p[1]] = p[0]

ss = {}
if os.path.isdir(a.subcons):
    for fn in sorted(os.listdir(a.subcons)):
        if fn.endswith(('.fa', '.fas', '.fasta')):
            name = fn.rsplit('.', 1)[0]
            recs = read_fa(os.path.join(a.subcons, fn))
            if len(recs) != 1:   # one sub-subfamily per file; several records were concatenated and halved as one (#96)
                sys.exit('%s: a sub-subfamily file must hold exactly one record, found %d' % (fn, len(recs)))
            seq = next(iter(recs.values()))
            ss[name] = seq[:len(seq) // 2] if a.already_doubled else seq
else:
    ss = read_fa(a.subcons)
    if a.already_doubled:
        ss = {k: v[:len(v) // 2] for k, v in ss.items()}
print('sub-subfamily consensuses: %s' % {k: len(v) for k, v in ss.items()})
for name in ss:
    if name not in FAMOF:
        print('  WARNING: %s has no family in %s; ignored' % (name, a.members), file=sys.stderr)

# doubled subjects
for name, seq in ss.items():
    with open('%s/ss_%s.fa' % (OUT, name), 'w') as f:
        f.write('>%s\n%s\n' % (name, seq + seq))

data = json.load(open('%s/candidates_summary.json' % OUT))
summary = data['summary']

for cid, s in sorted(summary.items()):
    sp, scaf, x, y = s['sp'], s['scf'], s['start'] + 1, s['end']
    fam = s['call']
    mono = MONO[fam]
    win = min(2 * mono, y - x + 1)
    q = fetch(sp, scaf, x, x + win - 1)
    with open('%s/_ssq.fa' % OUT, 'w') as f:
        f.write('>q\n%s\n' % q)
    scores = {}
    for name in ss:
        if FAMOF.get(name) != fam:
            continue
        o = run([SS, '-n', '-z', '3', '-m', '8C', '-E', '1000', '-b', '1', '%s/_ssq.fa' % OUT, '%s/ss_%s.fa' % (OUT, name)]).stdout
        best = 0.0
        pid = 0.0
        for line in o.splitlines():
            if line.startswith('#'):
                continue
            p = line.split('\t')
            if len(p) >= 12:
                try:
                    if float(p[11]) > best:
                        best = float(p[11])
                        pid = float(p[2])
                except ValueError:
                    pass
        scores[name] = (best, pid)
    if not scores:
        continue
    order = sorted(scores.items(), key=lambda kv: -kv[1][0])
    top, (tb, tp) = order[0]
    second, sb = (order[1][0], order[1][1][0]) if len(order) > 1 else (None, 0.0)
    delta = round(tb - sb, 1)
    s['subfamily'] = top
    s['ss_bits'] = tb
    s['ss_pid'] = tp
    s['ss_second'] = second
    s['ss_delta'] = delta
    s['ss_confident'] = bool(delta >= a.confident_delta or second is None)
    print('  %-3s %-7s -> %-9s %6.1f bits  pid=%.1f  second=%-9s delta=%5.1f  confident=%s'
          % (cid, fam, top, tb, tp, second, delta, s['ss_confident']))

with open('%s/candidates_summary.json' % OUT, 'w') as f:
    json.dump(data, f, indent=1)
print('\nupdated candidates_summary.json')
