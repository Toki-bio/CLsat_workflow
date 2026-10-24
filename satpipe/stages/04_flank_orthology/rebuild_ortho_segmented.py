# -*- coding: utf-8 -*-
"""Flank orthology on the segmented catalogue (L-ids), so that flank orthology, sub-subfamily
assignment, Table 1 and the viewer are all keyed to one consistent catalogue/id scheme.
Core method matches clean_ortho_pairs.py (600 bp flanks, overlap guard, residual-satellite drop,
ssearch36 >=80% id / >=200 bp, single best hit per target assembly, no side restriction). It is
not byte-for-byte identical: the loci source is the segmented catalogue rather than rebuilt
raw-merged loci, pairs are keyed by locus id rather than by flank id, the ssearch36 -b cutoff
differs (4000 here), and the output format differs. The resulting 53-pair / 18-both-flank-group
matrix is identical either way for the Darevskia catalogue (verified against the raw/merged
539-locus run).

Usage:
  rebuild_ortho_segmented.py --loci locus_master.tsv --genomes genomes.tsv --families families.tsv
                             --out OUT_DIR [--flank 600] [--min-aln 200] [--min-pid 80] [--guard 0]
                             [--key-pairs unp:nai,unp:val]
  locus_master.tsv  lid  scaffold  start0  end  species  family[,family...]  [mono]
  genomes.tsv       species_code  genome.fa   (order = row/column order of the matrix)
  families.tsv      family  mono  consensus.fa
Tools: ssearch36 and samtools must be on PATH.
Outputs: flanks.json, flank_orthologs.json, flank_class.tsv, ortho_pairs.tsv, ortho_matrix.txt
"""
import argparse
import collections
import json
import os
import shutil
import subprocess
import sys


def read_tsv(path, min_cols):
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < min_cols:
                sys.exit('%s: need %d columns, got %r' % (path, min_cols, line.strip()))
            rows.append(p)
    return rows


def read_fasta_seq(path):
    return ''.join(l.strip() for l in open(path) if not l.startswith('>')).upper()


ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--loci', required=True)
ap.add_argument('--genomes', required=True)
ap.add_argument('--families', required=True)
ap.add_argument('--out', required=True)
ap.add_argument('--flank', type=int, default=600)
ap.add_argument('--min-aln', type=int, default=200)
ap.add_argument('--min-pid', type=float, default=80.0)
ap.add_argument('--guard', type=int, default=0)
ap.add_argument('--key-pairs', default='', help='comma list of species pairs a:b to print as a KEY line')
a = ap.parse_args()

OUT = a.out
os.makedirs(OUT, exist_ok=True)
SS = shutil.which('ssearch36') or sys.exit('ssearch36 not on PATH')
SAM = shutil.which('samtools') or sys.exit('samtools not on PATH')

FLANK = a.flank
MIN_ALN = a.min_aln
MIN_PID = a.min_pid
GUARD = a.guard

GEN = {}
ASM = []
for p in read_tsv(a.genomes, 2):
    GEN[p[0]] = p[1]
    ASM.append(p[0])
MONO = {}
CONS = {}
for p in read_tsv(a.families, 3):
    MONO[p[0]] = int(p[1])
    CONS[p[0]] = read_fasta_seq(p[2])
with open(OUT + '/allcons.x2.fa', 'w') as f:
    for fam, s in CONS.items():
        f.write('>%s\n%s%s\n' % (fam, s, s))
SL = {}
for sp, g in GEN.items():
    for line in open(g + '.fai'):
        p = line.split('\t')
        SL[(sp, p[0])] = int(p[1])
TR = str.maketrans('ACGTN', 'TGCAN')
rc = lambda s: s.translate(TR)[::-1]

# ---- 1. load loci from the segmented catalogue, keep L-ids ----
loci = []       # (sp, scf, s, e, fams, lid, mono)
for p in read_tsv(a.loci, 6):
    lid, scf, s, e, sp, fam = p[0], p[1], int(p[2]), int(p[3]), p[4], p[5]
    fam0 = fam.split(',')[0]
    if fam0 not in MONO:
        sys.exit('locus %s: family %r is not in %s' % (lid, fam0, a.families))
    if sp not in GEN:
        sys.exit('locus %s: species %r is not in %s' % (lid, sp, a.genomes))
    loci.append((sp, scf, s, e, fam, lid, MONO[fam0]))
loci.sort(key=lambda x: (x[0], x[1], x[2]))
lid_of = {i: loci[i][5] for i in range(len(loci))}
print('loci: %d (%.2f Mb)' % (len(loci), sum(e - s for _, _, s, e, *_ in loci) / 1e6), flush=True)
lidx = collections.defaultdict(list)
for i, (sp, scf, s, e, fam, lid, mono) in enumerate(loci):
    lidx[(sp, scf)].append((s, e, i))


def faidx(sp, scf, s, e):
    s = max(1, s); e = min(SL[(sp, scf)], e)
    if e < s:
        return ''
    r = subprocess.run([SAM, 'faidx', GEN[sp], '%s:%d-%d' % (scf, s, e)], capture_output=True, text=True)
    seq = ''.join(l.strip() for l in r.stdout.splitlines() if not l.startswith('>')).upper()
    # a failed fetch of a region inside the contig is an error; it used to be filed as a 'truncated' flank,
    # silently losing an anchor (GLM audit 2026-10-24)
    if r.returncode or not seq:
        sys.exit('samtools faidx %s %s:%d-%d failed (exit %d): %s' % (GEN[sp], scf, s, e, r.returncode, r.stderr.strip()[-300:]))
    return seq


def overlaps_other(sp, scf, fs, fe, self_i):
    # flank [fs, fe] is 1-based inclusive; loci [s2, e2) are 0-based half-open.
    for (s2, e2, j) in lidx[(sp, scf)]:
        if j != self_i and (fs - 1) - GUARD < e2 and fe + GUARD > s2:
            return True
    return False


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flank_rules  # noqa: E402  (shared truncated / telomeric / lowcomplex / unique rules)


def classify(seq):
    return flank_rules.classify(seq, MIN_ALN)


# ---- 2/3. flanks with overlap guard + classification ----
flanks = {}          # (i, side) -> (sp, oriented_seq, (scf,fs,fe))
fclass = {}
flankseq = {}        # (i, side) -> raw genomic seq (for viewer flanks.json, not oriented)
for i, (sp, scf, s, e, fam, lid, mono) in enumerate(loci):
    for side in ('L', 'R'):
        if side == 'L':
            fs, fe = s - FLANK + 1, s
        else:
            fs, fe = e + 1, e + FLANK
        if fs < 1 or fe > SL[(sp, scf)]:
            fclass[(i, side)] = 'truncated'; continue
        if overlaps_other(sp, scf, fs, fe, i):
            fclass[(i, side)] = 'satellite'; continue
        seq = faidx(sp, scf, fs, fe)
        if not seq:
            fclass[(i, side)] = 'truncated'; continue
        oriented = seq if side == 'R' else rc(seq)
        kl = classify(seq)
        fclass[(i, side)] = kl
        flankseq[(i, side)] = (seq, fs, fe)
        if kl == 'unique':
            flanks[(i, side)] = (sp, oriented, (scf, fs, fe))
cc = collections.Counter(fclass.values())
print('flank classes: %s' % dict(cc), flush=True)
print('unique anchor flanks (pre-residual-drop): %d' % len(flanks), flush=True)

# residual-satellite drop
with open(OUT + '/uflanks.fa', 'w') as f:
    for (i, side), (sp, seq, g) in flanks.items():
        f.write('>%d_%s\n%s\n' % (i, side, seq))
o = subprocess.run([SS, '-n', '-z', '3', '-m', '8C', '-E', '1e-2', '-b', '1', '-d', '1',
                    OUT + '/uflanks.fa', OUT + '/allcons.x2.fa'], capture_output=True, text=True).stdout
sat_drop = set()
for line in o.splitlines():
    if line.startswith('#'):
        continue
    p = line.split('\t')
    if len(p) >= 12:
        try:
            alen, pid = int(p[3]), float(p[2])
        except ValueError:
            continue
        if alen >= 60 and pid >= 75:
            sat_drop.add(p[0])
for k in list(flanks):
    if '%d_%s' % (k[0], k[1]) in sat_drop:
        del flanks[k]; fclass[k] = 'satellite'
print('after residual-satellite drop: %d unique anchors' % len(flanks), flush=True)

# ---- 4. cross-assembly orthology at three thresholds (for viewer 70/80/90) ----
with open(OUT + '/uflanks.fa', 'w') as f:
    for (i, side), (sp, seq, g) in flanks.items():
        f.write('>%d_%s\n%s\n' % (i, side, seq))
o = subprocess.run([SS, '-n', '-z', '3', '-m', '8C', '-E', '1e-3', '-b', '4000', OUT + '/uflanks.fa',
                    OUT + '/uflanks.fa'], capture_output=True, text=True).stdout

# collect all cross-assembly hits with alen>=MIN_ALN (threshold on pid applied later per-level)
allhits = []   # (q, t, pid, alen, bits, side)
for line in o.splitlines():
    if line.startswith('#'):
        continue
    p = line.split('\t')
    if len(p) < 12:
        continue
    q, t = p[0], p[1]
    try:
        pid, alen, bits = float(p[2]), int(p[3]), float(p[11])
    except ValueError:
        continue
    if q == t or alen < MIN_ALN:
        continue
    qi, ti = int(q.split('_')[0]), int(t.split('_')[0])
    if loci[qi][0] == loci[ti][0]:
        continue
    side = q.split('_')[1]
    # NOTE: no side restriction -- matches clean_ortho_pairs.py exactly.
    # Flanks are oriented array->outward, so an L flank legitimately matches an R flank
    # when two orthologous arrays are in opposite orientations (inversion).
    allhits.append((q, t, pid, alen, bits, side))

# MIN_PID pairs (the manuscript/Fig-3 set), single best per target assembly per query
best = collections.defaultdict(dict)
for q, t, pid, alen, bits, side in allhits:
    if pid < MIN_PID:
        continue
    tsp = loci[int(t.split('_')[0])][0]
    if tsp not in best[q] or bits > best[q][tsp][0]:
        best[q][tsp] = (bits, t, pid, alen)
pairs = set()
detail = {}
for q, d in best.items():
    for tsp, (bits, t, pid, alen) in d.items():
        # undirected pair keyed by locus id (strip _side)
        qL, tL = lid_of[int(q.split('_')[0])], lid_of[int(t.split('_')[0])]
        key = tuple(sorted((qL, tL)))
        pairs.add(key)
        detail[key] = (pid, alen, bits)

# ---- matrix ----
M = collections.Counter()
lid2sp = {loci[i][5]: loci[i][0] for i in range(len(loci))}
for (x, y) in pairs:
    M[tuple(sorted((lid2sp[x], lid2sp[y])))] += 1
with open(OUT + '/ortho_matrix.txt', 'w') as f:
    hdr = '      ' + ''.join('%6s' % x for x in ASM)
    f.write(hdr + '\n'); print('\n=== ORTHOLOGY MATRIX (%d segmented loci) ===' % len(loci)); print(hdr)
    for x in ASM:
        row = '%5s ' % x + ''.join('     -' if x == y else '%6d' % M[tuple(sorted((x, y)))] for y in ASM)
        f.write(row + '\n'); print(row)
    tot = 'total unique undirected pairs: %d' % len(pairs)
    f.write('\n' + tot + '\n'); print('\n' + tot)
    if a.key_pairs:
        key = 'KEY: ' + '  '.join('%s-%s=%d' % (x, y, M[tuple(sorted((x, y)))])
                                  for x, y in (kp.split(':') for kp in a.key_pairs.split(',')))
        f.write(key + '\n'); print(key)

# ---- dump machine-readable artifacts for the viewer builder ----
# flank_class.tsv (L-id, side, class)
with open(OUT + '/flank_class.tsv', 'w') as f:
    f.write('locus\tside\tclass\n')
    for (i, side), kl in fclass.items():
        f.write('%s\t%s\t%s\n' % (lid_of[i], side, kl))
# per-locus flank sequences + truncation (for flanks.json); flank start/end are 1-based inclusive
flankdump = {}
for i, (sp, scf, s, e, fam, lid, mono) in enumerate(loci):
    rec = {'sp': sp, 'scf': scf, 'start': s, 'end': e, 'subtype': fam.split(',')[0],
           'mono': mono, 'length': e - s, 'n_monomers': round((e - s) / mono)}
    for side, key in (('LF', 'L'), ('RF', 'R')):
        fq = flankseq.get((i, key))
        cls = fclass.get((i, key), 'truncated')
        if fq:
            seq, fs, fe = fq
            rec[side] = {'start': fs, 'end': fe, 'len': len(seq), 'seq': seq,
                         'truncated': cls == 'truncated', 'class': cls}
        else:
            rec[side] = {'start': None, 'end': None, 'len': 0, 'seq': '',
                         'truncated': True, 'class': cls}
    flankdump[lid] = rec
json.dump(flankdump, open(OUT + '/flanks.json', 'w'))

# ortho at 3 thresholds keyed LF/RF -> thr -> {locus:[{locus,sp,pid,aln,bits}]}
ortho_json = {'LF': {'70': {}, '80': {}, '90': {}}, 'RF': {'70': {}, '80': {}, '90': {}},
              'min_aln': MIN_ALN, 'thresholds': [70, 80, 90]}
for q, t, pid, alen, bits, side in allhits:
    sidekey = 'LF' if side == 'L' else 'RF'
    qL = lid_of[int(q.split('_')[0])]
    tL = lid_of[int(t.split('_')[0])]
    tsp = loci[int(t.split('_')[0])][0]
    for thr in (70, 80, 90):
        if pid >= thr:
            ortho_json[sidekey][str(thr)].setdefault(qL, []).append(
                {'locus': tL, 'sp': tsp, 'pid': round(pid, 1), 'aln': alen, 'bits': round(bits, 1)})
json.dump(ortho_json, open(OUT + '/flank_orthologs.json', 'w'))

# save the MIN_PID undirected pairs list
with open(OUT + '/ortho_pairs.tsv', 'w') as f:
    f.write('locusA\tlocusB\tspA\tspB\tpid\taln\tbits\n')
    for (x, y) in sorted(pairs):
        pid, alen, bits = detail[(x, y)]
        f.write('%s\t%s\t%s\t%s\t%.1f\t%d\t%.1f\n' % (x, y, lid2sp[x], lid2sp[y], pid, alen, bits))

print('\nwrote flanks.json, flank_orthologs.json, flank_class.tsv, ortho_pairs.tsv, ortho_matrix.txt in %s' % OUT)
print('DONE_ORTHO')
