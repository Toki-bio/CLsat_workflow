# -*- coding: utf-8 -*-
"""CLEAN flank-orthology from the expanded loci, with per-pair FASTA alignments for manual audit.

Chain:
  1. loci = the expansion BEDs of every family, merged per species+scaffold
  2. for each locus, take flanks OUTSIDE the array; OVERLAP GUARD: a flank may not overlap any
     other locus (this is the bug that contaminated 35% of the old 'unique' flanks)
  3. classify each flank: satellite (matches a family consensus) / low-complexity/telomeric /
     unique. Only UNIQUE flanks are anchors.
  4. cross-assembly orthology: ssearch36 (NOT blastn), single best cross-assembly hit per
     target assembly per query, >=200 bp aligned, >=80% identity
  5. DELIVERABLE: for every orthologous pair, write a FASTA alignment of the two flanks
     (array->outward; UPPERCASE = residual satellite, lowercase = unique) so the shared unique
     sequence can be eyeballed. Plus the orthology matrix.

Usage:
  clean_ortho_pairs.py --expansions expansions.tsv --genomes genomes.tsv --families families.tsv
                       --out OUT_DIR [--flank 600] [--min-aln 200] [--min-pid 80] [--guard 0]
                       [--priority-pairs unp:nai,unp:val] [--max-alignments 120]
  expansions.tsv  family  expanded_merged.bed   (stage-11 output per family; col4 = species|...)
  genomes.tsv     species_code  genome.fa        (order = matrix order)
  families.tsv    family  mono  consensus.fa
Tools: ssearch36, samtools, needle (EMBOSS) on PATH.
Outputs in OUT_DIR: loci.bed  flank_class.tsv  ortho_pairs.tsv  ortho_matrix.txt
                    ortho_pairs_alignments.fasta   <-- the manual-inspection file
"""
import argparse
import collections
import os
import shutil
import subprocess
import sys
import tempfile


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


ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--expansions', required=True)
ap.add_argument('--genomes', required=True)
ap.add_argument('--families', required=True)
ap.add_argument('--out', required=True)
ap.add_argument('--flank', type=int, default=600)
ap.add_argument('--min-aln', type=int, default=200)
ap.add_argument('--min-pid', type=float, default=80.0)
ap.add_argument('--guard', type=int, default=0)
ap.add_argument('--priority-pairs', default='', help='species pairs a:b,... listed first in the alignment file')
ap.add_argument('--max-alignments', type=int, default=120)
a = ap.parse_args()

OUT = a.out
os.makedirs(OUT, exist_ok=True)
SS = shutil.which('ssearch36') or sys.exit('ssearch36 not on PATH')
SAM = shutil.which('samtools') or sys.exit('samtools not on PATH')
NEEDLE = shutil.which('needle') or sys.exit('needle (EMBOSS) not on PATH')

FLANK = a.flank           # bp of flank to test/align
MIN_ALN = a.min_aln
MIN_PID = a.min_pid
GUARD = a.guard           # a flank may not come within GUARD bp of another locus

GEN = {}
ASM = []
for p in read_tsv(a.genomes, 2):
    GEN[p[0]] = p[1]
    ASM.append(p[0])
CONS = {p[0]: ''.join(l.strip() for l in open(p[2]) if not l.startswith('>')).upper()
        for p in read_tsv(a.families, 3)}
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

# ---- 1. build the loci ----
raw = []
for fam, bed in ((p[0], p[1]) for p in read_tsv(a.expansions, 2)):
    for line in open(bed):
        p = line.rstrip('\n').split('\t')
        if len(p) >= 4:
            raw.append((p[3].split('|')[0], p[0], int(p[1]), int(p[2]), fam))
byscaf = collections.defaultdict(list)
for sp, scf, s, e, fam in raw:
    byscaf[(sp, scf)].append((s, e, fam))
loci = []
for (sp, scf), rs in byscaf.items():
    rs.sort()
    cs, ce, fams = rs[0][0], rs[0][1], {rs[0][2]}
    for s, e, fam in rs[1:]:
        if s <= ce:
            ce = max(ce, e); fams.add(fam)
        else:
            loci.append((sp, scf, cs, ce, fams)); cs, ce, fams = s, e, {fam}
    loci.append((sp, scf, cs, ce, fams))
loci.sort(key=lambda x: (x[0], x[1], x[2]))
with open(OUT + '/loci.bed', 'w') as f:
    for i, (sp, scf, s, e, fams) in enumerate(loci, 1):
        f.write('%s\t%d\t%d\tR%d\t%s\t%s\n' % (scf, s, e, i, sp, ','.join(sorted(fams))))
print('loci: %d (%.2f Mb)' % (len(loci), sum(e - s for _, _, s, e, _ in loci) / 1e6), flush=True)
lidx = collections.defaultdict(list)
for i, (sp, scf, s, e, fams) in enumerate(loci):
    lidx[(sp, scf)].append((s, e, i))


def faidx(sp, scf, s, e):
    s = max(1, s); e = min(SL[(sp, scf)], e)
    if e < s:
        return ''
    r = subprocess.run([SAM, 'faidx', GEN[sp], '%s:%d-%d' % (scf, s, e)], capture_output=True, text=True)
    return ''.join(l.strip() for l in r.stdout.splitlines() if not l.startswith('>')).upper() if not r.returncode else ''


# ---- 2/3. flanks with overlap guard + classification ----
def overlaps_other(sp, scf, fs, fe, self_i):
    # flank [fs, fe] is 1-based inclusive; loci [s2, e2) are 0-based half-open.
    # Convert the flank to 0-based half-open [fs-1, fe) before the interval test.
    for (s2, e2, j) in lidx[(sp, scf)]:
        if j != self_i and (fs - 1) - GUARD < e2 and fe + GUARD > s2:
            return True
    return False


def classify(seq):
    if len(seq) < MIN_ALN:
        return 'truncated'
    fu = seq.upper()
    if max(fu.count('CCCTAA'), fu.count('TTAGGG')) >= 8:
        return 'telomeric'
    if len(set(fu[i:i+4] for i in range(len(fu)-3))) / max(1, len(fu)-3) < 0.35:
        return 'lowcomplex'
    return 'unique'


flanks = {}          # (locus_i, side) -> (sp, seq_oriented_array->outward, genomic (scf,s,e))
fclass = {}
for i, (sp, scf, s, e, fams) in enumerate(loci):
    for side in ('L', 'R'):
        if side == 'L':
            fs, fe = s - FLANK + 1, s
        else:
            fs, fe = e + 1, e + FLANK
        if fs < 1 or fe > SL[(sp, scf)]:
            fclass[(i, side)] = 'truncated'; continue
        if overlaps_other(sp, scf, fs, fe, i):
            fclass[(i, side)] = 'satellite'; continue      # the guard: flank inside another locus
        seq = faidx(sp, scf, fs, fe)
        if not seq:
            fclass[(i, side)] = 'truncated'; continue
        # also drop if it matches a family consensus (residual satellite not caught by coords)
        oriented = seq if side == 'R' else rc(seq)     # array->outward
        kl = classify(seq)
        fclass[(i, side)] = kl
        if kl == 'unique':
            flanks[(i, side)] = (sp, oriented, (scf, fs, fe))
cc = collections.Counter(fclass.values())
print('flank classes: %s' % dict(cc), flush=True)
print('unique anchor flanks: %d' % len(flanks), flush=True)

# residual satellite check on the 'unique' set (ssearch36 vs consensus, drop if strong hit)
with open(OUT + '/uflanks.fa', 'w') as f:
    for (i, side), (sp, seq, g) in flanks.items():
        f.write('>%d_%s\n%s\n' % (i, side, seq))
# ssearch36 (NOT blastn) vs the doubled consensuses; drop a flank with any strong satellite hit
o = subprocess.run([SS, '-n', '-m', '8C', '-E', '1e-2', '-b', '1', '-d', '1',
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

# ---- 4. cross-assembly orthology (ssearch36), best hit per target assembly per query ----
with open(OUT + '/uflanks.fa', 'w') as f:
    for (i, side), (sp, seq, g) in flanks.items():
        f.write('>%d_%s\n%s\n' % (i, side, seq))
# all-vs-all with ssearch36 (both strands). library = the same file.
o = subprocess.run([SS, '-n', '-m', '8C', '-E', '1e-3', '-b', '2000', OUT + '/uflanks.fa',
                    OUT + '/uflanks.fa'], capture_output=True, text=True).stdout
best = collections.defaultdict(dict)     # query -> target_assembly -> (bits, target_id, pid, alen)
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
    if q == t or alen < MIN_ALN or pid < MIN_PID:
        continue
    qi = int(q.split('_')[0]); ti = int(t.split('_')[0])
    qsp, tsp = loci[qi][0], loci[ti][0]
    if qsp == tsp:
        continue
    if tsp not in best[q] or bits > best[q][tsp][0]:
        best[q][tsp] = (bits, t, pid, alen)
pairs = set()
detail = {}
for q, d in best.items():
    for tsp, (bits, t, pid, alen) in d.items():
        key = tuple(sorted((q, t)))
        pairs.add(key)
        detail[key] = (pid, alen, bits)

# ---- matrix ----
M = collections.Counter()
for (x, y) in pairs:
    M[tuple(sorted((loci[int(x.split('_')[0])][0], loci[int(y.split('_')[0])][0])))] += 1
with open(OUT + '/ortho_matrix.txt', 'w') as f:
    hdr = '      ' + ''.join('%6s' % x for x in ASM)
    f.write(hdr + '\n'); print('\n=== CLEAN ORTHOLOGY MATRIX ==='); print(hdr)
    for x in ASM:
        row = '%5s ' % x + ''.join('     -' if x == y else '%6d' % M[tuple(sorted((x, y)))] for y in ASM)
        f.write(row + '\n'); print(row)
    tot = 'total unique undirected pairs: %d' % len(pairs)
    f.write('\n' + tot + '\n'); print('\n' + tot)
    if a.priority_pairs:
        key = 'KEY: ' + '  '.join('%s-%s=%d' % (x, y, M[tuple(sorted((x, y)))])
                                  for x, y in (kp.split(':') for kp in a.priority_pairs.split(',')))
        f.write(key + '\n'); print(key)


# ---- 5. per-pair FASTA alignments for manual audit ----
def sw_align(x, y):
    """global alignment via EMBOSS needle; return the two aligned strings for display."""
    fa = tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False); fa.write('>a\n%s\n' % x); fa.close()
    fb = tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False); fb.write('>b\n%s\n' % y); fb.close()
    o = subprocess.run([NEEDLE, '-asequence', fa.name, '-bsequence', fb.name,
                        '-gapopen', '10', '-gapextend', '0.5', '-aformat', 'fasta',
                        '-stdout', '-auto'], capture_output=True, text=True).stdout
    os.unlink(fa.name); os.unlink(fb.name)
    seqs = []
    cur = []
    for line in o.splitlines():
        if line.startswith('>'):
            if cur:
                seqs.append(''.join(cur)); cur = []
        else:
            cur.append(line.strip())
    if cur:
        seqs.append(''.join(cur))
    return (seqs[0], seqs[1]) if len(seqs) >= 2 else (x, y)


PRIO = set(tuple(sorted(kp.split(':'))) for kp in a.priority_pairs.split(',') if kp)
pl = sorted(pairs)


def pr(k):
    sps = tuple(sorted((loci[int(k[0].split('_')[0])][0], loci[int(k[1].split('_')[0])][0])))
    return 0 if sps in PRIO else 1


pl.sort(key=pr)
with open(OUT + '/ortho_pairs_alignments.fasta', 'w') as f:
    f.write('; CLEAN flank-orthology pairs for manual validation.\n')
    f.write('; Each block: two orthologous UNIQUE flanks, aligned (needle). UPPERCASE input was\n')
    f.write('; the flank running array->outward; these are the unique-flank anchors (satellite\n')
    f.write('; and overlapping flanks already excluded). If the two sequences share long unique\n')
    f.write('; stretches, the orthology is real; if they only match in low-complexity, it is not.\n\n')
    shown = 0
    for k in pl:
        if shown >= a.max_alignments:
            break
        qi, qs = int(k[0].split('_')[0]), k[0].split('_')[1]
        ti, ts = int(k[1].split('_')[0]), k[1].split('_')[1]
        sa, sb = flanks.get((qi, qs)), flanks.get((ti, ts))
        if not sa or not sb:
            continue
        pid, alen, bits = detail[k]
        A, B = sw_align(sa[1], sb[1])
        mstr = ''.join('|' if x == y and x != '-' else ' ' for x, y in zip(A, B))
        f.write('>>> pair R%d(%s,%s) <-> R%d(%s,%s)  |  %.1f%% id over %d bp, bits %.0f\n'
                % (qi + 1, sa[0], qs, ti + 1, sb[0], ts, pid, alen, bits))
        f.write('    %s %s:%d-%d\n    %s %s:%d-%d\n' % (sa[0], sa[2][0], sa[2][1], sa[2][2],
                                                        sb[0], sb[2][0], sb[2][1], sb[2][2]))
        for x in range(0, len(A), 100):
            f.write('  %s  %s\n' % (sa[0][:3], A[x:x+100]))
            f.write('       %s\n' % mstr[x:x+100])
            f.write('  %s  %s\n\n' % (sb[0][:3], B[x:x+100]))
        shown += 1
print('\nwrote %s/ortho_pairs_alignments.fasta (%d pairs)' % (OUT, shown))
with open(OUT + '/flank_class.tsv', 'w') as f:
    f.write('locus\tside\tclass\n')
    for (i, side), kl in fclass.items():
        f.write('R%d\t%s\t%s\n' % (i + 1, side, kl))
print('wrote %s/ortho_matrix.txt, loci.bed, flank_class.tsv' % OUT)
