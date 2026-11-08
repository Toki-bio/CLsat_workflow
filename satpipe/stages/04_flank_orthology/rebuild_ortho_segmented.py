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
# single-copy anchors (flank_copy.py; author 2026-11-08). All off by default = the earlier behaviour exactly.
ap.add_argument('--trf-dat', action='append', default=[], metavar='CODE=PATH',
                help='genome-wide TRF .dat of a genome (repeatable; several per genome allowed): mask tandem letters')
ap.add_argument('--single-copy', action='store_true', help='anchor flanks must be single copy in their own genome')
ap.add_argument('--copy-min-pid', type=float, default=80.0, help='identity of a second copy (%%)')
ap.add_argument('--copy-min-aln', type=int, default=0, help='aligned letters of a second copy (0 = --min-aln)')
ap.add_argument('--clear-margin', type=float, default=0.0,
                help='best hit must beat the second-best locus: second < margin x best bits (0 = off, old best-hit rule)')
ap.add_argument('--reciprocal', action='store_true', help='keep a pair only if both loci choose each other')
ap.add_argument('--threads', type=int, default=8)
# refinements (author 2026-11-08, after the rules left molly with 84 and Darevskia with 7 anchors):
ap.add_argument('--repeat-mask', action='store_true',
                help='mask (N) the pieces of a flank that have another copy in the own genome (>= copy-min-pid, '
                     '>= mask-min-aln) and the pieces matching a catalogue family, instead of rejecting the whole '
                     'flank; the flank stays an anchor if >= min-aln letters are left (needs --single-copy)')
ap.add_argument('--mask-min-aln', type=int, default=100, help='shortest repeated piece that is masked')
ap.add_argument('--cut-at-locus', action='store_true',
                help='a flank ends before the next catalogue locus (and at the contig end) instead of being '
                     'rejected; kept if >= min-aln letters remain (for long --flank)')
a = ap.parse_args()
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flank_copy as FC  # noqa: E402
TRF_DAT = collections.defaultdict(list)
for x in a.trf_dat:
    if '=' not in x:
        sys.exit('--trf-dat wants CODE=PATH, got %r' % x)
    c, pth = x.split('=', 1)
    if not os.path.isfile(pth):
        sys.exit('--trf-dat %s: no such file' % pth)
    TRF_DAT[c].append(pth)
SINGLE = a.single_copy
COPY_MIN_PID = a.copy_min_pid
COPY_MIN_ALN = a.copy_min_aln or a.min_aln
CLEAR_MARGIN = a.clear_margin
RECIPROCAL = a.reciprocal
THREADS = a.threads
if RECIPROCAL and CLEAR_MARGIN <= 0:
    sys.exit('--reciprocal needs --clear-margin > 0')
REPMASK = a.repeat_mask
MASK_MIN_ALN = a.mask_min_aln
CUT = a.cut_at_locus
if REPMASK and not SINGLE:
    sys.exit('--repeat-mask needs --single-copy')

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
TR = str.maketrans('ACGTNRYKMSWBDHV', 'TGCANYRMKSWVHDB')   # IUPAC-aware (GLM audit #134)
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
        if CUT:
            fs, fe, why = FC.cut_flank(side, fs, fe, s, e, SL[(sp, scf)], lidx[(sp, scf)], i, GUARD)
            if fe - fs + 1 < MIN_ALN:
                fclass[(i, side)] = 'satellite' if why == 'locus' else 'truncated'
                continue
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
sat_masked = {}
if REPMASK:
    # mask the catalogue-family pieces instead of dropping the flank; ssearch36 gives one alignment per query x
    # family, so repeat until no piece is left (a flank can hold several separate pieces)
    for rnd_ in range(30):
        hit = collections.defaultdict(list)
        for line in o.splitlines():
            p = line.split('\t')
            if line.startswith('#') or len(p) < 12:
                continue
            try:
                alen, pid, qs, qe = int(p[3]), float(p[2]), int(p[6]), int(p[7])
            except ValueError:
                continue
            if alen >= 60 and pid >= 75:
                hit[p[0]].append((min(qs, qe), max(qs, qe)))
        if not hit:
            break
        for k in list(flanks):
            nm = '%d_%s' % k
            if nm in hit:
                sp, seq, g = flanks[k]
                mseq, n = FC.mask(seq, 1, FC.merge(hit[nm]))
                sat_masked[k] = sat_masked.get(k, 0) + n
                flanks[k] = (sp, mseq, g)
        with open(OUT + '/uflanks.fa', 'w') as f:
            for (i, side), (sp, seq, g) in flanks.items():
                f.write('>%d_%s\n%s\n' % (i, side, seq))
        o = subprocess.run([SS, '-n', '-z', '3', '-m', '8C', '-E', '1e-2', '-b', '4000', '-d', '1',
                            OUT + '/uflanks.fa', OUT + '/allcons.x2.fa'], capture_output=True, text=True).stdout
    sat_drop = set('%d_%s' % k for k, v in flanks.items() if FC.unmasked(v[1]) < MIN_ALN)
for k in list(flanks):
    if '%d_%s' % (k[0], k[1]) in sat_drop:
        del flanks[k]; fclass[k] = 'satellite'
print('after residual-satellite %s: %d unique anchors' % ('mask' if REPMASK else 'drop', len(flanks)), flush=True)

# ---- 3b. tandem mask from the run's genome-wide TRF calls (flank_copy.py rule 1) ----
masked_bp = {}
if TRF_DAT:
    need = collections.defaultdict(set)
    for (i, side), (sp, seq, (scf, fs, fe)) in flanks.items():
        need[sp].add(scf)
    TIV = {}
    for sp, paths in TRF_DAT.items():
        acc = collections.defaultdict(list)
        for pth in paths:
            for scf, iv in FC.parse_trf_dat(pth, need.get(sp, set())).items():
                acc[scf].extend(iv)
        TIV[sp] = {scf: FC.merge(iv) for scf, iv in acc.items()}
    miss = sorted(sp for sp in need if sp not in TRF_DAT)
    if miss:
        sys.exit('--trf-dat given but no TRF calls for genome(s): %s' % ', '.join(miss))
    for k in list(flanks):
        sp, oseq, (scf, fs, fe) = flanks[k]
        raw = flankseq[k][0]
        mseq, n = FC.mask(raw, fs, TIV[sp].get(scf, []))
        masked_bp[k] = n
        omseq = FC.n_union(mseq if k[1] == 'R' else rc(mseq), oseq)   # keep the satellite mask (--repeat-mask)
        if FC.unmasked(omseq) < MIN_ALN:
            del flanks[k]; fclass[k] = 'tandem'; continue
        flanks[k] = (sp, omseq, (scf, fs, fe))
    print('after tandem mask: %d anchors (%d flanks mostly tandem repeat)'
          % (len(flanks), sum(1 for v in fclass.values() if v == 'tandem')), flush=True)

# ---- 3c. repeat mask: pieces with another copy in the own genome become N (--repeat-mask) ----
rep_masked = {}
if REPMASK:
    MM = shutil.which('minimap2') or sys.exit('--repeat-mask needs minimap2 on PATH')
    by_sp = collections.defaultdict(list)
    for k, (sp, seq, g) in flanks.items():
        by_sp[sp].append(k)
    for sp in ASM:
        ks = by_sp.get(sp)
        if not ks:
            continue
        qfa = OUT + '/repmask_%s.fa' % sp
        place = {}
        with open(qfa, 'w') as f:
            for k in ks:
                nm = '%d_%s' % k
                f.write('>%s\n%s\n' % (nm, flanks[k][1]))
                place[nm] = flanks[k][2]
        r = subprocess.run([MM, '-k13', '-w4', '-s40', '-m25', '-A2', '-B4', '-O4,24', '-E2,1', '-c', '-N', '200', '-p', '0',
                            '-t', str(THREADS), GEN[sp], qfa], capture_output=True, text=True)
        if r.returncode:
            sys.exit('minimap2 failed on %s (exit %d): %s' % (sp, r.returncode, r.stderr.strip()[-300:]))
        piv = FC.paf_other_pieces(r.stdout.splitlines(), place, COPY_MIN_PID, MASK_MIN_ALN)
        for k in ks:
            nm = '%d_%s' % k
            sp_, seq, g = flanks[k]
            mseq, n = FC.mask(seq, 1, FC.merge(piv.get(nm, [])))
            rep_masked[k] = n
            if FC.unmasked(mseq) < MIN_ALN:
                del flanks[k]; fclass[k] = 'repeated'; continue
            flanks[k] = (sp_, mseq, g)
    print('after repeat mask: %d anchors (%d flanks mostly repeated pieces)'
          % (len(flanks), sum(1 for v in fclass.values() if v == 'repeated')), flush=True)

# ---- 4. cross-assembly orthology at three thresholds (for viewer 70/80/90) ----
with open(OUT + '/uflanks.fa', 'w') as f:
    for (i, side), (sp, seq, g) in flanks.items():
        f.write('>%d_%s\n%s\n' % (i, side, seq))
# threads only with --repeat-mask (long flanks); the published runs keep the single-thread call unchanged
o = subprocess.run([SS, '-n', '-z', '3', '-m', '8C', '-E', '1e-3', '-b', '4000']
                   + (['-T', str(THREADS)] if REPMASK else []) + [OUT + '/uflanks.fa', OUT + '/uflanks.fa'],
                   capture_output=True, text=True).stdout

# collect all cross-assembly hits with alen>=MIN_ALN (threshold on pid applied later per-level)
allhits = []   # (q, t, pid, alen, bits, side)
own_hits = collections.defaultdict(set)   # flank -> other loci of the SAME genome it matches (copy test, rule 2)
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
    qi, ti = int(q.split('_')[0]), int(t.split('_')[0])
    if loci[qi][0] == loci[ti][0]:
        # same genome: never an ortholog; with --single-copy it is evidence that the flank is not single copy
        if SINGLE and qi != ti and alen >= COPY_MIN_ALN and pid >= COPY_MIN_PID:
            # two flanks that share genome letters (loci < 2 x flank apart) are one place, not two copies
            gq, gt = flanks[(qi, q.split('_')[1])][2], flanks[(ti, t.split('_')[1])][2]
            if not (gq[0] == gt[0] and gq[1] <= gt[2] and gt[1] <= gq[2]):
                own_hits[q].add(ti)
        continue
    if q == t or alen < MIN_ALN:
        continue
    side = q.split('_')[1]
    # NOTE: no side restriction -- matches clean_ortho_pairs.py exactly.
    # Flanks are oriented array->outward, so an L flank legitimately matches an R flank
    # when two orthologous arrays are in opposite orientations (inversion).
    allhits.append((q, t, pid, alen, bits, side))

# ---- 4b. single copy in the own whole genome (flank_copy.py rule 2) ----
copies = {}    # (i, side) -> (other loci hit, other genome places, status)
if SINGLE:
    MM = shutil.which('minimap2') or sys.exit('--single-copy needs minimap2 on PATH')
    for k in list(flanks):
        nm = '%d_%s' % k
        copies[k] = [len(own_hits.get(nm, ())), None, '']
    by_sp = collections.defaultdict(list)
    for k, (sp, seq, g) in flanks.items():
        if not copies[k][0]:                      # already repeated among locus flanks: no genome search needed
            by_sp[sp].append(k)
    for sp in ASM:
        ks = by_sp.get(sp)
        if not ks:
            continue
        qfa = OUT + '/copytest_%s.fa' % sp
        place = {}
        with open(qfa, 'w') as f:
            for k in ks:
                nm = '%d_%s' % k
                f.write('>%s\n%s\n' % (nm, flanks[k][1]))
                scf, fs, fe = flanks[k][2]
                place[nm] = (scf, fs, fe)
        # sensitive seeding, every secondary kept: bench 2026-11-08 (tests/flank_copy_toy/bench_copy_search.py, 600-letter
        # flanks, copies with substitutions + indels): this finds 20/20 copies at 90, 85 and 82 %, 14/20 at 80 %, 0/20
        # false at 75 %; '-x asm20 -p 0.2' found 14/20 at 90 % and 1/20 at 85 %.
        r = subprocess.run([MM, '-k13', '-w4', '-s40', '-m25', '-A2', '-B4', '-O4,24', '-E2,1', '-c', '-N', '200', '-p', '0',
                            '-t', str(THREADS), GEN[sp], qfa], capture_output=True, text=True)
        if r.returncode:
            sys.exit('minimap2 failed on %s (exit %d): %s' % (sp, r.returncode, r.stderr.strip()[-300:]))
        n = FC.paf_other_copies(r.stdout.splitlines(), place, COPY_MIN_PID, COPY_MIN_ALN)
        for k in ks:
            copies[k][1] = n.get('%d_%s' % k, 0)
    for k in list(flanks):
        c = copies[k]
        if c[0] or (c[1] or 0):
            c[2] = 'repeated'
            del flanks[k]; fclass[k] = 'repeated'
        else:
            c[2] = 'single'
    anchor = set('%d_%s' % k for k in flanks)
    allhits = [h for h in allhits if h[0] in anchor and h[1] in anchor]
    print('after single-copy test: %d anchors (%d repeated in their own genome)'
          % (len(flanks), sum(1 for v in fclass.values() if v == 'repeated')), flush=True)

# MIN_PID pairs (the manuscript/Fig-3 set), single best per target assembly per query
pairs = set()
detail = {}
oneway = set()
if CLEAR_MARGIN > 0:
    # rules 3 and 4: clear best locus per query flank and target genome, then reciprocal
    cand = collections.defaultdict(list)
    for q, t, pid, alen, bits, side in allhits:
        if pid < MIN_PID:
            continue
        tsp = loci[int(t.split('_')[0])][0]
        cand[(q, tsp)].append((bits, lid_of[int(t.split('_')[0])], pid, alen))
    directed = {}
    for (q, tsp), hits in sorted(cand.items()):
        b = FC.clear_best(hits, CLEAR_MARGIN)
        if b is None:
            continue
        qL = lid_of[int(q.split('_')[0])]
        d = (qL, b[1])
        if d not in directed or b[0] > directed[d][0]:
            directed[d] = (b[0], b[2], b[3])
    if RECIPROCAL:
        keep, oneway = FC.reciprocal_pairs(set(directed))
    else:
        keep = set(tuple(sorted(d)) for d in directed)
    for key in keep:
        v = max((directed[d] for d in (key, key[::-1]) if d in directed), key=lambda x: x[0])
        pairs.add(key)
        detail[key] = (v[1], v[2], v[0])
else:
    best = collections.defaultdict(dict)
    for q, t, pid, alen, bits, side in allhits:
        if pid < MIN_PID:
            continue
        tsp = loci[int(t.split('_')[0])][0]
        if tsp not in best[q] or bits > best[q][tsp][0]:
            best[q][tsp] = (bits, t, pid, alen)
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

if TRF_DAT or SINGLE:
    # per anchor-candidate flank: tandem letters masked, other loci / other genome places hit, verdict
    with open(OUT + '/flank_copies.tsv', 'w') as f:
        f.write('locus\tside\tgenome\tmasked_tandem_bp\tother_locus_flanks\tother_genome_places\tclass'
                + ('\tflank_bp\tmasked_satellite_bp\tmasked_repeat_bp\tleft_bp' if REPMASK else '') + '\n')
        for (i, side) in sorted(set(masked_bp) | set(copies) | set(rep_masked), key=lambda k: (lid_of[k[0]], k[1])):
            c = copies.get((i, side), ['', '', ''])
            f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s' % (lid_of[i], side, loci[i][0], masked_bp.get((i, side), ''),
                                                c[0], '' if c[1] is None else c[1], fclass[(i, side)]))
            if REPMASK:
                fq = flankseq.get((i, side))
                left = FC.unmasked(flanks[(i, side)][1]) if (i, side) in flanks else 0
                f.write('\t%s\t%s\t%s\t%s' % (len(fq[0]) if fq else '', sat_masked.get((i, side), ''),
                                               rep_masked.get((i, side), ''), left))
            f.write('\n')
if RECIPROCAL:
    with open(OUT + '/ortho_oneway.tsv', 'w') as f:
        f.write('locusA\tlocusB\tspA\tspB\n')
        for x, y in sorted(oneway):
            f.write('%s\t%s\t%s\t%s\n' % (x, y, lid2sp[x], lid2sp[y]))
    print('one-way (not reciprocal) clear best hits, dropped: %d' % len(oneway))

print('\nwrote flanks.json, flank_orthologs.json, flank_class.tsv, ortho_pairs.tsv, ortho_matrix.txt in %s' % OUT)
print('DONE_ORTHO')
