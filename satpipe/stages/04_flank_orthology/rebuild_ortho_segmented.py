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
ap.add_argument('--extend-steps', default='',
                help='comma list of window lengths, e.g. 600,3000,5000: a flank not resolved at one length is tried at '
                     'the next (SINE_orth_loc adapt_flanks idea); needs --repeat-mask, implies --cut-at-locus')
ap.add_argument('--joint', action='store_true',
                help='second tier for loci left without a pair: both flanks of a locus must hit opposite flanks of the '
                     'same target locus (summed bits, clear best, reciprocal), and that left+right combination must '
                     'occur at no other locus of its own genome (SINE_orth_loc both_flanks idea)')
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
STEPS = [int(x) for x in a.extend_steps.split(',') if x.strip()] if a.extend_steps else [a.flank]
if a.extend_steps:
    if not REPMASK:
        sys.exit('--extend-steps needs --repeat-mask')
    if STEPS != sorted(set(STEPS)):
        sys.exit('--extend-steps must increase: %s' % a.extend_steps)
    CUT = True
JOINT = a.joint
if JOINT and (CLEAR_MARGIN <= 0 or not SINGLE):
    sys.exit('--joint needs --single-copy and --clear-margin')

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
sat_masked = {}      # letters masked as catalogue-family pieces (--repeat-mask)
masked_bp = {}       # letters masked as tandem repeat (TRF)
rep_masked = {}      # letters masked as pieces with another own-genome copy (--repeat-mask)
pre_rep = {}         # flank after the tandem mask, before the repeat mask (--joint uses these)
blocked = set()      # flanks a longer window cannot help: contig end, next locus
flank_len = {}       # window length used per flank


def prep(FLANK, keys):
    """Steps 2-3c for the flanks `keys` at window length FLANK. Called once (published runs) or once per
    --extend-steps length for the flanks not yet resolved. Only touches the given keys."""
    def cur():
        return [k for k in keys if k in flanks]
    for k in keys:                      # a retried flank starts clean
        for d in (flanks, flankseq, sat_masked, masked_bp, rep_masked, pre_rep):
            d.pop(k, None)
    for (i, side) in keys:
        sp, scf, s, e, fam, lid, mono = loci[i]
        flank_len[(i, side)] = FLANK
        if side == 'L':
            fs, fe = s - FLANK + 1, s
        else:
            fs, fe = e + 1, e + FLANK
        if CUT:
            fs, fe, why = FC.cut_flank(side, fs, fe, s, e, SL[(sp, scf)], lidx[(sp, scf)], i, GUARD)
            if fe - fs + 1 < MIN_ALN:
                fclass[(i, side)] = 'satellite' if why == 'locus' else 'truncated'
                blocked.add((i, side))
                continue
        if fs < 1 or fe > SL[(sp, scf)]:
            fclass[(i, side)] = 'truncated'; blocked.add((i, side)); continue
        if overlaps_other(sp, scf, fs, fe, i):
            fclass[(i, side)] = 'satellite'; continue
        mx = RAWMAX.get((i, side))
        if mx and mx[0] <= fs and fe <= mx[1]:
            seq = mx[2][fs - mx[0]:fe - mx[0] + 1]        # shorter windows are part of the longest one
        else:
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

    # residual-satellite drop (published runs) or mask (--repeat-mask)
    if REPMASK:
        sat_mask_windows(cur())
        sat_drop = set('%d_%s' % k for k in cur() if FC.unmasked(flanks[k][1]) < MIN_ALN)
    else:
        with open(OUT + '/uflanks.fa', 'w') as f:
            for (i, side), (sp, seq, g) in ((k, flanks[k]) for k in cur()):
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
    for k in cur():
        if '%d_%s' % (k[0], k[1]) in sat_drop:
            del flanks[k]; fclass[k] = 'satellite'
    print('after residual-satellite %s: %d unique anchors' % ('mask' if REPMASK else 'drop', len(flanks)), flush=True)

    # ---- 3b. tandem mask from the run's genome-wide TRF calls (flank_copy.py rule 1) ----
    if TRF_DAT:
        TIV = trf_intervals()
        for k in cur():
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

    for k in cur():
        pre_rep[k] = flanks[k][1]

    # ---- 3c. repeat mask: pieces with another copy in the own genome become N (--repeat-mask) ----
    if REPMASK:
        for k in cur():
            sp_, seq, g = flanks[k]
            iv = [(a_, min(b_, len(seq))) for a_, b_ in PIECES.get(k, ()) if a_ <= len(seq)]
            mseq, n = FC.mask(seq, 1, FC.merge(iv))
            rep_masked[k] = n
            if FC.unmasked(mseq) < MIN_ALN:
                del flanks[k]; fclass[k] = 'masked'; continue
            flanks[k] = (sp_, mseq, g)
        print('after repeat mask: %d anchors (%d flanks mostly repeated pieces)'
              % (len(flanks), sum(1 for v in fclass.values() if v == 'masked')), flush=True)


_TIV = {}


def trf_intervals():
    """TRF calls per genome and scaffold (scaffolds holding loci), parsed once; calls lying mostly inside catalogue
    arrays are the arrays themselves (or bridges between two arrays) and are not used."""
    if not _TIV:
        need = collections.defaultdict(set)
        for (sp, scf) in lidx:
            need[sp].add(scf)
        miss = sorted(sp for sp in need if sp not in TRF_DAT)
        if miss:
            sys.exit('--trf-dat given but no TRF calls for genome(s): %s' % ', '.join(miss))
        for sp, paths in TRF_DAT.items():
            acc = collections.defaultdict(list)
            for pth in paths:
                for scf, iv in FC.parse_trf_dat(pth, need.get(sp, set()), merged=False).items():
                    acc[scf].extend(iv)
            _TIV[sp] = {scf: FC.merge(FC.drop_calls_in_loci(iv, [(s2, e2) for s2, e2, _ in lidx[(sp, scf)]]))
                        for scf, iv in acc.items()}
    return _TIV


def sat_mask_windows(keys):
    """--repeat-mask: mask catalogue-family pieces (>= 60 letters at >= 75 %, the drop rule) in one search over
    overlapping windows. ssearch36 reports one alignment per query x family, so a whole long flank would need one
    search per piece. A second and third pass only for windows that still hit."""
    W = max(300, 2 * max(MONO.values()))
    step = W // 2
    todo = list(keys)
    for pass_ in range(3):
        n_w = 0
        with open(OUT + '/satwin.fa', 'w') as f:
            for k in todo:
                seq = flanks[k][1]
                for off in range(0, max(1, len(seq) - step), step):
                    w = seq[off:off + W]
                    if FC.unmasked(w) >= 60:
                        f.write('>%d_%s|%d\n%s\n' % (k[0], k[1], off, w))
                        n_w += 1
        if not n_w:
            return
        o = subprocess.run([SS, '-n', '-z', '3', '-m', '8C', '-E', '1e-2', '-b', '4000', '-d', '1', '-T', str(THREADS),
                            OUT + '/satwin.fa', OUT + '/allcons.x2.fa'], capture_output=True, text=True).stdout
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
                nm, off = p[0].rsplit('|', 1)
                i_, sd = nm.split('_')
                hit[(int(i_), sd)].append((min(qs, qe) + int(off), max(qs, qe) + int(off)))
        if not hit:
            return
        for k, iv in hit.items():
            sp, seq, g = flanks[k]
            mseq, n = FC.mask(seq, 1, FC.merge(iv))
            sat_masked[k] = sat_masked.get(k, 0) + n
            flanks[k] = (sp, mseq, g)
        todo = sorted(hit)


def window(k, L):
    """Genomic window of flank k at length L, as prep() cuts it; None if prep() would reject it on length."""
    i, side = k
    sp, scf, s, e, fam, lid, mono = loci[i]
    fs, fe = (s - L + 1, s) if side == 'L' else (e + 1, e + L)
    if CUT:
        fs, fe, why = FC.cut_flank(side, fs, fe, s, e, SL[(sp, scf)], lidx[(sp, scf)], i, GUARD)
        if fe - fs + 1 < MIN_ALN:
            return None
    if fs < 1 or fe > SL[(sp, scf)]:
        return None
    return sp, scf, fs, fe


RAWMAX = {}     # flank -> (fs, fe, raw seq) of its longest window (shorter windows are cut from it)
PIECES = {}     # flank -> pieces (1-based, oriented array->outward) with another own-genome copy
if REPMASK:
    # one minimap2 search per genome, on the longest windows (tandem letters masked): the pieces of every shorter
    # window are the same pieces cut at its length. Replaces one search per extension round; the whole-genome
    # copy test on the masked flanks (4b) still runs once at the end.
    MM = shutil.which('minimap2') or sys.exit('--repeat-mask needs minimap2 on PATH')
    TV = trf_intervals() if TRF_DAT else {}
    by_sp = collections.defaultdict(list)
    for i in range(len(loci)):
        for side in ('L', 'R'):
            w = window((i, side), STEPS[-1])
            if w:
                sp, scf, fs, fe = w
                RAWMAX[(i, side)] = (fs, fe, faidx(sp, scf, fs, fe))
                by_sp[sp].append((i, side))
    for sp in ASM:
        ks = by_sp.get(sp)
        if not ks:
            continue
        place = {}
        with open(OUT + '/repmask_%s.fa' % sp, 'w') as f:
            for k in ks:
                fs, fe, raw = RAWMAX[k]
                scf = loci[k[0]][1]
                m, _ = FC.mask(raw, fs, TV.get(sp, {}).get(scf, []))
                f.write('>%d_%s\n%s\n' % (k[0], k[1], m if k[1] == 'R' else rc(m)))
                place['%d_%s' % k] = (scf, fs, fe)
        r = subprocess.run([MM, '-k13', '-w4', '-s40', '-m25', '-A2', '-B4', '-O4,24', '-E2,1', '-c', '-N', '200', '-p', '0',
                            '-t', str(THREADS), GEN[sp], OUT + '/repmask_%s.fa' % sp], capture_output=True, text=True)
        if r.returncode:
            sys.exit('minimap2 failed on %s (exit %d): %s' % (sp, r.returncode, r.stderr.strip()[-300:]))
        for nm, iv in FC.paf_other_pieces(r.stdout.splitlines(), place, COPY_MIN_PID, MASK_MIN_ALN).items():
            i_, sd = nm.split('_')
            PIECES[(int(i_), sd)] = FC.merge(iv)
    print('repeat pieces: one search per genome on %d windows of up to %d letters' % (len(RAWMAX), STEPS[-1]),
          flush=True)

# one window length (published runs), or extend on demand (--extend-steps): a flank that is not resolved at one
# length (mostly tandem / masked / satellite pieces / simple) is tried again at the next, longer length; a flank
# stopped by the contig end or by the next locus is not (a longer window cannot help)
RETRY = {'tandem', 'masked', 'satellite', 'lowcomplex', 'telomeric'}
pending = [(i, side) for i in range(len(loci)) for side in ('L', 'R')]
for n_, L_ in enumerate(STEPS):
    if n_:
        print('extend to %d letters: %d flanks' % (L_, len(pending)), flush=True)
    prep(L_, pending)
    pending = [k for k in pending if k not in flanks and k not in blocked and fclass.get(k) in RETRY]

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
        # kept with --repeat-mask too: minimap2 reports at most -N 200 alignments per query, so on an unmasked
        # flank a piece with hundreds of copies hides the copies of the other pieces; the masked flank shows them
        # (molly 3,000 letters: 55 flanks passed without this test, 2026-11-08)
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

# ---- 5. joint tier (--joint): both flanks of a locus placed together ----
basis = {k: 'anchor' for k in pairs}
if JOINT:
    sp_of = {lid_of[i]: loci[i][0] for i in range(len(loci))}
    jk = [k for k in sorted(pre_rep) if (k[0], 'L') in pre_rep and (k[0], 'R') in pre_rep]
    with open(OUT + '/jflanks.fa', 'w') as f:
        for k in jk:
            f.write('>%d_%s\n%s\n' % (k[0], k[1], pre_rep[k]))
    o = subprocess.run([SS, '-n', '-z', '3', '-m', '8C', '-E', '1e-3', '-b', '4000', '-T', str(THREADS),
                        OUT + '/jflanks.fa', OUT + '/jflanks.fa'], capture_output=True, text=True).stdout
    cross = collections.defaultdict(lambda: collections.defaultdict(dict))   # (qi, qside) -> tlid -> tside -> bits
    same = collections.defaultdict(lambda: collections.defaultdict(dict))
    for line in o.splitlines():
        p = line.split('\t')
        if line.startswith('#') or len(p) < 12:
            continue
        try:
            pid, alen, bits = float(p[2]), int(p[3]), float(p[11])
        except ValueError:
            continue
        qi, qs_ = int(p[0].split('_')[0]), p[0].split('_')[1]
        ti, ts_ = int(p[1].split('_')[0]), p[1].split('_')[1]
        if qi == ti:
            continue
        if loci[qi][0] == loci[ti][0]:
            if alen >= COPY_MIN_ALN and pid >= COPY_MIN_PID:
                d = same[(qi, qs_)][ti]
                d[ts_] = max(d.get(ts_, 0), bits)
        elif alen >= MIN_ALN and pid >= MIN_PID:
            d = cross[(qi, qs_)][lid_of[ti]]
            d[ts_] = max(d.get(ts_, 0), bits)
    paired_to = collections.defaultdict(set)          # locus -> genomes it already has a pair in
    for x, y in pairs:
        paired_to[x].add(sp_of[y]); paired_to[y].add(sp_of[x])
    jdir = {}
    for i in sorted(set(k[0] for k in jk)):
        # the left+right combination must be unique among the loci of its own genome
        if FC.joint_scores(same.get((i, 'L'), {}), same.get((i, 'R'), {})):
            continue
        hl, hr = cross.get((i, 'L'), {}), cross.get((i, 'R'), {})
        sc = FC.joint_scores(hl, hr)
        by_g = collections.defaultdict(list)
        for t, b in sc.items():
            by_g[sp_of[t]].append((b, t))
        for g, hits in sorted(by_g.items()):
            if g in paired_to[lid_of[i]]:
                continue
            b = FC.clear_best(hits, CLEAR_MARGIN)
            if b is None or loci[i][0] in paired_to[b[1]]:
                continue
            jdir[(lid_of[i], b[1])] = b[0]
    jkeep, _ = FC.reciprocal_pairs(set(jdir))
    for key in sorted(jkeep):
        pairs.add(key)
        basis[key] = 'joint'
        detail[key] = (0.0, 0, max(jdir.get(key, 0), jdir.get(key[::-1], 0)))
    print('joint tier: %d loci with both flanks testable; %d pairs added' % (len(set(k[0] for k in jk)), len(jkeep)),
          flush=True)

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
    f.write('locusA\tlocusB\tspA\tspB\tpid\taln\tbits' + ('\tbasis' if JOINT else '') + '\n')
    for (x, y) in sorted(pairs):
        pid, alen, bits = detail[(x, y)]
        f.write('%s\t%s\t%s\t%s\t%.1f\t%d\t%.1f' % (x, y, lid2sp[x], lid2sp[y], pid, alen, bits)
                + ('\t' + basis[(x, y)] if JOINT else '') + '\n')

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
