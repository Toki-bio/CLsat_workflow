# -*- coding: utf-8 -*-
"""Re-run the flank-orthology pipeline on the FINAL 627-segment catalogue (L-ids), so that
flank orthology, sub-subfamily assignment, Table 1, and the viewer are all keyed to one
consistent catalogue/id scheme. Core method matches _clean_ortho_pairs.py (600 bp flanks,
overlap guard, residual-satellite drop, ssearch36 >=80% id / >=200 bp, single best hit per
target assembly, no side restriction). It is not byte-for-byte identical: the loci source is
the segmented catalogue rather than rebuilt raw-merged loci, pairs are keyed by locus id
rather than by flank id, the ssearch36 -b cutoff differs (4000 here), and the output format
differs. The resulting 53-pair / 18-both-flank-group matrix is identical either way
(catalogue-invariant; verified against the raw/merged 539-locus run)."""
import subprocess, os, collections, tempfile, json

BIN = '/staging/conda/envs/bioinfo/bin'
SD = ('/staging/Sergei/Lacertidae/Darevskia/Satellites/motifs_output/clean/'
      'cross_clusters/selected/subclusters/seeds')
CAT = '/staging/tmp/clsat_v6/viewer_rebuild/segmented_catalog.tsv'
OUT = '/staging/tmp/clsat_v6/viewer_rebuild/ortho627_noside'
os.makedirs(OUT, exist_ok=True)
SS = BIN + '/ssearch36'
SAM = BIN + '/samtools'

FLANK = 600
MIN_ALN = 200
MIN_PID = 80.0
GUARD = 0

GEN = {}
for line in open(SD + '/genomes.tsv'):
    p = line.rstrip('\n').split('\t')
    if len(p) >= 2:
        GEN[p[0]] = p[1]
CONS = {f: ''.join(l.strip() for l in open(SD + '/%s.cons.fa' % f) if not l.startswith('>')).upper()
        for f in ('CLsat1', 'CLsat2', 'CLsat3', 'CLsat4')}
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
MONO = {'CLsat1': 146, 'CLsat2': 145, 'CLsat3': 146, 'CLsat4': 116}
for line in open(CAT):
    p = line.rstrip('\n').split('\t')
    lid, scf, s, e, sp, fam = p[0], p[1], int(p[2]), int(p[3]), p[4], p[5]
    loci.append((sp, scf, s, e, fam, lid, MONO.get(fam.split(',')[0], 146)))
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
    return ''.join(l.strip() for l in r.stdout.splitlines() if not l.startswith('>')).upper() if not r.returncode else ''

def overlaps_other(sp, scf, fs, fe, self_i):
    # flank [fs, fe] is 1-based inclusive; loci [s2, e2) are 0-based half-open.
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

# ---- 4. cross-assembly orthology at three thresholds (for viewer 70/80/90) ----
with open(OUT + '/uflanks.fa', 'w') as f:
    for (i, side), (sp, seq, g) in flanks.items():
        f.write('>%d_%s\n%s\n' % (i, side, seq))
o = subprocess.run([SS, '-n', '-m', '8C', '-E', '1e-3', '-b', '4000', OUT + '/uflanks.fa',
                    OUT + '/uflanks.fa'], capture_output=True, text=True).stdout

# collect all cross-assembly hits with alen>=200 (threshold on pid applied later per-level)
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
    # NOTE: no side restriction -- matches the original _clean_ortho_pairs.py exactly.
    # Flanks are oriented array->outward, so an L flank legitimately matches an R flank
    # when two orthologous arrays are in opposite orientations (inversion).
    allhits.append((q, t, pid, alen, bits, side))

# 80% pairs (the manuscript/Fig-3 set), single best per target assembly per query
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
for (a, b) in pairs:
    M[tuple(sorted((lid2sp[a], lid2sp[b])))] += 1
ASM = ['mix', 'arm', 'val', 'unp', 'unm', 'nai']
with open(OUT + '/ortho_matrix.txt', 'w') as f:
    hdr = '      ' + ''.join('%6s' % x for x in ASM)
    f.write(hdr + '\n'); print('\n=== ORTHOLOGY MATRIX (627 segmented catalogue) ==='); print(hdr)
    for a in ASM:
        row = '%5s ' % a + ''.join('     -' if a == b else '%6d' % M[tuple(sorted((a, b)))] for b in ASM)
        f.write(row + '\n'); print(row)
    tot = 'total unique undirected pairs: %d' % len(pairs)
    f.write('\n' + tot + '\n'); print('\n' + tot)
    key = 'KEY: unp-nai=%d  unp-val=%d  unm-nai=%d  unm-val=%d' % (
        M[('nai', 'unp')], M[('unp', 'val')], M[('nai', 'unm')], M[('unm', 'val')])
    f.write(key + '\n'); print(key)

# ---- dump machine-readable artifacts for the viewer builder ----
# flank_class.tsv (L-id, side, class)
with open(OUT + '/flank_class.tsv', 'w') as f:
    f.write('locus\tside\tclass\n')
    for (i, side), kl in fclass.items():
        f.write('%s\t%s\t%s\n' % (lid_of[i], side, kl))
# per-locus flank sequences + truncation (for flanks.json)
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
json.dump(flankdump, open(OUT + '/flanks_627.json', 'w'))

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
json.dump(ortho_json, open(OUT + '/flank_orthologs_627.json', 'w'))

# save the 80% undirected pairs list
with open(OUT + '/ortho_pairs.tsv', 'w') as f:
    f.write('locusA\tlocusB\tspA\tspB\tpid\taln\tbits\n')
    for (a, b) in sorted(pairs):
        pid, alen, bits = detail[(a, b)]
        f.write('%s\t%s\t%s\t%s\t%.1f\t%d\t%.1f\n' % (a, b, lid2sp[a], lid2sp[b], pid, alen, bits))

print('\nwrote flanks_627.json, flank_orthologs_627.json, flank_class.tsv, ortho_pairs.tsv, ortho_matrix.txt')
print('DONE_ORTHO627')
