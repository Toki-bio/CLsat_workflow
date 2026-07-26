# -*- coding: utf-8 -*-
"""CLEAN flank-orthology from the fixed loci, with per-pair FASTA alignments for manual audit.

Replaces the contaminated published orthology entirely (per user: "old data is AI-generated
crap"). Chain:
  1. loci = v6 CLsat1 (MIN_ALN fixed) + published CLsat2/3/4, merged per species+scaffold
  2. for each locus, take flanks OUTSIDE the array; OVERLAP GUARD: a flank may not overlap any
     other locus (this is the bug that contaminated 35% of the old 'unique' flanks)
  3. classify each flank: satellite (matches a CLsat consensus) / low-complexity/telomeric /
     unique. Only UNIQUE flanks are anchors.
  4. cross-assembly orthology: ssearch36 (NOT blastn), single best cross-assembly hit per
     target assembly per query, >=200 bp aligned, >=80% identity
  5. DELIVERABLE: for every orthologous pair, write a FASTA alignment of the two flanks
     (array->outward; UPPERCASE = residual satellite, lowercase = unique) so the shared unique
     sequence can be eyeballed. Plus the orthology matrix.

Outputs in /staging/tmp/clsat_v6/ortho/:
  loci.bed  flanks.fa  flank_class.tsv  ortho_pairs.tsv  ortho_matrix.txt
  ortho_pairs_alignments.fasta   <-- the manual-inspection file
"""
import subprocess, os, collections, tempfile, random

BIN = '/staging/conda/envs/bioinfo/bin'
W = '/staging/tmp/mash_parentcheck'
V6 = '/staging/tmp/clsat_v6/full'
SD = ('/staging/Sergei/Lacertidae/Darevskia/Satellites/motifs_output/clean/'
      'cross_clusters/selected/subclusters/seeds')
PUB = '/staging/Sergei/Lacertidae/Darevskia/Satellites/manuscript/CLsat_full_clusters_6genome'
OUT = '/staging/tmp/clsat_v6/ortho'
os.makedirs(OUT, exist_ok=True)
SS = BIN + '/ssearch36'
SAM = BIN + '/samtools'

FLANK = 600           # bp of flank to test/align
MIN_ALN = 200
MIN_PID = 80.0
GUARD = 0             # a flank may not come within GUARD bp of another locus

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

# ---- 1. build the loci ----
raw = []
for line in open(V6 + '/v6full_CLsat1.expanded_merged.bed'):
    p = line.rstrip('\n').split('\t')
    if len(p) >= 4:
        raw.append((p[3].split('|')[0], p[0], int(p[1]), int(p[2]), 'CLsat1'))
for fam in ('CLsat2', 'CLsat3', 'CLsat4'):
    fn = SD + '/out_%s_v4.expanded_merged.bed' % fam
    for line in open(fn):
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
    for (s2, e2, j) in lidx[(sp, scf)]:
        if j != self_i and fs - GUARD < e2 and fe + GUARD > s2:
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
        # also drop if it matches a CLsat consensus (residual satellite not caught by coords)
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
for (a, b) in pairs:
    M[tuple(sorted((loci[int(a.split('_')[0])][0], loci[int(b.split('_')[0])][0])))] += 1
ASM = ['mix', 'arm', 'val', 'unp', 'unm', 'nai']
with open(OUT + '/ortho_matrix.txt', 'w') as f:
    hdr = '      ' + ''.join('%6s' % x for x in ASM)
    f.write(hdr + '\n'); print('\n=== CLEAN ORTHOLOGY MATRIX ==='); print(hdr)
    for a in ASM:
        row = '%5s ' % a + ''.join('     -' if a == b else '%6d' % M[tuple(sorted((a, b)))] for b in ASM)
        f.write(row + '\n'); print(row)
    tot = 'total unique undirected pairs: %d' % len(pairs)
    f.write('\n' + tot + '\n'); print('\n' + tot)
    key = 'KEY: unp-nai=%d  unp-val=%d  unm-nai=%d  unm-val=%d' % (
        M[('nai', 'unp')], M[('unp', 'val')], M[('nai', 'unm')], M[('unm', 'val')])
    f.write(key + '\n'); print(key)

# ---- 5. per-pair FASTA alignments for manual audit ----
def sw_align(a, b):
    """gapped Needleman-ish local via ssearch; return the two aligned strings for display."""
    fa = tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False); fa.write('>a\n%s\n' % a); fa.close()
    fb = tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False); fb.write('>b\n%s\n' % b); fb.close()
    o = subprocess.run([BIN + '/needle', '-asequence', fa.name, '-bsequence', fb.name,
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
    return (seqs[0], seqs[1]) if len(seqs) >= 2 else (a, b)

random.seed(1)
pl = sorted(pairs)
# prioritise the cross-haplotype pairs the paper cares about
def pr(k):
    sps = tuple(sorted((loci[int(k[0].split('_')[0])][0], loci[int(k[1].split('_')[0])][0])))
    return 0 if sps in (('nai', 'unp'), ('unp', 'val'), ('nai', 'unm'), ('unm', 'val')) else 1
pl.sort(key=pr)
with open(OUT + '/ortho_pairs_alignments.fasta', 'w') as f:
    f.write('; CLEAN flank-orthology pairs for manual validation.\n')
    f.write('; Each block: two orthologous UNIQUE flanks, aligned (needle). UPPERCASE input was\n')
    f.write('; the flank running array->outward; these are the unique-flank anchors (satellite\n')
    f.write('; and overlapping flanks already excluded). If the two sequences share long unique\n')
    f.write('; stretches, the orthology is real; if they only match in low-complexity, it is not.\n\n')
    shown = 0
    for k in pl:
        if shown >= 120:
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
        for a in range(0, len(A), 100):
            f.write('  %s  %s\n' % (sa[0][:3], A[a:a+100]))
            f.write('       %s\n' % mstr[a:a+100])
            f.write('  %s  %s\n\n' % (sb[0][:3], B[a:a+100]))
        shown += 1
print('\nwrote %s/ortho_pairs_alignments.fasta (%d pairs)' % (OUT, shown))
print('wrote %s/ortho_matrix.txt, loci.bed, flank_class.tsv-equivalent' % OUT)
with open(OUT + '/flank_class.tsv', 'w') as f:
    f.write('locus\tside\tclass\n')
    for (i, side), kl in fclass.items():
        f.write('R%d\t%s\t%s\n' % (i + 1, side, kl))
