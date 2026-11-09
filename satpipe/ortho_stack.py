"""Ortholog stacks for the viewer (owner's request 2026-11-09).

Every group of decided ortholog pairs (a connected component of ortho_pairs.tsv) is drawn as rows, one locus per
row, all turned to read the same way. Two products per group:

  * schematic (stacks.json): flanks cut into windows of WIN alignment columns counted from the array edge, each
    with the identity of every row to every other row (the page colours a row by its identity to the reference
    row the reader picks) and the masked share of the row's letters; monomers as boxes in columns.
  * text (stack/<group>.fa, .bed, .mask.json): the aligned rows as FASTA for ViewAlign, a BED track in alignment
    columns (flanks, monomer columns, left-shifted stretches) and a block mask with each row's masked letters.

Monomers: each array is cut at its family consensus phase (position 1 of the consensus), from ssearch36 hits of
the doubled consensus on overlapping windows of the array. A monomer of a row is placed under a monomer of the
reference row (the group's first row) only when the two are each other's clear best match (second < MARGIN x
best bits); these anchors are chained collinearly. Every other monomer is left-shifted: placed in the next free
column after the previous anchor, and labelled so. Arrays with more than MAX_MONO monomers are not compared
monomer by monomer (left-shifted whole).

Tools: samtools (flank and array letters), ssearch36 (orientation, monomers), mafft (flanks, anchored columns).
"""
from __future__ import annotations

import collections
import json
import os
import subprocess
import tempfile

WIN = 100            # flank window, alignment columns
MAX_MONO = 500       # larger arrays are left-shifted whole
COMP = str.maketrans('ACGTNacgtn-', 'TGCANtgcan-')


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


# ----------------------------------------------------------------------------- groups

def groups(pairs, sp_of, sp_order):
    """Connected components of the pairs. Rows ordered by species order, then locus id (numeric part)."""
    par = {}

    def f(x):
        while par.setdefault(x, x) != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for a, b in pairs:
        par[f(a)] = f(b)
    comp = collections.defaultdict(list)
    for x in list(par):
        comp[f(x)].append(x)
    rank = {s: i for i, s in enumerate(sp_order)}
    key = lambda x: (rank.get(sp_of[x], 99), int(x[1:]) if x[1:].isdigit() else 0, x)
    out = [sorted(v, key=key) for v in comp.values()]
    out.sort(key=lambda g: key(g[0]))
    return out


def orient(group, pairs, rel):
    """Row orientation (+1 / -1) relative to the group's first row. rel[(a, b)] = +1 when b reads the same way as
    a, -1 when reversed, None when unknown; walked breadth-first over the pairs. Unknown -> +1, reported."""
    adj = collections.defaultdict(list)
    gs = set(group)
    for a, b in pairs:
        if a in gs and b in gs:
            adj[a].append(b); adj[b].append(a)
    o = {group[0]: 1}
    unknown = []
    q = [group[0]]
    while q:
        a = q.pop(0)
        for b in sorted(adj[a]):
            if b in o:
                continue
            r = rel.get((a, b))
            if r is None:
                r = rel.get((b, a))
            if r is None:
                unknown.append((a, b)); r = 1
            o[b] = o[a] * r
            q.append(b)
    for x in group:
        o.setdefault(x, 1)
    return o, unknown


def pair_relation(hits):
    """hits: {(sideA, sideB): bits} of A's flanks against B's flanks, both oriented array -> outward, forward
    strand only. Same orientation = L~L / R~R; reversed = L~R / R~L. None when nothing matched."""
    same = hits.get(('L', 'L'), 0) + hits.get(('R', 'R'), 0)
    rev = hits.get(('L', 'R'), 0) + hits.get(('R', 'L'), 0)
    if not same and not rev:
        return None
    return 1 if same >= rev else -1


# ----------------------------------------------------------------------------- monomers

def boundaries_from_hits(hits, mono, n):
    """hits: [(win_start0, qs, qe, ss, se, bits)] of the doubled consensus (1-based qs..qe) on windows of the array
    (window letters 1-based ss..se, window starting at array index win_start0). A monomer starts where the
    consensus position is 1 (mod mono). Estimates from overlapping windows are merged (within mono/3) by median.
    Returns sorted 0-based array indices in (0, n)."""
    est = []
    for w0, qs, qe, ss, se, bits in hits:
        p = ss + ((1 - qs) % mono)
        while p <= se:
            a = w0 + p - 1
            if 0 < a < n:
                est.append(a)
            p += mono
    est.sort()
    out, cl = [], []
    for a in est:
        if cl and a - cl[0] > mono / 3:
            out.append(cl[len(cl) // 2]); cl = []
        cl.append(a)
    if cl:
        out.append(cl[len(cl) // 2])
    # two merged points closer than mono/2 are one boundary seen twice
    res = []
    for a in out:
        if res and a - res[-1] < mono / 2:
            continue
        res.append(a)
    return res


def segments(n, bounds):
    """Array of n letters cut at bounds -> [(a, b)] 0-based half-open, in order."""
    pts = [0] + [b for b in bounds if 0 < b < n] + [n]
    return [(pts[i], pts[i + 1]) for i in range(len(pts) - 1) if pts[i + 1] > pts[i]]


# ----------------------------------------------------------------------------- monomer anchors

def clear_anchors(bits, margin):
    """bits: {(x, y): score} between row monomers x and reference monomers y. Anchor = reciprocal clear best
    (second < margin x best on both sides). Then the collinear chain (x and y both increasing) with the largest
    summed score. Returns [(x, y)] sorted."""
    bx = collections.defaultdict(list)
    by = collections.defaultdict(list)
    for (x, y), s in bits.items():
        bx[x].append((s, y)); by[y].append((s, x))

    def clear(lst):
        lst = sorted(lst, reverse=True)
        if not lst:
            return None
        if len(lst) > 1 and lst[1][0] >= margin * lst[0][0]:
            return None
        return lst[0][1]
    best_x = {x: clear(v) for x, v in bx.items()}
    best_y = {y: clear(v) for y, v in by.items()}
    cand = sorted((x, y) for x, y in best_x.items() if y is not None and best_y.get(y) == x)
    if not cand:
        return []
    sc = [bits[c] for c in cand]
    dp, prev = sc[:], [-1] * len(cand)
    for i in range(len(cand)):
        for j in range(i):
            if cand[j][0] < cand[i][0] and cand[j][1] < cand[i][1] and dp[j] + sc[i] > dp[i]:
                dp[i], prev[i] = dp[j] + sc[i], j
    i = max(range(len(cand)), key=lambda k: dp[k])
    chain = []
    while i >= 0:
        chain.append(cand[i]); i = prev[i]
    return chain[::-1]


def layout(n_ref, rows_n, anchors):
    """Column layout. n_ref: monomers of the reference row; rows_n: {row: n monomers} (other rows);
    anchors: {row: [(x, y)]}. Column keys (j, k): j = reference monomer index (-1 = before the first), k = 0 the
    reference column itself, k >= 1 an inserted column after it. Returns (keys sorted, place) with
    place[row] = [(key, anchored)] per monomer of that row; the reference row is under key 'ref'."""
    place = {'ref': [((j, 0), True) for j in range(n_ref)]}
    keys = set(k for k, _ in place['ref'])
    for r, n in rows_n.items():
        anc = dict(anchors.get(r, []))
        pl = [None] * n
        for x, y in anc.items():
            pl[x] = ((y, 0), True)
        x = 0
        while x < n:
            if pl[x] is not None:
                x += 1
                continue
            a = x - 1                                     # previous anchored monomer of this row
            ya = anc[a] if a >= 0 else -1
            z = x
            while z < n and pl[z] is None:
                z += 1
            yb = anc[z] if z < n else n_ref               # next anchor's reference column
            slots = [(j, 0) for j in range(ya + 1, yb)]
            last = yb - 1 if yb - 1 > ya else ya
            t = 1
            while len(slots) < z - x:
                slots.append((last, t)); t += 1
            for i in range(x, z):
                pl[i] = (slots[i - x], False)
            x = z
        place[r] = pl
        keys.update(k for k, _ in pl)
    return sorted(keys), place


# ----------------------------------------------------------------------------- flank windows

def flank_windows(aln, masks, side, win=WIN):
    """aln: {row: aligned flank string}; masks: {row: set of 0-based letter indices of that row's flank (row
    orientation) that are masked}; side 'L' (array to the right of the flank) or 'R'. Windows of `win` columns
    counted from the array edge. Returns (n_windows, per-row {masked share, letters}, identity
    {rowA|rowB: [pid or None per window]}), windows listed from the array edge outward."""
    rows = list(aln)
    L = len(next(iter(aln.values()))) if aln else 0
    cols = list(range(L - 1, -1, -1)) if side == 'L' else list(range(L))
    nwin = (L + win - 1) // win
    letter_idx = {}
    for r in rows:
        s = aln[r]
        idx, k = [None] * L, 0
        it = range(L - 1, -1, -1) if side == 'L' else range(L)
        # letters are counted from the array edge too, so mask indices (from the edge) match
        for c in it:
            if s[c] != '-':
                idx[c] = k; k += 1
        letter_idx[r] = idx
    per_row = {r: {'masked': [], 'letters': []} for r in rows}
    ident = {}
    for w in range(nwin):
        cs = cols[w * win:(w + 1) * win]
        for r in rows:
            li = [letter_idx[r][c] for c in cs if letter_idx[r][c] is not None]
            m = sum(1 for i in li if i in masks.get(r, ()))
            per_row[r]['letters'].append(len(li))
            per_row[r]['masked'].append(round(m / len(li), 2) if li else None)
        for i, a in enumerate(rows):
            for b in rows[i + 1:]:
                same = cmp = 0
                for c in cs:
                    x, y = aln[a][c], aln[b][c]
                    if x == '-' and y == '-':
                        continue
                    cmp += 1
                    same += (x == y and x != '-' and x != 'N')
                ident.setdefault(a + '|' + b, []).append(round(100.0 * same / cmp, 1) if cmp >= win // 4 else None)
    return nwin, per_row, ident


# ----------------------------------------------------------------------------- external tools

def _run(cmd, **kw):
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, **kw)
    if r.returncode:
        raise SystemExit('ortho_stack: %s failed (exit %d): %s' % (cmd[0], r.returncode, r.stderr.strip()[-400:]))
    return r.stdout


def read_fa(text):
    out, name = {}, None
    for line in text.splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]; out[name] = []
        elif name is not None:
            out[name].append(line.strip())
    return {k: ''.join(v) for k, v in out.items()}


def write_fa(path, recs):
    with open(path, 'w') as f:
        for k, s in recs:
            f.write('>%s\n%s\n' % (k, s))


def ssearch_bits(ss, query, lib, tmp, forward_only=True, nlib=None):
    """All query x library best local alignments, -m 8C rows -> {(q, t): (bits, pid, alen, qs, qe, ss, se)}."""
    qf, lf = os.path.join(tmp, 'q.fa'), os.path.join(tmp, 'l.fa')
    write_fa(qf, query); write_fa(lf, lib)
    cmd = [ss, '-n', '-z', '3', '-m', '8C', '-E', '1000', '-d', '0', '-b', str(nlib or max(1, len(lib)))]
    if forward_only:
        cmd.append('-3')
    out = {}
    for line in _run(cmd + [qf, lf]).splitlines():
        p = line.split('\t')
        if line.startswith('#') or len(p) < 12:
            continue
        k = (p[0], p[1])
        v = (float(p[11]), float(p[2]), int(p[3]), int(p[6]), int(p[7]), int(p[8]), int(p[9]))
        if k not in out or v[0] > out[k][0]:
            out[k] = v
    return out


def mafft(mf, recs, tmp):
    """recs [(name, seq)] -> aligned {name: seq} (input order kept). One sequence or all empty: padded."""
    live = [(k, s) for k, s in recs if s]
    if len(live) < 2:
        L = max((len(s) for _, s in recs), default=0)
        return {k: s.ljust(L, '-') for k, s in recs}
    f = os.path.join(tmp, 'm.fa')
    write_fa(f, live)
    a = read_fa(_run([mf, '--quiet', '--auto', '--preservecase', f]))
    L = len(next(iter(a.values())))
    return {k: a.get(k, '-' * L).upper() for k, _ in recs}


# ----------------------------------------------------------------------------- build

def cut_monomers(ss, arr, cons, mono, tmp):
    """Monomer segments of an array (row orientation) at the consensus phase. Returns ([(a, b)], strand of the
    array relative to the consensus in row orientation)."""
    n = len(arr)
    if n < mono * 1.5:
        return [(0, n)], '+'
    step, wl = mono, 2 * mono
    starts = list(range(0, max(1, n - wl + 1), step))
    if starts[-1] + wl < n:
        starts.append(max(0, n - wl))
    dc = cons + cons
    best = {}
    for strand, seq in (('+', arr), ('-', rc(arr))):
        lib = [('w%d' % w, seq[w:w + wl]) for w in starts]
        h = ssearch_bits(ss, [('c', dc)], lib, tmp, forward_only=True, nlib=len(lib))
        best[strand] = [(int(t[1:]), v[3], v[4], v[5], v[6], v[0]) for (q, t), v in h.items()
                        if v[2] >= mono // 2 and v[1] >= 60]
    strand = '+' if sum(x[5] for x in best['+']) >= sum(x[5] for x in best['-']) else '-'
    b = boundaries_from_hits(best[strand], mono, n)
    if strand == '-':
        b = sorted(n - x for x in b)
    return segments(n, b), strand


def build_stacks(loci, pairs, sp_order, fam_cons, fam_colour, genomes, out_dir, tools, margin=0.9,
                 flank_view=5000, win=WIN, max_mono=MAX_MONO, say=print):
    """loci: flanks.json dict (lid -> record with sp, scf, start (0-based), end, subtype, mono, LF/RF with class,
    mask). pairs: [(a, b)]. genomes: code -> fasta (indexed). tools: {'samtools','ssearch36','mafft'} paths.
    Writes out_dir/stacks.json and out_dir/stack/G###.{fa,bed,mask.json}. Returns the stacks dict."""
    sam, ss, mf = tools['samtools'], tools['ssearch36'], tools['mafft']
    os.makedirs(os.path.join(out_dir, 'stack'), exist_ok=True)
    sp_of = {k: v['sp'] for k, v in loci.items()}
    lens = {}
    for code, g in genomes.items():
        for line in open(g + '.fai'):
            p = line.split('\t'); lens[(code, p[0])] = int(p[1])

    def fetch(sp, scf, a, b):          # 0-based half-open, clipped to the contig
        a, b = max(0, a), min(lens[(sp, scf)], b)
        if b <= a:
            return ''
        return ''.join(read_fa(_run([sam, 'faidx', genomes[sp], '%s:%d-%d' % (scf, a + 1, b)])).values()).upper()

    def masked_genomic(rec, side):
        out = set()
        for runs in ((rec.get(side) or {}).get('mask') or {}).values():
            for a, b in runs:
                out.update(range(a, b + 1))        # genomic 1-based
        return out

    result = {'window': win, 'max_mono': max_mono, 'margin': margin, 'groups': []}
    with tempfile.TemporaryDirectory() as tmp:
        for gi, grp in enumerate(groups(pairs, sp_of, sp_order), 1):
            gid = 'G%03d' % gi
            gp = [(a, b) for a, b in pairs if a in grp and b in grp]
            # --- orientation from the flanks (oriented array -> outward, forward strand only)
            rel = {}
            for a, b in gp:
                q, lb = [], []
                for x, lst in ((a, q), (b, lb)):
                    r = loci[x]
                    for side, key in (('L', 'LF'), ('R', 'RF')):
                        s = (r.get(key) or {}).get('seq') or ''
                        if s:
                            lst.append((side, rc(s) if side == 'L' else s))
                hits = {}
                if q and lb:
                    for (qs_, ts_), v in ssearch_bits(ss, q, lb, tmp).items():
                        if v[2] >= 100 and v[1] >= 70:
                            hits[(qs_, ts_)] = max(hits.get((qs_, ts_), 0), v[0])
                rel[(a, b)] = pair_relation(hits)
            ori, unknown = orient(grp, gp, rel)
            # --- row letters in row orientation
            rows = []
            for x in grp:
                r = loci[x]
                sp, scf, s, e = r['sp'], r['scf'], r['start'], r['end']
                lf, arr, rf = fetch(sp, scf, s - flank_view, s), fetch(sp, scf, s, e), fetch(sp, scf, e, e + flank_view)
                mL, mR = masked_genomic(r, 'LF'), masked_genomic(r, 'RF')
                # masked letters as indices counted from the array edge
                lf_m = {s - g for g in mL if s - flank_view < g <= s}           # genomic 1-based g; letter s is the edge
                rf_m = {g - e - 1 for g in mR if e < g <= e + flank_view}
                if ori[x] == 1:
                    row = {'lf': lf, 'arr': arr, 'rf': rf, 'lf_m': lf_m, 'rf_m': rf_m,
                           'lcls': r['LF']['class'], 'rcls': r['RF']['class']}
                else:
                    row = {'lf': rc(rf), 'arr': rc(arr), 'rf': rc(lf), 'lf_m': rf_m, 'rf_m': lf_m,
                           'lcls': r['RF']['class'], 'rcls': r['LF']['class']}
                fam = r['subtype']
                segs, cstrand = cut_monomers(ss, row['arr'], fam_cons[fam], r['mono'], tmp)
                row.update(id=x, sp=sp, scf=scf, start=s, end=e, ori='+' if ori[x] == 1 else '-', fam=fam,
                           mono=r['mono'], segs=segs, cstrand=cstrand)
                rows.append(row)
            # --- monomer anchors against the reference (first) row
            ref = rows[0]
            anchors, skipped = {}, []
            for row in rows[1:]:
                if len(ref['segs']) > max_mono or len(row['segs']) > max_mono:
                    skipped.append(row['id']); anchors[row['id']] = []
                    continue
                q = [('%d' % i, row['arr'][a:b]) for i, (a, b) in enumerate(row['segs'])]
                lb = [('%d' % i, ref['arr'][a:b]) for i, (a, b) in enumerate(ref['segs'])]
                h = ssearch_bits(ss, q, lb, tmp, nlib=len(lb))
                bits = {(int(k[0]), int(k[1])): v[0] for k, v in h.items() if v[1] >= 70}
                anchors[row['id']] = clear_anchors(bits, margin)
            keys, place = layout(len(ref['segs']), {r['id']: len(r['segs']) for r in rows[1:]}, anchors)
            place[ref['id']] = place.pop('ref')
            col = {k: i for i, k in enumerate(keys)}
            # --- flanks aligned (row orientation; L flank ends at the array)
            alnL = mafft(mf, [(r['id'], r['lf']) for r in rows], tmp)
            alnR = mafft(mf, [(r['id'], r['rf']) for r in rows], tmp)
            nL, prL, idL = flank_windows(alnL, {r['id']: r['lf_m'] for r in rows}, 'L', win)
            nR, prR, idR = flank_windows(alnR, {r['id']: r['rf_m'] for r in rows}, 'R', win)
            # --- monomer columns: anchored columns aligned by mafft, the rest left-justified
            percol = collections.defaultdict(list)
            for r in rows:
                for i, (k, anc) in enumerate(place[r['id']]):
                    a, b = r['segs'][i]
                    percol[col[k]].append((r['id'], r['arr'][a:b], anc))
            colaln = {}
            for c in range(len(keys)):
                items = percol.get(c, [])
                anchored = sum(1 for _, _, anc in items if anc) >= 2      # the reference + at least one row
                if anchored:
                    colaln[c] = mafft(mf, [(rid, s) for rid, s, _ in items], tmp)
                else:
                    L = max((len(s) for _, s, _ in items), default=0)
                    colaln[c] = {rid: s.ljust(L, '-') for rid, s, _ in items}
            colw = [len(next(iter(colaln[c].values()))) if colaln[c] else 0 for c in range(len(keys))]
            # --- FASTA, BED, mask
            names = {r['id']: '%s_%s_%s:%d-%d%s' % (r['id'], r['sp'], r['scf'], r['start'] + 1, r['end'], r['ori'])
                     for r in rows}
            fa, maskblocks = [], []
            LL, LR = len(next(iter(alnL.values()))), len(next(iter(alnR.values())))
            arr_off = LL
            for ri, r in enumerate(rows):
                parts = [alnL[r['id']]] + [colaln[c].get(r['id'], '-' * colw[c]) for c in range(len(keys))] + [alnR[r['id']]]
                fa.append((names[r['id']], ''.join(parts)))
                # masked letters -> 0-based alignment columns
                for side, aln, off, m in (('L', alnL[r['id']], 0, r['lf_m']), ('R', alnR[r['id']], LL + sum(colw), r['rf_m'])):
                    k, cols = 0, []
                    it = range(len(aln) - 1, -1, -1) if side == 'L' else range(len(aln))
                    for c in it:
                        if aln[c] != '-':
                            if k in m:
                                cols.append(off + c)
                            k += 1
                    cols.sort()
                    for c in cols:
                        if maskblocks and maskblocks[-1]['rows'] == [ri] and maskblocks[-1]['col_end'] == c - 1:
                            maskblocks[-1]['col_end'] = c
                        else:
                            maskblocks.append({'type': 'DIVERGENT', 'col_start': c, 'col_end': c, 'rows': [ri]})
            write_fa(os.path.join(out_dir, 'stack', gid + '.fa'), fa)
            json.dump({'row_headers': [n for n, _ in fa], 'blocks': maskblocks},
                      open(os.path.join(out_dir, 'stack', gid + '.mask.json'), 'w'), separators=(',', ':'))
            bed = ['track name="%s"' % gid, 'alignment\t0\t%d\tleft flank\t0\t.\t0\t0\t148,163,184' % LL]
            c0 = LL
            shifted = None
            for c in range(len(keys)):
                items = percol.get(c, [])
                anc = sum(1 for _, _, a in items if a) >= 2
                if anc:
                    if shifted:
                        bed.append('alignment\t%d\t%d\tnot anchored (left-shifted)\t0\t.\t0\t0\t203,213,225' % shifted); shifted = None
                    fam = rows[0]['fam']
                    rgb = ','.join(str(int(fam_colour.get(fam, '#60a5fa')[i:i + 2], 16)) for i in (1, 3, 5))
                    bed.append('alignment\t%d\t%d\tm%d\t0\t.\t0\t0\t%s' % (c0, c0 + colw[c], keys[c][0] + 1, rgb))
                elif colw[c]:
                    shifted = (shifted[0], c0 + colw[c]) if shifted else (c0, c0 + colw[c])
                c0 += colw[c]
            if shifted:
                bed.append('alignment\t%d\t%d\tnot anchored (left-shifted)\t0\t.\t0\t0\t203,213,225' % shifted)
            bed.append('alignment\t%d\t%d\tright flank\t0\t.\t0\t0\t148,163,184' % (c0, c0 + LR))
            open(os.path.join(out_dir, 'stack', gid + '.bed'), 'w').write('\n'.join(bed) + '\n')
            # --- schematic record
            spc = collections.Counter(r['sp'] for r in rows)
            warn = []
            if any(v > 1 for v in spc.values()):
                warn.append('two or more loci of one genome (pairs chain through each other)')
            if unknown:
                warn.append('orientation unknown for %d pair(s); drawn as found' % len(unknown))
            if skipped:
                warn.append('array over %d monomers: not compared monomer by monomer (left-shifted)' % max_mono)
            result['groups'].append({
                'id': gid, 'warn': warn, 'ncol': len(keys),
                'col_ref': [keys[c][0] if keys[c][1] == 0 else None for c in range(len(keys))],
                'pairs': [list(p) for p in gp],
                'flank': {'L': {'n': nL, 'ident': idL}, 'R': {'n': nR, 'ident': idR}},
                'rows': [{'id': r['id'], 'sp': r['sp'], 'name': names[r['id']], 'ori': r['ori'], 'fam': r['fam'],
                          'mono': r['mono'], 'len': len(r['arr']), 'cstrand': r['cstrand'],
                          'lcls': r['lcls'], 'rcls': r['rcls'], 'lf_len': len(r['lf']), 'rf_len': len(r['rf']),
                          'L': prL[r['id']], 'R': prR[r['id']],
                          'mono_boxes': [[col[k], r['segs'][i][1] - r['segs'][i][0], 1 if anc else 0]
                                         for i, (k, anc) in enumerate(place[r['id']])]}
                         for r in rows]})
            say('  stack %s: %d rows, %d monomer columns, %d anchored%s' % (
                gid, len(rows), len(keys), sum(len(v) for v in anchors.values()), ('; ' + '; '.join(warn)) if warn else ''))
    json.dump(result, open(os.path.join(out_dir, 'stacks.json'), 'w'), separators=(',', ':'))
    return result
