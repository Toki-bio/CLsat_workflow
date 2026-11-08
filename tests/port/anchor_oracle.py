"""WP9 differential oracle for anchor.awk: generate cases (gen) and evaluate them with satpipe/anchor.py (run).

  python tests/port/anchor_oracle.py gen DIR N SEED     -> DIR/cases.tsv (+ DIR/fa/*.fa for read_fasta cases)
  python tests/port/anchor_oracle.py run DIR/cases.tsv  -> one result line per case on stdout (tests/port/anchor_diff.awk
                                                           prints the same format from anchor.awk)
Case rows (TAB-separated; records 'name=seq,name=seq'):
  A  id q t                                         local_align
  F  id anchor motif max_mm                         anchor_frame
  R  id seq anchor mi min_aligned motif max_mm max_snap   rotate_to_anchor
  C  id seq motif max_mm                            check_frame
  P  id anchor_seq motif max_mm mi records          anchor_frame, then anchor_records and verify_records (as stage 21 / cli)
  RF id path                                        read_fasta
"""
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from satpipe import anchor as A  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
CONS = os.path.join(ROOT, 'satpipe', 'stages', '02_seed_expansion_scan', 'consensus')


def esc(s):
    return s.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t').replace('\r', '\\r')


def err(e):
    t = type(e).__name__
    if t in ('AnchorError', 'ValueError'):
        return 'ERR\t%s\t%s' % (t, esc(str(e)))
    if t in ('FileNotFoundError', 'IsADirectoryError', 'PermissionError'):
        t = 'OSError'
    return 'ERR\t%s\t' % t


def recs_str(d):
    return esc(','.join('%s=%s' % kv for kv in d.items()))


def parse_recs(f):
    out = {}
    for x in f.split(',') if f else []:
        n, s = x.split('=', 1)
        out[n] = s
    return out


# ------------------------------------------------------------------ generator
def cons():
    return {f: ''.join(A.read_fasta(os.path.join(CONS, f + '.cons.fa')).values()) for f in ('CLsat1', 'CLsat2', 'CLsat3', 'CLsat4')}


def rc(s):
    return s[::-1].translate(str.maketrans('ACGTacgt', 'TGCAtgca'))


def rnd(r, n, alpha='ACGT'):
    return ''.join(r.choice(alpha) for _ in range(n))


def mutate(r, s, rate):
    out = []
    for c in s:
        x = r.random()
        if x < rate * 0.6:
            out.append(r.choice('ACGT'))
        elif x < rate * 0.8:
            continue  # deletion
        elif x < rate:
            out.append(c + rnd(r, r.randint(1, 3)))
        else:
            out.append(c)
    return ''.join(out)


def decorate(r, s):
    """lower case, U for T, N runs, gaps and spaces: everything clean() and the comparisons see"""
    k = r.random()
    if k < 0.15:
        s = s.lower()
    elif k < 0.25:
        s = s.replace('T', 'U')
    elif k < 0.35 and len(s) > 10:
        p = r.randrange(len(s) - 5)
        s = s[:p] + 'N' * r.randint(1, 8) + s[p + 5:]
    elif k < 0.45 and len(s) > 10:
        p = r.randrange(len(s))
        s = s[:p] + r.choice(['-', '--', ' ', '- -']) + s[p:]
    elif k < 0.5:
        s = ''.join(c.lower() if r.random() < 0.3 else c for c in s)
    return s


def variant(r, base, C, good=False):
    k = r.random() if not good else 0.45 + 0.55 * r.random()
    s = base
    if k < 0.08:
        return ''                                                  # empty
    if k < 0.15:
        return rnd(r, r.randint(1, 12))                            # tiny
    if k < 0.22:
        return rnd(r, r.randint(40, 300))                          # junk
    if k < 0.27:
        return r.choice(['A', 'AC', 'AAAGAAG', 'ACG']) * r.randint(5, 80)   # periodic: alignment ties
    if k < 0.35:
        s = rc(s)                                                  # reverse complement
    if k < 0.40:
        s = s + r.choice(list(C.values()))                         # dimer of two families
    s = mutate(r, s, r.choice([0, 0, 0.02, 0.05, 0.1, 0.2, 0.35] if not good else [0, 0, 0.02, 0.05, 0.1]))
    if s and r.random() < (0.15 if not good else 0.03):
        a, b = sorted(r.randrange(len(s) + 1) for _ in range(2))
        s = s[a:b] if b - a > 5 else s                             # fragment
    if s:
        k2 = r.randrange(len(s))
        s = s[k2:] + s[:k2]                                        # any rotation
    return decorate(r, s)


def gen(d, n, seed):
    r = random.Random(seed)
    C = cons()
    os.makedirs(os.path.join(d, 'fa'), exist_ok=True)
    rows = []
    kinds = ['P'] * 45 + ['R'] * 18 + ['A'] * 12 + ['F'] * 10 + ['C'] * 5 + ['RF'] * 10
    motifs = ['AAAGAAG', 'AAAGAAG', 'aaagaag', '', 'CCCCCCCCCC', 'AAAG', 'A', 'GAAGAAAGAAGAAAG']
    for i in range(n):
        k = r.choice(kinds)
        cid = 'c%04d' % i
        fam = r.choice(list(C))
        if k == 'P':
            ak = r.random()
            anchor = C['CLsat1'] if ak < 0.5 else C[fam] if ak < 0.8 else variant(r, C[fam], C)
            mk = r.random()
            m = ('AAAGAAG' if mk < 0.4 else 'aaagaag' if mk < 0.5 else '' if mk < 0.7 else
                 r.choice(motifs + [anchor[j:j + r.randint(3, 10)] for j in [r.randrange(max(1, len(anchor)))]] + [rnd(r, 5)]))
            nrec = r.randint(1, 4)
            if r.random() < 0.12:                                  # tiny anchor and records shorter than the motif
                anchor, m = 'AAAGAAG' + rnd(r, r.randint(0, 3)), 'AAAGAAG'
                recs = ['r%d=%s' % (j, ('AAAGAAG' * 2)[r.randrange(7):][:r.randint(2, 8)] if r.random() < 0.5
                                    else rnd(r, r.randint(2, 6), 'AG')) for j in range(nrec)]
                rows.append([k, cid, anchor, m, str(r.choice([0, 1, 2, 3])), r.choice(['0.0', '0.3', '0.6']), ','.join(recs)])
                continue
            good = r.random() < 0.6                                # all records placeable-looking, so the success path is exercised
            recs = ['r%d=%s' % (j, variant(r, r.choice(list(C.values())), C, good)) for j in range(nrec)]
            rows.append([k, cid, anchor, m, str(r.choice([0, 1, 1, 1, 2])), r.choice(['0.6', '0.6', '0.0', '0.5', '0.9', '1.0']), ','.join(recs)])
        elif k == 'R':
            anchor = r.choice([A.anchor_frame(C['CLsat1'], 'AAAGAAG'), C[fam], rnd(r, r.randint(0, 60)), 'AAAGAAG' * r.randint(1, 30)])
            seq = variant(r, r.choice(list(C.values())), C)
            if r.random() < 0.2:                                   # motif copies at equal distance around the start: snap ties
                core = rnd(r, r.randint(20, 80))
                seq = core[:10] + 'AAAGAAG' + core[10:] + 'AAAGAAC' + core[:5]
            rows.append([k, cid, seq, anchor, r.choice(['0.6', '0.0', '0.3', '0.9']), str(r.choice([40, 40, 0, 5, 200])),
                         r.choice(motifs + [rnd(r, 3)]), str(r.choice([0, 1, 2, 3])), str(r.choice([40, 40, 0, 3, 1000]))])
        elif k == 'A':
            al = 'ACGT' if r.random() < 0.5 else r.choice(['AC', 'A', 'ACN', 'AT'])
            q = rnd(r, r.randint(0, 40), al)
            t = rnd(r, r.randint(0, 60), al) if r.random() < 0.6 else (q[3:] + q[:3]) * 2
            rows.append([k, cid, q, t])
        elif k == 'F':
            fk = r.random()
            if fk < 0.3:
                a = C[fam]
            elif fk < 0.5:
                a = rnd(r, 30) + 'AAAGAAG' + rnd(r, 30) + 'AAAGAAG' + rnd(r, 10)          # ambiguous
            elif fk < 0.6:
                a = ''
            elif fk < 0.7:
                a = rnd(r, r.randint(1, 6))                                             # motif longer than anchor
            else:
                a = variant(r, C[fam], C)
            rows.append([k, cid, a, r.choice(motifs + [rnd(r, 4)]), str(r.choice([0, 1, 2, 3]))])
        elif k == 'C':
            rows.append([k, cid, variant(r, C[fam], C), r.choice(motifs), str(r.choice([0, 1, 2]))])
        else:
            p = os.path.join(d, 'fa', cid + '.fa')
            nl = r.choice(['\n', '\r\n', '\r', '\n'])
            lines = []
            if r.random() < 0.2:
                lines.append('stray line before any header')
            names = []
            for j in range(r.randint(0, 4)):
                nm = r.choice(['s%d' % j, 's%d' % j, 's0', 'CLsat1.%d' % j])
                names.append(nm)
                hk = r.random()
                head = '>' + (nm if hk < 0.5 else ' ' + nm + ' desc text' if hk < 0.7 else nm + '\tx\x0by' if hk < 0.8
                              else '' if hk < 0.85 else '\x0c' + nm)
                lines.append(head)
                for _ in range(r.randint(0, 3)):
                    lines.append(r.choice(['', ' ', '\t']) + decorate(r, rnd(r, r.randint(0, 30))) + r.choice(['', ' ', '\t ', '\x0b']))
            data = nl.join(lines) + (nl if r.random() < 0.8 else '')
            with open(p, 'w', encoding='utf-8', newline='') as fh:
                fh.write(data)
            rows.append([k, cid, p.replace('\\', '/')])
    with open(os.path.join(d, 'cases.tsv'), 'w', encoding='utf-8', newline='\n') as fh:
        for x in rows:
            assert not any('\t' in f or '\n' in f for f in x), x
            fh.write('\t'.join(x) + '\n')


# ------------------------------------------------------------------ evaluator
def run_case(f):
    k, cid = f[0], f[1]
    out = [k, cid]
    try:
        if k == 'A':
            out.append(' '.join(str(v) for v in A.local_align(f[2], f[3])))
        elif k == 'F':
            out.append('OK\t' + A.anchor_frame(f[2], f[3], int(f[4])))
        elif k == 'R':
            rot, shift, ident, cols, snap = A.rotate_to_anchor(f[2], f[3], min_identity=float(f[4]), min_aligned=int(f[5]),
                                                               motif=f[6], max_mm=int(f[7]), max_snap=int(f[8]))
            out.append('OK\t%s\t%d\t%.17g\t%d\t%d' % (rot, shift, ident, cols, snap))
        elif k == 'C':
            out.append('1' if A.check_frame(f[2], f[3], int(f[4])) else '0')
        elif k == 'RF':
            out.append('OK\t' + recs_str(A.read_fasta(f[2])))
        elif k == 'P':
            mm, mi, recs = int(f[4]), float(f[5]), parse_recs(f[6])
            frame = A.anchor_frame(f[2], f[3], mm)
            out.append('FRAME\t' + frame)
            out.append('VERIFY\t' + esc('|'.join(A.verify_records(recs, frame, f[3], mm, mi))))
            try:
                rot, rows = A.anchor_records(recs, frame, f[3], mm, mi)
                out.append('OUT\t' + recs_str(rot))
                out.append('ROWS\t' + esc('|'.join('\t'.join(str(v) for v in x) for x in rows)))
                out.append('REVERIFY\t' + esc('|'.join(A.verify_records(rot, frame, f[3], mm, mi))))
            except Exception as e:  # noqa: BLE001
                out.append(err(e))
    except Exception as e:  # noqa: BLE001
        out.append(err(e))
    return '\t'.join(out)


def run(path):
    with open(path, encoding='utf-8', newline='\n') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            sys.stdout.write(run_case(f) + '\n')


if __name__ == '__main__':
    if sys.argv[1] == 'gen':
        gen(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
    else:
        sys.stdout.reconfigure(newline='\n')
        run(sys.argv[2])
