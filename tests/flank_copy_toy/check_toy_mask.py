"""Check the refinement toy (make_toy_mask.py). Usage: check_toy_mask.py expected.json rules_dir mask600_dir joint_dir long_dir ext_dir"""
import json
import sys

exp = json.load(open(sys.argv[1]))
bad = 0
for name, d in zip(('rules', 'mask600', 'joint', 'long', 'ext'), sys.argv[2:7]):
    got = sorted(sorted(l.split('\t')[:2]) for l in open(d + '/ortho_pairs.tsv').read().splitlines()[1:])
    want = sorted(sorted(x) for x in exp[name])
    print('%-8s pairs %s' % (name, ['-'.join(x) for x in got]))
    if got != want:
        print('FAIL %s: want %s' % (name, ['-'.join(x) for x in want]))
        bad += 1
# long run: no flank reaches into another locus, no genome letter serves two flanks of one genome
fl = json.load(open(sys.argv[5] + '/flanks.json'))
ivs = []
for lid, r in fl.items():
    for side in ('LF', 'RF'):
        f = r[side]
        if f['start'] is None:
            continue
        ivs.append((r['sp'], f['start'], f['end'], lid + side))
        for oid, (s2, e2) in exp['loci'].items():
            if oid != lid and fl[oid]['sp'] == r['sp'] and f['start'] - 1 < e2 and f['end'] > s2:
                print('FAIL flank %s %s reaches into locus %s' % (lid, side, oid))
                bad += 1
ivs.sort()
for x, y in zip(ivs, ivs[1:]):
    if x[0] == y[0] and y[1] <= x[2]:
        print('FAIL flanks %s and %s share letters' % (x[3], y[3]))
        bad += 1
print(open(sys.argv[5] + '/flank_copies.tsv').read())
# extend on demand: flanks resolved at the first length keep it
ln = {'%s|%s' % (p[0], p[1]): p[7] for p in (l.split('	') for l in open(sys.argv[6] + '/flank_copies.tsv').read().splitlines()[1:])}
for k, v in exp['ext_len'].items():
    if ln.get(k) != str(v):
        print('FAIL ext flank length %s: want %s got %s' % (k, v, ln.get(k)))
        bad += 1
# the joint pair is marked as such
jb = {tuple(sorted(l.split('	')[:2])): l.split('	')[-1] for l in open(sys.argv[4] + '/ortho_pairs.tsv').read().splitlines()[1:]}
if jb.get(('Ja', 'Jb')) != 'joint':
    print('FAIL joint pair Ja-Jb not marked joint: %s' % jb.get(('Ja', 'Jb')))
    bad += 1
# flanks.json keeps the masked letters (viewer): Fa's left flank (600) = 300 unique letters, then 300 letters of the
# dispersed repeat next to the array -> the repeat mask lies in its array-side half and the unique half is used
fe = json.load(open(sys.argv[6] + '/flanks.json'))['Fa']['LF']
rm = (fe.get('mask') or {}).get('repeat', [])
n_rep = sum(b - a + 1 for a, b in rm)
if not (250 <= n_rep <= 320 and rm and min(a for a, _ in rm) > fe['start'] + 250 and fe.get('left', 0) >= 250):
    print('FAIL flanks.json mask of Fa LF: %s, left %s (flank %s-%s)' % (rm, fe.get('left'), fe['start'], fe['end']))
    bad += 1
else:
    print('Fa LF mask: repeat %s, %s letters left' % (rm, fe['left']))
if bad:
    sys.exit('TOY (mask) FAILED: %d problems' % bad)
print('PASS toy (mask): every known answer reproduced')
