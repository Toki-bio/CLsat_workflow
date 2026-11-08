"""Check the refinement toy (make_toy_mask.py). Usage: check_toy_mask.py expected.json rules_dir mask600_dir long_dir"""
import json
import sys

exp = json.load(open(sys.argv[1]))
bad = 0
for name, d in zip(('rules', 'mask600', 'long'), sys.argv[2:5]):
    got = sorted(sorted(l.split('\t')[:2]) for l in open(d + '/ortho_pairs.tsv').read().splitlines()[1:])
    want = sorted(sorted(x) for x in exp[name])
    print('%-8s pairs %s' % (name, ['-'.join(x) for x in got]))
    if got != want:
        print('FAIL %s: want %s' % (name, ['-'.join(x) for x in want]))
        bad += 1
# long run: no flank reaches into another locus, no genome letter serves two flanks of one genome
fl = json.load(open(sys.argv[4] + '/flanks.json'))
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
print(open(sys.argv[4] + '/flank_copies.tsv').read())
if bad:
    sys.exit('TOY (mask) FAILED: %d problems' % bad)
print('PASS toy (mask): every known answer reproduced')
