"""Check satpipe.ortho_stack on the toy (make_toy.py). Usage: check_toy.py <toy_dir> <out_dir>"""
import json
import sys

T, O = sys.argv[1], sys.argv[2]
exp = json.load(open(T + '/expected.json'))
st = {g['id']: g for g in json.load(open(O + '/stacks.json'))['groups']}
bad = []


def ok(cond, msg):
    print(('PASS ' if cond else 'FAIL ') + msg)
    if not cond:
        bad.append(msg)


g, e = st['G001'], exp['G001']
ok([r['id'] for r in g['rows']] == e['rows'], 'G001 rows %s' % [r['id'] for r in g['rows']])
ok([r['ori'] for r in g['rows']] == e['ori'], 'G001 orientation %s (C1 reversed)' % [r['ori'] for r in g['rows']])
ref, rb, rc_ = g['rows']
ok(len(ref['mono_boxes']) == e['n_ref'], 'ref monomers %d' % len(ref['mono_boxes']))
lens = [b[1] for b in ref['mono_boxes']]
bounds = [sum(lens[:i + 1]) for i in range(len(lens) - 1)]
ok(bounds == e['ref_bounds'], 'ref boundaries at consensus phase %s' % bounds)
refcol = {b[0]: i for i, b in enumerate(ref['mono_boxes'])}
bcols = [refcol.get(b[0]) for b in rb['mono_boxes']]
ok(bcols == e['b_cols'] and all(b[2] for b in rb['mono_boxes']), 'B1 monomers under ref %s, all anchored (v6 v7 gap)' % bcols)
ccols = [refcol.get(b[0]) for b in rc_['mono_boxes']]
ok(ccols == e['c_cols'] and all(b[2] for b in rc_['mono_boxes']), 'C1 (reversed) monomers under ref %s' % ccols)
idL = g['flank']['L']['ident']['L1|L2']
ok(idL[0] is not None and idL[0] >= 95, 'A1|B1 left flank edge window identity %s' % idL[0])
ok(ref['L']['masked'][e['mask_window']] == 1.0 and ref['L']['masked'][0] == 0, 'A1 masked window %s' % ref['L']['masked'][:4])
idC = g['flank']['L']['ident']['L1|L3']
ok(sum(1 for x in idC if x is not None and x >= 95) >= 45, 'A1|C1 left flank: %d of %d windows >= 95 %% despite the insertion'
   % (sum(1 for x in idC if x is not None and x >= 95), len(idC)))
g2, e2 = st['G002'], exp['G002']
ok([r['id'] for r in g2['rows']] == e2['rows'], 'G002 rows')
b2 = g2['rows'][1]['mono_boxes']
ok(len(b2) == e2['b_n'] and sum(b[2] for b in b2) == e2['anchored_b'], 'G002 B monomers %d, anchored %d (left-shifted)'
   % (len(b2), sum(b[2] for b in b2)))
ok([b[0] for b in b2] == list(range(len(b2))), 'G002 left-shifted from the left edge %s' % [b[0] for b in b2])
g3, e3 = st['G003'], exp['G003']
b3 = g3['rows'][1]['mono_boxes']
ok([b[2] for b in b3] == e3['b_anchored'], 'G003 near-identical monomers: paired %s (3rd and 7th have two exact partners)' % [b[2] for b in b3])
ok(all(b[4] == 0 and b[5] is not None and b[5] >= 2 for b in b3 if b[2]), 'G003 paired monomers: 0 letters differ, closest other >= 2 %s'
   % [(b[4], b[5]) for b in b3])
fa = open(O + '/stack/G001.fa').read().split('>')[1:]
L = set(len(''.join(x.split('\n')[1:])) for x in fa)
ok(len(L) == 1, 'G001.fa rows of equal length %s' % L)
seqs = {x.split('\n')[0]: ''.join(x.split('\n')[1:]) for x in fa}
ok(all(len(s.replace('-', '')) > 0 for s in seqs.values()), 'G001.fa headers %s' % list(seqs))
m = json.load(open(O + '/stack/G001.mask.json'))
ok(any(b['rows'] == [0] and b['col_end'] - b['col_start'] + 1 >= 100 for b in m['blocks']), 'G001 mask: row 0 masked block of 100 letters')
bed = open(O + '/stack/G001.bed').read()
ok('left flank' in bed and 'right flank' in bed and '	m1	' in bed and '	m10	' in bed, 'G001 BED: flanks and monomer columns m1..m10')
ok('not anchored' in open(O + '/stack/G002.bed').read(), 'G002 BED: not anchored stretch')
print('TOY %s' % ('FAILED: ' + '; '.join(bad) if bad else 'ALL PASS'))
