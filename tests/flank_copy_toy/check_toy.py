"""Compare toy outputs with the known answers (make_toy.py). Usage: check_toy.py expected.json on_dir off_dir"""
import json
import sys

exp, on, off = json.load(open(sys.argv[1])), sys.argv[2], sys.argv[3]


def pairs(d):
    return sorted(sorted(l.split('\t')[:2]) for l in open(d + '/ortho_pairs.tsv').read().splitlines()[1:])


def classes(d):
    return {f'{p[0]}|{p[1]}': p[2] for p in (l.split('\t') for l in open(d + '/flank_class.tsv').read().splitlines()[1:])}


bad = 0
got = pairs(on)
want = sorted(sorted(x) for x in exp['on']['pairs'])
print('rules on : pairs', got)
if got != want:
    print('FAIL pairs: want', want); bad += 1
c = classes(on)
for k, v in exp['on']['class'].items():
    if c.get(k) != v:
        print('FAIL class %s: want %s got %s' % (k, v, c.get(k))); bad += 1
got_off = pairs(off)
print('rules off: pairs', got_off)
for x in exp['off_must_include']:
    if sorted(x) not in got_off:
        print('FAIL rules off lost pair', x); bad += 1
if any('C' in a and 'C' in b for a, b in got_off):
    print('note: rules off pairs a C locus (tie decided by input order), as expected of the old rule')
print('flank_copies.tsv:')
print(open(on + '/flank_copies.tsv').read())
if bad:
    sys.exit('TOY FAILED: %d problems' % bad)
print('PASS toy: every known answer reproduced')
