"""Test-only oracle: write tests/port/fixtures/pyjson_strings.tsv (random strings, Python's json.dumps and repr).

Line: flat_esc(s) TAB json.dumps(s) TAB repr(s) (repr column '-' when s holds a code point whose printability the
gawk side only approximates). Deterministic (seed 20261108). Run: python tests/port/gen_pyjson_cases.py [N]
"""
import json
import os
import random
import sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 400
rng = random.Random(20261108)

ASCII = [chr(c) for c in range(0x20, 0x7f)]
SPECIAL = ['"', '\\', '/', "'", '\t', '\n', '\r', '\b', '\f', '\x00', '\x01', '\x1f', '\x7f', ' ']
LATIN1 = [chr(c) for c in range(0x80, 0x100)]
BMP = ['\u00e9', '\u0416', '\u4e2d', '\u2028', '\u2029', '\ufeff', '\u00ad', '\u200b', '\u3000', '\uffff', '\u0100', '\u07ff', '\u0800']
ASTRAL = ['\U0001f600', '\U00010000', '\U0010fffd', '\U0001d11e']
# code points where the gawk repr follows Python exactly (ASCII, Latin-1 and the listed BMP/astral ones)
REPR_EXACT = set(ASCII + SPECIAL + LATIN1 + BMP + ASTRAL) - {'\uffff'}


def flat_esc(s):
    return s.replace('\\', '\\\\').replace('\t', '\\t').replace('\n', '\\n').replace('\r', '\\r')


def rand_string():
    kind = rng.random()
    n = rng.randint(0, 40)
    out = []
    for _ in range(n):
        r = rng.random()
        if r < 0.45:
            out.append(rng.choice(ASCII))
        elif r < 0.65:
            out.append(rng.choice(SPECIAL))
        elif r < 0.75:
            out.append(rng.choice(LATIN1))
        elif r < 0.85:
            out.append(rng.choice(BMP))
        elif r < 0.92:
            out.append(rng.choice(ASTRAL))
        else:
            c = rng.randint(0x100, 0x10ffff)
            if 0xd800 <= c <= 0xdfff:
                c = 0x1f601
            out.append(chr(c))
    if kind < 0.1:   # pure quote mixes, to exercise repr's quote choice
        out = [rng.choice(["'", '"', 'a', '\\']) for _ in range(rng.randint(1, 8))]
    return ''.join(out)


cases = ['', 'a', "'", '"', '\'"', 'it\'s', 'say "hi"', '\\', '/', 'C:\\work\\x', '\x7f', '\u00a0', '\U0001f600',
         'runs-v3', 'bin-majority-v1', '#1f77b4']
while len(cases) < N:
    cases.append(rand_string())

here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, 'fixtures', 'pyjson_strings.tsv'), 'w', encoding='utf-8', newline='\n') as fh:
    for s in cases:
        rep = repr(s) if all(c in REPR_EXACT for c in s) else '-'
        fh.write('%s\t%s\t%s\n' % (flat_esc(s), json.dumps(s), flat_esc(rep)))
print('wrote', len(cases), 'cases')
