# -*- coding: utf-8 -*-
"""Reconstruct chunks/chunk_44.dat from the pieces produced during the pathology hunt.
Offsets are chunk_44-local (0-based), derived from the exact split scripts used:
  split_chunk44.sh:      sub44_4 offset within chunk_44   = 1,500,000
  split_sub44_4.sh:      sub44_4_4 offset within sub44_4  =   360,000
  (manual) prefix = sub44_4_4[1-66500], offset 0 (excludes the ~57kb telomeric tail)
  (manual) half1 = prefix[1-33000] offset 0; half2 = prefix[33001-66500] offset 33,000
  (manual) q1 = half1[1-16500] offset 0; q2 = half1[16501-33000] offset 16,500
"""
import glob

BASE = '/staging/tmp/clsat_v6/scan85/chunks'

pieces = [
    (f'{BASE}/sub44/sub44_1.log', 0),
    (f'{BASE}/sub44/sub44_2.log', 500000),
    (f'{BASE}/sub44/sub44_3.log', 1000000),
    (f'{BASE}/sub44/sub44_4_split/q1.log', 1500000 + 360000 + 0),
    (f'{BASE}/sub44/sub44_4_split/q2.log', 1500000 + 360000 + 16500),
    (f'{BASE}/sub44/sub44_4_split/half2.log', 1500000 + 360000 + 33000),
]

TELOMERE_START = 1500000 + 360000 + 66500   # 1,926,500 -- excluded, documented not fabricated
TELOMERE_END = 1982912

rows = []
for path, offset in pieces:
    n = 0
    for line in open(path):
        if line.startswith('@'):
            continue
        p = line.split()
        if len(p) < 2 or not p[0].isdigit() or not p[1].isdigit():
            continue
        s, e = int(p[0]) + offset, int(p[1]) + offset
        p[0], p[1] = str(s), str(e)
        rows.append((s, e, ' '.join(p) + '\n'))
        n += 1
    print(f'{path}: {n} rows, offset {offset}')

rows.sort(key=lambda x: (x[0], x[1]))
with open(f'{BASE}/chunk_44.dat', 'w') as out:
    out.write('@chunk_44\n')
    for s, e, line in rows:
        out.write(line)

print(f'wrote chunk_44.dat: {len(rows)} rows total')
print(f'EXCLUDED telomeric region (not TRF-scanned, documented not fabricated): '
      f'chunk_44-local {TELOMERE_START}-{TELOMERE_END} '
      f'({TELOMERE_END - TELOMERE_START} bp, ~100% TTAGGG-family repeat, verified by direct '
      f'k-mer counting on 2026-07-15)')
