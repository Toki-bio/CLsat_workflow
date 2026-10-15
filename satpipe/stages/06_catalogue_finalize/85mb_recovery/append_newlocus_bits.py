# -*- coding: utf-8 -*-
"""Append the new R540 locus's bits rows (all 4 families) to the existing combined
scan_CLsat{1-4}.bits.tsv (which already cover the 539-locus catalog), producing complete
bits files covering all 540 loci."""

V = '/staging/tmp/clsat_v6'
NEW = '/staging/tmp/clsat_v6/scan85'

sources = {
    'CLsat1': f'{NEW}/new_candidates_out.bits.tsv',
    'CLsat2': f'{NEW}/newlocus_CLsat2.bits.tsv',
    'CLsat3': f'{NEW}/newlocus_CLsat3.bits.tsv',
    'CLsat4': f'{NEW}/newlocus_CLsat4.bits.tsv',
}

for fam, src in sources.items():
    existing = f'{V}/pass3/scan_{fam}.bits.tsv'
    added = 0
    with open(src) as f:
        header = f.readline()
        rows = f.readlines()
    with open(existing, 'a') as out:
        for line in rows:
            out.write(line)
            added += 1
    print(f'{fam}: appended {added} rows from {src} to {existing}')
