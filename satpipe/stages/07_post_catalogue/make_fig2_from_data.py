# -*- coding: utf-8 -*-
"""Regenerate Figure 2 from the deposited catalogue instead of hardcoded counts.

The published Figure 2 was rendered by figures_v6/render_fig2_connectors.py, which
carried the per-assembly / per-sub-subfamily matrix as a literal dict together with
`expected = {'mix': 128, ..., 'val': 80, ...}`. Two cells later drifted:

    CLsat3.1 / mix   2 -> 3   (the recovered mix CLsat3 locus, L541)
    CLsat1.5 / val   1 -> 0

so the figure showed mix = 128 and val = 80 while Table 2 and the Results give 129 and
79. Both versions total 627, which is why every total-based check passed.

This script derives the matrix from MASTER.loci.segmented.bed + locus_subfamily.tsv, so
the figure cannot drift from the catalogue again, and asserts the per-assembly totals
against an independent count of the BED. Design, colours, geometry, figsize and dpi are
unchanged from the original renderer, so the output is a drop-in replacement
(2090 x 1430 px).
"""
import collections
import csv
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

B = r'C:/work/CLsat/SUBMISSION/05_data_deposit'
OUT = r'C:/work/CLsat/manuscript/fig5_options/Figure_2_regenerated.png'

PREFIX = {'JAWWNF': 'arm', 'JAWWNG': 'val', 'JAWWNH': 'mix',
          'JAWWNI': 'nai', 'JAWDKN': 'unp', 'JAWDKO': 'unm'}
species = ['mix', 'arm', 'val', 'unp', 'unm', 'nai']
species_label = ['mix', 'arm', 'val', 'unp (pat.)', 'unm (mat.)', 'nai']
order = ['CLsat1.1', 'CLsat1.2', 'CLsat1.3', 'CLsat1.4', 'CLsat1.5',
         'CLsat2', 'CLsat3.1', 'CLsat3.2', 'CLsat4']
colors = {'CLsat1.1': '#08306b', 'CLsat1.2': '#2171b5', 'CLsat1.3': '#4292c6',
          'CLsat1.4': '#6baed6', 'CLsat1.5': '#c6dbef', 'CLsat2': '#e6550d',
          'CLsat3.1': '#00441b', 'CLsat3.2': '#41ab5d', 'CLsat4': '#6a51a3'}

# ---- derive the matrix from the deposit ----------------------------------
asm, bed_total = {}, collections.Counter()
for r in csv.reader(open(os.path.join(B, 'locus_catalogue', 'MASTER.loci.segmented.bed'),
                         encoding='utf-8'), delimiter='\t'):
    if len(r) >= 4:
        a = PREFIX.get(r[0][:6], '?')
        asm[r[3]] = a
        bed_total[a] += 1

sub = {r['locus']: r['subfamily'] for r in csv.DictReader(
    open(os.path.join(B, 'subfamily_assignment', 'locus_subfamily.tsv'), encoding='utf-8'),
    delimiter='\t')}

data = {f: {s: 0 for s in species} for f in order}
for locus, fam in sub.items():
    a = asm.get(locus)
    if a in data[fam] if fam in data else False:
        data[fam][a] += 1
    elif fam in data and a:
        data[fam][a] = data[fam].get(a, 0) + 1

totals = {sp: sum(data[f][sp] for f in order) for sp in species}
print('per-assembly totals from subfamily table :', totals)
print('per-assembly totals from the BED         :', dict(bed_total))
for sp in species:
    assert totals[sp] == bed_total[sp], 'mismatch for %s: %d vs %d' % (sp, totals[sp], bed_total[sp])
assert sum(totals.values()) == 627, 'grand total %d != 627' % sum(totals.values())
print('cross-check OK, grand total 627\n')

# ---- render (geometry identical to the original) -------------------------
BAR_W, MAX_LW = 0.65, 4.5

bounds = {sp: {} for sp in species}
for sp in species:
    bottom = 0
    for fam in order:
        bounds[sp][fam] = (bottom, bottom + data[fam][sp])
        bottom += data[fam][sp]


def smoothstep(x0, x1, y0, y1, n=40):
    t = np.linspace(0, 1, n)
    s = t * t * (3 - 2 * t)
    return x0 + (x1 - x0) * t, y0 + (y1 - y0) * s


fig, ax = plt.subplots(figsize=(9.5, 6.5))
x = np.arange(len(species))
bottom = np.zeros(len(species))
for fam in order:
    vals = np.array([data[fam][sp] for sp in species])
    ax.bar(x, vals, bottom=bottom, label=fam, color=colors[fam],
           edgecolor='white', linewidth=0.6, width=BAR_W, zorder=3)
    bottom += vals
for i, sp in enumerate(species):
    ax.text(i, bottom[i] + 3, str(totals[sp]), ha='center', va='bottom',
            fontsize=13, fontweight='bold', zorder=4)

max_seg = max(v for fam in order for v in data[fam].values())
for i in range(len(species) - 1):
    sl, sr = species[i], species[i + 1]
    xl, xr = i + BAR_W / 2, (i + 1) - BAR_W / 2
    for fam in order:
        bl0, bl1 = bounds[sl][fam]
        br0, br1 = bounds[sr][fam]
        if (bl1 - bl0) == 0 or (br1 - br0) == 0:
            continue
        lw = 0.6 + (((data[fam][sl] + data[fam][sr]) / 2) / max_seg) * MAX_LW
        xs, ys = smoothstep(xl, xr, (bl0 + bl1) / 2, (br0 + br1) / 2)
        ax.plot(xs, ys, color=colors[fam], linewidth=lw, alpha=0.85, zorder=2,
                solid_capstyle='round')

ax.set_xticks(x)
ax.set_xticklabels(species_label, rotation=20, ha='right', fontsize=11)
ax.set_ylabel('Segmented CLsat loci', fontsize=12)
ax.set_ylim(0, max(totals.values()) * 1.15)
ax.set_xlim(-0.6, len(species) - 1 + 0.6)
h, l = ax.get_legend_handles_labels()
ax.legend(h, l, title='sub-subfamily', loc='upper left', ncol=3, fontsize=9.5,
          title_fontsize=11, frameon=True)
for spine in ('top', 'right'):
    ax.spines[spine].set_visible(False)
plt.tight_layout()
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=220)
plt.close(fig)

from PIL import Image
im = Image.open(OUT)
print('wrote %s  (%dx%d)' % (OUT, im.width, im.height))
print('\nbar labels now: %s' % ' '.join('%s=%d' % (s, totals[s]) for s in species))
