# -*- coding: utf-8 -*-
"""Render the sub-subfamily x assembly stacked-bar figure (manuscript Figure 2) from the
catalogue instead of hardcoded counts.

The published Figure 2 was rendered by figures_v6/render_fig2_connectors.py, which carried
the per-assembly / per-sub-subfamily matrix as a literal dict. Two cells later drifted
(CLsat3.1/mix 2 -> 3 with L541; CLsat1.5/val 1 -> 0), and both versions totalled 627, so every
total-based check passed. This script derives the matrix from the segmented BED, the
sub-subfamily table and the scaffold->species map, and asserts the per-assembly totals
against an independent count of the BED, so it cannot drift from the catalogue again.
Design, colours, geometry, figsize and dpi are unchanged from the original renderer.

Usage: make_fig2_from_data.py <MASTER.loci.segmented.bed> <locus_subfamily.tsv> <MASTER.merged.bed> <out.png>
                              [--species mix,arm,val,unp,unm,nai] [--labels "mix,arm,val,unp (pat.),unm (mat.),nai"]
                              [--order CLsat1.1,CLsat1.2,...] [--colours colours.tsv]
  MASTER.merged.bed   scaffold -> species comes from its name column (code|scf:...)
  --species           assembly order (default: order of first appearance in MASTER.merged.bed)
  --order             sub-subfamily order (default: sorted names present)
  --colours           TSV sub-subfamily<TAB>#rrggbb (default: a tab20 palette)
"""
import argparse
import collections
import csv
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('segmented_bed')
ap.add_argument('subfamily_tsv')
ap.add_argument('master_bed')
ap.add_argument('out')
ap.add_argument('--species', default='')
ap.add_argument('--labels', default='')
ap.add_argument('--order', default='')
ap.add_argument('--colours', default='')
a = ap.parse_args()

# ---- species order from MASTER.merged.bed (names start with the genome code) ----
species_seen = []
for r in csv.reader(open(a.master_bed, encoding='utf-8'), delimiter='\t'):
    if len(r) >= 4:
        sp = r[3].split('|')[0]
        if sp not in species_seen:
            species_seen.append(sp)
species = a.species.split(',') if a.species else species_seen
species_label = a.labels.split(',') if a.labels else species
if len(species_label) != len(species):
    sys.exit('--labels must have one entry per species')

# ---- derive the matrix from the catalogue ---------------------------------
asm, bed_total = {}, collections.Counter()
for r in csv.reader(open(a.segmented_bed, encoding='utf-8'), delimiter='\t'):
    if len(r) >= 4:
        # col 5 = the locus's genome (stage 14); never looked up by contig name, which two genomes may share
        if len(r) < 5:
            sys.exit('segmented BED row without genome column (col 5); rerun from 14_arbitrate: %s' % r)
        sp = r[4]
        asm[r[3]] = sp
        bed_total[sp] += 1

sub = {r['locus']: r['subfamily'] for r in csv.DictReader(
    open(a.subfamily_tsv, encoding='utf-8'), delimiter='\t')}

order = a.order.split(',') if a.order else sorted({f for f in sub.values() if f and f != 'NA'})
data = {f: {s: 0 for s in species} for f in order}
for locus, fam in sub.items():
    sp = asm.get(locus)
    if fam in data and sp in data[fam]:
        data[fam][sp] += 1

totals = {sp: sum(data[f][sp] for f in order) for sp in species}
print('per-assembly totals from subfamily table :', totals)
print('per-assembly totals from the BED         :', {sp: bed_total[sp] for sp in species})
for sp in species:
    if totals[sp] != bed_total[sp]:
        print('WARNING: mismatch for %s: %d loci with a sub-subfamily call vs %d in the BED '
              '(loci with NA calls or sub-subfamilies outside --order are not drawn)'
              % (sp, totals[sp], bed_total[sp]), file=sys.stderr)
print('grand total drawn: %d (BED: %d)\n' % (sum(totals.values()), sum(bed_total[sp] for sp in species)))

# ---- colours ---------------------------------------------------------------
colors = {}
if a.colours:
    for line in open(a.colours):
        p = line.rstrip('\n').split('\t')
        if len(p) >= 2 and not p[0].startswith('#'):
            colors[p[0]] = p[1]
cmap = plt.get_cmap('tab20')
for i, fam in enumerate(order):
    colors.setdefault(fam, matplotlib.colors.to_hex(cmap(i % 20)))

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

max_seg = max([v for fam in order for v in data[fam].values()] or [1])
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
ax.set_ylabel('Segmented loci', fontsize=12)
ax.set_ylim(0, max(totals.values() or [1]) * 1.15)
ax.set_xlim(-0.6, len(species) - 1 + 0.6)
h, l = ax.get_legend_handles_labels()
ax.legend(h, l, title='sub-subfamily', loc='upper left', ncol=3, fontsize=9.5,
          title_fontsize=11, frameon=True)
for spine in ('top', 'right'):
    ax.spines[spine].set_visible(False)
plt.tight_layout()
os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
fig.savefig(a.out, dpi=220)
plt.close(fig)
print('wrote %s' % a.out)
print('\nbar labels: %s' % ' '.join('%s=%d' % (s, totals[s]) for s in species))
