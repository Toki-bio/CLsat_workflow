#!/usr/bin/env python3
"""Write viewer_config.json for an EXISTING viewer data directory (one that was deployed with the
constants embedded in index.html), from its shipped data files plus a small metadata JSON with
the things that used to live in the page: family order/monos/colours, species order and labels,
sub-subfamily consensuses, title. Used to migrate the live Darevskia viewer to the config-driven
page without re-running the pipeline.

Usage: make_config_from_viewer_dir.py <viewer_dir> <meta.json> [out.json]
meta.json: {"title": ..., "run_name": ..., "families": [{"name","mono","colour"}], "species": [{"code","label"}],
            "subfamilies": {"CLsat1.1": {"family": "CLsat1", "seq": "..."}}, "sat_bits": 50, "edge_tol": 50,
            "flank_view": 5000, "orthology": {"min_pid": 80, "min_aln": 200, "flank_bp": 600},
            "arbitration": {"delta": 20, "delta_ambig": 15}}
The consensus of each family is read from <viewer_dir>/<family>.cons.fa.
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from satpipe.viewer_lists import build_lists, load_viewer_dir  # noqa: E402

d, meta_path = sys.argv[1], sys.argv[2]
out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(d, 'viewer_config.json')
meta = json.load(open(meta_path))
master, assign, flanks, orth, sub = load_viewer_dir(d)
ol, fl, full = build_lists(master, assign, flanks, orth, sub)


def read_fa(p):
    return ''.join(l.strip() for l in open(p) if not l.startswith('>')).upper()


fams = []
for f in meta['families']:
    fams.append({'name': f['name'], 'mono': f['mono'], 'colour': f['colour'],
                 'consensus': read_fa(os.path.join(d, f['name'] + '.cons.fa')),
                 'bits_file': 'MASTER_%s.bits.tsv' % f['name']})
n_cat = sum(1 for r in master if not r[0].startswith('C'))
n_cand = sum(1 for r in master if r[0].startswith('C'))
n_flag = sum(1 for r in fl if not r['id'].startswith('C'))
n_contra = sum(1 for r in fl if r['cat'] == 'UNIQUE_SEED_CONTRADICTION')
o = meta.get('orthology', {'min_pid': 80, 'min_aln': 200, 'flank_bp': 600})
arb = meta.get('arbitration', {'delta': 20, 'delta_ambig': 15})
codes = [x['code'] for x in meta['species']]
cfg = {
    'generator': 'tests/viewer/make_config_from_viewer_dir.py',
    'title': meta.get('title', 'Satellite locus viewer'),
    'run_name': meta.get('run_name', ''),
    'families': fams,
    'species': meta['species'],
    'subfamilies': meta.get('subfamilies', {}),
    'window_step': max(f['mono'] for f in fams),
    'sat_bits': meta.get('sat_bits', 50), 'edge_tol': meta.get('edge_tol', 50), 'flank_view': meta.get('flank_view', 5000),
    'orthology': o,
    'counts': {'catalogue_loci': n_cat, 'candidates': n_cand, 'flagged': n_flag, 'contradictions': n_contra,
               'anchorable': len(ol), 'with_ortholog': sum(1 for r in ol if r['n'] > 0)},
    'lists': {'ortho': ol, 'flagged': fl, 'full': full},
    'notes': meta.get('notes') or {
        'ortho': 'Positional orthology from unique flanking sequence (>=%g%% identity, >=%d bp), single best hit per target '
                 'assembly. "Orthologous to" lists distinct partner assemblies then the total partner-locus count.' % (o['min_pid'], o['min_aln']),
        'flag': '%d catalogue loci annotated by arbitration%s. %d %s a seed/rescan contradiction (the rescan scored another '
                'family >=%g bits above the seed; the seed family is retained). The rest are close-margin calls (margin below '
                '%g bits) or loci without a seed.' % (n_flag, (', plus the %d uncatalogued candidates (not catalogue members)' % n_cand) if n_cand else '',
                                                   n_contra, 'is' if n_contra == 1 else 'are', arb['delta'], arb['delta_ambig']),
        'full': 'Every locus (%d, all %d assemblies) with its sub-subfamily call and both flanks\' cross-assembly orthology status '
                '(ssearch36, >=%g%% identity, >=%d bp, single best hit per target assembly). L/R badges show each flank\'s class; '
                'the double tick marks loci in a partner pair supported by both an L- and an R-sourced hit.' % (n_cat, len(codes), o['min_pid'], o['min_aln']),
    },
}
json.dump(cfg, open(out, 'w'), separators=(',', ':'))
print('wrote %s: %d families, %d species, lists %d/%d/%d' % (out, len(fams), len(codes), len(ol), len(fl), len(full)))
