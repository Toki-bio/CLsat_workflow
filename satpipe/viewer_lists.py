"""Build the viewer's three inspection lists (ORTHO_LOCI, FLAGGED_LOCI, FULL_TABLE) from the data
products instead of pasting them into index.html.

Inputs are the viewer data files themselves (so the lists can be rebuilt from any viewer
directory): the segmented catalogue BED with ids, the header-less assignment table, flanks.json,
flank_orthologs.json and subfamily.json.
"""
from __future__ import annotations

import json

FLAGGED = {'UNIQUE_SEED_CONTRADICTION', 'SEED_CONFLICT_RESOLVED_AMBIGUOUS', 'SEGMENT_RESOLVED_AMBIGUOUS',
           'NO_SEED_RESOLVED', 'NO_SEED_AMBIGUOUS', 'CANDIDATE_UNCATALOGUED'}


def read_master(path):
    """segmented catalogue rows: (id, sp, scf, start, end); id = 3rd pipe field of the name"""
    rows = []
    n = 0
    for line in open(path):
        p = line.rstrip('\n').split('\t')
        if len(p) < 4:
            continue
        n += 1
        parts = p[3].split('|')
        lid = parts[2] if len(parts) > 2 else 'L%d' % n
        rows.append((lid, parts[0], p[0], int(p[1]), int(p[2])))
    return rows


def read_assign(path):
    out = {}
    for line in open(path):
        p = line.rstrip('\n').split('\t')
        if len(p) < 7 or p[0].lower() in ('locus_id', 'locus'):
            continue
        out[p[0]] = {'fam': p[1], 'bits': p[2], 'other': p[3], 'other_bits': p[4], 'delta': p[5], 'flag': p[6],
                     'prov': p[7] if len(p) > 7 else ''}
    return out


def _fmt(x):
    try:
        return '%.1f' % float(x)
    except (TypeError, ValueError):
        return str(x)


def _best_partner(hits):
    """single best hit of one flank: highest bits (ties: higher pid, then id)"""
    return sorted(hits, key=lambda h: (-float(h.get('bits', 0)), -float(h.get('pid', 0)), h['locus']))[0]


def build_lists(master, assign, flanks, orthologs, subfam, threshold='80'):
    """master: list from read_master (catalogue order); assign: read_assign; flanks/orthologs/subfam:
    the loaded JSON products (may be None). Returns (ortho_loci, flagged_loci, full_table)."""
    flanks = flanks or {}
    orthologs = orthologs or {}
    subfam = subfam or {}
    fam_of = {lid: a['fam'] for lid, a in assign.items()}
    sp_of = {lid: sp for lid, sp, *_ in master}
    LF = (orthologs.get('LF') or {}).get(threshold, {})
    RF = (orthologs.get('RF') or {}).get(threshold, {})
    # orthology is undirected: a locus is a partner of every locus whose flank hit it, and of
    # every locus its own flanks hit. `sides` records whether the evidence came from an
    # L-sourced or an R-sourced hit list (own flank as query, or own locus as target).
    partners = {}     # lid -> {partner lid}
    pair_sides = {}   # frozenset{a,b} -> {'LF','RF'}: which hit lists (L- or R-sourced) support the pair
    forward = {}      # lid -> [partner lid] from its own flanks
    for key, tab in (('LF', LF), ('RF', RF)):
        for q, hits in tab.items():
            for h in hits:
                t = h['locus']
                partners.setdefault(q, set()).add(t)
                partners.setdefault(t, set()).add(q)
                forward.setdefault(q, []).append(t)
                pair_sides.setdefault(frozenset((q, t)), set()).add(key)

    full = []
    ortho = []
    flagged = []
    for lid, sp, scf, s, e in master:
        a = assign.get(lid)
        if a is None:
            continue
        fam = a['fam']
        if a['flag'] in FLAGGED:
            flagged.append({'id': lid, 'sp': sp, 'cat': a['flag'], 'f1': fam, 's1': a['bits'], 'f2': a['other'],
                            's2': a['other_bits'], 'd': a['delta'], 'fin': fam})
        if lid.startswith('C'):
            continue
        fl = flanks.get(lid) or {}
        lfc = ((fl.get('LF') or {}).get('class')) or 'none'
        rfc = ((fl.get('RF') or {}).get('class')) or 'none'
        row = {'id': lid, 'sp': sp, 'fam': fam}
        sf = subfam.get(lid) or {}
        row['sub'] = sf.get('subfamily')
        if sf.get('subfamily') and sf.get('confident'):
            row['subc'] = True
        row['lfc'] = lfc
        row['rfc'] = rfc
        for side, tab in (('lf', LF), ('rf', RF)):
            hits = tab.get(lid) or []
            if not hits:
                continue
            b = _best_partner(hits)
            row[side + 'p'] = b['locus']
            row[side + 'sp'] = b['sp']
            row[side + 'fam'] = fam_of.get(b['locus'], '')
            row[side + 'pid'] = _fmt(b['pid'])
            row[side + 'aln'] = str(b['aln'])
        mine = partners.get(lid, set())
        if mine:
            row['ortho'] = True
            # family agreement is judged on the best partner of each of the locus's OWN flanks
            # (the lfp/rfp columns); a locus that is only ever a target is judged on everything that hit it
            judged = [row[k] for k in ('lfp', 'rfp') if k in row] or sorted(mine)
            row['pfam'] = 'SAME' if all(fam_of.get(p) == fam for p in judged) else 'DIFFERENT'
        # both-sided = member of at least one partner pair "recovered from both flank sides"
        # (manuscript Methods 5): the pair is supported by an L-sourced AND an R-sourced hit, in
        # either direction. NOTE: the script that produced the published full_flank_ortho_table.tsv
        # (20 loci flagged, "18 of 53 pairs") was never deposited and its rule is not reproducible
        # from the shipped flank_orthologs.json (6 pairs / 11 loci under this definition); see
        # PROVENANCE.md. The viewer label says exactly what this flag means.
        if any(pair_sides.get(frozenset((lid, p)), set()) == {'LF', 'RF'} for p in mine):
            row['both'] = True
        full.append(row)
        nuf = (lfc == 'unique') + (rfc == 'unique')
        if nuf >= 1 or mine:
            by_sp = {}
            for p in mine:
                # a partner outside the catalogue means the orthology and catalogue files do not belong together;
                # filing it under '?' flipped the same/different-genome verdict (GLM audit #108)
                if p not in sp_of:
                    raise ValueError('ortholog partner %s of %s is not in the catalogue (flank_orthologs.json and '
                                     'MASTER.merged.bed are from different runs?)' % (p, lid))
                by_sp.setdefault(sp_of[p], set()).add(p)
            ortho.append({'id': lid, 'sp': sp, 'fam': fam, 'nuf': nuf, 'n': len(mine),
                          'partners': {k: sorted(v) for k, v in sorted(by_sp.items())}})
    # uncatalogued candidates first, then catalogue loci in catalogue order
    flagged.sort(key=lambda r: 0 if r['id'].startswith('C') else 1)
    return ortho, flagged, full


def load_viewer_dir(d):
    import os
    def j(name):
        p = os.path.join(d, name)
        return json.load(open(p)) if os.path.isfile(p) else None
    master = read_master(os.path.join(d, 'MASTER.merged.bed'))
    assign = read_assign(os.path.join(d, 'MASTER.final_assignment.clean.tsv'))
    return master, assign, j('flanks.json'), j('flank_orthologs.json'), j('subfamily.json')
