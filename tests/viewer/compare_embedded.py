#!/usr/bin/env python3
"""Rebuild the viewer's inspection lists from a viewer data directory and compare them with the
lists embedded in that directory's index.html (the v2.25 audit build). Development check for
satpipe.viewer_lists; exit 0 when all three lists are identical.

Usage: compare_embedded.py <viewer_dir>
"""
import json
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from satpipe.viewer_lists import build_lists, load_viewer_dir  # noqa: E402

d = sys.argv[1]
src = open(os.path.join(d, 'index.html'), encoding='utf-8').read()
emb = {}
for name in ('ORTHO_LOCI', 'FLAGGED_LOCI', 'FULL_TABLE'):
    m = re.search(r'^const %s = (.*?);\s*$' % name, src, re.M)
    emb[name] = json.loads(m.group(1))

master, assign, flanks, orth, sub = load_viewer_dir(d)
ortho, flagged, full = build_lists(master, assign, flanks, orth, sub)
# the viewer sorts FLAGGED_LOCI (contradictions first) only at load time; compare as sets of rows
ok = True
for name, mine in (('ORTHO_LOCI', ortho), ('FLAGGED_LOCI', flagged), ('FULL_TABLE', full)):
    theirs = emb[name]
    a = {r['id']: r for r in theirs}
    b = {r['id']: r for r in mine}
    if len(a) != len(theirs) or len(b) != len(mine):
        print(name, 'duplicate ids')
    missing = sorted(set(a) - set(b)); extra = sorted(set(b) - set(a))
    diffs = [(k, a[k], b[k]) for k in a if k in b and a[k] != b[k]]
    order_ok = [r['id'] for r in theirs] == [r['id'] for r in mine]
    # FULL_TABLE row order is not part of the contract (the viewer filters and renders rows as a set)
    same = not missing and not extra and not diffs and (order_ok or name == 'FULL_TABLE')
    print('%-13s embedded %d  rebuilt %d  %s' % (name, len(theirs), len(mine), 'IDENTICAL' if same else 'DIFF'))
    if not same:
        ok = False
        if missing: print('   missing:', missing[:10])
        if extra: print('   extra:', extra[:10])
        ia, ib = [r['id'] for r in theirs], [r['id'] for r in mine]
        if ia != ib and not missing and not extra:
            k = next(i for i in range(min(len(ia), len(ib))) if ia[i] != ib[i])
            print('   order differs at index %d: embedded %s rebuilt %s' % (k, ia[k:k+5], ib[k:k+5]))
        for k, x, y in diffs[:6]:
            print('   %s\n     emb %s\n     new %s' % (k, json.dumps(x, sort_keys=True), json.dumps(y, sort_keys=True)))
        if diffs:
            keys = {}
            for k, x, y in diffs:
                for f in set(x) | set(y):
                    if x.get(f) != y.get(f): keys[f] = keys.get(f, 0) + 1
            print('   differing fields:', keys)
sys.exit(0 if ok else 1)
