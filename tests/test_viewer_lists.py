"""viewer_lists.build_lists: with the decided pairs, a raw hit that is not a pair is no ortholog (owner 2026-11-08)."""
from satpipe.viewer_lists import build_lists, read_pairs

MASTER = [('L1', 'a', 's1', 0, 100), ('L2', 'b', 's1', 0, 100), ('L3', 'b', 's2', 0, 100)]
ASSIGN = {x: {'fam': 'F', 'bits': '1', 'other': '', 'other_bits': '', 'delta': '', 'flag': 'OK', 'prov': ''}
          for x in ('L1', 'L2', 'L3')}
FLANKS = {x: {'LF': {'class': 'unique'}, 'RF': {'class': 'unique'}} for x in ('L1', 'L2', 'L3')}
# L1's left flank hits L2 (best) and L3 (a tie or one-way hit that stage 19 did not keep)
HITS = {'LF': {'80': {'L1': [{'locus': 'L2', 'sp': 'b', 'pid': 99, 'aln': 500, 'bits': 900},
                             {'locus': 'L3', 'sp': 'b', 'pid': 98, 'aln': 500, 'bits': 880}]}}, 'RF': {'80': {}}}


def test_raw_hits_without_pairs_are_partners():
    ortho, _, _ = build_lists(MASTER, ASSIGN, FLANKS, HITS, {})
    assert {r['id']: r['n'] for r in ortho} == {'L1': 2, 'L2': 1, 'L3': 1}


def test_decided_pairs_only():
    ortho, _, full = build_lists(MASTER, ASSIGN, FLANKS, HITS, {}, pairs=[('L1', 'L2')])
    assert {r['id']: r['n'] for r in ortho} == {'L1': 1, 'L2': 1, 'L3': 0}
    row = {r['id']: r for r in full}
    assert row['L1']['lfp'] == 'L2' and 'ortho' not in row['L3']


def test_read_pairs(tmp_path):
    p = tmp_path / 'ortho_pairs.tsv'
    p.write_text('locusA\tlocusB\tspA\tspB\tpid\taln\tbits\nL1\tL2\ta\tb\t99.0\t500\t900.0\n')
    assert read_pairs(str(p)) == [('L1', 'L2')]
