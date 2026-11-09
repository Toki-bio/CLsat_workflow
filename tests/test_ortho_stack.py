"""Pure parts of satpipe.ortho_stack: groups, orientation, monomer boundaries, anchors, layout, flank windows."""
from satpipe import ortho_stack as S


def test_groups_components_and_order():
    sp = {'L1': 'Flat', 'L2': 'lat', 'L3': 'mex', 'L9': 'lat', 'L10': 'Fmex'}
    g = S.groups([('L1', 'L2'), ('L3', 'L1'), ('L9', 'L10')], sp, ['lat', 'Flat', 'Fmex', 'mex'])
    assert g == [['L2', 'L1', 'L3'], ['L9', 'L10']]


def test_pair_relation_and_orient():
    assert S.pair_relation({('L', 'L'): 300}) == 1
    assert S.pair_relation({('L', 'R'): 300, ('R', 'L'): 200}) == -1
    assert S.pair_relation({}) is None
    o, unk = S.orient(['A', 'B', 'C'], [('A', 'B'), ('C', 'B')], {('A', 'B'): -1, ('C', 'B'): -1})
    assert o == {'A': 1, 'B': -1, 'C': 1} and unk == []
    o, unk = S.orient(['A', 'B'], [('A', 'B')], {('A', 'B'): None})
    assert o == {'A': 1, 'B': 1} and unk == [('A', 'B')]


def test_boundaries_from_hits_phase():
    mono, n = 100, 1000
    # windows of 200 at step 100; array starts at consensus position 31 -> monomers start at 70, 170, ...
    hits = []
    for w in range(0, 801, 100):
        qs = 31 + w % mono          # consensus position at window letter 1
        hits.append((w, qs, qs + 199, 1, 200, 300.0))
    b = S.boundaries_from_hits(hits, mono, n)
    assert b == list(range(70, 1000, 100))
    assert S.segments(n, b)[0] == (0, 70) and S.segments(n, b)[-1] == (970, 1000)


def test_clear_anchors_reciprocal_and_collinear():
    bits = {(0, 0): 100, (1, 1): 100, (1, 2): 95,      # x1: not clear (95 >= 0.9 x 100)
            (2, 3): 100, (3, 2): 50, (4, 1): 100}      # (4,1): clear for x4 but y1's best... crosses the chain
    a = S.clear_anchors(bits, 0.9)
    assert (1, 1) not in a and a[0] == (0, 0) and (2, 3) in a
    assert all(a[i][0] < a[i + 1][0] and a[i][1] < a[i + 1][1] for i in range(len(a) - 1))


def test_layout_gap_and_leftshift():
    # ref 6 monomers; row B lost ref monomers 2,3 (anchors 0-0, 1-1, 2-4, 3-5); row C: nothing clear (3 monomers)
    keys, place = S.layout(6, {'B': 4, 'C': 3}, {'B': [(0, 0), (1, 1), (2, 4), (3, 5)], 'C': []})
    assert [k for k, _ in place['B']] == [(0, 0), (1, 0), (4, 0), (5, 0)]
    assert all(a for _, a in place['B'])
    assert place['C'] == [((0, 0), False), ((1, 0), False), ((2, 0), False)]
    # an extra monomer between adjacent anchors gets an inserted column
    keys, place = S.layout(2, {'D': 3}, {'D': [(0, 0), (2, 1)]})
    assert [k for k, _ in place['D']] == [(0, 0), (0, 1), (1, 0)] and keys == [(0, 0), (0, 1), (1, 0)]


def test_flank_windows_from_edge():
    # L flanks: the array is on the right; windows count from the right end
    aln = {'a': 'A' * 50 + 'C' * 100, 'b': '-' * 50 + 'C' * 50 + 'G' * 50}
    n, pr, idt = S.flank_windows(aln, {'a': set(range(10)), 'b': set()}, 'L', win=100)
    assert n == 2
    assert idt['a|b'][0] == 50.0                  # edge window: 50 C=C, 50 C/G
    assert pr['a']['masked'][0] == 0.1 and pr['b']['letters'] == [100, 0]
    assert idt['a|b'][1] == 0.0                   # 50 letters of a against gaps
