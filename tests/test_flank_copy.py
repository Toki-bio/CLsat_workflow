"""Single-copy anchor rules of stage 19 (flank_copy.py; author's decisions 2026-11-08). Pure functions, no tools."""
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    'flank_copy', os.path.join(HERE, '..', 'satpipe', 'stages', '04_flank_orthology', 'flank_copy.py'))
FC = importlib.util.module_from_spec(spec)
spec.loader.exec_module(FC)

DAT = """Tandem Repeats Finder Program written by:

Sequence: scfA



Parameters: 2 5 7 80 10 50 2000


100 199 10 10.0 10 95 0 180 25 25 25 25 2.0 ACGTACGTAC ACGTACGTAC
150 260 10 11.0 10 95 0 190 25 25 25 25 2.0 ACGTACGTAC ACGTACGTAC
1000 1049 25 2.0 25 90 0 80 25 25 25 25 2.0 ACGT ACGT
Sequence: scfB



Parameters: 2 5 7 80 10 50 2000


5 50 15 3.0 15 88 0 70 25 25 25 25 2.0 ACG ACG
"""


def test_parse_trf_dat_merges_and_keeps_only_wanted(tmp_path):
    p = tmp_path / 'g.dat'
    p.write_text(DAT)
    iv = FC.parse_trf_dat(str(p), {'scfA'})
    assert iv == {'scfA': [(100, 260), (1000, 1049)]}          # overlapping calls merged, scfB not stored
    assert FC.parse_trf_dat(str(p), set()) == {}


def test_merge_adjacent_and_nested():
    assert FC.merge([(10, 20), (21, 30), (12, 15), (40, 50)]) == [(10, 30), (40, 50)]


def test_mask_counts_only_letters_inside_flank():
    seq = 'A' * 100                      # flank at 1-based 201..300
    m, n = FC.mask(seq, 201, [(150, 210), (290, 400), (500, 600)])
    assert n == 10 + 11                  # 201..210 and 290..300
    assert m[:10] == 'N' * 10 and m[10] == 'A' and m[-11:] == 'N' * 11 and FC.unmasked(m) == 79


def test_mask_no_intervals_is_identity():
    assert FC.mask('ACGT', 5, []) == ('ACGT', 0)


def paf(q, t, ts, te, match, block):
    return '\t'.join(map(str, [q, 600, 0, 600, '+', t, 10**6, ts, te, match, block, 60]))


def test_paf_self_place_not_counted_second_copy_counted():
    place = {'f1': ('chr1', 1001, 1600)}           # 1-based inclusive genomic place of flank f1
    lines = [paf('f1', 'chr1', 1000, 1600, 600, 600),     # itself (0-based half-open 1000..1600)
             paf('f1', 'chr2', 5000, 5590, 500, 590),     # second copy: 84.7 % over 590
             paf('f1', 'chr1', 900000, 900400, 290, 400)]  # 72.5 %: below 80 %, not a copy
    n = FC.paf_other_copies(lines, place, 80.0, 200)
    assert n['f1'] == 1


def test_paf_short_or_weak_hits_do_not_count():
    place = {'f1': ('chr1', 1, 600)}
    lines = [paf('f1', 'chr3', 0, 150, 150, 150),          # 100 % but only 150 letters
             paf('f1', 'chr3', 0, 300, 230, 300)]          # 76.7 %
    assert FC.paf_other_copies(lines, place, 80.0, 200)['f1'] == 0


def test_paf_copy_on_same_scaffold_elsewhere_counts():
    place = {'f1': ('chr1', 1001, 1600)}
    lines = [paf('f1', 'chr1', 20000, 20600, 590, 600)]
    assert FC.paf_other_copies(lines, place, 80.0, 200)['f1'] == 1


def test_clear_best_single_locus():
    assert FC.clear_best([(500.0, 'L2', 95.0, 600)], 0.9)[1] == 'L2'


def test_clear_best_exact_tie_undecided():
    assert FC.clear_best([(500.0, 'L2'), (500.0, 'L3')], 0.9) is None


def test_clear_best_near_tie_undecided_and_clear_win():
    assert FC.clear_best([(500.0, 'L2'), (460.0, 'L3')], 0.9) is None       # 0.92 x best
    assert FC.clear_best([(500.0, 'L2'), (440.0, 'L3')], 0.9)[1] == 'L2'    # 0.88 x best


def test_clear_best_both_flanks_of_one_locus_are_not_a_tie():
    # a query can hit both flanks of the same target locus; that is one locus, not two competitors
    assert FC.clear_best([(500.0, 'L2'), (498.0, 'L2'), (100.0, 'L9')], 0.9)[1] == 'L2'


def test_clear_best_order_independent():
    a = [(300.0, 'L5'), (500.0, 'L2'), (100.0, 'L3')]
    assert FC.clear_best(a, 0.9) == FC.clear_best(list(reversed(a)), 0.9)


def test_reciprocal_pairs():
    keep, one = FC.reciprocal_pairs({('A', 'B'), ('B', 'A'), ('A', 'C'), ('D', 'E')})
    assert keep == {('A', 'B')}
    assert one == {('A', 'C'), ('D', 'E')}
