import importlib.util
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    'flank_rules', os.path.join(HERE, '..', 'satpipe', 'stages', '04_flank_orthology', 'flank_rules.py'))
fr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fr)


def rnd(n, seed):
    r = random.Random(seed)
    return ''.join(r.choice('ACGT') for _ in range(n))


def test_ordinary_dna_is_unique():
    # the old rule called most ordinary 600-letter flanks low-complexity
    calls = [fr.classify(rnd(600, s), 200) for s in range(200)]
    assert calls.count('unique') >= 195


def test_telomere_is_telomeric_not_lowcomplex():
    assert fr.classify('TTAGGG' * 100, 200) == 'telomeric'
    assert fr.classify('CCCTAA' * 100, 200) == 'telomeric'


def test_worn_telomere_in_simple_flank_is_telomeric():
    seq = 'CCCTAA' * 4 + 'AT' * 288
    assert fr.classify(seq, 200) == 'telomeric'


def test_telomere_motifs_in_ordinary_dna_stay_unique():
    seq = rnd(300, 1) + 'TTAGGG' * 4 + rnd(276, 2)
    assert fr.classify(seq, 200) == 'unique'


def test_simple_repeat_is_lowcomplex():
    assert fr.classify('AT' * 300, 200) == 'lowcomplex'
    assert fr.classify('AAAAC' * 120, 200) == 'lowcomplex'


def test_short_is_truncated():
    assert fr.classify(rnd(150, 3), 200) == 'truncated'
