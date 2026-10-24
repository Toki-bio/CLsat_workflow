import importlib.util
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, '..', 'satpipe', 'stages', '01_discovery_clustering', 'Step3b_monomerize.py')
spec = importlib.util.spec_from_file_location('step3b', SCRIPT)
m3b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m3b)
CONS = os.path.join(HERE, '..', 'satpipe', 'stages', '02_seed_expansion_scan', 'consensus')


def rnd(n, r, alphabet='ACGT'):
    return ''.join(r.choice(alphabet) for _ in range(n))


def mutate(s, r, subs=0.05, indels=2):
    s = list(s)
    for i in range(len(s)):
        if r.random() < subs:
            s[i] = r.choice([c for c in 'ACGT' if c != s[i]])
    for _ in range(indels):
        i = r.randrange(len(s))
        if r.random() < 0.5:
            del s[i]
        else:
            s.insert(i, r.choice('ACGT'))
    return ''.join(s)


def cons(name):
    return ''.join(l.strip() for l in open(os.path.join(CONS, name + '.cons.fa')) if not l.startswith('>')).upper()


def k_of(seq, pmin=80):
    return m3b.unit_count(seq, pmin, 0.70, 12)[0]


def test_single_unit_stays():
    r = random.Random(1)
    for n in (80, 146, 220, 300):
        assert k_of(rnd(n, r)) == 1


def test_random_and_at_rich_never_called():
    r = random.Random(2)
    calls = sum(k_of(rnd(r.randint(160, 300), r, 'AATTGC')) > 1 for _ in range(300))
    calls += sum(k_of(rnd(r.randint(160, 300), r)) > 1 for _ in range(300))
    assert calls == 0


def test_identical_dimer():
    u = rnd(146, random.Random(3))
    k, idv = m3b.unit_count(u + u, 80, 0.70, 12)
    assert k == 2 and idv > 0.99


def test_diverged_dimer_with_indels():
    r = random.Random(4)
    for _ in range(50):
        u = rnd(146, r)
        assert k_of(u + mutate(u, r, subs=0.05, indels=2)) == 2


def test_trimer_and_short_unit_dimer():
    r = random.Random(5)
    u = rnd(95, r)
    assert k_of(u + mutate(u, r, 0.03, 1) + mutate(u, r, 0.03, 1)) == 3
    v = rnd(84, r)
    assert k_of(v + mutate(v, r, 0.04, 1)) == 2


def test_exact_tetramer_and_hexamer_cut_to_the_unit():
    # a shift by 2 units matches as well as a shift by 1 unit; the unit (largest k) must win (GLM audit #31)
    r = random.Random(12)
    u = rnd(60, r)
    assert m3b.unit_count(u * 4, 50, 0.70, 12)[0] == 4
    assert m3b.unit_count(u * 6, 50, 0.70, 12)[0] == 6
    v = rnd(80, r)
    # mutated tetramer (indels make it a little shorter than 4 x 80, so period_min 75 lets k = 4 be tried)
    assert m3b.unit_count(v + mutate(v, r, 0.02, 1) + mutate(v, r, 0.02, 1) + mutate(v, r, 0.02, 1), 75, 0.70, 12)[0] == 4


def test_real_consensus_dimers():
    for fam in ('CLsat1', 'CLsat2', 'CLsat3', 'CLsat4'):
        c = cons(fam)
        k, idv = m3b.unit_count(c + c, 80, 0.70, 12)
        assert k == 2 and idv > 0.99, fam
        assert k_of(c) == 1, fam


def test_two_different_units_joined_stay():
    r = random.Random(6)
    assert k_of(cons('CLsat1') + rnd(146, r)) == 1


def test_unit_never_below_period_min():
    # (AATGG)n is periodic, but its unit (5) is far below period_min: no k is even tried
    s = ('AATGG' * 40)[:200]
    assert k_of(s, pmin=80) == 2 or k_of(s, pmin=80) == 1   # k <= 200 // 80 = 2: a 100-letter "unit" is allowed
    assert all(len(s) // k >= 80 for k in range(2, len(s) // 80 + 1))


def p_of(seq, pmin=80):
    return m3b.partial_period(seq, pmin, 0.80, 60, 12)[0]


def test_partial_unit_cut_to_one_unit():
    # the mixta case: a 146-letter CLsat1 unit plus its first 80 letters again (TRF period 215, 226 letters)
    c = cons('CLsat1')
    assert p_of(c + c[:80]) == len(c)
    r = random.Random(7)
    for _ in range(50):
        u = rnd(146, r)
        p = p_of(u + mutate(u[:r.randint(60, 120)], r, subs=0.05, indels=1))
        assert p is not None and abs(p - 146) <= 2


def test_partial_never_on_single_units_or_random():
    for fam in ('CLsat1', 'CLsat2', 'CLsat3', 'CLsat4'):
        assert p_of(cons(fam)) is None, fam
    r = random.Random(8)
    calls = sum(p_of(rnd(r.randint(140, 300), r)) is not None for _ in range(300))
    calls += sum(p_of(rnd(r.randint(140, 300), r, 'AATTGC')) is not None for _ in range(300))
    assert calls == 0


def test_partial_in_file_round_trip(tmp_path):
    c = cons('CLsat1')
    src = tmp_path / 'm.fa'
    src.write_text('>mix|ptg000337l_pilon:105428-128449|len=215|copies=105|locuslen=23022\n%s\n' % (c + c[:80]))
    subprocess.run([sys.executable, SCRIPT, str(src)], check=True, capture_output=True)
    lines = src.read_text().split('\n')
    assert lines[1] == c
    assert '|len=146|' in lines[0] and '|monomerized=partial:226:' in lines[0]


def test_header_rewrite():
    h = 'armf|h1tg000102l:2989387-3590591|len=292|copies=2059.7|locuslen=601205'
    out = m3b.monomerize_header(h, 2, 146, 292, 0.99)
    assert out == 'armf|h1tg000102l:2989387-3590591|len=146|copies=4119.4|locuslen=601205|monomerized=2:292:0.99'


def test_seeds_of_a_bridging_trf_interval():
    # the real case: one period-292 TRF row over 387 kb spanning two CLsat1 arrays and the 4.75 kb gap between them
    h = 'unm|JAWDKO010000018.1:156359-543662|len=146|copies=2649.6|locuslen=387304|monomerized=2:292:1.00'
    s, e = m3b.seed_headers(h, 146, 10)
    assert s == 'unm|JAWDKO010000018.1:156359-157818|len=146|copies=2649.6|locuslen=1460|monomerized=2:292:1.00|seed=start'
    assert e == 'unm|JAWDKO010000018.1:542203-543662|len=146|copies=2649.6|locuslen=1460|monomerized=2:292:1.00|seed=end'
    # no seed reaches into the gap (270231-274981)
    for x in (s, e):
        a, b = (int(v) for v in x.split('|')[1].split(':')[1].split('-'))
        assert b < 270231 or a > 274981


def test_scan_runs_splits_a_bridging_interval(tmp_path):
    import shutil
    if not (shutil.which('samtools') and shutil.which('ssearch36')):
        import pytest
        pytest.skip('needs samtools and ssearch36')
    r = random.Random(11)
    c = cons('CLsat1')
    arr = lambda n: ''.join(mutate(c, r, 0.03, 0) for _ in range(n))
    a1, gap, a2 = arr(40), rnd(4000, r), arr(50)
    seq = rnd(3000, r) + a1 + gap + a2 + rnd(3000, r)
    g = tmp_path / 'g.fna'
    g.write_text('>scf\n' + seq + '\n')
    subprocess.run(['samtools', 'faidx', str(g)], check=True)
    s, e = 3001, 3000 + len(a1) + len(gap) + len(a2)
    cf = os.path.join(CONS, 'CLsat1.cons.fa')
    runs = m3b.scan_runs(str(g), 'scf', s, e, [cf])
    assert len(runs) == 2
    (b1, e1), (b2, e2) = runs
    assert abs(b1 - s) <= 300 and abs(e1 - (s + len(a1) - 1)) <= 300
    assert abs(b2 - (s + len(a1) + len(gap))) <= 300 and abs(e2 - e) <= 300
    assert b2 - e1 > 3000          # the gap is not inside any run


def test_short_interval_keeps_one_seed():
    h = 'sp|scf:1001-3000|len=146|copies=13.7|locuslen=2000|monomerized=2:292:0.99'
    assert m3b.seed_headers(h, 146, 10) == [h]


def test_file_round_trip(tmp_path):
    r = random.Random(7)
    u = rnd(146, r)
    dimer = u + mutate(u, r)
    unit_len = round(len(dimer) / 2)
    fa = tmp_path / 'sp.fa'
    fa.write_text('>sp|scf1:1-5000|len=%d|copies=17.1|locuslen=5000\n%s\n>sp|scf2:1-3000|len=146|copies=20|locuslen=3000\n%s\n'
                  % (len(dimer), dimer, rnd(146, r)))
    subprocess.run([sys.executable, SCRIPT, str(fa)], check=True, capture_output=True,
                   env=dict(os.environ, PERIOD_MIN='80'))
    lines = fa.read_text().split('\n')
    # the 5000-letter TRF interval is longer than 2 x 10 units: one start seed and one end seed, same unit
    assert '|monomerized=2:%d:' % len(dimer) in lines[0] and '|len=%d|' % unit_len in lines[0] and '|copies=34.2|' in lines[0]
    w = 10 * unit_len
    assert lines[0].startswith('>sp|scf1:1-%d|' % w) and lines[0].endswith('|seed=start')
    assert lines[2].startswith('>sp|scf1:%d-5000|' % (5000 - w + 1)) and lines[2].endswith('|seed=end')
    assert lines[1] == lines[3] == dimer[:unit_len]
    assert lines[4] == '>sp|scf2:1-3000|len=146|copies=20|locuslen=3000'
    assert (tmp_path / 'sp.fa.before_monomerize').exists()
    rep = (tmp_path / 'sp.fa.monomerize.tsv').read_text().strip().split('\n')
    assert len(rep) == 2 and rep[1].split('\t')[2] == '2'
