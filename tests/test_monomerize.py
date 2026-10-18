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


def test_header_rewrite():
    h = 'armf|h1tg000102l:2989387-3590591|len=292|copies=2059.7|locuslen=601205'
    out = m3b.monomerize_header(h, 2, 146, 292, 0.99)
    assert out == 'armf|h1tg000102l:2989387-3590591|len=146|copies=4119.4|locuslen=601205|monomerized=2:292:0.99'


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
    assert '|monomerized=2:%d:' % len(dimer) in lines[0] and '|len=%d|' % unit_len in lines[0] and '|copies=34.2|' in lines[0]
    assert lines[1] == dimer[:unit_len]
    assert lines[2] == '>sp|scf2:1-3000|len=146|copies=20|locuslen=3000'
    assert (tmp_path / 'sp.fa.before_monomerize').exists()
    rep = (tmp_path / 'sp.fa.monomerize.tsv').read_text().strip().split('\n')
    assert len(rep) == 2 and rep[1].split('\t')[2] == '2'
