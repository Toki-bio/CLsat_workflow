"""Test-only: `satpipe run` (Python driver) with every stage except 21_anchor replaced by a stub (WP5 acceptance).
The bash twin is tests/port/stub_stages.sh; both must stay line-for-line equivalent. Usage: stub_run.py run CFG [options]."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from satpipe import cli, run as R  # noqa: E402

# writes $1 = cat of the inputs + "<sid> <STUB_TAG>", logs 20 lines, fails when STUB_FAIL == sid
SCRIPT = ('o=$1; shift; cat "$@" > "$o"; printf "%s %s\\n" "$0" "${STUB_TAG:-none}" >> "$o"; i=0; '
          'while [ $i -lt 20 ]; do i=$((i + 1)); echo "log line $i of $0"; done; [ "$0" != "${STUB_FAIL:-}" ]')
DEPS = {'00_tables': '', '02_trf': '00_tables', '02b_large': '02_trf', '03_motifs': '00_tables', '04_filter': '03_motifs',
        '05_rotate': '04_filter', '06_cluster': '05_rotate', '07_mafft': '06_cluster', '08_cross': '07_mafft',
        '09_assign': '08_cross', '10_seeds': '09_assign', '11_expand': '10_seeds 02b_large', '12_master': '11_expand',
        '13_rescan': '12_master', '14_arbitrate': '13_rescan', '15_catalogue': '14_arbitrate', '16_regions': '15_catalogue',
        '17_terminal': '16_regions', '18_subfamily': '17_terminal', '19_orthology': '16_regions', '20_candidates': '15_catalogue 19_orthology'}


def deps(ctx, sid):
    return [p for p in (os.path.join(ctx.run, d, 'stub.out') for d in DEPS[sid].split()) if os.path.isfile(p)]


def gen(ctx, sid, extra=(), cwd=False, env=None, log_name=None):
    d = ctx.d(sid)
    o = os.path.join(d, 'stub.out')
    ctx.sh(['bash', '-c', SCRIPT, sid, 'stub.out' if cwd else o] + deps(ctx, sid) + list(extra), sid,
           cwd=d if cwd else None, env=env, log_name=log_name)
    return [o, d]


def g(ctx, k):
    return ctx.cfg.data[k]


def body(sid, ctx):
    data = ctx.cfg.data
    if sid == '00_tables':
        return gen(ctx, sid, R._cons_files(ctx))
    if sid == '02_trf':
        return gen(ctx, sid, R._genome_files(ctx), log_name='02_trf_split')
    if sid == '03_motifs':
        return gen(ctx, sid, cwd=True, env=R._disc_env(ctx))
    if sid == '07_mafft':
        gen(ctx, sid, cwd=True, log_name='07_mafft_a')
        return gen(ctx, sid, cwd=True, log_name='07_mafft_b', env=ctx.env(STUB_TAG='b'))
    if sid == '09_assign':
        d = ctx.d(sid)
        r = R._traced_run(sid, [ctx.tool('ssearch36'), '-n', '-z', '3'] + deps(ctx, sid), capture_output=True, text=True,
                          env=ctx.env())
        o = os.path.join(d, 'stub.out')
        with open(o, 'w', newline='\n') as f:
            f.write(r.stdout)
        n = sum(1 for _ in open(o))
        ctx.manifest['centroids_n'] = n
        ctx.manifest['centroids_md5'] = R.md5sum(o)
        ctx.say('  centroids: %d' % n)
        return [o]
    if sid == '10_seeds':
        d = ctx.d(sid)
        o = os.path.join(d, 'stub.out')
        with open(o, 'w') as f:
            R._traced_run(sid, ['bash', '-c', 'echo seeds; cat "$@"', 'mk'] + deps(ctx, sid), stdout=f, env=ctx.env(), check=True)
        ctx.manifest['seed_counts'] = {fam: i + 1 for i, fam in enumerate(ctx.families)}
        return [o]
    if sid == '11_expand':
        e, sc = data['expansion'], data['scan']
        for fam in ctx.families:
            fd = ctx.cfg.families[fam]
            extra = ['--min-bits', e['min_bits'], '--scan-win-bp', sc['win_mono'] * fd['mono'], '--jobs', e['jobs'],
                     '--cluster-id', data['discovery']['cluster_id']]
            if e.get('to_scaffold_ends'):
                extra.append('--to-scaffold-ends')
            ctx.sh(['bash', '-c', 'printf "%s\\n" "$@" >> "$0"', os.path.join(ctx.d(sid), 'argv.txt')] + extra, sid,
                   log_name='11_expand_' + fam)
        return gen(ctx, sid, [os.path.join(ctx.d(sid), 'argv.txt')])
    if sid == '12_master':
        ctx.manifest['master_loci'] = 627
        return gen(ctx, sid)
    if sid == '14_arbitrate':
        a = data['arbitration']
        bmb = a['bin_min_bits'] if a.get('bin_min_bits') is not None else data['expansion']['min_bits']
        return gen(ctx, sid, env=ctx.env(BIN_MIN_BITS=bmb, DELTA=a['delta'], STUB_TAG='arb'))
    if sid == '15_catalogue':
        ctx.manifest['segmented_loci'] = 650
        ctx.manifest['catalogue_loci'] = 651
        ctx.manifest['manual_loci'] = [{'locus': 'L1', 'bed': ['spa_scf1', '0', '10', 'x'], 'family': 'CLsat1',
                                        'bits': {'CLsat1': 12.5, 'CLsat2': 0.0}}]
        return gen(ctx, sid, cwd=True)
    if sid == '19_orthology':
        ctx.manifest['ortholog_pairs'] = 53
        return gen(ctx, sid)
    return gen(ctx, sid)


def plan(sid, ctx):
    data = ctx.cfg.data
    ins = deps(ctx, sid)
    if sid == '00_tables':
        rot = (data.get('references') or {}).get('rotated_consensus')
        return R._cons_files(ctx) + ([rot] if rot else []), {'genomes': ctx.cfg.genomes, 'families': ctx.cfg.families,
                                                             'members': (data.get('subfamily') or {}).get('members')}
    if sid == '02_trf':
        return ins + R._genome_files(ctx), {'trf': data['discovery']['trf'], 'large': int(data['finalisation']['large_scaffold_bp']),
                                            'genomes': ctx.cfg.genomes}
    if sid == '03_motifs':
        return ins, {'discovery': data['discovery'], 'step3b': 'runs-v3'}
    if sid == '09_assign':
        return [data['family_assignment']] + ins + R._cons_files(ctx), {}
    if sid == '11_expand':
        return ins + R._cons_files(ctx) + R._genome_files(ctx), {'expansion': data['expansion'], 'scan': data['scan'],
                                                                 'mono': {f: fd['mono'] for f, fd in ctx.cfg.families.items()}}
    if sid == '14_arbitrate':
        return ins, {'arbitration': data['arbitration'], 'min_bits': data['expansion']['min_bits'], 'rule': 'bin-majority-v1'}
    if sid == '16_regions':
        p = {'fin': data['finalisation']}
        hd = data.get('hapdup') or {}
        if hd.get('pairs'):
            p['hapdup'] = hd
        return ins, p
    if sid == '18_subfamily':
        return ins + R._subx_inputs(ctx), {'subfamily': data.get('subfamily'), 'anchor': data.get('anchor')}
    if sid == '19_orthology':
        return ins, {'orthology': data['orthology']}
    if sid == '20_candidates':
        return ins, {'min_bits': data['expansion']['min_bits'], 'flank': data['orthology']['flank_bp']}
    return ins, {}


def make(sid, needs):
    def fn(ctx):
        return body(sid, ctx)
    fn.sid, fn.needs = sid, needs
    fn.plan = lambda ctx: plan(sid, ctx)
    return fn


for i, f in enumerate(R.STAGES):
    if f.sid != '21_anchor':
        R.STAGES[i] = make(f.sid, f.needs)

sys.exit(cli.main(sys.argv[1:]))
