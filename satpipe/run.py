"""The satpipe driver: a plain stage runner.

Each stage is a function that takes the run context, calls the published stage script with
explicit arguments and environment, and writes <run>/<stage>/DONE holding the fingerprints of
its inputs and the hash of the config section it depends on. A stage is skipped when its DONE
file matches; `--force` re-runs it. Everything a stage needs comes from the config: no path,
species code, family name or monomer length lives in this file.

Run layout (<out_dir>):
  tables/        genomes.tsv families.tsv members.tsv refs.rot (derived or copied)
  trf/           <code>.dat  (+ <code>.split/ when TRF runs here)  large_scaffolds.tsv
  large/<code>/<scaffold>/   chunked TRF of scaffolds above large_scaffold_bp + recovered seeds
  motifs/        Step 3-7 (motifs_output layout: <code>.fa, .clean.fas, clean/*.rot, clean/clusters_by_file/)
  cross/         Step 8 (centroids, cross_clusters/) + cluster_report.tsv (stage 09 validation)
  seeds/         <family>.members.fa  <family>.bed
  expand/        out_<family>.*  expansions.tsv
  master/        MASTER.merged.bed
  rescan/        MASTER_<family>.bits.tsv  arbitration_inputs.tsv
  arbitration/   stage-14 outputs
  catalogue/     final catalogue = arbitration + manual loci
  regions/       locus_master.tsv manifest.tsv all_regions.fa terminal/phase/subfamily tables + JSON
  ortho/         flank orthology
  logs/          one log per stage;  manifest.json;  DONE markers inside each stage directory
"""
from __future__ import annotations

import concurrent.futures
import datetime as dt
import hashlib
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field

from . import __version__
from .config import Config
from .stages import scripts_dir
from .tools import probe_all

SMALL = 64 * 1024 * 1024   # files up to this size are sha256-hashed; larger ones are fingerprinted by size+mtime


class StageError(Exception):
    pass


# ----------------------------------------------------------------------------- helpers

def fingerprint(path: str) -> str:
    st = os.stat(path)
    if st.st_size <= SMALL:
        h = hashlib.sha256()
        with open(path, 'rb') as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b''):
                h.update(chunk)
        return 'sha256:' + h.hexdigest()
    return 'size-mtime:%d:%d' % (st.st_size, int(st.st_mtime))


def md5sum(path: str) -> str:
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read_fasta_names(path: str):
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                yield line[1:].split()[0]


def minrot(s: str) -> str:
    d = s + s
    return min(d[i:i + len(s)] for i in range(len(s)))


def revcomp(s: str) -> str:
    return s.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]


def read_tsv(path, min_cols=1):
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) >= min_cols:
                rows.append(p)
    return rows


@dataclass
class Ctx:
    cfg: Config
    run: str
    force: bool = False
    log: object = None
    manifest: dict = field(default_factory=dict)

    # directories ----------------------------------------------------------
    def d(self, *parts) -> str:
        p = os.path.join(self.run, *parts)
        os.makedirs(p, exist_ok=True)
        return p

    @property
    def scripts(self) -> str:
        return scripts_dir()

    @property
    def families(self) -> list:
        return list(self.cfg.families)

    @property
    def codes(self) -> list:
        return list(self.cfg.genomes)

    def env(self, **extra) -> dict:
        e = dict(os.environ)
        e['LC_ALL'] = 'C'
        bd = self.cfg.bin_dir
        if bd:
            e['PATH'] = bd + os.pathsep + e.get('PATH', '')
        for name, path in (self.cfg.data.get('tools') or {}).items():
            if name != 'bin_dir':
                e['SATPIPE_TOOL_' + name.upper()] = path
        for k, v in extra.items():
            e[k] = str(v)
        return e

    def tool(self, name: str) -> str:
        return (self.cfg.data.get('tools') or {}).get(name) or name

    # process ----------------------------------------------------------------
    def sh(self, cmd: list, stage: str, cwd: str | None = None, env: dict | None = None, log_name: str | None = None):
        logf = os.path.join(self.d('logs'), (log_name or stage) + '.log')
        self.say('  $ ' + ' '.join(str(c) for c in cmd) + (f'   (cwd {cwd})' if cwd else ''))
        with open(logf, 'a') as lf:
            lf.write('\n### %s  %s\n$ %s\n' % (dt.datetime.now().isoformat(timespec='seconds'), cwd or '', ' '.join(map(str, cmd))))
            lf.flush()
            r = subprocess.run([str(c) for c in cmd], cwd=cwd, env=env or self.env(), stdout=lf, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            tail = ''.join(open(logf).readlines()[-15:])
            raise StageError('stage %s: command failed (exit %d): %s\n--- last lines of %s ---\n%s'
                             % (stage, r.returncode, ' '.join(map(str, cmd)), logf, tail))

    def say(self, msg: str):
        print(msg, flush=True)
        if self.log:
            self.log.write(msg + '\n'); self.log.flush()

    # DONE markers -------------------------------------------------------------
    def done_path(self, stage: str) -> str:
        return os.path.join(self.run, stage, 'DONE')

    def is_done(self, stage: str, inputs: list, params: dict) -> bool:
        p = self.done_path(stage)
        if self.force or not os.path.isfile(p):
            return False
        try:
            rec = json.load(open(p))
        except Exception:
            return False
        if rec.get('params_hash') != self.params_hash(params):
            return False
        for f, fp in rec.get('inputs', {}).items():
            if not os.path.exists(f) or fingerprint(f) != fp:
                return False
        return set(rec.get('inputs', {})) == set(inputs)

    def mark_done(self, stage: str, inputs: list, params: dict, outputs: list):
        rec = {'stage': stage, 'finished': dt.datetime.now().isoformat(timespec='seconds'),
               'params_hash': self.params_hash(params), 'params': params,
               'inputs': {f: fingerprint(f) for f in inputs},
               'outputs': {f: fingerprint(f) for f in outputs if os.path.isfile(f)}}
        os.makedirs(os.path.dirname(self.done_path(stage)), exist_ok=True)
        json.dump(rec, open(self.done_path(stage), 'w'), indent=1)
        self.manifest.setdefault('stages', {})[stage] = rec

    @staticmethod
    def params_hash(params: dict) -> str:
        return hashlib.sha256(json.dumps(params, sort_keys=True).encode()).hexdigest()[:16]


# ----------------------------------------------------------------------------- stage registry

STAGES: list = []


def stage(sid: str, needs: tuple = ()):
    def deco(fn):
        fn.sid = sid
        fn.needs = needs
        STAGES.append(fn)
        return fn
    return deco


def run_stage(fn, ctx: Ctx):
    sid = fn.sid
    inputs, params = fn.plan(ctx) if hasattr(fn, 'plan') else ([], {})
    if ctx.is_done(sid, inputs, params):
        ctx.say('[%s] up to date, skipped' % sid)
        return
    ctx.say('[%s] running' % sid)
    outputs = fn(ctx) or []
    ctx.mark_done(sid, inputs, params, outputs)
    ctx.say('[%s] done' % sid)


# ----------------------------------------------------------------------------- tables

@stage('00_tables')
def st_tables(ctx: Ctx):
    t = ctx.d('tables')
    with open(os.path.join(t, 'genomes.tsv'), 'w') as f:
        for code, g in ctx.cfg.genomes.items():
            f.write('%s\t%s\n' % (code, g))
    with open(os.path.join(t, 'families.tsv'), 'w') as f:
        for fam, fd in ctx.cfg.families.items():
            f.write('%s\t%d\t%s\n' % (fam, fd['mono'], fd['consensus']))
    sub = ctx.cfg.data.get('subfamily') or {}
    with open(os.path.join(t, 'members.tsv'), 'w') as f:
        for fam, members in (sub.get('members') or {}).items():
            for m in members:
                f.write('%s\t%s\n' % (fam, m))
    # Step-8 references: the phase-rotated consensuses
    rot = (ctx.cfg.data.get('references') or {}).get('rotated_consensus')
    refs = os.path.join(t, 'refs.rot')
    if rot:
        shutil.copyfile(rot, refs)
        ctx.manifest['refs_rot'] = {'source': rot, 'derived': False}
    else:
        with open(refs, 'w') as f:
            for fam, fd in ctx.cfg.families.items():
                seq = ''.join(l.strip() for l in open(fd['consensus']) if not l.startswith('>')).upper()
                f.write('>%s\n%s\n' % (fam, minrot(seq)))
        ctx.manifest['refs_rot'] = {'source': None, 'derived': True,
                                    'note': 'Step-5 (lexicographically minimal) rotation of each family consensus; '
                                            'the published Darevskia run used the file CLsat_cons.rot'}
    # genome indexes
    for code, g in ctx.cfg.genomes.items():
        if not os.path.isfile(g + '.fai'):
            ctx.sh([ctx.tool('samtools'), 'faidx', g], '00_tables')
    return [os.path.join(t, n) for n in ('genomes.tsv', 'families.tsv', 'members.tsv', 'refs.rot')]


def _plan_tables(ctx: Ctx):
    ins = [fd['consensus'] for fd in ctx.cfg.families.values()]
    rot = (ctx.cfg.data.get('references') or {}).get('rotated_consensus')
    if rot:
        ins.append(rot)
    return ins, {'genomes': ctx.cfg.genomes, 'families': ctx.cfg.families, 'members': (ctx.cfg.data.get('subfamily') or {}).get('members')}


st_tables.plan = _plan_tables


# ----------------------------------------------------------------------------- 02 TRF

def _large_threshold(ctx: Ctx) -> int:
    return int(ctx.cfg.data['finalisation']['large_scaffold_bp'])


def _scaffold_lengths(genome: str):
    for line in open(genome + '.fai'):
        p = line.split('\t')
        yield p[0], int(p[1])


@stage('02_trf', needs=('trf', 'parallel'))
def st_trf(ctx: Ctx):
    trf_dir = ctx.d('trf')
    disc = ctx.cfg.data['discovery']
    params = disc['trf']
    given = (ctx.cfg.data.get('inputs') or {}).get('trf_dat') or {}
    big_rows = []
    outs = []
    for code, genome in ctx.cfg.genomes.items():
        dat = os.path.join(trf_dir, code + '.dat')
        thr = _large_threshold(ctx)
        large = [(scf, ln) for scf, ln in _scaffold_lengths(genome) if ln > thr]
        for scf, ln in large:
            big_rows.append((code, scf, ln))
        if code in given:
            if os.path.islink(dat) or os.path.exists(dat):
                os.remove(dat)
            os.symlink(os.path.abspath(given[code]), dat)
            ctx.say('  %s: using given TRF output %s' % (code, given[code]))
            outs.append(dat)
            continue
        # run TRF per scaffold in parallel, scaffolds above the threshold excluded (stage 02b handles them)
        split = ctx.d('trf', code + '.split')
        names = []
        with open(genome) as fh:
            out = None
            skip = False
            for line in fh:
                if line.startswith('>'):
                    name = line[1:].split()[0]
                    skip = any(name == s for s, _ in large)
                    if out:
                        out.close()
                    out = None
                    if not skip:
                        names.append(name)
                        out = open(os.path.join(split, name + '.fa'), 'w')
                        out.write('>%s\n' % name)
                elif out and not skip:
                    out.write(line)
            if out:
                out.close()
        joblist = os.path.join(split, 'jobs.txt')
        with open(joblist, 'w') as f:
            for name in names:
                f.write(name + '.fa\n')
        suffix = '.' + params.replace(' ', '.') + '.dat'
        jobs = int(ctx.cfg.data['expansion']['jobs'])
        ctx.sh(['bash', '-c', 'cat jobs.txt | %s --jobs %d "%s {} %s -h > /dev/null 2>&1 || true"'
                % (ctx.tool('parallel'), jobs, ctx.tool('trf'), params)], '02_trf', cwd=split, log_name='02_trf_' + code)
        with open(dat, 'w') as out:
            missing = []
            for name in names:
                part = os.path.join(split, name + '.fa' + suffix)
                if not os.path.isfile(part):
                    missing.append(name)
                    continue
                with open(part) as fh:
                    out.write(fh.read())
            if missing:
                raise StageError('02_trf %s: TRF produced no output for %d scaffolds, e.g. %s' % (code, len(missing), missing[:5]))
        shutil.rmtree(split, ignore_errors=True)
        ctx.say('  %s: TRF on %d scaffolds -> %s' % (code, len(names), dat))
        outs.append(dat)
    with open(os.path.join(trf_dir, 'large_scaffolds.tsv'), 'w') as f:
        f.write('# scaffolds longer than finalisation.large_scaffold_bp=%d; scanned in chunks by stage 02b\n' % _large_threshold(ctx))
        for code, scf, ln in big_rows:
            f.write('%s\t%s\t%d\n' % (code, scf, ln))
    return outs + [os.path.join(trf_dir, 'large_scaffolds.tsv')]


def _plan_trf(ctx: Ctx):
    given = (ctx.cfg.data.get('inputs') or {}).get('trf_dat') or {}
    ins = list(given.values()) + [g + '.fai' for g in ctx.cfg.genomes.values()]
    return ins, {'trf': ctx.cfg.data['discovery']['trf'], 'large': _large_threshold(ctx), 'genomes': ctx.cfg.genomes}


st_trf.plan = _plan_trf


@stage('02b_large', needs=('trf', 'parallel', 'samtools'))
def st_large(ctx: Ctx):
    """Chunked TRF of the scaffolds stage 02 excluded, classified against the family consensuses;
    the kept candidates become extra seeds (published: the 85 Mb unp scaffold -> locus 540)."""
    rows = read_tsv(os.path.join(ctx.run, 'trf', 'large_scaffolds.tsv'), 3)
    fin = ctx.cfg.data['finalisation']
    disc = ctx.cfg.data['discovery']
    outs = []
    seeds_extra = {fam: [] for fam in ctx.families}
    for code, scf, ln in rows:
        d = ctx.d('large', code, scf)
        env = ctx.env(CHUNK_BP=fin['chunk_bp'], OVERLAP_BP=fin['overlap_bp'], TRF_PARAMS=disc['trf'],
                      TRF_TIMEOUT=fin['trf_timeout'], MIN_PIECE_BP=fin['min_piece_bp'],
                      JOBS=ctx.cfg.data['expansion']['jobs'])
        ctx.sh(['bash', os.path.join(ctx.scripts, '06_catalogue_finalize/85mb_recovery/run_chunked.sh'),
                ctx.cfg.genomes[code], scf, d], '02b_large', env=env, log_name='02b_large_%s_%s' % (code, scf))
        cand = os.path.join(d, 'candidates.tsv')
        ctx.sh(['python3', os.path.join(ctx.scripts, '06_catalogue_finalize/85mb_recovery/classify_candidates.py'),
                os.path.join(d, 'scaf.trf.dedup.dat'), os.path.join(ctx.run, 'tables', 'families.tsv'), cand,
                '--min-copies', str(disc['min_copies'])], '02b_large', log_name='02b_large_%s_%s' % (code, scf))
        for p in read_tsv(cand, 6):
            if p[0] == 'start':
                continue
            s, e, fam = int(p[0]), int(p[1]), p[4]
            seeds_extra[fam].append('%s\t%d\t%d\t%s|%s:%d-%d|RECOVERED\t0\t.\n' % (scf, s - 1, e, code, scf, s, e))
        outs.append(cand)
        uns = os.path.join(d, 'unscanned.tsv')
        if os.path.isfile(uns) and os.path.getsize(uns):
            ctx.manifest.setdefault('unscanned', []).extend(
                {'code': code, 'scaffold': scf, 'piece': p[0], 'start': int(p[1]), 'end': int(p[2]), 'note': p[4] if len(p) > 4 else ''}
                for p in read_tsv(uns, 3))
            ctx.say('  %s %s: %d piece(s) TRF could not finish (listed in %s)' % (code, scf, sum(1 for _ in read_tsv(uns, 3)), uns))
    for fam, lines in seeds_extra.items():
        p = os.path.join(ctx.d('large'), fam + '.recovered.bed')
        with open(p, 'w') as f:
            f.writelines(lines)
        outs.append(p)
    ctx.manifest['large_scaffolds'] = [{'code': c, 'scaffold': s, 'length': int(l)} for c, s, l in rows]
    return outs


def _plan_large(ctx: Ctx):
    return [os.path.join(ctx.run, 'trf', 'large_scaffolds.tsv'), os.path.join(ctx.run, 'tables', 'families.tsv')], \
        {'fin': ctx.cfg.data['finalisation'], 'trf': ctx.cfg.data['discovery']['trf'], 'min_copies': ctx.cfg.data['discovery']['min_copies']}


st_large.plan = _plan_large


# ----------------------------------------------------------------------------- 03-07 discovery

def _disc_env(ctx: Ctx):
    d = ctx.cfg.data['discovery']
    return ctx.env(PERIOD_MIN=d['period'][0], PERIOD_MAX=d['period'][1], MIN_COPIES=d['min_copies'],
                   TRF_PARAMS=d['trf'], SIMPLE_REPEAT_COV=d['simple_repeat_cov'], STRAND=d['strand'],
                   CLUSTER_ID=d['cluster_id'], CROSS_CLUSTER_ID=d['cross_cluster_id'])


def _dats(ctx):
    return [os.path.join(ctx.run, 'trf', c + '.dat') for c in ctx.codes]


@stage('03_motifs', needs=('gawk',))
def st_motifs(ctx: Ctx):
    m = ctx.d('motifs')
    ctx.sh(['bash', os.path.join(ctx.scripts, '01_discovery_clustering/Step3.sh'), os.path.join(ctx.run, 'trf'), m],
           '03_motifs', env=_disc_env(ctx))
    return [os.path.join(m, c + '.fa') for c in ctx.codes]


st_motifs.plan = lambda ctx: (_dats(ctx), {'period': ctx.cfg.data['discovery']['period'], 'min_copies': ctx.cfg.data['discovery']['min_copies']})


@stage('04_filter', needs=('trf', 'gawk'))
def st_filter(ctx: Ctx):
    m = ctx.d('motifs')
    for c in ctx.codes:
        ctx.sh(['bash', os.path.join(ctx.scripts, '01_discovery_clustering/Step4.sh'), c + '.fa'], '04_filter', cwd=m, env=_disc_env(ctx))
    return [os.path.join(m, c + '.fa.clean.fas') for c in ctx.codes]


st_filter.plan = lambda ctx: ([os.path.join(ctx.run, 'motifs', c + '.fa') for c in ctx.codes],
                              {'trf': ctx.cfg.data['discovery']['trf'], 'cov': ctx.cfg.data['discovery']['simple_repeat_cov']})


@stage('05_rotate')
def st_rotate(ctx: Ctx):
    m = ctx.d('motifs')
    ctx.sh(['bash', os.path.join(ctx.scripts, '01_discovery_clustering/Step5.sh'), m], '05_rotate', env=_disc_env(ctx))
    return [os.path.join(m, 'clean', c + '.fa.clean.fas.rot') for c in ctx.codes]


st_rotate.plan = lambda ctx: ([os.path.join(ctx.run, 'motifs', c + '.fa.clean.fas') for c in ctx.codes],
                              {'strand': ctx.cfg.data['discovery']['strand']})


@stage('06_cluster', needs=('vsearch',))
def st_cluster(ctx: Ctx):
    clean = os.path.join(ctx.run, 'motifs', 'clean')
    ctx.sh(['bash', os.path.join(ctx.scripts, '01_discovery_clustering/Step6.sh'), clean, os.path.join(clean, 'clusters_by_file')],
           '06_cluster', env=_disc_env(ctx))
    return []


st_cluster.plan = lambda ctx: ([os.path.join(ctx.run, 'motifs', 'clean', c + '.fa.clean.fas.rot') for c in ctx.codes],
                               {'id': ctx.cfg.data['discovery']['cluster_id'], 'strand': ctx.cfg.data['discovery']['strand']})


@stage('07_mafft', needs=('mafft',))
def st_mafft(ctx: Ctx):
    cb = os.path.join(ctx.run, 'motifs', 'clean', 'clusters_by_file')
    dirs = sorted(os.path.join(cb, x) for x in os.listdir(cb) if os.path.isdir(os.path.join(cb, x)))
    script = os.path.join(ctx.scripts, '01_discovery_clustering/Step7.sh')
    jobs = min(len(dirs), int(ctx.cfg.data['expansion']['jobs']))

    def one(d):
        ctx.sh(['bash', script], '07_mafft', cwd=d, log_name='07_mafft_' + os.path.basename(d))

    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, jobs)) as ex:
        list(ex.map(one, dirs))
    return []


st_mafft.plan = lambda ctx: ([os.path.join(ctx.run, 'motifs', 'clean', c + '.fa.clean.fas.rot') for c in ctx.codes], {})


# ----------------------------------------------------------------------------- 08-10 seeds

@stage('08_cross', needs=('vsearch',))
def st_cross(ctx: Ctx):
    x = ctx.d('cross')
    ctx.sh(['bash', os.path.join(ctx.scripts, '02_seed_expansion_scan/Step8_cross_cluster.sh'),
            os.path.join(ctx.run, 'motifs', 'clean', 'clusters_by_file'), os.path.join(ctx.run, 'tables', 'refs.rot'), x],
           '08_cross', env=_disc_env(ctx))
    cent = os.path.join(x, 'cross_species_centroids.fa')
    ctx.manifest['centroids_md5'] = md5sum(cent)
    ctx.manifest['centroids_n'] = sum(1 for _ in read_fasta_names(cent))
    ctx.say('  centroids: %d, md5 %s' % (ctx.manifest['centroids_n'], ctx.manifest['centroids_md5']))
    return [cent]


st_cross.plan = lambda ctx: ([os.path.join(ctx.run, 'tables', 'refs.rot')] + [os.path.join(ctx.run, 'motifs', 'clean', c + '.fa.clean.fas.rot') for c in ctx.codes],
                             {'id': ctx.cfg.data['discovery']['cross_cluster_id'], 'strand': ctx.cfg.data['discovery']['strand']})


@stage('09_assign', needs=('ssearch36',))
def st_assign(ctx: Ctx):
    """Validate the curated cluster -> family table against the Step-8 clusters, write a report
    of the unassigned clusters for the curator, then build the member FASTA per family."""
    x = os.path.join(ctx.run, 'cross')
    cc = os.path.join(x, 'cross_clusters')
    table = ctx.cfg.data['family_assignment']
    clusters = {}
    for fn in os.listdir(cc):
        if fn.startswith('clust_'):
            clusters[fn[6:]] = os.path.join(cc, fn)
    assigned = {}
    errors = []
    for p in read_tsv(table, 2):
        clu, fam = p[0].strip(), p[1].strip()
        if clu not in clusters:
            errors.append('cluster %s (family %s) does not exist in %s' % (clu, fam, cc))
        if fam not in ctx.cfg.families:
            errors.append('family %s is not in the config' % fam)
        assigned[clu] = fam
    for fam in ctx.families:
        if fam not in assigned.values():
            errors.append('family %s has no cluster in %s' % (fam, table))
    # report: every cluster, size, assigned family, best consensus hit of its centroid
    refs = {n: None for n in read_fasta_names(os.path.join(ctx.run, 'tables', 'refs.rot'))}
    cent_seq = {}
    name = None
    for line in open(os.path.join(x, 'cross_species_centroids.fa')):
        if line.startswith('>'):
            name = line[1:].split()[0]; cent_seq[name] = ''
        else:
            cent_seq[name] += line.strip()
    # which centroid belongs to which cluster: from the .uc file (S records)
    cent_of = {}
    for line in open(os.path.join(x, 'cross_species.uc')):
        p = line.rstrip('\n').split('\t')
        if p[0] == 'S':
            cent_of[p[1]] = p[8]
    sizes = {clu: sum(1 for _ in read_fasta_names(path)) for clu, path in clusters.items()}
    # best family bits per centroid: one ssearch36 run of all centroids vs doubled consensuses
    q = os.path.join(x, 'centroids_query.fa')
    lib = os.path.join(x, 'consensus_x2.fa')
    with open(q, 'w') as f:
        for clu, c in cent_of.items():
            if c in cent_seq:
                f.write('>%s\n%s\n' % (clu, cent_seq[c]))
    with open(lib, 'w') as f:
        for fam, fd in ctx.cfg.families.items():
            s = ''.join(l.strip() for l in open(fd['consensus']) if not l.startswith('>')).upper()
            f.write('>%s\n%s%s\n' % (fam, s, s))
    r = subprocess.run([ctx.tool('ssearch36'), '-n', '-m', '8C', '-E', '1000', '-b', '1', '-d', '1', q, lib],
                       capture_output=True, text=True, env=ctx.env())
    best = {}
    for line in r.stdout.splitlines():
        if line.startswith('#'):
            continue
        p = line.split('\t')
        if len(p) >= 12:
            try:
                bits = float(p[11])
            except ValueError:
                continue
            if p[0] not in best or bits > best[p[0]][1]:
                best[p[0]] = (p[1], bits, p[2])
    rep = os.path.join(x, 'cluster_report.tsv')
    with open(rep, 'w') as f:
        f.write('cluster\tsize\tassigned_family\tcentroid\tis_reference\tbest_family\tbest_bits\tbest_pid\n')
        for clu in sorted(clusters, key=int):
            c = cent_of.get(clu, '')
            b = best.get(clu, ('', 0.0, ''))
            f.write('%s\t%d\t%s\t%s\t%s\t%s\t%.1f\t%s\n' % (clu, sizes[clu], assigned.get(clu, ''), c,
                                                        'yes' if c in refs else 'no', b[0], b[1], b[2]))
    ctx.say('  cluster report: %s (%d clusters, %d assigned)' % (rep, len(clusters), len(assigned)))
    if errors:
        raise StageError('09_assign: the curated table does not match the clusters:\n  ' + '\n  '.join(errors)
                         + '\n  see %s' % rep)
    s = ctx.d('seeds')
    ctx.sh(['bash', os.path.join(ctx.scripts, '02_seed_expansion_scan/Step9_family_fasta.sh'), cc, table, s,
            os.path.join(ctx.run, 'tables', 'refs.rot')], '09_assign')
    return [os.path.join(s, fam + '.members.fa') for fam in ctx.families] + [rep]


st_assign.plan = lambda ctx: ([ctx.cfg.data['family_assignment'], os.path.join(ctx.run, 'cross', 'cross_species_centroids.fa')], {})


@stage('10_seeds', needs=('gawk',))
def st_seeds(ctx: Ctx):
    s = ctx.d('seeds')
    given = (ctx.cfg.data.get('inputs') or {}).get('seeds') or {}
    outs = []
    counts = {}
    for fam in ctx.families:
        bed = os.path.join(s, fam + '.bed')
        if fam in given:
            shutil.copyfile(given[fam], bed)
        else:
            with open(bed, 'w') as f:
                subprocess.run(['bash', os.path.join(ctx.scripts, '02_seed_expansion_scan/make_bed_from_clsats_headers.sh'),
                                os.path.join(s, fam + '.members.fa')], stdout=f, env=ctx.env(), check=True)
        rec = os.path.join(ctx.run, 'large', fam + '.recovered.bed')
        if os.path.isfile(rec) and os.path.getsize(rec):
            with open(bed, 'a') as f, open(rec) as r:
                f.write(r.read())
        counts[fam] = sum(1 for _ in open(bed))
        outs.append(bed)
    ctx.manifest['seed_counts'] = counts
    ctx.say('  seeds: ' + ', '.join('%s %d' % kv for kv in counts.items()))
    return outs


st_seeds.plan = lambda ctx: ([os.path.join(ctx.run, 'seeds', fam + '.members.fa') for fam in ctx.families
                              if fam not in ((ctx.cfg.data.get('inputs') or {}).get('seeds') or {})]
                             + list(((ctx.cfg.data.get('inputs') or {}).get('seeds') or {}).values())
                             + [p for p in (os.path.join(ctx.run, 'large', f + '.recovered.bed') for f in ctx.families) if os.path.isfile(p)], {})


# ----------------------------------------------------------------------------- 11-14 catalogue

def _expand_cmd(ctx: Ctx, bed: str, fam: str, out_prefix: str, scan_only: bool) -> list:
    e = ctx.cfg.data['expansion']
    sc = ctx.cfg.data['scan']
    fd = ctx.cfg.families[fam]
    cmd = ['bash', os.path.join(ctx.scripts, '02_seed_expansion_scan/clisat_expand_then_scan_merged.sh'),
           bed, os.path.join(ctx.run, 'tables', 'genomes.tsv'), fd['consensus'], str(fd['mono']), out_prefix,
           '--step-monomers', e['step_mono'], '--win-monomers', e['win_mono'], '--query-repeat', e['query_repeat'],
           '--min-aln-frac', e['min_aln_frac'], '--min-bits', e['min_bits'], '--max-steps', e['max_steps'],
           '--probe-monomers', e['probe_mono'], '--probe-min-step', e['probe_min_step'],
           '--scan-win-bp', sc['win_mono'] * fd['mono'], '--scan-step-bp', sc['step_mono'] * fd['mono'],
           '--scan-min-bits', sc['min_bits'], '--jobs', e['jobs']]
    if e.get('to_scaffold_ends'):
        cmd.append('--to-scaffold-ends')
    if scan_only:
        cmd.append('--scan-only')
    return cmd


@stage('11_expand', needs=('ssearch36', 'samtools', 'bedtools', 'parallel'))
def st_expand(ctx: Ctx):
    x = ctx.d('expand')
    outs = []
    with open(os.path.join(x, 'expansions.tsv'), 'w') as t:
        for fam in ctx.families:
            pre = os.path.join(x, 'out_' + fam)
            ctx.sh(_expand_cmd(ctx, os.path.join(ctx.run, 'seeds', fam + '.bed'), fam, pre, False), '11_expand', log_name='11_expand_' + fam)
            t.write('%s\t%s\n' % (fam, pre + '.expanded_merged.bed'))
            outs.append(pre + '.expanded_merged.bed')
    return outs + [os.path.join(x, 'expansions.tsv')]


st_expand.plan = lambda ctx: ([os.path.join(ctx.run, 'seeds', f + '.bed') for f in ctx.families] + [fd['consensus'] for fd in ctx.cfg.families.values()],
                              {'expansion': ctx.cfg.data['expansion'], 'scan': ctx.cfg.data['scan'], 'mono': {f: fd['mono'] for f, fd in ctx.cfg.families.items()}})


@stage('12_master', needs=('bedtools',))
def st_master(ctx: Ctx):
    m = ctx.d('master')
    ctx.sh(['bash', os.path.join(ctx.scripts, '03_master_arbitration/build_master.sh'), os.path.join(ctx.run, 'expand', 'expansions.tsv'), m], '12_master')
    n = sum(1 for _ in open(os.path.join(m, 'MASTER.merged.bed')))
    ctx.manifest['master_loci'] = n
    ctx.say('  MASTER loci: %d' % n)
    return [os.path.join(m, 'MASTER.merged.bed')]


st_master.plan = lambda ctx: ([os.path.join(ctx.run, 'expand', 'out_' + f + '.expanded_merged.bed') for f in ctx.families], {})


@stage('13_rescan', needs=('ssearch36', 'samtools', 'bedtools', 'parallel'))
def st_rescan(ctx: Ctx):
    r = ctx.d('rescan')
    master = os.path.join(ctx.run, 'master', 'MASTER.merged.bed')
    outs = []
    with open(os.path.join(r, 'arbitration_inputs.tsv'), 'w') as t:
        for fam in ctx.families:
            pre = os.path.join(r, 'MASTER_' + fam)
            ctx.sh(_expand_cmd(ctx, master, fam, pre, True), '13_rescan', log_name='13_rescan_' + fam)
            t.write('%s\t%d\t%s\t%s\n' % (fam, ctx.cfg.families[fam]['mono'], pre + '.bits.tsv', os.path.join(ctx.run, 'seeds', fam + '.bed')))
            outs.append(pre + '.bits.tsv')
    return outs + [os.path.join(r, 'arbitration_inputs.tsv')]


st_rescan.plan = lambda ctx: ([os.path.join(ctx.run, 'master', 'MASTER.merged.bed')] + [fd['consensus'] for fd in ctx.cfg.families.values()],
                              {'scan': ctx.cfg.data['scan'], 'jobs': ctx.cfg.data['expansion']['jobs']})


@stage('14_arbitrate', needs=('bedtools', 'gawk', 'python3'))
def st_arbitrate(ctx: Ctx):
    a = ctx.cfg.data['arbitration']
    out = ctx.d('arbitration')
    env = ctx.env(DELTA=a['delta'], DELTA_AMBIG=a['delta_ambig'], SEG_BIN=a['seg_bin'], SEG_MIN_RUN=a['seg_min_run'],
                  SEG_MIN_LOCUS=a['seg_min_locus'], SEG_SHORT=a['short_segments'])
    ctx.sh(['bash', os.path.join(ctx.scripts, '03_master_arbitration/arbitration_patched.sh'),
            os.path.join(ctx.run, 'rescan', 'arbitration_inputs.tsv'), os.path.join(ctx.run, 'master', 'MASTER.merged.bed'), out],
           '14_arbitrate', env=env)
    n = sum(1 for _ in open(os.path.join(out, 'MASTER.loci.segmented.bed')))
    ctx.manifest['segmented_loci'] = n
    ctx.say('  segmented loci: %d' % n)
    return [os.path.join(out, f) for f in ('MASTER.loci.segmented.bed', 'MASTER.final_assignment.tsv', 'MASTER.locus_bits.tsv')]


st_arbitrate.plan = lambda ctx: ([os.path.join(ctx.run, 'rescan', 'arbitration_inputs.tsv'), os.path.join(ctx.run, 'master', 'MASTER.merged.bed')]
                                 + [os.path.join(ctx.run, 'rescan', 'MASTER_' + f + '.bits.tsv') for f in ctx.families],
                                 {'arbitration': ctx.cfg.data['arbitration']})


@stage('15_catalogue', needs=('ssearch36', 'samtools', 'bedtools', 'parallel'))
def st_catalogue(ctx: Ctx):
    """Final catalogue = arbitration output + manual loci (tagged MANUAL, numbered after the
    pipeline loci, scanned with every consensus for their bits)."""
    cat = ctx.d('catalogue')
    arb = os.path.join(ctx.run, 'arbitration')
    for f in ('MASTER.loci.segmented.bed', 'MASTER.final_assignment.tsv', 'MASTER.locus_bits.tsv', 'MASTER.loci.bed'):
        shutil.copyfile(os.path.join(arb, f), os.path.join(cat, f))
    shutil.copyfile(os.path.join(ctx.run, 'master', 'MASTER.merged.bed'), os.path.join(cat, 'MASTER.merged.bed'))
    for fam in ctx.families:
        shutil.copyfile(os.path.join(ctx.run, 'rescan', 'MASTER_%s.bits.tsv' % fam), os.path.join(cat, 'MASTER_%s.bits.tsv' % fam))
    manual = ctx.cfg.data.get('manual_loci')
    added = []
    if manual and os.path.getsize(manual):
        n_master = sum(1 for _ in open(os.path.join(cat, 'MASTER.loci.bed')))
        rows = [p for p in read_tsv(manual, 3) if not p[0].startswith(('track', 'browser'))]
        mbed = os.path.join(cat, 'manual.bed')
        with open(mbed, 'w') as f:
            for p in rows:
                f.write('\t'.join(p[:4] if len(p) >= 4 else p[:3] + ['MANUAL']) + '\n')
        # bits of the manual loci against every family (scan-only), max per locus
        maxbits = {i: {} for i in range(len(rows))}
        for fam in ctx.families:
            pre = os.path.join(cat, 'manual_' + fam)
            ctx.sh(_expand_cmd(ctx, mbed, fam, pre, True), '15_catalogue', log_name='15_catalogue_' + fam)
            with open(pre + '.bits.tsv') as fh, open(os.path.join(cat, 'MASTER_%s.bits.tsv' % fam), 'a') as out:
                fh.readline()
                for line in fh:
                    out.write(line)
                    p = line.rstrip('\n').split('\t')
                    for i, r in enumerate(rows):
                        if p[0] == r[0] and int(r[1]) < int(p[1]) <= int(r[2]):
                            maxbits[i][fam] = max(maxbits[i].get(fam, 0.0), float(p[2]))
        with open(os.path.join(cat, 'MASTER.merged.bed'), 'a') as mb, open(os.path.join(cat, 'MASTER.loci.bed'), 'a') as lb, \
             open(os.path.join(cat, 'MASTER.loci.segmented.bed'), 'a') as sb, open(os.path.join(cat, 'MASTER.locus_bits.tsv'), 'a') as bb, \
             open(os.path.join(cat, 'MASTER.final_assignment.tsv'), 'a') as fa:
            for i, p in enumerate(rows):
                lid = 'L%d' % (n_master + i + 1)
                name = p[3] if len(p) >= 4 else 'MANUAL'
                fam = p[6] if len(p) >= 7 and p[6] in ctx.cfg.families else max(maxbits[i], key=maxbits[i].get, default=ctx.families[0])
                bits = maxbits[i]
                others = sorted(((b, f) for f, b in bits.items() if f != fam), reverse=True)
                best_other, best_bits = (others[0][1], others[0][0]) if others else ('', 0.0)
                mb.write('%s\t%s\t%s\t%s\tMANUAL\n' % (p[0], p[1], p[2], name))
                lb.write('%s\t%s\t%s\t%s\n' % (p[0], p[1], p[2], lid))
                sb.write('%s\t%s\t%s\t%s\n' % (p[0], p[1], p[2], lid))
                bb.write(lid + ''.join('\t%s' % ('%g' % bits.get(f, 0.0)) for f in ctx.families) + '\n')
                fa.write('%s\t%s\t%g\t%s\t%g\t%g\tMANUAL\tMANUAL\n' % (lid, fam, bits.get(fam, 0.0), best_other, best_bits, bits.get(fam, 0.0) - best_bits))
                added.append({'locus': lid, 'bed': p[:4], 'family': fam, 'bits': bits})
    ctx.manifest['manual_loci'] = added
    ctx.manifest['catalogue_loci'] = sum(1 for _ in open(os.path.join(cat, 'MASTER.loci.segmented.bed')))
    ctx.say('  catalogue: %d segmented loci (%d manual)' % (ctx.manifest['catalogue_loci'], len(added)))
    # a header-less copy of the assignment for the viewer (its name in the viewer repo)
    shutil.copyfile(os.path.join(cat, 'MASTER.final_assignment.tsv'), os.path.join(cat, 'MASTER.final_assignment.clean.tsv'))
    return [os.path.join(cat, f) for f in ('MASTER.merged.bed', 'MASTER.loci.segmented.bed', 'MASTER.final_assignment.tsv', 'MASTER.locus_bits.tsv')]


st_catalogue.plan = lambda ctx: ([os.path.join(ctx.run, 'arbitration', 'MASTER.final_assignment.tsv')] + ([ctx.cfg.data['manual_loci']] if ctx.cfg.data.get('manual_loci') else []), {})


# ----------------------------------------------------------------------------- 16-19 regions, terminal, subfamily, orthology

@stage('16_regions', needs=('samtools',))
def st_regions(ctx: Ctx):
    reg = ctx.d('regions')
    cat = os.path.join(ctx.run, 'catalogue')
    ctx.sh(['bash', os.path.join(ctx.scripts, '05_subfamily_terminal_boundary/build_locus_master.sh'),
            os.path.join(cat, 'MASTER.merged.bed'), cat, reg, os.path.join(ctx.run, 'tables', 'families.tsv')], '16_regions')
    ctx.sh(['bash', os.path.join(ctx.scripts, '06_catalogue_finalize/subfamily_rerun/run_extract.sh'), reg,
            os.path.join(ctx.run, 'tables', 'genomes.tsv')], '16_regions')
    return [os.path.join(reg, f) for f in ('locus_master.tsv', 'manifest.tsv', 'all_regions.fa')]


st_regions.plan = lambda ctx: ([os.path.join(ctx.run, 'catalogue', f) for f in ('MASTER.merged.bed', 'MASTER.loci.segmented.bed', 'MASTER.final_assignment.tsv')], {})


@stage('17_terminal', needs=('ssearch36', 'samtools'))
def st_terminal(ctx: Ctx):
    reg = os.path.join(ctx.run, 'regions')
    t = ctx.cfg.data['terminal']
    env = ctx.env(WINDOW_MONO=t['window_mono'], EDGE_TOL=t['edge_tol'])
    ctx.sh(['bash', os.path.join(ctx.scripts, '05_subfamily_terminal_boundary/terminal_monomers.sh'), reg,
            os.path.join(ctx.run, 'tables', 'genomes.tsv'), os.path.join(ctx.run, 'tables', 'families.tsv')], '17_terminal', env=env)
    ctx.sh(['bash', os.path.join(ctx.scripts, '05_subfamily_terminal_boundary/boundary_phase.sh'), reg,
            os.path.join(ctx.run, 'tables', 'families.tsv')], '17_terminal')
    ctx.sh(['python3', os.path.join(ctx.scripts, '05_subfamily_terminal_boundary/build_json.py'), reg], '17_terminal')
    return [os.path.join(reg, f) for f in ('terminal_monomers.tsv', 'locus_boundary_phase.tsv', 'boundary_phase.json')]


st_terminal.plan = lambda ctx: ([os.path.join(ctx.run, 'regions', f) for f in ('locus_master.tsv', 'manifest.tsv')], {'terminal': ctx.cfg.data['terminal']})


@stage('18_subfamily', needs=('ssearch36', 'samtools'))
def st_subfamily(ctx: Ctx):
    reg = os.path.join(ctx.run, 'regions')
    sub = ctx.cfg.data.get('subfamily') or {}
    cons = sub.get('consensuses')
    if cons and os.path.isdir(cons):
        ctx.sh(['bash', os.path.join(ctx.scripts, '05_subfamily_terminal_boundary/assign_subfamily.sh'), reg, cons,
                os.path.join(ctx.run, 'tables', 'members.tsv')], '18_subfamily')
    else:
        ctx.say('  no sub-subfamily consensus directory: writing an empty locus_subfamily.tsv')
        with open(os.path.join(reg, 'locus_subfamily.tsv'), 'w') as f:
            f.write('locus\tmajor\tsubfamily\tpct_id\taln_len\tbits\tsecond\tsecond_bits\tdelta\n')
    ctx.sh(['python3', os.path.join(ctx.scripts, '05_subfamily_terminal_boundary/build_subfam_term_json.py'), reg,
            os.path.join(ctx.run, 'tables', 'families.tsv')], '18_subfamily', env=ctx.env(CONFIDENT_DELTA=sub.get('confident_delta', 10)))
    return [os.path.join(reg, f) for f in ('locus_subfamily.tsv', 'subfamily.json', 'terminal_monomers.json')]


st_subfamily.plan = lambda ctx: ([os.path.join(ctx.run, 'regions', f) for f in ('locus_master.tsv', 'terminal_monomers.tsv')],
                                 {'subfamily': ctx.cfg.data.get('subfamily')})


@stage('19_orthology', needs=('ssearch36', 'samtools'))
def st_orthology(ctx: Ctx):
    o = ctx.cfg.data['orthology']
    out = ctx.d('ortho')
    ctx.sh(['python3', os.path.join(ctx.scripts, '04_flank_orthology/rebuild_ortho_segmented.py'),
            '--loci', os.path.join(ctx.run, 'regions', 'locus_master.tsv'), '--genomes', os.path.join(ctx.run, 'tables', 'genomes.tsv'),
            '--families', os.path.join(ctx.run, 'tables', 'families.tsv'), '--out', out,
            '--flank', o['flank_bp'], '--min-aln', o['min_aln'], '--min-pid', o['min_pid'], '--guard', o['guard_bp']], '19_orthology')
    pairs = sum(1 for _ in open(os.path.join(out, 'ortho_pairs.tsv'))) - 1
    ctx.manifest['ortholog_pairs'] = pairs
    ctx.say('  ortholog locus pairs: %d' % pairs)
    return [os.path.join(out, f) for f in ('flanks.json', 'flank_orthologs.json', 'ortho_pairs.tsv')]


st_orthology.plan = lambda ctx: ([os.path.join(ctx.run, 'regions', 'locus_master.tsv')], {'orthology': ctx.cfg.data['orthology']})


@stage('20_candidates', needs=('ssearch36', 'samtools', 'bedtools'))
def st_candidates(ctx: Ctx):
    """Post-catalogue singleton screen (manuscript Methods 3.1): every singleton per-species
    cluster is scored against the family references; hits above the floor that lie outside the
    catalogue become inspection candidates C1.. (scanned, flanked, optionally sub-subfamily-
    assigned) for the viewer. Skipped when discovery did not run here (seeds given)."""
    cb = os.path.join(ctx.run, 'motifs', 'clean', 'clusters_by_file')
    out = ctx.d('candidates')
    if not os.path.isdir(cb):
        ctx.say('  no per-species clusters in this run (seeds were given): no singleton screen')
        return []
    cat_bed = os.path.join(ctx.run, 'catalogue', 'MASTER.merged.bed')
    env = ctx.env(MIN_BITS=ctx.cfg.data['expansion']['min_bits'])
    ctx.sh(['bash', os.path.join(ctx.scripts, '07_post_catalogue/singleton_screen.sh'), cb,
            os.path.join(ctx.run, 'tables', 'refs.rot'), out, cat_bed], '20_candidates', env=env)
    # candidates.tsv: id species scaffold start1 end1, ordered by best bits (descending)
    best = {}
    for p in read_tsv(os.path.join(out, 'sing_best.tsv'), 3):
        best[p[0]] = float(p[2])
    absent = [l.strip() for l in open(os.path.join(out, 'candidates_absent.txt')) if l.strip()] \
        if os.path.isfile(os.path.join(out, 'candidates_absent.txt')) else []
    absent.sort(key=lambda h: (-best.get(h, 0.0), h))
    ctab = os.path.join(out, 'candidates.tsv')
    with open(ctab, 'w') as f:
        for i, h in enumerate(absent, 1):
            sp, reg = h.split('|')[0], h.split('|')[1]
            scf, coords = reg.rsplit(':', 1)
            s, e = coords.split('-')
            f.write('C%d\t%s\t%s\t%s\t%s\n' % (i, sp, scf, s, e))
    ctx.manifest['candidates'] = len(absent)
    ctx.say('  singleton loci above the floor and outside the catalogue: %d' % len(absent))
    if not absent:
        return [ctab]
    tables = os.path.join(ctx.run, 'tables')
    ctx.sh(['python3', os.path.join(ctx.scripts, '07_post_catalogue/candidate_scan.py'), ctab,
            os.path.join(tables, 'genomes.tsv'), os.path.join(tables, 'families.tsv'), out,
            '--flank', ctx.cfg.data['orthology']['flank_bp']], '20_candidates')
    ctx.sh(['python3', os.path.join(ctx.scripts, '07_post_catalogue/candidate_flank5k.py'), out,
            os.path.join(tables, 'genomes.tsv'), os.path.join(tables, 'families.tsv'),
            '--sat-bits', ctx.cfg.data['expansion']['min_bits']], '20_candidates')
    sub = ctx.cfg.data.get('subfamily') or {}
    if sub.get('consensuses') and os.path.isdir(sub['consensuses']) and sub.get('members'):
        ctx.sh(['python3', os.path.join(ctx.scripts, '07_post_catalogue/candidate_subsubfamily.py'), out,
                os.path.join(tables, 'genomes.tsv'), os.path.join(tables, 'families.tsv'), sub['consensuses'],
                os.path.join(tables, 'members.tsv'), '--already-doubled',
                '--confident-delta', sub.get('confident_delta', 10)], '20_candidates')
    return [ctab, os.path.join(out, 'candidates_summary.json')]


st_candidates.plan = lambda ctx: ([os.path.join(ctx.run, 'catalogue', 'MASTER.merged.bed'), os.path.join(ctx.run, 'tables', 'refs.rot')],
                                  {'min_bits': ctx.cfg.data['expansion']['min_bits'], 'flank': ctx.cfg.data['orthology']['flank_bp']})


# ----------------------------------------------------------------------------- driver

def stage_ids() -> list:
    return [fn.sid for fn in STAGES]


def run(cfg: Config, from_stage: str | None = None, to_stage: str | None = None, force: bool = False, only: list | None = None) -> dict:
    run_dir = cfg.data['out_dir']
    os.makedirs(run_dir, exist_ok=True)
    ctx = Ctx(cfg=cfg, run=run_dir, force=force)
    ctx.log = open(os.path.join(ctx.d('logs'), 'satpipe.log'), 'a')
    ctx.say('=== satpipe %s run %s  %s' % (__version__, cfg.data['run_name'], dt.datetime.now().isoformat(timespec='seconds')))
    # manifest: config copy, tool versions, code version
    mpath = os.path.join(run_dir, 'manifest.json')
    if os.path.isfile(mpath):
        try:
            ctx.manifest = json.load(open(mpath))
        except Exception:
            ctx.manifest = {}
    ctx.manifest['satpipe_version'] = __version__
    ctx.manifest['config_file'] = cfg.path
    ctx.manifest['config'] = cfg.data
    ctx.manifest['tools'] = {p.tool.name: {'path': p.path, 'version': p.version, 'error': p.error} for p in probe_all(cfg.data.get('tools') or {})}
    try:
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        ctx.manifest['code_commit'] = subprocess.run(['git', '-C', here, 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip() or None
    except Exception:
        ctx.manifest['code_commit'] = None
    ctx.manifest.setdefault('runs', []).append(dt.datetime.now().isoformat(timespec='seconds'))
    ids = stage_ids()
    sel = ids
    if from_stage:
        if from_stage not in ids:
            raise StageError('unknown stage %r; stages: %s' % (from_stage, ' '.join(ids)))
        sel = sel[sel.index(from_stage):]
    if to_stage:
        if to_stage not in ids:
            raise StageError('unknown stage %r; stages: %s' % (to_stage, ' '.join(ids)))
        sel = sel[:sel.index(to_stage) + 1]
    if only:
        sel = [s for s in sel if s in only]
    # stage contracts as files: a run that starts from given seeds skips discovery (03-09)
    given_seeds = (cfg.data.get('inputs') or {}).get('seeds') or {}
    if given_seeds and all(f in given_seeds for f in cfg.families):
        skip = {'03_motifs', '04_filter', '05_rotate', '06_cluster', '07_mafft', '08_cross', '09_assign'}
        ctx.say('seeds given for every family: skipping ' + ' '.join(sorted(skip)))
        sel = [s for s in sel if s not in skip]
    try:
        for fn in STAGES:
            if fn.sid in sel:
                run_stage(fn, ctx)
    finally:
        json.dump(ctx.manifest, open(mpath, 'w'), indent=1, default=str)
        ctx.log.close()
    return ctx.manifest
