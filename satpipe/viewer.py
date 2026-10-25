"""`satpipe viewer`: assemble the viewer directory for a run.

Writes <run>/viewer/ with the data files the CLsat viewer reads, in the formats it expects
(see viewer/README.md in CLsat_workflow and the viewer repository), plus `viewer_config.json`,
which carries everything that used to be pasted into index.html: family names, monomer lengths,
colours and consensus sequences, species order and display names, optional sub-subfamily
consensuses, and the three inspection lists (ORTHO_LOCI, FLAGGED_LOCI, FULL_TABLE) generated
from the JSON products by satpipe.viewer_lists.
"""
from __future__ import annotations

import gzip
import json
import os
import shutil
import subprocess
import sys

from . import __version__
from .config import Config
from .run import Ctx, _expand_cmd, read_tsv
from .viewer_lists import build_lists, read_assign, read_master

VENDORED_HTML = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'viewer_html', 'index.html')
PALETTE = ['#60a5fa', '#34d399', '#f472b6', '#fbbf24', '#a78bfa', '#f87171', '#2dd4bf', '#fb923c', '#c084fc', '#4ade80']


def _gzip(src: str, dst: str):
    with open(src, 'rb') as fi, gzip.GzipFile(dst, 'wb', compresslevel=9, mtime=0) as fo:
        shutil.copyfileobj(fi, fo)


def _read_fasta(path: str) -> dict:
    out, name = {}, None
    for line in open(path):
        if line.startswith('>'):
            name = line[1:].split()[0]; out[name] = ''
        elif name:
            out[name] += line.strip().upper()
    return out


def load_run_config(run_dir: str) -> Config:
    m = json.load(open(os.path.join(run_dir, 'manifest.json')))
    return Config(path=m.get('config_file', ''), base_dir=os.path.dirname(m.get('config_file', '')) or run_dir, data=m['config'])


def build(run_dir: str, out_dir: str | None = None, viewer_html: str | None = None, flank_scan: bool = True) -> dict:
    run_dir = os.path.abspath(run_dir)
    cfg = load_run_config(run_dir)
    out = out_dir or os.path.join(run_dir, 'viewer')
    os.makedirs(out, exist_ok=True)
    ctx = Ctx(cfg=cfg, run=run_dir)
    cat = os.path.join(run_dir, 'catalogue')
    reg = os.path.join(run_dir, 'regions')
    ortho = os.path.join(run_dir, 'ortho')
    fams = list(cfg.families)
    codes = list(cfg.genomes)
    say = lambda s: print(s, flush=True)

    # ---- 1. catalogue BED in viewer format: scf start end sp|scf:start1-end|lid --------------
    # the genome of each locus is col 5 of the segmented catalogue (never looked up by contig name:
    # genomes may share contig names)
    rows = []
    seg = read_tsv(os.path.join(cat, 'MASTER.loci.segmented.bed'))
    if any(len(p) < 5 for p in seg):
        raise SystemExit('viewer: MASTER.loci.segmented.bed has no genome column (col 5); it was written by an older '
                         'satpipe. Rerun from 14_arbitrate.')
    for p in seg:
        rows.append((p[0], int(p[1]), int(p[2]), p[4], p[3]))
    # candidates from stage 20 (optional)
    cand = None
    cpath = os.path.join(run_dir, 'candidates', 'candidates_summary.json')
    if os.path.isfile(cpath):
        cand = json.load(open(cpath))
        for cid, s in sorted(cand['summary'].items(), key=lambda kv: int(kv[0][1:]) if kv[0][1:].isdigit() else 0):
            rows.append((s['scf'], int(s['start']), int(s['end']), s['sp'], cid))
    with open(os.path.join(out, 'MASTER.merged.bed'), 'w') as f:
        for scf, s, e, sp, lid in rows:
            f.write('%s\t%d\t%d\t%s|%s:%d-%d|%s\n' % (scf, s, e, sp, scf, s + 1, e, lid))

    # ---- 2. assignment (no header) + 3. bits matrix ---------------------------------------
    shutil.copyfile(os.path.join(cat, 'MASTER.final_assignment.tsv'), os.path.join(out, 'MASTER.final_assignment.clean.tsv'))
    shutil.copyfile(os.path.join(cat, 'MASTER.locus_bits.tsv'), os.path.join(out, 'MASTER.locus_bits.tsv'))
    if cand:
        with open(os.path.join(out, 'MASTER.final_assignment.clean.tsv'), 'a') as fa, open(os.path.join(out, 'MASTER.locus_bits.tsv'), 'a') as fb:
            for cid, s in sorted(cand['summary'].items(), key=lambda kv: int(kv[0][1:]) if kv[0][1:].isdigit() else 0):
                fa.write('%s\t%s\t%s\t%s\t%s\t%s\tCANDIDATE_UNCATALOGUED\t%s\n' % (cid, s['call'], s['call_bits'], s['comp'] or '', s['comp_bits'], s['delta'], s['call']))
                fb.write(cid + ''.join('\t%s' % s['bits'].get(fam, 0) for fam in fams) + '\n')

    # ---- 4. per-family bits tracks (gzip) ---------------------------------------------------
    for fam in fams:
        src = os.path.join(cat, 'MASTER_%s.bits.tsv' % fam)
        tmp = os.path.join(out, 'MASTER_%s.bits.tsv' % fam)
        shutil.copyfile(src, tmp)
        cb = os.path.join(run_dir, 'candidates', 'CAND_%s.bits.tsv' % fam)
        if cand and os.path.isfile(cb):
            with open(tmp, 'a') as fo, open(cb) as fi:
                fo.write(fi.read())
        _gzip(tmp, tmp + '.gz')
        os.remove(tmp)

    # ---- 5. flank windows (background track) ----------------------------------------------
    flank_view = 5000
    if flank_scan:
        fdir = os.path.join(run_dir, 'viewer_flankscan')
        os.makedirs(fdir, exist_ok=True)
        bed = os.path.join(fdir, 'flanks.bed')
        lengths = {}
        for code, g in cfg.genomes.items():
            for line in open(g + '.fai'):
                p = line.split('\t'); lengths[(code, p[0])] = int(p[1])
        with open(bed, 'w') as f:
            for scf, s, e, sp, lid in rows:
                L = lengths.get((sp, scf), e)
                a, b = max(0, s - flank_view), s
                if b > a:
                    f.write('%s\t%d\t%d\t%s|%s:%d-%d|FLANK\n' % (scf, a, b, sp, scf, a + 1, b))
                a, b = e, min(L, e + flank_view)
                if b > a:
                    f.write('%s\t%d\t%d\t%s|%s:%d-%d|FLANK\n' % (scf, a, b, sp, scf, a + 1, b))
        # earlier flank scans are reused only for the same flank windows: a rebuilt catalogue changes flanks.bed, and
        # stale scans would then lose the new loci's flanks (GLM audit #56)
        import hashlib
        sig = hashlib.sha256(open(bed, 'rb').read() + repr(sorted(cfg.families)).encode()).hexdigest()
        sigf = os.path.join(fdir, 'flanks.bed.sha256')
        if not os.path.isfile(sigf) or open(sigf).read().strip() != sig:
            for fam in fams:
                for f in os.listdir(fdir):
                    if f.startswith('flank_' + fam + '.'):
                        os.remove(os.path.join(fdir, f))
            open(sigf, 'w').write(sig + '\n')
        tmp = os.path.join(out, 'MASTER_flank.bits.tsv')
        with open(tmp, 'w') as fo:
            fo.write('scaffold\tpos\tbits\tclsat\tspecies\n')   # species: the genome (prefix of the scan name)
            for fam in fams:
                pre = os.path.join(fdir, 'flank_' + fam)
                if not os.path.isfile(pre + '.bits.tsv'):
                    say('  flank scan %s' % fam)
                    ctx.sh(_expand_cmd(ctx, bed, fam, pre, True), 'viewer_flankscan', log_name='viewer_flankscan_' + fam)
                with open(pre + '.bits.tsv') as fi:
                    fi.readline()
                    for line in fi:
                        p = line.rstrip('\n').split('\t')
                        if len(p) >= 4:
                            fo.write('%s\t%s\t%s\t%s\t%s\n' % (p[0], p[1], p[2], fam, p[3].split('|', 1)[0]))
            cf = os.path.join(run_dir, 'candidates', 'CAND_FLANK.bits.tsv')
            if cand and os.path.isfile(cf):
                # candidate_flank5k.py scans each candidate's flanks with its called family;
                # rows: scaffold pos0 bits "sp|scf:start1-end|CANDFLANK" FLANK sp
                call_of = {'%s|%s:%d-%d' % (s['sp'], s['scf'], s['start'] + 1, s['end']): s['call'] for s in cand['summary'].values()}
                for p in read_tsv(cf, 4):
                    fam0 = call_of.get(p[3].rsplit('|', 1)[0], fams[0])
                    fo.write('%s\t%d\t%s\t%s\t%s\n' % (p[0], int(p[1]) + 1, p[2], fam0, p[3].split('|', 1)[0]))
        _gzip(tmp, tmp + '.gz')
        os.remove(tmp)

    # ---- 6. JSON products -------------------------------------------------------------------
    flanks = json.load(open(os.path.join(ortho, 'flanks.json')))
    if cand:
        for cid, fl in cand['flanks'].items():
            s = cand['summary'][cid]
            fam0 = s['call']
            mono = cfg.families[fam0]['mono']
            flanks[cid] = {'sp': s['sp'], 'scf': s['scf'], 'start': s['start'], 'end': s['end'], 'subtype': fam0, 'mono': mono,
                           'length': s['end'] - s['start'], 'n_monomers': round((s['end'] - s['start']) / mono),
                           'LF': dict(fl['LF'], **{'class': 'candidate'}), 'RF': dict(fl['RF'], **{'class': 'candidate'})}
    json.dump(flanks, open(os.path.join(out, 'flanks.json'), 'w'), separators=(',', ':'))
    shutil.copyfile(os.path.join(ortho, 'flank_orthologs.json'), os.path.join(out, 'flank_orthologs.json'))
    subfam = {}
    for name in ('subfamily.json', 'terminal_monomers.json'):
        p = os.path.join(reg, name)
        if os.path.isfile(p):
            shutil.copyfile(p, os.path.join(out, name))
    if os.path.isfile(os.path.join(out, 'subfamily.json')):
        subfam = json.load(open(os.path.join(out, 'subfamily.json')))
        if cand:
            for cid, s in cand['summary'].items():
                if 'subfamily' in s:
                    subfam[cid] = {'major': s['call'], 'subfamily': s['subfamily'], 'pid': s.get('ss_pid'), 'bits': s.get('ss_bits'),
                                   'second': s.get('ss_second'), 'delta': s.get('ss_delta'), 'confident': s.get('ss_confident')}
            json.dump(subfam, open(os.path.join(out, 'subfamily.json'), 'w'), separators=(',', ':'))

    # ---- 7. scaffold lengths, consensuses ---------------------------------------------------
    with open(os.path.join(out, 'scaffold_lengths.tsv'), 'w') as f:
        f.write('scaffold\tlength\tspecies\n')
        want = {(sp, scf) for scf, s, e, sp, lid in rows}   # one row per genome + contig holding a locus
        for code, g in cfg.genomes.items():
            for line in open(g + '.fai'):
                p = line.split('\t')
                if (code, p[0]) in want:
                    f.write('%s\t%s\t%s\n' % (p[0], p[1], code))
    cons = {}
    for fam, fd in cfg.families.items():
        shutil.copyfile(fd['consensus'], os.path.join(out, fam + '.cons.fa'))
        cons[fam] = ''.join(_read_fasta(fd['consensus']).values())

    # ---- 8. viewer_config.json --------------------------------------------------------------
    master = read_master(os.path.join(out, 'MASTER.merged.bed'))
    assign = read_assign(os.path.join(out, 'MASTER.final_assignment.clean.tsv'))
    ortho_json = json.load(open(os.path.join(out, 'flank_orthologs.json')))
    ol, fl_list, full = build_lists(master, assign, flanks, ortho_json, subfam)
    sub = cfg.data.get('subfamily') or {}
    subfamilies = {}
    cdir = sub.get('consensuses')
    if cdir and os.path.isdir(cdir):
        for fam, members in (sub.get('members') or {}).items():
            for m in members:
                p = os.path.join(cdir, m + '.fa')
                if os.path.isfile(p):
                    seq = ''.join(_read_fasta(p).values())
                    # the Darevskia subx/ files are DOUBLED consensuses
                    half = seq[:len(seq) // 2]
                    if len(seq) % 2 == 0 and half == seq[len(seq) // 2:]:
                        seq = half
                    subfamilies[m] = {'family': fam, 'seq': seq}
    view = cfg.data.get('viewer') or {}
    labels = view.get('species_names') or {}
    n_cat = sum(1 for r in master if not r[0].startswith('C'))
    n_cand = sum(1 for r in master if r[0].startswith('C'))
    n_flag = sum(1 for r in fl_list if not r['id'].startswith('C'))
    n_contra = sum(1 for r in fl_list if r['cat'] == 'UNIQUE_SEED_CONTRADICTION')
    o = cfg.data['orthology']
    config = {
        'generator': 'satpipe %s' % __version__,
        'title': view.get('title') or '%s locus viewer' % cfg.data['run_name'],
        'run_name': cfg.data['run_name'],
        'families': [{'name': fam, 'mono': cfg.families[fam]['mono'],
                      'colour': cfg.families[fam].get('colour') or PALETTE[i % len(PALETTE)],
                      'consensus': cons[fam], 'bits_file': 'MASTER_%s.bits.tsv' % fam} for i, fam in enumerate(fams)],
        'species': [{'code': c, 'label': labels.get(c, c)} for c in codes],
        'subfamilies': subfamilies,
        'window_step': max(cfg.families[f]['mono'] for f in fams),
        'sat_bits': cfg.data['expansion']['min_bits'],
        'edge_tol': cfg.data['finalisation']['edge_tol'],
        'flank_view': flank_view,
        'orthology': {'min_pid': o['min_pid'], 'min_aln': o['min_aln'], 'flank_bp': o['flank_bp']},
        'counts': {'catalogue_loci': n_cat, 'candidates': n_cand, 'flagged': n_flag, 'contradictions': n_contra,
                   'anchorable': len(ol), 'with_ortholog': sum(1 for r in ol if r['n'] > 0)},
        'lists': {'ortho': ol, 'flagged': fl_list, 'full': full},
        'notes': {
            'ortho': 'Positional orthology from unique flanking sequence (>=%g%% identity, >=%d bp), single best hit per target '
                     'assembly. "Orthologous to" lists distinct partner assemblies then the total partner-locus count.'
                     % (o['min_pid'], o['min_aln']),
            'flag': '%d catalogue loci annotated by arbitration%s. Only %d %s a seed/rescan contradiction (the rescan scored another '
                    'family >=%g bits above the seed; the seed family is retained). The rest are close-margin calls (segment or '
                    'seed-conflict resolved with a margin below %g bits) or loci without a seed.'
                    % (n_flag, (', plus the %d uncatalogued candidates (not catalogue members)' % n_cand) if n_cand else '',
                       n_contra, 'is' if n_contra == 1 else 'are', cfg.data['arbitration']['delta'], cfg.data['arbitration']['delta_ambig']),
            'full': 'Every locus (%d, all %d assemblies) with its sub-subfamily call and both flanks\' cross-assembly orthology '
                    'status (ssearch36, >=%g%% identity, >=%d bp, single best hit per target assembly). L/R badges show each '
                    'flank\'s class.' % (n_cat, len(codes), o['min_pid'], o['min_aln']),
        },
    }
    json.dump(config, open(os.path.join(out, 'viewer_config.json'), 'w'), separators=(',', ':'))

    # ---- 9. the page --------------------------------------------------------------------------
    html = viewer_html or VENDORED_HTML
    if os.path.isfile(html):
        shutil.copyfile(html, os.path.join(out, 'index.html'))
    else:
        say('  WARNING: viewer page not found (%s); data files and viewer_config.json written without index.html' % html)
    open(os.path.join(out, '.nojekyll'), 'w').close()
    say('viewer: %d catalogue loci, %d candidates, %d anchorable, %d flagged -> %s' % (n_cat, n_cand, len(ol), n_flag, out))
    return config
