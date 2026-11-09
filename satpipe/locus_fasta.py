"""Per-locus FASTA files for the viewer (viewer.locus_fasta; owner 2026-11-09).

Each catalogue locus (and stage-20 candidate) is cut from its own assembly with `flank` letters on each side,
clipped to the contig, and written as <out>/seq/<locus>.fa.gz with the header

    >L297 Fmex h2tg000025l:2525015-2829895

(1-based inclusive genomic range). The viewer's "Get FASTA" and "Analyze Monomers" slice the view from this file, so
loci on assemblies that are not in NCBI get their letters too. gzip is written with mtime 0: two runs are
byte-identical.

Standalone use on an existing viewer folder (adds seq/ and the `locus_fasta` key in viewer_config.json):

    python -m satpipe.locus_fasta --run-dir RUN --viewer-dir RUN/viewer [--flank 5000] [--samtools samtools]
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import tempfile

SEQ_DIR = 'seq'


def read_viewer_bed(path: str) -> list:
    """MASTER.merged.bed of a viewer folder -> [(scf, start0, end, sp, lid)]."""
    rows = []
    for line in open(path):
        p = line.rstrip('\n').split('\t')
        if len(p) < 4:
            continue
        sp, _, lid = p[3].split('|')
        rows.append((p[0], int(p[1]), int(p[2]), sp, lid))
    return rows


def write_locus_fasta(rows: list, genomes: dict, out_dir: str, flank: int = 5000, samtools: str = 'samtools',
                      wrap: int = 60) -> dict:
    """Cut every row (scf, start0, end, sp, lid) +/- flank from genomes[sp]; returns {'n', 'bytes', 'flank'}."""
    sd = os.path.join(out_dir, SEQ_DIR)
    os.makedirs(sd, exist_ok=True)
    lengths = {}
    for sp in {r[3] for r in rows}:
        for line in open(genomes[sp] + '.fai'):
            p = line.split('\t')
            lengths[(sp, p[0])] = int(p[1])
    by_sp = {}
    for scf, s, e, sp, lid in rows:
        L = lengths[(sp, scf)]
        a, b = max(0, s - flank), min(L, e + flank)
        by_sp.setdefault(sp, []).append((lid, scf, a + 1, b))
    n = nbytes = 0
    for sp, items in by_sp.items():
        with tempfile.TemporaryDirectory() as tmp:
            reg = os.path.join(tmp, 'regions.txt')
            with open(reg, 'w') as f:
                for lid, scf, a1, b in items:
                    f.write('%s:%d-%d\n' % (scf, a1, b))
            res = subprocess.run([samtools, 'faidx', '-r', reg, genomes[sp]],
                                 capture_output=True, text=True)
            if res.returncode != 0:
                raise SystemExit('locus_fasta: samtools faidx failed for %s: %s' % (sp, res.stderr.strip()))
        recs = res.stdout.split('>')[1:]
        if len(recs) != len(items):
            raise SystemExit('locus_fasta: %s: asked %d regions, samtools gave %d' % (sp, len(items), len(recs)))
        for (lid, scf, a1, b), rec in zip(items, recs):
            name, _, seq = rec.partition('\n')
            seq = seq.replace('\n', '')
            if name.strip() != '%s:%d-%d' % (scf, a1, b) or len(seq) != b - a1 + 1:
                raise SystemExit('locus_fasta: %s: region %s:%d-%d came back as %s (%d letters)'
                                 % (lid, scf, a1, b, name.strip(), len(seq)))
            txt = '>%s %s %s:%d-%d\n' % (lid, sp, scf, a1, b) + \
                  ''.join(seq[i:i + wrap] + '\n' for i in range(0, len(seq), wrap))
            path = os.path.join(sd, lid + '.fa.gz')
            with open(path, 'wb') as fo, gzip.GzipFile(filename='', mode='wb', fileobj=fo, mtime=0, compresslevel=9) as gz:
                gz.write(txt.encode())
            n += 1
            nbytes += os.path.getsize(path)
    return {'n': n, 'bytes': nbytes, 'flank': flank}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--run-dir', required=True)
    ap.add_argument('--viewer-dir', required=True)
    ap.add_argument('--flank', type=int, default=5000)
    ap.add_argument('--samtools', default='samtools')
    a = ap.parse_args(argv)
    from .viewer import load_run_config
    cfg = load_run_config(os.path.abspath(a.run_dir))
    rows = read_viewer_bed(os.path.join(a.viewer_dir, 'MASTER.merged.bed'))
    res = write_locus_fasta(rows, dict(cfg.genomes), a.viewer_dir, a.flank, a.samtools)
    cp = os.path.join(a.viewer_dir, 'viewer_config.json')
    if os.path.isfile(cp):
        c = json.load(open(cp))
        c['locus_fasta'] = {'dir': SEQ_DIR, 'flank': res['flank'], 'n': res['n']}
        json.dump(c, open(cp, 'w'), separators=(',', ':'))
    print('locus_fasta: %d loci, %.1f MB -> %s' % (res['n'], res['bytes'] / 1e6, os.path.join(a.viewer_dir, SEQ_DIR)))


if __name__ == '__main__':
    main()
