"""satpipe.locus_fasta: each locus +/- flank cut from its own assembly, clipped at contig ends, byte-identical reruns."""
import gzip
import os
import random
import shutil
import subprocess

import pytest

from satpipe.locus_fasta import write_locus_fasta, read_viewer_bed

SAMTOOLS = shutil.which('samtools')


def _genome(path, contigs):
    with open(path, 'w') as f:
        for name, seq in contigs.items():
            f.write('>%s\n' % name + ''.join(seq[i:i + 70] + '\n' for i in range(0, len(seq), 70)))
    subprocess.run([SAMTOOLS, 'faidx', path], check=True)


@pytest.mark.skipif(SAMTOOLS is None, reason='samtools not on PATH')
def test_cut_clip_and_same_contig_names(tmp_path):
    rnd = random.Random(1)
    a = {'ctg1': ''.join(rnd.choice('ACGT') for _ in range(3000))}
    b = {'ctg1': ''.join(rnd.choice('ACGT') for _ in range(2000))}   # same contig name, other assembly
    ga, gb = str(tmp_path / 'a.fa'), str(tmp_path / 'b.fa')
    _genome(ga, a)
    _genome(gb, b)
    bed = tmp_path / 'MASTER.merged.bed'
    bed.write_text('ctg1\t1000\t1200\tA|ctg1:1001-1200|L1\n'      # flank fits on both sides
                   'ctg1\t100\t1900\tB|ctg1:101-1900|L2\n')        # flank clipped at both contig ends
    rows = read_viewer_bed(str(bed))
    out = tmp_path / 'v'
    res = write_locus_fasta(rows, {'A': ga, 'B': gb}, str(out), flank=500, samtools=SAMTOOLS)
    assert res['n'] == 2

    def read(lid):
        t = gzip.open(str(out / 'seq' / (lid + '.fa.gz')), 'rt').read()
        head, _, body = t.partition('\n')
        return head, body.replace('\n', '')

    h, s = read('L1')
    assert h == '>L1 A ctg1:501-1700' and s == a['ctg1'][500:1700]
    h, s = read('L2')
    assert h == '>L2 B ctg1:1-2000' and s == b['ctg1']
    first = (out / 'seq' / 'L1.fa.gz').read_bytes()
    write_locus_fasta(rows, {'A': ga, 'B': gb}, str(out), flank=500, samtools=SAMTOOLS)
    assert (out / 'seq' / 'L1.fa.gz').read_bytes() == first
