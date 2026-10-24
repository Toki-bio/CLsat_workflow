import textwrap

import pytest


def write(p, s):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(textwrap.dedent(s), encoding="utf-8")
    return p


@pytest.fixture(autouse=True)
def _long_tmp_paths_allowed(monkeypatch):
    """pytest's temp folders are longer than the out_dir limit that protects ssearch36 (satpipe.check.MAX_BASE_PATH);
    the fixtures never run ssearch36, so the limit is lifted here and tested on its own in test_check.py."""
    import satpipe.check
    monkeypatch.setattr(satpipe.check, "MAX_BASE_PATH", 10000)


@pytest.fixture
def project(tmp_path):
    """A tiny synthetic project: two genomes, their TRF .dat, one family, curated table."""
    d = tmp_path
    for code in ("spa", "spb"):
        write(d / "genomes" / f"{code}.fna", f">{code}_scf1\n" + "ACGT" * 50 + "\n")
        write(d / "trf" / f"{code}.dat", "Sequence: scf1\n\n1 500 146 3.4 146 95 2 300 ACGT ACGTACGT\n")
    write(d / "ref" / "SatX.cons.fa", ">SatX\n" + "ACGTTGCA" * 18 + "ACGT\n")   # 148 bp
    write(d / "ref" / "cluster_family.tsv", "# cluster\tfamily\n0\tSatX\n5\tSatX\tnote\n")
    write(d / "ref" / "manual.bed", "spa_scf1\t0\t1000\tspa|spa_scf1:1-1000|MANUAL\t0\t.\n")
    cfg = write(d / "configs" / "run.yaml", """
        run_name: toy
        out_dir: ../runs/toy
        genomes:
          spa: ../genomes/spa.fna
          spb: ../genomes/spb.fna
        inputs:
          trf_dat:
            spa: ../trf/spa.dat
            spb: ../trf/spb.dat
        families:
          SatX: {consensus: ../ref/SatX.cons.fa, mono: 148}
        family_assignment: ../ref/cluster_family.tsv
        manual_loci: ../ref/manual.bed
        """)
    return d, cfg
