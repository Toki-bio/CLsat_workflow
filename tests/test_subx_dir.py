import types

from satpipe.run import subx_dir


def ctx_for(tmp_path, cons):
    return types.SimpleNamespace(run=str(tmp_path / 'run'), cfg=types.SimpleNamespace(data={'subfamily': {'consensuses': cons}}))


def test_fasta_becomes_doubled_files(tmp_path):
    # the final-panel config gives the 9 sub-subfamily consensuses as one FASTA (fig1b_rotated.fasta);
    # stage 18 once skipped it silently and wrote an empty table
    fa = tmp_path / 'subs.fasta'
    fa.write_text('>CLsat1.1 x\nAAAG\nCT\n>CLsat3.2\nGGTTA\n')
    d = subx_dir(ctx_for(tmp_path, str(fa)))
    assert d == str(tmp_path / 'run' / 'tables' / 'subx')
    assert (tmp_path / 'run' / 'tables' / 'subx' / 'CLsat1.1.fa').read_text() == '>CLsat1.1\nAAAGCTAAAGCT\n'
    assert (tmp_path / 'run' / 'tables' / 'subx' / 'CLsat3.2.fa').read_text() == '>CLsat3.2\nGGTTAGGTTA\n'


def test_directory_used_as_given(tmp_path):
    d = tmp_path / 'subx'
    d.mkdir()
    assert subx_dir(ctx_for(tmp_path, str(d))) == str(d)


def test_missing_or_absent(tmp_path):
    assert subx_dir(ctx_for(tmp_path, str(tmp_path / 'nope.fa'))) is None
    assert subx_dir(ctx_for(tmp_path, None)) is None
