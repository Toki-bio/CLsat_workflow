import pytest
import os

from satpipe import cli
from satpipe.check import FAIL, OK, WARN, run_check
from satpipe.tools import TOOLS, probe_all


def levels(rep, section=None):
    return [(i.level, i.message) for i in rep.items if section is None or i.section == section]


def test_fixture_passes_without_tools(project):
    d, cfg = project
    rep = run_check(str(cfg), do_tools=False)
    assert not rep.failed, rep.to_text()
    msgs = " ".join(m for _, m in levels(rep))
    assert "148 bp = mono" in msgs
    assert "2 cluster->family rows" in msgs
    assert "manual_loci: 1 intervals" in msgs
    assert "stage script references present" in msgs


def test_mono_mismatch_is_warning(project):
    d, cfg = project
    cfg.write_text(cfg.read_text().replace("mono: 148", "mono: 146"))
    rep = run_check(str(cfg), do_tools=False)
    assert not rep.failed
    assert any(lv == WARN and "consensus length 148 != mono 146" in m for lv, m in levels(rep, "families"))


def test_unknown_family_in_table_fails(project):
    d, cfg = project
    (d / "ref" / "cluster_family.tsv").write_text("0\tSatX\n7\tSatQ\n")
    rep = run_check(str(cfg), do_tools=False)
    assert rep.failed
    assert any("family 'SatQ' is not in" in m for _, m in levels(rep, "family_assignment"))


def test_family_without_cluster_fails(project):
    d, cfg = project
    cfg.write_text(cfg.read_text().replace(
        "families:\n", "families:\n  SatY: {consensus: ../ref/SatX.cons.fa, mono: 148}\n"))
    rep = run_check(str(cfg), do_tools=False)
    assert any(lv == FAIL and "family SatY has no cluster" in m for lv, m in levels(rep, "family_assignment"))


def test_missing_genome_fails(project):
    d, cfg = project
    os.remove(d / "genomes" / "spb.fna")
    rep = run_check(str(cfg), do_tools=False)
    assert rep.failed
    assert any(lv == FAIL and "spb" in m and "not found" in m for lv, m in levels(rep, "genomes"))


def test_bad_bed_is_reported(project):
    d, cfg = project
    (d / "ref" / "manual.bed").write_text("scf\t10\t5\tspa|x\n")
    rep = run_check(str(cfg), do_tools=False)
    assert any(lv == FAIL and "0 <= start < end" in m for lv, m in levels(rep, "manual_loci"))


def test_tools_missing_are_reported(project, monkeypatch, tmp_path):
    (tmp_path / "empty").mkdir()
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    rep = run_check(str(project[1]), do_files=False, do_tools=True)
    assert rep.failed
    names = {m.split(":")[0] for lv, m in levels(rep, "tools")}
    assert {"trf", "vsearch", "ssearch36", "gawk", "samtools", "bedtools", "parallel", "mafft"} <= names
    required = {t.name for t in TOOLS if t.required}
    failed = {m.split(":")[0] for lv, m in levels(rep, "tools") if lv == FAIL}
    assert required <= failed


def test_tool_override_and_bin_dir(project, tmp_path):
    d, cfg = project
    fake = tmp_path / "bin"
    fake.mkdir()
    exe = fake / "vsearch"
    exe.write_text("#!/bin/sh\necho 'vsearch v2.30.0_linux_x86_64' >&2\n")
    exe.chmod(0o755)
    cfg.write_text(cfg.read_text() + f"\ntools:\n  bin_dir: {fake.as_posix()}\n")
    rep = run_check(str(cfg), do_files=False, do_tools=True)
    msgs = dict((m.split(":")[0].split(" ")[0], (lv, m)) for lv, m in levels(rep, "tools"))
    assert msgs["bin_dir"][0] == OK
    # the fake is located through bin_dir on every platform (it may not be runnable on Windows)
    assert "vsearch" in msgs and "not found" not in msgs["vsearch"][1]
    assert msgs["vsearch"][1].split(": ")[1].startswith(fake.as_posix()) or "v2.30.0" in msgs["vsearch"][1]


def test_cli_exit_codes(project, capsys):
    d, cfg = project
    assert cli.main(["check", str(cfg), "--no-tools"]) == 0
    assert cli.main(["check", str(cfg), "--no-tools", "--json"]) == 0
    out = capsys.readouterr().out
    assert '"failed": false' in out
    (d / "ref" / "cluster_family.tsv").write_text("x\tSatX\n")
    assert cli.main(["check", str(cfg), "--no-tools"]) == 1
    assert cli.main(["run", "--list-stages", str(cfg)]) == 0
    assert "02_trf" in capsys.readouterr().out
    assert cli.main(["stages"]) == 0
    assert cli.main(["config", "--schema-doc"]) == 0
    assert cli.main(["config", str(cfg)]) == 0
    assert cli.main(["config"]) == 2


def test_probe_all_runs_without_crashing():
    probes = probe_all({})
    assert len(probes) == len(TOOLS)


def test_long_out_dir_fails(project, monkeypatch):
    # ssearch36 cannot open a file whose path is 120+ characters; a run folder that leaves no room fails the check
    import satpipe.check
    d, cfg = project
    monkeypatch.setattr(satpipe.check, "MAX_BASE_PATH", 20)
    assert cli.main(["check", str(cfg), "--no-tools"]) == 1


def test_unknown_anchor_family_and_tool_key_are_config_errors(project):
    # both passed the check and failed only later (anchor) or were silently ignored (tools) (GLM audit #125, #86)
    from satpipe.config import ConfigError, load_config
    d, cfg = project
    base = cfg.read_text()
    cfg.write_text(base + "\nanchor: {family: NoSuchFamily}\n")
    with pytest.raises(ConfigError) as e:
        load_config(str(cfg))
    assert "anchor.family" in str(e.value)
    cfg.write_text(base + "\ntools: {ssearch_36: /opt/x/ssearch36}\n")
    with pytest.raises(ConfigError) as e:
        load_config(str(cfg))
    assert "tools.ssearch_36" in str(e.value)
