import os

import pytest

from satpipe.config import ConfigError, load_config
from satpipe.schema import SCHEMA, to_markdown, validate


def test_defaults_are_published_values(project):
    d, cfg = project
    x = load_config(str(cfg)).data
    assert x["discovery"]["trf"] == "2 5 7 80 10 50 2000"
    assert x["discovery"]["period"] == [80, 300]
    assert x["discovery"]["strand"] == "plus"
    assert x["discovery"]["cluster_id"] == 0.90 and x["discovery"]["cross_cluster_id"] == 0.90
    assert x["expansion"]["min_bits"] == 50 and x["expansion"]["min_aln_frac"] == 0.67
    assert x["expansion"]["max_steps"] == 500 and x["expansion"]["probe_mono"] == 100
    assert x["arbitration"] == {"delta": 20.0, "delta_ambig": 15.0, "seg_bin": 1000, "seg_min_run": 10,
                                "seg_min_locus": 2000, "short_segments": "drop",
                                "bin_min_bits": None, "max_gap_monomers": 100, "bin_share": 0.75}   # bin-majority family: not in the published run (2026-10-25)
    assert x["orthology"] == {"flank_bp": 600, "guard_bp": 0, "min_aln": 200, "min_pid": 80.0}
    assert x["subfamily"]["confident_delta"] == 10.0
    assert x["scan"] == {"win_mono": 2, "step_mono": 1, "min_bits": 0.0}


def test_paths_resolve_relative_to_config(project):
    d, cfg = project
    c = load_config(str(cfg))
    assert c.data["genomes"]["spa"] == (d / "genomes" / "spa.fna").as_posix()
    assert c.data["out_dir"] == (d / "runs" / "toy").as_posix()
    assert c.data["families"]["SatX"]["consensus"].endswith("ref/SatX.cons.fa")


def test_unknown_key_is_error(project):
    d, cfg = project
    cfg.write_text(cfg.read_text() + "\nexpansion:\n  min_bitz: 3\n")
    with pytest.raises(ConfigError) as e:
        load_config(str(cfg))
    assert any("expansion.min_bitz: unknown key" in m for m in e.value.errors)


def test_missing_required_and_bad_types_all_reported():
    errors = []
    validate({"run_name": "x y", "genomes": {}, "families": {"F": {"mono": "146"}},
              "discovery": {"strand": "minus", "cluster_id": 0}}, SCHEMA, "", errors, None)
    joined = "\n".join(errors)
    assert "out_dir: required" in joined
    assert "run_name: 'x y' does not match" in joined
    assert "genomes: expected at least 1 entries" in joined
    assert "families.F.consensus: required" in joined
    assert "families.F.mono: expected an integer" in joined
    assert "discovery.strand: must be one of" in joined
    assert "discovery.cluster_id: must be > 0.0" in joined
    assert "family_assignment: required" in joined


def test_cross_checks(project):
    d, cfg = project
    cfg.write_text(cfg.read_text() + "\nsubfamily:\n  members:\n    SatY: [a]\n  consensuses: ../ref\n"
                   "viewer:\n  species_names:\n    zzz: Z\n")
    with pytest.raises(ConfigError) as e:
        load_config(str(cfg))
    j = "\n".join(e.value.errors)
    assert "subfamily.members.SatY: not a family" in j
    assert "viewer.species_names.zzz: not a code" in j


def test_trf_dat_code_must_be_a_genome(project):
    d, cfg = project
    cfg.write_text(cfg.read_text().replace("spb: ../trf/spb.dat", "spc: ../trf/spb.dat"))
    with pytest.raises(ConfigError) as e:
        load_config(str(cfg))
    assert any("inputs.trf_dat.spc: not a code" in m for m in e.value.errors)
    # partial coverage is a warning, not an error
    cfg.write_text(cfg.read_text().replace("    spc: ../trf/spb.dat\n", ""))
    c = load_config(str(cfg))
    assert any("covers only part" in w for w in c.warnings)


def test_schema_doc_lists_every_leaf():
    md = to_markdown()
    for key in ("`run_name`", "`genomes`", "`families.<key>.mono`", "`expansion.min_aln_frac`",
                "`arbitration.short_segments`", "`references.rotated_consensus`", "`tools`"):
        assert key in md, key


def test_darevskia_config_is_schema_valid():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    c = load_config(os.path.join(here, "configs", "darevskia_clsat.yaml"))
    assert set(c.genomes) == {"arm", "mix", "nai", "unm", "unp", "val"}
    assert [c.families[f]["mono"] for f in ("CLsat1", "CLsat2", "CLsat3", "CLsat4")] == [146, 145, 146, 116]
    assert c.data["arbitration"]["short_segments"] == "drop"
    assert c.data["tools"]["bin_dir"] == "/staging/conda/envs/bioinfo/bin"
    assert c.data["family_assignment"].endswith("/satpipe/stages/02_seed_expansion_scan/cluster_family.tsv")
    assert os.path.isfile(c.data["family_assignment"])
    assert all(os.path.isfile(f["consensus"]) for f in c.families.values())
    assert os.path.isfile(c.data["manual_loci"])
    assert not c.warnings


def test_example_minimal_is_schema_valid():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    c = load_config(os.path.join(here, "configs", "example_minimal.yaml"))
    assert c.data["expansion"]["jobs"] == 16
