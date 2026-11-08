#!/usr/bin/env python3
"""Write the config variants of the WP3 golden test (PORT_SPEC 5a) into tests/port/bad/ (one YAML file each).

    python3 tests/port/make_bad_configs.py tests/port/bad

Each variant is one edit of BASE. Most are invalid (schema errors of schema.py:252-335, cross-checks of
config.py:42-82); a few are valid on purpose (warning only, path forms, number forms, duplicate key, `yes`).
tests/port/golden_config.sh compares the bash loader with tools/py_oracle.py on every file."""
import os
import sys

BASE = """run_name: t1
out_dir: out
genomes:
  spa: g/a.fa
  spb: g/b.fa
families:
  SatX:
    consensus: ref/x.fa
    mono: 146
  SatY: {consensus: ref/y.fa, mono: 80}
family_assignment: ref/assign.tsv
"""


def drop(line):
    assert line in BASE, line
    return BASE.replace(line, "")


def sub(a, b):
    assert a in BASE, a
    return BASE.replace(a, b, 1)


def add(text):
    return BASE + text


V = {
    # unknown keys at every level
    "unknown_top": add("foo: 1\n"),
    "unknown_nested": add("expansion:\n  min_bitz: 3\n"),
    "unknown_family_key": sub("    mono: 146\n", "    mono: 146\n    colur: '#aabbcc'\n"),
    "unknown_inputs_key": add("inputs: {foo: 1}\n"),
    "unknown_several": add("zzz: 1\naaa: 2\nhapdup: {pairz: 1, edge_tol: -1}\n"),
    # missing required
    "missing_run_name": drop("run_name: t1\n"),
    "missing_out_dir": drop("out_dir: out\n"),
    "missing_genomes": drop("genomes:\n  spa: g/a.fa\n  spb: g/b.fa\n"),
    "missing_families": drop("families:\n  SatX:\n    consensus: ref/x.fa\n    mono: 146\n  SatY: {consensus: ref/y.fa, mono: 80}\n"),
    "missing_family_assignment": drop("family_assignment: ref/assign.tsv\n"),
    "null_run_name": sub("run_name: t1", "run_name:"),
    "family_no_consensus": sub("    consensus: ref/x.fa\n", ""),
    "family_null_block": sub("  SatY: {consensus: ref/y.fa, mono: 80}\n", "  SatY:\n"),
    "empty_file": "",
    "comments_only": "# nothing here\n",
    # map / min_items / key patterns
    "genomes_empty": sub("genomes:\n  spa: g/a.fa\n  spb: g/b.fa\n", "genomes: {}\n"),
    "genomes_list": sub("genomes:\n  spa: g/a.fa\n  spb: g/b.fa\n", "genomes: [g/a.fa]\n"),
    "code_with_space": sub("  spb: g/b.fa\n", "  'sp b': g/b.fa\n"),
    "code_with_dash": sub("  spb: g/b.fa\n", "  sp-b: g/b.fa\n"),
    "family_name_bad": sub("  SatY: {", "  Sat/Y: {"),
    "families_value_int": sub("  SatY: {consensus: ref/y.fa, mono: 80}\n", "  SatY: 5\n"),
    # types
    "mono_str": sub("mono: 146", "mono: '146'"),
    "mono_float": sub("mono: 146", "mono: 1.5"),
    "mono_bool": sub("mono: 146", "mono: true"),
    "mono_zero": sub("mono: 146", "mono: 0"),
    "mono_negative_big": sub("mono: 146", "mono: -123456789012345678901234567890"),
    "mono_list": sub("mono: 146", "mono: [1, 'a', 2.5, null, true]"),
    "mono_map": sub("mono: 146", "mono: {a: 1, b: x}"),
    "colour_bad": sub("mono: 146\n", "mono: 146\n    colour: '#xyz'\n"),
    "colour_int": sub("mono: 146\n", "mono: 146\n    colour: 123\n"),
    "run_name_space": sub("run_name: t1", "run_name: my run"),
    "run_name_int": sub("run_name: t1", "run_name: 5"),
    "run_name_newline": sub("run_name: t1", 'run_name: "a\\nb\\tc"'),
    "run_name_quote": sub("run_name: t1", "run_name: \"it's\""),
    "trf_bad": add("discovery:\n  trf: '2 5 7'\n"),
    "trf_int": add("discovery:\n  trf: 2\n"),
    "period_one": add("discovery:\n  period: [1]\n"),
    "period_reversed": add("discovery:\n  period: [300, 80]\n"),
    "period_items_str": add("discovery:\n  period: [a, b]\n"),
    "period_zero_three": add("discovery:\n  period: [0, 5, 7]\n"),
    "period_str": add("discovery:\n  period: '80 300'\n"),
    "period_null_item": add("discovery:\n  period: [300, null]\n"),
    "strand_minus": add("discovery:\n  strand: minus\n"),
    "strand_int": add("discovery:\n  strand: 1\n"),
    "strand_list": add("discovery:\n  strand: [plus]\n"),
    "cluster_id_0": add("discovery:\n  cluster_id: 0\n"),
    "cluster_id_1.5": add("discovery:\n  cluster_id: 1.5\n"),
    "cluster_id_str": add("discovery:\n  cluster_id: '0.9'\n"),
    "monomerize_quoted": add("discovery:\n  monomerize: 'yes'\n"),
    "monomerize_int": add("discovery:\n  monomerize: 1\n"),
    "monomer_min_id_excl": add("discovery:\n  monomer_min_id: 0.0\n  min_copies: 0\n  simple_repeat_cov: -1\n"),
    "min_bits_bool": add("expansion:\n  min_bits: true\n"),
    "min_bits_negative": add("expansion:\n  min_bits: -0.5\n  min_aln_frac: 1.5\n"),
    "jobs_zero": add("expansion:\n  jobs: 0\n"),
    "expansion_int": add("expansion: 5\n"),
    "expansion_list": add("expansion: [1]\n"),
    "short_segments_keep": add("arbitration:\n  short_segments: keep\n"),
    "bin_share_2": add("arbitration:\n  bin_share: 2\n  bin_min_bits: -1\n"),
    "window_mono_zero": add("terminal:\n  window_mono: 0\n"),
    "window_mono_negative": add("terminal:\n  window_mono: -1.5\n  edge_tol: 1.0\n"),
    "min_pid_101": add("orthology:\n  min_pid: 101\n"),
    "trf_timeout_59": add("finalisation:\n  trf_timeout: 59\n  min_piece_bp: 999\n"),
    "anchor_motif_bad": add("anchor:\n  motif: ACGN\n  min_identity: -0.1\n"),
    "anchor_family_bad_pattern": add("anchor:\n  family: 'a b'\n"),
    # hapdup
    "pairs_short": add("hapdup:\n  pairs: [[spa]]\n"),
    "pairs_int_item": add("hapdup:\n  pairs: [[spa, 5]]\n"),
    "pairs_str": add("hapdup:\n  pairs: spa\n"),
    "pairs_flat": add("hapdup:\n  pairs: [spa, spb]\n"),
    "pairs_bad_code": add("hapdup:\n  pairs: [[spa, sp-b], [spa, spb, spc]]\n"),
    # paths
    "path_empty": sub("out_dir: out", "out_dir: ''"),
    "path_blank": sub("out_dir: out", "out_dir: '  '"),
    "path_int": sub("out_dir: out", "out_dir: 5"),
    "path_list": sub("out_dir: out", "out_dir: [a]"),
    "genome_path_null": sub("  spb: g/b.fa\n", "  spb:\n"),
    # cross-checks (config.py:42-82)
    "members_unknown_family": add("subfamily:\n  consensuses: subx\n  members:\n    SatZ: [a]\n"),
    "members_duplicate": add("subfamily:\n  consensuses: subx\n  members:\n    SatX: [a, b, a]\n"),
    "members_other_family": add("subfamily:\n  consensuses: subx\n  members:\n    SatX: [a, SatY, SatX, SatY]\n"),
    "members_no_consensuses": add("subfamily:\n  members:\n    SatX: [a]\n"),
    "members_empty_list": add("subfamily:\n  consensuses: subx\n  members:\n    SatX: []\n"),
    "members_bad_name": add("subfamily:\n  consensuses: subx\n  members:\n    SatX: ['a b']\n"),
    "members_null_value": add("subfamily:\n  consensuses: subx\n  members:\n    SatX:\n"),
    "members_empty_map": add("subfamily:\n  members: {}\n"),
    "species_names_unknown": add("viewer:\n  species_names: {spz: Z, spa: A}\n"),
    "species_names_int": add("viewer:\n  species_names: {spa: 5}\n"),
    "anchor_family_unknown": add("anchor:\n  family: NoSuchFamily\n"),
    "tools_unknown": add("tools:\n  ssearch_36: /opt/x\n  bin_dir: /opt/bin\n  zzz: /x\n"),
    "trf_dat_unknown_code": add("inputs:\n  trf_dat:\n    spc: t/c.dat\n"),
    "trf_dat_partial_warning": add("inputs:\n  trf_dat:\n    spb: t/b.dat\n"),
    "seeds_unknown_family": add("inputs:\n  seeds:\n    SatZ: s/z.bed\n    SatX: s/x.bed\n"),
    "cross_many": add("inputs:\n  trf_dat: {spc: t.dat}\n  seeds: {SatQ: q.bed}\n"
                      "subfamily:\n  members: {SatQ: [SatX]}\nviewer:\n  species_names: {q: Q}\n"
                      "anchor: {family: Q}\ntools: {nope: /x}\ndiscovery: {period: [9, 8]}\n"),
    "schema_and_cross": add("anchor: {family: NoSuchFamily}\nexpansion: {jobs: 0}\n"),
    # top level
    "top_list": "[a, b]\n",
    "top_scalar": "hello\n",
    # valid on purpose
    "ok_home_paths": sub("out_dir: out", "out_dir: ~/x/../y") + "manual_loci: ~\nreferences: {rotated_consensus: '~root/r.fa'}\n",
    "ok_unknown_user": sub("out_dir: out", "out_dir: ~nosuchuser_satpipe/x"),
    "ok_dotdot_paths": sub("out_dir: out", "out_dir: ../../../../../../../../x/./y//z/..")
                       + "manual_loci: a\\b/../c.bed\n",
    "ok_abs_paths": sub("out_dir: out", "out_dir: //srv//x/") + "manual_loci: ///m/./n\n",
    "ok_float_int_written": add("orthology:\n  min_pid: 80\nexpansion:\n  min_bits: 50\n  min_aln_frac: 0.670\n"
                                "scan: {min_bits: 0}\narbitration: {delta: 1.0e+3, bin_min_bits: 12345678901234567890}\n"
                                "terminal: {window_mono: .5}\n"),
    "float_1e3_is_string": add("subfamily: {confident_delta: 1e3}\n"),
    "ok_duplicate_key": BASE + "run_name: t2\n",
    "ok_yes_bool": add("discovery:\n  monomerize: no\nexpansion:\n  to_scaffold_ends: yes\nanchor: {independent: On}\n"),
    "ok_yes_for_str": sub("run_name: t1", "run_name: yes"),
    "ok_nulls": add("inputs:\nexpansion:\ndiscovery: {period: null, trf: ~}\nhapdup: {pairs: [[spa, spb]], edge_tol: }\n"),
    "ok_everything": add("inputs:\n  trf_dat: {spa: t/a.dat, spb: t/b.dat}\n  seeds: {SatX: s/x.bed}\n"
                         "references: {rotated_consensus: r.fa}\n"
                         "discovery: {trf: '2 7 7 80 10 50 500', period: [50, 400], strand: both, monomerize: false}\n"
                         "subfamily:\n  consensuses: subx/\n  members:\n    SatX: [x1, x2]\n    SatY: [y1]\n"
                         "anchor: {family: SatY, motif: AAAGAAG}\ntools: {bin_dir: /opt/bin, trf: /opt/trf}\n"
                         "viewer: {title: \"T: x # y\", species_names: {spa: 'A. a', spb: \"B b'\"}}\n"
                         "manual_loci: m.bed\n"),
    "ok_unicode": sub("run_name: t1", "run_name: t1") + "viewer: {title: \"Darevskia é 中\", "
                  "species_names: {spa: \"\\u00e9\\tq\"}}\n",
    "ok_newline_key": sub("  spb: g/b.fa\n", "  \"spb\\n\": g/b.fa\n"),
    # expected differences (PORT_SPEC 3: refused instead of retyped; tests/port/golden_config.expected)
    "yaml_tab_indent": sub("  spb: g/b.fa\n", "\tspb: g/b.fa\n"),
    "xfail_octal": sub("mono: 146", "mono: 010"),
    "xfail_underscore": add("finalisation:\n  chunk_bp: 1_000\n"),
}


def main(out):
    os.makedirs(out, exist_ok=True)
    for name, text in V.items():
        with open(os.path.join(out, name + ".yaml"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
    print(f"{len(V)} variants written to {out}")


if __name__ == "__main__":
    main(sys.argv[1])
