#!/usr/bin/env python3
"""Write the test-generated satpipe configs as YAML fixtures for the yaml.awk differential (WP2).

    python3 tests/port/make_yaml_fixtures.py OUTDIR

Sources (read from this checkout, nothing hand-typed except the mutation strings copied from the tests):
  conftest.py `project` fixture config (dedented, as the fixture writes it) and the edits that test_schema.py and
  test_check.py apply to it; the heredoc configs of tests/toy/run_driver_toy.sh and tests/molly/run_molly.sh,
  expanded by bash with placeholder variables (defaults and the variants the spec names: MONOMERIZE, STRAND, LARGE)."""
import os
import re
import subprocess
import sys
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))


def read(p):
    with open(os.path.join(REPO, p), encoding="utf-8") as fh:
        return fh.read()


def write(out, name, text):
    with open(os.path.join(out, name), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)


def main(out):
    os.makedirs(out, exist_ok=True)
    m = re.search(r'write\(d / "configs" / "run.yaml", """(.*?)"""\)', read("tests/conftest.py"), re.S)
    base = textwrap.dedent(m.group(1))
    write(out, "conftest_base.yaml", base)
    # edits copied from tests/test_schema.py and tests/test_check.py (line numbers of this branch)
    edits = {
        "schema_36_unknown_nested.yaml": base + "\nexpansion:\n  min_bitz: 3\n",
        "schema_59_cross.yaml": base + "\nsubfamily:\n  members:\n    SatY: [a]\n  consensuses: ../ref\n"
                                       "viewer:\n  species_names:\n    zzz: Z\n",
        "schema_70_trfdat_code.yaml": base.replace("spb: ../trf/spb.dat", "spc: ../trf/spb.dat"),
        "schema_115_hapdup_same.yaml": base + "\nhapdup:\n  pairs: [[spa, spa]]\n",
        "schema_124_hapdup_short.yaml": base + "\nhapdup:\n  pairs: [[spa]]\n",
        "check_26_mono.yaml": base.replace("mono: 148", "mono: 146"),
        "check_42_family.yaml": base.replace("families:\n", "families:\n  SatY: {consensus: ../ref/SatX.cons.fa, mono: 148}\n"),
        "check_82_bindir.yaml": base + "\ntools:\n  bin_dir: /tmp/pytest-of-x/fakebin\n",
        "check_125_anchor.yaml": base + "\nanchor: {family: NoSuchFamily}\n",
        "check_129_tools.yaml": base + "\ntools: {ssearch_36: /opt/x/ssearch36}\n",
    }
    edits["schema_75_trfdat_partial.yaml"] = edits["schema_70_trfdat_code.yaml"].replace("    spc: ../trf/spb.dat\n", "")
    for k, v in edits.items():
        write(out, k, v)

    env = dict(os.environ, W="/w/toy work", TOY="/data/toy", REPO="/repo/satpipe", TAB="/w/cluster_family.tsv",
               M="/data/molly", JOBS="12")
    for script, name, variants in (
            ("tests/toy/run_driver_toy.sh", "toy", [{}, {"MONOMERIZE": "false", "STRAND": "both", "LARGE": "200000"}]),
            ("tests/molly/run_molly.sh", "molly", [{}])):
        src = read(script)
        hd = re.search(r'cat > "\$W/\w+\.yaml" <<EOF\n(.*?)\nEOF\n', src, re.S).group(1)
        for i, var in enumerate(variants):
            sh = "mono() { echo 146; }\ncat <<EOF\n" + hd + "\nEOF\n"
            text = subprocess.run(["bash", "-c", sh], env=dict(env, **var), capture_output=True, text=True,
                                  encoding="utf-8", check=True).stdout
            write(out, f"{name}{i}.yaml", text)


if __name__ == "__main__":
    main(sys.argv[1])
