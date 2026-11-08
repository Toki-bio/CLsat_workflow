#!/usr/bin/env python3
"""Test helper (WP1): rebuild the schema spec from schema.tsv and compare it field by field with SCHEMA.

  python tests/port/check_schema_tsv.py SCHEMA_TSV [SCHEMA_PY]      # compare; exit 1 on any difference
  python tests/port/check_schema_tsv.py --markdown [SCHEMA_PY]      # print the Python to_markdown() of SCHEMA_PY
  python tests/port/check_schema_tsv.py --patterns SCHEMA_TSV       # print 'pattern<TAB>ere<TAB>string<TAB>0|1'
                                                                    #   (Python re.match result) for the ERE test

Independent of tools/gen_schema_tsv.py: it parses the TSV with its own code. Comparison is type-strict
(50 and 50.0 differ, True and 1 differ), so the float-default trap of PORT_SPEC 2.1 is checked.
"""
from __future__ import annotations

import ast
import importlib.util
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def load_module(path):
    if path is None:
        path = os.path.join(ROOT, "satpipe", "schema.py")
    spec = importlib.util.spec_from_file_location("_schema_under_test", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        text = fh.read()
    if "\r" in text:
        raise SystemExit(f"{path}: contains CR")
    lines = text.split("\n")
    if lines[-1] != "":
        raise SystemExit(f"{path}: no final newline")
    head = lines[0].lstrip("#").split("\t")
    return [dict(zip(head, ln.split("\t"), strict=True)) for ln in lines[1:-1]]


def rebuild(rows, problems):
    nodes = {}
    root = None
    for r in rows:
        p = r["path"]
        d = {"type": r["type"], "required": {"True": True, "False": False}[r["required"]],
             "default": json.loads(r["default_json"]), "doc": r["doc"]}
        dflt = d["default"]
        # the three default texts must agree with each other, type-strictly
        if ast.literal_eval(r["default_repr"]) != dflt or type(ast.literal_eval(r["default_repr"])) is not type(dflt):
            problems.append(f"{p}: default_repr {r['default_repr']!r} != default_json {r['default_json']!r}")
        if r["default_str"] != str(dflt):
            problems.append(f"{p}: default_str {r['default_str']!r} != str(default) {str(dflt)!r}")
        tn = {type(None): "null", bool: "bool", int: "int", float: "float", str: "str", list: "list"}[type(dflt)]
        if r["default_type"] != tn:
            problems.append(f"{p}: default_type {r['default_type']} != {tn}")
        for k in ("min", "max"):
            if r[k] != "":
                d[k] = ast.literal_eval(r[k])
        if r["exclusive_min"] == "True":
            d["exclusive_min"] = True
        elif r["exclusive_min"] != "":
            problems.append(f"{p}: exclusive_min {r['exclusive_min']!r}")
        for k in ("pattern", "key_pattern"):
            if r[k] != "":
                d[k] = r[k]
            if (r[k] == "") != (r[k + "_ere"] == ""):
                problems.append(f"{p}: {k} and {k}_ere disagree on presence")
        if r["choices"] != "":
            d["choices"] = ast.literal_eval(r["choices"])
        for k in ("length", "min_items"):
            if r[k] != "":
                d[k] = int(r[k])
        exp_origin = "derived" if p == "anchor/family" else "default"
        if r["absent_origin"] != exp_origin:
            problems.append(f"{p}: absent_origin {r['absent_origin']} != {exp_origin}")
        if d["type"] == "block":
            d["fields"] = {}
        nodes[p] = (d, r["elem_type"])
        if p == ".":
            root = d
            continue
        parent, _, last = p.rpartition("/")
        parent = parent or "."
        if parent not in nodes:
            problems.append(f"{p}: parent row {parent} missing or after the child")
            continue
        pd = nodes[parent][0]
        if last == "<key>":
            pd["value"] = d
        elif last == "<item>":
            pd["item"] = d
        else:
            pd["fields"][last] = d
    for p, (d, elem) in nodes.items():
        want = d["item"]["type"] if d["type"] == "list" else d["value"]["type"] if d["type"] == "map" else ""
        if elem != want:
            problems.append(f"{p}: elem_type {elem!r} != {want!r}")
    return root


def strict_diff(a, b, path, problems):
    """Type-strict recursive comparison; dict key ORDER also compared (fields are walked in order)."""
    if type(a) is not type(b):
        problems.append(f"{path}: type {type(a).__name__} ({a!r}) != {type(b).__name__} ({b!r})")
        return
    if isinstance(a, dict):
        if list(a) != list(b):
            problems.append(f"{path}: keys {list(a)} != {list(b)}")
        for k in a:
            if k in b:
                strict_diff(a[k], b[k], f"{path}.{k}", problems)
    elif isinstance(a, (list, tuple)):
        if len(a) != len(b):
            problems.append(f"{path}: length {len(a)} != {len(b)}")
        for i, (x, y) in enumerate(zip(a, b)):
            strict_diff(x, y, f"{path}[{i}]", problems)
    elif a != b:
        problems.append(f"{path}: {a!r} != {b!r}")


def canon(spec):
    """Python SCHEMA with dict keys in a fixed order (F() order differs from the rebuild's), for strict_diff."""
    order = ["type", "required", "default", "doc", "min", "max", "exclusive_min", "pattern", "key_pattern",
             "choices", "length", "min_items", "fields", "value", "item"]
    out = {}
    for k in order:
        if k in spec:
            v = spec[k]
            if k == "fields":
                v = {fk: canon(fv) for fk, fv in v.items()}
            elif k in ("value", "item"):
                v = canon(v)
            out[k] = v
    extra = set(spec) - set(order)
    if extra:
        raise SystemExit(f"spec keys not covered by the test: {sorted(extra)}")
    return out


PATTERN_STRINGS = ["", "a", "abc", "A_1", "a-b", "a.b", "a b", "a|b", "a:b", "a/b", "-", ".", "_", "<key>", "x\\y",
                   "#12ab3F", "#12AB3f", "#12ab3G", "#12ab3", "12ab3F", "#12ab3F0", "2 5 7 80 10 50 2000",
                   "2 5 7 80 10 50", "2 5 7 80 10 50 2000 1", "2  5 7 80 10 50 2000", " 2 5 7 80 10 50 2000",
                   "ACGT", "acgtACGT", "ACGTN", "AAAGAAG", "CLsat1", "CLsat1.2", "unif", "trf", "bin_dir", "mafft-7"]


def main(argv):
    if argv and argv[0] == "--markdown":
        sys.stdout.buffer.write(load_module(argv[1] if len(argv) > 1 else None).to_markdown().encode("utf-8"))
        return 0
    if argv and argv[0] == "--patterns":
        pats = set()
        for r in read_tsv(argv[1]):
            for k in ("pattern", "key_pattern"):
                if r[k]:
                    pats.add((r[k], r[k + "_ere"]))
        out = []
        for pat, ere in sorted(pats):
            for s in PATTERN_STRINGS:
                out.append(f"{pat}\t{ere}\t{s}\t{1 if re.match(pat, s) else 0}")
        sys.stdout.buffer.write(("\n".join(out) + "\n").encode("utf-8"))
        return 0
    tsv, schema_py = argv[0], (argv[1] if len(argv) > 1 else None)
    problems = []
    rows = read_tsv(tsv)
    root = rebuild(rows, problems)
    strict_diff(canon(root), canon(load_module(schema_py).SCHEMA), "SCHEMA", problems)
    n_leaf = sum(1 for r in rows if r["type"] not in ("block",))
    for p in problems:
        print("DIFF", p)
    print(f"schema.tsv rows {len(rows)} (non-block {n_leaf}); differences {len(problems)}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
