#!/usr/bin/env python3
"""Generate lib/satpipe/schema.tsv (and the two stage tables) from the Python driver.

One-time generator for the bash port (PORT_SPEC.md section 1, WP1). schema.tsv is GENERATED: never edit it by
hand while satpipe/schema.py exists; rerun this script whenever schema.py changes.

  python tools/gen_schema_tsv.py schema   [--schema PATH] [-o FILE]   # default -o lib/satpipe/schema.tsv
  python tools/gen_schema_tsv.py stages   [-o FILE]                   # default -o lib/satpipe/stages_registry.tsv
  python tools/gen_schema_tsv.py runorder [-o FILE]                   # default -o lib/satpipe/run_stages.tsv

--schema loads SCHEMA from any schema.py file (e.g. another branch's, extracted with `git show`), so the
generator can be checked against schema versions that are not checked out.

schema.tsv: one header line starting with '#', then one row per schema node in definition order (depth first,
the order validate() and to_markdown() walk). Columns (tab separated, LF line ends, ASCII):

  path           '.' for the root block; else components joined by '/'; '<key>' = the value spec of a map,
                 '<item>' = the item spec of a list (e.g. families/<key>/mono, hapdup/pairs/<item>/<item>)
  type           block map list str int float bool enum path
  required       True / False
  elem_type      list: item type; map: value type; else empty
  default_type   Python type of the default OBJECT: null int float str bool list (a float field with default 50
                 has default_type int: absent from the YAML, validate() returns the int untouched, schema.py:270)
  default_json   json.dumps(default)              (null when no default)
  default_str    str(default)                     (None when no default; '' default -> empty cell)
  default_repr   repr(default)
  min max        repr() of the bound ('0' vs '0.0' matters in messages); empty when absent
  exclusive_min  True when set; empty when absent
  pattern        the Python regex as written (error messages quote it)
  pattern_ere    the same regex for gawk: \\d -> [0-9], '\\-' inside a bracket -> '-' last; no backslash left
  key_pattern    map key regex as written;  key_pattern_ere as above
  choices        Python list repr (['plus', 'both'])
  length         list length;  min_items: map/list minimum
  absent_origin  origin of the TYPED-flat row when the key is absent (PORT_SPEC 2.1): default, or derived for the
                 fields the driver fills at run time and the schema leaves null (anchor/family, run.py:1046)
  doc            the doc string (asserted free of tab, CR, LF)
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

COLUMNS = ["path", "type", "required", "elem_type", "default_type", "default_json", "default_str", "default_repr",
           "min", "max", "exclusive_min", "pattern", "pattern_ere", "key_pattern", "key_pattern_ere", "choices",
           "length", "min_items", "absent_origin", "doc"]

# PORT_SPEC 2.1: anchor/family is not filled by the schema; run.py:1046 falls back to the first family at run time.
DERIVED = {"anchor/family"}

KNOWN_KEYS = {"type", "required", "default", "doc", "choices", "min", "max", "pattern", "item", "value", "fields",
              "key_pattern", "length", "min_items", "exclusive_min"}


def load_schema(path: str | None):
    if path is None:
        sys.path.insert(0, ROOT)
        from satpipe import schema as mod
        return mod.SCHEMA
    spec = importlib.util.spec_from_file_location("_satpipe_schema_src", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.SCHEMA


def type_name(v) -> str:
    if v is None:
        return "null"
    if isinstance(v, bool):
        return "bool"
    for t, n in ((int, "int"), (float, "float"), (str, "str"), (list, "list")):
        if isinstance(v, t):
            return n
    raise SystemExit(f"gen_schema_tsv: unsupported default type {type(v).__name__}")


def to_ere(pat: str) -> str:
    """Python regex (the subset schema.py uses) -> gawk ERE with the same meaning."""
    e = pat.replace(r"\d", "[0-9]")
    e = re.sub(r"\[([^\]]*)\\-\]", r"[\1-]", e)     # [..\-] -> [..-]  (dash last = literal in ERE)
    if "\\" in e:
        raise SystemExit(f"gen_schema_tsv: pattern {pat!r} needs a manual ERE translation")
    return e


def cell(s: str) -> str:
    if any(c in s for c in "\t\r\n"):
        raise SystemExit(f"gen_schema_tsv: tab/CR/LF in a cell: {s!r}")
    if any(ord(c) > 126 for c in s):
        raise SystemExit(f"gen_schema_tsv: non-ASCII in a cell: {s!r}")
    return s


def rows_of(spec: dict, path: str, out: list):
    unknown = set(spec) - KNOWN_KEYS
    if unknown:
        raise SystemExit(f"gen_schema_tsv: unknown spec keys {sorted(unknown)} at {path}")
    t = spec["type"]
    d = spec.get("default")
    if "exclusive_min" in spec and spec["exclusive_min"] is not True:
        raise SystemExit(f"gen_schema_tsv: exclusive_min must be True when given ({path})")
    elem = spec["item"]["type"] if t == "list" else spec["value"]["type"] if t == "map" else ""
    r = {
        "path": path or ".",
        "type": t,
        "required": repr(bool(spec.get("required"))),
        "elem_type": elem,
        "default_type": type_name(d),
        "default_json": json.dumps(d),
        "default_str": str(d),
        "default_repr": repr(d),
        "min": repr(spec["min"]) if spec.get("min") is not None else "",
        "max": repr(spec["max"]) if spec.get("max") is not None else "",
        "exclusive_min": "True" if spec.get("exclusive_min") else "",
        "pattern": spec.get("pattern") or "",
        "pattern_ere": to_ere(spec["pattern"]) if spec.get("pattern") else "",
        "key_pattern": spec.get("key_pattern") or "",
        "key_pattern_ere": to_ere(spec["key_pattern"]) if spec.get("key_pattern") else "",
        "choices": repr(spec["choices"]) if "choices" in spec else "",
        "length": str(spec["length"]) if "length" in spec else "",
        "min_items": str(spec["min_items"]) if "min_items" in spec else "",
        "absent_origin": "derived" if path in DERIVED else "default",
        "doc": spec.get("doc", ""),
    }
    out.append("\t".join(cell(r[c]) for c in COLUMNS))
    sub = (lambda k: f"{path}/{k}" if path else k)
    if t == "block":
        for k, f in spec["fields"].items():
            if "/" in k or k.startswith("<"):
                raise SystemExit(f"gen_schema_tsv: field name {k!r} cannot be a path component")
            rows_of(f, sub(k), out)
    elif t == "map":
        rows_of(spec["value"], sub("<key>"), out)
    elif t == "list":
        rows_of(spec["item"], sub("<item>"), out)


def gen_schema(schema_path: str | None) -> str:
    rows: list = []
    rows_of(load_schema(schema_path), "", rows)
    return "#" + "\t".join(COLUMNS) + "\n" + "\n".join(rows) + "\n"


def gen_stages() -> str:
    """The static contract table of `satpipe stages` (satpipe/stages/__init__.py STAGES, cli.py:42-56)."""
    sys.path.insert(0, ROOT)
    from satpipe.stages import STAGES
    cols = ["id", "name", "wired", "optional", "hand", "needs", "scripts"]
    out = ["#" + "\t".join(cols)]
    for st in STAGES:
        vals = [st.id, st.name, repr(st.wired), repr(st.optional), st.hand, ",".join(st.needs), ",".join(st.scripts)]
        if any("," in s for s in st.needs + st.scripts):
            raise SystemExit("gen_schema_tsv: ',' inside a needs/scripts item")
        out.append("\t".join(cell(v) for v in vals))
    return "\n".join(out) + "\n"


def gen_runorder() -> str:
    """The run order of `satpipe run` (run.py @stage decorators: ids in run order, driver-declared needs)."""
    sys.path.insert(0, ROOT)
    from satpipe.run import STAGES
    out = ["#id\tneeds"]
    for fn in STAGES:
        out.append(cell(fn.sid) + "\t" + cell(",".join(fn.needs)))
    return "\n".join(out) + "\n"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("what", choices=["schema", "stages", "runorder"])
    p.add_argument("--schema", help="schema.py to load [satpipe/schema.py of this checkout]")
    p.add_argument("-o", "--out", help="output file ('-' = stdout)")
    a = p.parse_args(argv)
    if a.what == "schema":
        text, dflt = gen_schema(a.schema), "lib/satpipe/schema.tsv"
    elif a.what == "stages":
        text, dflt = gen_stages(), "lib/satpipe/stages_registry.tsv"
    else:
        text, dflt = gen_runorder(), "lib/satpipe/run_stages.tsv"
    out = a.out or os.path.join(ROOT, dflt)
    if out == "-":
        sys.stdout.buffer.write(text.encode("ascii"))
    else:
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
        with open(out, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
