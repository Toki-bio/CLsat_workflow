#!/usr/bin/env python3
"""Python oracle for the bash port's config loader (PORT_SPEC 5a, WP3). Test-only; deleted at switch-over.

    python3 tools/py_oracle.py flat CONFIG.yaml       # = config_stream (lib/satpipe/config.sh)
    python3 tools/py_oracle.py reload CONFIG.yaml DUMP.yaml   # exit 0 when yaml.safe_load(DUMP) == load_config(CONFIG).data
    python3 tools/py_oracle.py tools                  # names of satpipe/tools.py TOOLS, space separated
    env SATPIPE_SCHEMA_PY=path/to/schema.py           # validate with another schema.py (e.g. another branch's)

flat: load_config(CONFIG) -> TYPED flat (PORT_SPEC 2.1) on stdout, then one `!warning<TAB>msg` line per warning;
exit 0. On ConfigError: one `!error<TAB>msg` line per error, exit 1. OSError: `!error<TAB>str(e)`, exit 1.
  type   path for schema path fields, str for str/enum, map for map/block, else the value's Python type
         (int bool float null list); a float field absent from the YAML stays the int default (schema.py:270)
  text   str(value) (float repr), containers: len; text and path components escaped \\ \\t \\n \\r
  origin file when the key is in the YAML with a non-null value, else the schema's absent origin
         (derived for anchor/family, default otherwise)"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import yaml  # noqa: E402

from satpipe import config as C  # noqa: E402

DERIVED = {"anchor/family"}
MISSING = object()

if os.environ.get("SATPIPE_SCHEMA_PY"):
    _spec = importlib.util.spec_from_file_location("_oracle_schema", os.environ["SATPIPE_SCHEMA_PY"])
    _mod = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_mod)
    C.SCHEMA = _mod.SCHEMA
    C.validate = _mod.validate


def esc(s):
    return s.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n").replace("\r", "\\r")


def tname(v):
    if v is None:
        return "null"
    for t, n in ((bool, "bool"), (int, "int"), (float, "float"), (str, "str"), (list, "list"), (dict, "map")):
        if isinstance(v, t):
            return n
    raise TypeError(type(v).__name__)


def walk(val, raw, spec, sp, path, out):
    origin = "file" if raw is not MISSING and raw is not None else ("derived" if sp in DERIVED else "default")
    t = tname(val)
    if spec["type"] == "path" and t == "str":
        t = "path"
    if t in ("map", "list"):
        out.append(f"{path}\t{t}\t{len(val)}\t{origin}")
        raw_ok = raw is not MISSING and raw is not None
        if t == "list":
            for i, v in enumerate(val):
                r = raw[i] if raw_ok and isinstance(raw, list) and i < len(raw) else MISSING
                if spec["type"] == "list":
                    walk(v, r, spec["item"], f"{sp}/<item>", f"{path}/{i}", out)
                else:   # a default list object (discovery.period)
                    walk(v, MISSING, {"type": tname(v)}, sp, f"{path}/{i}", out)
            return
        rmap = {}
        if raw_ok and isinstance(raw, dict):
            for k, v in raw.items():
                rmap[str(k)] = v
        for k, v in val.items():
            child = f"{path}/{esc(k)}" if path else esc(k)
            if spec["type"] == "block":
                walk(v, rmap.get(k, MISSING), spec["fields"][k], f"{sp}/{k}" if sp else k, child, out)
            else:
                walk(v, rmap.get(k, MISSING), spec["value"], f"{sp}/<key>", child, out)
        return
    text = repr(val) if t == "float" else str(val)
    out.append(f"{path}\t{t}\t{esc(text)}\t{origin}")


def flat(cfg_path):
    try:
        cfg = C.load_config(cfg_path)
    except C.ConfigError as e:
        return ["!error\t" + esc(m) for m in e.errors], 1
    except OSError as e:
        return ["!error\t" + esc(str(e))], 1
    with open(cfg.path, encoding="utf-8") as fh:
        raw = yaml.safe_load(fh) or {}
    out = []
    walk(cfg.data, raw, C.SCHEMA, "", "", out)
    out = out[1:]                     # the root row is implicit
    out += ["!warning\t" + esc(w) for w in cfg.warnings]
    return out, 0


def typed_equal(a, b):
    if type(a) is not type(b):
        return False
    if isinstance(a, dict):
        return list(a) == list(b) and all(typed_equal(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(typed_equal(x, y) for x, y in zip(a, b))
    return a == b or (a != a and b != b)


def main(argv):
    if argv[0] == "flat":
        out, rc = flat(argv[1])
        sys.stdout.buffer.write("".join(x + "\n" for x in out).encode("utf-8", "surrogateescape"))
        return rc
    if argv[0] == "reload":
        cfg = C.load_config(argv[1])
        with open(argv[2], encoding="utf-8") as fh:
            back = yaml.safe_load(fh)
        if typed_equal(back, cfg.data):
            return 0
        print(f"reload differs:\n  loaded: {cfg.data!r}\n  dumped: {back!r}")
        return 1
    if argv[0] == "tools":
        from satpipe.tools import TOOLS
        print(" ".join(t.name for t in TOOLS))
        return 0
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
