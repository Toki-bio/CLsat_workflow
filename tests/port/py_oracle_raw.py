#!/usr/bin/env python3
"""PyYAML oracle for lib/satpipe/yaml.awk (PORT_SPEC.md section 6, WP2): print the RAW flat form of a YAML file
exactly as yaml.awk documents it (header of lib/satpipe/yaml.awk), using yaml.safe_load as config.py:93 does.

    python3 tests/port/py_oracle_raw.py CONFIG.yaml

Errors: one line '!error<TAB>...' and exit 1. 'YAML parse error: ...' for PyYAML errors (config.py:95 wording),
'top level must be a mapping' (config.py:99), 'unsupported value type X' for values the RAW form cannot carry
(datetime from a timestamp; the reader refuses those)."""
import sys

import yaml


def esc_text(s):
    return s.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n").replace("\r", "\\r")


def esc_key(s):
    return esc_text(s).replace("/", "\\x2f")


def rows(node, prefix, out):
    items = node.items() if isinstance(node, dict) else enumerate(node)
    for k, v in items:
        comp = esc_key(str(k)) if isinstance(node, dict) else str(k)
        p = f"{prefix}/{comp}" if prefix is not None else comp
        if isinstance(v, dict):
            out.append(f"{p}\tmap\t{len(v)}")
            rows(v, p, out)
        elif isinstance(v, list):
            out.append(f"{p}\tlist\t{len(v)}")
            rows(v, p, out)
        elif isinstance(v, bool):
            out.append(f"{p}\tbool\t{v}")
        elif isinstance(v, int):
            out.append(f"{p}\tint\t{v}")
        elif isinstance(v, float):
            out.append(f"{p}\tfloat\t{v!r}")
        elif v is None:
            out.append(f"{p}\tnull\tNone")
        elif isinstance(v, str):
            out.append(f"{p}\tstr\t{esc_text(v)}")
        else:
            raise TypeError(type(v).__name__)


def main(path):
    out = []
    try:
        with open(path, encoding="utf-8") as fh:
            raw = yaml.safe_load(fh)
        if raw is None:
            raw = {}
        if not isinstance(raw, dict):
            out = ["!error\ttop level must be a mapping"]
        else:
            rows(raw, None, out)
    except yaml.YAMLError as e:
        out = ["!error\tYAML parse error: " + " ".join(str(e).split())]
    except TypeError as e:
        out = [f"!error\tunsupported value type {e}"]
    except UnicodeDecodeError as e:
        out = [f"!error\tUnicodeDecodeError {e}"]
    sys.stdout.buffer.write(("".join(line + "\n" for line in out)).encode("utf-8", "backslashreplace"))
    return 1 if out and out[0].startswith("!error") else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
