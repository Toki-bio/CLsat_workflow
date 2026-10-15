"""satpipe command line: check | run | viewer | stages | config."""
from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .check import run_check
from .config import ConfigError, dump_config, load_config
from .schema import to_markdown
from .stages import STAGES

NOT_YET = ("not implemented yet: the driver is work-plan step 5 (HANDOFF.md). "
           "`satpipe check` and `satpipe config` work; `satpipe stages` lists the contracts.")


def cmd_check(a: argparse.Namespace) -> int:
    rep = run_check(a.config, do_files=not a.no_files, do_tools=not a.no_tools)
    if a.json:
        print(json.dumps(rep.to_json(), indent=2))
    else:
        print(rep.to_text())
    return 1 if rep.failed else 0


def cmd_config(a: argparse.Namespace) -> int:
    if a.schema_doc:
        sys.stdout.write(to_markdown())
        return 0
    if not a.config:
        print("error: config path required unless --schema-doc", file=sys.stderr)
        return 2
    try:
        cfg = load_config(a.config)
    except ConfigError as e:
        for err in e.errors:
            print(f"error: {err}", file=sys.stderr)
        return 1
    sys.stdout.write(dump_config(cfg))
    return 0


def cmd_stages(a: argparse.Namespace) -> int:
    print(f"{'stage':12} {'wired':5} {'name':46} scripts")
    for st in STAGES:
        scripts = ", ".join(st.scripts) or "(command in docs/files.md)"
        print(f"{st.id:12} {'yes' if st.wired else 'no':5} {st.name:46} {scripts}")
        notes = []
        if st.optional:
            notes.append("optional")
        if st.hand:
            notes.append(f"hand-curated input: {st.hand}")
        if st.needs:
            notes.append("needs " + ", ".join(st.needs))
        if notes:
            print(f"{'':12} {'':5} {'':46} [{'; '.join(notes)}]")
    return 0


def cmd_not_yet(a: argparse.Namespace) -> int:
    print(f"satpipe {a.cmd}: {NOT_YET}", file=sys.stderr)
    return 2


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="satpipe", description="Reusable satellite-DNA locus pipeline.")
    p.add_argument("--version", action="version", version=f"satpipe {__version__}")
    sp = p.add_subparsers(dest="cmd", required=True)

    c = sp.add_parser("check", help="validate a config: schema, input files, external tools")
    c.add_argument("config")
    c.add_argument("--no-files", action="store_true", help="skip input-file checks (config authoring on another machine)")
    c.add_argument("--no-tools", action="store_true", help="skip external-tool checks")
    c.add_argument("--json", action="store_true")
    c.set_defaults(func=cmd_check)

    g = sp.add_parser("config", help="print the normalised config (defaults filled, paths resolved) or the schema docs")
    g.add_argument("config", nargs="?")
    g.add_argument("--schema-doc", action="store_true", help="print docs/config.md content")
    g.set_defaults(func=cmd_config)

    s = sp.add_parser("stages", help="list stage contracts and their wiring status")
    s.set_defaults(func=cmd_stages)

    r = sp.add_parser("run", help="run stages (not implemented yet)")
    r.add_argument("config")
    r.add_argument("--from", dest="from_stage")
    r.add_argument("--to", dest="to_stage")
    r.set_defaults(func=cmd_not_yet)

    v = sp.add_parser("viewer", help="assemble the viewer directory for a run (not implemented yet)")
    v.add_argument("run_dir")
    v.set_defaults(func=cmd_not_yet)
    return p


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
