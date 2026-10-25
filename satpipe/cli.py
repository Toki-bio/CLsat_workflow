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


def cmd_viewer(a: argparse.Namespace) -> int:
    from .run import StageError
    from .viewer import build
    try:
        build(a.run_dir, out_dir=a.out, viewer_html=a.viewer_html, flank_scan=not a.no_flank_scan)
    except (StageError, FileNotFoundError) as e:
        print(f"satpipe viewer: {e}", file=sys.stderr)
        return 1
    return 0


def cmd_run(a: argparse.Namespace) -> int:
    from .run import StageError, run, stage_ids
    if a.list_stages:
        print("\n".join(stage_ids()))
        return 0
    rep = run_check(a.config, do_files=True, do_tools=not a.skip_tool_check)
    if rep.failed:
        print(rep.to_text(), file=sys.stderr)
        print("satpipe run: refusing to start, fix the check failures above (or --skip-tool-check)", file=sys.stderr)
        return 1
    try:
        manifest = run(rep.config, from_stage=a.from_stage, to_stage=a.to_stage, force=a.force,
                       only=a.only.split(",") if a.only else None)
    except StageError as e:
        print(f"satpipe run: {e}", file=sys.stderr)
        return 1
    keys = ("centroids_n", "centroids_md5", "seed_counts", "master_loci", "segmented_loci", "catalogue_loci", "ortholog_pairs")
    print("summary: " + "; ".join(f"{k}={manifest[k]}" for k in keys if k in manifest))
    return 0


def cmd_anchor(a: argparse.Namespace) -> int:
    from . import anchor as A
    rep = run_check(a.config, do_files=True, do_tools=False)
    if rep.failed:
        print(rep.to_text(), file=sys.stderr)
        return 1
    d = rep.config.data
    ac = d.get("anchor") or {}
    fams = d["families"]
    afam = ac.get("family") or next(iter(fams))
    motif, mm, mi = ac.get("motif") or "", ac.get("max_motif_mismatch", 1), ac.get("min_identity", 0.6)
    try:
        frame = A.anchor_frame("".join(A.read_fasta(fams[afam]["consensus"]).values()), motif, mm)
        bad = 0
        for path in a.fasta:
            for fp in A.fasta_files(path):
                recs = A.read_fasta(fp)
                if a.check:
                    probs = A.verify_records(recs, frame, motif, mm, mi)
                    for p in probs:
                        print(f"{fp}: {p}", file=sys.stderr)
                    bad += len(probs)
                    print(f"{fp}: {len(recs) - len(probs)}/{len(recs)} records in the anchor frame")
                else:
                    if any("-" in s for s in recs.values()):
                        print(f"{fp}: gapped sequences (an alignment): rotate the unaligned sequences and re-align", file=sys.stderr)
                        return 1
                    rot, _ = A.anchor_records(recs, frame, motif, mm, mi)
                    out = a.out if (a.out and len(a.fasta) == 1 and not __import__("os").path.isdir(path)) else fp + ".anchored.fa"
                    A.write_fasta(out, rot)
                    print(f"{fp} -> {out}")
    except A.AnchorError as e:
        print(f"satpipe anchor: {e}", file=sys.stderr)
        return 1
    return 1 if bad else 0


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

    r = sp.add_parser("run", help="run the stages (resumes: stages whose DONE marker matches are skipped)")
    r.add_argument("config")
    r.add_argument("--from", dest="from_stage", help="first stage to run (see --list-stages)")
    r.add_argument("--to", dest="to_stage", help="last stage to run")
    r.add_argument("--only", help="comma-separated stage ids to run")
    r.add_argument("--force", action="store_true", help="re-run selected stages even if their DONE marker matches")
    r.add_argument("--skip-tool-check", action="store_true")
    r.add_argument("--list-stages", action="store_true")
    r.set_defaults(func=cmd_run)

    v = sp.add_parser("viewer", help="assemble the viewer directory (data files + viewer_config.json + page) for a finished run")
    v.add_argument("run_dir")
    v.add_argument("--out", help="output directory [<run_dir>/viewer]")
    v.add_argument("--viewer-html", help="index.html to copy in [the copy vendored with satpipe]")
    v.add_argument("--no-flank-scan", action="store_true", help="skip the +-5 kb flank background scan (needs ssearch36)")
    v.set_defaults(func=cmd_viewer)

    n = sp.add_parser("anchor", help="rotate FASTA files into the run's anchor frame, or verify that they are in it")
    n.add_argument("config")
    n.add_argument("fasta", nargs="+", help="FASTA file(s) or directories of them")
    n.add_argument("--check", action="store_true", help="only verify; exit 1 if any record is out of frame (works on alignments)")
    n.add_argument("-o", "--out", help="output file (single input only) [<input>.anchored.fa]")
    n.set_defaults(func=cmd_anchor)
    return p


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
