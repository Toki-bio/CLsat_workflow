"""Config schema for satpipe: definition, validator, and Markdown rendering.

The schema is a plain nested dict so that it can be rendered to docs/config.md and
validated without any dependency beyond PyYAML. Every parameter that the stage scripts
currently hard-code appears here with the published Darevskia value as its default, so a
config that sets only the required keys reproduces the published behaviour, except where a later correction
was made the default and is said so in its field doc (discovery.monomerize; ssearch36 now always runs with fixed
statistics, -z 3, which the published run did not).
"""
from __future__ import annotations

import os
import re
from typing import Any

# ----------------------------------------------------------------------------- helpers


def F(type: str, *, required: bool = False, default: Any = None, doc: str = "", **kw) -> dict:
    """Build a field spec. Extra keys: choices, min, max, pattern, item, value, fields,
    key_pattern, length, min_items, exclusive_min."""
    d = {"type": type, "required": required, "default": default, "doc": doc}
    d.update(kw)
    return d


CODE_PATTERN = r"^[A-Za-z0-9_]+$"      # species codes end up in FASTA headers as "code|scf:s-e"
NAME_PATTERN = r"^[A-Za-z0-9_.\-]+$"   # family names end up in file names and BED names
TOOL_PATTERN = r"^[A-Za-z0-9_\-]+$"

# ----------------------------------------------------------------------------- schema

SCHEMA: dict = F("block", doc="satpipe run configuration", fields={
    "run_name": F("str", required=True, pattern=NAME_PATTERN,
                  doc="Name of the run; used in the manifest and in output file names."),
    "out_dir": F("path", required=True,
                 doc="Run directory. Every stage writes to <out_dir>/<stage>/."),
    "genomes": F("map", required=True, key_pattern=CODE_PATTERN, min_items=1, value=F("path"),
                 doc="Species/assembly code -> genome FASTA. The code is used in every locus "
                     "name (`code|scaffold:start-end`), so keep it short and free of `|`, `:` "
                     "and `-`."),
    "inputs": F("block", doc="Pre-computed inputs that let a run start after a stage.", fields={
        "trf_dat": F("map", key_pattern=CODE_PATTERN, value=F("path"),
                     doc="Code -> TRF `.dat` file for that genome. When given for every genome, "
                         "stage 02 (TRF scan) is skipped. The published Darevskia run starts here."),
        "seeds": F("map", key_pattern=NAME_PATTERN, value=F("path"),
                   doc="Family -> seed BED6 (0-based). When given for every family, stages 03-10 "
                       "are skipped and the run starts at expansion (stage 11)."),
    }),
    "families": F("map", required=True, key_pattern=NAME_PATTERN, min_items=1,
                  doc="One block per satellite family.",
                  value=F("block", fields={
                      "consensus": F("path", required=True,
                                     doc="Single-record FASTA with one monomer of the family consensus."),
                      "mono": F("int", required=True, min=1,
                                doc="Monomer length in bp. Drives every window size (3xMONO walk, "
                                    "2xMONO scan window, D3 widths, viewer step). Should equal the "
                                    "consensus length; `satpipe check` warns if it does not."),
                      "colour": F("str", pattern=r"^#[0-9A-Fa-f]{6}$",
                                  doc="Display colour for the viewer (hex). Optional."),
                  })),
    "references": F("block", doc="Reference files that feed specific stages.", fields={
        "rotated_consensus": F("path",
                               doc="Phase-rotated family consensuses prepended in Step 8 "
                                   "(`CLsat_cons.rot` in the Darevskia run). Prepending makes each "
                                   "reference win centroid election for its family. If absent the "
                                   "driver derives it by applying the Step-5 rotation to each "
                                   "`families.*.consensus` and records that deviation in the manifest."),
    }),
    "discovery": F("block", doc="Stages 02-08: TRF scan, motif extraction, filtering, clustering.", fields={
        "trf": F("str", default="2 5 7 80 10 50 2000", pattern=r"^\d+( \d+){6}$",
                 doc="TRF parameters: match mismatch delta PM PI minscore maxperiod."),
        "period": F("list", default=[80, 300], item=F("int", min=1), length=2,
                    doc="Motif period range [min, max] accepted by Step 3."),
        "min_copies": F("int", default=10, min=1, doc="Minimum TRF copy number for a motif (Step 3)."),
        "monomerize": F("bool", default=True,
                        doc="Step 3b: replace a motif that is k tandem copies of a shorter unit (TRF reporting "
                            "an array at a multiple of its unit, e.g. 292 for a 146-letter unit) by one unit. "
                            "Test: the motif against itself rotated by len/k, banded global alignment, identity "
                            ">= monomer_min_id; units never shorter than period[0]. Added 2026-10-18; the "
                            "published run had no such step (set false to reproduce it exactly)."),
        "monomer_min_id": F("float", default=0.70, min=0.0, max=1.0, exclusive_min=True,
                            doc="Step 3b: minimum self-identity (motif vs itself rotated by one unit) to call "
                                "a motif k units. Random and AT-rich random sequence of motif lengths stays "
                                "far below (0 of 2 000 called in the 2026-10-18 test)."),
        "simple_repeat_cov": F("int", default=50, min=0,
                               doc="Step 4: motifs whose own TRF re-scan covers more than this many bp "
                                   "are dropped as simple repeats."),
        "strand": F("enum", default="plus", choices=["plus", "both"],
                    doc="Step 5 canonical rotation and vsearch clustering strand. The published run "
                        "is plus-only (deliberate; it is what left L541 unclustered)."),
        "cluster_id": F("float", default=0.90, min=0.0, max=1.0, exclusive_min=True,
                        doc="vsearch identity for per-species clustering (Step 6)."),
        "cross_cluster_id": F("float", default=0.90, min=0.0, max=1.0, exclusive_min=True,
                              doc="vsearch identity for cross-species clustering (Step 8)."),
    }),
    "family_assignment": F("path", required=True,
                           doc="Curated TSV: `cluster<TAB>family[<TAB>note]`, cluster = vsearch cluster "
                               "number from Step 8. Validated at stage 09: every listed cluster must "
                               "exist, every family must have at least one cluster; unassigned clusters "
                               "are reported with sizes and best-consensus bits for the curator."),
    "expansion": F("block", doc="Stage 11: seed expansion (`clisat_expand_then_scan_merged.sh`).", fields={
        "min_bits": F("float", default=50, min=0, doc="--min-bits: minimum ssearch36 bits to accept a window."),
        "min_aln_frac": F("float", default=0.67, min=0.0, max=1.0,
                          doc="--min-aln-frac: fraction of the window (win_mono x MONO) that must align."),
        "step_mono": F("int", default=3, min=1, doc="--step-monomers: dense-walk step, in monomers."),
        "win_mono": F("int", default=3, min=1, doc="--win-monomers: probe window width, in monomers."),
        "query_repeat": F("int", default=2, min=1,
                          doc="--query-repeat N: the ssearch36 query is N+1 tandem copies of the consensus."),
        "probe_mono": F("int", default=100, min=1, doc="--probe-monomers: coarse-probe initial step, in monomers."),
        "probe_min_step": F("int", default=50, min=1, doc="--probe-min-step: coarse-probe minimum step, in monomers."),
        "max_steps": F("int", default=500, min=0, doc="--max-steps: dense-walk cap per direction (0 = unlimited)."),
        "to_scaffold_ends": F("bool", default=False, doc="--to-scaffold-ends."),
        "jobs": F("int", default=16, min=1, doc="--jobs: GNU parallel workers."),
    }),
    "scan": F("block", doc="Stage 13: symmetric rescan of MASTER loci with every family consensus.", fields={
        "win_mono": F("int", default=2, min=1, doc="Scan window width in monomers (--scan-win-bp = win_mono x MONO)."),
        "step_mono": F("int", default=1, min=1, doc="Scan step in monomers (--scan-step-bp = step_mono x MONO)."),
        "min_bits": F("float", default=0, min=0, doc="--scan-min-bits: minimum bits to emit a scan row."),
    }),
    "arbitration": F("block", doc="Stage 14: segmentation and family arbitration.", fields={
        "delta": F("float", default=20, min=0,
                   doc="DELTA: seed/rescan contradiction threshold (UNIQUE_SEED_CONTRADICTION at delta <= -DELTA)."),
        "delta_ambig": F("float", default=15, min=0, doc="DELTA_AMBIG: ambiguity flag threshold."),
        "seg_bin": F("int", default=1000, min=1, doc="SEG_BIN: D3b bin width in bp."),
        "seg_min_run": F("int", default=10, min=1, doc="SEG_MIN_RUN: minimum run of bins to start a segment."),
        "seg_min_locus": F("int", default=2000, min=1,
                           doc="SEG_MIN_LOCUS: segments shorter than this are handled per `short_segments`."),
        "short_segments": F("enum", default="drop", choices=["drop", "merge"],
                            doc="SEG_SHORT: `drop` is the published behaviour; `merge` folds short segments "
                                "into the longer neighbour."),
    }),
    "orthology": F("block", doc="Flank orthology (`rebuild_ortho_segmented.py`).", fields={
        "flank_bp": F("int", default=600, min=1, doc="FLANK: flank length extracted outward from each array edge."),
        "guard_bp": F("int", default=0, min=0, doc="GUARD: extra distance a flank must keep from another locus."),
        "min_aln": F("int", default=200, min=1, doc="Minimum aligned length for an ortholog hit."),
        "min_pid": F("float", default=80.0, min=0, max=100, doc="Minimum percent identity for an ortholog hit."),
    }),
    "subfamily": F("block", doc="Sub-subfamily assignment at array boundaries.", fields={
        "consensuses": F("path",
                         doc="Directory of per-sub-subfamily consensus FASTA files (`subx/` in the Darevskia "
                             "run) or one multi-record FASTA. Optional: without it the stage is skipped."),
        "members": F("map", key_pattern=NAME_PATTERN,
                     value=F("list", item=F("str", pattern=NAME_PATTERN), min_items=1),
                     doc="Family -> list of its sub-subfamily names (the FAMSUBS table). Each family listed "
                         "must exist in `families`."),
        "confident_delta": F("float", default=10, min=0,
                             doc="Bits margin over the second-best sub-subfamily for a confident call."),
    }),
    "terminal": F("block", doc="Terminal-monomer and boundary-phase measurement.", fields={
        "window_mono": F("float", default=1.5, exclusive_min=True, min=0,
                         doc="Window into the array from each terminus, in monomers."),
        "edge_tol": F("int", default=10, min=0, doc="EDGE_TOL: tolerated gap between alignment and array edge, bp."),
    }),
    "finalisation": F("block", doc="Optional finalisation stages (large-scaffold recovery, edge refinement, re-arbitration).",
                      fields={
        "large_scaffold_bp": F("int", default=50000000, min=1,
                               doc="Scaffolds longer than this are TRF-scanned in overlapping chunks "
                                   "(the published run needed this for one 85.8 Mb unp scaffold)."),
        "chunk_bp": F("int", default=2000000, min=1, doc="Chunk size for the chunked TRF scan."),
        "overlap_bp": F("int", default=50000, min=0, doc="Chunk overlap for the chunked TRF scan."),
        "trf_timeout": F("int", default=1800, min=60,
                         doc="Seconds after which a chunk's TRF run is stopped and the chunk bisected (the published "
                             "run's chunk 44 of the 85 Mb scaffold never finished: a telomeric tail)."),
        "min_piece_bp": F("int", default=20000, min=1000,
                          doc="Pieces shorter than this are no longer bisected; if TRF still times out on them they "
                              "are listed in unscanned.tsv and left out (never imputed)."),
        "edge_tol": F("int", default=50, min=0, doc="EDGE_TOL used by edge refinement."),
        "rearbitrate": F("bool", default=True, doc="No effect (kept so old configs load): stage 14 always arbitrates the final catalogue; there is no separate re-arbitration pass."),
    }),
    "anchor": F("block", doc="Anchor phasing (stage 21): every consensus the run writes is rotated into the frame of "
                             "the anchor family, so figures, alignments and deposited FASTA share one start.", fields={
        "family": F("str", pattern=NAME_PATTERN,
                    doc="Family whose consensus defines the frame. Default: the first family under `families`."),
        "motif": F("str", default="", pattern=r"^[ACGTacgt]*$",
                   doc="Motif the anchor frame must start with (Darevskia CLsat: AAAGAAG). Empty: the anchor "
                       "family's consensus is used exactly as given."),
        "independent": F("bool", default=False,
                         doc="Families are unrelated (no common frame exists, e.g. different satellites of one genome): "
                             "each family consensus is its own frame and only sub-subfamilies are rotated, onto their "
                             "family. Default false: all families must be placeable on the anchor."),
        "max_motif_mismatch": F("int", default=1, min=0,
                                doc="Mismatches tolerated between the motif and the start of a rotated sequence "
                                    "(diverged families, e.g. CLsat4, may differ by one base)."),
        "min_identity": F("float", default=0.6, min=0.0, max=1.0,
                          doc="Minimum identity of the local alignment to the anchor for a sequence to be placed; "
                              "below it the stage fails instead of guessing a frame."),
    }),
    "manual_loci": F("path",
                     doc="Optional BED of hand-added loci (L541 in the Darevskia catalogue): chrom start end "
                         "`<genome>|...` [score strand family]; col 4 must start with a genome code, col 7 (optional) "
                         "names the family, else the best-scoring family is used. Entries are tagged MANUAL in "
                         "provenance and listed in the manifest."),
    "tools": F("map", key_pattern=TOOL_PATTERN, value=F("path"),
               doc="Optional executable overrides: tool name -> path. The special key `bin_dir` is a "
                   "directory searched before PATH for every tool."),
    "viewer": F("block", doc="Viewer configuration written by `satpipe viewer`.", fields={
        "title": F("str", doc="Page title."),
        "species_names": F("map", key_pattern=CODE_PATTERN, value=F("str"),
                           doc="Code -> display name. Keys must exist in `genomes`."),
    }),
})


# ----------------------------------------------------------------------------- validator


class SchemaError(Exception):
    pass


def _norm_path(value: str, base_dir: str | None) -> str:
    v = os.path.expanduser(str(value))
    if base_dir and not os.path.isabs(v):
        v = os.path.join(base_dir, v)
    v = os.path.normpath(v)
    return v.replace("\\", "/")


def _check_number(val, spec, path, errors):
    lo, hi = spec.get("min"), spec.get("max")
    if lo is not None:
        if spec.get("exclusive_min") and not val > lo:
            errors.append(f"{path}: must be > {lo}, got {val}")
        elif not spec.get("exclusive_min") and val < lo:
            errors.append(f"{path}: must be >= {lo}, got {val}")
    if hi is not None and val > hi:
        errors.append(f"{path}: must be <= {hi}, got {val}")


def validate(node: Any, spec: dict, path: str, errors: list, base_dir: str | None):
    """Validate `node` against `spec`; return the normalised value (defaults filled, paths
    resolved). Errors are appended to `errors` as human-readable strings with dotted paths."""
    t = spec["type"]
    if node is None:
        if spec.get("required"):
            errors.append(f"{path}: required")
            return None
        if t == "block":
            node = {}
        else:
            return spec.get("default")

    if t == "str":
        if not isinstance(node, str):
            errors.append(f"{path}: expected a string, got {type(node).__name__}")
            return node
        pat = spec.get("pattern")
        if pat and not re.match(pat, node):
            errors.append(f"{path}: '{node}' does not match {pat}")
        return node
    if t == "enum":
        if node not in spec["choices"]:
            errors.append(f"{path}: must be one of {spec['choices']}, got {node!r}")
        return node
    if t == "bool":
        if not isinstance(node, bool):
            errors.append(f"{path}: expected true/false, got {node!r}")
        return node
    if t == "int":
        if isinstance(node, bool) or not isinstance(node, int):
            errors.append(f"{path}: expected an integer, got {node!r}")
            return node
        _check_number(node, spec, path, errors)
        return node
    if t == "float":
        if isinstance(node, bool) or not isinstance(node, (int, float)):
            errors.append(f"{path}: expected a number, got {node!r}")
            return node
        _check_number(node, spec, path, errors)
        return float(node)
    if t == "path":
        if not isinstance(node, str) or not node.strip():
            errors.append(f"{path}: expected a file path, got {node!r}")
            return node
        return _norm_path(node, base_dir)
    if t == "list":
        if not isinstance(node, list):
            errors.append(f"{path}: expected a list, got {type(node).__name__}")
            return node
        if "length" in spec and len(node) != spec["length"]:
            errors.append(f"{path}: expected exactly {spec['length']} items, got {len(node)}")
        if "min_items" in spec and len(node) < spec["min_items"]:
            errors.append(f"{path}: expected at least {spec['min_items']} items, got {len(node)}")
        return [validate(v, spec["item"], f"{path}[{i}]", errors, base_dir) for i, v in enumerate(node)]
    if t == "map":
        if not isinstance(node, dict):
            errors.append(f"{path}: expected a mapping, got {type(node).__name__}")
            return node
        if "min_items" in spec and len(node) < spec["min_items"]:
            errors.append(f"{path}: expected at least {spec['min_items']} entries, got {len(node)}")
        out = {}
        kp = spec.get("key_pattern")
        for k, v in node.items():
            k = str(k)
            if kp and not re.match(kp, k):
                errors.append(f"{path}.{k}: key does not match {kp}")
            out[k] = validate(v, spec["value"], f"{path}.{k}", errors, base_dir)
        return out
    if t == "block":
        if not isinstance(node, dict):
            errors.append(f"{path}: expected a mapping, got {type(node).__name__}")
            return node
        fields = spec["fields"]
        for k in node:
            if k not in fields:
                errors.append(f"{path + '.' if path else ''}{k}: unknown key (allowed: {', '.join(fields)})")
        out = {}
        for k, fspec in fields.items():
            out[k] = validate(node.get(k), fspec, f"{path}.{k}" if path else k, errors, base_dir)
        return out
    raise SchemaError(f"bad spec type {t!r} at {path}")


# ----------------------------------------------------------------------------- docs


def _type_str(spec: dict) -> str:
    t = spec["type"]
    if t == "enum":
        return "one of " + ", ".join(f"`{c}`" for c in spec["choices"])
    if t == "list":
        n = f" x{spec['length']}" if "length" in spec else ""
        return f"list of {spec['item']['type']}{n}"
    if t == "map":
        return f"map -> {spec['value']['type']}"
    return t


def to_markdown(spec: dict = SCHEMA) -> str:
    """Render the schema as the reference table in docs/config.md."""
    rows: list = []

    def walk(s: dict, path: str):
        t = s["type"]
        if t == "block":
            for k, f in s["fields"].items():
                walk(f, f"{path}.{k}" if path else k)
            return
        if t == "map" and s["value"]["type"] == "block":
            req = "yes" if s.get("required") else ""
            rows.append(f"| `{path}.<key>` | map of blocks | {req} | | {s['doc']} |")
            walk(s["value"], f"{path}.<key>")
            return
        req = "yes" if s.get("required") else ""
        dflt = "" if s.get("default") is None else f"`{s['default']}`"
        rows.append(f"| `{path}` | {_type_str(s)} | {req} | {dflt} | {s['doc']} |")

    walk(spec, "")
    head = ("# satpipe configuration reference\n\n"
            "Generated from `satpipe/schema.py` by `satpipe config --schema-doc`; do not edit by hand.\n\n"
            "Relative paths are resolved against the directory that holds the config file. Unknown keys "
            "are errors. Defaults are the values of the published Darevskia run.\n\n"
            "| key | type | required | default | meaning |\n|---|---|---|---|---|\n")
    return head + "\n".join(rows) + "\n"
