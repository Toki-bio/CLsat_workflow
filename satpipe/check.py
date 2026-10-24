"""`satpipe check`: config schema, input files, external tools."""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Optional

from .config import Config, ConfigError, load_config
from .stages import STAGES, missing_scripts
from .tools import probe_all

OK, WARN, FAIL = "OK", "WARN", "FAIL"


@dataclass
class Item:
    level: str
    section: str
    message: str


@dataclass
class Report:
    items: list = field(default_factory=list)
    config: Optional[Config] = None

    def add(self, level: str, section: str, message: str) -> None:
        self.items.append(Item(level, section, message))

    @property
    def failed(self) -> bool:
        return any(i.level == FAIL for i in self.items)

    def counts(self) -> dict:
        return {lv: sum(1 for i in self.items if i.level == lv) for lv in (OK, WARN, FAIL)}

    def to_text(self) -> str:
        lines = [f"[{it.level:4}] {it.section}: {it.message}" for it in self.items]
        c = self.counts()
        verdict = "FAILED" if self.failed else "passed"
        lines.append(f"check {verdict}: {c[OK]} ok, {c[WARN]} warnings, {c[FAIL]} failures")
        return "\n".join(lines)

    def to_json(self) -> dict:
        return {"failed": self.failed, "counts": self.counts(),
                "items": [{"level": i.level, "section": i.section, "message": i.message} for i in self.items]}


# ----------------------------------------------------------------------------- file helpers


def read_fasta(path: str, max_records: int = 100000) -> list:
    """Return [(name, length)] without holding sequences; stops after max_records."""
    out: list = []
    name, n = None, 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, n))
                    if len(out) >= max_records:
                        return out
                hdr = line[1:].strip()
                name, n = (hdr.split()[0] if hdr else ""), 0
            elif name is not None:
                n += len(line.strip())
    if name is not None:
        out.append((name, n))
    return out


def _exists(report: Report, section: str, label: str, path: str, kind: str = "file", fatal: bool = True) -> bool:
    ok = os.path.isdir(path) if kind == "dir" else os.path.isfile(path)
    if ok:
        report.add(OK, section, f"{label}: {path}")
        return True
    report.add(FAIL if fatal else WARN, section, f"{label}: {kind} not found: {path}")
    return False


def check_family_assignment(path: str, families: set, report: Report) -> None:
    sec = "family_assignment"
    if not _exists(report, sec, "table", path):
        return
    seen = {f: set() for f in families}
    clusters: dict = {}
    nrows = 0
    with open(path, encoding="utf-8") as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                report.add(FAIL, sec, f"line {ln}: need `cluster<TAB>family`, got {line.strip()!r}")
                continue
            clu, fam = p[0].strip(), p[1].strip()
            if not re.fullmatch(r"\d+", clu):
                report.add(FAIL, sec, f"line {ln}: cluster must be a non-negative integer, got {clu!r}")
                continue
            if fam not in families:
                report.add(FAIL, sec, f"line {ln}: family {fam!r} is not in `families` ({', '.join(sorted(families))})")
                continue
            if clu in clusters:
                report.add(FAIL, sec, f"line {ln}: cluster {clu} listed twice ({clusters[clu]} and {fam})")
                continue
            clusters[clu] = fam
            seen[fam].add(clu)
            nrows += 1
    for fam, cl in seen.items():
        if not cl:
            report.add(FAIL, sec, f"family {fam} has no cluster in the table")
    if nrows:
        report.add(OK, sec, f"{nrows} cluster->family rows; " + ", ".join(f"{f}: {len(seen[f])}" for f in sorted(seen)))
    report.add(OK, sec, "cluster existence is checked at stage 09 against the Step-8 output (needs a run)")


def check_bed(path: str, label: str, codes: set, report: Report, section: str) -> None:
    if not _exists(report, section, label, path):
        return
    n = 0
    with open(path, encoding="utf-8") as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                report.add(FAIL, section, f"{label} line {ln}: fewer than 3 columns")
                continue
            try:
                s, e = int(p[1]), int(p[2])
            except ValueError:
                report.add(FAIL, section, f"{label} line {ln}: start/end not integers")
                continue
            if not 0 <= s < e:
                report.add(FAIL, section, f"{label} line {ln}: need 0 <= start < end (BED is 0-based half-open)")
            if len(p) >= 4 and p[3].split("|")[0] not in codes:
                report.add(WARN, section, f"{label} line {ln}: name {p[3]!r} does not start with a known species code")
            n += 1
    report.add(OK, section, f"{label}: {n} intervals")


# ----------------------------------------------------------------------------- main


def check_files(cfg: Config, report: Report) -> None:
    d = cfg.data
    codes = set(d["genomes"])
    for code, g in sorted(d["genomes"].items()):
        if _exists(report, "genomes", code, g):
            fai = g + ".fai"
            if os.path.isfile(fai):
                report.add(OK, "genomes", f"{code}: index present ({fai})")
            else:
                report.add(WARN, "genomes", f"{code}: no .fai index; stage 11 will run `samtools faidx` "
                                             f"(needs write access next to the genome)")
    inputs = d.get("inputs") or {}
    for code, p in sorted((inputs.get("trf_dat") or {}).items()):
        _exists(report, "inputs.trf_dat", code, p)
    for fam, p in sorted((inputs.get("seeds") or {}).items()):
        check_bed(p, f"seeds {fam}", codes, report, "inputs.seeds")
    for fam, fd in sorted(d["families"].items()):
        if _exists(report, "families", f"{fam} consensus", fd["consensus"]):
            recs = read_fasta(fd["consensus"], 5)
            if len(recs) != 1:
                report.add(FAIL, "families", f"{fam}: consensus must hold exactly one record, found {len(recs)}")
            else:
                name, ln = recs[0]
                if ln != fd["mono"]:
                    report.add(WARN, "families", f"{fam}: consensus length {ln} != mono {fd['mono']} "
                                                  f"(mono drives every window size; check which is intended)")
                else:
                    report.add(OK, "families", f"{fam}: consensus {name}, {ln} bp = mono")
    rot = (d.get("references") or {}).get("rotated_consensus")
    if rot:
        if _exists(report, "references", "rotated_consensus", rot):
            names = {n for n, _ in read_fasta(rot, 1000)}
            missing = [f for f in d["families"] if f not in names]
            if missing:
                report.add(WARN, "references",
                           f"rotated_consensus lacks records named like families: {', '.join(missing)} "
                           f"(Step 9 identifies references by header; the Darevskia file uses CLsat1..4)")
    else:
        report.add(WARN, "references", "rotated_consensus not given: the driver will derive it by Step-5 rotation "
                                        "of each consensus and record the deviation (the published run used CLsat_cons.rot)")
    check_family_assignment(d["family_assignment"], set(d["families"]), report)
    sub = d.get("subfamily") or {}
    if sub.get("consensuses"):
        p = sub["consensuses"]
        names: set = set()
        if os.path.isdir(p):
            fas = [f for f in os.listdir(p) if re.search(r"\.(fa|fas|fasta|fna)$", f)]
            if fas:
                report.add(OK, "subfamily", f"consensuses: directory with {len(fas)} FASTA files")
                names = {re.sub(r"\.(fa|fas|fasta|fna)$", "", f) for f in fas}
            else:
                report.add(FAIL, "subfamily", f"consensuses: directory has no FASTA files: {p}")
        elif os.path.isfile(p):
            names = {n for n, _ in read_fasta(p)}
            report.add(OK, "subfamily", f"consensuses: FASTA with {len(names)} records")
        else:
            report.add(WARN, "subfamily", f"consensuses: not found: {p} (the sub-subfamily stage will write an empty table)")
        if names:
            for fam, members in (sub.get("members") or {}).items():
                lost = [m for m in members if m not in names]
                if lost:
                    report.add(WARN, "subfamily", f"members of {fam} without a consensus record/file: {', '.join(lost)}")
    else:
        report.add(OK, "subfamily", "no sub-subfamily consensuses: stage skipped")
    if d.get("manual_loci"):
        check_bed(d["manual_loci"], "manual_loci", codes, report, "manual_loci")
    out = d["out_dir"]
    if os.path.isdir(out):
        extra = " (non-empty; stages with DONE markers would be skipped)" if os.listdir(out) else ""
        report.add(OK, "out_dir", f"exists: {out}{extra}")
    elif os.path.isdir(os.path.dirname(out.rstrip("/")) or "."):
        report.add(OK, "out_dir", f"will be created: {out}")
    else:
        report.add(WARN, "out_dir", f"parent directory does not exist: {out}")
    # FASTA36 (ssearch36) cannot open a file whose path is 120+ characters. The expansion/scan script passes long
    # paths through short symlinks; the other stages put files at most ~45 characters below out_dir or TMPDIR.
    for what, p in (("out_dir", os.path.abspath(out)), ("TMPDIR", os.path.abspath(os.environ.get("TMPDIR") or "/tmp"))):
        if len(p) > MAX_BASE_PATH:
            report.add(FAIL, "paths", f"{what} is {len(p)} characters ({p}); keep it at most {MAX_BASE_PATH}, "
                                      "or ssearch36 cannot open the files below it (path limit 119)")


MAX_BASE_PATH = 70


def check_tools(cfg: Optional[Config], report: Report) -> None:
    overrides = dict(cfg.data.get("tools") or {}) if cfg else {}
    bin_dir = overrides.get("bin_dir")
    if bin_dir:
        if os.path.isdir(bin_dir):
            report.add(OK, "tools", f"bin_dir: {bin_dir}")
        else:
            report.add(FAIL, "tools", f"bin_dir not found: {bin_dir}")
    for pr in probe_all(overrides):
        t = pr.tool
        if not pr.found:
            report.add(FAIL if t.required else WARN, "tools",
                       f"{t.name}: not found (stages {t.stages})" + (f"; {t.note}" if t.note else ""))
        elif pr.error:
            report.add(FAIL if t.required else WARN, "tools", f"{t.name}: {pr.path}: {pr.error}")
        else:
            report.add(OK, "tools", f"{t.name} {pr.version}: {pr.path}")


def run_check(config_path: str, do_files: bool = True, do_tools: bool = True) -> Report:
    report = Report()
    try:
        cfg = load_config(config_path)
    except ConfigError as e:
        for err in e.errors:
            report.add(FAIL, "config", err)
        if do_tools:
            check_tools(None, report)
        return report
    except OSError as e:
        report.add(FAIL, "config", str(e))
        return report
    report.config = cfg
    report.add(OK, "config", f"{cfg.path}: schema valid; run_name={cfg.data['run_name']}, "
                             f"{len(cfg.genomes)} genomes, {len(cfg.families)} families")
    for w in cfg.warnings:
        report.add(WARN, "config", w)
    miss = missing_scripts()
    if miss:
        report.add(FAIL, "package", "stage scripts missing from the installed package: " + ", ".join(miss))
    else:
        report.add(OK, "package", f"{sum(len(s.scripts) for s in STAGES)} stage script references present")
    if do_files:
        check_files(cfg, report)
    if do_tools:
        check_tools(cfg, report)
    return report
