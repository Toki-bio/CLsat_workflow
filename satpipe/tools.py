"""External tool registry and version probing."""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class Tool:
    name: str
    version_args: tuple
    version_re: str
    required: bool
    stages: str
    note: str = ""
    min_version: Optional[tuple] = None


# Versions of the published run (PIPELINE.md): TRF 4.10.0-rc.2, VSEARCH 2.30.0, MAFFT 7.526,
# FASTA 36.3.8e, EMBOSS 6.6.0, BLAST+ 2.17.0, samtools 1.22.1, bedtools 2.31.1, Python 3.12.11.
TOOLS = [
    Tool("bash", ("--version",), r"version (\d+\.\d+[\w.()-]*)", True, "all",
         "bash >= 4 (associative arrays)", (4,)),
    Tool("gawk", ("--version",), r"GNU Awk (\d+\.\d+\.\d+)", True, "04, 10, 14",
         "must be GNU awk: PROCINFO[\"sorted_in\"], asorti, 3-argument match; mawk/busybox awk break Steps 4 and 10",
         (4,)),
    Tool("python3", ("--version",), r"Python (\d+\.\d+\.\d+)", True, "driver, 14, orthology, subfamily",
         "python >= 3.10", (3, 10)),
    Tool("trf", ("-v",), r"Version (\d+\.\d+\.\d+[\w.-]*)", True, "02, 04, finalisation",
         "exit status is not 0 on success in every build; scripts use `|| true`"),
    Tool("vsearch", ("--version",), r"vsearch v?(\d+\.\d+\.\d+)", True, "06, 08"),
    Tool("mafft", ("--version",), r"v(\d+\.\d+)", True, "07"),
    Tool("ssearch36", ("-h",), r"version[: ]+(\d+\.\d+\.\d+\w*)", True, "11, 13, orthology, subfamily, terminal",
         "FASTA package (fasta3 / fasta36)"),
    Tool("samtools", ("--version",), r"samtools (\d+\.\d+[\w.]*)", True, "11+", "", (1, 12)),
    Tool("bedtools", ("--version",), r"bedtools v?(\d+\.\d+\.\d+)", True, "11, 12, 13, 14, 15, 20"),
    Tool("parallel", ("--version",), r"GNU parallel (\d+)", True, "11, 13, finalisation", "GNU parallel"),
    Tool("cons", ("-version",), r"EMBOSS:(\d+\.\d+\.\d+[\w.]*)", False, "finalisation (consensus deliverable)", "EMBOSS"),
    Tool("water", ("-version",), r"EMBOSS:(\d+\.\d+\.\d+[\w.]*)", False, "finalisation", "EMBOSS"),
    Tool("blastn", ("-version",), r"blastn: (\d+\.\d+\.\d+\+?)", False, "07 (IMO-TaqI comparison)", "BLAST+"),
    Tool("makeblastdb", ("-version",), r"makeblastdb: (\d+\.\d+\.\d+\+?)", False, "07", "BLAST+"),
]


@dataclass
class Probe:
    tool: Tool
    path: Optional[str]
    version: Optional[str]
    error: Optional[str]

    @property
    def found(self) -> bool:
        return self.path is not None


def _parse_version(s: str) -> tuple:
    return tuple(int(x) for x in re.findall(r"\d+", s)[:3])


def _exe(p: str) -> bool:
    return os.path.isfile(p) and os.access(p, os.X_OK)


def locate(name: str, override: Optional[str] = None, bin_dir: Optional[str] = None) -> Optional[str]:
    if override:
        return override.replace("\\", "/") if _exe(override) else None
    if bin_dir:
        cand = os.path.join(bin_dir, name)
        if _exe(cand):
            return cand.replace("\\", "/")
    found = shutil.which(name)
    if found is None and name == "python3":
        found = shutil.which("python")
    return found.replace("\\", "/") if found else None


def probe(tool: Tool, override: Optional[str] = None, bin_dir: Optional[str] = None, timeout: int = 30) -> Probe:
    path = locate(tool.name, override, bin_dir)
    if path is None:
        return Probe(tool, None, None, "not found" + (f" (override {override})" if override else ""))
    try:
        cp = subprocess.run([path, *tool.version_args], capture_output=True, text=True, timeout=timeout,
                            stdin=subprocess.DEVNULL, errors="replace")
        out = (cp.stdout or "") + "\n" + (cp.stderr or "")
    except (OSError, subprocess.TimeoutExpired) as e:
        return Probe(tool, path, None, f"could not run: {e}")
    m = re.search(tool.version_re, out)
    version = m.group(1) if m else None
    err = None
    if version is None:
        first = (out.strip().splitlines() or ["(no output)"])[0][:80]
        err = "version not recognised: " + first
    elif tool.min_version and _parse_version(version) < tool.min_version:
        err = f"version {version} < required {'.'.join(map(str, tool.min_version))}"
    if tool.name == "gawk" and version and "GNU Awk" not in out:
        err = "not GNU awk"
    return Probe(tool, path, version, err)


def probe_all(overrides: Optional[dict] = None) -> list:
    overrides = overrides or {}
    bin_dir = overrides.get("bin_dir")
    return [probe(t, overrides.get(t.name), bin_dir) for t in TOOLS]
