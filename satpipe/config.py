"""Load and validate a satpipe config file."""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any

import yaml

from .schema import SCHEMA, validate


class ConfigError(Exception):
    def __init__(self, errors: list):
        super().__init__("\n".join(errors))
        self.errors = errors


@dataclass
class Config:
    path: str
    base_dir: str
    data: dict
    warnings: list = field(default_factory=list)

    @property
    def genomes(self) -> dict:
        return self.data["genomes"]

    @property
    def families(self) -> dict:
        return self.data["families"]

    def tool_override(self, name: str):
        return (self.data.get("tools") or {}).get(name)

    @property
    def bin_dir(self):
        return (self.data.get("tools") or {}).get("bin_dir")


def _cross_checks(d: dict, errors: list, warnings: list) -> None:
    """Checks that need more than one key."""
    codes = set(d["genomes"] or {})
    fams = set(d["families"] or {})
    inputs = d.get("inputs") or {}
    trf = inputs.get("trf_dat") or {}
    for code in trf:
        if code not in codes:
            errors.append(f"inputs.trf_dat.{code}: not a code in `genomes`")
    if trf and set(trf) <= codes and set(trf) != codes:
        warnings.append("inputs.trf_dat covers only part of `genomes`; stage 02 (TRF) will run for the rest: "
                        + ", ".join(sorted(codes - set(trf))))
    for fam in (inputs.get("seeds") or {}):
        if fam not in fams:
            errors.append(f"inputs.seeds.{fam}: not a family in `families`")
    sub = d.get("subfamily") or {}
    for fam, members in (sub.get("members") or {}).items():
        if fam not in fams:
            errors.append(f"subfamily.members.{fam}: not a family in `families`")
        if members and len(set(members)) != len(members):
            errors.append(f"subfamily.members.{fam}: duplicate names")
        for m in members or []:
            if m in fams and m != fam:
                errors.append(f"subfamily.members.{fam}: '{m}' is also a family name")
    if sub.get("members") and not sub.get("consensuses"):
        errors.append("subfamily.members given but subfamily.consensuses is missing")
    for code in ((d.get("viewer") or {}).get("species_names") or {}):
        if code not in codes:
            errors.append(f"viewer.species_names.{code}: not a code in `genomes`")
    per = (d.get("discovery") or {}).get("period")
    if isinstance(per, list) and len(per) == 2 and all(isinstance(x, int) for x in per) and per[0] > per[1]:
        errors.append(f"discovery.period: min {per[0]} > max {per[1]}")


def load_config(path: str) -> Config:
    """Parse YAML, validate against the schema, fill defaults, resolve paths.
    Raises ConfigError with every problem found (not just the first)."""
    path = os.path.abspath(path)
    base_dir = os.path.dirname(path)
    with open(path, encoding="utf-8") as fh:
        try:
            raw = yaml.safe_load(fh)
        except yaml.YAMLError as e:
            raise ConfigError([f"YAML parse error: {e}"])
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise ConfigError(["top level must be a mapping"])
    errors: list = []
    warnings: list = []
    data = validate(raw, SCHEMA, "", errors, base_dir)
    if not errors:
        _cross_checks(data, errors, warnings)
    if errors:
        raise ConfigError(errors)
    return Config(path=path.replace("\\", "/"), base_dir=base_dir.replace("\\", "/"),
                  data=data, warnings=warnings)


def dump_config(cfg: Config) -> str:
    return yaml.safe_dump(cfg.data, sort_keys=False, default_flow_style=False, allow_unicode=True)
