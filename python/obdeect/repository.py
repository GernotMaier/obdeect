"""Read a telescope's parameter records from a simulation-models checkout."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from obdeect.formats import ModelError
from obdeect.propagation import ambient_group_index


@dataclass(frozen=True)
class TelescopeRecords:
    """Parameter records of one telescope model and the files they name."""

    model: str
    version: str
    site: str | None
    records: dict[str, dict[str, Any]]
    files: dict[str, Path]
    propagation_group_index: float

    def read_file(self, name: str) -> str:
        """Text of the file named by parameter `name`."""
        if name not in self.files:
            raise ModelError(f"parameter {name} has no file in the simulation-models checkout")
        return self.files[name].read_text(encoding="utf-8")

    @property
    def provenance(self) -> dict[str, Any]:
        return {"model": self.model, "model_version": self.version, "site": self.site}


def _component(text: object) -> str:
    if not isinstance(text, str) or text in ("", ".", "..") or any(c in text for c in "/\\\0"):
        raise ModelError(f"invalid name in simulation-models: {text!r}")
    return text


def _load(path: Path) -> dict[str, Any]:
    try:
        content = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ModelError(f"cannot read {path}: {error}") from error
    if not isinstance(content, dict):
        raise ModelError(f"{path}: expected a JSON object")
    return content


def _parameter_path(root: Path, instrument: str, name: str, version: str) -> Path:
    return root / "model_parameters" / instrument / name / f"{name}-{version}.json"


def _production_versions(root: Path, version: str) -> list[str]:
    """The production and, for a patch production, the productions it patches
    (its model_version_history), oldest first, as simtools merges them."""
    versions = {version}
    info = root / "productions" / version / "info.yml"
    text = info.read_text(encoding="utf-8") if info.is_file() else ""
    if re.search(r"^model_update:\s*[\"']?patch_update", text, re.MULTILINE):
        history = re.search(
            r"^model_version_history:(.*?)(?=^\S|\Z)", text, re.MULTILINE | re.DOTALL
        )
        versions |= set(re.findall(r"\d+\.\d+\.\d+", history[1] if history else ""))
    return sorted(versions, key=lambda v: tuple(int(part) for part in v.split(".")))


def _production_table(
    root: Path, versions: list[str], instrument: str
) -> tuple[dict[str, str], str | None]:
    """Parameter versions and design model of `instrument`, merged over `versions`."""
    table: dict[str, str] = {}
    parent = None
    deprecated: list[str] = []
    found = False
    for version in versions:
        path = root / "productions" / version / f"{instrument}.json"
        if not path.is_file():
            continue
        manifest = _load(path)
        parameters = manifest.get("parameters", {}).get(instrument)
        if manifest.get("model_version") != version or not isinstance(parameters, dict):
            raise ModelError(f"invalid production table {instrument} {version}")
        found = True
        table |= parameters
        parent = manifest.get("design_model", {}).get(instrument, parent)
        deprecated = manifest.get("deprecated_parameters", deprecated)
    if not found:
        raise ModelError(f"no production table {instrument} in {versions[-1]}")
    return {k: v for k, v in table.items() if k not in deprecated}, parent


def read_telescope(root: Path, model: str, version: str) -> TelescopeRecords:
    """Records of `model` in production `version` (a full or a patch
    production), following design-model inheritance, plus the site's air
    group index at the observation level."""
    root = Path(root).resolve()
    if (
        not (root / "productions").is_dir()
        and (root / "simulation-models" / "productions").is_dir()
    ):
        root = root / "simulation-models"
    versions = _production_versions(root, _component(version))
    chosen: dict[str, tuple[str, str]] = {}

    def inherit(instrument: str, chain: tuple[str, ...] = ()) -> None:
        if _component(instrument) in chain:
            raise ModelError("cyclic design-model inheritance")
        table, parent = _production_table(root, versions, instrument)
        if parent is not None:
            inherit(parent, (*chain, instrument))
        for name, parameter_version in table.items():
            chosen[_component(name)] = (instrument, _component(parameter_version))

    inherit(model)
    records: dict[str, dict[str, Any]] = {}
    files: dict[str, Path] = {}
    sites = set()
    for name, (instrument, parameter_version) in sorted(chosen.items()):
        path = _parameter_path(root, instrument, name, parameter_version)
        record = _load(path)
        if record.get("parameter") != name or record.get("parameter_version") != parameter_version:
            raise ModelError(f"{path}: identity does not match the production table")
        records[name] = record
        if isinstance(record.get("site"), str):
            sites.add(record["site"])
        value = record.get("value")
        if record.get("file") and isinstance(value, str) and value.lower() != "none":
            filename = _component(value.partition("#rpol:")[0])
            for candidate in (
                path.parent / filename,
                root / "model_parameters" / "Files" / filename,
            ):
                if candidate.is_file():
                    files[name] = candidate
                    break
    if len(sites) > 1:
        raise ModelError(f"{model}: parameters disagree on the site")
    site = next(iter(sites), None)
    return TelescopeRecords(
        model, version, site, records, files, _group_index(root, versions, site)
    )


def _group_index(root: Path, versions: list[str], site: str | None) -> float:
    """Air group index at the observation level, from the site's atmosphere."""
    if site is None:
        return 1.0
    instrument = f"OBS-{_component(site)}"
    table = _production_table(root, versions, instrument)[0]
    profile = _load(
        _parameter_path(root, instrument, "atmospheric_profile", table["atmospheric_profile"])
    )
    level = _load(
        _parameter_path(
            root, instrument, "corsika_observation_level", table["corsika_observation_level"]
        )
    )
    if level.get("unit") != "m":
        raise ModelError("corsika_observation_level must be in metres")
    path = (
        root
        / "model_parameters"
        / instrument
        / "atmospheric_profile"
        / _component(profile["value"])
    )
    return ambient_group_index(path.read_text(encoding="utf-8"), float(level["value"]))
