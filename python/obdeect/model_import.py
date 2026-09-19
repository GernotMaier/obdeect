"""Import a selected simulation-models production into auditable scene IR.

The importer deliberately uses the standard library only.  It records the
complete selected parameter records and hashes declared assets; compiling
those records into optical geometry remains an explicit later stage.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


class ImportError(ValueError):
    """An input violates the simulation-models manifest contract."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def component(value: str) -> bool:
    return isinstance(value, str) and value not in ("", ".", "..") and "/" not in value and "\\" not in value


def within_root(path: Path, root: Path) -> Path:
    if not path.resolve().is_relative_to(root):
        raise ImportError(f"source path escapes checkout: {path}")
    return path


def load_json(path: Path) -> dict[str, Any]:
    try:
        content = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ImportError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(content, dict):
        raise ImportError(f"JSON root must be an object: {path}")
    return content


def record(path: Path, root: Path) -> dict[str, str]:
    return {"path": path.relative_to(root).as_posix(), "sha256": sha256(path)}


def resolve_model(root: Path, model: str, version: str) -> dict[str, Any]:
    """Return complete JSON-serializable provenance IR or raise ImportError."""
    root = root.resolve()
    # Released simulation-models checkouts commonly contain the data package
    # in a nested ``simulation-models/`` directory, while exported data trees
    # use that directory itself as their root.  Accept both layouts without
    # searching arbitrary descendants.
    nested_root = root / "simulation-models"
    if not (root / "productions").is_dir() and (nested_root / "productions").is_dir():
        root = nested_root
    if not component(model) or not component(version):
        raise ImportError("model and version must be simple path components")
    manifest_path = root / "productions" / version / f"{model}.json"
    manifest = load_json(within_root(manifest_path, root))
    if manifest.get("model_version") != version:
        raise ImportError(f"manifest version mismatch in {manifest_path}")
    if manifest.get("production_table_name") != model:
        raise ImportError(f"manifest table name mismatch in {manifest_path}")
    tables = manifest.get("parameters")
    if not isinstance(tables, dict) or set(tables) != {model} or not isinstance(tables[model], dict):
        raise ImportError("production manifest must contain exactly the requested parameter table")

    parameters: dict[str, Any] = {}
    records: dict[str, dict[str, str]] = {"production_manifest": record(manifest_path, root)}
    assets: dict[str, dict[str, str]] = {}
    for name, parameter_version in sorted(tables[model].items()):
        if not component(name) or not component(parameter_version):
            raise ImportError("parameter names and versions must be strings")
        parameter_path = root / "model_parameters" / model / name / f"{name}-{parameter_version}.json"
        parameter = load_json(within_root(parameter_path, root))
        if parameter.get("instrument") != model or parameter.get("parameter") != name:
            raise ImportError(f"parameter identity mismatch in {parameter_path}")
        if parameter.get("parameter_version") != parameter_version:
            raise ImportError(f"parameter version mismatch in {parameter_path}")
        if not all(key in parameter for key in ("type", "value", "unit", "file")):
            raise ImportError(f"parameter record is incomplete: {parameter_path}")
        if not isinstance(parameter["file"], bool):
            raise ImportError(f"file flag must be boolean: {parameter_path}")
        parameters[name] = parameter
        records[f"parameter:{name}"] = record(parameter_path, root)
        if parameter["file"] and parameter["value"] is not None:
            value = parameter["value"]
            if not component(value):
                raise ImportError(f"unsafe or invalid asset name in {parameter_path}")
            asset_path = within_root(root / "model_parameters" / "Files" / value, root)
            if not asset_path.is_file():
                raise ImportError(f"declared model asset is missing: {asset_path}")
            assets[name] = record(asset_path, root)
    return {
        "format": "obdeect.simulation-models-ir.v1",
        "model": model,
        "model_version": version,
        "source_root": root.name,
        "input_records": records,
        "assets": assets,
        "parameters": parameters,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Import a pinned simulation-models production manifest.")
    parser.add_argument("root", type=Path, help="path to the simulation-models repository root")
    parser.add_argument("model", help="production table, e.g. LSTN-design")
    parser.add_argument("--version", default="6.3.0", help="production model version")
    parser.add_argument("--output", type=Path, required=True, help="destination scene-IR JSON")
    args = parser.parse_args()
    try:
        scene = resolve_model(args.root, args.model, args.version)
    except ImportError as error:
        raise SystemExit(f"import failed: {error}") from error
    args.output.write_text(json.dumps(scene, indent=2, sort_keys=True) + "\n", encoding="utf-8")
