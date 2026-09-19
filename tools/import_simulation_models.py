#!/usr/bin/env python3
"""Resolve one simulation-models production into a provenance-bearing scene IR.

This importer intentionally has no dependency beyond the Python standard
library. It does not interpret sim_telarray's geometry files yet; it freezes
the complete parameter selection and all declared file assets so that a later
C++ geometry compiler has an auditable, deterministic input.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


class ImportError(ValueError):
    """An input is missing or violates the simulation-models manifest contract."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
    """Return a complete, JSON-serializable model IR or raise ImportError."""
    root = root.resolve()
    if Path(model).name != model or Path(version).name != version:
        raise ImportError("model and version must be simple path components")
    manifest_path = root / "productions" / version / f"{model}.json"
    manifest = load_json(manifest_path)
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
        if not isinstance(name, str) or not isinstance(parameter_version, str):
            raise ImportError("parameter names and versions must be strings")
        parameter_path = root / "model_parameters" / model / name / f"{name}-{parameter_version}.json"
        parameter = load_json(parameter_path)
        if parameter.get("instrument") != model or parameter.get("parameter") != name:
            raise ImportError(f"parameter identity mismatch in {parameter_path}")
        if parameter.get("parameter_version") != parameter_version:
            raise ImportError(f"parameter version mismatch in {parameter_path}")
        if "value" not in parameter or "unit" not in parameter or "file" not in parameter:
            raise ImportError(f"parameter record is incomplete: {parameter_path}")
        parameters[name] = {key: parameter[key] for key in ("type", "unit", "value", "file")}
        records[f"parameter:{name}"] = record(parameter_path, root)
        if parameter["file"]:
            value = parameter["value"]
            # A file-capable parameter may intentionally be disabled by a
            # null value. Preserve that selection in ``parameters`` without
            # manufacturing an asset dependency.
            if value is None:
                continue
            if not isinstance(value, str) or Path(value).name != value:
                raise ImportError(f"unsafe or invalid asset name in {parameter_path}")
            asset_path = root / "model_parameters" / "Files" / value
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


if __name__ == "__main__":
    main()
