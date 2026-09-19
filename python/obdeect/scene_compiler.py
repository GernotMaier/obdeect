"""Compile a provenance-checked simulation-models IR into generic scene data.

This adapter deliberately does not invent dish sag, panel normals, camera
geometry, or material behaviour from incomplete catalogue records. It turns
the mirror-list asset into unit-normalised facet records and reports every
other imported parameter as deferred. The result is a stable handoff to the
native scene compiler, not a claim of a trace-ready production scene.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


class SceneCompileError(ValueError):
    """The provenance IR cannot be compiled without guessing optical data."""


_SHAPES = {0: "circle", 1: "hexagon_flat_y", 2: "square", 3: "hexagon_flat_x"}
_UNIT_TO_M = {"m": 1.0, "cm": 0.01, "mm": 0.001}
_CONSUMED = {"mirror_list", "focal_length", "mirror_focal_length"}


def _number(value: object, context: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise SceneCompileError(f"{context} must be a finite number")
    return float(value)


def _length_m(parameter: dict[str, Any], name: str, *, allow_zero: bool = False) -> float:
    unit = parameter.get("unit")
    if unit not in _UNIT_TO_M:
        raise SceneCompileError(f"{name} has unsupported length unit {unit!r}")
    value = _number(parameter.get("value"), name) * _UNIT_TO_M[unit]
    if value < 0.0 or (value == 0.0 and not allow_zero):
        qualifier = "non-negative" if allow_zero else "positive"
        raise SceneCompileError(f"{name} must be {qualifier}")
    return value


def parse_simtel_mirror_list(
    contents: str, *, fallback_focal_length_m: float | None
) -> list[dict[str, Any]]:
    """Parse the documented sim_telarray mirror-list columns, in centimetres.

    The sim_telarray reader consumes five required columns (x, y, diameter,
    focal length, shape) and an optional z position.  It ignores text after
    those fields, which model files commonly use for comments and panel IDs.
    This parser retains only the documented geometric values: trailing data is
    deliberately not interpreted as a normal, rotation, or alignment.

    Focal length zero is legal only when an explicit catalogue fallback is
    supplied.  It is not a prescription for a facet normal.
    """
    facets: list[dict[str, Any]] = []
    for line_number, source_line in enumerate(contents.splitlines(), start=1):
        line = source_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 5:
            raise SceneCompileError(f"mirror list line {line_number}: expected at least 5 columns")
        try:
            x_cm, y_cm, diameter_cm, focal_cm = (float(item) for item in fields[:4])
            shape_code = int(fields[4])
            # sim_telarray parses a sixth floating-point field when present.
            # A comment immediately after the required fields means no z was
            # supplied, rather than an invalid optical datum.
            z_cm = 0.0
            if len(fields) >= 6:
                try:
                    z_cm = float(fields[5])
                except ValueError:
                    pass
        except ValueError as error:
            raise SceneCompileError(
                f"mirror list line {line_number}: invalid numeric field"
            ) from error
        if not all(math.isfinite(item) for item in (x_cm, y_cm, diameter_cm, focal_cm, z_cm)):
            raise SceneCompileError(f"mirror list line {line_number}: non-finite value")
        if diameter_cm <= 0.0:
            raise SceneCompileError(f"mirror list line {line_number}: diameter must be positive")
        if shape_code not in _SHAPES:
            raise SceneCompileError(
                f"mirror list line {line_number}: unsupported shape {shape_code}"
            )
        focal_m = focal_cm * 0.01
        if focal_m == 0.0:
            if fallback_focal_length_m is None:
                raise SceneCompileError(
                    f"mirror list line {line_number}: zero focal length has no fallback"
                )
            focal_m = fallback_focal_length_m
        if focal_m <= 0.0:
            raise SceneCompileError(
                f"mirror list line {line_number}: focal length must be positive"
            )
        facet: dict[str, Any] = {
            "id": len(facets),
            "centre_m": [x_cm * 0.01, y_cm * 0.01, z_cm * 0.01],
            "shape": _SHAPES[shape_code],
            "diameter_m": diameter_cm * 0.01,
            "focal_length_m": focal_m,
        }
        facets.append(facet)
    if not facets:
        raise SceneCompileError("mirror list contains no facets")
    return facets


def compile_scene(ir: dict[str, Any], source_root: Path) -> dict[str, Any]:
    """Compile ``obdeect.simulation-models-ir.v1`` into generic scene data."""
    if ir.get("format") != "obdeect.simulation-models-ir.v1":
        raise SceneCompileError("expected obdeect.simulation-models-ir.v1")
    model = ir.get("model")
    version = ir.get("model_version")
    parameters = ir.get("parameters")
    assets = ir.get("assets")
    if not isinstance(model, str) or not model or not isinstance(version, str) or not version:
        raise SceneCompileError("IR model identity is invalid")
    if not isinstance(parameters, dict) or not isinstance(assets, dict):
        raise SceneCompileError("IR parameters and assets must be objects")
    mirror = parameters.get("mirror_list")
    asset = assets.get("mirror_list")
    if (
        not isinstance(mirror, dict)
        or mirror.get("file") is not True
        or not isinstance(asset, dict)
    ):
        raise SceneCompileError("IR does not provide a tracked mirror_list asset")
    relative_asset = asset.get("path")
    expected_hash = asset.get("sha256")
    if not isinstance(relative_asset, str) or not isinstance(expected_hash, str):
        raise SceneCompileError("mirror-list asset record is invalid")
    root = source_root.resolve()
    # Keep the compiler's source-root contract aligned with the importer: a
    # released checkout may contain its data package in simulation-models/.
    nested_root = root / "simulation-models"
    if not (root / "model_parameters").is_dir() and (nested_root / "model_parameters").is_dir():
        root = nested_root
    path = (root / relative_asset).resolve()
    if not path.is_relative_to(root) or not path.is_file():
        raise SceneCompileError("mirror-list asset is missing or escapes source root")
    actual_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual_hash != expected_hash:
        raise SceneCompileError("mirror-list asset hash does not match IR provenance")
    fallback = None
    if isinstance(parameters.get("mirror_focal_length"), dict):
        fallback = _length_m(
            parameters["mirror_focal_length"], "mirror_focal_length", allow_zero=True
        )
        # A zero catalogue fallback means that every panel must provide its
        # own focal length.  It is not an optical value to manufacture.
        if fallback == 0.0:
            fallback = None
    facets = parse_simtel_mirror_list(
        path.read_text(encoding="utf-8"), fallback_focal_length_m=fallback
    )
    explicit_z_count = sum(facet["centre_m"][2] != 0.0 for facet in facets)
    report = {
        "consumed": sorted(name for name in parameters if name in _CONSUMED),
        "deferred": sorted(name for name in parameters if name not in _CONSUMED),
        "unsupported": [],
        "trace_blockers": [
            "facet normals and tip/tilt alignment are unavailable: sim_telarray mirror-list "
            "columns encode centres, footprint, focal length, shape, and optional z only",
            "camera, obstructions, and material response remain deferred",
        ],
        "facet_geometry_evidence": {
            "source_format": "sim_telarray mirror-list (x, y, diameter, focal_length, shape[, z])",
            "facet_count": len(facets),
            "explicit_nonzero_z_count": explicit_z_count,
            "normal_status": "unavailable",
            "in_plane_orientation_status": "unavailable",
            "alignment_status": "unavailable",
            "interpretation": "No normals, rotations, or alignment "
            "are inferred from facet centres, "
            "focal lengths, shape codes, or z positions.",
        },
    }
    compiled = {
        "format": "obdeect.compiled-scene.v1",
        "provenance": {
            "model": model,
            "model_version": version,
            "input_records": ir.get("input_records", {}),
            "assets": assets,
        },
        "primary": {"kind": "segmented_mirror", "facets": facets},
        "report": report,
    }
    canonical = json.dumps(compiled, sort_keys=True, separators=(",", ":")).encode()
    compiled["scene_sha256"] = hashlib.sha256(canonical).hexdigest()
    return compiled


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compile a simulation-models IR into generic scene data."
    )
    parser.add_argument("ir", type=Path, help="provenance-checked scene IR JSON")
    parser.add_argument(
        "--source-root", type=Path, required=True, help="root used to resolve IR asset paths"
    )
    parser.add_argument("--output", type=Path, required=True, help="compiled generic-scene JSON")
    args = parser.parse_args()
    try:
        ir = json.loads(args.ir.read_text(encoding="utf-8"))
        if not isinstance(ir, dict):
            raise SceneCompileError("IR root must be an object")
        scene = compile_scene(ir, args.source_root)
    except (OSError, json.JSONDecodeError, SceneCompileError) as error:
        raise SystemExit(f"scene compilation failed: {error}") from error
    args.output.write_text(json.dumps(scene, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
