"""Compile a provenance-checked simulation-models IR into generic scene data.

This adapter verifies the selected production and extracts documented mirror
footprints and camera pixel layouts. It does not invent dish sag, panel
normals, detector surfaces, or material behaviour. The output is an audited
handoff, not a trace-ready production scene.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from obdeect.camera_config import CameraConfigError, parse_camera_layout
from obdeect.model_import import ImportError as ModelImportError
from obdeect.model_import import component, record, resolve_model


class SceneCompileError(ValueError):
    """The provenance IR cannot be compiled without guessing optical data."""


_SHAPES = {0: "circle", 1: "hexagon_flat_y", 2: "square", 3: "hexagon_flat_x"}
_UNIT_TO_M = {"m": 1.0, "cm": 0.01, "mm": 0.001}


def parse_simtel_segmentation(contents: str) -> list[dict[str, Any]]:
    """Parse explicit sim_telarray hex and ring footprint records, in cm/deg.

    Ring groups are expanded to stable per-segment IDs. The recorded gap is
    retained but no physical mask or surface normal is inferred from it.
    """
    segments: list[dict[str, Any]] = []
    for line_number, source_line in enumerate(contents.splitlines(), start=1):
        fields = source_line.split("#", 1)[0].split()
        if not fields:
            continue
        kind = fields[0].lower()
        if kind not in {"hex", "yhex", "ring"}:
            raise SceneCompileError(
                f"segmentation line {line_number}: unsupported type {fields[0]}"
            )
        expected_fields = {5, 6} if kind in {"hex", "yhex"} else {5, 6, 7}
        if len(fields) not in expected_fields:
            raise SceneCompileError(f"segmentation line {line_number}: wrong field count")
        try:
            count = int(fields[1])
            values = [float(value) for value in fields[2:]]
        except ValueError as error:
            raise SceneCompileError(f"segmentation line {line_number}: invalid number") from error
        if count < 1 or not all(math.isfinite(value) for value in values):
            raise SceneCompileError(f"segmentation line {line_number}: invalid count or value")
        if kind in {"hex", "yhex"}:
            x_cm, y_cm, diameter_cm, rotation_deg = (*values, 0.0)[0:4]
            if count != 1 or diameter_cm <= 0.0:
                raise SceneCompileError(f"segmentation line {line_number}: invalid hex footprint")
            segments.append({
                "id": len(segments),
                "shape": "hexagon",
                "centre_xy_m": [x_cm * 0.01, y_cm * 0.01],
                "diameter_m": diameter_cm * 0.01,
                "rotation_deg": rotation_deg,
            })
        else:
            inner_cm, outer_cm, span_deg, start_deg, gap_cm = (*values, 0.0, 0.0)[0:5]
            if inner_cm < 0.0 or outer_cm <= inner_cm or span_deg <= 0.0 or gap_cm < 0.0:
                raise SceneCompileError(f"segmentation line {line_number}: invalid ring footprint")
            if count * span_deg > 360.0 + 1e-9:
                raise SceneCompileError(f"segmentation line {line_number}: ring exceeds full turn")
            for index in range(count):
                segments.append({
                    "id": len(segments),
                    "shape": "annular_sector",
                    "inner_radius_m": inner_cm * 0.01,
                    "outer_radius_m": outer_cm * 0.01,
                    "start_deg": start_deg + index * 360.0 / count,
                    "span_deg": span_deg,
                    "gap_m": gap_cm * 0.01,
                })
    if not segments:
        raise SceneCompileError("segmentation contains no segments")
    return segments


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


def _signed_length_m(parameter: dict[str, Any], name: str) -> float:
    unit = parameter.get("unit")
    if unit not in _UNIT_TO_M:
        raise SceneCompileError(f"{name} has unsupported length unit {unit!r}")
    return _number(parameter.get("value"), name) * _UNIT_TO_M[unit]


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
    if fallback_focal_length_m is not None and (
        not math.isfinite(fallback_focal_length_m) or fallback_focal_length_m <= 0.0
    ):
        raise SceneCompileError("fallback focal length must be finite and positive")
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
            if len(fields) >= 6 and not fields[5].startswith("#"):
                z_cm = float(fields[5])
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


def derive_nominal_single_reflector(
    facets: list[dict[str, Any]], parameters: dict[str, Any]
) -> None:
    """Attach the unperturbed sim_telarray panel centres and normals.

    This follows ``tel_setup_primary`` in sim_imaging.c. Random distance,
    focal-length and alignment draws remain separate unresolved run inputs.
    """
    focal_length = _length_m(parameters["focal_length"], "focal_length")
    dish_length = _length_m(parameters["dish_shape_length"], "dish_shape_length", allow_zero=True)
    if dish_length == 0:
        dish_length = focal_length
    offset = _signed_length_m(parameters["mirror_offset"], "mirror_offset")
    parabolic = parameters["parabolic_dish"].get("value")
    if not isinstance(parabolic, bool):
        raise SceneCompileError("parabolic_dish must be boolean")
    for facet in facets:
        x, y, file_z = facet["centre_m"]
        radius = math.hypot(x, y)
        if file_z:
            height = abs(file_z)
        elif parabolic:
            height = radius * radius / (4 * dish_length)
        else:
            if radius > dish_length:
                raise SceneCompileError("panel centre lies outside Davies-Cotton dish radius")
            height = dish_length - math.sqrt(dish_length * dish_length - radius * radius)
        distance = math.hypot(focal_length - height, radius)
        if file_z > 0:
            z = file_z - offset
        else:
            z = focal_length - math.sqrt(distance * distance - radius * radius) - offset
        inclination = 0.5 * math.asin(radius / distance)
        if radius:
            nx = -math.sin(inclination) * x / radius
            ny = -math.sin(inclination) * y / radius
        else:
            nx = ny = 0.0
        facet["nominal_centre_m"] = [x, y, z]
        facet["nominal_normal"] = [nx, ny, math.cos(inclination)]


def compile_scene(
    ir: dict[str, Any], source_root: Path, *, simtel_root: Path | None = None
) -> dict[str, Any]:
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
    root = source_root.resolve()
    # Keep the compiler's source-root contract aligned with the importer: a
    # released checkout may contain its data package in simulation-models/.
    nested_root = root / "simulation-models"
    if not (root / "model_parameters").is_dir() and (nested_root / "model_parameters").is_dir():
        root = nested_root
    records = ir.get("input_records", {})
    if not isinstance(records, dict):
        raise SceneCompileError("IR input_records must be an object")
    try:
        source_ir = resolve_model(source_root, model, version)
    except ModelImportError as error:
        raise SceneCompileError(f"cannot verify source production: {error}") from error
    for key in ("source_root", "input_records", "assets", "parameters"):
        if ir.get(key) != source_ir[key]:
            raise SceneCompileError(f"IR {key} differs from the verified source production")
    verified_assets = {name: root / entry["path"] for name, entry in assets.items()}

    fallback = None
    if isinstance(parameters.get("mirror_focal_length"), dict):
        fallback = _length_m(
            parameters["mirror_focal_length"], "mirror_focal_length", allow_zero=True
        )
        # A zero catalogue fallback means that every panel must provide its
        # own focal length.  It is not an optical value to manufacture.
        if fallback == 0.0:
            fallback = None
    consumed = set()
    nominal_geometry = False
    if "mirror_list" in verified_assets:
        facets = parse_simtel_mirror_list(
            verified_assets["mirror_list"].read_text(encoding="utf-8"),
            fallback_focal_length_m=fallback,
        )
        primary = {"kind": "segmented_mirror", "facets": facets}
        consumed.add("mirror_list")
        if "mirror_focal_length" in parameters:
            consumed.add("mirror_focal_length")
        nominal_fields = {"focal_length", "dish_shape_length", "mirror_offset", "parabolic_dish"}
        nominal_geometry = parameters.get("mirror_class", {}).get(
            "value"
        ) == 0 and nominal_fields.issubset(parameters)
        if nominal_geometry:
            derive_nominal_single_reflector(facets, parameters)
            consumed.update(nominal_fields)
            consumed.add("mirror_class")
        evidence = {
            "source_format": "sim_telarray mirror-list (x, y, diameter, focal_length, shape[, z])",
            "facet_count": len(facets),
            "explicit_nonzero_z_count": sum(facet["centre_m"][2] != 0.0 for facet in facets),
            "normal_status": "nominal_unperturbed" if nominal_geometry else "unavailable",
            "in_plane_orientation_status": "unavailable",
            "alignment_status": "unavailable",
            "interpretation": (
                "Nominal panel centres and normals follow sim_telarray tel_setup_primary; "
                "run-specific random alignment and distance remain unresolved."
                if nominal_geometry
                else "No normals, rotations, or alignment are inferred from facet centres."
            ),
        }
    elif "primary_mirror_segmentation" in verified_assets:
        segments = parse_simtel_segmentation(
            verified_assets["primary_mirror_segmentation"].read_text(encoding="utf-8")
        )
        primary = {"kind": "segmented_footprints", "segments": segments}
        consumed.add("primary_mirror_segmentation")
        evidence = {
            "source_format": "sim_telarray segmentation (hex/ring)",
            "facet_count": len(segments),
            "normal_status": "unavailable",
            "alignment_status": "unavailable",
            "interpretation": "Footprints parsed; surface sag and normals unresolved.",
        }
    else:
        raise SceneCompileError("IR has no tracked primary mirror geometry asset")
    secondary = None
    if "secondary_mirror_segmentation" in verified_assets:
        secondary = {
            "kind": "segmented_footprints",
            "segments": parse_simtel_segmentation(
                verified_assets["secondary_mirror_segmentation"].read_text(encoding="utf-8")
            ),
        }
        consumed.add("secondary_mirror_segmentation")
    camera = None
    nested_assets = {}
    unresolved_references = []
    if "camera_config_file" in verified_assets:
        try:
            camera = parse_camera_layout(
                verified_assets["camera_config_file"].read_text(encoding="utf-8")
            )
        except CameraConfigError as error:
            raise SceneCompileError(str(error)) from error
        consumed.add("camera_config_file")
        declared_pixels = parameters.get("camera_pixels", {}).get("value")
        if isinstance(declared_pixels, bool) or not isinstance(declared_pixels, int):
            raise SceneCompileError("camera_pixels must be an integer count")
        if len(camera["pixels"]) != declared_pixels:
            raise SceneCompileError("camera pixel count differs from production record")
        consumed.add("camera_pixels")
        for pixel_type in camera["pixel_types"]:
            for filename in pixel_type["response_files"]:
                if not component(filename):
                    raise SceneCompileError(f"unsafe camera response filename: {filename}")
                path = root / "model_parameters" / "Files" / filename
                if not path.resolve().is_relative_to(root):
                    raise SceneCompileError(f"camera response path escapes source root: {filename}")
                if path.is_file():
                    nested_assets[filename] = record(path, root)
                else:
                    # sim_telarray fileopen() also searches its compiled-in
                    # cfg/CTA path. Require that root explicitly so the
                    # fallback cannot vary unnoticed between installations.
                    fallback = (
                        simtel_root.resolve() / "cfg" / "CTA" / filename
                        if simtel_root is not None
                        else None
                    )
                    if (
                        fallback is not None
                        and fallback.is_file()
                        and fallback.resolve().is_relative_to(simtel_root.resolve())
                    ):
                        nested_assets[filename] = {
                            **record(fallback, simtel_root.resolve()),
                            "source_root": "sim_telarray",
                        }
                    else:
                        unresolved_references.append(filename)
    deferred = sorted(set(parameters) - consumed)
    for name in deferred:
        if parameters[name].get("required_for_trace") is True:
            raise SceneCompileError(f"required field {name} is not supported by the compiler")
    report = {
        "native_trace_ready": False,
        "consumed": sorted(consumed),
        "deferred": deferred,
        "unsupported": [],
        "field_coverage": {
            name: {
                "disposition": "consumed" if name in consumed else "deferred",
                "reason": "parsed into geometry"
                if name in consumed
                else "not compiled into an optical scene",
            }
            for name in sorted(parameters)
        },
        "trace_blockers": [
            "run-specific panel alignment and distance are not compiled"
            if nominal_geometry
            else "facet surface normals and alignment are not compiled",
            "physical detector surfaces, obstructions, and materials remain deferred",
            *[f"camera response asset is missing: {name}" for name in unresolved_references],
        ],
        "facet_geometry_evidence": evidence,
    }
    if camera is not None:
        xs = [pixel["centre_xy_m"][0] for pixel in camera["pixels"]]
        ys = [pixel["centre_xy_m"][1] for pixel in camera["pixels"]]
        report["camera_layout_evidence"] = {
            "pixel_count": len(camera["pixels"]),
            "pixel_type_count": len(camera["pixel_types"]),
            "centre_bounds_xy_m": [[min(xs), min(ys)], [max(xs), max(ys)]],
            "unresolved_response_files": sorted(set(unresolved_references)),
            "surface_status": "unavailable",
        }
    compiled = {
        "format": "obdeect.compiled-scene.v1",
        "provenance": {
            "model": model,
            "model_version": version,
            "input_records": records,
            "assets": assets,
            "nested_assets": nested_assets,
        },
        "primary": primary,
        "report": report,
    }
    if secondary is not None:
        compiled["secondary"] = secondary
    if camera is not None:
        compiled["camera"] = camera
    canonical = json.dumps(compiled, sort_keys=True, separators=(",", ":")).encode()
    compiled["scene_sha256"] = hashlib.sha256(canonical).hexdigest()
    return compiled


def require_trace_ready(scene: dict[str, Any]) -> None:
    """Reject a compiled scene that cannot be passed to a production tracer."""
    report = scene.get("report", {})
    blockers = report.get("trace_blockers", [])
    if not report.get("native_trace_ready", False):
        blockers = [*blockers, "native production scene binding is unavailable"]
    if blockers:
        raise SceneCompileError("scene is not trace-ready: " + "; ".join(blockers))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compile a simulation-models IR into generic scene data."
    )
    parser.add_argument("ir", type=Path, help="provenance-checked scene IR JSON")
    parser.add_argument(
        "--source-root", type=Path, required=True, help="root used to resolve IR asset paths"
    )
    parser.add_argument(
        "--simtel-root",
        type=Path,
        help="explicit sim_telarray installation for cfg/CTA camera tables",
    )
    parser.add_argument("--output", type=Path, required=True, help="compiled generic-scene JSON")
    parser.add_argument(
        "--require-trace-ready",
        action="store_true",
        help="fail if unresolved optical geometry or response prevents production tracing",
    )
    args = parser.parse_args()
    try:
        ir = json.loads(args.ir.read_text(encoding="utf-8"))
        if not isinstance(ir, dict):
            raise SceneCompileError("IR root must be an object")
        scene = compile_scene(ir, args.source_root, simtel_root=args.simtel_root)
        if args.require_trace_ready:
            require_trace_ready(scene)
    except (OSError, json.JSONDecodeError, SceneCompileError) as error:
        raise SystemExit(f"scene compilation failed: {error}") from error
    args.output.write_text(json.dumps(scene, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
