"""Compile a provenance-checked simulation-models IR into generic optical_model data.

This adapter verifies the selected production and extracts documented mirror
footprints and focal-plane layouts. It does not invent dish sag, panel
normals, detector surfaces, or material behaviour. The output is an audited
handoff, not a trace-ready production optical_model.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from obdeect.camera_config import CameraConfigError, parse_camera_layout, parse_camera_layout_ecsv
from obdeect.model_import import ImportError as ModelImportError
from obdeect.model_import import component, record, resolve_model


class OpticalModelCompileError(ValueError):
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
            raise OpticalModelCompileError(
                f"segmentation line {line_number}: unsupported type {fields[0]}"
            )
        expected_fields = {5, 6} if kind in {"hex", "yhex"} else {5, 6, 7}
        if len(fields) not in expected_fields:
            raise OpticalModelCompileError(f"segmentation line {line_number}: wrong field count")
        try:
            count = int(fields[1])
            values = [float(value) for value in fields[2:]]
        except ValueError as error:
            raise OpticalModelCompileError(
                f"segmentation line {line_number}: invalid number"
            ) from error
        if count < 1 or not all(math.isfinite(value) for value in values):
            raise OpticalModelCompileError(
                f"segmentation line {line_number}: invalid count or value"
            )
        if kind in {"hex", "yhex"}:
            x_cm, y_cm, diameter_cm, rotation_deg = (*values, 0.0)[0:4]
            if count != 1 or diameter_cm <= 0.0:
                raise OpticalModelCompileError(
                    f"segmentation line {line_number}: invalid hex footprint"
                )
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
                raise OpticalModelCompileError(
                    f"segmentation line {line_number}: invalid ring footprint"
                )
            if count * span_deg > 360.0 + 1e-9:
                raise OpticalModelCompileError(
                    f"segmentation line {line_number}: ring exceeds full turn"
                )
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
        raise OpticalModelCompileError("segmentation contains no segments")
    return segments


def _number(value: object, context: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise OpticalModelCompileError(f"{context} must be a finite number")
    return float(value)


def _length_m(parameter: dict[str, Any], name: str, *, allow_zero: bool = False) -> float:
    unit = parameter.get("unit")
    if unit not in _UNIT_TO_M:
        raise OpticalModelCompileError(f"{name} has unsupported length unit {unit!r}")
    value = _number(parameter.get("value"), name) * _UNIT_TO_M[unit]
    if value < 0.0 or (value == 0.0 and not allow_zero):
        qualifier = "non-negative" if allow_zero else "positive"
        raise OpticalModelCompileError(f"{name} must be {qualifier}")
    return value


def _signed_length_m(parameter: dict[str, Any], name: str) -> float:
    unit = parameter.get("unit")
    if unit not in _UNIT_TO_M:
        raise OpticalModelCompileError(f"{name} has unsupported length unit {unit!r}")
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
        raise OpticalModelCompileError("fallback focal length must be finite and positive")
    facets: list[dict[str, Any]] = []
    for line_number, source_line in enumerate(contents.splitlines(), start=1):
        line = source_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        # Astropy ECSV keeps its column names in one un-commented header row.
        # The numerical columns that follow have the same documented layout as
        # sim_telarray mirror lists, so skip only this exact header form.
        if fields[:5] == ["mirror_x", "mirror_y", "mirror_diameter", "focal_length", "shape_type"]:
            continue
        if len(fields) < 5:
            raise OpticalModelCompileError(
                f"mirror list line {line_number}: expected at least 5 columns"
            )
        try:
            x_cm, y_cm, diameter_cm, focal_cm = (float(item) for item in fields[:4])
            shape_value = float(fields[4])
            shape_code = int(shape_value)
            # sim_telarray parses a sixth floating-point field when present.
            # A comment immediately after the required fields means no z was
            # supplied, rather than an invalid optical datum.
            z_cm = 0.0
            if len(fields) >= 6 and not fields[5].startswith("#"):
                z_cm = float(fields[5])
        except ValueError as error:
            raise OpticalModelCompileError(
                f"mirror list line {line_number}: invalid numeric field"
            ) from error
        if not all(
            math.isfinite(item) for item in (x_cm, y_cm, diameter_cm, focal_cm, shape_value, z_cm)
        ):
            raise OpticalModelCompileError(f"mirror list line {line_number}: non-finite value")
        if shape_value != shape_code:
            raise OpticalModelCompileError(
                f"mirror list line {line_number}: shape must be an integer"
            )
        if diameter_cm <= 0.0:
            raise OpticalModelCompileError(
                f"mirror list line {line_number}: diameter must be positive"
            )
        if shape_code not in _SHAPES:
            raise OpticalModelCompileError(
                f"mirror list line {line_number}: unsupported shape {shape_code}"
            )
        focal_m = focal_cm * 0.01
        if focal_m == 0.0:
            if fallback_focal_length_m is None:
                raise OpticalModelCompileError(
                    f"mirror list line {line_number}: zero focal length has no fallback"
                )
            focal_m = fallback_focal_length_m
        if focal_m <= 0.0:
            raise OpticalModelCompileError(
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
        raise OpticalModelCompileError("mirror list contains no facets")
    return facets


def parse_obscuration_cylinders(contents: str) -> list[dict[str, Any]]:
    """Parse finite opaque cylinders in the documented primary-mirror frame."""
    header: list[str] | None = None
    cylinders: list[dict[str, Any]] = []
    required = {"id", "x1", "y1", "z1", "x2", "y2", "z2", "diameter"}
    for line_number, source_line in enumerate(contents.splitlines(), start=1):
        line = source_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if header is None:
            header = fields
            if not required.issubset(header):
                raise OpticalModelCompileError("obscuration cylinders lack required columns")
            continue
        if len(fields) != len(header):
            raise OpticalModelCompileError(
                f"obscuration cylinder line {line_number}: wrong column count"
            )
        row = dict(zip(header, fields, strict=True))
        try:
            values = [float(row[name]) for name in ("x1", "y1", "z1", "x2", "y2", "z2", "diameter")]
        except ValueError as error:
            raise OpticalModelCompileError(
                f"obscuration cylinder line {line_number}: invalid number"
            ) from error
        if not all(math.isfinite(value) for value in values) or values[-1] <= 0:
            raise OpticalModelCompileError(
                f"obscuration cylinder line {line_number}: invalid geometry"
            )
        if math.dist(values[:3], values[3:6]) <= 1.0e-12:
            raise OpticalModelCompileError(f"obscuration cylinder line {line_number}: zero length")
        cylinders.append({
            "id": row["id"],
            "first_endpoint_m": values[:3],
            "second_endpoint_m": values[3:6],
            "diameter_m": values[6],
        })
    if header is None or not cylinders:
        raise OpticalModelCompileError("obscuration cylinder table is empty")
    return cylinders


def parse_wavelength_response(contents: str, name: str) -> list[dict[str, float]]:
    """Read a one-dimensional model response table without discarding angle data."""
    header: list[str] | None = None
    response: list[dict[str, float]] = []
    for line_number, source_line in enumerate(contents.splitlines(), start=1):
        line = source_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if header is None:
            header = fields
            if "wavelength" not in header or "reflectivity" not in header:
                raise OpticalModelCompileError(f"{name} lacks wavelength and reflectivity columns")
            if "incidence_angle" in header:
                raise OpticalModelCompileError(f"{name} has incidence-dependent response")
            continue
        if len(fields) != len(header):
            raise OpticalModelCompileError(f"{name} line {line_number}: wrong column count")
        row = dict(zip(header, fields, strict=True))
        try:
            wavelength, value = float(row["wavelength"]), float(row["reflectivity"])
        except ValueError as error:
            raise OpticalModelCompileError(
                f"{name} line {line_number}: invalid response"
            ) from error
        if (
            not math.isfinite(wavelength)
            or not math.isfinite(value)
            or wavelength <= 0
            or not 0 <= value <= 1
        ):
            raise OpticalModelCompileError(f"{name} line {line_number}: invalid response")
        response.append({"wavelength_nm": wavelength, "response": value})
    if len(response) < 2 or any(
        left["wavelength_nm"] > right["wavelength_nm"]
        for left, right in zip(response, response[1:])
    ):
        raise OpticalModelCompileError(f"{name} must have strictly increasing wavelengths")
    # Some released response tables repeat a wavelength for independent
    # measurements. A table-driven transport kernel needs one ordinate per
    # abscissa, so retain their arithmetic mean instead of depending on parser
    # row order.
    reduced: list[dict[str, float]] = []
    samples: list[int] = []
    for entry in response:
        if reduced and entry["wavelength_nm"] == reduced[-1]["wavelength_nm"]:
            samples[-1] += 1
            reduced[-1]["response"] += (entry["response"] - reduced[-1]["response"]) / samples[-1]
        else:
            reduced.append(entry)
            samples.append(1)
    return reduced


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
        raise OpticalModelCompileError("parabolic_dish must be boolean")
    for facet in facets:
        x, y, file_z = facet["centre_m"]
        radius = math.hypot(x, y)
        if file_z:
            height = abs(file_z)
        elif parabolic:
            height = radius * radius / (4 * dish_length)
        else:
            if radius > dish_length:
                raise OpticalModelCompileError(
                    "panel centre lies outside Davies-Cotton dish radius"
                )
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


def _even_polynomial_surface(parameter: dict[str, Any], name: str) -> list[float]:
    """Convert a physical-centimetre even polynomial into SI coefficients."""
    value, unit = parameter.get("value"), parameter.get("unit")
    if (
        not isinstance(value, list)
        or not isinstance(unit, list)
        or len(value) != len(unit)
        or not 1 <= len(value) <= 20
        or any(entry != "cm" for entry in unit)
    ):
        raise OpticalModelCompileError(f"{name} must be a physical-centimetre coefficient list")
    coefficients = []
    for index, coefficient in enumerate(value):
        if isinstance(coefficient, bool) or not isinstance(coefficient, (int, float)):
            raise OpticalModelCompileError(f"{name} has a non-numeric coefficient")
        converted = float(coefficient) * 0.01 ** (1 - 2 * index)
        if not math.isfinite(converted):
            raise OpticalModelCompileError(f"{name} has a non-finite coefficient")
        coefficients.append(converted)
    if any(value != 0.0 for value in coefficients[13:]):
        raise OpticalModelCompileError(f"{name} exceeds the native 13-coefficient surface contract")
    return coefficients[:13] + [0.0] * (13 - len(coefficients[:13]))


def _dual_reflector_surfaces(parameters: dict[str, Any]) -> dict[str, dict[str, Any]] | None:
    """Return imported SSTS-like aspheres without inferring absent geometry."""
    required = {
        "primary_mirror_parameters",
        "primary_mirror_diameter",
        "primary_mirror_hole_diameter",
        "secondary_mirror_parameters",
        "secondary_mirror_diameter",
        "secondary_mirror_hole_diameter",
        "focal_surface_parameters",
    }
    present = required & set(parameters)
    if not present:
        return None
    if present != required:
        raise OpticalModelCompileError(
            "dual-reflector optical prescription is incomplete: "
            + ", ".join(sorted(required - present))
        )

    def radial_scale(prefix: str) -> float | None:
        parameter = parameters.get(f"{prefix}_ref_radius")
        if parameter is None:
            return None
        return _length_m(parameter, f"{prefix}_ref_radius")

    def polynomial(prefix: str) -> tuple[list[float], float]:
        """Return SI sag coefficients and their explicit radial scale.

        sim_telarray accepts the SC prescription as ``z = R sum(a_i (r/R)^2i)``.
        The SST records use R=1 cm, which is algebraically identical to the
        physical-centimetre convention; SCT uses its 558.63-cm reference
        radius and must therefore retain the scale rather than expanding the
        coefficients as physical centimetre powers.
        """
        scale = radial_scale(prefix)
        if scale is None:
            return _even_polynomial_surface(
                parameters[f"{prefix}_parameters"], f"{prefix}_parameters"
            ), 1.0
        values = parameters[f"{prefix}_parameters"].get("value")
        units = parameters[f"{prefix}_parameters"].get("unit")
        if (
            not isinstance(values, list)
            or not isinstance(units, list)
            or len(values) != len(units)
            or not 1 <= len(values) <= 20
            or any(unit != "cm" for unit in units)
        ):
            raise OpticalModelCompileError(
                f"{prefix}_parameters must be a centimetre coefficient list"
            )
        if any(isinstance(value, bool) or not isinstance(value, (int, float)) for value in values):
            raise OpticalModelCompileError(f"{prefix}_parameters has a non-numeric coefficient")
        coefficient_m = [float(value) * scale for value in values]
        if not all(math.isfinite(value) for value in coefficient_m):
            raise OpticalModelCompileError(f"{prefix}_parameters has a non-finite coefficient")
        if any(value != 0.0 for value in coefficient_m[13:]):
            raise OpticalModelCompileError(
                f"{prefix}_parameters exceeds the native 13-coefficient surface contract"
            )
        return coefficient_m[:13] + [0.0] * (13 - len(coefficient_m[:13])), scale

    def surface(prefix: str) -> dict[str, Any]:
        outer = _length_m(parameters[f"{prefix}_diameter"], f"{prefix}_diameter") * 0.5
        inner = (
            _length_m(
                parameters[f"{prefix}_hole_diameter"], f"{prefix}_hole_diameter", allow_zero=True
            )
            * 0.5
        )
        if inner >= outer:
            raise OpticalModelCompileError(f"{prefix} hole must be smaller than its aperture")
        coefficient_m, scale = polynomial(prefix)
        return {
            "coefficient_m": coefficient_m,
            "radial_scale_m": scale,
            "inner_radius_m": inner,
            "outer_radius_m": outer,
        }

    primary = surface("primary_mirror")
    secondary = surface("secondary_mirror")
    focal_value, focal_scale = polynomial("focal_surface")
    return {
        "primary": primary,
        "secondary": secondary,
        "focal_surface": {
            "coefficient_m": focal_value,
            "radial_scale_m": focal_scale,
            "inner_radius_m": 0.0,
            "outer_radius_m": None,
        },
    }


def compile_optical_model(
    ir: dict[str, Any], source_root: Path, *, simtel_root: Path | None = None
) -> dict[str, Any]:
    """Compile ``obdeect.simulation-models-optical-model-ir.v1`` into generic optical model data."""
    if ir.get("format") != "obdeect.simulation-models-optical-model-ir.v1":
        raise OpticalModelCompileError("expected obdeect.simulation-models-optical-model-ir.v1")
    model = ir.get("model")
    version = ir.get("model_version")
    parameters = ir.get("parameters")
    assets = ir.get("assets")
    if not isinstance(model, str) or not model or not isinstance(version, str) or not version:
        raise OpticalModelCompileError("IR model identity is invalid")
    if not isinstance(parameters, dict) or not isinstance(assets, dict):
        raise OpticalModelCompileError("IR parameters and assets must be objects")
    root = source_root.resolve()
    # Keep the compiler's source-root contract aligned with the importer: a
    # released checkout may contain its data package in simulation-models/.
    nested_root = root / "simulation-models"
    if not (root / "model_parameters").is_dir() and (nested_root / "model_parameters").is_dir():
        root = nested_root
    records = ir.get("input_records", {})
    if not isinstance(records, dict):
        raise OpticalModelCompileError("IR input_records must be an object")
    try:
        source_ir = resolve_model(source_root, model, version)
    except ModelImportError as error:
        raise OpticalModelCompileError(f"cannot verify source production: {error}") from error
    for key in ("source_root", "input_records", "assets", "parameters"):
        if ir.get(key) != source_ir[key]:
            raise OpticalModelCompileError(f"IR {key} differs from the verified source production")
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
    # Some dual-mirror catalogues set every panel focal-length column to zero
    # and provide the telescope focal length as the documented common value.
    # Use it only as an explicit fallback, never as an inferred prescription.
    if fallback is None and isinstance(parameters.get("focal_length"), dict):
        fallback = _length_m(parameters["focal_length"], "focal_length")
    consumed = set()
    dual_surfaces = _dual_reflector_surfaces(parameters)
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
        if fallback is not None and "focal_length" in parameters:
            consumed.add("focal_length")
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
        raise OpticalModelCompileError("IR has no tracked primary mirror geometry asset")
    if "telescope_obscuration_cylinders" in verified_assets:
        primary["cylinder_obscurers"] = parse_obscuration_cylinders(
            verified_assets["telescope_obscuration_cylinders"].read_text(encoding="utf-8")
        )
        consumed.add("telescope_obscuration_cylinders")
    if "mirror_reflectivity" in verified_assets:
        try:
            primary["reflectivity"] = parse_wavelength_response(
                verified_assets["mirror_reflectivity"].read_text(encoding="utf-8"),
                "mirror_reflectivity",
            )
            consumed.add("mirror_reflectivity")
        except OpticalModelCompileError as error:
            if "incidence-dependent response" not in str(error):
                raise
    secondary = None
    if "secondary_mirror_segmentation" in verified_assets:
        secondary = {
            "kind": "segmented_footprints",
            "segments": parse_simtel_segmentation(
                verified_assets["secondary_mirror_segmentation"].read_text(encoding="utf-8")
            ),
        }
        consumed.add("secondary_mirror_segmentation")
    if dual_surfaces is not None:
        consumed.update({
            "primary_mirror_parameters",
            "primary_mirror_diameter",
            "primary_mirror_hole_diameter",
            "secondary_mirror_parameters",
            "secondary_mirror_diameter",
            "secondary_mirror_hole_diameter",
            "focal_surface_parameters",
        })
        consumed.update(
            name
            for name in (
                "primary_mirror_ref_radius",
                "secondary_mirror_ref_radius",
                "focal_surface_ref_radius",
            )
            if name in parameters
        )
        primary["aspheric_surface"] = dual_surfaces["primary"]
        secondary = {"kind": "aspheric_mirror", **dual_surfaces["secondary"]}
        if "secondary_mirror_reflectivity" in verified_assets:
            try:
                secondary["reflectivity"] = parse_wavelength_response(
                    verified_assets["secondary_mirror_reflectivity"].read_text(encoding="utf-8"),
                    "secondary_mirror_reflectivity",
                )
                consumed.add("secondary_mirror_reflectivity")
            except OpticalModelCompileError as error:
                if "incidence-dependent response" not in str(error):
                    raise
    camera = None
    nested_assets = {}
    unresolved_references = []
    camera_asset_name = next(
        (name for name in ("camera_config_file", "camera_pixel_layout") if name in verified_assets),
        None,
    )
    if camera_asset_name is not None:
        try:
            contents = verified_assets[camera_asset_name].read_text(encoding="utf-8")
            camera = (
                parse_camera_layout_ecsv(contents)
                if camera_asset_name == "camera_pixel_layout"
                else parse_camera_layout(contents)
            )
        except CameraConfigError as error:
            raise OpticalModelCompileError(str(error)) from error
        consumed.add(camera_asset_name)
        declared_pixels = parameters.get("camera_pixels", {}).get("value")
        if isinstance(declared_pixels, bool) or not isinstance(declared_pixels, int):
            raise OpticalModelCompileError("camera_pixels must be an integer count")
        if len(camera["pixels"]) != declared_pixels:
            raise OpticalModelCompileError(
                "pixel count differs: focal-plane element count differs from production record"
            )
        consumed.add("camera_pixels")
        for pixel_type in camera["pixel_types"]:
            for filename in pixel_type["response_files"]:
                if not component(filename):
                    raise OpticalModelCompileError(f"unsafe camera response filename: {filename}")
                path = root / "model_parameters" / "Files" / filename
                if not path.resolve().is_relative_to(root):
                    raise OpticalModelCompileError(
                        f"camera response path escapes source root: {filename}"
                    )
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
            raise OpticalModelCompileError(
                f"required field {name} is not supported by the compiler"
            )
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
                else "not compiled into an optical optical model",
            }
            for name in sorted(parameters)
        },
        "trace_blockers": [
            "run-specific panel alignment and distance are not compiled"
            if nominal_geometry
            else "facet surface normals and alignment are not compiled",
            "physical detector surfaces and materials remain deferred",
            *(
                ["primary incidence-dependent reflectivity remains deferred"]
                if "mirror_reflectivity" in parameters and "mirror_reflectivity" not in consumed
                else []
            ),
            *(
                ["quadrilateral obscurers remain deferred"]
                if "telescope_obscuration_quadrilaterals" in parameters
                else []
            ),
            *[f"camera response asset is missing: {name}" for name in unresolved_references],
            *(
                ["cylinder obscurers are not supported by axisymmetric native export"]
                if dual_surfaces is not None and primary.get("cylinder_obscurers")
                else []
            ),
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
    compiled_focal_surface = None
    if dual_surfaces is not None:
        compiled_focal_surface = dict(dual_surfaces["focal_surface"])
        # The camera layout bounds the only model-derived focal aperture.  A
        # native curved detector is intentionally not fabricated until its
        # coordinate transform and active boundary are represented by the
        # native optical_model format.
        if camera is not None:
            compiled_focal_surface["outer_radius_m"] = _camera_extent_m(camera)
    compiled = {
        "format": "obdeect.compiled-optical-model.v1",
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
    if isinstance(parameters.get("focal_length"), dict):
        compiled["focal_length_m"] = _length_m(parameters["focal_length"], "focal_length")
    if secondary is not None:
        compiled["secondary"] = secondary
    if compiled_focal_surface is not None:
        compiled["focal_surface"] = compiled_focal_surface
    if camera is not None:
        compiled["camera"] = camera
    canonical = json.dumps(compiled, sort_keys=True, separators=(",", ":")).encode()
    compiled["optical_model_sha256"] = hashlib.sha256(canonical).hexdigest()
    return compiled


def require_trace_ready(optical_model: dict[str, Any]) -> None:
    """Reject a compiled optical model that cannot be passed to a production tracer."""
    report = optical_model.get("report", {})
    blockers = report.get("trace_blockers", [])
    if not report.get("native_trace_ready", False):
        blockers = [*blockers, "native production optical model binding is unavailable"]
    if blockers:
        raise OpticalModelCompileError("optical model is not trace-ready: " + "; ".join(blockers))


def _camera_extent_m(camera: dict[str, Any]) -> float:
    """Bound every known entrance footprint, including its half diameter."""
    try:
        radii = {
            pixel_type["id"]: float(pixel_type["funnel_diameter_m"]) * 0.5
            for pixel_type in camera.get("pixel_types", [])
        }
        if any(not math.isfinite(radius) or radius <= 0 for radius in radii.values()):
            raise ValueError("invalid pixel entrance radius")
        extent = 0.0
        for pixel in camera["pixels"]:
            x, y = pixel["centre_xy_m"]
            if radii and pixel.get("type_id") not in radii:
                raise ValueError("pixel entrance type is undefined")
            extent = max(
                extent,
                math.hypot(float(x), float(y)) + radii.get(pixel.get("type_id"), 0.0),
            )
    except (KeyError, IndexError, TypeError, ValueError) as error:
        raise OpticalModelCompileError(
            "native export found an invalid focal-plane layout"
        ) from error
    if not math.isfinite(extent) or extent <= 0.0:
        raise OpticalModelCompileError("native export found an invalid focal-plane extent")
    return extent


def native_surface_rows(
    optical_model: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[tuple[float, float]]]:
    """Return model-derived planar surfaces for the native optical_model adapter.

    The exporter only accepts facets with explicit nominal centres and normals.
    It never reconstructs a missing telescope prescription. The focal surface
    is represented by the imported focal-plane extent at the selected focal
    length, so the native tracer can preserve detector-boundary losses.
    """
    report = optical_model.get("report", {})
    if report.get("facet_geometry_evidence", {}).get("normal_status") != "nominal_unperturbed":
        raise OpticalModelCompileError("native export requires explicit facet normals")
    facets = optical_model.get("primary", {}).get("facets", [])
    if not facets:
        raise OpticalModelCompileError("native export has no primary facets")
    rows: list[dict[str, Any]] = []
    for facet in facets:
        centre = facet.get("nominal_centre_m")
        normal = facet.get("nominal_normal")
        if not isinstance(centre, list) or not isinstance(normal, list):
            raise OpticalModelCompileError("native export found a facet without nominal placement")
        rows.append({
            "surface_id": int(facet["id"]),
            "role": "mirror",
            "shape": facet["shape"],
            "centre_m": centre,
            "normal": normal,
            "diameter_m": facet["diameter_m"],
            "focal_length_m": facet["focal_length_m"],
        })
    camera = optical_model.get("camera", {})
    # ``parse_simtel_mirror_list`` retains the focal length for every facet;
    # use that explicit catalogue value for the detector plane.  A compiled
    # optical_model intentionally does not copy the whole parameter table, so this is
    # the only value the native exporter is allowed to consume here.
    focal_length = optical_model.get("focal_length_m")
    if (
        not isinstance(focal_length, (int, float))
        or not math.isfinite(float(focal_length))
        or focal_length <= 0
    ):
        focal_lengths = [float(facet["focal_length_m"]) for facet in facets]
        focal_length = focal_lengths[0]
        if any(
            not math.isclose(value, focal_length, rel_tol=1e-9, abs_tol=1e-12)
            for value in focal_lengths
        ):
            raise OpticalModelCompileError("native export requires an explicit focal length")
    pixels = camera.get("pixels", [])
    if not pixels:
        raise OpticalModelCompileError("native export requires focal-plane layout")
    extent = _camera_extent_m(camera)
    rows.append({
        "surface_id": max(row["surface_id"] for row in rows) + 1,
        "role": "detector",
        "shape": "circle",
        "centre_m": [0.0, 0.0, float(focal_length)],
        "normal": [0.0, 0.0, 1.0],
        "diameter_m": 2.0 * extent,
    })
    obscurers = optical_model.get("primary", {}).get("cylinder_obscurers", [])
    if not isinstance(obscurers, list):
        raise OpticalModelCompileError("native export found invalid cylinder obscurers")
    next_surface_id = max(row["surface_id"] for row in rows) + 1
    native_obscurers = []
    for obscurer in obscurers:
        if not isinstance(obscurer, dict):
            raise OpticalModelCompileError("native export found invalid cylinder obscurer")
        try:
            first = [float(value) for value in obscurer["first_endpoint_m"]]
            second = [float(value) for value in obscurer["second_endpoint_m"]]
            diameter = float(obscurer["diameter_m"])
        except (KeyError, TypeError, ValueError) as error:
            raise OpticalModelCompileError(
                "native export found invalid cylinder obscurer"
            ) from error
        if (
            len(first) != 3
            or len(second) != 3
            or not all(math.isfinite(value) for value in (*first, *second, diameter))
            or diameter <= 0
            or math.dist(first, second) <= 1e-12
        ):
            raise OpticalModelCompileError("native export found invalid cylinder obscurer")
        native_obscurers.append({
            "surface_id": next_surface_id,
            "first": first,
            "second": second,
            "diameter_m": diameter,
        })
        next_surface_id += 1
    for row in rows:
        normal = row["normal"]
        # Stable in-plane orientation for polygonal apertures.  Project the
        # telescope x axis into the panel plane; use y when the panel normal
        # is parallel to x.  Circular surfaces do not consume this field.
        candidate = [1.0, 0.0, 0.0]
        projection = sum(candidate[index] * normal[index] for index in range(3))
        tangent = [candidate[index] - projection * normal[index] for index in range(3)]
        tangent_norm = math.sqrt(sum(value * value for value in tangent))
        if tangent_norm <= 1e-12:
            candidate = [0.0, 1.0, 0.0]
            projection = sum(candidate[index] * normal[index] for index in range(3))
            tangent = [candidate[index] - projection * normal[index] for index in range(3)]
            tangent_norm = math.sqrt(sum(value * value for value in tangent))
        row["tangent"] = [value / tangent_norm for value in tangent]
    reflectivity = optical_model.get("primary", {}).get("reflectivity", [])
    if not isinstance(reflectivity, list):
        raise OpticalModelCompileError("native export found invalid primary reflectivity")
    native_reflectivity = []
    for entry in reflectivity:
        if not isinstance(entry, dict):
            raise OpticalModelCompileError("native export found invalid primary reflectivity")
        try:
            wavelength = float(entry["wavelength_nm"])
            value = float(entry["response"])
        except (KeyError, TypeError, ValueError) as error:
            raise OpticalModelCompileError(
                "native export found invalid primary reflectivity"
            ) from error
        if (
            not math.isfinite(wavelength)
            or not math.isfinite(value)
            or wavelength <= 0
            or not 0 <= value <= 1
        ):
            raise OpticalModelCompileError("native export found invalid primary reflectivity")
        native_reflectivity.append((wavelength, value))
    if any(
        left[0] >= right[0] for left, right in zip(native_reflectivity, native_reflectivity[1:])
    ):
        raise OpticalModelCompileError("native export found unordered primary reflectivity")
    return rows, native_obscurers, native_reflectivity


def write_native_optical_model(optical_model: dict[str, Any], output: Path) -> None:
    """Write the strict, dependency-free native optical model surface table."""
    primary = optical_model.get("primary", {})
    secondary = optical_model.get("secondary", {})
    focal = optical_model.get("focal_surface", {})
    if (
        isinstance(primary, dict)
        and isinstance(secondary, dict)
        and isinstance(focal, dict)
        and isinstance(primary.get("aspheric_surface"), dict)
        and secondary.get("kind") == "aspheric_mirror"
    ):
        _write_native_axisymmetric_optical_model(optical_model, output)
        return
    rows, obscurers, reflectivity = native_surface_rows(optical_model)
    lines = [
        "obdeect-optical-model-v1",
        "surface_id,role,shape,cx_m,cy_m,cz_m,nx,ny,nz,tx,ty,tz,diameter_m,focal_length_m",
    ]
    provenance = optical_model.get("provenance", {})
    model = provenance.get("model")
    version = provenance.get("model_version")
    records = provenance.get("input_records", {})
    if (
        not isinstance(model, str)
        or not isinstance(version, str)
        or "," in model
        or "," in version
        or not isinstance(records, dict)
    ):
        raise OpticalModelCompileError("native export requires model provenance")
    content_hash = optical_model.get("optical_model_sha256")
    if not isinstance(content_hash, str) or len(content_hash) != 64:
        raise OpticalModelCompileError("native export requires optical_model_sha256")
    lines.insert(1, f"provenance,{model},{version},{content_hash}")
    for row in rows:
        centre = row["centre_m"]
        normal = row["normal"]
        tangent = row["tangent"]
        lines.append(
            ",".join([
                str(row["surface_id"]),
                row["role"],
                row["shape"],
                *(f"{value:.17g}" for value in (*centre, *normal, *tangent, row["diameter_m"])),
                f"{row.get('focal_length_m', 0.0):.17g}",
            ])
        )
    for obscurer in obscurers:
        lines.append(
            ",".join([
                "obscurer_cylinder",
                str(obscurer["surface_id"]),
                *(
                    f"{value:.17g}"
                    for value in (*obscurer["first"], *obscurer["second"], obscurer["diameter_m"])
                ),
            ])
        )
    for wavelength, value in reflectivity:
        lines.append(f"primary_reflectivity,{wavelength:.17g},{value:.17g}")
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_native_axisymmetric_optical_model(optical_model: dict[str, Any], output: Path) -> None:
    """Write an exact rotationally symmetric SST/SCT optical prescription.

    The source model defines these surfaces as an even polynomial around the
    telescope optical axis.  This format carries the radius convention and
    never replaces them with a faceted approximation.
    """
    provenance = optical_model.get("provenance", {})
    model, version = provenance.get("model"), provenance.get("model_version")
    content_hash = optical_model.get("optical_model_sha256")
    if (
        not isinstance(model, str)
        or not isinstance(version, str)
        or not isinstance(content_hash, str)
    ):
        raise OpticalModelCompileError("native export requires model provenance")
    if len(content_hash) != 64 or "," in model or "," in version:
        raise OpticalModelCompileError("native export requires valid model provenance")

    def row(role: str, surface: dict[str, Any]) -> str:
        coefficients = surface.get("coefficient_m")
        outer = surface.get("outer_radius_m")
        inner = surface.get("inner_radius_m", 0.0)
        scale = surface.get("radial_scale_m", 1.0)
        if (
            not isinstance(coefficients, list)
            or len(coefficients) != 13
            or not all(
                isinstance(value, (int, float)) and math.isfinite(value) for value in coefficients
            )
            or not all(
                isinstance(value, (int, float)) and math.isfinite(value)
                for value in (inner, outer, scale)
            )
            or inner < 0
            or outer <= inner
            or scale <= 0
        ):
            raise OpticalModelCompileError(f"invalid {role} aspheric surface")
        # The native intersection is z = vertex + polynomial(r), so retain
        # the physical vertex separately from the zeroed polynomial constant.
        local = [float(value) for value in coefficients]
        vertex = local[0]
        local[0] = 0.0
        return ",".join([
            role,
            *(f"{value:.17g}" for value in (vertex, inner, outer, scale, *local)),
        ])

    primary_source = optical_model["primary"]
    if primary_source.get("cylinder_obscurers"):
        raise OpticalModelCompileError(
            "native axisymmetric export cannot represent cylinder obscurers"
        )
    primary = primary_source["aspheric_surface"]
    secondary = optical_model["secondary"]
    focal = optical_model.get("focal_surface")
    if not isinstance(focal, dict) or focal.get("outer_radius_m") is None:
        raise OpticalModelCompileError("native dual-mirror export requires a bounded focal surface")
    lines = [
        "obdeect-axisymmetric-optical-model-v1",
        f"provenance,{model},{version},{content_hash}",
        "role,vertex_z_m,inner_radius_m,outer_radius_m,radial_scale_m,c0_m,c1_m,c2_m,c3_m,c4_m,c5_m,c6_m,c7_m,c8_m,c9_m,c10_m,c11_m,c12_m",
        row("primary", primary),
        row("secondary", secondary),
        row("detector", focal),
    ]
    for role, surface, response_source in (
        ("primary", primary, primary_source),
        ("secondary", secondary, secondary),
    ):
        response = response_source.get("reflectivity", [])
        if not isinstance(response, list):
            raise OpticalModelCompileError(f"invalid {role} reflectivity")
        for entry in response:
            if not isinstance(entry, dict):
                raise OpticalModelCompileError(f"invalid {role} reflectivity")
            try:
                wavelength = float(entry["wavelength_nm"])
                value = float(entry["response"])
            except (KeyError, TypeError, ValueError) as error:
                raise OpticalModelCompileError(f"invalid {role} reflectivity") from error
            if (
                not math.isfinite(wavelength)
                or not math.isfinite(value)
                or wavelength <= 0
                or not 0 <= value <= 1
            ):
                raise OpticalModelCompileError(f"invalid {role} reflectivity")
            lines.append(f"{role}_reflectivity,{wavelength:.17g},{value:.17g}")
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compile a simulation-models IR into generic optical model data."
    )
    parser.add_argument(
        "--input", type=Path, required=True, help="provenance-checked optical model IR JSON"
    )
    parser.add_argument(
        "--source-root", type=Path, required=True, help="root used to resolve IR asset paths"
    )
    parser.add_argument(
        "--simtel-root",
        type=Path,
        help="explicit sim_telarray installation for cfg/CTA camera tables",
    )
    parser.add_argument(
        "--output", type=Path, required=True, help="compiled generic-optical model JSON"
    )
    parser.add_argument("--native-output", type=Path, help="optional native surface table")
    parser.add_argument(
        "--require-trace-ready",
        action="store_true",
        help="fail if unresolved optical geometry or response prevents production tracing",
    )
    args = parser.parse_args()
    try:
        ir = json.loads(args.input.read_text(encoding="utf-8"))
        if not isinstance(ir, dict):
            raise OpticalModelCompileError("IR root must be an object")
        optical_model = compile_optical_model(ir, args.source_root, simtel_root=args.simtel_root)
        if args.require_trace_ready:
            require_trace_ready(optical_model)
        if args.native_output:
            write_native_optical_model(optical_model, args.native_output)
    except (OSError, json.JSONDecodeError, OpticalModelCompileError) as error:
        raise SystemExit(f"optical model compilation failed: {error}") from error
    args.output.write_text(
        json.dumps(optical_model, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Compiled {optical_model['provenance']['model']} geometry into {args.output}")
    if args.native_output:
        print(f"Wrote native surface table to {args.native_output}")


if __name__ == "__main__":
    main()
