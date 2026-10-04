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

from obdeect.camera_config import (
    CameraConfigError,
    parse_camera_layout,
    parse_camera_layout_ecsv,
    parse_camera_pixel_types,
)
from obdeect.camera_surfaces import compile_camera_surfaces
from obdeect.model_import import ImportError as ModelImportError
from obdeect.model_import import component, record, resolve_model


class OpticalModelCompileError(ValueError):
    """The provenance IR cannot be compiled without guessing optical data."""


_SHAPES = {0: "circle", 1: "hexagon_flat_y", 2: "square", 3: "hexagon_flat_x"}
_UNIT_TO_M = {"m": 1.0, "cm": 0.01, "mm": 0.001}
_RESPONSE_METADATA = {"reflectivity_rms", "reflectivity_min", "reflectivity_max"}


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
                "rotation_deg": rotation_deg + (90.0 if kind == "yhex" else 0.0),
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
                    "start_deg": start_deg + index * span_deg,
                    "span_deg": span_deg,
                    "gap_m": gap_cm * 0.01,
                })
    if not segments:
        raise OpticalModelCompileError("segmentation contains no segments")
    return segments


def parse_model_segmentation(parameter: dict[str, Any]) -> list[dict[str, Any]]:
    """Convert structured production footprints through the documented simtel parser."""
    groups = parameter.get("value")
    if not isinstance(groups, list) or not groups:
        raise OpticalModelCompileError("structured segmentation must contain footprint groups")
    lines = []
    for group in groups:
        if not isinstance(group, dict):
            raise OpticalModelCompileError("structured segmentation group must be an object")
        kind = group.get("kind")
        required = {"kind", "count"} | (
            {"x", "y", "diameter"} if kind in {"hex", "yhex"} else {"r_min", "r_max", "dphi"}
        )
        optional = {"rotation"} if kind in {"hex", "yhex"} else {"phi0", "gap"}
        count = group.get("count")
        if (
            kind not in {"hex", "yhex", "ring"}
            or not required.issubset(group)
            or set(group) - required - optional
            or isinstance(count, bool)
            or not isinstance(count, int)
        ):
            raise OpticalModelCompileError("invalid structured segmentation fields")

        def quantity(name: str, *, angle: bool = False) -> float:
            value = group.get(name)
            if value is None and name in optional:
                return 0.0
            if not isinstance(value, dict) or set(value) != {"value", "unit"}:
                raise OpticalModelCompileError(f"segmentation {name} requires value and unit")
            if angle:
                if value["unit"] != "deg":
                    raise OpticalModelCompileError(f"segmentation {name} requires degrees")
                return _number(value["value"], f"segmentation {name}")
            return _signed_length_m(value, f"segmentation {name}") * 100

        names = (
            ("x", "y", "diameter", "rotation")
            if kind != "ring"
            else ("r_min", "r_max", "dphi", "phi0", "gap")
        )
        values = [quantity(name, angle=name in {"rotation", "dphi", "phi0"}) for name in names]
        lines.append(" ".join([kind, str(count), *(format(value, ".17g") for value in values)]))
    return parse_simtel_segmentation("\n".join(lines))


def _secondary_segment_frame(segments: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Apply simtel secondary's pi rotation about y to finite local footprints."""
    transformed = []
    for segment in segments:
        if segment["shape"] == "hexagon":
            transformed.append({
                **segment,
                "centre_xy_m": [-segment["centre_xy_m"][0], segment["centre_xy_m"][1]],
                "rotation_deg": 180 - segment["rotation_deg"],
            })
        else:
            transformed.append({
                **segment,
                "start_deg": 180 - segment["start_deg"] - segment["span_deg"],
                "gap_at_start": True,
            })
    return transformed


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


def _validated_mirror_response(values: object, name: str) -> list[dict[str, float]]:
    """Validate spectral knots or a rectangular spectral/incidence grid."""
    if not isinstance(values, list) or len(values) < 2:
        raise OpticalModelCompileError(f"trace model requires {name} reflectivity")
    angular = any(isinstance(entry, dict) and "incidence_angle_deg" in entry for entry in values)
    result = []
    for entry in values:
        required = {"wavelength_nm", "response"} | ({"incidence_angle_deg"} if angular else set())
        if (
            not isinstance(entry, dict)
            or not required.issubset(entry)
            or set(entry) - required - _RESPONSE_METADATA
        ):
            raise OpticalModelCompileError(f"invalid {name} reflectivity columns")
        parsed = {key: _number(value, f"{name} {key}") for key, value in entry.items()}
        if (
            parsed["wavelength_nm"] <= 0
            or not 0 <= parsed["response"] <= 1
            or (angular and not 0 <= parsed["incidence_angle_deg"] <= 90)
            or any(parsed[key] < 0 for key in _RESPONSE_METADATA & parsed.keys())
        ):
            raise OpticalModelCompileError(f"invalid {name} reflectivity range")
        result.append(parsed)
    if angular:
        wavelengths = sorted({entry["wavelength_nm"] for entry in result})
        angles = sorted({entry["incidence_angle_deg"] for entry in result})
        keys = {(entry["wavelength_nm"], entry["incidence_angle_deg"]) for entry in result}
        if (
            len(wavelengths) < 2
            or len(angles) < 2
            or len(keys) != len(result)
            or len(result) != len(wavelengths) * len(angles)
        ):
            raise OpticalModelCompileError(
                f"{name} requires a complete unique incidence-angle grid"
            )
        result.sort(key=lambda entry: (entry["wavelength_nm"], entry["incidence_angle_deg"]))
    elif any(
        first["wavelength_nm"] >= second["wavelength_nm"]
        for first, second in zip(result, result[1:])
    ):
        raise OpticalModelCompileError(f"trace model found unordered {name} reflectivity")
    return result


def parse_wavelength_response(
    contents: str, name: str, value_column: str = "reflectivity"
) -> list[dict[str, float]]:
    """Read spectral response, preserving a complete incidence-angle grid when supplied."""
    header: list[str] | None = None
    response: list[dict[str, float]] = []
    for line_number, source_line in enumerate(contents.splitlines(), start=1):
        line = source_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if header is None:
            header = fields
            if set(header) - _RESPONSE_METADATA not in (
                {"wavelength", value_column},
                {"wavelength", "incidence_angle", value_column},
            ) or len(header) != len(set(header)):
                raise OpticalModelCompileError(
                    f"{name} lacks wavelength and {value_column} columns"
                )
            continue
        if len(fields) != len(header):
            raise OpticalModelCompileError(f"{name} line {line_number}: wrong column count")
        row = dict(zip(header, fields, strict=True))
        try:
            wavelength, value = float(row["wavelength"]), float(row[value_column])
            metadata = {key: float(row[key]) for key in _RESPONSE_METADATA & row.keys()}
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
        entry = {"wavelength_nm": wavelength, "response": value}
        for key, measured in metadata.items():
            entry[key] = _number(measured, f"{name} {key}")
            if entry[key] < 0:
                raise OpticalModelCompileError(f"{name} line {line_number}: invalid {key}")
        if "incidence_angle" in row:
            try:
                entry["incidence_angle_deg"] = float(row["incidence_angle"])
            except ValueError as error:
                raise OpticalModelCompileError(
                    f"{name} line {line_number}: invalid angle"
                ) from error
        response.append(entry)
    if header is not None and "incidence_angle" in header:
        return _validated_mirror_response(response, name)
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
            # Metadata means describe the averaged rows, not propagated errors.
            for key in entry.keys() - {"wavelength_nm"}:
                reduced[-1][key] += (entry[key] - reduced[-1][key]) / samples[-1]
        else:
            reduced.append(entry)
            samples.append(1)
    return reduced


def parse_angular_response(contents: str, name: str) -> list[dict[str, float]]:
    """Read measured lightguide response without adding geometric losses twice."""
    header = None
    result = []
    for line in contents.splitlines():
        fields = line.split("#", 1)[0].split()
        if not fields:
            continue
        if header is None:
            header = fields
            if set(header) not in (
                {"incidence_angle", "efficiency"},
                {"incidence_angle", "efficiency", "inverse_cosine"},
            ):
                raise OpticalModelCompileError(f"{name} has unsupported angular columns")
            continue
        if len(fields) != len(header):
            raise OpticalModelCompileError(f"{name} has wrong angular column count")
        try:
            values = dict(zip(header, map(float, fields), strict=True))
        except ValueError as error:
            raise OpticalModelCompileError(f"{name} has invalid angular response") from error
        angle, response = values["incidence_angle"], values["efficiency"]
        if (
            not all(math.isfinite(value) for value in values.values())
            or not 0 <= angle <= 90
            or not 0 <= response <= 1
        ):
            raise OpticalModelCompileError(f"{name} has invalid angular range")
        entry = {"incidence_angle_deg": angle, "response": response}
        if "inverse_cosine" in values:
            if angle == 90 or not math.isclose(
                values["inverse_cosine"], 1 / math.cos(math.radians(angle)), rel_tol=1e-7
            ):
                raise OpticalModelCompileError(f"{name} inverse_cosine metadata differs from angle")
            entry["inverse_cosine"] = values["inverse_cosine"]
        result.append(entry)
    if len(result) < 2 or any(
        a["incidence_angle_deg"] >= b["incidence_angle_deg"] for a, b in zip(result, result[1:])
    ):
        raise OpticalModelCompileError(f"{name} requires increasing incidence angles")
    return result


def _apply_mirror_degradation(source: dict[str, Any], parameter: dict[str, Any], name: str) -> None:
    factor = _number(parameter.get("value"), name)
    if parameter.get("unit") not in (None, "", "null") or not 0 <= factor <= 1:
        raise OpticalModelCompileError(f"{name} requires a dimensionless response fraction")
    if "reflectivity" not in source:
        raise OpticalModelCompileError(f"{name} lacks a mirror response to scale")
    source["reflectivity"] = [
        {
            key: value * factor if key == "response" or key in _RESPONSE_METADATA else value
            for key, value in entry.items()
        }
        for entry in source["reflectivity"]
    ]
    source["degradation_applied"] = {
        "factor": factor,
        "source_parameter": name,
        "semantics": "multiplicative_response_once",
    }


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
        # A nominal facet normal is the bisector of the reverse incident
        # direction (+z for a star) and the direction from the facet centre
        # to the focal point.  Derive it from the final placement, including
        # mirror_offset, so every chief ray reaches the compiled focal plane.
        focal_distance = math.hypot(radius, focal_length - offset - z)
        if not math.isfinite(focal_distance) or focal_distance <= 0.0:
            raise OpticalModelCompileError("facet has no finite focal-point direction")
        outgoing_x = -x / focal_distance
        outgoing_y = -y / focal_distance
        outgoing_z = (focal_length - offset - z) / focal_distance
        normal_length = math.sqrt(
            outgoing_x * outgoing_x + outgoing_y * outgoing_y + (outgoing_z + 1.0) ** 2
        )
        if not math.isfinite(normal_length) or normal_length <= 0.0:
            raise OpticalModelCompileError("facet has no finite nominal normal")
        nx = outgoing_x / normal_length
        ny = outgoing_y / normal_length
        facet["nominal_centre_m"] = [x, y, z]
        facet["nominal_normal"] = [nx, ny, (outgoing_z + 1.0) / normal_length]


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
    # sim_imaging.c::tel_setup_secondary rotates the local mirror frame by pi.
    # Export telescope-frame sag, retaining the vertex offset from sim_config.c.
    secondary_coefficients = secondary["coefficient_m"]
    if secondary_coefficients[0] < primary["coefficient_m"][0]:
        secondary_coefficients[0] *= -1
    secondary_coefficients[1:] = [-value for value in secondary_coefficients[1:]]
    focal_value, focal_scale = polynomial("focal_surface")
    overall_offset = (
        _signed_length_m(parameters["mirror_offset"], "mirror_offset")
        if "mirror_offset" in parameters
        else 0.0
    )
    for coefficients in (primary["coefficient_m"], secondary_coefficients, focal_value):
        coefficients[0] -= overall_offset
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


def compile_optical_model(ir: dict[str, Any], source_root: Path) -> dict[str, Any]:
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
    for key in (
        "source_root",
        "input_records",
        "assets",
        "parameters",
        "source_parameter_coverage",
    ):
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
    if "mirror_list" in verified_assets and dual_surfaces is None:
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
    elif "primary_mirror_segmentation" in parameters:
        if "primary_mirror_segmentation" in verified_assets:
            segments = parse_simtel_segmentation(
                verified_assets["primary_mirror_segmentation"].read_text(encoding="utf-8")
            )
        else:
            segments = parse_model_segmentation(parameters["primary_mirror_segmentation"])
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
    if "mirror_degraded_reflection" in parameters:
        _apply_mirror_degradation(
            primary, parameters["mirror_degraded_reflection"], "mirror_degraded_reflection"
        )
        consumed.add("mirror_degraded_reflection")
    secondary = None
    if "secondary_mirror_segmentation" in parameters:
        segmentation = parameters["secondary_mirror_segmentation"]
        if "secondary_mirror_segmentation" in verified_assets:
            segments = parse_simtel_segmentation(
                verified_assets["secondary_mirror_segmentation"].read_text(encoding="utf-8")
            )
        elif segmentation["value"] is None:
            segments = []
        else:
            segments = parse_model_segmentation(segmentation)
        secondary = {"kind": "segmented_footprints", "segments": _secondary_segment_frame(segments)}
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
        secondary = {
            "kind": "aspheric_mirror",
            **dual_surfaces["secondary"],
            "segments": secondary.get("segments", []) if secondary is not None else [],
        }
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
    if secondary is not None and "secondary_mirror_degraded_reflection" in parameters:
        _apply_mirror_degradation(
            secondary,
            parameters["secondary_mirror_degraded_reflection"],
            "secondary_mirror_degraded_reflection",
        )
        consumed.add("secondary_mirror_degraded_reflection")
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
        if "camera_pixel_types" in parameters:
            try:
                camera["pixel_types"] = parse_camera_pixel_types(parameters["camera_pixel_types"])
            except CameraConfigError as error:
                raise OpticalModelCompileError(str(error)) from error
            type_ids = {entry["id"] for entry in camera["pixel_types"]}
            if any(pixel["type_id"] not in type_ids for pixel in camera["pixels"]):
                raise OpticalModelCompileError(
                    "camera layout references undefined physical pixel type"
                )
            consumed.add("camera_pixel_types")
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
                    unresolved_references.append(filename)
    if camera is not None:
        if "camera_filter" in verified_assets:
            camera["filter_response"] = parse_wavelength_response(
                verified_assets["camera_filter"].read_text(encoding="utf-8"),
                "camera_filter",
                "transmission",
            )
            consumed.add("camera_filter")
        if "lightguide_efficiency_vs_incidence_angle" in verified_assets:
            camera["lightguide_response"] = parse_angular_response(
                verified_assets["lightguide_efficiency_vs_incidence_angle"].read_text(
                    encoding="utf-8"
                ),
                "lightguide_efficiency_vs_incidence_angle",
            )
            consumed.add("lightguide_efficiency_vs_incidence_angle")
        if "camera_transmission" in parameters:
            fraction = _number(parameters["camera_transmission"]["value"], "camera_transmission")
            if not 0 <= fraction <= 1:
                raise OpticalModelCompileError("camera_transmission must be a response fraction")
            camera["transmission"] = fraction
            consumed.add("camera_transmission")
    deferred = sorted(set(parameters) - consumed)
    for name in deferred:
        if parameters[name].get("required_for_trace") is True:
            raise OpticalModelCompileError(
                f"required field {name} is not supported by the compiler"
            )
    report = {
        "trace_ready": False,
        "native_trace_ready": False,
        "production_trace_ready": False,
        "consumed": sorted(consumed),
        "deferred": deferred,
        "unsupported": sorted(
            name
            for name, entry in ir.get("source_parameter_coverage", {}).items()
            if entry["disposition"] == "unsupported"
        ),
        "source_parameter_coverage": ir.get("source_parameter_coverage", {}),
        "field_coverage": {
            name: {
                "disposition": "consumed" if name in consumed else "deferred",
                "reason": "parsed into geometry"
                if name in consumed
                else "not applied by the optical tracer",
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
        # embedded trace-model representation.
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
    if camera is not None and camera.get("pixel_types") and "camera_pixel_types" in consumed:
        focus = parameters.get("focus_offset")
        if focus is None:
            report["trace_blockers"].append("physical camera placement lacks explicit focus offset")
        else:
            values, units = focus.get("value"), focus.get("unit")
            if (
                not isinstance(values, list)
                or len(values) != 4
                or not isinstance(units, list)
                or len(units) != 4
                or units[0] not in _UNIT_TO_M
                or any(_number(value, "focus_offset") != 0 for value in values[2:])
            ):
                raise OpticalModelCompileError(
                    "camera placement requires a resolved static focus offset"
                )
            offset = _number(values[0], "focus_offset") * _UNIT_TO_M[units[0]]
            if compiled_focal_surface is not None:
                compiled_focal_surface["coefficient_m"][0] -= offset
                focal = compiled_focal_surface
                mode = parameters.get("pixels_parallel", {}).get("value")
            else:
                overall = _signed_length_m(parameters["mirror_offset"], "mirror_offset")
                focal = {
                    "coefficient_m": [compiled["focal_length_m"] + offset - overall],
                    "radial_scale_m": 1.0,
                }
                mode = 1
                compiled["detector_vertex_z_m"] = focal["coefficient_m"][0]
            try:
                camera.update(
                    compile_camera_surfaces(camera, focal, mode, reflected=dual_surfaces is None)
                )
            except CameraConfigError as error:
                raise OpticalModelCompileError(str(error)) from error
            consumed.add("focus_offset")
            if "pixels_parallel" in parameters:
                consumed.add("pixels_parallel")
            camera["physical_geometry_convention"] = {
                "source": "sim_telarray/common/sim_imaging.c::camera_setup_inclined_pixels",
                "orientation_mode": mode,
                "prime_focus_pi_rotation": dual_surfaces is None,
                "boundary": "physical_entrance",
                "guide_wall_profile": "unavailable",
            }
            report["consumed"] = sorted(consumed)
            report["deferred"] = sorted(set(parameters) - consumed)
            for name in ("focus_offset", "pixels_parallel"):
                if name in report["field_coverage"] and name in consumed:
                    report["field_coverage"][name] = {
                        "disposition": "consumed",
                        "reason": "physical camera placement",
                    }
            report["camera_layout_evidence"]["surface_status"] = "physical_entrance_planes"
            report["trace_blockers"] = [
                "concentrator wall profiles, windows and detector material response remain deferred"
                if item == "physical detector surfaces and materials remain deferred"
                else item
                for item in report["trace_blockers"]
            ]
    try:
        compiled["trace_model"] = build_trace_model(compiled)
        compiled["report"]["trace_ready"] = True
        compiled["report"]["native_trace_ready"] = True
        compiled["plot_geometry"] = build_plot_geometry(compiled["trace_model"], report)
    except OpticalModelCompileError as error:
        compiled["report"]["trace_blockers"].append(str(error))
    canonical = json.dumps(
        compiled, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode()
    compiled["optical_model_sha256"] = hashlib.sha256(canonical).hexdigest()
    return compiled


def build_plot_geometry(trace_model: dict[str, Any], report: dict[str, Any]) -> dict[str, Any]:
    """Index only the finite geometry used by native tracing for diagnostic plots."""
    components = []
    if trace_model["kind"] == "segmented":
        for role, field in (
            ("primary", "primary_facets"),
            ("detector", "detector_surfaces"),
            ("opaque_cylinder", "cylinder_obscurers"),
        ):
            components.extend(
                {"id": item["id"], "role": role, "source": f"trace_model.{field}"}
                for item in trace_model.get(field, [])
            )
    else:
        components = [
            {
                "id": trace_model.get(f"{role}_surface_id", index),
                "role": role,
                "source": f"trace_model.{role}",
            }
            for index, role in enumerate(("primary", "secondary", "detector"))
            if role != "detector" or not trace_model.get("detector_surfaces")
        ]
        components.extend(
            {"id": item["id"], "role": "detector", "source": "trace_model.detector_surfaces"}
            for item in trace_model.get("detector_surfaces", [])
        )
    unavailable = []
    if any("physical detector surfaces and materials" in item for item in report["trace_blockers"]):
        unavailable.append("physical_pixel_boundaries")
    if any("quadrilateral obscurers" in item for item in report["trace_blockers"]):
        unavailable.append("quadrilateral_obscurers")
    return {
        "schema_version": 1,
        "frame": {"origin": "telescope_optical_reference", "axes": "+x,+y,+z", "unit": "m"},
        "components": components,
        "unavailable_roles": unavailable,
    }


def require_trace_ready(optical_model: dict[str, Any]) -> None:
    """Reject a compiled optical model that cannot be passed to a production tracer."""
    report = optical_model.get("report", {})
    blockers = report.get("trace_blockers", [])
    if not report.get("trace_ready", False):
        blockers = [*blockers, "trace-model binding is unavailable"]
    if blockers:
        raise OpticalModelCompileError("optical model is not trace-ready: " + "; ".join(blockers))


def _camera_extent_m(camera: dict[str, Any]) -> float:
    """Bound every known entrance footprint, including its half diameter."""
    try:
        radii = {
            pixel_type["id"]: float(pixel_type["funnel_diameter_m"])
            * 0.5
            * {0: 1.0, 1: 2 / math.sqrt(3), 2: math.sqrt(2), 3: 2 / math.sqrt(3)}.get(
                pixel_type.get("funnel_shape_code", 0), 1.0
            )
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
        raise OpticalModelCompileError("trace model found an invalid focal-plane layout") from error
    if not math.isfinite(extent) or extent <= 0.0:
        raise OpticalModelCompileError("trace model found an invalid focal-plane extent")
    return extent


def trace_surface_rows(
    optical_model: dict[str, Any], *, include_response: bool = True
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[tuple[float, float]]]:
    """Return model-derived planar surfaces for the embedded trace model.

    The exporter only accepts facets with explicit nominal centres and normals.
    It never reconstructs a missing telescope prescription. The focal surface
    is represented by the imported focal-plane extent at the selected focal
    length, so the tracer can preserve detector-boundary losses.
    """
    report = optical_model.get("report", {})
    if report.get("facet_geometry_evidence", {}).get("normal_status") != "nominal_unperturbed":
        raise OpticalModelCompileError("trace model requires explicit facet normals")
    facets = optical_model.get("primary", {}).get("facets", [])
    if not facets:
        raise OpticalModelCompileError("trace model has no primary facets")
    rows: list[dict[str, Any]] = []
    for facet in facets:
        centre = facet.get("nominal_centre_m")
        normal = facet.get("nominal_normal")
        if not isinstance(centre, list) or not isinstance(normal, list):
            raise OpticalModelCompileError("trace model found a facet without nominal placement")
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
    # the only value the trace-model builder is allowed to consume here.
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
            raise OpticalModelCompileError("trace model requires an explicit focal length")
    pixels = camera.get("pixels", [])
    if not pixels:
        raise OpticalModelCompileError("trace model requires focal-plane layout")
    detector_id = max(row["surface_id"] for row in rows) + 1
    physical = camera.get("entrance_surfaces", [])
    if physical:
        rows.extend(
            {**plane, "surface_id": detector_id + index, "role": "detector"}
            for index, plane in enumerate(physical)
        )
    else:
        extent = _camera_extent_m(camera)
        rows.append({
            "surface_id": detector_id,
            "role": "detector",
            "shape": "circle",
            "centre_m": [0.0, 0.0, float(optical_model.get("detector_vertex_z_m", focal_length))],
            "normal": [0.0, 0.0, 1.0],
            "diameter_m": 2.0 * extent,
        })
    obscurers = optical_model.get("primary", {}).get("cylinder_obscurers", [])
    if not isinstance(obscurers, list):
        raise OpticalModelCompileError("trace model found invalid cylinder obscurers")
    next_surface_id = max(row["surface_id"] for row in rows) + 1
    trace_obscurers = []
    for obscurer in obscurers:
        if not isinstance(obscurer, dict):
            raise OpticalModelCompileError("trace model found invalid cylinder obscurer")
        try:
            first = [float(value) for value in obscurer["first_endpoint_m"]]
            second = [float(value) for value in obscurer["second_endpoint_m"]]
            diameter = float(obscurer["diameter_m"])
        except (KeyError, TypeError, ValueError) as error:
            raise OpticalModelCompileError("trace model found invalid cylinder obscurer") from error
        if (
            len(first) != 3
            or len(second) != 3
            or not all(math.isfinite(value) for value in (*first, *second, diameter))
            or diameter <= 0
            or math.dist(first, second) <= 1e-12
        ):
            raise OpticalModelCompileError("trace model found invalid cylinder obscurer")
        trace_obscurers.append({
            "surface_id": next_surface_id,
            "first": first,
            "second": second,
            "diameter_m": diameter,
        })
        next_surface_id += 1
    for row in rows:
        if "tangent" in row:
            continue
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
    trace_reflectivity = []
    if include_response and not isinstance(reflectivity, list):
        raise OpticalModelCompileError("trace model found invalid primary reflectivity")
    if include_response and reflectivity:
        validated = _validated_mirror_response(reflectivity, "primary")
        if any("incidence_angle_deg" in entry for entry in validated):
            raise OpticalModelCompileError(
                "CSV export cannot represent incidence-dependent response"
            )
        trace_reflectivity = [(entry["wavelength_nm"], entry["response"]) for entry in validated]
    return rows, trace_obscurers, trace_reflectivity


def _trace_response(response: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Export applied values; retain measurement metadata in the imported optical model."""
    return [
        {
            key: entry[key]
            for key in ("wavelength_nm", "response", "incidence_angle_deg")
            if key in entry
        }
        for entry in response
    ]


def _camera_response_fields(optical_model: dict[str, Any]) -> dict[str, Any]:
    camera = optical_model.get("camera", {})
    if not any(key in camera for key in ("filter_response", "lightguide_response", "transmission")):
        return {}
    response = {
        # Identity adds no extra loss when no scalar response was supplied.
        "camera_transmission": camera.get("transmission", 1.0),
        "semantics": "measured_complete_response",
    }
    if "filter_response" in camera:
        response["camera_filter"] = _trace_response(camera["filter_response"])
    if "lightguide_response" in camera:
        response["lightguide_efficiency"] = _trace_response(camera["lightguide_response"])
    return {"camera_response": response}


def build_trace_model(optical_model: dict[str, Any]) -> dict[str, Any]:
    """Build the finite surfaces consumed directly from the compiled JSON.

    This is deliberately embedded in the canonical optical-model document:
    there is no second, flattened transport file to keep in sync.
    """
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
        return _axisymmetric_trace_model(optical_model)
    rows, obscurers, _ = trace_surface_rows(optical_model, include_response=False)
    response = primary.get("reflectivity", [])
    if not isinstance(response, list):
        raise OpticalModelCompileError("trace model found invalid primary reflectivity")
    reflectivity = _validated_mirror_response(response, "primary") if response else []
    facets = []
    detectors = []
    for row in rows:
        surface = {
            "id": row["surface_id"],
            "shape": row["shape"],
            "centre_m": row["centre_m"],
            "normal": row["normal"],
            "tangent": row["tangent"],
            "diameter_m": row["diameter_m"],
        }
        if row["role"] == "mirror":
            facets.append({**surface, "focal_length_m": row["focal_length_m"]})
        else:
            detectors.append(surface)
    return {
        "kind": "segmented",
        **_camera_response_fields(optical_model),
        "primary_facets": facets,
        "detector_surfaces": detectors,
        "cylinder_obscurers": [
            {
                "id": item["surface_id"],
                "first_endpoint_m": item["first"],
                "second_endpoint_m": item["second"],
                "diameter_m": item["diameter_m"],
            }
            for item in obscurers
        ],
        "primary_reflectivity": _trace_response(reflectivity),
    }


def native_surface_rows(
    optical_model: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[tuple[float, float]]]:
    """Return the planar surfaces used by the dependency-free native writer."""
    return trace_surface_rows(optical_model)


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
    provenance = optical_model.get("provenance", {})
    model = provenance.get("model")
    version = provenance.get("model_version")
    content_hash = optical_model.get("optical_model_sha256")
    if (
        not isinstance(model, str)
        or not isinstance(version, str)
        or "," in model
        or "," in version
        or not isinstance(content_hash, str)
        or len(content_hash) != 64
    ):
        raise OpticalModelCompileError("native export requires model provenance")
    lines = [
        "obdeect-optical-model-v1",
        f"provenance,{model},{version},{content_hash}",
        "surface_id,role,shape,cx_m,cy_m,cz_m,nx,ny,nz,tx,ty,tz,diameter_m,focal_length_m",
    ]
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
                    for value in (
                        *obscurer["first"],
                        *obscurer["second"],
                        obscurer["diameter_m"],
                    )
                ),
            ])
        )
    for wavelength, value in reflectivity:
        lines.append(f"primary_reflectivity,{wavelength:.17g},{value:.17g}")
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_native_axisymmetric_optical_model(optical_model: dict[str, Any], output: Path) -> None:
    """Write an exact rotationally symmetric SST/SCT optical prescription."""
    provenance = optical_model.get("provenance", {})
    model = provenance.get("model")
    version = provenance.get("model_version")
    content_hash = optical_model.get("optical_model_sha256")
    if (
        not isinstance(model, str)
        or not isinstance(version, str)
        or not isinstance(content_hash, str)
        or len(content_hash) != 64
        or "," in model
        or "," in version
    ):
        raise OpticalModelCompileError("native export requires valid model provenance")

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
        local = [float(value) for value in coefficients]
        vertex = local[0]
        local[0] = 0.0
        return ",".join([
            role,
            *(f"{value:.17g}" for value in (vertex, inner, outer, scale, *local)),
        ])

    lines = [
        "obdeect-axisymmetric-optical-model-v1",
        f"provenance,{model},{version},{content_hash}",
        "role,vertex_z_m,inner_radius_m,outer_radius_m,radial_scale_m,c0_m,c1_m,c2_m,c3_m,c4_m,c5_m,c6_m,c7_m,c8_m,c9_m,c10_m,c11_m,c12_m",
        row("primary", primary),
        row("secondary", secondary),
        row("detector", focal),
    ]
    for role, response_source in (("primary", primary_source), ("secondary", secondary)):
        response = response_source.get("reflectivity", [])
        if not isinstance(response, list):
            raise OpticalModelCompileError(f"invalid {role} reflectivity")
        response = _validated_mirror_response(response, role) if response else []
        if any("incidence_angle_deg" in entry for entry in response):
            raise OpticalModelCompileError(
                "CSV export cannot represent incidence-dependent response"
            )
        for entry in response:
            wavelength, value = entry["wavelength_nm"], entry["response"]
            lines.append(f"{role}_reflectivity,{wavelength:.17g},{value:.17g}")
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _axisymmetric_trace_model(optical_model: dict[str, Any]) -> dict[str, Any]:
    """Preserve exact rotationally symmetric surfaces inside compiled JSON."""
    primary = optical_model["primary"]["aspheric_surface"]
    secondary = optical_model["secondary"]
    focal = optical_model.get("focal_surface")
    if not isinstance(focal, dict) or focal.get("outer_radius_m") is None:
        raise OpticalModelCompileError("trace model requires a bounded focal surface")

    def surface(source: dict[str, Any], role: str) -> dict[str, Any]:
        coefficients = source.get("coefficient_m")
        outer = source.get("outer_radius_m")
        inner = source.get("inner_radius_m", 0.0)
        scale = source.get("radial_scale_m", 1.0)
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
        local = [float(value) for value in coefficients]
        vertex = local[0]
        local[0] = 0.0
        return {
            "vertex_z_m": vertex,
            "inner_radius_m": float(inner),
            "outer_radius_m": float(outer),
            "radial_scale_m": float(scale),
            "coefficient_m": local,
        }

    primary_segments = optical_model["primary"].get("segments", [])
    secondary_segments = optical_model["secondary"].get("segments", [])
    next_id = 0
    masks = {}
    for role, segments in (("primary", primary_segments), ("secondary", secondary_segments)):
        masks[role] = []
        for segment in segments:
            masks[role].append({**segment, "id": next_id})
            next_id += 1
    return {
        "primary_segments": masks["primary"],
        "secondary_segments": masks["secondary"],
        "primary_surface_id": next_id,
        "secondary_surface_id": next_id + 1,
        "detector_surface_id": next_id + 2,
        "detector_surfaces": [
            {"id": next_id + 3 + index, **plane}
            for index, plane in enumerate(
                optical_model.get("camera", {}).get("entrance_surfaces", [])
            )
        ],
        "block_incoming_secondary": True,
        "kind": "axisymmetric",
        **_camera_response_fields(optical_model),
        "primary": surface(primary, "primary"),
        "secondary": surface(secondary, "secondary"),
        "detector": surface(focal, "detector"),
        "primary_reflectivity": _trace_response(
            _validated_mirror_response(optical_model["primary"].get("reflectivity"), "primary")
        ),
        "secondary_reflectivity": _trace_response(
            _validated_mirror_response(secondary.get("reflectivity"), "secondary")
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compile a simulation-models production into an optical model."
    )
    parser.add_argument(
        "--source-root", type=Path, required=True, help="simulation-models repository root"
    )
    parser.add_argument("--model", required=True, help="production table, e.g. LSTN-design")
    parser.add_argument("--version", required=True, help="production model version")
    parser.add_argument("--output", type=Path, required=True, help="compiled optical model JSON")
    parser.add_argument(
        "--require-trace-ready",
        action="store_true",
        help="fail if unresolved optical geometry or response prevents production tracing",
    )
    args = parser.parse_args()
    try:
        ir = resolve_model(args.source_root, args.model, args.version)
        optical_model = compile_optical_model(ir, args.source_root)
        if args.require_trace_ready:
            require_trace_ready(optical_model)
    except (OSError, ModelImportError, OpticalModelCompileError) as error:
        raise SystemExit(f"optical model compilation failed: {error}") from error
    args.output.write_text(
        json.dumps(optical_model, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"Compiled {optical_model['provenance']['model']} geometry into {args.output}")


if __name__ == "__main__":
    main()
