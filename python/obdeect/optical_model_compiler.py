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
import random
from pathlib import Path
from typing import Any

from obdeect.camera_config import (
    PIXEL_APERTURE_SHAPES,
    CameraConfigError,
    parse_camera_layout,
    parse_camera_layout_ecsv,
    parse_camera_pixel_types,
)
from obdeect.camera_surfaces import compile_camera_surfaces
from obdeect.model_import import REFERENCE_DIAGNOSTIC_PARAMETERS, component, record, resolve_model
from obdeect.model_import import ImportError as ModelImportError
from obdeect.pixel_response_compiler import bind_pixel_responses, measured_pixel_response
from obdeect.simtel_tables import TableImportError, parse_rpol_table


class OpticalModelCompileError(ValueError):
    """The provenance IR cannot be compiled without guessing optical data."""


_SHAPES = {0: "circle", 1: "hexagon_flat_x", 2: "square", 3: "hexagon_flat_y"}
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
        token = fields[0].lower()
        aliases = {
            "hexsegments": "hex",
            "yhexsegments": "yhex",
            "ringsegments": "ring",
            "squaresegments": "square",
            "circularsegments": "circle",
            "polygonsegments": "polygon",
        }
        kind = next((value for name, value in aliases.items() if name.startswith(token)), token)
        if kind == "polygon":
            try:
                count = int(fields[1])
                angle = math.radians(float(fields[2]))
                values = [float(v) for v in " ".join(fields[3:]).replace(",", " ").split()]
            except (ValueError, IndexError) as error:
                raise OpticalModelCompileError("invalid polygon segment") from error
            if (
                count != 1
                or len(values) < 6
                or len(values) % 2
                or not all(math.isfinite(v) for v in values)
                or not math.isfinite(angle)
            ):
                raise OpticalModelCompileError("polygon segment requires one finite polygon")
            vertices = [
                [
                    0.01 * (x * math.cos(angle) + y * math.sin(angle)),
                    0.01 * (-x * math.sin(angle) + y * math.cos(angle)),
                ]
                for x, y in zip(values[::2], values[1::2])
            ]
            if vertices[0] == vertices[-1]:
                vertices.pop()
            turns = [
                (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
                for a, b, c in zip(
                    vertices, vertices[1:] + vertices[:1], vertices[2:] + vertices[:2]
                )
            ]
            if len(vertices) < 3 or not (all(v > 0 for v in turns) or all(v < 0 for v in turns)):
                raise OpticalModelCompileError("polygon segment must be strictly convex")
            segments.append({"id": len(segments), "shape": "polygon", "vertices_xy_m": vertices})
            continue
        if kind not in {"hex", "yhex", "ring", "square", "circle"}:
            raise OpticalModelCompileError(
                f"segmentation line {line_number}: unsupported type {fields[0]}"
            )
        expected_fields = {5, 6} if kind != "ring" else {5, 6, 7}
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
        if kind in {"hex", "yhex", "square", "circle"}:
            x_cm, y_cm, diameter_cm, rotation_deg = (*values, 0.0)[0:4]
            if count != 1 or diameter_cm <= 0.0:
                raise OpticalModelCompileError(
                    f"segmentation line {line_number}: invalid hex footprint"
                )
            segments.append({
                "id": len(segments),
                "shape": {"square": "square", "circle": "circle"}.get(kind, "hexagon"),
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
            {"x", "y", "diameter"}
            if kind in {"hex", "yhex", "square", "circle"}
            else {"r_min", "r_max", "dphi"}
        )
        optional = {"rotation"} if kind in {"hex", "yhex", "square", "circle"} else {"phi0", "gap"}
        count = group.get("count")
        if (
            kind not in {"hex", "yhex", "ring", "square", "circle"}
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
        if segment["shape"] == "polygon":
            transformed.append({
                **segment,
                "vertices_xy_m": [[-x, y] for x, y in segment["vertices_xy_m"]],
            })
        elif segment["shape"] in {"hexagon", "square", "circle"}:
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
    contents: str, *, fallback_focal_length_m: float | None, allow_automatic: bool = False
) -> list[dict[str, Any]]:
    """Parse the documented sim_telarray mirror-list columns, in centimetres.

    The sim_telarray reader consumes four required columns (x, y, diameter,
    focal length), an optional shape (default circle), and an optional z position.
    It ignores text after
    those fields, which model files commonly use for comments and panel IDs.
    This parser retains only the documented geometric values: trailing data is
    deliberately not interpreted as a normal, rotation, or alignment.

    Focal length zero requires an explicit catalogue fallback or automatic
    prescription resolution. Negative values request a random focal-length error.
    """
    if fallback_focal_length_m is not None and (
        not math.isfinite(fallback_focal_length_m) or fallback_focal_length_m <= 0.0
    ):
        raise OpticalModelCompileError("fallback focal length must be finite and positive")
    facets: list[dict[str, Any]] = []
    for line_number, source_line in enumerate(contents.splitlines(), start=1):
        line = source_line.strip()
        if not line or line.startswith(("#", "%")):
            continue
        fields = line.split("#", 1)[0].split("%", 1)[0].split()
        # Astropy ECSV keeps its column names in one un-commented header row.
        # The numerical columns that follow have the same documented layout as
        # sim_telarray mirror lists, so skip only this exact header form.
        if fields[:5] == ["mirror_x", "mirror_y", "mirror_diameter", "focal_length", "shape_type"]:
            continue
        if len(fields) < 4:
            raise OpticalModelCompileError(
                f"mirror list line {line_number}: expected at least 4 columns"
            )
        try:
            x_cm, y_cm, diameter_cm, focal_cm = (float(item) for item in fields[:4])
            shape_value = float(fields[4]) if len(fields) >= 5 else 0.0
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
        catalogue_focal_m = focal_cm * 0.01
        focal_m = abs(catalogue_focal_m)
        if focal_m == 0.0:
            if fallback_focal_length_m is None and not allow_automatic:
                raise OpticalModelCompileError(
                    f"mirror list line {line_number}: zero focal length has no fallback"
                )
            focal_m = fallback_focal_length_m or 0.0
        if focal_m < 0.0:
            raise OpticalModelCompileError(
                f"mirror list line {line_number}: focal length must be positive"
            )
        facet: dict[str, Any] = {
            "id": len(facets),
            "centre_m": [x_cm * 0.01, y_cm * 0.01, z_cm * 0.01],
            "shape": _SHAPES[shape_code],
            "diameter_m": diameter_cm * 0.01,
            "focal_length_m": focal_m,
            "catalogue_focal_length_m": catalogue_focal_m,
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


def parse_rectangular_response(
    contents: str, name: str, *, filename_options: str = ""
) -> dict[str, Any]:
    """Normalize RPOL formats 2/3 and all declared grid interpolation options."""
    try:
        table = parse_rpol_table(contents, default_dimension=3, filename_options=filename_options)
    except TableImportError as error:
        raise OpticalModelCompileError(f"{name}: {error}") from error
    if not table["y"]:
        raise OpticalModelCompileError(f"{name} requires a two-dimensional table")
    if any(not 0 <= value <= 1 for value in table["response"]):
        raise OpticalModelCompileError(f"{name} response must be a fraction")
    return {**table, "clip": table["interpolation"]["boundary"] == "zero"}


def parse_spatial_response(
    contents: str, name: str, *, filename_options: str = ""
) -> dict[str, Any]:
    table = parse_rectangular_response(contents, name, filename_options=filename_options)
    return {
        "x_m": [v * 0.01 for v in table["x"]],
        "y_m": [v * 0.01 for v in table["y"]],
        "response": table["response"],
        "clip": table["clip"],
        "interpolation": table["interpolation"],
    }


def parse_obscuration_quadrilaterals(contents: str) -> list[dict[str, Any]]:
    """Read the production ECSV vertices in metres, retaining convex finite plates."""
    header = None
    result = []
    names = [f"{axis}{index}" for index in range(1, 5) for axis in "xyz"]
    for source in contents.splitlines():
        line = source.split("#", 1)[0].strip()
        if not line:
            continue
        fields = line.split()
        if header is None:
            header = fields
            if not set(names).issubset(header):
                raise OpticalModelCompileError("quadrilateral obscurers lack vertex columns")
            continue
        if len(fields) != len(header):
            raise OpticalModelCompileError("quadrilateral obscurer has wrong column count")
        row = dict(zip(header, fields, strict=True))
        try:
            vertices = [[float(row[f"{axis}{i}"]) for axis in "xyz"] for i in range(1, 5)]
        except ValueError as error:
            raise OpticalModelCompileError(
                "quadrilateral obscurer has nonnumeric vertices"
            ) from error
        if any(not math.isfinite(v) for point in vertices for v in point):
            raise OpticalModelCompileError("quadrilateral obscurer has nonfinite vertices")
        edges = [
            [b - a for a, b in zip(vertices[i], vertices[(i + 1) % 4], strict=True)]
            for i in range(4)
        ]

        def cross(a, b):
            return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]

        normal = cross(edges[0], edges[1])
        magnitude = math.sqrt(sum(v * v for v in normal))
        if magnitude <= 1e-12:
            raise OpticalModelCompileError("quadrilateral obscurer is degenerate")
        normal = [v / magnitude for v in normal]
        if abs(sum(normal[j] * (vertices[3][j] - vertices[0][j]) for j in range(3))) > 1e-12:
            raise OpticalModelCompileError("quadrilateral obscurer must be planar")
        if any(
            sum(a * b for a, b in zip(cross(edges[i], edges[(i + 1) % 4]), normal, strict=True))
            <= 1e-12
            for i in range(4)
        ):
            raise OpticalModelCompileError("quadrilateral obscurer must be convex")
        result.append({"shape": "quadrilateral", "vertices_m": vertices})
    return result


def compile_secondary_baffle(
    parameter: dict[str, Any], overall_offset_m: float
) -> dict[str, Any] | None:
    """Manual 12.3.3: open baffle with optional radial thickness and end radius."""
    values, units = parameter.get("value"), parameter.get("unit")
    if not isinstance(values, list) or not 3 <= len(values) <= 5:
        raise OpticalModelCompileError("secondary baffle requires 3 to 5 lengths")
    if isinstance(units, str):
        units = [units] * len(values)
    if (
        not isinstance(units, list)
        or len(units) != len(values)
        or any(u not in _UNIT_TO_M for u in units)
    ):
        raise OpticalModelCompileError("secondary baffle requires length units")
    z1, z2, r1, thickness, r2 = (
        [_number(v, "secondary baffle") * _UNIT_TO_M[u] for v, u in zip(values, units)] + [0.0, 0.0]
    )[:5]
    if r1 == 0:
        return None
    r2 = r2 or r1
    if r1 < 0 or r2 < 0 or thickness < 0 or z1 == z2:
        raise OpticalModelCompileError("invalid secondary baffle")
    return {
        "shape": "hollow_frustum",
        "first_endpoint_m": [0.0, 0.0, z1 - overall_offset_m],
        "second_endpoint_m": [0.0, 0.0, z2 - overall_offset_m],
        "first_radius_m": r1,
        "second_radius_m": r2,
        "thickness_m": thickness,
    }


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
            or set(entry) - required - _RESPONSE_METADATA - {"interpolation"}
        ):
            raise OpticalModelCompileError(f"invalid {name} reflectivity columns")
        parsed = {
            key: _number(value, f"{name} {key}")
            for key, value in entry.items()
            if key != "interpolation"
        }
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
    if values and "interpolation" in values[0]:
        result[0]["interpolation"] = values[0]["interpolation"]
    return result


def parse_wavelength_response(
    contents: str, name: str, value_column: str = "reflectivity", *, filename_options: str = ""
) -> list[dict[str, float]]:
    """Read spectral response, preserving a complete incidence-angle grid when supplied."""
    first = next((line.strip() for line in contents.splitlines() if line.strip()), "")
    numeric = next(
        (
            line.strip()
            for line in contents.splitlines()
            if line.strip() and not line.lstrip().startswith(("#", "%"))
        ),
        "",
    )
    try:
        float(numeric.split()[0])
        plain = True
    except ValueError, IndexError:
        plain = False
    if first.startswith("#@RPOL@") or plain:
        try:
            table = parse_rpol_table(contents, filename_options=filename_options)
        except TableImportError as error:
            raise OpticalModelCompileError(f"{name}: {error}") from error
        rows = []
        for i, wavelength in enumerate(table["x"]):
            if table["y"]:
                rows.extend(
                    {
                        "wavelength_nm": wavelength,
                        "incidence_angle_deg": angle,
                        "response": table["response"][i * len(table["y"]) + j],
                    }
                    for j, angle in enumerate(table["y"])
                )
            else:
                rows.append({"wavelength_nm": wavelength, "response": table["response"][i]})
        rows[0]["interpolation"] = table["interpolation"]
        return _validated_mirror_response(rows, name)

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


def parse_angular_response(
    contents: str, name: str, *, filename_options: str = ""
) -> list[dict[str, float]]:
    """Read measured lightguide response without adding geometric losses twice."""
    numeric = next(
        (
            line.strip()
            for line in contents.splitlines()
            if line.strip() and not line.lstrip().startswith(("#", "%"))
        ),
        "",
    )
    try:
        float(numeric.split()[0])
        plain = True
    except ValueError, IndexError:
        plain = False
    if "#@RPOL@" in contents or plain:
        try:
            table = parse_rpol_table(contents, filename_options=filename_options)
        except TableImportError as error:
            raise OpticalModelCompileError(f"{name}: {error}") from error
        if (
            table["y"]
            or any(not 0 <= a <= 90 for a in table["x"])
            or any(not 0 <= v <= 1 for v in table["response"])
        ):
            raise OpticalModelCompileError(f"invalid {name} angular response")
        rows = [
            {"incidence_angle_deg": x, "response": y} for x, y in zip(table["x"], table["response"])
        ]
        rows[0]["interpolation"] = table["interpolation"]
        return rows

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
        distance = facet.get("resolved_focus_distance_m", math.hypot(focal_length - height, radius))
        if file_z > 0:
            z = file_z - offset
        else:
            z = focal_length - math.sqrt(distance * distance - radius * radius) - offset
        # A nominal facet normal is the bisector of the reverse incident
        # direction (+z for a star) and the direction from the facet centre
        # to the focal point.  Derive it from the final placement, including
        # mirror_offset, so every chief ray reaches the compiled focal plane.
        focal_distance = distance
        if not math.isfinite(focal_distance) or focal_distance <= 0.0:
            raise OpticalModelCompileError("facet has no finite focal-point direction")
        opt = parameters.get("mirror_opt", {}).get("value", [0, 0, 0])
        inclination = (
            0.5 * math.asin(radius / focal_distance) * (1 + opt[0] * radius**2 / focal_length**2)
        )
        phi = math.atan2(y, x)
        facet["nominal_centre_m"] = [x, y, z]
        facet["nominal_normal"] = [
            -math.cos(phi) * math.sin(inclination),
            -math.sin(phi) * math.sin(inclination),
            math.cos(inclination),
        ]
    apply_panel_alignment(facets, {}, 0, None)


def resolve_focus_offset(parameter: dict[str, Any], zenith_deg: float | None) -> float:
    """Manual section 12.3.3: linear cosine/sine camera focus deformation."""
    values, units = parameter["value"], parameter["unit"]
    base = _number(values[0], "focus_offset") * _UNIT_TO_M[units[0]]
    if not values[2] and not values[3]:
        return base
    if zenith_deg is None:
        raise OpticalModelCompileError("zenith-dependent focus requires alignment zenith")
    theta, reference = math.radians(_number(zenith_deg, "zenith")), math.radians(values[1])
    return (
        base
        + _number(values[2], "focus_offset")
        * _UNIT_TO_M[units[0] if units[2] == "null" else units[2]]
        * (math.cos(theta) - math.cos(reference))
        + _number(values[3], "focus_offset")
        * _UNIT_TO_M[units[0] if units[3] == "null" else units[3]]
        * (math.sin(theta) - math.sin(reference))
    )


def resolve_panel_prescription(
    facets: list[dict[str, Any]], parameters: dict[str, Any], seed: int
) -> None:
    """Manual 11.7/12.3.3, with the exact tel_setup_primary selection rules.

    A positive catalogue focal length is fixed. Zero is automatic; negative
    fixes the nominal magnitude but requests a manufacturing-error draw.
    Positive catalogue heights prohibit distance-error draws.
    """
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed < 2**64:
        raise OpticalModelCompileError("panel prescription requires a uint64 seed")
    flip = parameters.get("flip_mirrors", {}).get("value", False)
    if flip not in (False, True, 0, 1):
        raise OpticalModelCompileError("flip_mirrors must be boolean")
    if flip:
        for facet in facets:
            facet["centre_m"][:2] = facet["centre_m"][:2][::-1]
            facet["shape"] = {
                "hexagon_flat_x": "hexagon_flat_y",
                "hexagon_flat_y": "hexagon_flat_x",
            }.get(facet["shape"], facet["shape"])
    focal = _length_m(parameters["focal_length"], "focal_length")
    dish = _length_m(parameters["dish_shape_length"], "dish_shape_length", allow_zero=True) or focal
    radii = [math.hypot(*facet["centre_m"][:2]) for facet in facets]
    minimum, maximum = min(radii), max(radii)
    grading = (
        _signed_length_m(parameters["grading_of_focal_length"], "grading_of_focal_length")
        if "grading_of_focal_length" in parameters
        else 0.0
    )
    distance_sigma = (
        _length_m(
            parameters["mirror_align_random_distance"],
            "mirror_align_random_distance",
            allow_zero=True,
        )
        if "mirror_align_random_distance" in parameters
        else 0.0
    )
    error = parameters.get("random_focal_length", {"value": [0, 0], "unit": ["m", "m"]})
    values, units = error["value"], error["unit"]
    if (
        not isinstance(values, list)
        or len(values) != 2
        or not isinstance(units, list)
        or len(units) != 2
        or any(unit not in _UNIT_TO_M for unit in units)
    ):
        raise OpticalModelCompileError(
            "random_focal_length requires RMS and clipping width with length units"
        )
    sigma, limit = [
        _number(v, "random_focal_length") * _UNIT_TO_M[u] for v, u in zip(values, units)
    ]
    if sigma < 0 or limit < 0:
        raise OpticalModelCompileError("random_focal_length widths must be nonnegative")
    scales = parameters.get("mirror_f_scale", {}).get("value", [1] * len(facets))
    if not isinstance(scales, list):
        scales = [scales] * len(facets)
    if len(scales) != len(facets) or any(_number(v, "mirror_f_scale") <= 0 for v in scales):
        raise OpticalModelCompileError("mirror_f_scale must give one positive scale per facet")
    opt = parameters.get("mirror_opt", {}).get("value", [0, 0, 0])
    if not isinstance(opt, list) or len(opt) != 3:
        raise OpticalModelCompileError("mirror_opt requires three coefficients")
    opt = [_number(v, "mirror_opt") for v in opt]
    for facet, radius, scale in zip(facets, radii, scales):
        z = facet["centre_m"][2]
        if z:
            height = abs(z)
        elif parameters["parabolic_dish"]["value"]:
            height = radius**2 / (4 * dish)
        else:
            if radius > dish:
                raise OpticalModelCompileError("facet lies outside dish sphere")
            height = dish - math.sqrt(dish * dish - radius * radius)
        nominal_distance = math.hypot(focal - height, radius)
        distance = nominal_distance * (1 + opt[2] + opt[1] * radius * radius / focal**2)
        rng = random.Random(f"obdeect-panel-distance-v1:{seed}:{facet['id']}")
        if z <= 0 and distance_sigma:
            distance += rng.gauss(0, distance_sigma)
        if distance <= radius or not math.isfinite(distance):
            raise OpticalModelCompileError("resolved facet distance must exceed radial offset")
        facet["resolved_focus_distance_m"] = distance
        catalogue = facet.get("catalogue_focal_length_m", facet["focal_length_m"])
        if catalogue > 0:
            facet["focal_length_m"] = catalogue * scale
            continue
        nominal = (
            abs(catalogue) * scale if catalogue < 0 else facet["focal_length_m"] or nominal_distance
        )
        if catalogue == 0 and maximum > minimum:
            nominal += grading * (radius - 0.5 * (minimum + maximum)) / (maximum - minimum)
        rng = random.Random(f"obdeect-panel-focal-v1:{seed}:{facet['id']}")
        if limit > 100 * sigma:
            error_draw = limit * (rng.random() - 0.5)
        elif sigma:
            # Exact clipped-normal draw by inverse CDF avoids unbounded rejection.
            from statistics import NormalDist

            normal = NormalDist()
            lo = normal.cdf(-limit / sigma) if limit else 0.0
            hi = normal.cdf(limit / sigma) if limit else 1.0
            probability = max(
                math.nextafter(0.0, 1.0),
                min(math.nextafter(1.0, 0.0), lo + (hi - lo) * rng.random()),
            )
            error_draw = sigma * normal.inv_cdf(probability)
        else:
            error_draw = 0.0
        facet["focal_length_m"] = nominal + error_draw
        if facet["focal_length_m"] <= 0 or not math.isfinite(facet["focal_length_m"]):
            raise OpticalModelCompileError("resolved facet focal length must be positive")


def apply_panel_alignment(
    facets: list[dict[str, Any]], parameters: dict[str, Any], seed: int, zenith_deg: float | None
) -> None:
    """Resolve panel errors once; retain sim_telarray's Euler aperture frame.

    Per-panel random streams are keyed by catalogue ID, independent of ordering.
    Angular widths are normal-distribution standard deviations, not FWHM.
    """
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed < 2**64:
        raise OpticalModelCompileError("alignment requires a uint64 seed")
    widths = []
    for name in ("mirror_align_random_horizontal", "mirror_align_random_vertical"):
        record = parameters.get(name)
        if record is None:
            widths.append(0.0)
            continue
        values = record.get("value")
        if (
            not isinstance(values, list)
            or len(values) != 4
            or record.get("unit") not in (["deg", "deg", "null", "null"], ["deg"] * 4)
        ):
            raise OpticalModelCompileError(f"{name} requires four degree components")
        rms, reference, cosine, sine = [_number(v, name) for v in values]
        if rms < 0:
            raise OpticalModelCompileError(f"{name} RMS must be nonnegative")
        if (cosine or sine) and zenith_deg is None:
            raise OpticalModelCompileError(
                "zenith-dependent panel alignment requires alignment zenith"
            )
        theta = math.radians(
            _number(zenith_deg, "alignment zenith") if zenith_deg is not None else reference
        )
        reference = math.radians(reference)
        widths.append(
            math.radians(
                math.sqrt(
                    rms * rms
                    + cosine * cosine * (math.cos(theta) - math.cos(reference)) ** 2
                    + sine * sine * (math.sin(theta) - math.sin(reference)) ** 2
                )
            )
        )
    for facet in facets:
        normal = facet["nominal_normal"]
        x, y, _ = facet["nominal_centre_m"]
        phi = math.atan2(y, x)
        inclination = math.acos(max(-1.0, min(1.0, normal[2])))
        rng = random.Random(f"obdeect-panel-alignment-v1:{seed}:{facet.get('id', (x, y))}")
        horizontal, vertical = (rng.gauss(0.0, sigma) for sigma in widths)
        inclination += -horizontal * math.sin(phi) + vertical * math.cos(phi)
        gamma = horizontal * math.cos(phi) + vertical * math.sin(phi)
        c, t, ci, si, cg, sg = (
            math.cos(phi),
            math.sin(phi),
            math.cos(inclination),
            math.sin(inclination),
            math.cos(gamma),
            math.sin(gamma),
        )
        facet["nominal_normal"] = [-c * si * cg - t * sg, -t * si * cg + c * sg, ci * cg]
        u = [c * ci, t * ci, si]
        v = [-t * cg + c * si * sg, c * cg + t * si * sg, -ci * sg]
        facet["nominal_tangent"] = [c * a - t * b for a, b in zip(u, v)]
        facet["response_basis_u"] = u
        facet["response_basis_v"] = v


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


def compile_optical_model(
    ir: dict[str, Any],
    source_root: Path,
    *,
    scatter_seed: int | None = 0,
    alignment_seed: int = 0,
    alignment_zenith_deg: float | None = None,
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

    def table_options(name):
        value = parameters.get(name, {}).get("value")
        return value.partition("#rpol:")[2] if isinstance(value, str) else ""

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
    if (
        fallback is None
        and parameters.get("mirror_class", {}).get("value") != 0
        and isinstance(parameters.get("focal_length"), dict)
    ):
        fallback = _length_m(parameters["focal_length"], "focal_length")
    consumed = set()
    dual_surfaces = _dual_reflector_surfaces(parameters)
    nominal_geometry = False
    if "mirror_list" in verified_assets and dual_surfaces is None:
        facets = parse_simtel_mirror_list(
            verified_assets["mirror_list"].read_text(encoding="utf-8"),
            fallback_focal_length_m=fallback,
            allow_automatic=parameters.get("mirror_class", {}).get("value") == 0,
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
            resolve_panel_prescription(facets, parameters, alignment_seed)
            derive_nominal_single_reflector(facets, parameters)
            apply_panel_alignment(facets, parameters, alignment_seed, alignment_zenith_deg)
            consumed.update(
                name
                for name in (
                    "mirror_align_random_horizontal",
                    "mirror_align_random_vertical",
                    "mirror_align_random_distance",
                    "random_focal_length",
                    "grading_of_focal_length",
                    "mirror_f_scale",
                    "flip_mirrors",
                    "mirror_opt",
                )
                if name in parameters
            )
            primary["alignment"] = {
                "seed": alignment_seed,
                "zenith_angle_deg": alignment_zenith_deg,
                "method": "sim_telarray_euler_frame",
            }
            consumed.update(nominal_fields)
            consumed.add("mirror_class")
        evidence = {
            "source_format": "sim_telarray mirror-list (x, y, diameter, focal_length, shape[, z])",
            "facet_count": len(facets),
            "explicit_nonzero_z_count": sum(facet["centre_m"][2] != 0.0 for facet in facets),
            "normal_status": "nominal_unperturbed" if nominal_geometry else "unavailable",
            "in_plane_orientation_status": "sim_telarray_euler_frame"
            if nominal_geometry
            else "unavailable",
            "alignment_status": "compiled_seeded" if nominal_geometry else "unavailable",
            "interpretation": (
                "Nominal panel centres and normals follow sim_telarray tel_setup_primary; "
                "seeded panel alignment and exact aperture frames are compiled."
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
    if "telescope_obscuration_quadrilaterals" in verified_assets:
        primary["opaque_obscurers"] = parse_obscuration_quadrilaterals(
            verified_assets["telescope_obscuration_quadrilaterals"].read_text(encoding="utf-8")
        )
        consumed.add("telescope_obscuration_quadrilaterals")
    if "mirror_reflectivity" in verified_assets:
        try:
            primary["reflectivity"] = parse_wavelength_response(
                verified_assets["mirror_reflectivity"].read_text(encoding="utf-8"),
                "mirror_reflectivity",
                filename_options=table_options("mirror_reflectivity"),
            )
            consumed.add("mirror_reflectivity")
        except OpticalModelCompileError as error:
            if "incidence-dependent response" not in str(error):
                raise
    for name in ("primary_mirror_degraded_map", "primary_degraded_map"):
        if name in verified_assets:
            primary["degradation_map"] = parse_spatial_response(
                verified_assets[name].read_text(encoding="utf-8"),
                name,
                filename_options=table_options(name),
            )
            consumed.add(name)
        elif name in parameters and parameters[name].get("value") in (None, "none"):
            consumed.add(name)
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
        shadow_fields = {"secondary_mirror_shadow_diameter", "secondary_mirror_shadow_offset"}
        if shadow_fields.issubset(parameters):
            diameter = _signed_length_m(
                parameters["secondary_mirror_shadow_diameter"],
                "secondary_mirror_shadow_diameter",
            )
            if diameter < 0:
                diameter = 2 * secondary["outer_radius_m"]
            offset = _signed_length_m(
                parameters["secondary_mirror_shadow_offset"], "secondary_mirror_shadow_offset"
            )
            if offset <= 0:
                radius = diameter * 0.5 / secondary["radial_scale_m"]
                offset = secondary["coefficient_m"][0]
                # Only concave secondaries use their edge for the shadow plane;
                # telescope-frame coefficients have the opposite local sag sign.
                if secondary["coefficient_m"][1] < 0:
                    offset = sum(
                        value * radius ** (2 * index)
                        for index, value in enumerate(secondary["coefficient_m"])
                    )
            elif "mirror_offset" in parameters:
                offset -= _signed_length_m(parameters["mirror_offset"], "mirror_offset")
            secondary["incoming_shadow"] = {"diameter_m": diameter, "z_m": offset}
            consumed.update(shadow_fields)
        if "secondary_mirror_reflectivity" in verified_assets:
            try:
                secondary["reflectivity"] = parse_wavelength_response(
                    verified_assets["secondary_mirror_reflectivity"].read_text(encoding="utf-8"),
                    "secondary_mirror_reflectivity",
                    filename_options=table_options("secondary_mirror_reflectivity"),
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
    if "mirror_reflection_random_angle" in parameters:
        if scatter_seed is None:
            scatter_seed = 0
        parameter = parameters["mirror_reflection_random_angle"]
        values, units = parameter.get("value"), parameter.get("unit")
        if (
            not isinstance(values, list)
            or len(values) != 3
            or units != ["deg", "null", "deg"]
            or isinstance(scatter_seed, bool)
            or not isinstance(scatter_seed, int)
            or not 0 <= scatter_seed < 2**64
        ):
            raise OpticalModelCompileError(
                "mirror scatter requires three components and uint64 seed"
            )
        first, fraction, second = [
            _number(value, "mirror_reflection_random_angle") for value in values
        ]
        if not 0 <= first < 90 or not 0 <= second < 90 or not 0 <= fraction < 1:
            raise OpticalModelCompileError("invalid mirror scatter widths or mixture fraction")
        scatter = {
            "sigma1_rad": math.radians(first),
            "fraction2": fraction,
            "sigma2_rad": math.radians(second),
            "method": "surface_slopes" if dual_surfaces is not None else "outgoing_angles",
            "seed": scatter_seed,
        }
        primary["scatter"] = scatter
        if dual_surfaces is not None:
            secondary["scatter"] = dict(scatter)
        consumed.add("mirror_reflection_random_angle")
    if dual_surfaces is not None and "secondary_mirror_baffle" in parameters:
        overall = (
            _signed_length_m(parameters["mirror_offset"], "mirror_offset")
            if "mirror_offset" in parameters
            else 0.0
        )
        baffle = compile_secondary_baffle(parameters["secondary_mirror_baffle"], overall)
        if baffle is not None:
            primary.setdefault("opaque_obscurers", []).append(baffle)
        consumed.add("secondary_mirror_baffle")
    for name in ("secondary_mirror_degraded_map", "secondary_degraded_map"):
        if name in verified_assets and secondary is not None:
            secondary["degradation_map"] = parse_spatial_response(
                verified_assets[name].read_text(encoding="utf-8"),
                name,
                filename_options=table_options(name),
            )
            consumed.add(name)
        elif name in parameters and parameters[name].get("value") in (None, "none"):
            consumed.add(name)
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
        if "camera_scale_factor" in parameters:
            scale = _number(parameters["camera_scale_factor"].get("value"), "camera_scale_factor")
            if scale <= 0:
                raise OpticalModelCompileError("camera_scale_factor must be positive")
            for pixel in camera["pixels"]:
                pixel["centre_xy_m"] = [value * scale for value in pixel["centre_xy_m"]]
                if "z_offset_m" in pixel:
                    pixel["z_offset_m"] *= scale
            for pixel_type in camera["pixel_types"]:
                for key in ("cathode_diameter_m", "funnel_diameter_m", "funnel_depth_m"):
                    pixel_type[key] *= scale
            camera["scale_factor"] = scale
            consumed.add("camera_scale_factor")
        if "camera_rotate" in parameters:
            rotation = parameters["camera_rotate"]
            if rotation.get("unit") != "deg":
                raise OpticalModelCompileError("camera_rotate must be in degrees")
            camera["rotation_deg"] = _number(rotation.get("value"), "camera_rotate")
            consumed.add("camera_rotate")
        if "camera_pixel_types" in parameters:
            layout_types = {entry["id"]: entry for entry in camera["pixel_types"]}
            try:
                camera["pixel_types"] = parse_camera_pixel_types(parameters["camera_pixel_types"])
            except CameraConfigError as error:
                raise OpticalModelCompileError(str(error)) from error
            for pixel_type in camera["pixel_types"]:
                pixel_type["response_files"] = list(
                    layout_types.get(pixel_type["id"], {}).get("response_files", [])
                )
            type_ids = {entry["id"] for entry in camera["pixel_types"]}
            if any(pixel["type_id"] not in type_ids for pixel in camera["pixels"]):
                raise OpticalModelCompileError(
                    "camera layout references undefined physical pixel type"
                )
            consumed.add("camera_pixel_types")
            for pixel_type in camera["pixel_types"]:
                for key in ("cathode_diameter_m", "funnel_diameter_m", "funnel_depth_m"):
                    pixel_type[key] *= camera.get("scale_factor", 1.0)
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
        if "camera_degraded_map" in verified_assets:
            camera["degradation_map"] = parse_spatial_response(
                verified_assets["camera_degraded_map"].read_text(encoding="utf-8"),
                "camera_degraded_map",
                filename_options=table_options("camera_degraded_map"),
            )
            consumed.add("camera_degraded_map")
        elif "camera_degraded_map" in parameters and parameters["camera_degraded_map"].get(
            "value"
        ) in (None, "none"):
            consumed.add("camera_degraded_map")
        housing_fields = {"camera_body_diameter", "camera_body_shape", "camera_depth"}
        if (dual_surfaces is not None or nominal_geometry) and housing_fields.issubset(parameters):
            shape = parameters["camera_body_shape"].get("value")
            if (
                isinstance(shape, bool)
                or not isinstance(shape, int)
                or shape not in PIXEL_APERTURE_SHAPES
            ):
                raise OpticalModelCompileError("unsupported camera housing footprint")
            camera["housing"] = {
                "shape": _SHAPES[shape],
                "diameter_m": _length_m(
                    parameters["camera_body_diameter"], "camera_body_diameter", allow_zero=True
                ),
                "depth_m": _length_m(parameters["camera_depth"], "camera_depth", allow_zero=True),
                "front_z_m": (
                    dual_surfaces["focal_surface"]["coefficient_m"][0]
                    if dual_surfaces is not None
                    else _length_m(parameters["focal_length"], "focal_length")
                    - _signed_length_m(parameters["mirror_offset"], "mirror_offset")
                ),
                "sidewall_convention": "sim_telarray_cylindrical_sidewalls",
            }
            consumed.update(housing_fields)
        if "camera_filter" in verified_assets:
            camera["filter_response"] = parse_wavelength_response(
                verified_assets["camera_filter"].read_text(encoding="utf-8"),
                "camera_filter",
                "transmission",
                filename_options=table_options("camera_filter"),
            )
            consumed.add("camera_filter")
        if "lightguide_efficiency_vs_incidence_angle" in verified_assets:
            camera["lightguide_response"] = parse_angular_response(
                verified_assets["lightguide_efficiency_vs_incidence_angle"].read_text(
                    encoding="utf-8"
                ),
                "lightguide_efficiency_vs_incidence_angle",
                filename_options=table_options("lightguide_efficiency_vs_incidence_angle"),
            )
            consumed.add("lightguide_efficiency_vs_incidence_angle")
        if "camera_transmission" in parameters:
            fraction = _number(parameters["camera_transmission"]["value"], "camera_transmission")
            if not 0 <= fraction <= 1:
                raise OpticalModelCompileError("camera_transmission must be a response fraction")
            camera["transmission"] = fraction
            consumed.add("camera_transmission")
    if camera is not None and "camera_degraded_efficiency" in parameters:
        factor = _number(
            parameters["camera_degraded_efficiency"].get("value"), "camera_degraded_efficiency"
        )
        if not 0 <= factor <= 1:
            raise OpticalModelCompileError("camera_degraded_efficiency must be a fraction")
        camera["transmission"] = camera.get("transmission", 1.0) * factor
        consumed.add("camera_degraded_efficiency")
    if camera is not None and camera.get("pixel_types"):
        tables = []
        for pixel_type in camera["pixel_types"]:
            fields = pixel_type.get("source_fields", {})
            references = pixel_type.get("response_files", [])
            angle_name = fields.get("lightguide_angle_parameter")
            wavelength_name = fields.get("lightguide_wavelength_parameter")
            angular_path = verified_assets.get(angle_name) if angle_name else None
            wavelength_path = verified_assets.get(wavelength_name) if wavelength_name else None

            def response_file(filename):
                base = filename.partition("#rpol:")[0]
                if not component(base):
                    raise OpticalModelCompileError("unsafe pixel response filename")
                known = {path.resolve() for path in verified_assets.values() if path.name == base}
                if not known:
                    candidates = (
                        root / "model_parameters" / "Files" / base,
                        verified_assets[camera_asset_name].parent / base,
                    )
                    known = {
                        path.resolve()
                        for path in candidates
                        if path.is_file() and path.resolve().is_relative_to(root)
                    }
                if len(known) != 1:
                    raise OpticalModelCompileError(
                        f"pixel response file is missing or ambiguous: {filename}"
                    )
                path = known.pop()
                nested_assets[filename] = record(path, root)
                return path

            if references:
                angular_path = response_file(references[0])
                wavelength_path = response_file(references[1]) if len(references) > 1 else None
            if angular_path is None and angle_name is not None:
                raise OpticalModelCompileError(
                    f"pixel type references an unresolved angular response: {angle_name}"
                )
            if wavelength_path is None and wavelength_name is not None:
                raise OpticalModelCompileError(
                    f"pixel type references an unresolved spectral response: {wavelength_name}"
                )
            if angular_path is not None:
                try:
                    tables.append(
                        measured_pixel_response(
                            pixel_type["id"],
                            angular_path.read_text(encoding="utf-8"),
                            wavelength_path.read_text(encoding="utf-8")
                            if wavelength_path is not None
                            else None,
                        )
                    )
                except TableImportError as error:
                    raise OpticalModelCompileError(str(error)) from error
                if angle_name:
                    consumed.add(angle_name)
                if wavelength_name:
                    consumed.add(wavelength_name)
            else:
                transparency = fields.get(
                    "funnel_transparency", pixel_type.get("funnel_transparency")
                )
                reflection = fields.get(
                    "funnel_wall_reflectivity", pixel_type.get("funnel_wall_reflectivity")
                )
                if transparency is not None and reflection is not None:
                    transparency = _number(transparency, "funnel transparency")
                    reflection = _number(reflection, "funnel wall reflectivity")
                    if not 0 <= transparency <= 1 or not 0 <= reflection <= 1:
                        raise OpticalModelCompileError("numeric funnel response must be fractional")
                    tables.append({
                        "id": pixel_type["id"],
                        "method": "single_reflection",
                        "transparency": transparency,
                        "wall_reflectivity": reflection,
                    })
        if tables:
            unresolved_references = [
                name for name in unresolved_references if name not in nested_assets
            ]
            if len(tables) != len(camera["pixel_types"]):
                raise OpticalModelCompileError(
                    "every physical pixel type needs its optical response"
                )
            camera["pixel_optical_response_tables"] = tables
            camera.pop("lightguide_response", None)  # Already applied by the actual pixel type.
    diagnostics = set(parameters) & REFERENCE_DIAGNOSTIC_PARAMETERS
    deferred = sorted(set(parameters) - consumed - diagnostics)
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
        "reference_diagnostics": {
            name: {
                "reason": "derived reference data; not an individual-ray optical prescription",
                "parameter_record": records[f"parameter:{name}"],
                "asset": assets.get(name),
            }
            for name in sorted(diagnostics)
        },
        "field_coverage": {
            name: {
                "disposition": "reference_diagnostic"
                if name in diagnostics
                else "consumed"
                if name in consumed
                else "deferred",
                "reason": "parsed into geometry"
                if name in consumed
                else "derived reference data retained for comparison"
                if name in diagnostics
                else "not applied by the optical tracer",
            }
            for name in sorted(parameters)
        },
        "trace_blockers": [
            *(
                ["mirror scatter requires an explicit --scatter-seed"]
                if "mirror_reflection_random_angle" in parameters
                and "mirror_reflection_random_angle" not in consumed
                else []
            ),
            *([] if nominal_geometry else ["facet surface normals and alignment are not compiled"]),
            "physical detector surfaces and materials remain deferred",
            *(
                ["primary incidence-dependent reflectivity remains deferred"]
                if "mirror_reflectivity" in parameters and "mirror_reflectivity" not in consumed
                else []
            ),
            *(
                ["quadrilateral obscurers remain deferred"]
                if "telescope_obscuration_quadrilaterals" in parameters
                and "telescope_obscuration_quadrilaterals" not in consumed
                else []
            ),
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
    if (
        camera is not None
        and camera.get("pixel_types")
        and ("camera_pixel_types" in consumed or "camera_config_file" in consumed)
    ):
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
                or units[1] != "deg"
                or any(unit not in _UNIT_TO_M and unit != "null" for unit in units[2:])
            ):
                raise OpticalModelCompileError(
                    "camera placement requires a resolved static focus offset"
                )
            offset = resolve_focus_offset(focus, alignment_zenith_deg)
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
            report["deferred"] = sorted(set(parameters) - consumed - diagnostics)
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
        trace = compiled["trace_model"]
        identifiers = [
            entry["id"]
            for field, value in trace.items()
            if isinstance(value, list)
            for entry in value
            if isinstance(entry, dict) and "id" in entry
        ]
        identifiers += [
            trace[name]
            for name in ("primary_surface_id", "secondary_surface_id", "detector_surface_id")
            if name in trace
        ]
        next_identifier = max(identifiers, default=-1) + 1
        trace["opaque_obscurers"] = [
            dict(surface, id=next_identifier + index)
            for index, surface in enumerate(primary.get("opaque_obscurers", []))
        ]
        if trace["kind"] == "axisymmetric":
            for cylinder in primary.get("cylinder_obscurers", []):
                radius = cylinder["diameter_m"] * 0.5
                trace["opaque_obscurers"].append({
                    "id": next_identifier + len(trace["opaque_obscurers"]),
                    "shape": "solid_frustum",
                    "first_endpoint_m": cylinder["first_endpoint_m"],
                    "second_endpoint_m": cylinder["second_endpoint_m"],
                    "first_radius_m": radius,
                    "second_radius_m": radius,
                    "thickness_m": 0.0,
                })

        for role, source in (("primary", primary), ("secondary", secondary), ("camera", camera)):
            if source is not None and "degradation_map" in source:
                trace[f"{role}_degradation"] = source["degradation_map"]
                if role == "primary" and trace["kind"] == "segmented":
                    trace["primary_degradation_in_facet_frame"] = True
                if role == "camera":
                    angle = math.radians(camera.get("rotation_deg", 0))
                    sign = -1 if trace["kind"] == "segmented" else 1
                    trace["camera_degradation"] = {
                        **source["degradation_map"],
                        "x_basis": [sign * math.cos(angle), -sign * math.sin(angle), 0.0],
                        "y_basis": [math.sin(angle), math.cos(angle), 0.0],
                    }
                if role == "secondary":
                    trace[f"{role}_degradation"] = {
                        **source["degradation_map"],
                        "x_basis": [-1.0, 0.0, 0.0],
                        "y_basis": [0.0, 1.0, 0.0],
                    }

        if (
            camera is not None
            and camera.get("pixel_optical_response_tables")
            and camera.get("cathode_surfaces")
        ):
            try:
                pixel_response = bind_pixel_responses(camera, trace)
            except TableImportError as error:
                raise OpticalModelCompileError(str(error)) from error
            if pixel_response is not None:
                trace["pixel_responses"] = pixel_response
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
            ("obscurer", "incoming_obscurer_planes"),
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
        for field in ("primary_to_secondary_planes", "incoming_obscurer_planes"):
            components.extend(
                {"id": item["id"], "role": "obscurer", "source": f"trace_model.{field}"}
                for item in trace_model.get(field, [])
            )
        components.extend(
            {
                "id": item["id"],
                "role": "opaque_cylinder",
                "source": "trace_model.primary_to_secondary_cylinders",
            }
            for item in trace_model.get("primary_to_secondary_cylinders", [])
        )
    components.extend(
        {"id": item["id"], "role": "opaque_surface", "source": "trace_model.opaque_obscurers"}
        for item in trace_model.get("opaque_obscurers", [])
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
        if facet["shape"] != "circle" and "nominal_tangent" not in facet:
            raise OpticalModelCompileError("polygonal facets require an explicit aperture frame")
        rows.append({
            "surface_id": int(facet["id"]),
            "role": "mirror",
            "shape": facet["shape"],
            "centre_m": centre,
            "normal": normal,
            **({"tangent": facet["nominal_tangent"]} if "nominal_tangent" in facet else {}),
            "diameter_m": facet["diameter_m"],
            "focal_length_m": facet["focal_length_m"],
            **(
                {
                    "response_basis_u": facet["response_basis_u"],
                    "response_basis_v": facet["response_basis_v"],
                }
                if "response_basis_u" in facet
                else {}
            ),
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
    """Retain measured values and explicit simtel interpolation semantics."""
    rows = [
        {
            key: entry[key]
            for key in ("wavelength_nm", "response", "incidence_angle_deg")
            if key in entry
        }
        for entry in response
    ]
    if rows:
        rows[0]["interpolation"] = response[0].get(
            "interpolation",
            {
                "boundary": "clamp",
                "scheme": 1,
                "x_log": False,
                "y_log": False,
                "value_log": False,
                "coefficients": [],
            },
        )
    return rows


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
            facets.append({
                **surface,
                "focal_length_m": row["focal_length_m"],
                **(
                    {
                        "response_basis_u": row["response_basis_u"],
                        "response_basis_v": row["response_basis_v"],
                    }
                    if "response_basis_u" in row
                    else {}
                ),
            })
        else:
            detectors.append({
                **surface,
                **{
                    key: row[key]
                    for key in ("source_pixel_id", "source_type_id", "enabled")
                    if key in row
                },
            })
    incoming = []
    housing = optical_model.get("camera", {}).get("housing", {})
    if housing.get("diameter_m", 0) > 0:
        identifier = (
            max(
                [surface["id"] for surface in facets + detectors]
                + [item["surface_id"] for item in obscurers]
            )
            + 1
        )
        front, depth = housing["front_z_m"], housing["depth_m"]
        for z in [front] if depth == 0 else [front, front + depth]:
            incoming.append({
                "id": identifier,
                "shape": housing["shape"],
                "centre_m": [0.0, 0.0, z],
                "normal": [0.0, 0.0, 1.0],
                "tangent": [1.0, 0.0, 0.0],
                "diameter_m": housing["diameter_m"],
            })
            identifier += 1
    return {
        "kind": "segmented",
        "incoming_obscurer_planes": incoming,
        **({"primary_scatter": primary["scatter"]} if "scatter" in primary else {}),
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

    if (
        primary.get("scatter")
        or primary.get("opaque_obscurers")
        or optical_model.get("camera", {}).get("housing")
    ):
        raise OpticalModelCompileError(
            "surface-table export cannot represent mirror scatter or camera housing; use JSON"
        )
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
    if (
        primary_source.get("scatter")
        or optical_model["secondary"].get("scatter")
        or optical_model["secondary"].get("incoming_shadow")
        or optical_model.get("camera", {}).get("housing")
    ):
        raise OpticalModelCompileError(
            "CSV export cannot represent scatter or flight-specific obscurers; use compiled JSON"
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
    camera = optical_model.get("camera", {})
    housing = camera.get("housing", {})
    obscurer_planes, obscurer_cylinders, incoming_planes = [], [], []
    identifier = next_id + 3 + len(camera.get("entrance_surfaces", []))
    if housing.get("diameter_m", 0) > 0:
        front = housing["front_z_m"]
        back = front - housing["depth_m"]
        obscurer_planes.append({
            "id": identifier,
            "shape": housing["shape"],
            "centre_m": [0.0, 0.0, back],
            "normal": [0.0, 0.0, 1.0],
            "tangent": [1.0, 0.0, 0.0],
            "diameter_m": housing["diameter_m"],
        })
        if housing["depth_m"] > 0:
            # sim_telarray's usual build checks the rear face and cylindrical
            # sidewall. The front-face check is behind the disabled
            # RAYTRACING_CAMERA_FRONT_AFTER_PRIMARY build flag.
            obscurer_cylinders.append({
                "id": identifier + 1,
                "first_endpoint_m": [0.0, 0.0, back],
                "second_endpoint_m": [0.0, 0.0, front],
                "diameter_m": housing["diameter_m"],
            })
    shadow = secondary.get("incoming_shadow", {})
    if shadow.get("diameter_m", 0) > 0:
        incoming_planes.append({
            "id": identifier + len(obscurer_planes) + len(obscurer_cylinders),
            "shape": "circle",
            "centre_m": [0.0, 0.0, shadow["z_m"]],
            "normal": [0.0, 0.0, 1.0],
            "tangent": [1.0, 0.0, 0.0],
            "diameter_m": shadow["diameter_m"],
        })
    return {
        "primary_segments": masks["primary"],
        "secondary_segments": masks["secondary"],
        "primary_to_secondary_planes": obscurer_planes,
        "primary_to_secondary_cylinders": obscurer_cylinders,
        "incoming_obscurer_planes": incoming_planes,
        "primary_surface_id": next_id,
        "secondary_surface_id": next_id + 1,
        "detector_surface_id": next_id + 2,
        "detector_surfaces": [
            {"id": next_id + 3 + index, **plane}
            for index, plane in enumerate(
                optical_model.get("camera", {}).get("entrance_surfaces", [])
            )
        ],
        "block_incoming_secondary": "incoming_shadow" not in secondary,
        "kind": "axisymmetric",
        **(
            {"primary_scatter": optical_model["primary"]["scatter"]}
            if "scatter" in optical_model["primary"]
            else {}
        ),
        **({"secondary_scatter": secondary["scatter"]} if "scatter" in secondary else {}),
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
        "--scatter-seed",
        type=int,
        default=0,
        help="explicit uint64 seed for measured mirror scatter",
    )
    parser.add_argument("--alignment-seed", type=int, default=0)
    parser.add_argument("--alignment-zenith-deg", type=float)
    parser.add_argument(
        "--require-trace-ready",
        action="store_true",
        help="fail if unresolved optical geometry or response prevents production tracing",
    )
    args = parser.parse_args()
    try:
        ir = resolve_model(args.source_root, args.model, args.version)
        optical_model = compile_optical_model(
            ir,
            args.source_root,
            scatter_seed=args.scatter_seed,
            alignment_seed=args.alignment_seed,
            alignment_zenith_deg=args.alignment_zenith_deg,
        )
        if args.require_trace_ready:
            require_trace_ready(optical_model)
    except (OSError, ModelImportError, OpticalModelCompileError) as error:
        raise SystemExit(f"optical model compilation failed: {error}") from error
    args.output.write_text(
        json.dumps(optical_model, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"Compiled {optical_model['provenance']['model']} optical model into {args.output}")


if __name__ == "__main__":
    main()
