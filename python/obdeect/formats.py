"""Readers for simulation-models parameter records and sim_telarray data files.

Every reader converts to SI units and rejects input it cannot represent
exactly; none of them places optics in a telescope (see ``compiler``).
"""

from __future__ import annotations

import math
import re
from typing import Any

from obdeect.simtel_tables import TableImportError, parse_rpol_table


class ModelError(ValueError):
    """Model data cannot be converted without guessing."""


_UNIT_TO_M = {"m": 1.0, "cm": 0.01, "mm": 0.001}
# sim_telarray mirror and camera-body shape codes.
SHAPES = {0: "circle", 1: "hexagon_flat_x", 2: "square", 3: "hexagon_flat_y"}


def number(value: object, name: str) -> float:
    """Return a finite number or raise."""
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ModelError(f"{name} must be a finite number")
    return float(value)


def length_unit_m(unit: object, name: str) -> float:
    """Metres per length `unit`."""
    if unit not in _UNIT_TO_M:
        raise ModelError(f"{name} has unsupported length unit {unit!r}")
    return _UNIT_TO_M[unit]


def signed_length_m(record: dict[str, Any], name: str) -> float:
    """Return a parameter record's length value in metres."""
    return number(record.get("value"), name) * length_unit_m(record.get("unit"), name)


def length_m(record: dict[str, Any], name: str, *, allow_zero: bool = False) -> float:
    """Return a positive (or non-negative) length in metres."""
    value = signed_length_m(record, name)
    if value < 0.0 or (value == 0.0 and not allow_zero):
        raise ModelError(f"{name} must be {'non-negative' if allow_zero else 'positive'}")
    return value


def components(record: dict[str, Any], name: str, count: int | None = None) -> tuple[list, list]:
    """Values and per-value units of a list parameter; a single unit applies to all values."""
    values, units = record.get("value"), record.get("unit")
    if not isinstance(values, list) or not values or count not in (None, len(values)):
        raise ModelError(f"{name} requires {count or 'a list of'} components")
    if units is None or isinstance(units, str):
        units = [units] * len(values)
    if not isinstance(units, list) or len(units) != len(values):
        raise ModelError(f"{name}: units do not match its values")
    return values, units


def lengths_m(record: dict[str, Any], name: str, count: int) -> list[float]:
    """Return a list of `count` lengths in metres."""
    values, units = components(record, name, count)
    return [number(v, name) * length_unit_m(u, name) for v, u in zip(values, units, strict=True)]


def _data_lines(contents: str, comments: tuple[str, ...] = ("#", "%")) -> list[list[str]]:
    rows = []
    for line in contents.splitlines():
        for marker in comments:
            line = line.split(marker, 1)[0]
        if line.strip():
            rows.append(line.split())
    return rows


def read_mirror_list(contents: str) -> list[dict[str, Any]]:
    """Read a sim_telarray mirror list: x, y, diameter, focal length [cm], shape, [z].

    Text after the documented columns is ignored, as sim_telarray does. A focal
    length of zero asks for the automatic prescription; a negative one fixes its
    magnitude but asks for a random manufacturing error.
    """
    facets = []
    for fields in _data_lines(contents):
        if fields[:4] == ["mirror_x", "mirror_y", "mirror_diameter", "focal_length"]:
            continue  # ECSV column header
        if len(fields) < 4:
            raise ModelError(f"mirror list row {len(facets)}: expected at least 4 columns")
        try:
            x, y, diameter, focal = (float(v) for v in fields[:4])
            shape = float(fields[4]) if len(fields) > 4 else 0.0
            z = float(fields[5]) if len(fields) > 5 else 0.0
        except ValueError as error:
            raise ModelError(f"mirror list row {len(facets)}: invalid number") from error
        if not all(math.isfinite(v) for v in (x, y, diameter, focal, shape, z)):
            raise ModelError(f"mirror list row {len(facets)}: non-finite value")
        if shape not in SHAPES or diameter <= 0.0:
            raise ModelError(f"mirror list row {len(facets)}: invalid shape or diameter")
        facets.append({
            "id": len(facets),
            "centre_m": [0.01 * x, 0.01 * y, 0.01 * z],
            "shape": SHAPES[int(shape)],
            "diameter_m": 0.01 * diameter,
            "catalogue_focal_length_m": 0.01 * focal,
        })
    if not facets:
        raise ModelError("mirror list contains no mirrors")
    return facets


def read_segmentation(contents: str) -> list[dict[str, Any]]:
    """Read sim_telarray segmentation records (hex, yhex, square, circle, ring, polygon)."""
    names = {
        "hexsegments": "hex",
        "yhexsegments": "yhex",
        "ringsegments": "ring",
        "squaresegments": "square",
        "circularsegments": "circle",
        "polygonsegments": "polygon",
    }
    segments: list[dict[str, Any]] = []
    for row, fields in enumerate(_data_lines(contents, ("#",)), start=1):
        token = fields[0].lower()
        kind = next((value for name, value in names.items() if name.startswith(token)), token)
        try:
            count = int(fields[1])
            values = [float(v) for v in " ".join(fields[2:]).replace(",", " ").split()]
        except (ValueError, IndexError) as error:
            raise ModelError(f"segmentation row {row}: invalid number") from error
        if count < 1 or not all(math.isfinite(v) for v in values):
            raise ModelError(f"segmentation row {row}: invalid count or value")
        if kind == "polygon":
            if count != 1 or len(values) < 7 or len(values) % 2 == 0:
                raise ModelError(f"segmentation row {row}: polygon needs an angle and points")
            angle = math.radians(values[0])
            vertices = [
                [
                    0.01 * (x * math.cos(angle) + y * math.sin(angle)),
                    0.01 * (-x * math.sin(angle) + y * math.cos(angle)),
                ]
                for x, y in zip(values[1::2], values[2::2])
            ]
            if vertices[0] == vertices[-1]:
                vertices.pop()
            turns = [
                (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
                for a, b, c in zip(
                    vertices, vertices[1:] + vertices[:1], vertices[2:] + vertices[:2]
                )
            ]
            if len(vertices) < 3 or not (all(t > 0 for t in turns) or all(t < 0 for t in turns)):
                raise ModelError(f"segmentation row {row}: polygon must be strictly convex")
            segments.append({"id": len(segments), "shape": "polygon", "vertices_xy_m": vertices})
        elif kind in ("hex", "yhex", "square", "circle"):
            if count != 1 or len(values) not in (3, 4) or values[2] <= 0:
                raise ModelError(f"segmentation row {row}: invalid footprint")
            x, y, diameter, rotation = (*values, 0.0)[:4]
            segments.append({
                "id": len(segments),
                "shape": {"square": "square", "circle": "circle"}.get(kind, "hexagon"),
                "centre_xy_m": [0.01 * x, 0.01 * y],
                "diameter_m": 0.01 * diameter,
                "rotation_deg": rotation + (90.0 if kind == "yhex" else 0.0),
            })
        elif kind == "ring":
            if len(values) not in (3, 4, 5):
                raise ModelError(f"segmentation row {row}: invalid ring")
            inner, outer, span, start, gap = (*values, 0.0, 0.0)[:5]
            if inner < 0 or outer <= inner or span <= 0 or gap < 0 or count * span > 360 + 1e-9:
                raise ModelError(f"segmentation row {row}: invalid ring")
            segments.extend(
                {
                    "id": len(segments) + i,
                    "shape": "annular_sector",
                    "inner_radius_m": 0.01 * inner,
                    "outer_radius_m": 0.01 * outer,
                    "start_deg": start + i * span,
                    "span_deg": span,
                    "gap_m": 0.01 * gap,
                    "gap_at_start": False,
                }
                for i in range(count)
            )
        else:
            raise ModelError(f"segmentation row {row}: unsupported type {fields[0]}")
    if not segments:
        raise ModelError("segmentation contains no segments")
    return segments


def segmentation_from_record(record: dict[str, Any]) -> list[dict[str, Any]]:
    """Convert a structured simulation-models segmentation into sim_telarray rows."""
    groups = record.get("value")
    if not isinstance(groups, list) or not groups:
        raise ModelError("structured segmentation must list footprint groups")
    lines = []
    for group in groups:
        kind = group.get("kind") if isinstance(group, dict) else None
        footprint = kind in ("hex", "yhex", "square", "circle")
        names = (
            ("x", "y", "diameter", "rotation")
            if footprint
            else ("r_min", "r_max", "dphi", "phi0", "gap")
        )
        optional = {"rotation"} if footprint else {"phi0", "gap"}
        if kind not in ("hex", "yhex", "square", "circle", "ring") or set(group) - {
            "kind",
            "count",
            *names,
        }:
            raise ModelError("invalid structured segmentation fields")

        def quantity(name: str, group: dict = group, optional: set = optional) -> float:
            value = group.get(name)
            if value is None and name in optional:
                return 0.0
            if not isinstance(value, dict) or set(value) != {"value", "unit"}:
                raise ModelError(f"segmentation {name} requires value and unit")
            if name in ("rotation", "dphi", "phi0"):
                if value["unit"] != "deg":
                    raise ModelError(f"segmentation {name} requires degrees")
                return number(value["value"], name)
            return signed_length_m(value, name) * 100

        if isinstance(group.get("count"), bool) or not isinstance(group.get("count"), int):
            raise ModelError("segmentation count must be an integer")
        values = " ".join(format(quantity(name), ".17g") for name in names)
        lines.append(f"{kind} {group['count']} {values}")
    return read_segmentation("\n".join(lines))


def secondary_segment_frame(segments: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Apply sim_telarray's pi rotation about y of the secondary mirror frame."""
    result = []
    for segment in segments:
        if segment["shape"] == "polygon":
            result.append({
                **segment,
                "vertices_xy_m": [[-x, y] for x, y in segment["vertices_xy_m"]],
            })
        elif segment["shape"] == "annular_sector":
            start = 180 - segment["start_deg"] - segment["span_deg"]
            result.append({**segment, "start_deg": start, "gap_at_start": True})
        else:
            x, y = segment["centre_xy_m"]
            result.append({
                **segment,
                "centre_xy_m": [-x, y],
                "rotation_deg": 180 - segment["rotation_deg"],
            })
    return result


def read_response(
    contents: str, name: str, value_column: str = "reflectivity", *, filename_options: str = ""
) -> dict[str, Any]:
    """Read a response versus wavelength, optionally versus incidence angle.

    Returns ``{"x": wavelengths_nm, "y": incidence_deg (empty for 1-D),
    "values": [x][y] row-major, "interpolation": RPOL options}``.
    """
    rows = _data_lines(contents)
    try:
        float(rows[0][0])
        numeric_first = True
    except ValueError, IndexError:
        numeric_first = False
    if contents.lstrip().startswith("#@RPOL@") or numeric_first:
        try:
            table = parse_rpol_table(
                contents, api_options="clip=0,yscale=deg2rad", filename_options=filename_options
            )
        except TableImportError as error:
            raise ModelError(f"{name}: {error}") from error
        result = {
            "x": table["x"],
            "y": [math.degrees(v) for v in table["y"]],
            "values": table["response"],
            "interpolation": table["interpolation"],
        }
    else:
        header = rows[0]
        columns = {"wavelength", value_column} | (
            {"incidence_angle"} if "incidence_angle" in header else set()
        )
        if not columns.issubset(header) or len(header) != len(set(header)):
            raise ModelError(f"{name} lacks wavelength and {value_column} columns")
        samples: dict[tuple[float, float], list[float]] = {}
        for row in rows[1:]:
            if len(row) != len(header):
                raise ModelError(f"{name}: wrong column count")
            entry = dict(zip(header, row, strict=True))
            try:
                key = (float(entry["wavelength"]), float(entry.get("incidence_angle", 0.0)))
                samples.setdefault(key, []).append(float(entry[value_column]))
            except ValueError as error:
                raise ModelError(f"{name}: invalid number") from error
        x = sorted({k[0] for k in samples})
        y = sorted({k[1] for k in samples}) if "incidence_angle" in header else []
        if len(samples) != len(x) * max(1, len(y)):
            raise ModelError(f"{name} requires a complete wavelength-angle grid")
        # Repeated measurements of one wavelength are averaged.
        values = [sum(samples[xx, yy]) / len(samples[xx, yy]) for xx in x for yy in (y or [0.0])]
        result = {
            "x": x,
            "y": y,
            "values": values,
            "interpolation": {
                "boundary": "clamp",
                "scheme": 1,
                "x_log": False,
                "value_log": False,
                "coefficients": [],
            },
        }
    if len(result["x"]) < 2 or any(not 0.0 <= v <= 1.0 for v in result["values"]):
        raise ModelError(f"{name} needs at least two wavelengths and values in [0, 1]")
    return result


def read_spatial_map(contents: str, name: str, *, filename_options: str = "") -> dict[str, Any]:
    """Read a 2-D relative response map in cm coordinates; return metres."""
    try:
        table = parse_rpol_table(contents, default_dimension=3, filename_options=filename_options)
    except TableImportError as error:
        raise ModelError(f"{name}: {error}") from error
    if not table["y"] or any(not 0.0 <= v <= 1.0 for v in table["response"]):
        raise ModelError(f"{name} must be a 2-D table of fractions")
    ny = len(table["y"])
    return {
        "x_m": [0.01 * v for v in table["x"]],
        "y_m": [0.01 * v for v in table["y"]],
        "values": [table["response"][i * ny : (i + 1) * ny] for i in range(len(table["x"]))],
        "outside": "zero" if table["interpolation"]["boundary"] == "zero" else "clamp",
    }


def read_structure_table(contents: str, name: str) -> dict[str, Any]:
    """Read an ECSV structure table; check units (m) and the coordinate frame.

    Returns ``{"frame": str, "rows": [{column: value}]}``; numeric columns are floats.
    """
    header = None
    units: dict[str, str] = {}
    frame = None
    rows = []
    for line in contents.splitlines():
        stripped = line.strip()
        if stripped.startswith("#"):
            column = re.match(
                r"#\s*-\s*\{name:\s*(\w+),\s*(?:unit:\s*(\w+),\s*)?datatype", stripped
            )
            if column and column[2]:
                units[column[1]] = column[2]
            if match := re.match(r"#\s*coordinate_frame:\s*(\S+)", stripped):
                frame = match[1]
            continue
        if not stripped:
            continue
        fields = stripped.split()
        if header is None:
            header = fields
            continue
        if len(fields) != len(header):
            raise ModelError(f"{name} row {len(rows)}: wrong column count")
        row: dict[str, Any] = {}
        for column, field in zip(header, fields, strict=True):
            if column in ("id", "group"):
                row[column] = field
                continue
            try:
                row[column] = float(field)
            except ValueError as error:
                raise ModelError(f"{name} row {len(rows)}: invalid {column}") from error
            if not math.isfinite(row[column]):
                raise ModelError(f"{name} row {len(rows)}: non-finite {column}")
        rows.append(row)
    if header is None or not rows:
        raise ModelError(f"{name} is empty")
    for column in header:
        if column not in ("id", "group") and units.get(column) != "m":
            raise ModelError(f"{name}: column {column} must be declared in metres")
    if frame is None:
        raise ModelError(f"{name} does not declare a coordinate_frame")
    return {"frame": frame, "columns": header, "rows": rows}
