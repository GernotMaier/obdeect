"""Parse the documented geometry subset of a sim_telarray camera file."""

from __future__ import annotations

import math
import shlex
from collections import Counter
from typing import Any


class CameraConfigError(ValueError):
    """Camera geometry is malformed or uses an unsupported directive."""


def _integer(value: str, line: int) -> int:
    try:
        number = int(value)
    except ValueError as error:
        raise CameraConfigError(f"camera line {line}: invalid integer") from error
    if number < 0:
        raise CameraConfigError(f"camera line {line}: negative identifier")
    return number


def _number(value: str, line: int) -> float:
    try:
        number = float(value)
    except ValueError as error:
        raise CameraConfigError(f"camera line {line}: invalid number") from error
    if not math.isfinite(number):
        raise CameraConfigError(f"camera line {line}: non-finite number")
    return number


def parse_camera_layout(contents: str) -> dict[str, Any]:
    """Return pixel centres/types and retain directives outside geometry scope."""
    pixel_types: dict[int, dict[str, Any]] = {}
    pixels: list[dict[str, Any]] = []
    seen_pixels: set[int] = set()
    deferred: Counter[str] = Counter()
    rotation_deg = 0.0
    for line_number, line in enumerate(contents.splitlines(), start=1):
        try:
            fields = shlex.split(line, comments=True)
        except ValueError as error:
            raise CameraConfigError(f"camera line {line_number}: invalid quoting") from error
        if not fields:
            continue
        directive = fields[0]
        if directive == "PixType":
            if len(fields) < 9:
                raise CameraConfigError(f"camera line {line_number}: short PixType")
            identifier = _integer(fields[1], line_number)
            if identifier in pixel_types:
                raise CameraConfigError(f"camera line {line_number}: duplicate PixType")
            shape = _integer(fields[3], line_number)
            funnel_shape = _integer(fields[5], line_number)
            diameter = _number(fields[4], line_number)
            funnel_diameter = _number(fields[6], line_number)
            depth = _number(fields[7], line_number)
            if shape > 3 or funnel_shape > 3 or diameter <= 0 or funnel_diameter <= 0 or depth < 0:
                raise CameraConfigError(f"camera line {line_number}: invalid pixel footprint")
            references = []
            for field in fields[8:10]:
                try:
                    float(field)
                except ValueError:
                    references.append(field)
            pixel_types[identifier] = {
                "id": identifier,
                "cathode_shape_code": shape,
                "cathode_diameter_m": diameter * 0.01,
                "funnel_shape_code": funnel_shape,
                "funnel_diameter_m": funnel_diameter * 0.01,
                "funnel_depth_m": depth * 0.01,
                "response_files": references,
                "source_columns": fields[2:],
            }
        elif directive == "Pixel":
            if len(fields) < 5:
                raise CameraConfigError(f"camera line {line_number}: short Pixel")
            identifier = _integer(fields[1], line_number)
            if identifier in seen_pixels:
                raise CameraConfigError(f"camera line {line_number}: duplicate Pixel")
            seen_pixels.add(identifier)
            pixels.append({
                "id": identifier,
                "type_id": _integer(fields[2], line_number),
                "centre_xy_m": [
                    _number(fields[3], line_number) * 0.01,
                    _number(fields[4], line_number) * 0.01,
                ],
                "source_columns": fields[5:],
            })
        elif directive == "Rotate":
            if len(fields) != 2:
                raise CameraConfigError(f"camera line {line_number}: invalid Rotate")
            rotation_deg += _number(fields[1], line_number)
        elif directive in {"AnalogSumTrigger", "DigitalSumTrigger", "MajorityTrigger", "Trigger"}:
            deferred[directive] += 1
        else:
            raise CameraConfigError(f"camera line {line_number}: unsupported directive {directive}")
    if not pixels or not pixel_types:
        raise CameraConfigError("camera has no pixels or pixel types")
    if any(pixel["type_id"] not in pixel_types for pixel in pixels):
        raise CameraConfigError("pixel refers to undefined PixType")
    return {
        "kind": "pixel_layout",
        "pixel_types": [pixel_types[key] for key in sorted(pixel_types)],
        "pixels": pixels,
        "rotation_deg": rotation_deg,
        "deferred_directives": dict(sorted(deferred.items())),
    }


def parse_camera_layout_ecsv(contents: str) -> dict[str, Any]:
    """Read the geometric subset of a simulation-models camera-layout ECSV.

    Only ID and focal-plane coordinates are consumed. The table describes the
    physical layout and is not interpreted as readout or response behaviour.
    """
    header: list[str] | None = None
    entries: list[dict[str, Any]] = []
    identifiers: set[int] = set()
    for line_number, source_line in enumerate(contents.splitlines(), start=1):
        line = source_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if header is None:
            header = fields
            required = {"pixel_id", "type_id", "x_cm", "y_cm"}
            if not required.issubset(header):
                raise CameraConfigError("camera ECSV lacks geometric columns")
            continue
        if len(fields) != len(header):
            raise CameraConfigError(f"camera ECSV line {line_number}: wrong column count")
        values = dict(zip(header, fields, strict=True))
        identifier = _integer(values["pixel_id"], line_number)
        if identifier in identifiers:
            raise CameraConfigError(f"camera ECSV line {line_number}: duplicate identifier")
        identifiers.add(identifier)
        entries.append({
            "id": identifier,
            "type_id": _integer(values["type_id"], line_number),
            "centre_xy_m": [
                _number(values["x_cm"], line_number) * 0.01,
                _number(values["y_cm"], line_number) * 0.01,
            ],
        })
    if header is None or not entries:
        raise CameraConfigError("camera ECSV has no layout entries")
    return {
        "kind": "focal_plane_layout",
        "pixel_types": [],
        "pixels": entries,
        "rotation_deg": 0.0,
        "deferred_directives": {},
    }
