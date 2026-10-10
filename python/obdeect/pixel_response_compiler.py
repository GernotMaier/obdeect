"""Resolve sim_telarray's measured and single-reflection pixel response contracts."""

from __future__ import annotations

import bisect
import math
import re
from typing import Any

from obdeect.simtel_tables import TableImportError, parse_rpol_table


def _curve(contents: str) -> dict[str, Any]:
    numeric = []
    for line in contents.splitlines():
        if not line.strip() or line.lstrip().startswith(("#", "%")):
            continue
        try:
            float(line.split()[0])
        except ValueError:
            continue  # A named ECSV column header; all data rows remain explicit.
        numeric.append(line)
    if len(numeric) == 1:
        fields = numeric[0].split("#", 1)[0].split("%", 1)[0].split()
        if len(fields) != 2:
            raise TableImportError("constant pixel response requires two columns")
        coordinate, response = map(float, fields)
        if not math.isfinite(coordinate) or not math.isfinite(response):
            raise TableImportError("pixel response values must be finite")
        return {"x": [coordinate], "y": [], "response": [response]}
    table = parse_rpol_table("\n".join(numeric))
    if table["y"]:
        raise TableImportError("pixel response requires a one-dimensional curve")
    return table


def _linear(table, coordinate):
    axis, values = table["x"], table["response"]
    if coordinate <= axis[0]:
        return values[0]
    if coordinate >= axis[-1]:
        return values[-1]
    high = bisect.bisect_right(axis, coordinate)
    fraction = (coordinate - axis[high - 1]) / (axis[high] - axis[high - 1])
    return (1 - fraction) * values[high - 1] + fraction * values[high]


def measured_pixel_response(
    type_id: int, angular_contents: str, wavelength_contents: str | None = None
) -> dict[str, Any]:
    """Freeze the reference PixType tangent and integer-nanometre lookup tables.

    Sampling is the sim_telarray convention in sim_imaging.c camera reading;
    the generic native kernel receives its bin widths and arrays explicitly.
    Wavelength curves are absolute average efficiencies, divided by the mean
    angular efficiency over the file's #THETA range (100 sine-weighted samples).
    """
    angle = _curve(angular_contents)
    if any(not 0 <= value <= 1 for value in angle["response"]):
        raise TableImportError("angular pixel efficiency must be a fraction")
    angular_bins, tangent_span = 1000, 10.0
    table = {
        "id": type_id,
        "method": "measured",
        "tangent_bin_width": tangent_span / angular_bins,
        "angular_efficiency": [
            _linear(angle, math.degrees(math.atan(tangent_span * (i + 0.5) / angular_bins)))
            for i in range(angular_bins)
        ],
    }
    if wavelength_contents is not None:
        wavelength = _curve(wavelength_contents)
        match = re.search(
            r"^#THETA:\s*([+\-\d.eE]+)(?:\s+([+\-\d.eE]+))?", wavelength_contents, re.M | re.I
        )
        lower = float(match[1]) if match else 0.0
        upper = float(match[2] or match[1]) if match else 0.0
        angles = [lower + (i + 0.5) * (upper - lower) / 100 for i in range(100)]
        weights = [math.sin(math.radians(abs(a))) if a else 1e-6 for a in angles]
        total = sum(weights)
        average = sum(w * _linear(angle, a) for w, a in zip(weights, angles)) / total
        if average == 0:
            average = 1e30  # The reference's explicitly zero-efficiency convention.
        table.update(
            wavelength_bin_width_nm=1.0,
            wavelength_bin_origin_nm=0.5,
            spectral_correction=[_linear(wavelength, i + 0.5) / average for i in range(1000)],
        )
    return table


def bind_pixel_responses(camera: dict[str, Any], trace: dict[str, Any]) -> dict[str, Any] | None:
    tables = camera.get("pixel_optical_response_tables", [])
    if not tables:
        return None
    by_type = {entry["id"]: entry for entry in tables}
    cathodes = {entry["source_pixel_id"]: entry for entry in camera["cathode_surfaces"]}
    bindings = []
    for detector in trace["detector_surfaces"]:
        type_id = detector["source_type_id"]
        if type_id not in by_type:
            raise TableImportError(f"physical pixel type {type_id} has no optical response")
        binding = {"id": detector["id"], "table_id": type_id}
        if by_type[type_id]["method"] == "single_reflection":
            cathode = cathodes[detector["source_pixel_id"]]
            binding["cathode"] = {
                key: cathode[key]
                for key in ("shape", "centre_m", "normal", "tangent", "diameter_m")
            }
        bindings.append(binding)
    return {"tables": tables, "bindings": bindings}
