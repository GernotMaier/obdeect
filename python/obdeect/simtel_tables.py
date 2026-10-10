"""Compile the RPOL table conventions documented in sim_telarray manual 11.2–11.3."""

from __future__ import annotations

import math
import re
from typing import Any


class TableImportError(ValueError):
    """An RPOL table cannot be represented without discarding a declared option."""


def _polynomial_product(first, second):
    result = [0.0] * (len(first) + len(second) - 1)
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            result[i + j] += a * b
    return result


def _quadratic(x, values, origin, width):
    result = [0.0] * 3
    for i in range(3):
        polynomial = [values[i]]
        for j in range(3):
            if i != j:
                polynomial = _polynomial_product(
                    polynomial, [(origin - x[j]) / (x[i] - x[j]), width / (x[i] - x[j])]
                )
        for j, coefficient in enumerate(polynomial):
            result[j] += coefficient
    return result


def interval_coefficients(axis, values, scheme):
    """Normalized-interval cubic coefficients, resolved once before native tracing."""
    if scheme in (0, 1):
        return []
    count = len(axis)
    if count < 3:
        return [[a, b - a, 0.0, 0.0] for a, b in zip(values, values[1:])]
    if scheme == 2:
        coefficients = []
        for i in range(count - 1):
            left = max(0, i - 1)
            right = min(i, count - 3)
            width = axis[i + 1] - axis[i]
            a = _quadratic(axis[left : left + 3], values[left : left + 3], axis[i], width)
            b = _quadratic(axis[right : right + 3], values[right : right + 3], axis[i], width)
            coefficients.append([a[0], a[1] + b[0] - a[0], a[2] + b[1] - a[1], b[2] - a[2]])
        return coefficients
    # The reference falls back to its piecewise polynomial when fewer than
    # four supporting points are available for a cubic spline.
    if count < 4:
        return interval_coefficients(axis, values, 2)
    widths = [b - a for a, b in zip(axis, axis[1:])]
    rhs = [0.0] * count
    for i in range(1, count - 1):
        rhs[i] = 3 * (
            (values[i + 1] - values[i]) / widths[i] - (values[i] - values[i - 1]) / widths[i - 1]
        )
    diagonal, upper, solution = [0.0] * count, [0.0] * count, [0.0] * count
    if scheme == 4:
        rhs[0] = 3 * (values[1] - values[0]) / widths[0]
        rhs[-1] = -3 * (values[-1] - values[-2]) / widths[-1]
        diagonal[0], upper[0] = 2 * widths[0], 0.5
        solution[0] = rhs[0] / diagonal[0]
    else:
        diagonal[0] = 1.0
    for i in range(1, count - 1):
        diagonal[i] = 2 * (axis[i + 1] - axis[i - 1]) - widths[i - 1] * upper[i - 1]
        upper[i] = widths[i] / diagonal[i]
        solution[i] = (rhs[i] - widths[i - 1] * solution[i - 1]) / diagonal[i]
    second = [0.0] * count
    if scheme == 4:
        diagonal[-1] = widths[-1] * (2 - upper[-2])
        solution[-1] = (rhs[-1] - widths[-1] * solution[-2]) / diagonal[-1]
        second[-1] = solution[-1]
    coefficients = [None] * (count - 1)
    for i in range(count - 2, -1, -1):
        second[i] = solution[i] - upper[i] * second[i + 1]
        slope = (values[i + 1] - values[i]) / widths[i] - widths[i] * (
            second[i + 1] + 2 * second[i]
        ) / 3
        third = (second[i + 1] - second[i]) / (3 * widths[i])
        coefficients[i] = [
            values[i],
            slope * widths[i],
            second[i] * widths[i] ** 2,
            third * widths[i] ** 3,
        ]
    return coefficients


def parse_rpol_table(
    contents: str,
    *,
    default_dimension: int = 1,
    filename_options: str = "",
    api_options: str = "",
    default_boundary: str = "clamp",
) -> dict[str, Any]:
    """Read all three table formats; retain declared interpolation and boundary rules.

    API options override the header, and filename options override both.
    The reference scales logarithms after taking them.
    Column-count and allocation hints have no effect on the compiled kernel.
    """
    lines = contents.splitlines()
    header = next((line.strip() for line in lines if line.strip()), "")
    match = re.match(r"#@RPOL@(?:\[(.*?)\])?\s*([123])(?:\s|$)(.*)", header)
    dimension = int(match[2]) if match else default_dimension
    markup = match[1] if match else None
    text = " ".join((match[3] if match else "", api_options, filename_options))
    options = dict(boundary=default_boundary, scheme=1, x_log=False, y_log=False, value_log=False)
    columns, scales = dict(x=0, y=1, z=2), dict(x=1.0, y=1.0, z=1.0)
    for token in re.split(r"[,\s]+", text.replace("OPTIONS:", "").strip()):
        if not token:
            continue
        key, _, value = token.lower().partition("=")
        if key == "clip":
            options["boundary"] = "clamp" if value in ("0", "off", "no") else "zero"
        elif key == "noclip":
            options["boundary"] = "clamp"
        elif key == "scheme":
            try:
                options["scheme"] = int(value)
            except ValueError as error:
                raise TableImportError("invalid RPOL interpolation scheme") from error
            if not 0 <= options["scheme"] <= 4:
                raise TableImportError("unsupported RPOL interpolation scheme")
        elif key in ("xlog", "ylog", "zlog"):
            options[
                {
                    "xlog": "x_log",
                    "ylog": "y_log" if dimension != 1 else "value_log",
                    "zlog": "value_log",
                }[key]
            ] = True
        elif key in ("xcol", "ycol", "zcol"):
            try:
                columns[key[0]] = int(value) - 1
            except ValueError as error:
                raise TableImportError("invalid RPOL column index") from error
            if columns[key[0]] < 0:
                raise TableImportError("RPOL columns start at one")
            if dimension == 1 and key == "zcol":
                columns["y"] = columns["z"]
        elif key in ("xscale", "yscale", "zscale"):
            try:
                scales[key[0]] = {"deg2rad": math.pi / 180, "rad2deg": 180 / math.pi}.get(
                    value, float(value) if value not in ("deg2rad", "rad2deg") else 1
                )
            except ValueError as error:
                raise TableImportError("invalid RPOL scale") from error
            if not math.isfinite(scales[key[0]]):
                raise TableImportError("RPOL scale must be finite")
            scales[key[0]] = scales[key[0]] or 1.0
        elif key not in ("zxmax", "zxmin", "verbose", "rows", "cols", "columns"):
            raise TableImportError(f"unknown RPOL option: {token}")
    rows = []
    for line in lines:
        line = line.strip()
        if markup and line.startswith(markup):
            line = line[len(markup) :].strip()
        elif not line or line.startswith(("#", "%")):
            continue
        try:
            row = [float(v) for v in line.split("#", 1)[0].split("%", 1)[0].split()]
        except ValueError as error:
            raise TableImportError("RPOL rows must be numeric") from error
        if not row or any(not math.isfinite(v) for v in row):
            raise TableImportError("RPOL rows must be finite")
        rows.append(row)

    def scaled(value, coordinate):
        logarithmic = options[{"x": "x_log", "y": "y_log", "z": "value_log"}[coordinate]]
        if logarithmic:
            if value <= 0:
                raise TableImportError("logarithmic RPOL values must be positive")
            return math.exp(math.log(value) * scales[coordinate])
        return value * scales[coordinate]

    try:
        if dimension == 1:
            samples = [
                (scaled(row[columns["x"]], "x"), scaled(row[columns["y"]], "z")) for row in rows
            ]
            axis = [x for x, _ in samples]
            values = [v for _, v in samples]
            if len(axis) < 2 or any(a >= b for a, b in zip(axis, axis[1:])):
                raise TableImportError("RPOL axis must increase strictly")
            if options["x_log"] and min(axis) <= 0 or options["value_log"] and min(values) <= 0:
                raise TableImportError("logarithmic RPOL values must be positive")
            x = [math.log(v) for v in axis] if options["x_log"] else axis
            y = [math.log(v) for v in values] if options["value_log"] else values
            options["coefficients"] = interval_coefficients(x, y, options["scheme"])
            return dict(x=axis, y=[], response=values, interpolation=options)
        if dimension == 2:
            y = [scaled(v, "y") for v in rows[0]]
            if any(len(row) != len(y) + 1 for row in rows[1:]):
                raise TableImportError("RPOL matrix has invalid columns")
            samples = [
                (scaled(row[0], "x"), yy, scaled(value, "z"))
                for row in rows[1:]
                for yy, value in zip(y, row[1:])
            ]
        else:
            samples = [
                (
                    scaled(row[columns["x"]], "x"),
                    scaled(row[columns["y"]], "y"),
                    scaled(row[columns["z"]], "z"),
                )
                for row in rows
            ]
    except (IndexError, KeyError) as error:
        raise TableImportError("RPOL rows lack their requested columns") from error
    grid = {(x, y): value for x, y, value in samples}
    x, y = sorted({p[0] for p in samples}), sorted({p[1] for p in samples})
    if not grid or len(grid) != len(samples) or len(grid) != len(x) * len(y):
        raise TableImportError("RPOL table must be a complete grid with unique coordinates")
    if (
        options["x_log"]
        and min(x) <= 0
        or options["y_log"]
        and min(y) <= 0
        or options["value_log"]
        and min(grid.values()) <= 0
    ):
        raise TableImportError("logarithmic RPOL values must be positive")
    options["scheme"] = 1  # Reference implementation always uses bilinear 2-D interpolation.
    options["coefficients"] = []
    return dict(x=x, y=y, response=[grid[xx, yy] for xx in x for yy in y], interpolation=options)
