"""Fail-closed validation gate for model-derived CTAO optical traces.

The comparison inputs are deliberately normalised CSV tables, one row per
resolved photon.  Adapters for sim_telarray/ROBAST belong at the boundary; the
gate never guesses coordinate, time, or wavelength conventions from a native
simulator output file.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from obdeect.result_contract import TERMINAL_STATUSES


class ProductionValidationError(ValueError):
    """A production scene or comparison input does not meet the validation contract."""


_REQUIRED_COLUMNS = (
    "photon_id",
    "status",
    "focal_x_m",
    "focal_y_m",
    "path_length_m",
    "arrival_time_ns",
    "incidence_primary_deg",
    "incidence_secondary_deg",
    "incidence_focal_deg",
)


@dataclass(frozen=True)
class ComparisonRow:
    """One resolved input photon measured by one optical engine."""

    photon_id: int
    status: str
    focal_x_m: float
    focal_y_m: float
    path_length_m: float
    arrival_time_ns: float
    incidence_primary_deg: float
    incidence_secondary_deg: float
    incidence_focal_deg: float


def _number(row: dict[str, str], column: str, path: Path, line: int) -> float:
    try:
        value = float(row[column])
    except (KeyError, TypeError, ValueError) as error:
        raise ProductionValidationError(f"{path}:{line}: invalid {column}") from error
    if not math.isfinite(value):
        raise ProductionValidationError(f"{path}:{line}: non-finite {column}")
    return value


def read_comparison_table(path: Path) -> dict[int, ComparisonRow]:
    """Read a normalised, one-row-per-resolved-photon comparison table."""
    try:
        handle = path.open(newline="", encoding="utf-8")
    except OSError as error:
        raise ProductionValidationError(f"cannot read {path}: {error}") from error
    with handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or set(_REQUIRED_COLUMNS) - set(reader.fieldnames):
            raise ProductionValidationError(f"{path}: missing required comparison columns")
        rows: dict[int, ComparisonRow] = {}
        for line, row in enumerate(reader, 2):
            try:
                photon_id = int(row["photon_id"])
            except (KeyError, TypeError, ValueError) as error:
                raise ProductionValidationError(f"{path}:{line}: invalid photon_id") from error
            if photon_id < 0 or photon_id in rows:
                raise ProductionValidationError(f"{path}:{line}: non-unique photon_id")
            status = row["status"]
            if status not in TERMINAL_STATUSES:
                raise ProductionValidationError(f"{path}:{line}: invalid status")
            rows[photon_id] = ComparisonRow(
                photon_id,
                status,
                *(_number(row, column, path, line) for column in _REQUIRED_COLUMNS[2:]),
            )
    if not rows:
        raise ProductionValidationError(f"{path}: comparison table is empty")
    return rows


def _maximum(values: list[float]) -> float:
    return max(values) if values else 0.0


def validate_scene(scene: dict[str, Any], telescope_family: str) -> None:
    """Require the compiled scene to contain the physical model for its family."""
    report = scene.get("report")
    if not isinstance(report, dict):
        raise ProductionValidationError("compiled scene lacks report")
    if report.get("trace_blockers"):
        raise ProductionValidationError(
            "compiled scene is not production ready: " + "; ".join(report["trace_blockers"])
        )
    if report.get("native_trace_ready") is not True:
        raise ProductionValidationError("compiled scene has no production native binding")
    primary = scene.get("primary", {})
    if not isinstance(primary, dict) or not primary.get("facets"):
        raise ProductionValidationError("compiled scene lacks finite primary facets")
    if telescope_family == "SST":
        secondary = scene.get("secondary", {})
        if not isinstance(secondary, dict) or not secondary.get("facets"):
            raise ProductionValidationError("SST production scene lacks finite secondary facets")
    if not isinstance(scene.get("detector"), dict):
        raise ProductionValidationError("compiled scene lacks a physical detector surface")
    if not isinstance(scene.get("materials"), dict):
        raise ProductionValidationError("compiled scene lacks material bindings")


def compare(
    obdeect: dict[int, ComparisonRow],
    simtel: dict[int, ComparisonRow],
    tolerances: dict[str, float],
) -> dict[str, Any]:
    """Compare identities first, then terminal states and measured optical quantities."""
    if set(obdeect) != set(simtel):
        missing_obdeect = sorted(set(simtel) - set(obdeect))
        missing_simtel = sorted(set(obdeect) - set(simtel))
        raise ProductionValidationError(
            f"resolved photon identities differ: missing_obdeect={missing_obdeect[:5]}, "
            f"missing_simtel={missing_simtel[:5]}"
        )
    allowed = {"focal_position_m", "path_length_m", "arrival_time_ns", "incidence_deg"}
    if set(tolerances) != allowed or any(
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value < 0
        for value in tolerances.values()
    ):
        raise ProductionValidationError(
            "tolerances must define non-negative focal_position_m, path_length_m, "
            "arrival_time_ns, incidence_deg"
        )
    status_mismatches = []
    focal_residuals = []
    path_residuals = []
    time_residuals = []
    incidence_residuals = []
    first_divergence = None
    for photon_id in sorted(obdeect):
        left, right = obdeect[photon_id], simtel[photon_id]
        if left.status != right.status:
            status_mismatches.append(photon_id)
            if first_divergence is None:
                first_divergence = {
                    "photon_id": photon_id,
                    "field": "status",
                    "obdeect": left.status,
                    "sim_telarray": right.status,
                }
            continue
        # A terminal loss has no focal-plane or incidence measurement.  Its
        # path is useful diagnostic information but is not an equivalence
        # observable: different transport implementations can terminate at
        # different points on the same rejected trajectory.  Require its
        # stable category above, and reserve numerical optical checks for
        # photons that both engines actually detected.
        if left.status != "detected":
            continue
        focal = math.hypot(left.focal_x_m - right.focal_x_m, left.focal_y_m - right.focal_y_m)
        path = abs(left.path_length_m - right.path_length_m)
        arrival_time = abs(left.arrival_time_ns - right.arrival_time_ns)
        incidence = max(
            abs(left.incidence_primary_deg - right.incidence_primary_deg),
            abs(left.incidence_secondary_deg - right.incidence_secondary_deg),
            abs(left.incidence_focal_deg - right.incidence_focal_deg),
        )
        focal_residuals.append(focal)
        path_residuals.append(path)
        time_residuals.append(arrival_time)
        incidence_residuals.append(incidence)
        if first_divergence is None and (
            focal > tolerances["focal_position_m"]
            or path > tolerances["path_length_m"]
            or arrival_time > tolerances["arrival_time_ns"]
            or incidence > tolerances["incidence_deg"]
        ):
            first_divergence = {
                "photon_id": photon_id,
                "field": "optical_residual",
                "focal_position_m": focal,
                "path_length_m": path,
                "arrival_time_ns": arrival_time,
                "incidence_deg": incidence,
            }
    summary = {
        "photon_count": len(obdeect),
        "status_mismatch_count": len(status_mismatches),
        "max_focal_position_residual_m": _maximum(focal_residuals),
        "max_path_length_residual_m": _maximum(path_residuals),
        "max_arrival_time_residual_ns": _maximum(time_residuals),
        "max_incidence_residual_deg": _maximum(incidence_residuals),
        "first_divergence": first_divergence,
    }
    if status_mismatches or first_divergence is not None:
        raise ProductionValidationError("comparison failed: " + json.dumps(summary, sort_keys=True))
    return summary


def main() -> None:
    """Validate one production scene and a normalised sim_telarray comparison."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scene", type=Path, required=True, help="compiled obdeect scene JSON")
    parser.add_argument("--telescope-family", choices=("LST", "MST", "SST"), required=True)
    parser.add_argument(
        "--obdeect-arrivals", type=Path, required=True, help="normalised obdeect comparison CSV"
    )
    parser.add_argument(
        "--simtel-arrivals", type=Path, required=True, help="normalised sim_telarray comparison CSV"
    )
    parser.add_argument(
        "--tolerances",
        type=Path,
        required=True,
        help="JSON object of declared comparison tolerances",
    )
    parser.add_argument("--output", type=Path, required=True, help="validation summary JSON")
    args = parser.parse_args()
    try:
        scene = json.loads(args.scene.read_text(encoding="utf-8"))
        tolerances = json.loads(args.tolerances.read_text(encoding="utf-8"))
        if not isinstance(scene, dict) or not isinstance(tolerances, dict):
            raise ProductionValidationError("scene and tolerances must be JSON objects")
        validate_scene(scene, args.telescope_family)
        summary = compare(
            read_comparison_table(args.obdeect_arrivals),
            read_comparison_table(args.simtel_arrivals),
            tolerances,
        )
    except (OSError, json.JSONDecodeError, ProductionValidationError) as error:
        raise SystemExit(f"production validation failed: {error}") from error
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Validated {args.telescope_family} production trace: {args.output}")


if __name__ == "__main__":
    main()
