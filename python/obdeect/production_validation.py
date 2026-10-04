"""Fail-closed validation gate for model-derived CTAO optical traces.

The comparison inputs are deliberately normalised CSV tables, one row per
resolved photon.  Adapters for sim_telarray/ROBAST belong at the boundary; the
gate never guesses coordinate, time, or wavelength conventions from a native
simulator output file.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from obdeect.result_contract import TERMINAL_STATUSES


class ProductionValidationError(ValueError):
    """A production optical model or comparison input does not meet the validation contract."""


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
    wavelength_nm: float | None = None
    source_weight: float | None = None
    throughput: float | None = None
    terminal_surface_id: str | None = None
    optical_model_sha256: str | None = None
    source_sha256: str | None = None
    final_direction: tuple[float, float, float] | None = None
    interaction_surface_ids: tuple[int, ...] | None = None
    response_loss_fraction: float | None = None
    terminal_loss_fraction: float | None = None
    optical_path_m: float | None = None


def _number(row: dict[str, str], column: str, path: Path, line: int) -> float:
    try:
        value = float(row[column])
    except (KeyError, TypeError, ValueError) as error:
        raise ProductionValidationError(f"{path}:{line}: invalid {column}") from error
    if not math.isfinite(value):
        raise ProductionValidationError(f"{path}:{line}: non-finite {column}")
    return value


PhotonIdentity = int | tuple[int, ...]


def read_comparison_table(path: Path) -> dict[PhotonIdentity, ComparisonRow]:
    """Read a normalised, one-row-per-resolved-photon comparison table."""
    try:
        handle = path.open(newline="", encoding="utf-8")
    except OSError as error:
        raise ProductionValidationError(f"cannot read {path}: {error}") from error
    with handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or set(_REQUIRED_COLUMNS) - set(reader.fieldnames):
            raise ProductionValidationError(f"{path}: missing required comparison columns")
        rows: dict[PhotonIdentity, ComparisonRow] = {}
        for line, row in enumerate(reader, 2):
            try:
                photon_id = int(row["photon_id"])
            except (KeyError, TypeError, ValueError) as error:
                raise ProductionValidationError(f"{path}:{line}: invalid photon_id") from error
            context_names = ("run_id", "event_id", "array_id", "telescope_id", "bunch_id")
            identity: PhotonIdentity = photon_id
            if any(name in row for name in context_names):
                try:
                    context = tuple(int(row[name]) for name in context_names)
                except (KeyError, TypeError, ValueError) as error:
                    raise ProductionValidationError(
                        f"{path}:{line}: invalid identity context"
                    ) from error
                if any(value < 0 for value in context):
                    raise ProductionValidationError(f"{path}:{line}: negative identity context")
                identity = (*context, photon_id)
            if photon_id < 0 or identity in rows:
                raise ProductionValidationError(f"{path}:{line}: non-unique photon_id")
            status = row["status"]
            if status not in TERMINAL_STATUSES:
                raise ProductionValidationError(f"{path}:{line}: invalid status")
            values = tuple(_number(row, column, path, line) for column in _REQUIRED_COLUMNS[2:])
            extra = {}
            for column in ("wavelength_nm", "source_weight", "throughput", "optical_path_m"):
                if column in row:
                    extra[column] = _number(row, column, path, line)
            loss_fields = ("response_loss_fraction", "terminal_loss_fraction")
            if any(column in row for column in loss_fields):
                for column in loss_fields:
                    extra[column] = _number(row, column, path, line)
            for column in ("terminal_surface_id", "optical_model_sha256", "source_sha256"):
                if column in row:
                    extra[column] = row[column]
                    if not row[column]:
                        raise ProductionValidationError(f"{path}:{line}: empty {column}")
                    if column.endswith("sha256") and (
                        len(row[column]) != 64
                        or any(character not in "0123456789abcdef" for character in row[column])
                    ):
                        raise ProductionValidationError(f"{path}:{line}: invalid {column}")
            direction_columns = ("final_dx", "final_dy", "final_dz")
            if any(column in row for column in direction_columns):
                direction = tuple(_number(row, column, path, line) for column in direction_columns)
                extra["final_direction"] = direction
            if "interaction_surface_ids" in row:
                try:
                    surfaces = tuple(
                        int(value) for value in row["interaction_surface_ids"].split(";") if value
                    )
                except (AttributeError, ValueError) as error:
                    raise ProductionValidationError(
                        f"{path}:{line}: invalid interaction surface IDs"
                    ) from error
                extra["interaction_surface_ids"] = surfaces
            parsed = ComparisonRow(
                photon_id,
                status,
                *values,
                **extra,
            )
            try:
                _validate_comparison_row(parsed)
            except ProductionValidationError as error:
                raise ProductionValidationError(f"{path}:{line}: {error}") from error
            rows[identity] = parsed
    if not rows:
        raise ProductionValidationError(f"{path}: comparison table is empty")
    return rows


def _maximum(values: list[float]) -> float:
    return max(values) if values else 0.0


def _validate_comparison_row(row: ComparisonRow) -> None:
    """Apply numerical invariants to rows supplied directly by library adapters."""
    if isinstance(row.photon_id, bool) or not isinstance(row.photon_id, int) or row.photon_id < 0:
        raise ProductionValidationError("invalid photon identity")
    if row.status not in TERMINAL_STATUSES:
        raise ProductionValidationError("invalid terminal status")
    measurements = (
        row.focal_x_m,
        row.focal_y_m,
        row.path_length_m,
        row.arrival_time_ns,
        row.incidence_primary_deg,
        row.incidence_secondary_deg,
        row.incidence_focal_deg,
    )
    if (
        not all(math.isfinite(value) for value in measurements)
        or row.path_length_m < 0
        or any(not 0 <= value <= 90 for value in measurements[4:])
    ):
        raise ProductionValidationError("invalid optical measurement range")
    if row.source_weight is not None and (
        not math.isfinite(row.source_weight) or row.source_weight < 0
    ):
        raise ProductionValidationError("invalid source weight")
    if row.wavelength_nm is not None and (
        not math.isfinite(row.wavelength_nm) or row.wavelength_nm <= 0
    ):
        raise ProductionValidationError("invalid wavelength")
    if row.throughput is not None and (
        not 0 <= row.throughput <= 1 or (row.status != "detected" and row.throughput != 0)
    ):
        raise ProductionValidationError("invalid throughput")
    if row.final_direction is not None and (
        len(row.final_direction) != 3
        or not all(math.isfinite(value) for value in row.final_direction)
        or (
            row.status == "detected"
            and not math.isclose(
                sum(value * value for value in row.final_direction), 1, rel_tol=0, abs_tol=1e-9
            )
        )
    ):
        raise ProductionValidationError("invalid output direction")

    if row.optical_path_m is not None and (
        not math.isfinite(row.optical_path_m) or row.optical_path_m < 0
    ):
        raise ProductionValidationError("invalid phase optical path")
    losses = (row.response_loss_fraction, row.terminal_loss_fraction)
    if any(value is not None for value in losses):
        if (
            row.throughput is None
            or any(value is None or not math.isfinite(value) or value < 0 for value in losses)
            or not math.isclose(row.throughput + sum(losses), 1, rel_tol=0, abs_tol=1e-9)
        ):
            raise ProductionValidationError("optical loss fractions do not close")
        if row.status == "detected" and row.terminal_loss_fraction != 0:
            raise ProductionValidationError("detected photon has terminal loss")
    if row.interaction_surface_ids is not None and any(
        isinstance(value, bool) or not isinstance(value, int) or not 0 <= value < 4294967295
        for value in row.interaction_surface_ids
    ):
        raise ProductionValidationError("invalid interaction surface IDs")


def validate_optical_model(optical_model: dict[str, Any], telescope_family: str) -> None:
    """Require an intact current optical model with explicit production readiness."""
    report = optical_model.get("report")
    if not isinstance(report, dict):
        raise ProductionValidationError("compiled optical model lacks report")
    if report.get("production_trace_ready") is not True or report.get("trace_blockers"):
        raise ProductionValidationError("compiled optical model is not production ready")
    if optical_model.get("format") != "obdeect.compiled-optical-model.v1":
        raise ProductionValidationError("unsupported compiled optical model format")
    try:
        canonical = json.dumps(
            {key: value for key, value in optical_model.items() if key != "optical_model_sha256"},
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode()
    except (TypeError, ValueError) as error:
        raise ProductionValidationError(
            "compiled optical model has invalid canonical content"
        ) from error
    if optical_model.get("optical_model_sha256") != hashlib.sha256(canonical).hexdigest():
        raise ProductionValidationError("compiled optical model content hash differs")
    if report.get("native_trace_ready") is not True:
        raise ProductionValidationError("compiled optical model lacks native trace readiness")
    trace_model = optical_model.get("trace_model")
    if not isinstance(trace_model, dict):
        raise ProductionValidationError("compiled optical model lacks trace geometry")
    kind = trace_model.get("kind")
    if kind == "segmented":
        if not trace_model.get("primary_facets"):
            raise ProductionValidationError("compiled optical model lacks finite primary facets")
        if not trace_model.get("detector_surfaces"):
            raise ProductionValidationError("compiled optical model lacks detector surfaces")
    elif kind == "axisymmetric":
        for role in ("primary", "secondary", "detector"):
            if not isinstance(trace_model.get(role), dict):
                raise ProductionValidationError(
                    f"compiled optical model lacks axisymmetric {role} surface"
                )
    else:
        raise ProductionValidationError("compiled optical model has unsupported trace geometry")
    if telescope_family in {"SST", "SCT"} and kind != "axisymmetric":
        raise ProductionValidationError(
            f"{telescope_family} production optical model lacks axisymmetric secondary surfaces"
        )


def validate_comparison_fixture(fixture: dict[str, Any]) -> tuple[int, set[str]]:
    """Validate declared acceptance coverage before any simulator is executed."""
    if not isinstance(fixture, dict):
        raise ProductionValidationError("fixture must declare detection and surface coverage")
    minimum = fixture.get("minimum_detected")
    surfaces = fixture.get("required_surfaces")
    if (
        isinstance(minimum, bool)
        or not isinstance(minimum, int)
        or minimum < 0
        or (minimum == 0 and fixture.get("allow_all_loss") is not True)
        or not isinstance(surfaces, list)
        or any(not isinstance(surface, str) or not surface for surface in surfaces)
    ):
        raise ProductionValidationError("fixture must declare detection and surface coverage")
    for column in ("optical_model_sha256", "source_sha256"):
        expected = fixture.get(column)
        if (
            not isinstance(expected, str)
            or len(expected) != 64
            or any(character not in "0123456789abcdef" for character in expected)
        ):
            raise ProductionValidationError(f"fixture lacks valid {column}")
    return minimum, set(surfaces)


def compare(
    obdeect: dict[PhotonIdentity, ComparisonRow],
    simtel: dict[PhotonIdentity, ComparisonRow],
    tolerances: dict[str, float],
    fixture: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Compare identities first, then terminal states and measured optical quantities."""
    if not obdeect or not simtel:
        raise ProductionValidationError("comparison contains no resolved photons")
    for rows in (obdeect, simtel):
        contextual = isinstance(next(iter(rows)), tuple)
        for identity, row in rows.items():
            components = identity if isinstance(identity, tuple) else (identity,)
            if (
                isinstance(identity, tuple) != contextual
                or len(components) != (6 if contextual else 1)
                or any(
                    isinstance(value, bool) or not isinstance(value, int) or value < 0
                    for value in components
                )
                or components[-1] != row.photon_id
            ):
                raise ProductionValidationError("invalid photon identity binding")
            _validate_comparison_row(row)
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
    optical_path_residuals = []
    time_residuals = []
    incidence_residuals = []
    first_divergence = None
    minimum = 1
    required_surfaces: set[str] = set()
    if fixture is not None:
        minimum, required_surfaces = validate_comparison_fixture(fixture)
        for engine, rows in (("obdeect", obdeect), ("sim_telarray", simtel)):
            for photon in rows.values():
                for column in (
                    "wavelength_nm",
                    "source_weight",
                    "throughput",
                    "terminal_surface_id",
                    "final_direction",
                ):
                    if getattr(photon, column) is None:
                        raise ProductionValidationError(
                            f"{engine} photon {photon.photon_id} lacks {column}"
                        )
                for column in ("optical_model_sha256", "source_sha256"):
                    if getattr(photon, column) != fixture[column]:
                        raise ProductionValidationError(
                            f"{engine} photon {photon.photon_id}: {column} binding differs"
                        )
                if photon.status == "detected" and photon.terminal_surface_id in {
                    "unrecorded",
                    "4294967295",
                }:
                    raise ProductionValidationError(
                        f"{engine} photon {photon.photon_id}: missing detector surface identity"
                    )
    for engine, rows in (("obdeect", obdeect), ("sim_telarray", simtel)):
        detections = [row for row in rows.values() if row.status == "detected"]
        if len(detections) < minimum:
            raise ProductionValidationError(
                f"{engine}: detection coverage {len(detections)} < {minimum}"
            )
        seen = {row.terminal_surface_id for row in detections}
        if required_surfaces - seen:
            raise ProductionValidationError(f"{engine}: missing detector surface coverage")
    for photon_id in sorted(obdeect):
        left, right = obdeect[photon_id], simtel[photon_id]
        if (
            left.interaction_surface_ids != right.interaction_surface_ids
            and first_divergence is None
        ):
            left_ids = left.interaction_surface_ids or ()
            right_ids = right.interaction_surface_ids or ()
            interaction_index = next(
                (
                    index
                    for index, pair in enumerate(zip(left_ids, right_ids))
                    if pair[0] != pair[1]
                ),
                min(len(left_ids), len(right_ids)),
            )
            first_divergence = {
                "photon_id": photon_id,
                "field": "interaction_surface_ids",
                "interaction_index": interaction_index,
                "obdeect": left_ids,
                "sim_telarray": right_ids,
            }
        for column in (
            "wavelength_nm",
            "source_weight",
            "throughput",
            "terminal_surface_id",
            "response_loss_fraction",
            "terminal_loss_fraction",
        ):
            if getattr(left, column) != getattr(right, column) and first_divergence is None:
                first_divergence = {
                    "photon_id": photon_id,
                    "field": column,
                    "obdeect": getattr(left, column),
                    "sim_telarray": getattr(right, column),
                }
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
        if left.final_direction is not None or right.final_direction is not None:
            if left.final_direction is None or right.final_direction is None:
                raise ProductionValidationError(f"photon {photon_id}: missing detector direction")
            direction_residual = math.dist(left.final_direction, right.final_direction)
            # Incidence tolerance in degrees also bounds the outgoing angular residual.
            if (
                direction_residual
                > 2 * math.sin(math.radians(min(tolerances["incidence_deg"], 180)) / 2)
                and first_divergence is None
            ):
                first_divergence = {
                    "photon_id": photon_id,
                    "field": "final_direction",
                    "residual": direction_residual,
                }
        if left.optical_path_m is not None or right.optical_path_m is not None:
            if left.optical_path_m is None or right.optical_path_m is None:
                raise ProductionValidationError(f"photon {photon_id}: missing phase optical path")
            phase_path = abs(left.optical_path_m - right.optical_path_m)
            optical_path_residuals.append(phase_path)
            if phase_path > tolerances["path_length_m"] and first_divergence is None:
                first_divergence = {
                    "photon_id": photon_id,
                    "field": "optical_path_m",
                    "residual": phase_path,
                }
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
        "detected_count": sum(row.status == "detected" for row in obdeect.values()),
        "status_mismatch_count": len(status_mismatches),
        "max_focal_position_residual_m": _maximum(focal_residuals),
        "max_path_length_residual_m": _maximum(path_residuals),
        "max_optical_path_residual_m": _maximum(optical_path_residuals),
        "max_arrival_time_residual_ns": _maximum(time_residuals),
        "max_incidence_residual_deg": _maximum(incidence_residuals),
        "first_divergence": first_divergence,
    }
    if all(
        row.source_weight is not None and row.throughput is not None for row in obdeect.values()
    ):
        incident = sum(row.source_weight for row in obdeect.values())
        arriving = sum(row.source_weight * row.throughput for row in obdeect.values())
        summary["incident_weight"] = incident
        summary["arriving_weight"] = arriving
        summary["lost_optical_weight"] = incident - arriving
        if all(
            row.response_loss_fraction is not None and row.terminal_loss_fraction is not None
            for row in obdeect.values()
        ):
            summary["response_loss_weight"] = sum(
                row.source_weight * row.response_loss_fraction for row in obdeect.values()
            )
            summary["terminal_loss_weight"] = sum(
                row.source_weight * row.terminal_loss_fraction for row in obdeect.values()
            )
    if status_mismatches or first_divergence is not None:
        raise ProductionValidationError("comparison failed: " + json.dumps(summary, sort_keys=True))
    return summary


def main() -> None:
    """Validate one production optical model and a normalised sim_telarray comparison."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--optical-model", type=Path, required=True, help="compiled obdeect optical model JSON"
    )
    parser.add_argument("--telescope-family", choices=("LST", "MST", "SST", "SCT"), required=True)
    parser.add_argument(
        "--fixture", type=Path, required=True, help="declared coverage and input hashes JSON"
    )
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
        optical_model = json.loads(args.optical_model.read_text(encoding="utf-8"))
        tolerances = json.loads(args.tolerances.read_text(encoding="utf-8"))
        fixture = json.loads(args.fixture.read_text(encoding="utf-8"))
        if not isinstance(optical_model, dict) or not isinstance(tolerances, dict):
            raise ProductionValidationError("optical model and tolerances must be JSON objects")
        validate_optical_model(optical_model, args.telescope_family)
        if not isinstance(fixture, dict) or fixture.get(
            "optical_model_sha256"
        ) != optical_model.get("optical_model_sha256"):
            raise ProductionValidationError("fixture optical model hash binding differs")
        summary = compare(
            read_comparison_table(args.obdeect_arrivals),
            read_comparison_table(args.simtel_arrivals),
            tolerances,
            fixture,
        )
    except (OSError, json.JSONDecodeError, ProductionValidationError) as error:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps({"passed": False, "error": str(error)}, indent=2) + "\n", encoding="utf-8"
        )
        raise SystemExit(f"production validation failed: {error}") from error
    summary["passed"] = True
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Validated {args.telescope_family} production trace: {args.output}")


if __name__ == "__main__":
    main()
