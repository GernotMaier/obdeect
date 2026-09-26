"""Backend-neutral optical arrival records.

The contract deliberately stops at the optical detector surface.  It is small
enough for CSV interchange and strict enough for simtools analysis to reject
ambiguous or physically invalid records before calculating PSF observables.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path


class ArrivalContractError(ValueError):
    """An arrival file does not satisfy the optical result contract."""


@dataclass(frozen=True)
class OpticalArrival:
    """One resolved photon at the optical boundary, or one terminal loss."""

    photon_id: int
    source_kind: str
    wavelength_nm: float
    emission_time_ns: float
    source_weight: float
    throughput: float
    status: str
    path_length_m: float
    focal_x_m: float | None
    focal_y_m: float | None
    focal_z_m: float | None
    interaction_points_m: tuple[tuple[float, float, float], ...]
    incidence_primary_deg: float | None = None
    incidence_secondary_deg: float | None = None
    incidence_focal_deg: float | None = None

    @property
    def detected(self) -> bool:
        """Whether the photon reached the optical detector surface."""

        return self.status == "detected"

    @property
    def optical_weight(self) -> float:
        """Input weight after optical transmission."""

        return self.source_weight * self.throughput


_REQUIRED = {
    "contract_version",
    "photon_id",
    "source_kind",
    "wavelength_nm",
    "emission_time_ns",
    "source_weight",
    "throughput",
    "status",
    "point_count",
    "path_length_m",
}

_STATUSES = {
    "detected",
    "blocked_camera",
    "blocked_mast",
    "blocked_obscurer",
    "missed_primary",
    "missed_screen",
    "no_detector",
    "invalid_input",
    "escaped_scene",
    "interaction_limit",
}


def _number(row: dict[str, str], name: str, path: Path, line: int) -> float:
    try:
        value = float(row[name])
    except (KeyError, TypeError, ValueError) as error:
        raise ArrivalContractError(f"{path}:{line}: invalid {name}") from error
    if not math.isfinite(value):
        raise ArrivalContractError(f"{path}:{line}: non-finite {name}")
    return value


def _point(row: dict[str, str], index: int, path: Path, line: int) -> tuple[float, float, float]:
    return (
        _number(row, f"x{index}_m", path, line),
        _number(row, f"y{index}_m", path, line),
        _number(row, f"z{index}_m", path, line),
    )


def read_arrivals(path: Path) -> list[OpticalArrival]:
    """Read and validate an obdeect optical-arrival CSV."""

    path = Path(path)
    try:
        handle = path.open(newline="", encoding="utf-8")
    except OSError as error:
        raise ArrivalContractError(f"cannot read arrival file {path}: {error}") from error
    with handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or ())
        missing = _REQUIRED - fields
        if missing:
            raise ArrivalContractError(f"{path}: missing columns: {', '.join(sorted(missing))}")
        arrivals: list[OpticalArrival] = []
        for line, row in enumerate(reader, start=2):
            if row["contract_version"] != "obdeect-arrival-v1":
                raise ArrivalContractError(f"{path}:{line}: unsupported contract_version")
            try:
                photon_id = int(row["photon_id"])
                point_count = int(row["point_count"])
            except (TypeError, ValueError) as error:
                raise ArrivalContractError(f"{path}:{line}: invalid integer field") from error
            if photon_id < 0 or point_count < 1 or point_count > 4:
                raise ArrivalContractError(f"{path}:{line}: invalid photon_id or point_count")
            wavelength = _number(row, "wavelength_nm", path, line)
            source_kind = row["source_kind"]
            if source_kind not in {"star", "illuminator", "laser"}:
                raise ArrivalContractError(f"{path}:{line}: invalid source_kind")
            status = row["status"]
            if status not in _STATUSES:
                raise ArrivalContractError(f"{path}:{line}: invalid status")
            emission_time = _number(row, "emission_time_ns", path, line)
            source_weight = _number(row, "source_weight", path, line)
            throughput = _number(row, "throughput", path, line)
            path_length = _number(row, "path_length_m", path, line)
            incidence = tuple(
                _number(row, name, path, line)
                if name in fields and row.get(name) not in (None, "")
                else None
                for name in (
                    "incidence_primary_deg",
                    "incidence_secondary_deg",
                    "incidence_focal_deg",
                )
            )
            if wavelength <= 0 or source_weight < 0 or not 0 <= throughput <= 1 or path_length < 0:
                raise ArrivalContractError(f"{path}:{line}: invalid optical scalar")
            if any(value is not None and not 0 <= value <= 90 for value in incidence):
                raise ArrivalContractError(f"{path}:{line}: invalid incidence angle")
            points = tuple(_point(row, index, path, line) for index in range(point_count))
            focal = points[-1] if status == "detected" else (None, None, None)
            arrivals.append(
                OpticalArrival(
                    photon_id,
                    source_kind,
                    wavelength,
                    emission_time,
                    source_weight,
                    throughput,
                    status,
                    path_length,
                    *focal,
                    points,
                    *incidence,
                )
            )
    return arrivals
