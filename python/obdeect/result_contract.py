"""Backend-neutral optical arrival records.

The contract deliberately stops at the optical detector surface.  It is small
enough for CSV interchange and strict enough for simtools analysis to reject
ambiguous or physically invalid records before calculating PSF observables.
"""

from __future__ import annotations

import csv
import math
from collections.abc import Iterator
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
    run_id: int = 0
    event_id: int = 0
    array_id: int = 0
    telescope_id: int = 0
    bunch_id: int = 0
    arrival_time_ns: float | None = None
    terminal_surface_id: int | None = None
    final_direction: tuple[float, float, float] | None = None
    interaction_surface_ids: tuple[int, ...] | None = None
    response_loss_fraction: float | None = None
    terminal_loss_fraction: float | None = None

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

TERMINAL_STATUSES = frozenset({
    "detected",
    "blocked_camera",
    "blocked_mast",
    "blocked_obscurer",
    "missed_primary",
    "missed_secondary",
    "missed_screen",
    "no_detector",
    "invalid_input",
    "escaped_optical_model",
    "interaction_limit",
    "intersection_failure",
})


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


def iter_arrivals(path: Path) -> Iterator[OpticalArrival]:
    """Read one arrival-v1 batch with unique photon identities.

    A detected photon must include a detector vertex after its launch vertex.
    Zero response at that detector is legitimate and retains ``detected``;
    terminal losses always have zero arriving throughput.
    """

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
        identities: set[tuple[int, ...]] = set()
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
            context_names = ("run_id", "event_id", "array_id", "telescope_id", "bunch_id")
            if any(name in fields for name in context_names) and not all(
                name in fields for name in context_names
            ):
                raise ArrivalContractError(f"{path}:{line}: incomplete photon identity context")
            try:
                context = tuple(int(row[name]) if name in fields else 0 for name in context_names)
            except (TypeError, ValueError) as error:
                raise ArrivalContractError(
                    f"{path}:{line}: invalid photon identity context"
                ) from error
            identity = (*context, photon_id)
            if any(value < 0 for value in context) or identity in identities:
                raise ArrivalContractError(f"{path}:{line}: duplicate photon_id {photon_id}")
            identities.add(identity)
            wavelength = _number(row, "wavelength_nm", path, line)
            source_kind = row["source_kind"]
            if source_kind not in {"star", "illuminator", "laser", "replay"}:
                raise ArrivalContractError(f"{path}:{line}: invalid source_kind")
            status = row["status"]
            if status not in TERMINAL_STATUSES:
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
            if status == "detected" and point_count < 2:
                raise ArrivalContractError(f"{path}:{line}: detected row lacks detector vertex")
            if status != "detected" and throughput != 0:
                raise ArrivalContractError(f"{path}:{line}: lost photon has nonzero throughput")
            if any(value is not None and not 0 <= value <= 90 for value in incidence):
                raise ArrivalContractError(f"{path}:{line}: invalid incidence angle")
            points = tuple(_point(row, index, path, line) for index in range(point_count))
            loss_fields = ("response_loss_fraction", "terminal_loss_fraction")
            losses = (None, None)
            if any(name in fields for name in loss_fields):
                losses = tuple(_number(row, name, path, line) for name in loss_fields)
                if any(value < 0 for value in losses) or not math.isclose(
                    throughput + sum(losses), 1, rel_tol=0, abs_tol=1e-9
                ):
                    raise ArrivalContractError(
                        f"{path}:{line}: optical loss fractions do not close"
                    )
                if status == "detected" and losses[1] != 0:
                    raise ArrivalContractError(f"{path}:{line}: detected photon has terminal loss")
            arrival_time = (
                _number(row, "arrival_time_ns", path, line) if "arrival_time_ns" in fields else None
            )
            terminal_surface = None
            if "terminal_surface_id" in fields:
                try:
                    terminal_surface = int(row["terminal_surface_id"])
                except (TypeError, ValueError) as error:
                    raise ArrivalContractError(
                        f"{path}:{line}: invalid terminal_surface_id"
                    ) from error
                if not 0 <= terminal_surface <= 4294967295:
                    raise ArrivalContractError(f"{path}:{line}: invalid terminal_surface_id")
            direction_names = ("final_dx", "final_dy", "final_dz")
            direction = None
            if any(name in fields for name in direction_names):
                direction = tuple(_number(row, name, path, line) for name in direction_names)
                if status == "detected" and not math.isclose(
                    sum(value * value for value in direction), 1, rel_tol=0, abs_tol=1e-9
                ):
                    raise ArrivalContractError(f"{path}:{line}: non-unit detector direction")
            focal = points[-1] if status == "detected" else (None, None, None)
            surface_ids = None
            if "interaction_surface_ids" in fields:
                try:
                    surface_ids = tuple(
                        int(value) for value in row["interaction_surface_ids"].split(";") if value
                    )
                except (AttributeError, ValueError) as error:
                    raise ArrivalContractError(
                        f"{path}:{line}: invalid interaction surface IDs"
                    ) from error
                if len(surface_ids) != point_count - 1 or any(
                    not 0 <= value <= 4294967295 for value in surface_ids
                ):
                    raise ArrivalContractError(
                        f"{path}:{line}: inconsistent interaction surface IDs"
                    )
            yield OpticalArrival(
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
                *context,
                arrival_time,
                terminal_surface,
                direction,
                surface_ids,
                *losses,
            )


def read_arrivals(path: Path) -> list[OpticalArrival]:
    """Read and validate a complete arrival-v1 batch."""
    return list(iter_arrivals(path))
