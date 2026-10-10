"""sim_telarray mount transforms at the ground-to-telescope input boundary."""

from __future__ import annotations

import csv
import math
import os
import random
from dataclasses import dataclass
from pathlib import Path


def _z(vector, angle):
    x, y, z = vector
    c, s = math.cos(angle), math.sin(angle)
    return (c * x - s * y, s * x + c * y, z)


def _y(vector, angle):
    x, y, z = vector
    c, s = math.cos(angle), math.sin(angle)
    return (c * x + s * z, y, -s * x + c * z)


@dataclass(frozen=True)
class TelescopeFrame:
    """N/West/up ground coordinates; model origin is the elevation-axis position."""

    origin_m: tuple[float, float, float]
    az_alt_offset_m: float
    alt_optics_offset_m: float
    azimuth_rad: float
    zenith_rad: float
    known_azimuth_deg: float
    known_zenith_deg: float

    def to_local(self, position, direction):
        p = _z(tuple(a - b for a, b in zip(position, self.origin_m, strict=True)), self.azimuth_rad)
        p = _y((p[0] - self.az_alt_offset_m, p[1], p[2]), -self.zenith_rad)
        p = (p[0] + self.alt_optics_offset_m, p[1], p[2])
        return p, _y(_z(direction, self.azimuth_rad), -self.zenith_rad)

    def to_ground(self, position, direction):
        p = _y((position[0] - self.alt_optics_offset_m, position[1], position[2]), self.zenith_rad)
        p = _z((p[0] + self.az_alt_offset_m, p[1], p[2]), -self.azimuth_rad)
        p = tuple(a + b for a, b in zip(p, self.origin_m, strict=True))
        return p, _z(_y(direction, self.zenith_rad), -self.azimuth_rad)


def resolve_telescope_frame(
    context, azimuth_deg, zenith_deg, *, seed=0, telescope_id=None, pointing_errors=True
):
    """Known errors affect reported pointing; unknown errors affect the actual frame.

    In reference tel_newdir the unknown zenith draw has a minus sign. Keys are
    explicit and reproducible, without claiming the reference's sequential RNG.
    Optics diagnostics disable pointing errors in both backends.
    """
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed < 2**64:
        raise ValueError("pointing seed must be uint64")
    if any(not math.isfinite(v) for v in (azimuth_deg, zenith_deg)):
        raise ValueError("pointing angles must be finite")
    if context["axis_origin_m"] is None:
        raise ValueError("ground-coordinate replay requires an explicit elevation-axis origin")
    if len(context["axis_origin_m"]) != 3 or len(context["axes_offsets_m"]) != 2:
        raise ValueError("mount context requires three origin coordinates and two axis offsets")
    if any(
        not math.isfinite(v)
        for v in (
            *context["axis_origin_m"],
            *context["axes_offsets_m"],
            context["known_error_deg"],
            context["unknown_error_deg"],
        )
    ):
        raise ValueError("mount and pointing parameters must be finite")
    if context["known_error_deg"] < 0 or context["unknown_error_deg"] < 0:
        raise ValueError("pointing error RMS must be nonnegative")
    known = context["known_error_deg"] if pointing_errors else 0
    unknown = context["unknown_error_deg"] if pointing_errors else 0
    telescope_id = context.get("pointing_key", 0) if telescope_id is None else telescope_id
    rng = random.Random(f"obdeect-telescope-pointing-v1:{seed}:{telescope_id}")
    theta = zenith_deg + rng.gauss(0, known)
    phi = azimuth_deg + rng.gauss(0, known)
    actual_phi = phi + rng.gauss(0, unknown)
    actual_theta = theta - rng.gauss(0, unknown)
    return TelescopeFrame(
        tuple(context["axis_origin_m"]),
        *context["axes_offsets_m"],
        math.radians(actual_phi),
        math.radians(actual_theta),
        phi,
        theta,
    )


def single_telescope_identity(source: Path) -> str:
    """Return the sole telescope identity in a replay input, or reject mixed identities."""
    with Path(source).open(newline="") as stream:
        reader = csv.DictReader(stream)
        identities = {row.get("telescope_id") or "0" for row in reader}
    if not identities:
        raise ValueError("ground replay input contains no photons")
    if len(identities) != 1:
        raise ValueError("ground replay invocation requires one consistent telescope identity")
    return identities.pop()


def prepare_ground_photons(
    source: Path,
    output: Path,
    context,
    azimuth_deg,
    zenith_deg,
    *,
    seed=0,
    pointing_errors=True,
    telescope_id: str | None = None,
):
    """Preserve input identities/weights/time; transform at the boundary, never in analysis."""
    source = Path(source)
    output = Path(output)
    if source.resolve() == output.resolve() or (
        source.exists() and output.exists() and os.path.samefile(source, output)
    ):
        raise ValueError("ground photon input and output must be different files")
    with source.open(newline="") as stream, output.open("w", newline="") as destination:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames or []
        required = {"x_m", "y_m", "z_m", "dx", "dy", "dz"}
        if not required.issubset(fields):
            raise ValueError("ground photon input lacks position/direction columns")
        writer = csv.DictWriter(destination, fieldnames=fields)
        writer.writeheader()
        frames = {}
        for row in reader:
            local_context = dict(context)
            if all(name in fields for name in ("telescope_x_m", "telescope_y_m", "telescope_z_m")):
                local_context["axis_origin_m"] = [
                    float(row[f"telescope_{axis}_m"]) for axis in "xyz"
                ]
            row_identity = row.get("telescope_id") or "0"
            if telescope_id is not None and row_identity != telescope_id:
                raise ValueError(
                    "ground replay photon telescope identity changed during preparation"
                )
            identity = row_identity if telescope_id is None else telescope_id
            key = (identity, tuple(local_context["axis_origin_m"] or []))
            if key not in frames:
                frames[key] = resolve_telescope_frame(
                    local_context,
                    azimuth_deg,
                    zenith_deg,
                    seed=seed,
                    telescope_id=identity,
                    pointing_errors=pointing_errors,
                )
            frame = frames[key]
            position = tuple(float(row[f"{axis}_m"]) for axis in "xyz")
            direction = tuple(float(row[f"d{axis}"]) for axis in "xyz")
            if any(not math.isfinite(v) for v in (*position, *direction)):
                raise ValueError("ground photon coordinates must be finite")
            p, d = frame.to_local(position, direction)
            for axis, value, component in zip("xyz", p, d, strict=True):
                row[f"{axis}_m"] = format(value, ".17g")
                row[f"d{axis}"] = format(component, ".17g")
            writer.writerow(row)
