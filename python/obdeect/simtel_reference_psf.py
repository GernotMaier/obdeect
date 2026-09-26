"""Recalculate radial PSF references from archived sim_telarray imaging lists."""

import argparse
import bisect
import csv
import gzip
import hashlib
import io
import json
import math
import re
from pathlib import Path

FRACTIONS = (0.5, 0.68, 0.8, 0.9, 0.95, 0.99)
PROFILE_STEPS = 1000
_ROTATION = re.compile(r"^# Camera rotation angle = ([\d.+-]+) deg$")
_ILLUMINATION = re.compile(
    r"^# Telescope \d+ with (\d+) photons .* falling on an area of ([\d.+-]+) m\^2$"
)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def analyse_imaging_list(data: bytes) -> dict:
    """Use the same camera rotation and all focal crossings as simtools PSFImage."""
    angle = None
    launched = None
    illuminated_area = None
    points = []
    for line_number, line in enumerate(data.decode("ascii").splitlines(), 1):
        if match := _ROTATION.fullmatch(line):
            angle = math.radians(float(match[1]))
        elif match := _ILLUMINATION.fullmatch(line):
            launched = int(match[1])
            illuminated_area = float(match[2])
        elif line and not line.startswith("#"):
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f"imaging list line {line_number}: missing focal coordinates")
            x, y = float(fields[2]), float(fields[3])
            if not math.isfinite(x) or not math.isfinite(y):
                raise ValueError(f"imaging list line {line_number}: nonfinite coordinate")
            points.append((x, y))
    if angle is None or launched is None or illuminated_area is None or not points:
        raise ValueError("imaging list lacks rotation, illuminated area, or focal crossings")
    if launched < len(points) or illuminated_area <= 0:
        raise ValueError("imaging list has inconsistent photon count or area")
    cosine, sine = math.cos(angle), math.sin(angle)
    rotated = [(x * cosine - y * sine, y * cosine + x * sine) for x, y in points]
    count = len(rotated)
    centroid_x = math.fsum(x for x, _ in rotated) / count
    centroid_y = math.fsum(y for _, y in rotated) / count
    radii = sorted(math.hypot(x - centroid_x, y - centroid_y) for x, y in rotated)
    return {
        "focal_plane_crossings": count,
        "launched_photons": launched,
        "illuminated_area_m2": illuminated_area,
        "effective_area_m2": illuminated_area * count / launched,
        "camera_rotation_deg": math.degrees(angle),
        "centroid_x_m": centroid_x / 100,
        "centroid_y_m": centroid_y / 100,
        "radii_cm": radii,
    }


def integration_radius(radii_cm: list[float], fraction: float) -> float:
    """Smallest empirical radius containing at least the requested fraction."""
    if not radii_cm or not 0 < fraction <= 1:
        raise ValueError("radii must be nonempty and fraction must be in (0, 1]")
    return radii_cm[math.ceil(fraction * len(radii_cm)) - 1]


def profile_bytes(radii_cm: list[float]) -> bytes:
    """Encode the empirical cumulative distribution at 0.1% increments."""
    handle = io.StringIO(newline="")
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow((
        "integration_radius_cm",
        "enclosed_count",
        "enclosed_fraction",
        "target_containment_fraction",
    ))
    zero_count = bisect.bisect_right(radii_cm, 0.0)
    writer.writerow((0.0, zero_count, zero_count / len(radii_cm), 0.0))
    for step in range(1, PROFILE_STEPS + 1):
        fraction = step / PROFILE_STEPS
        radius = integration_radius(radii_cm, fraction)
        enclosed_count = bisect.bisect_right(radii_cm, radius)
        writer.writerow((radius, enclosed_count, enclosed_count / len(radii_cm), fraction))
    return handle.getvalue().encode("ascii")


def write_profile(path: Path, radii_cm: list[float]) -> None:
    """Write the empirical cumulative distribution at 0.1% increments."""
    path.write_bytes(profile_bytes(radii_cm))


def update_summary(summary_path: Path, root: Path, *, check: bool = False) -> None:
    """Verify archived photon bytes and attach numerical PSF products."""
    summary = json.loads(summary_path.read_text())
    for row in summary["rows"]:
        archive = root / row["imaging_list_archive"]
        compressed = archive.read_bytes()
        if _sha256(compressed) != row["sha256"]["imaging_list_archive"]:
            raise ValueError(f"compressed imaging-list hash mismatch: {archive}")
        raw = gzip.decompress(compressed)
        if _sha256(raw) != row["sha256"]["imaging_list"]:
            raise ValueError(f"imaging-list hash mismatch: {archive}")
        result = analyse_imaging_list(raw)
        for key in ("focal_plane_crossings", "centroid_x_m", "centroid_y_m"):
            if not math.isclose(result[key], row[key], abs_tol=1e-11):
                raise ValueError(f"{archive}: {key} differs from recorded reference")
        if not math.isclose(result["effective_area_m2"], row["effective_area_m2"], rel_tol=1e-12):
            raise ValueError(f"{archive}: effective area differs from recorded reference")
        radii = result.pop("radii_cm")
        integration_radii = {
            str(fraction): integration_radius(radii, fraction) for fraction in FRACTIONS
        }
        exact_d80 = 2 * integration_radii["0.8"] / 100
        profile = root / row["cumulative_profile"]
        expected_profile = profile_bytes(radii)
        if check:
            if (
                row.get("integration_radii_cm") != integration_radii
                or row.get("exact_d80_m") != exact_d80
                or profile.read_bytes() != expected_profile
                or row["sha256"].get("cumulative_profile") != _sha256(expected_profile)
            ):
                raise ValueError(f"radial PSF product differs from archived photons: {profile}")
        else:
            row["integration_radii_cm"] = integration_radii
            row["exact_d80_m"] = exact_d80
            profile.parent.mkdir(parents=True, exist_ok=True)
            profile.write_bytes(expected_profile)
            row["sha256"]["cumulative_profile"] = _sha256(expected_profile)
    if check:
        return
    summary["radial_analysis"] = {
        "centre": "mean of camera-rotation-corrected focal-plane crossings",
        "weights": "one per crossing, including pixel number -1",
        "integration_radius": (
            "smallest empirical centroid-centred circle containing at least the target fraction"
        ),
        "profile": "empirical cumulative distribution sampled at target fractions 0, 0.001, ..., 1",
        "position_unit": "cm in the focal plane",
        "d80_m": "simtools iterative result; exact_d80_m is twice the empirical 80% radius",
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_dir", type=Path)
    parser.add_argument("--check", action="store_true", help="verify without rewriting products")
    args = parser.parse_args()
    update_summary(args.reference_dir / "summary.json", args.reference_dir, check=args.check)


if __name__ == "__main__":
    main()
